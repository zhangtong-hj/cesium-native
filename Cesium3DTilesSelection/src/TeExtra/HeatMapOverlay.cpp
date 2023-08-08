#include "Cesium3DTilesSelection/TeExtra/HeatMapOverlay.h"

#include "Cesium3DTilesSelection/BoundingVolume.h"
#include "Cesium3DTilesSelection/RasterOverlayTileProvider.h"
#include "Cesium3DTilesSelection/spdlog-cesium.h"
#include "TileUtilities.h"

#include <CesiumAsync/AsyncSystem.h>
#include <CesiumAsync/IAssetAccessor.h>
#include <CesiumGeospatial/GlobeRectangle.h>
#include <CesiumUtility/IntrusivePointer.h>

#include <memory>
#include <string>
#include <vector>
#include "Cesium3DTilesSelection/Variogram.h"
using namespace CesiumGeometry;
using namespace CesiumGeospatial;
using namespace CesiumUtility;

namespace Cesium3DTilesSelection {
namespace {
bool InsidePolygon(glm::dvec2 point, std::vector<glm::dvec2> polygon) {
  bool inside = false;
  float x = point.x, y = point.y;
  int l = polygon.size();
  int i, j = l - 1;
  for (i = 0; i < l; j = i, i++) {
    float ax = polygon[j].x;
    float ay = polygon[j].y;
    float bx = polygon[i].x;
    float by = polygon[i].y;
    if ((ay >= y && by < y) || (by >= y && ay < y)) {
      inside ^= (ax + (bx - ax) * (y - ay) / (by - ay) > x);
    }
  }
  return inside;
}

void rasterizePolygons(
    LoadedRasterOverlayImage& loaded,
    const CesiumGeospatial::GlobeRectangle& rectangle,
    const glm::dvec2& textureSize,
    const std::vector<CartographicPolygon>& cartographicPolygons,
    const HeatMapDataSource& data,
    const HeatMapGridSetting& grid,
    bool invertSelection) {

  CesiumGltf::ImageCesium& image = loaded.image.emplace();

  std::byte insideColor;
  std::byte outsideColor;
  int t = (int)data.radius;
  insideColor = static_cast<std::byte>(0);
  outsideColor = static_cast<std::byte>(0);

  double minX = data.bnd[0].x, maxX = data.bnd[0].x, minY = data.bnd[0].y,
               maxY = data.bnd[0].y;
  // find max&min rectangle of bound
  for(const auto& point:data.bnd){
    minX = glm::min(minX, point.x);
    maxX = glm::max(maxX, point.x);
    minY = glm::min(minY, point.y);
    maxY = glm::max(maxY, point.y);
  }

  const CesiumGeospatial::GlobeRectangle coverRectangle(minX,minY,maxX,maxY);
  // create a 1x1 mask if the rectangle is completely outside all polygons
  if (!rectangle.computeIntersection(coverRectangle)) {
    loaded.moreDetailAvailable = false;
    image.width = 1;
    image.height = 1;
    image.channels = 4;
    image.bytesPerChannel = 1;
    image.pixelData.resize(4, outsideColor);
    return;
  }

  const double rectangleWidth = rectangle.computeWidth();
  const double rectangleHeight = rectangle.computeHeight();
  // create source image
  loaded.moreDetailAvailable = true;
  image.width = int32_t(glm::round(textureSize.x));
  image.height = int32_t(glm::round(textureSize.y));
  image.channels = 4;
  image.bytesPerChannel = 1;
  image.pixelData.resize(
      size_t(image.width * image.height * image.channels),
      outsideColor);

  // get grid params
  const double MIN_GRID_SIZE=0;
  double offsetX = 0, offsetY = 0, sizeX = grid.size.x, sizeY = grid.size.y,
         sampleNumX = grid.smplFreq.x, sampleNumY = grid.smplFreq.y,
         lineWidth = grid.lineW;
  int gridWidth=0,gridHeight=0,sampleMode=grid.smplMode;
  bool useGrid = false;

  std::vector<float> vecHeatValuePxl, vecHeatValue, vecOpacityValue, sampleNum;
  std::vector<bool> isLine, inPoly;
  isLine.resize(image.width * image.height, 0);
  inPoly.resize(image.width * image.height, 0);

  if (sizeX > MIN_GRID_SIZE && sizeY > MIN_GRID_SIZE){
    gridWidth = (int)(rectangleWidth / sizeX) + 2;
    gridHeight =  (int)(rectangleHeight / sizeY) + 2;
    useGrid = true;
    vecHeatValue.resize(gridWidth * gridHeight, 0);
    sampleNum.resize(gridWidth * gridHeight, 0);
    vecOpacityValue.resize(gridWidth * gridHeight, 0);
  } else {
    vecHeatValuePxl.resize(image.width * image.height, 0);
  }

  
  std::vector<CesiumGeospatial::GlobeRectangle> vecHpBoundingVolume;
  double radius = data.radius;
  size_t width = size_t(image.width);
  size_t height = size_t(image.height);

  if(useGrid){
    for (size_t j = 0; j < gridHeight; ++j)
      for (size_t i = 0; i < gridWidth; ++i)
        for (size_t sx = 0; sx < sampleNumX; ++sx)
          for (size_t sy = 0; sy < sampleNumY; ++sy) {
            const double pixelY = ((int)(rectangle.getSouth() / sizeY) + j + (sy + 0.5) / sampleNumY) * sizeY;
            const double pixelX = ((int)(rectangle.getWest() / sizeX) + i + (sx + 0.5) / sampleNumX) * sizeX;
            const glm::dvec2 v(pixelX, pixelY);
            double dSumDis = 0;
            float heatValue = 0;
            double gridN = ((int)(rectangle.getSouth() / sizeY) + j +1) *
                           sizeY,
                   gridS=((int)(rectangle.getSouth() / sizeY) + j) * sizeY,
                   gridW = ((int)(rectangle.getWest() / sizeX) + i) * sizeX,
                   gridE = ((int)(rectangle.getWest() / sizeX) + i+1) * sizeX;
            for (const auto& samplePoint : data.pos) {
              float pointValue = samplePoint.z;
              if (samplePoint.x > gridW && samplePoint.x <= gridE &&
                  samplePoint.y > gridS && samplePoint.y <= gridN) {
                if (sampleMode == 0) {
                  if (vecHeatValue[gridWidth * j + i] <= 0) {
                    vecHeatValue[gridWidth * j + i] = pointValue;
                  } else {
                    vecHeatValue[gridWidth * j + i] = glm::min(
                        vecHeatValue[gridWidth * j + i], pointValue);
                  }
                } else if (sampleMode == 1) {
                  vecHeatValue[gridWidth * j + i] =
                      glm::max(vecHeatValue[gridWidth * j + i], pointValue);
                }
              }
              glm::dvec2 sample(samplePoint.x, samplePoint.y);
              dSumDis += 1.0 / distance(sample, v);
            }
            for (const auto& samplePoint : data.pos) {
              glm::dvec2 sample(samplePoint.x, samplePoint.y);
              heatValue += 1.0 / distance(sample, v) * samplePoint.z /
                                dSumDis;
            }
            if (sampleMode == 0) {
              if (vecHeatValue[gridWidth * j + i] <= 0) {
                vecHeatValue[gridWidth * j + i] = heatValue;
              } else {
                vecHeatValue[gridWidth * j + i] =
                    glm::min(vecHeatValue[gridWidth * j + i], heatValue);
              }
            } else if (sampleMode == 1) {
              vecHeatValue[gridWidth * j + i] =
                  glm::max(vecHeatValue[gridWidth * j + i], heatValue);
            } else {
              vecHeatValue[gridWidth * j + i] +=
                  heatValue / sampleNumX / sampleNumY;
            }
          }
  }else{
    for (size_t j = 0; j < height; ++j) {
      const double pixelY =
          rectangle.getSouth() +
          rectangleHeight * (1.0 - (double(j) + 0.5) / double(height));
      for (size_t i = 0; i < width; ++i) {
        const double pixelX = rectangle.getWest() + rectangleWidth *
                                                        (double(i) + 0.5) /
                                                        double(width);
        const glm::dvec2 v(pixelX, pixelY);

        if (InsidePolygon(v, data.bnd)) {
          inPoly[image.width * j + i] = true;
          double dSumDis = 0;
          for (const auto& samplePoint : data.pos) {
            glm::dvec2 sample(samplePoint.x, samplePoint.y);
            dSumDis += 1.0 / distance(sample, v);
          }
          for (const auto& samplePoint : data.pos) {
            glm::dvec2 sample(samplePoint.x, samplePoint.y);
            vecHeatValuePxl[image.width * j + i] +=
                1.0 / distance(sample, v) * samplePoint.z / dSumDis;
          }
        }
      }
    }
  }

  for (size_t j = 0; j < height; ++j) {
    for (size_t i = 0; i < width; ++i) {
      float fPointVal = 0;
      if (useGrid) {
        const double pixelY =
            rectangle.getSouth() +
            rectangleHeight * (1.0 - (double(j) + 0.5) / double(height));
        const double pixelX = rectangle.getWest() + rectangleWidth *
                                                        (double(i) + 0.5) /
                                                        double(width);
        const glm::dvec2 v(pixelX, pixelY);
        if (!InsidePolygon(v, data.bnd)) {
          image.pixelData[(image.width * j + i) * 4 + 0] = (std::byte)0;
          image.pixelData[(image.width * j + i) * 4 + 1] = (std::byte)0;
          image.pixelData[(image.width * j + i) * 4 + 2] = (std::byte)0;
          image.pixelData[(image.width * j + i) * 4 + 3] = (std::byte)0;
          continue;
        }
        const double pixelLeft =
            rectangle.getWest() + rectangleWidth * double(i) / double(width),
                     pixelRight = rectangle.getWest() + rectangleWidth *
                                                          (double(i) + 1) /
                                                          double(width),
                     pixelBottom =
                         rectangle.getSouth() +
                              rectangleHeight *
                                  (1.0 - (double(j) + 1) / double(height)),
                     pixelTop =
                         rectangle.getSouth() +
                              rectangleHeight *
                                  (1.0 - double(j) / double(height));
        const double dislineX = std::fmod(pixelX, sizeX),
                     dislineY = std::fmod(pixelY, sizeY);
        if (dislineX < lineWidth || dislineY < lineWidth ||
            ((int)(pixelLeft / sizeX) < (int)(pixelRight / sizeX)) ||
            ((int)(pixelBottom / sizeY) < (int)(pixelTop / sizeY))) {
          image.pixelData[(image.width * j + i) * 4 + 0] = (std::byte)128;
          image.pixelData[(image.width * j + i) * 4 + 1] = (std::byte)128;
          image.pixelData[(image.width * j + i) * 4 + 2] = (std::byte)128;
          image.pixelData[(image.width * j + i) * 4 + 3] = (std::byte)255;
          continue;
        }
        int gridY = (int)(pixelY / sizeY) - (int)(rectangle.getSouth() / sizeY),
            gridX = (int)(pixelX / sizeX) - (int)(rectangle.getWest() / sizeX);
        fPointVal = vecHeatValue[gridWidth * gridY + gridX];
      }else if (!inPoly[image.width * j + i]) {
        image.pixelData[(image.width * j + i) * 4 + 0] = (std::byte)0;
        image.pixelData[(image.width * j + i) * 4 + 1] = (std::byte)0;
        image.pixelData[(image.width * j + i) * 4 + 2] = (std::byte)0;
        image.pixelData[(image.width * j + i) * 4 + 3] = (std::byte)0;
        continue;
      } else {
        fPointVal = vecHeatValuePxl[image.width * j + i];
      }

      int k = 0;
      for (k = 0; k < data.val.size()-1 && fPointVal >= data.val[k + 1]; k++) {
      }
      if (k < 1) {
        image.pixelData[(image.width * j + i) * 4 + 0] = (std::byte)0;
        image.pixelData[(image.width * j + i) * 4 + 1] = (std::byte)0;
        image.pixelData[(image.width * j + i) * 4 + 2] = (std::byte)0;
        image.pixelData[(image.width * j + i) * 4 + 3] = (std::byte)0;
      }else if (k == data.val.size() - 1) {
        image.pixelData[(image.width * j + i) * 4 + 0] = (std::byte)data.col[k].x;
        image.pixelData[(image.width * j + i) * 4 + 1] =
            (std::byte)data.col[k].y;
        image.pixelData[(image.width * j + i) * 4 + 2] =
            (std::byte)data.col[k].z;
        image.pixelData[(image.width * j + i) * 4 + 3] = (std::byte)255;
      } else {
        glm::vec3 low = (glm::vec3)data.col[k];
        glm::vec3 high = (glm::vec3)data.col[k + 1];
        float a = (fPointVal - data.val[k]) / (data.val[k + 1] - data.val[k]);
        glm::vec3 mixed = (1 - a) * low + a * high;
        glm::uvec3 umixed = (glm::uvec3)((glm::ivec3)mixed);
        image.pixelData[(image.width * j + i) * 4 + 0] = (std::byte)umixed.x;
        image.pixelData[(image.width * j + i) * 4 + 1] = (std::byte)umixed.y;
        image.pixelData[(image.width * j + i) * 4 + 2] = (std::byte)umixed.z;
        image.pixelData[(image.width * j + i) * 4 + 3] = (std::byte)255;
      }
    }
  }
}

  //Make up a heat map;
  //1)calculate the position and radius of each heat point;
  //2)calculate the numerial value of pixel value;
  //  Test each heat point's bounding volume in/out the rectangle;
  //      For the heat point int the rectangle;
  //          Add up its heat value in the pixel;
  //3)Turn the numerical value to the pixel value;
  //
} // namespace

class CESIUM3DTILESSELECTION_API HeatMapTileProvider final
    : public RasterOverlayTileProvider {

private:
  std::vector<CartographicPolygon> _polygons;
  HeatMapDataSource _dataSource;
  HeatMapGridSetting _gridSetting;
  bool _invertSelection;

public:
  HeatMapTileProvider(
      const IntrusivePointer<const RasterOverlay>& pOwner,
      const CesiumAsync::AsyncSystem& asyncSystem,
      const std::shared_ptr<CesiumAsync::IAssetAccessor>& pAssetAccessor,
      const std::shared_ptr<IPrepareRendererResources>&
          pPrepareRendererResources,
      const std::shared_ptr<spdlog::logger>& pLogger,
      const CesiumGeospatial::Projection& projection,
      const std::vector<CartographicPolygon>& polygons,
      HeatMapDataSource dataSource,
      HeatMapGridSetting gridSetting,
      bool invertSelection)
      : RasterOverlayTileProvider(
            pOwner,
            asyncSystem,
            pAssetAccessor,
            std::nullopt,
            pPrepareRendererResources,
            pLogger,
            projection,
            // computeCoverageRectangle(projection, polygons)),
            projectRectangleSimple(
                projection,
                CesiumGeospatial::GlobeRectangle(
                    -CesiumUtility::Math::OnePi,
                    -CesiumUtility::Math::PiOverTwo,
                    CesiumUtility::Math::OnePi,
                    CesiumUtility::Math::PiOverTwo))),
        _polygons(polygons),
        _dataSource(dataSource),
        _gridSetting(gridSetting),
        _invertSelection(invertSelection) {}

  virtual CesiumAsync::Future<LoadedRasterOverlayImage>
  loadTileImage(RasterOverlayTile& overlayTile) override {
    // Choose the texture size according to the geometry screen size and raster
    // SSE, but no larger than the maximum texture size.
    const RasterOverlayOptions& options = this->getOwner().getOptions();
    glm::dvec2 textureSize = glm::min(
        overlayTile.getTargetScreenPixels() / options.maximumScreenSpaceError,
        glm::dvec2(options.maximumTextureSize));

    return this->getAsyncSystem().runInWorkerThread(
        [&polygons = this->_polygons,
         invertSelection = this->_invertSelection,
         projection = this->getProjection(),
         rectangle = overlayTile.getRectangle(),
         dataSource = this->_dataSource,
         gridSetting = this->_gridSetting,
         textureSize]() -> LoadedRasterOverlayImage {
          const CesiumGeospatial::GlobeRectangle tileRectangle =
              CesiumGeospatial::unprojectRectangleSimple(projection, rectangle);

          LoadedRasterOverlayImage result;
          result.rectangle = rectangle;

          rasterizePolygons(
              result,
              tileRectangle,
              textureSize,
              polygons,
              dataSource,
              gridSetting,
              invertSelection);

          return result;
        });
  }
};

HeatMapOverlay::HeatMapOverlay(
    const std::string& name,
    const std::vector<CesiumGeospatial::CartographicPolygon>& polygons,
    HeatMapDataSource dataSource,
    HeatMapGridSetting gridSetting,
    bool invertSelection,
    const CesiumGeospatial::Ellipsoid& ellipsoid,
    const CesiumGeospatial::Projection& projection,
    const RasterOverlayOptions& overlayOptions)
    : RasterOverlay(name, overlayOptions),
      _polygons(polygons),
      _dataSource(dataSource),
      _gridSetting(gridSetting),
      _invertSelection(invertSelection),
      _ellipsoid(ellipsoid),
      _projection(projection) {}

HeatMapOverlay::~HeatMapOverlay() {}

CesiumAsync::Future<RasterOverlay::CreateTileProviderResult>
HeatMapOverlay::createTileProvider(
    const CesiumAsync::AsyncSystem& asyncSystem,
    const std::shared_ptr<CesiumAsync::IAssetAccessor>& pAssetAccessor,
    const std::shared_ptr<CreditSystem>& /*pCreditSystem*/,
    const std::shared_ptr<IPrepareRendererResources>& pPrepareRendererResources,
    const std::shared_ptr<spdlog::logger>& pLogger,
    CesiumUtility::IntrusivePointer<const RasterOverlay> pOwner) const {

  pOwner = pOwner ? pOwner : this;

  return asyncSystem.createResolvedFuture<CreateTileProviderResult>(
      IntrusivePointer<RasterOverlayTileProvider>(
          new HeatMapTileProvider(
              pOwner,
              asyncSystem,
              pAssetAccessor,
              pPrepareRendererResources,
              pLogger,
              this->_projection,
              this->_polygons,
              this->_dataSource,
              this->_gridSetting,
              this->_invertSelection)));
}

} // namespace Cesium3DTilesSelection
