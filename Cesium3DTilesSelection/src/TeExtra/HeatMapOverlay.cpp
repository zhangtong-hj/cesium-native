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
void rasterizePolygons(
    LoadedRasterOverlayImage& loaded,
    const CesiumGeospatial::GlobeRectangle& rectangle,
    const glm::dvec2& textureSize,
    const std::vector<CartographicPolygon>& cartographicPolygons,
    const HeatMapDataSource& data,
    bool invertSelection) {

  CesiumGltf::ImageCesium& image = loaded.image.emplace();

  std::byte insideColor;
  std::byte outsideColor;
  int t = (int)data.radius;
  insideColor = static_cast<std::byte>(0);
  outsideColor = static_cast<std::byte>(0);

  // create a 1x1 mask if the rectangle is completely inside a polygon
  // if (Cesium3DTilesSelection::CesiumImpl::withinPolygons(
  //         rectangle,
  //         cartographicPolygons)) {
  //   loaded.moreDetailAvailable = false;
  //   image.width = 1;
  //   image.height = 1;
  //   image.channels = 1;
  //   image.bytesPerChannel = 1;
  //   image.pixelData.resize(1, insideColor);
  //   return;
  // }

  bool completelyOutsidePolygons = true;
  for (const CartographicPolygon& selection : cartographicPolygons) {
    const std::optional<CesiumGeospatial::GlobeRectangle>& boundingRectangle =
        selection.getBoundingRectangle();

    if (boundingRectangle &&
        rectangle.computeIntersection(*boundingRectangle)) {
      completelyOutsidePolygons = false;
      break;
    }
  }

  // create a 1x1 mask if the rectangle is completely outside all polygons
  if (completelyOutsidePolygons) {
    loaded.moreDetailAvailable = false;
    image.width = 1;
    image.height = 1;
    image.channels = 1;
    image.bytesPerChannel = 1;
    image.pixelData.resize(1, outsideColor);
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


  std::vector<float> vecHeatValue,vecOpacityValue;
  vecHeatValue.resize(image.width * image.height, 0);
  vecOpacityValue.resize(image.width * image.height, 0);
  std::vector<CesiumGeospatial::GlobeRectangle> vecHpBoundingVolume;
  double radius = data.radius;
  size_t width = size_t(image.width);
  size_t height = size_t(image.height);
  if (data.pos.size() > 0) {
    vecHpBoundingVolume.push_back(CesiumGeospatial::GlobeRectangle(
        data.pos[0].x - radius,
        data.pos[0].y - radius,
        data.pos[0].x + radius,
        data.pos[0].y + radius));
    for (auto pos : data.pos) {
      vecHpBoundingVolume[0] = vecHpBoundingVolume[0].computeUnion(CesiumGeospatial::GlobeRectangle(
              pos.x - radius,
              pos.y - radius,
              pos.x + radius,
              pos.y + radius));
    }
    for (int idx = 0; idx < vecHpBoundingVolume.size(); idx++) {
      if (!rectangle.computeIntersection(vecHpBoundingVolume[idx]))
        continue;
      int N = data.pos.size();
      glm::dvec2 center(data.pos[idx].x, data.pos[idx].y);
      for (size_t j = 0; j < height; ++j) {
        const double pixelY =
            rectangle.getSouth() +
            rectangleHeight * (1.0 - (double(j) + 0.5) / double(height));
        for (size_t i = 0; i < width; ++i) {
          const double pixelX = rectangle.getWest() + rectangleWidth *
                                                          (double(i) + 0.5) /
                                                          double(width);
          const glm::dvec2 v(pixelX, pixelY);
          const Cartographic cartPoint(pixelX, pixelY);
          if (!vecHpBoundingVolume[idx].contains(cartPoint))
            continue;
          // vector3 point(pixelX, pixelY, 0);
          // matrix D = matrix(N, 1);
          //
          // // Calculate Variogram Model for point i
          // matrix tmp =
          //     calculateVariogram(point, data.vecPoint, N, D, data.vmodel);
          // matrix W = data.mVInvt->I->multiply(D);
          // point.z = 0;
          // for (int n = 0; n < N; ++n) {
          //   point.z += W(n, 0) * data.vecPoint[n].z;
          // }
          //
          // vecHeatValue[image.width * j + i] = point.z;
          double dSumDis = 0;
          for (auto samplePoint :  data.pos) {
            glm::dvec2 sample(samplePoint.x, samplePoint.y);
            dSumDis += 1.0/distance(sample, v);
          }
          for (auto samplePoint : data.pos) {
            glm::dvec2 sample(samplePoint.x, samplePoint.y);
            vecHeatValue[image.width * j + i] +=
                1.0 / distance(sample, v) * samplePoint.z / dSumDis;
          }
          
        }
      }
    }
    for (size_t j = 0; j < height; ++j) {
      for (size_t i = 0; i < width; ++i) {
        int k = 0;
        float fPointVal = vecHeatValue[image.width * j + i];
        for (k = 0; k < data.val.size()-1 && fPointVal > data.val[k + 1]; k++) {
        }
        if (k == data.val.size() - 1) {
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

  //
  // for (int pixelX = 0; pixelX < image.width; pixelX++) {
  //   for (int pixelY = 0; pixelY < image.height; pixelY++) {
  //     int32_t t =glm::round(
  //         vecHeatValue[image.width * pixelY + pixelX] / 100 * 255);
  //     image.pixelData[(image.width * pixelY + pixelX) * 4+3] = (std::byte)t;
  //     if (t>0)
  //       image.pixelData[(image.width * pixelY + pixelX) * 4] = static_cast<std::byte>(0xff);
  //     
  //   }
  // }
  // TODO: this is naive approach, use line-triangle
  // intersections to rasterize one row at a time
  // NOTE: also completely ignores antimeridian (really these
  // calculations should be normalized to the first vertex)
  //for (const CartographicPolygon& polygon : cartographicPolygons) {
  //  const std::vector<glm::dvec2>& vertices = polygon.getVertices();
  //  const std::vector<uint32_t>& indices = polygon.getIndices();
  //  for (size_t triangle = 0; triangle < indices.size() / 3; ++triangle) {
  //    const glm::dvec2& a = vertices[indices[3 * triangle]];
  //    const glm::dvec2& b = vertices[indices[3 * triangle + 1]];
  //    const glm::dvec2& c = vertices[indices[3 * triangle + 2]];

      // TODO: deal with the corner cases here
  //    const double minX = glm::min(a.x, glm::min(b.x, c.x));
  //    const double minY = glm::min(a.y, glm::min(b.y, c.y));
  //    const double maxX = glm::max(a.x, glm::max(b.x, c.x));
  //    const double maxY = glm::max(a.y, glm::max(b.y, c.y));

  //    const CesiumGeospatial::GlobeRectangle triangleBounds(
  //        minX,
  //        minY,
  //        maxX,
  //        maxY);

      // skip this triangle if it is entirely outside the tile bounds
  //    if (!rectangle.computeIntersection(triangleBounds)) {
  //      continue;
  //    }

  //    const glm::dvec2 ab = b - a;
  //    const glm::dvec2 ab_perp(-ab.y, ab.x);
  //    const glm::dvec2 bc = c - b;
  //    const glm::dvec2 bc_perp(-bc.y, bc.x);
  //    const glm::dvec2 ca = a - c;
  //    const glm::dvec2 ca_perp(-ca.y, ca.x);

  //    size_t width = size_t(image.width);
  //    size_t height = size_t(image.height);

  //    for (size_t j = 0; j < height; ++j) {
  //      const double pixelY =
  //          rectangle.getSouth() +
  //          rectangleHeight * (1.0 - (double(j) + 0.5) / double(height));
  //      for (size_t i = 0; i < width; ++i) {
  //        const double pixelX = rectangle.getWest() + rectangleWidth *
  //                                                        (double(i) + 0.5) /
  //                                                        double(width);
  //        const glm::dvec2 v(pixelX, pixelY);

  //        const glm::dvec2 av = v - a;
  //        const glm::dvec2 cv = v - c;

  //        const double v_proj_ab_perp = glm::dot(av, ab_perp);
  //        const double v_proj_bc_perp = glm::dot(cv, bc_perp);
  //        const double v_proj_ca_perp = glm::dot(cv, ca_perp);

          // will determine in or out, irrespective of winding
  //        if ((v_proj_ab_perp >= 0.0 && v_proj_ca_perp >= 0.0 &&
  //             v_proj_bc_perp >= 0.0) ||
  //            (v_proj_ab_perp <= 0.0 && v_proj_ca_perp <= 0.0 &&
  //             v_proj_bc_perp <= 0.0)) {
  //         image.pixelData[(width * j + i)*4] = insideColor;
  //        }
  //      }
  //    }
  //  }
  //}
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
              invertSelection);

          return result;
        });
  }
};

HeatMapOverlay::HeatMapOverlay(
    const std::string& name,
    const std::vector<CesiumGeospatial::CartographicPolygon>& polygons,
    HeatMapDataSource dataSource,
    bool invertSelection,
    const CesiumGeospatial::Ellipsoid& ellipsoid,
    const CesiumGeospatial::Projection& projection,
    const RasterOverlayOptions& overlayOptions)
    : RasterOverlay(name, overlayOptions),
      _polygons(polygons),
      _dataSource(dataSource),
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
              this->_invertSelection)));
}

} // namespace Cesium3DTilesSelection
