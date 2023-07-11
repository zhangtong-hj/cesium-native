#pragma once

#include "Cesium3DTilesSelection/Library.h"
#include "Cesium3DTilesSelection/RasterOverlay.h"
#include "Cesium3DTilesSelection/RasterOverlayTileProvider.h"
#include "Cesium3DTilesSelection/Variogram.h"
#include <CesiumAsync/AsyncSystem.h>
#include <CesiumGeospatial/CartographicPolygon.h>
#include <CesiumGeospatial/Ellipsoid.h>
#include <CesiumGeospatial/Projection.h>

#include <spdlog/fwd.h>

#include <memory>
#include <string>
#include <vector>

namespace Cesium3DTilesSelection {

struct HeatMapDataSource {
  std::vector<glm::dvec3> pos;
  std::vector<glm::u8vec3> col;
  std::vector<float> val;
  float radius;
  VariogramModel vmodel;
  vector3* vecPoint;
  matrix* mVInvt;
};

class CESIUM3DTILESSELECTION_API HeatMapOverlay final : public RasterOverlay {

public:
  HeatMapOverlay(
      const std::string& name,
      const std::vector<CesiumGeospatial::CartographicPolygon>& polygons,
      HeatMapDataSource dataSource,
      bool invertSelection,
      const CesiumGeospatial::Ellipsoid& ellipsoid,
      const CesiumGeospatial::Projection& projection,
      const RasterOverlayOptions& overlayOptions = {});
  virtual ~HeatMapOverlay() override;

  virtual CesiumAsync::Future<CreateTileProviderResult> createTileProvider(
      const CesiumAsync::AsyncSystem& asyncSystem,
      const std::shared_ptr<CesiumAsync::IAssetAccessor>& pAssetAccessor,
      const std::shared_ptr<CreditSystem>& pCreditSystem,
      const std::shared_ptr<IPrepareRendererResources>&
          pPrepareRendererResources,
      const std::shared_ptr<spdlog::logger>& pLogger,
      CesiumUtility::IntrusivePointer<const RasterOverlay> pOwner)
      const override;

  const std::vector<CesiumGeospatial::CartographicPolygon>&
  getPolygons() const noexcept {
    return this->_polygons;
  }
  const HeatMapDataSource& getHeatMapDataSource() { return this->_dataSource; }

  bool getInvertSelection() const noexcept { return this->_invertSelection; }

private:
  std::vector<CesiumGeospatial::CartographicPolygon> _polygons;
  bool _invertSelection;
  CesiumGeospatial::Ellipsoid _ellipsoid;
  CesiumGeospatial::Projection _projection;
  HeatMapDataSource _dataSource;
};
} // namespace Cesium3DTilesSelection
