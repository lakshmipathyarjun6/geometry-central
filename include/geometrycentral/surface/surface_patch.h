#pragma once

#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_curve.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/vector2.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

struct PatchPointParams {
  size_t closestAxisPoint;                 // index of closest point on axis
  double axisPointDistance;                // distance away from axis point
  std::complex<double> axisPointDirection; // outgoing direction from axis point
};

class SurfacePatch {

public:
  // constructors
  SurfacePatch(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry, GeodesicAlgorithmExact* mmpSolver,
               VectorHeatMethodSolver* vectorHeatSolver);

  void createCustomAxis(std::vector<SurfacePoint>& axisPoints);

  void createDefaultAxis();

  void deformAxis(int index, std::complex<double> newDir);

  void getAxis(std::vector<SurfacePoint>& axis);

  void getPoints(std::vector<SurfacePoint>& points);

  void getParameterizedAxis(std::vector<CurvePointParams>& parameterizedAxis);

  void getParameterizedPoints(std::vector<PatchPointParams>& parameterizedPoints);

  std::complex<double> getDirAtAxisIndex(int index);

  Vector2 getInitDir();

  double getPatchSpreadCoefficient();

  void parameterizePatch(std::map<size_t, size_t> customClosestPointBinding = {});

  void reconstructPatch();

  void rotateAxis(Vector2 newDir);

  void saveAxisStartPointAndDirection();

  void setPatchSpreadCoefficient(double patchSpreadCoefficient);

  void transfer(SurfacePatch* target, const SurfacePoint& targetMeshStart, const SurfacePoint& targetMeshDirEndpoint);

  void translate(const SurfacePoint& newStartPoint);

  void setAxisPoints(std::vector<SurfacePoint>& axisPoints);

  void setPatchPoints(const std::vector<SurfacePoint>& points);

  void setParameterizedAxis(std::vector<CurvePointParams>& parameterizedAxisPoints);

  void setParameterizedPatch(std::vector<PatchPointParams>& parameterizedPatchPoints);

private:
  bool approxEqual(const SurfacePoint& pA, const SurfacePoint& pB);

  std::vector<SurfacePoint> connectPointsWithGeodesic(const SurfacePoint& pt1, const SurfacePoint& pt2,
                                                      double& distance);

  size_t indexOfClosestPointOnAxis(double diffusedVal,
                                   const std::vector<std::tuple<SurfacePoint, double>>& zippedDistances);

  Vector2 localDir(const SurfacePoint& pt1, const SurfacePoint& pt2);

  std::vector<SurfacePoint> pruneApproxEqualEntries(const std::vector<SurfacePoint>& source);

  std::unique_ptr<ManifoldSurfaceMesh> m_mesh;
  std::unique_ptr<VertexPositionGeometry> m_geometry;
  std::unique_ptr<VectorHeatMethodSolver> m_vectorHeatSolver;
  std::unique_ptr<GeodesicAlgorithmExact> m_mmpSolver;

  SurfaceCurve* m_axis;

  std::vector<PatchPointParams> m_parameterizedPoints;
  std::vector<SurfacePoint> m_points;

  double m_patchSpreadCoefficient;
};
