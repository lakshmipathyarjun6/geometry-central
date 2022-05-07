#pragma once

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/vector2.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

struct params {
  size_t cp;   // index of closest point on axis
  double dist; // distance away from axis
  std::complex<double> dir;
};

class SurfacePatch {

public:
  // constructors
  SurfacePatch(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry,
               HeatMethodDistanceSolver* distanceHeatSolver, VectorHeatMethodSolver* vectorHeatSolver);

  void createDefaultAxis();

  void get(std::vector<SurfacePoint>& axis, std::vector<SurfacePoint>& boundary);

  void rotate(double angle);

  void transfer(SurfacePatch* target, const Vertex& targetMeshStart, Vector2 initDir);

  void translate(const Vertex& newStartVertex);

private:
  std::complex<double> axisTangent(size_t idx, const std::vector<SurfacePoint>& axis);

  std::vector<params> closestPointAndDirection();

  void computeAxisAnglesAndDistances();

  std::vector<SurfacePoint> connectPointsWithGeodesic(const SurfacePoint& pt1, const SurfacePoint& pt2);

  void constructDenselySampledAxis();

  Vector2 evaluateLogMap(const VertexData<Vector2>& logMap, const SurfacePoint& pt);

  double evaluateVertexDataAtPoint(const VertexData<double>& u, const SurfacePoint& pt);

  size_t indexOfClosestPointOnAxis(double diffusedVal,
                                   const std::vector<std::tuple<SurfacePoint, double>>& zippedDistances);

  Vector2 localDir(const SurfacePoint& pt1, const SurfacePoint& pt2);

  std::unique_ptr<ManifoldSurfaceMesh> m_mesh;
  std::unique_ptr<VertexPositionGeometry> m_geometry;
  std::unique_ptr<HeatMethodDistanceSolver> m_distanceHeatSolver;
  std::unique_ptr<VectorHeatMethodSolver> m_vectorHeatSolver;

  std::vector<SurfacePoint> m_patchAxisSparse;
  std::vector<std::complex<double>> m_patchAxisSparseAngles;
  std::vector<double> m_patchAxisSparseDistances;
  std::vector<size_t> m_patchAxisSparseDenseIdx;

  std::vector<SurfacePoint> m_patchAxisDense;

  std::vector<SurfacePoint> m_patchBoundary;
};
