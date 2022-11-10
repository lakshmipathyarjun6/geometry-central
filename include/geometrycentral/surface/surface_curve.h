#pragma once

#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

struct CurvePointParams {
  double nextPointDistance;                // distance to next point
  std::complex<double> nextPointDirection; // outgoing direction relative to incoming direction
};

class SurfaceCurve {

public:
  SurfaceCurve(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry, GeodesicAlgorithmExact* mmpSolver);

  void createFromPoints(std::vector<SurfacePoint>& curvePoints);

  void deformAngle(int pointIndex, std::complex<double> newDir);

  void deformDistance(int pointIndex, double newDist);

  void get(std::vector<SurfacePoint>& points);

  void getParameterized(std::vector<CurvePointParams>& parameterizedPoints);

  std::complex<double> getDirAtIndex(int index);

  Vector2 getStartDir();

  void rotate(Vector2 newDir);

  void saveStartParams();

  void setParameterized(std::vector<CurvePointParams>& parameterizedPoints);

  void setPoints(const std::vector<SurfacePoint>& points);

  void transfer(SurfaceCurve* target, const SurfacePoint& targetMeshStart, const SurfacePoint& targetMeshDirEndpoint);

  void translate(const SurfacePoint& newStartPoint);

private:
  bool approxEqual(const SurfacePoint& pA, const SurfacePoint& pB);

  void computeAnglesAndDistances();

  void computeInitialDirection();

  std::vector<SurfacePoint> connectPointsWithGeodesic(const SurfacePoint& pt1, const SurfacePoint& pt2,
                                                      double& distance);

  Vector2 inTangentBasis(const SurfacePoint& pA, const SurfacePoint& pB, const SurfacePoint& p);

  Vector2 localDir(const SurfacePoint& pt1, const SurfacePoint& pt2);

  Vector2 parallelTransport(const SurfacePoint& startPoint, const SurfacePoint& endpoint, double initAngle);

  std::vector<SurfacePoint> pruneApproxEqualEntries(const std::vector<SurfacePoint>& source);

  void recompute();

  SurfacePoint m_startPoint;
  Vector2 m_initDir;

  SurfacePoint m_savedStartPoint;
  Vector2 m_savedInitDir;

  std::vector<SurfacePoint> m_points;
  std::vector<std::complex<double>> m_angles;
  std::vector<double> m_distances;

  std::unique_ptr<ManifoldSurfaceMesh> m_mesh;
  std::unique_ptr<VertexPositionGeometry> m_geometry;
  std::unique_ptr<GeodesicAlgorithmExact> m_mmpSolver;
};
