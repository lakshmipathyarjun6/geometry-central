#pragma once

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

#include <complex>
#include <memory>
#include <tuple>
#include <vector>

namespace geometrycentral {
namespace surface {

// Stateful class. Allows efficient repeated solves

class VectorHeatMethodSolver {

public:
  // === Constructor
  VectorHeatMethodSolver(IntrinsicGeometryInterface& geom, double tCoef = 1.0);


  // === Scalar Extension
  VertexData<double> extendScalar(const std::vector<std::tuple<Vertex, double>>& sources);
  VertexData<double> extendScalar(const std::vector<std::tuple<SurfacePoint, double>>& sources);


  // === Tangent Vector Extension
  VertexData<Vector2> transportTangentVector(Vertex sourceVert, Vector2 sourceVec);
  VertexData<Vector2> transportTangentVectors(const std::vector<std::tuple<Vertex, Vector2>>& sources);
  VertexData<Vector2> transportTangentVectors(const std::vector<std::tuple<SurfacePoint, Vector2>>& sources);


  // === The Logarithmic map
  VertexData<Vector2> computeLogMap(const Vertex& sourceVert, double vertexDistanceShift = 0., double vertAngleRad = 0.,
                                    bool invert = false);
  VertexData<Vector2> computeLogMap(const SurfacePoint& sourceP);

  VertexData<Vector2> computeLogMapIncrementalRadial(const Vertex& sourceVert, bool invert = false);
  VertexData<Vector2> computeLogMapIncrementalHorizontal(const Vertex& sourceVert, double vertAngleRad = 0.,
                                                         bool invert = false);

  VertexData<Vector2> getHorizontalTangentVectors();
  VertexData<Vector2> getRadialTangentVectors();


  // === Options and parameters
  const double tCoef; // the time parameter used for heat flow, measured as time = tCoef * mean_edge_length^2
                      // default: 1.0


private:
  // === Members

  // Basics
  // TODO FIXME this should probably become a manfiold surface mesh, at least for now
  SurfaceMesh& mesh;
  IntrinsicGeometryInterface& geom;

  // Parameters
  double shortTime; // the actual time used for heat flow computed from tCoef

  // Solvers
  std::unique_ptr<PositiveDefiniteSolver<double>> scalarHeatSolver;
  std::unique_ptr<LinearSolver<std::complex<double>>> vectorHeatSolver;
  std::unique_ptr<PositiveDefiniteSolver<double>> poissonSolver;
  SparseMatrix<double> massMat;

  Vector<std::complex<double>> radialSol;
  Vector<std::complex<double>> horizontalSol;

  VertexData<Vector2> radialTangentVecs;
  VertexData<Vector2> horizontalTangentVecs;
  Vector<double> distance;

  // Helpers
  void ensureHaveScalarHeatSolver();
  void ensureHaveVectorHeatSolver();
  void ensureHavePoissonSolver();

  void addVertexOutwardBall(Vertex v, Vector<std::complex<double>>& distGradRHS);
};


} // namespace surface
} // namespace geometrycentral
