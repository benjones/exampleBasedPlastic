#pragma once

#include <exampleGraph.h>
#include <Eigen/dense>

struct EGPosition{
  int simplex;
  std::vector<double> coords;
};

//move through an example graph
class EGTraverser{
public:
  EGTraverser(ExampleGraph* _graph = nullptr):
	graph{_graph},
	currentSimplex{-1}
  {}
  

  //update the current timestep
  void traverse();

  
  EGPosition restPosition;
  double restSpringStrength;
  
  EGPosition currentPosition;


  ExampleGraph* graph;
  int currentSimplex;
  
  //sp(i,j) is the simplex adjacent to i to go toward
  //to get to simplex j
  Eigen::MatrixXi shortestPaths;
  //edges between simplices
  std::vector<std::pair<int,int>> edges;
  
  void computeShortestPaths();
  void computeEdges();

  EGPosition interpolatePoints(const EGPosition& first,
							   const EGPosition& second,
							   double alpha);
  double distanceBetweenSimplexPoints(const EGPosition& first,
									  const EGPosition& second);
  
  //find the point on the boundary of first's simplex, and its neighbor
  void getCorrespondingPoint(const EGPosition& first,
							  int neighborSimplex,
							  EGPosition& out1,
							  EGPosition& out2);

  
  static void zeroWeightsSum(std::vector<double>& weights);
  static void zeroWeightsSum(Eigen::VectorXd& weights);


  static EGPosition addBcVector(const EGPosition& start, const std::vector<double>& direction);
  static EGPosition addBcVector(const EGPosition& start, const Eigen::VectorXd& direction);

};

