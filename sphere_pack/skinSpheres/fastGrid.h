#pragma once

#include "eigen/Dense"
#include <vector>


//takes in a tetmesh and answers "which tet is point p in" real fast
//uses a uniform grid


template<typename U1, typename U2>
class FastGrid{

  public:
  FastGrid(const Eigen::DenseBase<U1>& _tetVertices, 
		   const Eigen::DenseBase<U2>& _tetTets,
		   size_t _numSplits);

  template <typename U3, typename U4>
  size_t containingTet(
	  const Eigen::MatrixBase<U3>& point, 
	  U4&& barycentricCoords );


  Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> tetVertices;
  Eigen::Matrix<int, Eigen::Dynamic, 4, Eigen::RowMajor> tetTets;

  size_t numSplits;
  using ScalarType = float;
  using VecType = Eigen::Matrix<ScalarType, 3, 1>;

  VecType minCorner, maxCorner, delta;

  std::vector<std::vector<size_t>> grid;
  //grid[x,y,z] = grid[x*numSlplits^2 + y*numSplits + z]

  

};

template<typename U1, typename U2>
FastGrid<U1, U2> makeFastGrid(
	const Eigen::DenseBase<U1>& tetVertices,
	const Eigen::DenseBase<U2>& tetTets,
	size_t numSplits){
  return FastGrid<U1, U2>(tetVertices, tetTets, numSplits);
}

#include "fastGrid.cpp"
