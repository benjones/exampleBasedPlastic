#include "fastGrid.h"
#include <cassert>
#include <iostream>

template <typename U1, typename U2>
FastGrid<U1, U2>::FastGrid
(const Eigen::DenseBase<U1>& _tetVertices,
	const Eigen::DenseBase<U2>& _tetTets,
	size_t _numSplits)
  :
  tetVertices{_tetVertices},
  tetTets{_tetTets},
  numSplits{_numSplits},
  grid(numSplits*numSplits*numSplits){
	
	//compute max and min
	minCorner = tetVertices.colwise().minCoeff();
	maxCorner = tetVertices.colwise().maxCoeff();

	delta = (1.0/numSplits)*(maxCorner - minCorner);
	
	minCorner -= 0.5*delta;
	maxCorner += 0.5*delta;

	delta = (1.0/numSplits)*(maxCorner - minCorner);

	Eigen::Matrix<ScalarType,4,3> tet;
	for(auto tInd = 0; tInd < tetTets.rows(); ++tInd){
	  for(auto i : {0,1,2,3}){
		tet.row(i) = tetVertices.row(tetTets(tInd, i));
	  }
	  
	  VecType tetMin = tet.colwise().minCoeff();
	  tetMin.array() -= 0.01;
	  VecType tetMax = tet.colwise().maxCoeff();
	  tetMax.array() += 0.01; //fudge factor
	  
	  VecType nCells = (tetMax - tetMin);
	  nCells(0) /= delta(0);
	  nCells(1) /= delta(1);
	  nCells(2) /= delta(2);

	  VecType start = (tetMin - minCorner);
	  start(0) /= delta(0);
	  start(1) /= delta(1);
	  start(2) /= delta(2);


	  size_t count = 0;
	  for(int x = floor(start[0]); 
		  x <= floor(start[0]) + ceil(nCells[0]) && 
			x < numSplits; 
		  ++x){
		
		for(int y = floor(start[1]); 
			y <= floor(start[1]) + ceil(nCells[1]) &&
			  y < numSplits;
			++y){
		  for(int z = floor(start[2]);
			  z <= floor(start[2]) + ceil(nCells[2]) &&
				z < numSplits;
			  ++z){
			
			grid[x*numSplits*numSplits + y*numSplits + z].push_back(tInd);
			count++;
		  }
		}
	  }

	}
  }

template <typename U1, typename U2>
template <typename U3, typename U4>
size_t FastGrid<U1, U2>::containingTet(
	const Eigen::MatrixBase<U3>& point,
	U4&& barycentricCoords){
  
  
  size_t x = (point[0] - minCorner[0])/delta[0];
  size_t y = (point[1] - minCorner[1])/delta[1];
  size_t z = (point[2] - minCorner[2])/delta[2];
  
  int bestTet;
  ScalarType bestBcError = std::numeric_limits<ScalarType>::infinity();
  VecType bestBcs;

  for(auto & tInd : grid[x*numSplits*numSplits + y*numSplits + z]){

	Eigen::Matrix<ScalarType, 3, 3> T;
	T.col(0) = (tetVertices.row(tetTets(tInd, 0)) - 
				tetVertices.row(tetTets(tInd, 3))).transpose();
	T.col(1) = (tetVertices.row(tetTets(tInd, 1)) -
				tetVertices.row(tetTets(tInd, 3))).transpose();
	T.col(2) = (tetVertices.row(tetTets(tInd, 2)) -
				tetVertices.row(tetTets(tInd, 3))).transpose();
	
	VecType bcReduced = T.inverse()*
	  (point - 
		  tetVertices.row(tetTets(tInd, 3)).transpose());

	if(bcReduced(0) >= -0.01 &&
	   bcReduced(1) >= -0.01 &&
	   bcReduced(2) >= -0.01 &&
	   bcReduced.sum() <= 1.03){
	  barycentricCoords(0) = bcReduced(0);
	  barycentricCoords(1) = bcReduced(1);
	  barycentricCoords(2) = bcReduced(2);
	  barycentricCoords(3) = 1.0 - bcReduced.sum();
	  return tInd;
	} else {
	  auto mc = bcReduced.minCoeff();
	  auto s = bcReduced.sum();
	  auto thisError = std::max<ScalarType>(
		  mc < 0 ? fabs(mc) : 0, 
		  (s > 1) ? s - 1 : 0);

	  assert(thisError > 0);

	  if(thisError < bestBcError){
		bestBcError = thisError;
		bestTet = tInd;
		bestBcs = bcReduced;
	  }
							 
	}
  }
		
  std::cout << "checked : " << grid[x*numSplits*numSplits + y*numSplits + z].size() << " tets" << std::endl;
  std::cout << "point not in tet mesh.  This is bad" << std::endl;
  std::cout << point << std::endl;
  
  std::cout << "using best tet, number: " << bestTet 
			<< "\nwith bcs: " << bestBcs << std::endl;

  barycentricCoords(0) = bestBcs(0);
  barycentricCoords(1) = bestBcs(1);
  barycentricCoords(2) = bestBcs(2);
  barycentricCoords(3) = 1.0 - bestBcs.sum();
  return bestTet;

}
