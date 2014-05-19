#include "egTraverser.h"

#include <numeric> //I <3 accumulate + inner_product

#include "utils.h"
#include "cppitertools/range.hpp"
#include "cppitertools/enumerate.hpp"
using iter::range;
using iter::enumerate;



void EGTraverser::traverse(){

  
  currentPosition = interpolatePoints(currentPosition,
									  restPosition,
									  restSpringStrength);
  
  
}


EGPosition EGTraverser::interpolatePoints(const EGPosition& first,
										  const EGPosition& second,
										  double alpha){



  if(first.simplex == second.simplex){
	assert(first.coords.size() == second.coords.size());

	EGPosition ret;
	ret.simplex = first.simplex;
	ret.coords.resize(first.coords.size());
	//interpolate coords in this simplex
	std::transform(first.coords.begin(),
				   first.coords.end(),
				   second.coords.begin(),
				   ret.coords.begin(),
				   [alpha](double f, double s){
					 return (1.0 - alpha)*f + alpha*s;
				   });
				   
	
	return ret;
	
  } else {
	double totalDistance = distanceBetweenSimplexPoints(first, second);
	double distanceToTravel = alpha*totalDistance;
	
	int nextSimplex = shortestPaths(first.simplex, second.simplex);
	if(nextSimplex < 0){
	  //just set it since it's not possible to pull there
	  return alpha == 0 ? first : second;
	}

	EGPosition boundary1, boundary2;
	getCorrespondingPoint(first, nextSimplex,
						  boundary1, boundary2);
	
	double distanceToBoundary = distanceBetweenSimplexPoints(first, boundary1);
	
	if(distanceToTravel <= distanceToBoundary){
	  //we can stay in this simplex
	  return interpolatePoints(first, boundary1,
							   distanceToTravel/distanceToBoundary);
	}
	
	return interpolatePoints(boundary2, second,
							 (distanceToTravel - distanceToBoundary)/
							 (totalDistance - distanceToBoundary));
	

  }
  
  
}

void EGTraverser::computeShortestPaths(){
  assert(graph);

  computeEdges();

  shortestPaths.resize(graph->simplices.size(),
					   graph->simplices.size());

  shortestPaths.setConstant(-1);
  
  Eigen::MatrixXi distances(shortestPaths.rows(),
							shortestPaths.cols());
  distances.setConstant(9999);  

  for(auto i : range(shortestPaths.rows())){
	shortestPaths(i,i) = i;
	distances(i,i) = 0;
  }


  for(auto& edge : edges){
	distances(edge.first,edge.second) = 1;
	distances(edge.second,edge.first) = 1;
	shortestPaths(edge.first, edge.second) = edge.second;
	shortestPaths(edge.second, edge.first) = edge.first;
  }
  
  for(auto iter : range(graph->simplices.size())){
	for(auto& edge: edges){
	  //update all paths from i that are shorter via j
	  for(auto k : range(graph->simplices.size())){
		if((distances(edge.second, k) +1) < distances(edge.first, k)){
		  distances(edge.first, k) = distances(edge.second, k) + 1;
		  shortestPaths(edge.first, k) = edge.second;
		}
		if((distances(edge.first, k) +1) < distances(edge.second, k)){
		  distances(edge.second, k) = distances(edge.first, k) + 1;
		  shortestPaths(edge.second, k) = edge.first;
		}
	  }
	}
  }
}


void EGTraverser::computeEdges(){
  assert(graph);
  edges.clear();
  std::vector<int> samesies;
  for(auto i : range(graph->simplices.size())){
	for(auto j : range(i, graph->simplices.size())){
	  samesies.clear();
	  std::set_union(graph->simplices[i].begin(),
					 graph->simplices[i].end(),
					 graph->simplices[j].begin(),
					 graph->simplices[j].end(),
					 std::back_inserter(samesies));

	  //todo, do more with this ?
	  if(!samesies.empty()){
		//i and j are already sorted... but oh well
		edges.push_back(makeSortedPair(i,j));
	  }
	}
  }
  //the edges are all sorted
}


double EGTraverser::distanceBetweenSimplexPoints(const EGPosition& first,
												 const EGPosition& second){

  if(first.simplex == second.simplex){
	return sqrt(std::inner_product(first.coords.begin(), first.coords.end(),
								   second.coords.begin(), 0,
								   std::plus<double>{},
								   [](double a, double b){return (a-b)*(a-b);}));
  } else {
	//find next simplex
	int nextSimplex = shortestPaths(first.simplex, second.simplex);
	if(nextSimplex < 0){ return std::numeric_limits<double>::infinity();}
	
	//find the closest point on the neighboring simplex
	EGPosition nearestFirst, nearestSecond;
	getCorrespondingPoint(first, nextSimplex, nearestFirst, nearestSecond);
	
	return distanceBetweenSimplexPoints(first, nearestFirst) +
	  distanceBetweenSimplexPoints(nearestSecond, second);

  }

}

void EGTraverser::getCorrespondingPoint(const EGPosition& first,
										int nextSimplex,
										EGPosition& out1,
										EGPosition& out2){


  std::vector<std::pair<int, int>> weightCorrespondences;
  for(auto e1 : enumerate(graph->simplices[first.simplex])){
	for(auto e2 : enumerate(graph->simplices[nextSimplex])){
	  if(e1.element == e2.element){
		weightCorrespondences.emplace_back(e1.index, e2.index);
	  }
	}
  }

  auto sum = std::accumulate(weightCorrespondences.begin(),
							 weightCorrespondences.end(),
							 0,
							 [&first](double s, std::pair<int, int> p){
							   return s + first.coords[p.first];
							 });
  
  out1.simplex = first.simplex;
  out2.simplex = nextSimplex;
  out1.coords.assign(first.coords.size(), 0.0);
  out2.coords.assign(graph->simplices[nextSimplex].size(), 0.0);

  if(sum > 0){
	for(auto& wi : weightCorrespondences){
	  out1.coords[wi.first] = first.coords[wi.first]/sum;
	  out2.coords[wi.second] = out1.coords[wi.first];
	}
  } else {
	double factor = 1.0/weightCorrespondences.size();
	for(auto& wi : weightCorrespondences){
	  out1.coords[wi.first] = factor;
	  out2.coords[wi.second] = factor;
	}
  }
  
  
}


void EGTraverser::zeroWeightsSum(std::vector<double>& weights){

  auto sum = std::accumulate(weights.begin(), weights.end(), 0);
  for(auto& w : weights){
	w -= sum/weights.size();
  }
}
void EGTraverser::zeroWeightsSum(Eigen::VectorXd& weights){

  weights.array() -= weights.sum()/weights.rows();
  if(fabs(weights.sum()) > 1e-8){
	std::cout << "zeroWeightsSum failed" << std::endl;
	exit(1);
  }
}



EGPosition EGTraverser::addBcVector(const EGPosition& start, const std::vector<double>& direction){
  
  assert(direction.size() == start.coords.size());
  
  EGPosition ret;
  ret.simplex = start.simplex;
  ret.coords = start.coords;
  
  //find scale
  double s = 1;
  for(auto i : range(direction.size())){
	if(direction[i] < 0){ 
	  s = std::min(-ret.coords[i]/direction[i], s);
	}
	if(direction[i] > 0){
	  s = std::min((1 - ret.coords[i])/direction[i], s);
	}
  }
  
  std::transform(ret.coords.begin(), ret.coords.end(), direction.begin(),
				 ret.coords.begin(),
				 [s](double a, double b){
				   return a + s*b;
				 });
  if(std::any_of(ret.coords.begin(), ret.coords.end(),
			  [](double a){return a < 0 || a > 1;})){
	
	std::cout << "addBcVector screwed up" << std::endl;
	std::cout << start.coords << std::endl;
	std::cout << "dir: " << direction << std::endl;
	std::cout << "ret: " << ret.coords << std::endl;
	exit(1);
	
  }
  
  return ret;

}
EGPosition EGTraverser::addBcVector(const EGPosition& start, const Eigen::VectorXd& direction){
  return addBcVector(start, eigenToStd(direction));
  
}
