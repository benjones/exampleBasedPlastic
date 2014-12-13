#pragma once

#include <vector>
#include "Eigen/Dense"
//#include "Tform.h"
//#include "Widget.h"
#include "Quat.h"

//#include "Skeleton.cpp"
//class Skinning;

struct ExampleNode{
  
  Eigen::Vector2d renderPosition;
  Eigen::Vector3d color;
  
  //std::vector<Tform3> transforms;
  struct TransformInfo{
	Eigen::Vector3d translation;
	Quat rotation;
  };
  std::vector<TransformInfo> transformations;
  
};

//ripped out of teh fast skinning stuff.  Just need the save/load functionality

class ExampleGraph {
public:
  std::vector<ExampleNode> nodes;
  std::vector<std::vector<size_t>> simplices;

  //when he handles your mouse click...
  //just widget things <3 <3 <3
  /*
  void draw();
  bool down(int x, int y, bool right_click, 
			bool shift_down, bool control_down,
			bool meta_down);
  bool up(int x, int y, bool right_click, 
		  bool shift_down, bool control_down,
		  bool meta_down);

  bool drag(int x, int y, bool right_click,
			bool shift_down, bool control_down, 
			bool meta_down);
			
  bool move(int x, int y, 
			bool shift_down, bool control_down, 
			bool meta_down);

  bool inside(int x, int y);

  int whichNodeIsHere(int x, int y);

  int whichSimplexIsHere(int x, int y, std::vector<double>& barycentricCoords);

  void makeSimplex();

  void setTransformsToNode(int selectedNode);

  void setInterpolatedTransform(int selectedSimplex, 
								const std::vector<double>& barycentricCoords);
  */
  void save(const std::string& filename);
  void load(const std::string& filename);
  /*
  int selectedNode; //which node are we moving
  

  const int drawingWidth = 400;
  const int drawingHeight = 400;
  const double pointRadius = 0.02;
  */
  //Skeleton<Bone>* parent = nullptr;
  /*
  bool selectingSimplexNodes;
  std::vector<size_t> currentlySelectedSimplexNodes;

  

  //return distance to simplex
  double projectPointToLine(int x, int y, const std::vector<size_t>& simp, 
						  std::vector<double>& barycentricCoords);
  
  double projectPointToTriangle(int x, int y, const std::vector<size_t>& simp,
							  std::vector<double>& barycentricCoords);

  
  ExampleNode::TransformInfo 
  interpolateTransform(size_t boneIndex, const std::vector<double>& weights);
  */
};
