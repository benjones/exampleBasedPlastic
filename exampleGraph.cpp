#include "exampleGraph.h"
/*
#ifdef __APPLE__

#  include <OpenGL/gl.h>
#  include <GLUT/glut.h>

#elif defined(_WIN32)
//#define FREEGLUT_STATIC
#  include <GL/freeglut.h>
#else
#  define GL_GLEXT_PROTOTYPES
#  include <GL/gl.h>
#  include <GL/glext.h>
#  include <GL/glut.h>
#endif
*/
#include <iostream>
#include <fstream>
#include <numeric>
#include <functional>
//#include "Skinning.h"

/*
static inline void drawCircle(Eigen::Vector2d p, double r){
  
  const size_t numSegments = 10;
  glBegin(GL_POLYGON);
  for(size_t i = 0; i < numSegments; ++i){
	double theta = i*2*M_PI/numSegments;
	glVertex2d(p(0) + r*cos(theta),
			   p(1) + r*sin(theta));
							
  }
  glEnd();
  
}


void ExampleGraph::draw(){

  GLint oldViewport[4];
  glGetIntegerv(GL_VIEWPORT, oldViewport);

  //mask off area to draw on
  glViewport(0, 0, drawingWidth, drawingHeight);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  gluOrtho2D(0, 1, 0, 1);
  glColor3f(0,0,0);
  glBegin(GL_QUADS);
  glVertex2f(0, 0);
  glVertex2f(1, 0);
  glVertex2f(1, 1);
  glVertex2f(0, 1);
  glEnd();
  
  glColor3f(1,1,1);

  for(auto& simp : simplices){

	glBegin(GL_LINE_LOOP);
	for(auto n : simp){
	  glVertex2dv(nodes[n].renderPosition.data());
	}
	glEnd();

  }

  for(size_t i = 0; i <  simplices.size(); ++i){
	auto& simp = simplices[i];
	auto col = 0.5 + 0.5*i/simplices.size();
	glColor3d(col, col, col);
	glBegin(GL_TRIANGLE_FAN);
	for(auto n : simp){
	  glVertex2dv(nodes[n].renderPosition.data());
	}
	glEnd();
  }


  for(auto& node : nodes){
	//	std::cout << node.color << std::endl;
	//	std::cout << node.color.data() << std::endl;
	glColor3dv(node.color.data());
	drawCircle(node.renderPosition, pointRadius);
	
  }



  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  
  glViewport(oldViewport[0],
			 oldViewport[1],
			 oldViewport[2],
			 oldViewport[3]);

};

bool ExampleGraph::inside(int x, int y){

  return x < drawingWidth && y < drawingHeight;
  
}


bool ExampleGraph::down(int x, int y, bool, //right_click, 
						bool shift_down, bool,//control_down,
						bool){ //meta_down
  
  if(!inside(x, y)){
	is_down = false;
	is_selected = false;
	is_hover = false;
	return false;
  } else {

	if(shift_down){
	  selectingSimplexNodes = true;
	  
	  int here = whichNodeIsHere(x,y);
	  if(here >= 0){
		currentlySelectedSimplexNodes.push_back(here);
	  }

	} else {
	  selectingSimplexNodes = false;
	  makeSimplex();
	}


	if(is_down){
	  last_x = x;
	  last_y = y;

	} else {
	  is_down = true;
	  down_x = x;
	  down_y = y;
	  last_x = x;
	  last_y = y;
	  
	  std::vector<double> barycentricCoords;
	  int selectedSimplex = whichSimplexIsHere(x,y, barycentricCoords);
	  if(selectedSimplex >= 0){
		setInterpolatedTransform(selectedSimplex, barycentricCoords);
		
	  }

	  selectedNode = whichNodeIsHere(x,y);

	}
	return true;
  }

  return false;
}
bool ExampleGraph::up(int x, int y, bool right_click, 
					  bool shift_down, bool control_down,
					  bool meta_down){

  is_down = false;
  is_selected = false;
  is_hover = false;
  

  return inside(x,y);

   
}

bool ExampleGraph::drag(int x, int y, bool right_click,
						bool shift_down, bool control_down, 
						bool meta_down){

  if(!is_down || !inside(x,y)){
	return false;
  }
  if(selectedNode >= 0){
	nodes[selectedNode].renderPosition = 
	  Eigen::Vector2d(static_cast<double>(x)/drawingWidth,
					  static_cast<double>(y)/drawingHeight);
  }
  std::vector<double> barycentricCoords;
  int selectedSimplex = whichSimplexIsHere(x,y, barycentricCoords);
  if(selectedSimplex >= 0){
	setInterpolatedTransform(selectedSimplex, barycentricCoords);
  }

  return true;
}

bool ExampleGraph::move(int x, int y, 
						bool shift_down, bool control_down, 
						bool meta_down){
  return false;
}

int ExampleGraph::whichNodeIsHere(int x, int y){
  if(nodes.empty()) return -1;
  Eigen::Vector2d p(static_cast<double>(x)/drawingWidth,
				   static_cast<double>(y)/drawingHeight);

  double rSquared = pointRadius*pointRadius;
  auto minIter = std::min_element(nodes.begin(),
								  nodes.end(),
								  [=](const ExampleNode& a, const ExampleNode& b){
									return (a.renderPosition - p).squaredNorm() <
									(b.renderPosition - p).squaredNorm();
								  });
  double dist = (minIter->renderPosition - p).squaredNorm();
  return (dist < rSquared)  ? std::distance(nodes.begin(), minIter) : -1;

}

void ExampleGraph::setTransformsToNode(int selectedNode){
  auto& node = nodes[selectedNode];
  std::cout << "setting transforms" << std::endl;
  std::vector<Bone*> bones = gather_bones(parent->roots);
  for(size_t i = 0; i < bones.size(); ++i){
	bones[i]->translation = node.transformations[i].translation;
	bones[i]->rotation = node.transformations[i].rotation;
  }
  //  parent->dial_in_each_T;
}
*/
void ExampleGraph::save(std::string filename){
  
  std::ofstream out(filename);
  if(!out){
	std::cout << "couldn't save example graph:  " << filename << std::endl;
	return;
  }
  std::cout << "saving: " << filename << std::endl;
  //format: num nodes, then for each node

  out << nodes.size() << std::endl;
  
  for(auto& n : nodes){
	//render pos:
	out << n.renderPosition(0) << ' ' << n.renderPosition(1) << ' '
		<< n.color(0) << ' ' 
		<< n.color(1) << ' ' 
		<< n.color(2) << std::endl;
	//# of transformations
	out << n.transformations.size() << std::endl;
	for(auto& t : n.transformations){
	  //transformations in order 
	  out << t.translation(0) << ' ' 
		  << t.translation(1) << ' ' 
		  << t.translation(2) << ' ' 
		  << t.rotation.w() << ' ' 
		  << t.rotation.x() << ' ' 
		  << t.rotation.y() << ' ' 
		  << t.rotation.z() << std::endl;
	}
  }
  
  out << simplices.size() << std::endl;
  for(auto& simp : simplices){
	out << simp.size() << std::endl;
	for(auto n : simp){
	  out << n << ' ';
	}
	out << std::endl;
  }

}


void ExampleGraph::load(std::string filename){
  std::ifstream ins(filename);
  if(!ins){
	std::cout << "couldn't open example graph file: " << filename << std::endl;
	return;
  }
  std::cout << "loading: " << filename << std::endl;
  size_t numNodes;
  ins >> numNodes;

  nodes.resize(numNodes);
  
  for(auto& n : nodes){
	ins >> n.renderPosition(0) >> n.renderPosition(1)
		>> n.color(0) >> n.color(1) >> n.color(2);
	size_t numTransformations;
	ins >> numTransformations;
	n.transformations.resize(numTransformations);
	
	for(auto& t : n.transformations){
	  ins >> t.translation(0) >> t.translation(1) >> t.translation(2)
		  >> t.rotation.w() >> t.rotation.x() >> t.rotation.y() >> t.rotation.z();
	}
  }
  std::cout << "loaded: " << nodes.size() << " nodes to EG" << std::endl;
  
  size_t numSimplices;
  ins >> numSimplices;
  simplices.resize(numSimplices);
  for(auto& simp : simplices){
	size_t simplexSize;
	ins >> simplexSize;
	simp.resize(simplexSize);
	for(auto& n : simp){
	  ins >> n;
	}
  }
  std::cout << "read : " << simplices.size() << " simplices" << std::endl;
}

/*
void ExampleGraph::makeSimplex(){

  if(currentlySelectedSimplexNodes.empty()){
	return;
  }

  std::sort(currentlySelectedSimplexNodes.begin(),
			currentlySelectedSimplexNodes.end());
  
  //remove duplicates
  currentlySelectedSimplexNodes.resize(std::distance(currentlySelectedSimplexNodes.begin(),
													 std::unique(currentlySelectedSimplexNodes.begin(),
																 currentlySelectedSimplexNodes.end())));

  
  std::cout << "making simplex with: " << currentlySelectedSimplexNodes.size() << " nodes" << std::endl;
  

  //check if this is a party of any other simplex
  for(auto& simp : simplices){
	if(std::includes(simp.begin(), simp.end(),
					 currentlySelectedSimplexNodes.begin(),
					 currentlySelectedSimplexNodes.end())){
	  std::cout << "simplex is included in an existing simplex.  Ignoring" << std::endl;
	  currentlySelectedSimplexNodes.clear();
	  return;
	}
  }

  //remove all simplices that are fully contained in this one
  //stdlib algorithms FTW
  simplices.erase(std::remove_if(simplices.begin(),
								 simplices.end(),
								 [this](const std::vector<size_t>& simp){
								   return std::includes(currentlySelectedSimplexNodes.begin(),
														currentlySelectedSimplexNodes.end(),
														simp.begin(), simp.end());
								 }),
				  simplices.end());
  



  simplices.push_back(currentlySelectedSimplexNodes);
  currentlySelectedSimplexNodes.clear();




}


int ExampleGraph::whichSimplexIsHere(int x, int y, std::vector<double>& barycentricCoords){

  const double threshold = 0.05; //you can be a couple pixels away from a simplex I guess.

  int best = -1;
  double closestDistance = std::numeric_limits<double>::max();
  std::vector<double> theseBarycentricCoords;
  for(size_t i = 0; i < simplices.size(); ++i){
	auto& simp = simplices[i];
	double thisDistance;
	switch (simp.size()){
	case 2:
	  thisDistance = projectPointToLine(x, y, simp, theseBarycentricCoords);
	  if( thisDistance < closestDistance){
		closestDistance = thisDistance;
		best = i;
		barycentricCoords = theseBarycentricCoords;
	  }
	  break;
	case 3:
	  thisDistance = projectPointToTriangle(x,y, simp, theseBarycentricCoords);
	  if(thisDistance < closestDistance){
		closestDistance = thisDistance;
		best = i;
		barycentricCoords = theseBarycentricCoords;
	  }
	  break;
	case 1:
	  break; //ignore 1-node simplices
	default:
	  std::cout << "unsupported simplex: size " << simp.size() << std::endl;
	}
  }

  return closestDistance < threshold ? best : -1;

}

double ExampleGraph::projectPointToLine(int x, int y, 
										const std::vector<size_t>& simp,
										std::vector<double>& barycentricCoords){
  
  Eigen::Vector2d p(static_cast<double>(x)/drawingWidth,
					static_cast<double>(y)/drawingHeight);
  Eigen::Vector2d v = 
	nodes[simp[1]].renderPosition - 
	nodes[simp[0]].renderPosition;
  
  Eigen::Vector2d s = p - nodes[simp[0]].renderPosition;

  barycentricCoords.resize(2);
  //clamp it to the range [0,1]
  barycentricCoords[1] = std::max(std::min(s.dot(v)/v.squaredNorm(), 1.0), 0.0);
  barycentricCoords[0] = 1 - barycentricCoords[1];

  auto projectedPoint = 
	barycentricCoords[0]*nodes[simp[0]].renderPosition +
	barycentricCoords[1]*nodes[simp[1]].renderPosition;
  //this should handle nans gracefully
  return (p - projectedPoint).norm();
}


double ExampleGraph::projectPointToTriangle(int x, int y, const std::vector<size_t>& simp,
							std::vector<double>& barycentricCoords){
  Eigen::Vector2d p(static_cast<double>(x)/drawingWidth,
					static_cast<double>(y)/drawingHeight);

  Eigen::Vector2d p0 = nodes[simp[0]].renderPosition,
	p1 = nodes[simp[1]].renderPosition,
	p2 = nodes[simp[2]].renderPosition;
  
  
  Eigen::Matrix2d basis;
  basis.col(0) = p0 - p2;
  basis.col(1) = p1 - p2;

  Eigen::Vector2d lambda01 = basis.inverse()*(p - p2);
  barycentricCoords.resize(3);
  barycentricCoords[0] = lambda01(0);
  barycentricCoords[1] = lambda01(1);
  barycentricCoords[2] = 1.0 - barycentricCoords[0] - barycentricCoords[1];

  auto projectedPoint = barycentricCoords[0]*p0 +
	barycentricCoords[1]*p1 +
	barycentricCoords[2]*p2;

  //if the point's not in the triangle
  if(std::any_of(barycentricCoords.begin(), barycentricCoords.end(),
				 [](double b){ return b < 0 || b > 1;})){
	return std::numeric_limits<double>::infinity();
	//punt on finding how far we are for now.

  }

  return (p - projectedPoint).norm();
}

void ExampleGraph::setInterpolatedTransform(int selectedSimplex,
											const std::vector<double>& barycentricCoords){
  
  std::vector<Bone*> bones = gather_bones(parent->roots);
  std::vector<double> weights(nodes.size(), 0.0);
  for(size_t i = 0; i < simplices[selectedSimplex].size(); ++i){
	weights[simplices[selectedSimplex][i]] = barycentricCoords[i];
  }
  for(size_t i = 0; i < bones.size(); ++i){
	auto trans = interpolateTransform(i, weights);
	bones[i]->translation = trans.translation;
	bones[i]->rotation = trans.rotation;
	
  }

  

}

ExampleNode::TransformInfo ExampleGraph::interpolateTransform(size_t boneIndex, 
															  const std::vector<double>& weights){

  std::vector<Quat> quats(weights.size());
  std::transform(nodes.begin(),nodes.end(),
				 quats.begin(),
				 [this,boneIndex](const ExampleNode& node){
				   return node.transformations[boneIndex].rotation;});
  
  ExampleNode::TransformInfo ret;
  
  
  ret.translation = std::inner_product(nodes.begin(), nodes.end(),
									   weights.begin(), Vec3{0.0,0.0,0.0},
									   std::plus<Vec3>{}, //add the weighted vectors
									   [this,boneIndex](const ExampleNode& node, double w){
										 return w*node.transformations[boneIndex].translation;}
									   );
  ret.rotation = qlerp(quats, weights);
  
  return ret;
}
*/
