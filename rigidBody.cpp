
#include "rigidBody.h"

#include <fstream>
#include <vector>
#include "cppitertools/range.hpp"
#include <BulletCollision/CollisionShapes/btBoxShape.h>

using iter::range;

//write out an obj for this object
void RigidBody::dump(std::string filename){

  std::ofstream outs(filename.c_str());

  auto transform = bulletBody->getWorldTransform();
  
  if(rbType == RB_BOX){
    std::vector<btVector3> vertices(12);
    
    for(auto i : range(8)){
      dynamic_cast<btBoxShape*>(shape.get())->getVertex(i, vertices[i]);
      vertices[i] = transform(vertices[i]);
      outs << "v " << vertices[i][0] << ' ' << vertices[i][1] << ' ' << vertices[i][2] << std::endl;
    }

    outs << "f 1 2 3" << std::endl
	 << "f 2 4 3" << std::endl
	 << "f 1 6 2" << std::endl
	 << "f 1 5 6" << std::endl
	 << "f 1 3 7" << std::endl
	 << "f 1 7 5" << std::endl
	 << "f 2 6 4" << std::endl
	 << "f 4 6 8" << std::endl
	 << "f 3 4 7" << std::endl
	 << "f 4 8 7" << std::endl
	 << "f 5 7 6" << std::endl
	 << "f 6 7 8" << std::endl;

 
  } //todo, draw other RB shapes


}


