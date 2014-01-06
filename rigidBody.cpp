
#include "rigidBody.h"

#include <fstream>
#include <vector>
#include "cppitertools/range.hpp"
#include <BulletCollision/CollisionShapes/btBoxShape.h>
#include <BulletCollision/CollisionShapes/btStaticPlaneShape.h>

using iter::range;

//write out an obj for this object
void RigidBody::dump(std::string filename){

  std::ofstream outs(filename.c_str());

  auto transform = bulletBody->getWorldTransform();
  
  switch(rbType){
  case  RB_BOX:{
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

	break;
  } 
  case RB_PLANE:{
	auto* planeShape = dynamic_cast<btStaticPlaneShape*>(shape.get());
	auto normal = planeShape->getPlaneNormal();
	auto localPoint = -planeShape->getPlaneConstant()*normal;
	
	auto span1 = normal.cross(btVector3{1, 0, 0});
	if(span1.length2() < 1e-3){ //normal is nearly parallel to 1, 0, 0, so don't use it
	  span1 = normal.cross(btVector3{0, 0, 1});
	}
	auto span2 = span1.cross(normal);
	auto printPoint = [&](const btVector3& p){ outs << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << std::endl;};
	
	auto tmp = localPoint - 100*span1 - 100*span2;
	printPoint(tmp);
	tmp = localPoint + 100*span1 - 100*span2;
	printPoint(tmp);
	tmp = localPoint + 100*span1 + 100*span2;
	printPoint(tmp);
	tmp = localPoint - 100*span1 + 100*span2;
	printPoint(tmp);

	//check winding, which is probably unnecessary since we know the cross product ordering...
	if(normal.dot(span1.cross(span2)) > 0){
	  outs << "f 1 2 3\nf 1 3 4" << std::endl;
	} else {
	  outs << "f 1 3 2\nf 1 4 3" << std::endl;
	}

	break;
  }
  }

}


