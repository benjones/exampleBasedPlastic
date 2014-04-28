
#include "rigidBody.h"

#include <fstream>
#include <vector>
#include "cppitertools/range.hpp"
#include <BulletCollision/CollisionShapes/btCollisionShape.h>
#include <BulletDynamics/Dynamics/btRigidBody.h>
#include <BulletCollision/CollisionShapes/btBoxShape.h>
#include <BulletCollision/CollisionShapes/btStaticPlaneShape.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>
#include <BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h>
#include <iostream>

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
  case RB_TRIMESH:{
	
	for(auto i : range(meshVertices.size()/3)){
	  btVector3 pos{meshVertices[3*i],
		  meshVertices[3*i +1],
		  meshVertices[3*i +2]};
	  pos = transform(pos);
	  outs << "v " << pos.x() << ' '
		   << pos.y() << ' '
		   << pos.z() << std::endl;
	}
	for(auto i : range(meshIndices.size()/3)){
	  outs << "f " << (meshIndices[3*i]  + 1) << ' '
		   << (meshIndices[3*i + 1] + 1)<< ' '
		   << (meshIndices[3*i + 2] + 1) << std::endl;
	}
	break;
  }
  }

}


void RigidBody::loadTrimesh(std::string filename){
  std::cout << "reading: " << filename << std::endl;
  std::ifstream ins(filename);
  if(!ins){
	std::cout << "couldn't read: " << filename << std::endl;
	exit(1);
  }

	//just in case?
  meshVertices.clear();
  meshIndices.clear();  
  
  for(std::string line; std::getline(ins, line);){

	std::stringstream linestream(line);
	std::string lineType;
	linestream >> lineType;

	if(lineType == "v"){
	  double data;
	  linestream >> data;
	  meshVertices.push_back(data);
	  linestream >> data;
	  meshVertices.push_back(data);
	  linestream >> data;
	  meshVertices.push_back(data);

	} else if(lineType == "f"){
	  auto numSlashes = std::count(begin(line), end(line), '/');
	  if(numSlashes == 0){
		int data;
		for(auto i : range(3)){
		  linestream >> data;
		  meshIndices.push_back(data -1);
		}
	  } else if(numSlashes == 3){
		int data, trash;
		char slash;
		for(auto i : range(3)){
		  linestream >> data >> slash >> trash;
		  assert(slash == '/');
		  meshIndices.push_back(data -1);
		}
	  } else if(numSlashes == 6){
		if(line.find("//") != std::string::npos){
		  int data, trash;
		  char slash;
		  for(auto i : range(3)){
			linestream >> data >> slash;
			assert(slash == '/');
			linestream >> slash >> trash;
			assert(slash == '/');
			meshIndices.push_back(data -1);
		  }
		} else {
		  int data, trash1, trash2;
		  char slash1, slash2;
		  for(auto i : range(3)){
			linestream >> data >> slash1 >> trash1 >> slash2 >> trash2;
			assert(slash1 == '/');
			assert(slash2 == '/');
			meshIndices.push_back(data - 1);
		  }
		}
	  } else {
		std::cout << "ignoring malformed face line: " << line << std::endl;
	  }
	} else if(lineType == "vt"){
	  std::cout << "ignoring texture coord" << std::endl;
	} else if(lineType == "vn"){
	  std::cout << "ignoring vertex normal" << std::endl;
	} else {
	  std::cout << "unknown line, ignoring it: " << line << std::endl;
	}
  }

  triMesh = 
	std::unique_ptr<btTriangleIndexVertexArray>{
	new btTriangleIndexVertexArray{
	  static_cast<int>(meshIndices.size()/3),
	  meshIndices.data(),
	  3*sizeof(int),
	  static_cast<int>(meshVertices.size()/3),
	  meshVertices.data(),
	  3*sizeof(typename decltype(meshVertices)::value_type)}};
  
  shape = std::unique_ptr<btCollisionShape>{
	new btGImpactMeshShape{triMesh.get()}};
  dynamic_cast<btGImpactMeshShape*>(shape.get())->updateBound();
  std::cout << meshVertices.size()/3 << " vertices\n" 
			<< meshIndices.size()/3 << "triangles" << std::endl;

}


