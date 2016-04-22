
#include "rigidBody.h"

#include <fstream>
#include <vector>

#include <BulletCollision/CollisionShapes/btCollisionShape.h>
#include <BulletDynamics/Dynamics/btRigidBody.h>
#include <BulletCollision/CollisionShapes/btBoxShape.h>
#include <BulletCollision/CollisionShapes/btSphereShape.h>

#include <BulletCollision/Gimpact/btGImpactShape.h>
#include <BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h>
#include <iostream>
#include "plyIO.hpp"

//#include "cppitertools/range.hpp"
//using iter::range;
#include "range.hpp"
using benlib::range;

using RMMatrix3d = Eigen::Matrix<double,Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix3f = Eigen::Matrix<float ,Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix3i = Eigen::Matrix<int,   Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix4i = Eigen::Matrix<int,   Eigen::Dynamic, 4, Eigen::RowMajor>;



//write out an obj for this object
void RigidBody::dump(std::string filename){


  auto transform = bulletBody->getWorldTransform();
  
  switch(rbType){
  case  RB_BOX:{
	  std::ofstream outs(filename.c_str());

	//repeat vertices so texturing works better
	Eigen::Matrix<float, 14, 3, Eigen::RowMajor> vertices;
    

    for(auto i : range(8)){
	  btVector3 untransformedVertex;
      dynamic_cast<btBoxShape*>(shape.get())->getVertex(i, untransformedVertex);
      btVector3 transformedVertex = transform(untransformedVertex);
	  vertices(i,0) = transformedVertex[0];
	  vertices(i,1) = transformedVertex[1];
	  vertices(i,2) = transformedVertex[2];
      //outs << "v " << vertices[i][0] << ' ' << vertices[i][1] << ' ' << vertices[i][2] << std::endl;
    }
	vertices.row(8) = vertices.row(4);
	vertices.row(9) = vertices.row(5);
	vertices.row(10) = vertices.row(6);
	vertices.row(11) = vertices.row(7);
	vertices.row(12) = vertices.row(5);
	vertices.row(13) = vertices.row(4);

	Eigen::Matrix<float, 14, 1> us;
	us << 0.5f, 0.25f, 
	  0.5f, 0.25f, 
	  0.5f, 0.25f, 
	  0.5f, 0.25f, 
	  0.75f, 0.0f, 
	  0.75f, 0.0f, 
	  0.25f, 0.5f;
	Eigen::Matrix<float, 14, 1> vs;
	vs << 0.75f, 0.75f,
	  0.5f, 0.5f,
	  0.0f, 0.0f,
	  0.25f, 0.25f,
	  0.75f, 0.75f,
	  0.5f, 0.5f,
	  1.0f, 1.0f;


	RMMatrix3i triangles(12, 3);
	triangles << 
	  0, 1, 2,
	  1, 2, 3,
	  2, 3, 6,
	  3, 6, 7,
	  4, 5, 7,
	  4, 6, 7,
	  0, 8, 10,
	  0, 2, 10,
	  1, 3, 9,
	  9, 11, 3,
	  0, 1, 12,
	  0, 12, 13;
	/*    outs << 
	    "f 1 2 3" << std::endl
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
	*/
	writePLYWithTexture(outs, vertices, us, vs, triangles);
	break;
  } 
  case RB_PLANE:{
	  std::ofstream outs(filename.c_str());

	auto* planeShape = dynamic_cast<btStaticPlaneShape*>(shape.get());
	//auto printPoint = [&](const btVector3& p){ 
	//outs << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << std::endl;
	//};
	auto vertices = getPlaneVertices(planeShape);
	Eigen::Matrix<float, 4, 1> us;
	Eigen::Matrix<float, 4, 1> vs;
	us(0) = 0;
	vs(0) = 0;
	
	us(1) = 1;
	vs(1) = 0;

	us(2) = 1;
	vs(2) = 1;

	us(3) = 0;
	vs(3) = 1;

	RMMatrix3i triangles(2, 3);
	//fuck it... 2 sided
	triangles << 0, 1, 2, 0, 2, 3;

	writePLYWithTexture(outs, vertices, us, vs, triangles);
	break;
  }
  case RB_TRIMESH:{
	  std::ofstream outs(filename.c_str());

	std::vector<float> outputVertices(meshVertices.size());
	for(auto i : range(meshVertices.size()/3)){
	  btVector3 pos{meshVertices[3*i],
		  meshVertices[3*i +1],
		  meshVertices[3*i +2]};
	  pos = transform(pos);
	  outputVertices[3*i    ] = pos[0];
	  outputVertices[3*i + 1] = pos[1];
	  outputVertices[3*i + 2] = pos[2];
	  /*	  outs << "v " << pos.x() << ' '
		   << pos.y() << ' '
		   << pos.z() << std::endl;
	  */
	}
	/*	for(auto i : range(meshIndices.size()/3)){
	  outs << "f " << (meshIndices[3*i]  + 1) << ' '
		   << (meshIndices[3*i + 1] + 1)<< ' '
		   << (meshIndices[3*i + 2] + 1) << std::endl;
		   }*/

	writePLY(outs, outputVertices, meshIndices);

	break;
  }
  case RB_SPHERE:{
	std::ofstream outs(filename + ".sphere");
	outs << transform.getOrigin().x()
		 << ' ' << transform.getOrigin().y()
		 << ' ' << transform.getOrigin().z() << std::endl
		 << dynamic_cast<btSphereShape*>(shape.get())->getRadius() << std::endl;

	
  }
  }
  
}


void RigidBody::loadTrimesh(std::string filename, double scale){
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
	  meshVertices.push_back(scale*data);
	  std::cout<<meshVertices[meshVertices.size()-1]<<std::endl;
	  linestream >> data;
	  meshVertices.push_back(scale*data);
	  std::cout<<meshVertices[meshVertices.size()-1]<<std::endl;
	  linestream >> data;
	  meshVertices.push_back(scale*data);
	  std::cout<<meshVertices[meshVertices.size()-1]<<std::endl;
	} else if(lineType == "f"){
	  auto numSlashes = std::count(begin(line), end(line), '/');
	  if(numSlashes == 0){
		int data;
		for(auto i : range(3)){
		  (void)i;//unused
		  linestream >> data;
		  meshIndices.push_back(data -1);
		}
	  } else if(numSlashes == 3){
		int data, trash;
		char slash;
		for(auto i : range(3)){
		  (void)i; //unused
		  linestream >> data >> slash >> trash;
		  assert(slash == '/');
		  meshIndices.push_back(data -1);
		}
	  } else if(numSlashes == 6){
		if(line.find("//") != std::string::npos){
		  int data, trash;
		  char slash;
		  for(auto i : range(3)){
			(void) i; //unused
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
			(void)i; //unused
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

	//int meshIndicessizedividedbythree - meshIndices.size()/3;
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


