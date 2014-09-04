#pragma once

#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <sstream>



inline void writePlyHeader(std::ofstream& outs, size_t numVertices, size_t numFaces){

  outs << "ply\n"
	"format binary_little_endian 1.0\n"
	"element vertex " << numVertices
	   << "\n"
	"property float32 x\n"
	"property float32 y\n"
	"property float32 z\n"
	"element face " << numFaces
	   << "\n" 
	"property list uchar int32 vertex_indices\n"
	"end_header\n";
  

}

//should work for eigen type stuff
//has to be ROW MAJOR FLOATX3 Matrix
template<typename VType, typename TType>
void writePLY(std::ofstream& outs, const VType& vertices, const TType& triangles){

  assert(vertices.cols() == 3);
  assert(triangles.cols() == 3);
  writePlyHeader(outs, vertices.rows(), triangles.rows());
  
  outs.write(reinterpret_cast<const char*>(vertices.data()),
			 3*vertices.rows()*sizeof(float));
  const unsigned char three{3u};
  for(size_t i = 0; i < triangles.rows(); ++i){
	outs.write(reinterpret_cast<const char*>(&three), sizeof(three));
	outs.write(reinterpret_cast<const char*>(triangles.data() + 3*i),
			   3*sizeof(int));
  }
}

inline void writePLY(std::ofstream& outs, 
			  const std::vector<float>& vertices, 
			  const std::vector<int>& triangles){

  assert((vertices.size() % 3) == 0);
  assert((triangles.size() % 3) == 0);

  writePlyHeader(outs, vertices.size()/3, triangles.size()/3);
  
  outs.write(reinterpret_cast<const char*>(vertices.data()),
			 vertices.size()*sizeof(float));

  const unsigned char three{3u};
  for(size_t i = 0; i < triangles.size()/3; ++i){
	outs.write(reinterpret_cast<const char*>(&three), sizeof(three));
	outs.write(reinterpret_cast<const char*>(triangles.data()+3*i), 
			   3*sizeof(int));
  }

}



/*this will work for files that I wrote, but probably not much else.
Assumes only 3D vertices with floats, and triangles, with int indices
and nothing else
MATRICES NEED TO BE ROW MAJOR!!!
*/
template<typename VType, typename TType>
void readPLY(std::ifstream& ins, VType& vertices, TType& triangles){
  std::string magic;
  ins >> magic;
  assert(magic == "ply");
  std::string firstWord;
  //read until we get to the end of the header
  int numVertices{-1};
  int numFaces{-1};
  for(std::string line; std::getline(ins, line);){
	//std::cout << "read line: " << line << std::endl;
	std::stringstream lineStream(line);
	lineStream >> firstWord;
	if(firstWord == std::string("element")){
	  std::string secondWord;
	  lineStream >> secondWord;
	  //std::cout << "second word: " << secondWord << std::endl;
	  if(secondWord == "vertex"){
		lineStream >> numVertices;
	  } else if(secondWord == "face"){
		lineStream >> numFaces;
	  }
	} else if(firstWord == std::string("end_header")){
	  break;
	}
	//ignore anythign else
	
  }
  if(numVertices <= 0){ 
	vertices.resize(0,3);
	triangles.resize(0,3);
	return;
  }
  assert(numFaces > 0);
  vertices.resize(numVertices, 3);
  triangles.resize(numFaces, 3);
  ins.read(reinterpret_cast<char*>(vertices.data()),
		   numVertices*3*sizeof(float));
  unsigned char three;
  for(int i = 0; i < numFaces; ++i){
	ins.read(reinterpret_cast<char*>(&three), sizeof(three));
	assert(three == 3);
	ins.read(reinterpret_cast<char*>(triangles.data() + 3*i),
			 3*sizeof(int));
  }
  

}
