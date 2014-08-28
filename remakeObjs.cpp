#include <eigen/dense>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <iostream>

int main(int argc, char** argv){
  
  const std::string usage("./remakeObjs <format string, eg: frames/foo-%03d.%04d>");
  if(argc != 2){
	std::cout << usage << std::endl;
	exit(1);
  }
  
  const std::string formatString{argv[1]};
  const std::string objFormat = formatString + ".obj";
  const std::string binFormat = formatString + ".bin";
  char filename[1024];

  Eigen::MatrixXd binVertices;
  for(size_t objectIndex = 0; ; ++objectIndex){
	std::cout << "object: " << objectIndex << std::endl;
	Eigen::MatrixXd objVertices;
	Eigen::MatrixXi objFaces;
	
	sprintf(filename, objFormat.c_str(), objectIndex, 0);
	if(!igl::readOBJ(filename, objVertices, objFaces)){
	  std::cout << "no obj for index: " << objectIndex << std::endl;
	  break;
	}

	for(size_t frameIndex = 1; ; ++frameIndex){
	  sprintf(filename, binFormat.c_str(), objectIndex, frameIndex);
	  std::cout << "reading binary file: " << filename << std::endl;
	  std::ifstream ins(filename, std::ios_base::in | std::ios_base::binary);
	  if(!ins){
		std::cout << "no bin for index: " << frameIndex << std::endl;
		break;
	  }
	  size_t numVerts;
	  ins.read(reinterpret_cast<char*>(&numVerts),
			   sizeof(numVerts));
	  assert(numVerts = objVertices.rows());
	  binVertices.resize(numVerts, 3);
	  ins.read(reinterpret_cast<char*>(binVertices.data()),
			   3*numVerts*sizeof(double));
	  sprintf(filename, objFormat.c_str(), objectIndex, frameIndex);
	  std::cout << "writing: " << filename << std::endl;
	  igl::writeOBJ(filename, binVertices, objFaces);
	}
  }

  return 0;
}
