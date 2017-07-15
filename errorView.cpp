#include "plyIO.hpp"
#include "utils.h"

using RMMatrix3d = Eigen::Matrix<double,Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix3f = Eigen::Matrix<float ,Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix3i = Eigen::Matrix<int,   Eigen::Dynamic, 3, Eigen::RowMajor>;
using Vector3 = Eigen::Vector3f;

struct Sphere{
  Vector3 center;
  float radius;
};

float sphereDist(const Vector3& v, const Sphere& s){
  return (s.center - v).norm() - s.radius;
}

int main(int argc, char** argv){

  const auto usage = "errorView <pathToFrames> <outputFolder> <camera xml file>";
  if(argc != 4){
	std::cout << usage << std::endl;
	return 1;
  }

  std::string framePath(argv[1]);
  std::string outputPath(argv[2]);
  auto cameraInfo = readFile(argv[3]);

  
  int numObjects = 1000;
  char plyName[1024];
  std::string plyFormat = framePath + "/foo-%03d.%04d.ply";

  char sphereName[1024];
  std::string sphereFormat = framePath + "/foo-%03d.%04d.sph";

  char outName[1024];
  std::string outFormat = outputPath + "/foo-%03d.%04d.ply";

  char mitsubaName[1024];
  std::string mitsubaOutFormat = outputPath + "/mitsubaFrame%04d.xml";

  char meshNameOnly[1024];
  std::string meshOnlyFormat = "foo-%03d.%04d.ply";

  std::string boringShading =
	"<boolean name=\"faceNormals\" value=\"true\" />\n"
	"<bsdf type=\"twosided\">"
	"<bsdf type=\"roughdiffuse\">"
	"<float name=\"alpha\" value=\"0.7\" />"
	"<rgb name=\"reflectance\" value=\"#aaaaaa\" />"
	"</bsdf></bsdf></shape>";
  
  for(int frame = 0; ; ++frame){
	bool found = false;
	std::cout << "processing frame " << frame << std::endl;

	sprintf(mitsubaName, mitsubaOutFormat.c_str(), frame);
	std::ofstream mitsubaOut(mitsubaName);
	mitsubaOut <<   "<?xml version=\"1.0\" encoding =\"utf-8\" ?>\n"
	  "<scene version=\"0.5.0\">\n";
	mitsubaOut << cameraInfo;

	
	for(int i = 0; i < numObjects ; i++){
	  sprintf(plyName, plyFormat.c_str(), i, frame);
	  std::cout << "trying to read " << plyName << std::endl;
	  std::ifstream plyIn(plyName);
	  if(!plyIn.good()){

		if(frame == 0){
		  numObjects = i;
		  break;
		} else {
		  //copy from first frame

		  sprintf(meshNameOnly, meshOnlyFormat.c_str(), i, 0);
		  std::cout << "using frame 0 version... " << meshNameOnly << std::endl;
		  mitsubaOut << "<shape type=\"ply\" ><string name=\"filename\" value=\""
					 << meshNameOnly << "\" /> \n";
		  mitsubaOut << boringShading << std::endl;
		  continue;
		}
	  } else {
		found = true;
	  }
	  RMMatrix3f plyVertices;
	  RMMatrix3i plyTriangles;
	  readPLY(plyIn, plyVertices, plyTriangles);

	  sprintf(sphereName, sphereFormat.c_str(), i, frame);
	  std::ifstream sphereIn(sphereName);
	  std::cout << "reading spheres: " << sphereName << std::endl;
	  
	  sprintf(outName, outFormat.c_str(), i, frame);
	  std::ofstream plyOut(outName);

	  if(sphereIn.good()){
		std::vector<Sphere> spheres;
		int numSpheres;
		sphereIn >> numSpheres;
		while(true){
		  Sphere s;
		  sphereIn >> s.center.x() >> s.center.y() >> s.center.z() >> s.radius;
		  if(sphereIn.good()){
			spheres.push_back(s);
		  } else {
			break;
		  }
		}
		std::cout << "has " << spheres.size() << " spheres" << std::endl;

		RMMatrix3f colors = RMMatrix3f::Constant(plyVertices.rows(), 3, 1.0f);
		for(auto row = 0; row < plyVertices.rows(); ++row){
		  Vector3 v = plyVertices.row(row);
		  auto it = std::min_element(spheres.begin(), spheres.end(),
			  [&v](const Sphere& a, const Sphere& b){
				return std::abs(sphereDist(v, a))
				< std::abs(sphereDist(v, b));
			  });
		  
		  auto dist = sphereDist(v, *it);
		  if(row == 0){ 
			std::cout << "dist " << dist << std::endl;
		  }
		  const float scale = 0.5; //max dist
		  if(dist < 0){ //make this blue
			auto diff = std::min(1.0f, -dist/scale);
			colors(row, 0) -= diff;
			colors(row, 1) -= diff;
		  } else {
			auto diff = std::min(1.0f, dist/scale);
			colors(row, 1) -= diff;
			colors(row, 2) -= diff;
			
		  }
		  
		}
		std::cout << "writing " << outName  << " with colors" << std::endl;
		sprintf(meshNameOnly, meshOnlyFormat.c_str(), i, frame);
		writePLYWithColor(plyOut, plyVertices, colors, plyTriangles);
		mitsubaOut << "<shape type=\"ply\"><string name=\"filename\" value=\""
				   << meshNameOnly << "\" />\n"
				   << "<boolean name=\"faceNormals\" value=\"true\" />\n"
				   << "<bsdf type=\"twosided\">"
		  "<bsdf type=\"diffuse\" >"
		  "<texture name=\"reflectance\" type=\"vertexColors\" />"
		  "</bsdf></bsdf></shape>" << std::endl;
	  } else {
		std::cout << "writing " << outName  << " without colors" << std::endl;
		sprintf(meshNameOnly, meshOnlyFormat.c_str(), i, frame);
		writePLY(plyOut, plyVertices, plyTriangles);
		mitsubaOut << "<shape type=\"ply\"><string name=\"filename\" value=\""
				   << meshNameOnly << "\" />\n"
				   << boringShading  << std::endl;

		
	  }

	  
	  
	  
	  
	}
	mitsubaOut << "</scene>" << std::endl;
	if(!found){
	  break;
	}
  }
  
}
