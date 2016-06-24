#include "slVector.H"
#include "slUtil.H"
#include "plyIO.hpp"

using RMMatrix3f = Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix4f = Eigen::Matrix<float, Eigen::Dynamic, 4, Eigen::RowMajor>;
using RMMatrix3i = Eigen::Matrix<int  , Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix4i = Eigen::Matrix<int  , Eigen::Dynamic, 4, Eigen::RowMajor>;

void readSphereFile(char *fname, std::vector<std::pair<SlVector3, double> > &spheres) {
  int nspheres;
  SlVector3 c;
  double r;
  std::ifstream in(fname, std::ios::in);

  in>>nspheres;
  for (unsigned int i=0; i<nspheres; i++) {
	in>>c[0]>>c[1]>>c[2]>>r;
	spheres.push_back(std::pair<SlVector3, double>(c, r));
  }
}

int main(int argc, char *argv[]) {
  std::vector<std::pair<SlVector3, double> > spheres;
  std::vector<double> distances;
  double mind;

  if (argc != 4) {
	std::cout<<"usage: acg sphere_file mesh_file output"<<std::endl;
	exit(0);
  }
  readSphereFile(argv[1], spheres);

  RMMatrix3f vertices;
  RMMatrix3i trash;
  std::cout << "loading deformed triMesh: " << argv[2] << std::endl;
  {
	std::ifstream plyIn(argv[2]);
	readPLY(plyIn, vertices, trash);
  }
  std::cout<<vertices.rows()<<" "<<vertices.cols()<<std::endl;

  //readObjFile(argv[2], vertices, triangles);
  
  for (unsigned int i=0; i<vertices.rows(); i++) {
	double mind = std::numeric_limits<double>::max();
	for (unsigned int j=0; j<spheres.size(); j++) {
	  double d = mag (SlVector3(vertices(i,0),vertices(i,1),vertices(i,2)) - spheres[j].first) - spheres[j].second;
	  if (d < mind) mind = d;
	}
	distances.push_back(mind);
  }
  std::ofstream out(argv[3], std::ios::out);
  //out<<distances.size()<<std::endl;
  for (unsigned int i=0; i<distances.size(); i++) {
	out<<distances[i]<<std::endl;
  }
}
