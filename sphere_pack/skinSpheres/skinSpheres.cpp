
#include <fstream>

#include "readOBJ.h"
#include "readMESH.h"
#include "fastGrid.h"
#include "plyIO.hpp"
#include "range.hpp"
#include "slVector.H"


using RMMatrix3f = Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix4f = Eigen::Matrix<float, Eigen::Dynamic, 4, Eigen::RowMajor>;
using RMMatrix3i = Eigen::Matrix<int  , Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix4i = Eigen::Matrix<int  , Eigen::Dynamic, 4, Eigen::RowMajor>;


void readSphereFile(char *fname, std::vector<std::pair<Eigen::Vector3f, double> > &spheres) {
  int nspheres;
  float a,b,c;
  double r;
  std::ifstream in(fname, std::ios::in);

  in>>nspheres;
  for (unsigned int i=0; i<nspheres; i++) {
	in>>a>>b>>c>>r;
	spheres.push_back(std::pair<Eigen::Vector3f, double>(Eigen::Vector3f(a,b,c), r));
  }
}

int main(int argc, char **argv){

  const std::string usage = "./skinSpheres <pathToHighResMesh> <pathToRestSpaceTets> <pathToDeformedMesh.ply (should have the same vertices as the previous .mesh)> <outputFilename.ply>";

  if(argc != 5){
	std::cout << usage << std::endl;
	exit(1);
  }


  const std::string restSpaceTetFilename{argv[2]};
  const std::string deformedTriMeshFilename{argv[3]};


  std::vector<std::pair<Eigen::Vector3f, double> > spheres;
  readSphereFile(argv[1], spheres);

  RMMatrix3f restSpaceTetVertices;
  RMMatrix4i tetTets;
  RMMatrix3i tetTriangleTrash;
  std::cout << "loading rest space tet Mesh: " << restSpaceTetFilename << std::endl;
  igl::readMESH(restSpaceTetFilename, restSpaceTetVertices, tetTets, tetTriangleTrash);
  std::cout<<tetTets.rows()<<" "<<tetTets.cols()<<std::endl;
  std::cout<<restSpaceTetVertices.rows()<<" "<<restSpaceTetVertices.cols()<<std::endl;

  RMMatrix3f deformedTetVertices;
  RMMatrix3i defTriangleTrash;
  std::cout << "loading deformed triMesh: " << deformedTriMeshFilename << std::endl;
  {
	std::ifstream plyIn(deformedTriMeshFilename);
	readPLY(plyIn, deformedTetVertices, defTriangleTrash);
  }
  std::cout<<deformedTetVertices.rows()<<" "<<deformedTetVertices.cols()<<std::endl;

  std::cout << "making grid" << std::endl;
  auto grid = makeFastGrid(restSpaceTetVertices, tetTets, 16);


  RMMatrix4f barycentricCoordinates(spheres.size(), 4);
  RMMatrix4i containingTetVertices(spheres.size(), 4);
  for(auto i : benlib::range(spheres.size())){
	Eigen::Vector3f pos = spheres[i].first;

	auto whichTet =
	  grid.containingTet(pos, barycentricCoordinates.row(i));
	containingTetVertices.row(i) = tetTets.row(whichTet);
  }


  std::ofstream out(argv[4], std::ios::out);
  out<<spheres.size()<<std::endl;
  for(auto i : benlib::range(spheres.size())){
	Eigen::Vector3f pos =
	  barycentricCoordinates(i, 0)*deformedTetVertices.row(
		  containingTetVertices(i, 0)) +
	  barycentricCoordinates(i, 1)*deformedTetVertices.row(
		  containingTetVertices(i, 1)) +
	  barycentricCoordinates(i, 2)*deformedTetVertices.row(
		  containingTetVertices(i, 2)) +
	  barycentricCoordinates(i, 3)*deformedTetVertices.row(
		  containingTetVertices(i, 3));
	out<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" "<<spheres[i].second<<std::endl;
  }

  return 0;
}
