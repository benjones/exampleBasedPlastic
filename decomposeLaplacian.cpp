
#include <unsupported/Eigen/SparseExtra>
#include <igl/readMESH.h>
#include <igl/writeDMAT.h>
#include <igl/cotmatrix.h>
#include <iostream>

using RMMatrix3d = Eigen::Matrix<double,Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix3f = Eigen::Matrix<float ,Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix3i = Eigen::Matrix<int,   Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix4i = Eigen::Matrix<int,   Eigen::Dynamic, 4, Eigen::RowMajor>;


int main(int argc, char** argv){

  std::string usage{"decomposeLaplacian <directory>"};
  
  if(argc != 2){
	std::cout << usage << std::endl;
	return 1;
  }
  std::string directory{argv[1]};
  
  RMMatrix3d vertices;
  RMMatrix4i tets;
  RMMatrix3i triangles;

  if(!igl::readMESH(directory + "/mesh.mesh",
					vertices,
					tets,
					triangles)){
	std::cout << "couldn't open the file: " << directory + "/mesh.mesh" << std::endl;
	return 1;
  }


  Eigen::SparseMatrix<double> cotangentLaplacian;
  igl::cotmatrix(vertices, tets, cotangentLaplacian);
  
  Eigen::MatrixXd denseLaplacian{cotangentLaplacian};
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver{denseLaplacian};
  Eigen::MatrixXd eigenvectors = solver.eigenvectors();

  std::cout << "computed.  Writing evecs to " << directory << "/eigenvectors.dmat" << std::endl;
  igl::writeDMAT(directory + "/eigenvectors.dmat", eigenvectors);

  return 0;
}
