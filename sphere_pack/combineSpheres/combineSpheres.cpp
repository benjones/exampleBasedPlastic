#include "slVector.H"
#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>

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
  std::vector<std::vector<std::pair<SlVector3, double> > > spheres;

  if (argc != 4) {
	std::cout<<"usage: in_sphere_file_format nexamples output_filename"<<std::endl;
	exit(0);
  }

  int nexamples = atoi(argv[2]);
  spheres.resize(nexamples);

  for (int i=0; i<nexamples; i++) {
	char fname[80];
	sprintf (fname, argv[1], i);
	readSphereFile(fname, spheres[i]);
  }

  std::ofstream out(argv[3], std::ios::out);
  std::cout<<"writing "<<argv[3]<<std::endl;
  out<<nexamples<<std::endl<<spheres[0].size()<<std::endl;
  for (unsigned int i=0; i<nexamples; i++) {
	for (unsigned int j=0; j<spheres[i].size(); j++) {
	  std::pair<SlVector3, double> &s = spheres[i][j];
	  out<<s.first[0]<<" "<<s.first[1]<<" "<<s.first[2]<<" "<<s.second<<std::endl;
	}
  }
}
