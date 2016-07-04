#include <iostream>
#include <cstring>

int main(int argc, char *argv[]) {
  char command[300];
  
  if (argc != 8) {
	std::cout<<"usage: doAll path_to_example_files path_to_input_sphere_pack_file path_to_boneOptimizer path_to_destruction path_to_hiresSkinner nexample outputfile"<<std::endl;
  }
  
  sprintf (command, "cp %s/* .", argv[1]);
  std::cout<<command<<std::endl;
  system (command);
  
  sprintf (command, "cp %s .", argv[2]);
  std::cout<<command<<std::endl;
  system (command);
  
  sprintf (command, "%s/makeVolumePly mesh.mesh coarseWeights.dmat", argv[3]);
  std::cout<<command<<std::endl;
  system (command);

  for (unsigned int i=0; i<atoi(argv[6]); i++) {
	sprintf (command, "%s/boneOptimizer volumeMesh.ply volumeWeights.dmat -e exampleGraph.txt %i examplePose.%02d.ply", argv[3], i, i);
	std::cout<<command<<std::endl;
	system(command);

	sprintf (command, "%s/sphere_pack/skinSpheres/skinSpheres ./spheres.txt mesh.mesh examplePose.%02d.ply exampleSpheres.%02d.txt", argv[4], i, i);
	std::cout<<command<<std::endl;
	system(command);
	
	sprintf(command, "%s/skinSingleMesh mesh.obj mesh.mesh examplePose.%02d.ply skinnedExample.%02d.ply", argv[5], i, i);
	std::cout<<command<<std::endl;
	system(command);
	
	sprintf(command, "%s/sphere_pack/vtk/build/sphere_resize skinnedExample.%02d.ply exampleSpheres.%02d.txt exampleSpheresResized.%02d.txt", argv[4], i, i, i);
	std::cout<<command<<std::endl;
	system(command);
  }

  sprintf(command, "%s/sphere_pack/combineSpheres/combineSpheres exampleSpheresResized.%s.txt %i %s", argv[4], "%02d", atoi(argv[6]), argv[7]);
  std::cout<<command<<std::endl;
  system(command);
}
