#include <vtkImplicitPolyDataDistance.h>
#include <vtkSphere.h>
#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>
#include <vtkPLYReader.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCommand.h>
#include <vtkCallbackCommand.h>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>

#include "parse_json.h"

struct Sphere {
  double center[3];
  double radius;
};




void errorCallbackFunc ( vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
  std::cerr << "An error occurred in the vtkOBJReader: " << std::endl;
  std::cerr << reinterpret_cast<char *>(callData) << std::endl;
  
  exit(2);
}

void readSphereFile(char *fname, std::vector<Sphere> &spheres) {
  int nspheres;
  Sphere s;
  std::ifstream in(fname, std::ios::in);
  
  in>>nspheres;
  for (unsigned int i=0; i<nspheres; i++) {
	in>>s.center[0]>>s.center[1]>>s.center[2]>>s.radius;
	spheres.push_back(s);
  }
}

int main(int argc, char *argv[]) {
	if (argc != 4) {
      std::cerr << "Usage: " << argv[0] << " mesh.ply input_sphere_file output_sphere_file" << std::endl;
      exit(1);
	}

   vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
   reader->SetFileName(argv[1]);
   //set up error handling in case the obj file is bad
   vtkSmartPointer<vtkCallbackCommand> errorCallback = vtkSmartPointer<vtkCallbackCommand>::New();
   errorCallback->SetCallback ( errorCallbackFunc );
   reader->AddObserver ( vtkCommand::ErrorEvent, errorCallback );

   //read the obj_file
   reader->Update();

   //create implicitDistance filter
   vtkSmartPointer<vtkImplicitPolyDataDistance> implicitDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
   vtkSmartPointer<vtkPolyData> poly_data = reader->GetOutput();
   implicitDistance->SetInput(poly_data);

   std::vector<Sphere> spheres;
   readSphereFile(argv[2], spheres);
   for (int i=0; i<spheres.size(); i++) {
	 double d = implicitDistance->EvaluateFunction(spheres[i].center);
	 std::cout<<spheres[i].radius<<" "<<d<<std::endl;
	 spheres[i].radius = d;
   }

   //write out the spheres
   ofstream spheres_out;
   spheres_out.open(argv[3]);
   spheres_out << spheres.size() << std::endl;
   for (int i=0; i<spheres.size(); i++) {
	 spheres_out << spheres[i].center[0] << " ";
	 spheres_out << spheres[i].center[1] << " ";
	 spheres_out << spheres[i].center[2] << " ";
	 spheres_out << -spheres[i].radius << std::endl;
   }

   spheres_out.close();

   return 0;
}

