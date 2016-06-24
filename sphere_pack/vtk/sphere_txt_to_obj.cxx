#include <vtkSphere.h>
#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "vtkOBJWriter.h"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>

struct Sphere {
   double center[3];
   double radius;
};



int main(int argc, char** argv) {
   if (argc != 3) {
      std::cerr << "Usage: " << argv[0] << " spheres.txt spheres.obj" << std::endl;

      exit(1);
   }


   //write out t:he spheres
   ifstream spheres_in;
   std::vector<Sphere> spheres;
   spheres_in.open(argv[1]);

   int num_spheres;

   spheres_in >> num_spheres;

   for (int i=0; i<num_spheres; i++) {
      Sphere s;
      spheres_in >> s.center[0];
      spheres_in >> s.center[1];
      spheres_in >> s.center[2];
      spheres_in >> s.radius;
      spheres.push_back(s);

   }

   vtkSmartPointer<vtkAppendPolyData> sphere_poly_list = vtkSmartPointer<vtkAppendPolyData>::New();


   for (int i=0; i<spheres.size(); i++) {
      vtkSmartPointer<vtkSphereSource> sphere_source = vtkSmartPointer<vtkSphereSource>::New();
      sphere_source->SetRadius(spheres[i].radius);
      sphere_source->SetCenter(spheres[i].center);
      sphere_source->Update();
      vtkSmartPointer<vtkPolyData> sphere_poly = vtkSmartPointer<vtkPolyData>::New();
      sphere_poly->ShallowCopy(sphere_source->GetOutput());
      sphere_poly_list->AddInputData(sphere_poly);
      sphere_poly_list->Update();
   }

   vtkSmartPointer<vtkOBJWriter> writer = vtkSmartPointer<vtkOBJWriter>::New();
   writer->SetInputData(sphere_poly_list->GetOutput());
   writer->SetFileName(argv[2]);
   writer->Update();

   return 0;
}

