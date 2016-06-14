#include <vtkImplicitPolyDataDistance.h>
#include <vtkSphere.h>
#include <vtkOBJReader.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkXMLImageDataWriter.h>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>

#include "parse_json.h"

struct Sphere {
   double center[3];
   double radius;
};

int main(int argc, char** argv) {
   if (argc != 2) {
      std::cerr << "Usage: " << argv[0] << " config.json" << std::endl;

      exit(1);
   }


   json_info info;
   parse_json(argv[1], info);

   //reader in input
   vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
   reader->SetFileName(info.obj_file.c_str());
   reader->Update();

   //create implicitDistance filter
   vtkSmartPointer<vtkImplicitPolyDataDistance> implicitDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
   vtkSmartPointer<vtkPolyData> poly_data = reader->GetOutput();
   implicitDistance->SetInput(poly_data);


   //First, generate a list of positions
   vtkSmartPointer<vtkPoints> positions = vtkSmartPointer<vtkPoints>::New();

   //this is lazy -- just user the imagedata class to generate the points
   //after specifying the Origin and Spacing

   vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
   int grid_res[3] = {-1,-1,-1};
   double bounds[6];
   poly_data->ComputeBounds();
   poly_data->GetBounds(bounds);

   //expand the bounds by 10%
   double len = bounds[1]-bounds[0];
   bounds[0] -= 0.05*len;
   bounds[1] += 0.05*len;
   len = bounds[3]-bounds[2];
   bounds[2] -= 0.05*len;
   bounds[3] += 0.05*len;
   len = bounds[5]-bounds[4];
   bounds[4] -= 0.05*len;
   bounds[5] += 0.05*len;

   //figure out the dimensions, use uniform scaling to adjust the grid
   double stretch[3];
   stretch[0] = bounds[1]-bounds[0];
   stretch[1] = bounds[3]-bounds[2];
   stretch[2] = bounds[5]-bounds[4];
   double min_stretch = stretch[0];
   int min_dim = 0;
   if (stretch[1] < min_stretch) {
      min_stretch = stretch[1];
      min_dim = 1;
   }
   if (stretch[2] < min_stretch) {
      min_stretch = stretch[2];
      min_dim = 2;
   }

   //best set of integer units
   grid_res[min_dim] = info.grid_res;
   grid_res[(min_dim+1)%3] = int(stretch[(min_dim+1)%3]/stretch[min_dim]*grid_res[min_dim]);
   grid_res[(min_dim+2)%3] = int(stretch[(min_dim+2)%3]/stretch[min_dim]*grid_res[min_dim]);


   image->SetDimensions(grid_res);
   image->SetOrigin(bounds[0], bounds[2], bounds[4]);
   double spacing[3] = {0};
   double min_spacing = 9e99;
   for (int i=0; i<3; i++) {
      if (grid_res[i] > 1) {
         spacing[i] = (bounds[2*i+1]-bounds[2*i]) / (grid_res[i]-1);
      } else {
         spacing[i] = 1;
      }

      if (spacing[i] < min_spacing) {
         min_spacing = spacing[i];
      }
   }
   image->SetSpacing(spacing);

   for (int i=0; i<image->GetNumberOfPoints(); i++) {
      positions->InsertNextPoint(image->GetPoint(i));
   }

   //populate the grid
   vtkSmartPointer<vtkDoubleArray> signed_distance = vtkSmartPointer<vtkDoubleArray>::New();
   signed_distance->SetNumberOfComponents(1);
   signed_distance->SetNumberOfTuples(positions->GetNumberOfPoints());


   //compute the distance to the surface
   //also compute the first sphere to pack
   if (info.verbosity > 1) {
      std::cout << "Input model is: " << info.obj_file << std::endl;
      std::cout << "  Enlarged Bounding box dimensions: " << stretch[0] << " x " << stretch[1] << " x " << stretch[2] << std::endl;
      std::cout << "  Using grid resolution: " << grid_res[0] << " x " << grid_res[1] << " x " << grid_res[2] << std::endl;
      std::cout << std::endl;

      std::cout << "Building the signed distance field" << std::endl;
   }

   Sphere min_sphere;
   min_sphere.radius = 0;
   min_sphere.center[0] = 0;
   min_sphere.center[1] = 0;
   min_sphere.center[2] = 0;
   for (int i=0; i<positions->GetNumberOfPoints(); i++) {
      double pos[3];
      positions->GetPoint(i, pos);
      double v = implicitDistance->EvaluateFunction(pos);

      signed_distance->SetValue(i,v);

      if (v < min_sphere.radius) {
         min_sphere.radius = v;
         min_sphere.center[0] = pos[0];
         min_sphere.center[1] = pos[1];
         min_sphere.center[2] = pos[2];
      }
   }

   //attach the data to the image
   image->GetPointData()->SetScalars(signed_distance);

   //write the signed distance field image
   if (info.vti_file.length() > 0) {
      if (info.verbosity > 1) {
         std::cout << "Writing the signed distance field to: " << info.vti_file << std::endl;
         std::cout << std::endl;
      }

      vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
      writer->SetFileName(info.vti_file.c_str());
      writer->SetInputData(image);
      writer->Update();
   }

   if (info.verbosity > 1) {
      std::cout << "Computing the sphere packing" << std::endl;
   }

   //finally, do the sphere packing
   //use a negative threshold to only work with inside
   //scale by input parameter to adjust fewest number of spheres
   double THRESHOLD = -1.0 * info.min_radius_scale * min_spacing;
   int count = 0;

   std::vector<Sphere> spheres;

   while (min_sphere.radius < THRESHOLD && count < info.max_spheres) {
      //insert the min_sphere
      vtkSmartPointer<vtkSphere> sphereDistance = vtkSmartPointer<vtkSphere>::New();
      //negate the radius, account for OFFSET before inserting
      //default info.offset is 0, which does no extra scaling
      min_sphere.radius = (1.0+info.offset)*min_sphere.radius;
      spheres.push_back(min_sphere);

      sphereDistance->SetRadius(-1.0*min_sphere.radius);
      sphereDistance->SetCenter(min_sphere.center);


      if (info.verbosity > 1) {
        std::cout << "Sphere #" << count << ": " << min_sphere.center[0] << " " << min_sphere.center[1] << " " << min_sphere.center[2] << " " << -min_sphere.radius << std::endl; 
      }

      //go through each point do a min subtract
      min_sphere.radius = 0;
      min_sphere.center[0] = 0;
      min_sphere.center[1] = 0;
      min_sphere.center[2] = 0;
      for (int i=0; i<positions->GetNumberOfPoints(); i++) {
         double pos[3];
         positions->GetPoint(i, pos);
         double v = signed_distance->GetValue(i);
         double s_v = sphereDistance->EvaluateFunction(pos);
         double new_v = v;

         if (s_v <= 0) {
            //inside sphere, zero out
            signed_distance->SetValue(i,0);
         } else if (v >= 0) {
            //skip, outside
         } else {
            //replace with minimum
            signed_distance->SetValue(i,std::min(v,s_v));
         }

         if (signed_distance->GetValue(i) < min_sphere.radius) {
            min_sphere.radius = signed_distance->GetValue(i);
            min_sphere.center[0] = pos[0];
            min_sphere.center[1] = pos[1];
            min_sphere.center[2] = pos[2];
         }

      }

      //attach the data to the image
      image->GetPointData()->SetScalars(signed_distance);

      //write the intermediate images
      if (info.verbosity > 2) {
         vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
         char filename[1024];
         sprintf(filename, "sphere.%04d.vti", count);
         writer->SetFileName(filename);
         writer->SetInputData(image);
         writer->Update();
      }
      count++;

   }

   if (info.spheres_file.length() > 0) {
      if (info.verbosity > 0) {
         std::cout << "Writing out the sphere file: " << info.spheres_file << " with " << spheres.size() << " spheres." << std::endl;
      }

      //write out the spheres
      ofstream spheres_out;
      spheres_out.open(info.spheres_file.c_str());

      spheres_out << spheres.size() << std::endl;
      for (int i=0; i<spheres.size(); i++) {
         spheres_out << spheres[i].center[0] << " ";
         spheres_out << spheres[i].center[1] << " ";
         spheres_out << spheres[i].center[2] << " ";
         spheres_out << -spheres[i].radius << std::endl;
      }

      spheres_out.close();
   } else {
      if (info.verbosity > 1) {
         std::cout << "No spheres filename was specified, but we packed: " << spheres.size() << " spheres" << std::endl;
      }
   }

   return 0;
}

