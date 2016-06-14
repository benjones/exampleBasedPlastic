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


struct Sphere {
   double center[3];
   double radius;
};

int main(int argc, char** argv) {
   if (argc != 6) {
      std::cerr << "Usage: " << argv[0] << " [grid_res] [min_radius] [in.obj] [spheres.txt] [signed_dist.vti]" << std::endl;
      std::cerr << "[grid_res] is an integer dictating the number of voxels along the minimum dimension of the bounding box of the input" << std::endl;
      std::cerr << "[min_radius] is a scale parameter for the radius of the minimum sphere used, scaled by side length of a voxel, 1.0 = side of voxel, <1.0 might be dangerous, >1.0 uses larger spheres, recommended is 2.0" << std::endl;
      std::cerr << "[in.obj] is the input obj file" << std::endl;
      std::cerr << "[spheres.txt] is the name of the output spheres file" << std::endl;
      std::cerr << "[signed_dist.vti] is the name of VTK image file that stores the signed distance field" << std::endl;

      exit(1);
   }


   //reader in input
   vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
   reader->SetFileName(argv[3]);
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
   grid_res[min_dim] = atoi(argv[1]);
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
   std::cout << "Input model is: " << argv[3] << std::endl;
   std::cout << "  Enlarged Bounding box dimensions: " << stretch[0] << " x " << stretch[1] << " x " << stretch[2] << std::endl;
   std::cout << "  Using grid resolution: " << grid_res[0] << " x " << grid_res[1] << " x " << grid_res[2] << std::endl;
   std::cout << std::endl;

   std::cout << "Building the signed distance field" << std::endl;

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
   std::cout << "Writing the signed distance field to: " << argv[5] << std::endl;
   std::cout << std::endl;

   vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
   writer->SetFileName(argv[5]);
   writer->SetInputData(image);
   writer->Update();

   std::cout << "Computing the sphere packing" << std::endl;

   //finally, do the sphere packing
   //use a negative threshold to only work with inside
   double THRESHOLD = -1.0 * atof(argv[2]) * min_spacing;
   //OFFSET = 0 will just place the spheres based on input data
   //OFFSET > 0 will push them out based on radius, and cover parts of
   //the exterior
   double OFFSET = 0.0;
   int count = 0;

   std::vector<Sphere> spheres;

   while (min_sphere.radius < THRESHOLD) {
      //insert the min_sphere
      vtkSmartPointer<vtkSphere> sphereDistance = vtkSmartPointer<vtkSphere>::New();
      //negate the radius, account for OFFSET before inserting
      min_sphere.radius = (1.0+OFFSET)*min_sphere.radius;
      spheres.push_back(min_sphere);

      sphereDistance->SetRadius(-1.0*min_sphere.radius);
      sphereDistance->SetCenter(min_sphere.center);


      //std::cout << "Sphere #" << count << ": " << min_sphere.center[0] << " " << min_sphere.center[1] << " " << min_sphere.center[2] << " " << -min_sphere.radius << std::endl; 

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
//       vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
//       char filename[1024];
//       sprintf(filename, "sphere.%04d.vti", count);
//       writer->SetFileName(filename);
//       writer->SetInputData(image);
//       writer->Update();
      count++;

   }

   std::cout << "Writing out the sphere file: " << argv[4] << " with " << spheres.size() << " spheres." << std::endl;

   //write out the spheres
   ofstream spheres_out;
   spheres_out.open(argv[4]);

   spheres_out << spheres.size() << std::endl;
   for (int i=0; i<spheres.size(); i++) {
      spheres_out << spheres[i].center[0] << " ";
      spheres_out << spheres[i].center[1] << " ";
      spheres_out << spheres[i].center[2] << " ";
      spheres_out << spheres[i].radius << std::endl;
   }

   spheres_out.close();

   return 0;
}

