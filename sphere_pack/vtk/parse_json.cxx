#include "parse_json.h"

#include <fstream>
#include <iostream>
#include <json/json.h>
#include <climits>


void parse_json(char* filename, json_info& info) {
   std::ifstream ins(filename);
  
   Json::Value root;
   Json::Reader jReader;

   if(!jReader.parse(ins, root)){
      std::cout << "couldn't read input file: " << filename << '\n'
         << jReader.getFormattedErrorMessages() << std::endl;
      exit(1);
   }

   info.grid_res = root.get("grid_res", 10).asInt();
   info.min_radius_scale = root.get("min_radius_scale", 2.0).asDouble();
   info.offset = root.get("offset", 0.0).asDouble();
   info.max_spheres = root.get("max_spheres", INT_MAX).asInt();

   //0 = none
   //1 = debug messages
   //2 = each sphere
   //3 = write debug files
   info.verbosity = root.get("verbosity", 1).asInt();

   info.obj_file = root.get("obj_file", "").asString();
   info.spheres_file = root.get("spheres_file", "").asString();
   info.vti_file = root.get("vti_file", "").asString();
}


