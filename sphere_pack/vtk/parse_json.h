#ifndef PARSE_JSON_H
#define PARSE_JSON_H

#include <string>

class json_info {
   public:
      int grid_res;
      double min_radius_scale;
      double max_radius_scale;
      double offset;
      int max_spheres;
      int verbosity;

      std::string obj_file;
      std::string spheres_file;
      std::string vti_file;
};


void parse_json(char* filename, json_info& info);


#endif
