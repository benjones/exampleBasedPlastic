#include <iostream>
#include <limits>
#include <set>
#include <algorithm>
#include <cmath>
#include <chrono>


enum State { far, considered, accepted };

void fmm(double* dist, State* states, int N, std::set<int>& ids) {
   while (!ids.empty()) {
      //first get the min out of the id list;
      std::set<int>::iterator id_iter = ids.begin();
      double min_dist = std::numeric_limits<double>::max();
      int cid = -1;
      while (id_iter != ids.end()) {
         int hid = *id_iter;

         if (dist[hid] < min_dist) {
            min_dist = dist[hid];
            cid = hid;
         }

         id_iter++;
      }

      int cy = cid/N;
      int cx = cid%N;


      std::cout << cx << " " << cy << " " << dist[cid] << std::endl;

      ids.erase(cid);
   
      if (states[cid] == considered) {
         states[cid] = accepted;
   
         //check all neighbors
         int nx, ny;
         int nid;

         //first cx-1,cy
         nx = cx-1; ny = cy;
         nid = ny*N + nx;
         if (nx >= 0 && states[nid] != accepted) {
            //estimate from neighbors
            double dx = std::numeric_limits<double>::max();
            if (nx < N-1 && dist[nid+1] < dx) {
               dx = dist[nid+1];
            }
            if (nx > 0 && dist[nid-1] < dx) {
               dx = dist[nid-1];
            }
            double dy = std::numeric_limits<double>::max();
            if (ny < N-1 && dist[nid+N] < dy) {
               dy = dist[nid+N];
            }
            if (ny > 0 && dist[nid-N] < dy) {
               dy = dist[nid-N];
            }
            double delta = 2 - (dx-dy)*(dx-dy);
            double d_trial = std::min(dx+1,dy+1);
            if (delta >= 0) {
               d_trial = (dx + dy + sqrt(delta)) / 2;
            }
            if (d_trial < dist[nid]) {
               dist[nid] = d_trial;
            }
            if (states[nid] == far) {
               states[nid] = considered;
               ids.insert(ny*N + nx);
            }
         }

         //second cx+1,cy
         nx = cx+1; ny = cy;
         nid = ny*N + nx;
         if (nx < N && states[nid] != accepted) {
            //estimate from neighbors
            double dx = std::numeric_limits<double>::max();
            if (nx < N-1 && dist[nid+1] < dx) {
               dx = dist[nid+1];
            }
            if (nx > 0 && dist[nid-1] < dx) {
               dx = dist[nid-1];
            }
            double dy = std::numeric_limits<double>::max();
            if (ny < N-1 && dist[nid+N] < dy) {
               dy = dist[nid+N];
            }
            if (ny > 0 && dist[nid-N] < dy) {
               dy = dist[nid-N];
            }
            double delta = 2 - (dx-dy)*(dx-dy);
            double d_trial = std::min(dx+1,dy+1);
            if (delta >= 0) {
               d_trial = (dx + dy + sqrt(delta)) / 2;
            }
            if (d_trial < dist[nid]) {
               dist[nid] = d_trial;
            }
            if (states[nid] == far) {
               states[nid] = considered;
               ids.insert(ny*N + nx);
            }
         }


         //third cx,cy-1
         nx = cx; ny = cy-1;
         nid = ny*N + nx;
         if (ny >= 0 && states[nid] != accepted) {
            //estimate from neighbors
            double dx = std::numeric_limits<double>::max();
            if (nx < N-1 && dist[nid+1] < dx) {
               dx = dist[nid+1];
            }
            if (nx > 0 && dist[nid-1] < dx) {
               dx = dist[nid-1];
            }
            double dy = std::numeric_limits<double>::max();
            if (ny < N-1 && dist[nid+N] < dy) {
               dy = dist[nid+N];
            }
            if (ny > 0 && dist[nid-N] < dy) {
               dy = dist[nid-N];
            }
            double delta = 2 - (dx-dy)*(dx-dy);
            double d_trial = std::min(dx+1,dy+1);
            if (delta >= 0) {
               d_trial = (dx + dy + sqrt(delta)) / 2;
            }
            if (d_trial < dist[nid]) {
               dist[nid] = d_trial;
            }
            if (states[nid] == far) {
               states[nid] = considered;
               ids.insert(ny*N + nx);
            }
         }


         //fourth cx,cy+1
         nx = cx; ny = cy+1;
         nid = ny*N + nx;
         if (ny < N && states[nid] != accepted) {
            //estimate from neighbors
            double dx = std::numeric_limits<double>::max();
            if (nx < N-1 && dist[nid+1] < dx) {
               dx = dist[nid+1];
            }
            if (nx > 0 && dist[nid-1] < dx) {
               dx = dist[nid-1];
            }
            double dy = std::numeric_limits<double>::max();
            if (ny < N-1 && dist[nid+N] < dy) {
               dy = dist[nid+N];
            }
            if (ny > 0 && dist[nid-N] < dy) {
               dy = dist[nid-N];
            }
            double delta = 2 - (dx-dy)*(dx-dy);
            double d_trial = std::min(dx+1,dy+1);
            if (delta >= 0) {
               d_trial = (dx + dy + sqrt(delta)) / 2;
            }
            if (d_trial < dist[nid]) {
               dist[nid] = d_trial;
            }
            if (states[nid] == far) {
               states[nid] = considered;
               ids.insert(ny*N + nx);
            }
         }

         
      } else {
         //someone else updated this node
      }
   }
}

int main() {
   const int N = 11;
   double dist[N*N];
   State states[N*N];

   //initialize
   for (int y=0; y<N; y++) {
      for (int x=0; x<N; x++) {
         dist[y*N+x] = std::numeric_limits<double>::infinity();
         states[y*N+x] = far;
      }
   }

   //set boundary condition
   std::set< int > ids;
   
   for (int x=0; x<N; x++) {
      //layer outside of the surface
      int cy = N/2;
      int cx = x;
      int cid = N*cy + cx;
      dist[cid] = -0.5;
      states[cid] = accepted;

      //ids.insert(cid);

      //layer inside of the surface
      cy = cy+1;
      cid = N*cy + cx;
      dist[cid] = +0.5;
      states[cid] = considered;

      ids.insert(cid);

   }


   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   fmm(dist, states, N, ids);


   end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end-start;

   std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

   return 0;
}


