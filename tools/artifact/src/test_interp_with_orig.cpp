#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "critical_points.hpp"
#include "posterization.hpp"
#include <array>
#include <iostream>
#include <vector>
#include <data_preprocess.hpp>

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cout << "Usage: ./interp_detect <input_file> <N> <dim1 dim2 ...> "
              << std::endl;
    return 0;
  }
  int N = atoi(argv[3]);
  int num_elements = 1;
  std::vector<int> global_dimensions;
  for (int i = 0; i < N; i++) {
    global_dimensions.push_back(atoi(argv[4 + i]));
    num_elements *= global_dimensions[i];
  }

  // float *data = new float[num_elements];
  std::vector<float> odata(num_elements);
  std::vector<float> ddata(num_elements);
  SZ::readfile<float>(argv[1], num_elements, odata.data());
  SZ::readfile<float>(argv[2], num_elements, ddata.data());

  // normalize the data;
  float omin,omax;
    float dmin,dmax;
float orange = normalization(odata, omin, omax);
float dorange = normalization(ddata, dmin, dmax);

  // constrcut the critical point map;
  CriticalPointsCalculator odata_cp(odata.data(), N, global_dimensions.data());
  std::vector<int> odata_cp_map = odata_cp.get_critical_points_map();

  CriticalPointsCalculator ddata_cp(ddata.data(), N, global_dimensions.data());
  std::vector<int> ddata_cp_map = ddata_cp.get_critical_points_map();

odata_cp.set_local_range(1e-2);
ddata_cp.set_local_range(1e-2);

  std::vector<int> match_index;
  std::vector<int> xpaddings;
  std::vector<int> ypaddings;
  std::vector<int> zpaddings;

SZ::writefile("cp_map.dat", ddata_cp_map.data(),num_elements);
std::cout<<ddata_cp_map[0]<<std::endl;

  if (N == 2) {
    match_index.reserve(num_elements);
    xpaddings.reserve(num_elements);
    ypaddings.reserve(num_elements);
    int dxpadding;
    int dypadding;
    int oxpadding;
    int oypadding;
    int match_count = 0;
    for (int i = 0; i < num_elements; i++) {
      if (ddata_cp_map[i] == 4 || ddata_cp_map[i] == -4) {

        bool dmatch =
            ddata_cp.try_match2d(ddata_cp_map,i, dxpadding, dypadding);
        
        if (dmatch) {
          bool omatch =
              odata_cp.try_match2d(odata_cp_map, i, oxpadding, oypadding);
          if (!omatch) {
            match_index.push_back(i);
            xpaddings.push_back(dxpadding);
            ypaddings.push_back(dypadding);
            match_count++;
          }
        }
      }
    }
    std::cout << "match_count: " << match_count << std::endl;

  } else if (N == 3) {
    match_index.reserve(num_elements);
    xpaddings.reserve(num_elements);
    ypaddings.reserve(num_elements);
    zpaddings.reserve(num_elements);
    int dxpadding;
    int dypadding;
    int dzpadding;
    int oxpadding;
    int oypadding;
    int ozpadding;
    int match_count = 0;
    for (int i = 0; i < num_elements; i++) {
      if (ddata_cp_map[i] == 6 || ddata_cp_map[i] == -6) {

        bool dmatch = ddata_cp.try_match3d(ddata_cp_map,i, dxpadding,
                                           dypadding, dzpadding);
        
        std::cout << "dmatch: " << dmatch << std::endl;

        if (dmatch) {
          bool omatch = odata_cp.try_match3d(odata_cp_map ,i, oxpadding,
                                             oypadding, ozpadding);
          if (!omatch) {
            match_index.push_back(i);
            xpaddings.push_back(dxpadding);
            ypaddings.push_back(dypadding);
            zpaddings.push_back(dzpadding);
            match_count++;
          }
        }
      }
    }
    std::cout << "match_count: " << match_count << std::endl;
  }

  // posterization.evaluation( segmentation_map);

  // SZ::writefile("segmentation_map.dat", segmentation_map.data(),
  // segmentation_map.size());

  return 0;
}