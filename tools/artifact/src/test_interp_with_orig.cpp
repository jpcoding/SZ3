#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "critical_points.hpp"
#include "posterization.hpp"
#include <array>
#include <data_preprocess.hpp>
#include <iostream>
#include <ostream>
#include <vector>

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
  float omin, omax;
  float dmin, dmax;
  float orange = normalization(odata, omin, omax);

  float drange = normalization_range(ddata, omin, omax, orange);
  // float drange = normalization(ddata, dmin, dmax);

  SZ::Timer timer;

  // constrcut the critical point map;
  CriticalPointsCalculator odata_cp(odata.data(), N, global_dimensions.data());

  CriticalPointsCalculator ddata_cp(ddata.data(), N, global_dimensions.data());

  double detection_time = 0;
  timer.start();
  std::vector<int> odata_cp_map = odata_cp.get_critical_points_map();

  std::vector<int> ddata_cp_map = ddata_cp.get_critical_points_map();
  detection_time += timer.stop("get_critical_points_map");

  float local_range = 1e-3;

  if (argc == (4 + N + 1)) {
    local_range = atof(argv[argc - 1]);
  }
  odata_cp.set_local_range(local_range);
  ddata_cp.set_local_range(local_range);
  // local_range = 1e-5

  float interp_threshold = 1e-2;

  std::cout << "local_range = " << local_range << std::endl;

  std::vector<int> match_index;
  std::vector<int> xpaddings;
  std::vector<int> ypaddings;
  std::vector<int> zpaddings;
  std::vector<float> interp_errors;

  if (N == 2) {
    match_index.reserve(num_elements);
    xpaddings.reserve(num_elements);
    ypaddings.reserve(num_elements);
    int dxpadding;
    int dypadding;
    int oxpadding;
    int oypadding;
    int match_count = 0;
    interp_threshold = 1e-3;
    timer.start();
    for (int i = 0; i < num_elements; i++) {
      int idx = i % global_dimensions[0];
      int idy = i / global_dimensions[0];
      // if ((idx & 1) || (idy & 1)) {
      //   continue;
      // }
      if (ddata_cp_map[i] == 4 || ddata_cp_map[i] == -4) {

        float interp_error = 0;
        bool dmatch =
            ddata_cp.try_match2d_interp(ddata_cp_map, i, dxpadding, dypadding,
                                        interp_error, interp_threshold);
        //     bool dmatch =
        // ddata_cp.try_match2d(ddata_cp_map, i, dxpadding, dypadding);

        if (dmatch) {
          // bool omatch =
          //     odata_cp.try_match2d(odata_cp_map, i, oxpadding, oypadding);
          bool omatch = (ddata_cp_map[i] == odata_cp_map[i]);

          // if (i==3204420)
          // {
          //   std::cout<< "dmatch: " << dmatch << std::endl;
          //   std::cout<< "omatch: " << omatch << std::endl;
          // }
          if (!omatch) {
            match_index.push_back(i);
            xpaddings.push_back(dxpadding);
            ypaddings.push_back(dypadding);
            interp_errors.push_back(interp_error);
            match_count++;
          }
        }
      }
    }
    detection_time += timer.stop("pattern match");
    std::cout << "match_count: " << match_count << std::endl;
    std::cout << "detection_time: " << detection_time << std::endl;

    SZ::writefile("match_index.dat", match_index.data(), match_index.size());
    SZ::writefile("xpaddings.dat", xpaddings.data(), xpaddings.size());
    SZ::writefile("ypaddings.dat", ypaddings.data(), ypaddings.size());
    SZ::writefile("interp_errors.dat", interp_errors.data(),
                  interp_errors.size());

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
    interp_threshold = 1e-3;

    timer.start();

    for (int i = 0; i < num_elements; i++) {

      int idx = i % global_dimensions[0];
      int idy = (i / global_dimensions[0]) % global_dimensions[1];
      int idz = i / (global_dimensions[0] * global_dimensions[1]);
      // if ((idx & 1) || (idy & 1) || (idz & 1)) {
      //   continue;
      // }
      if (ddata_cp_map[i] == 6 || ddata_cp_map[i] == -6) {


        float interp_error=0;
        // bool dmatch = ddata_cp.try_match3d(ddata_cp_map, i, dxpadding,
        //  dypadding, dzpadding);
        bool dmatch = ddata_cp.try_match3d_interp(
            ddata_cp_map, i, dxpadding, dypadding, dzpadding, interp_error,
            interp_threshold);
            
        if (dmatch) {
          // bool omatch = odata_cp.try_match3d(odata_cp_map, i, oxpadding,
                                            //  oypadding, ozpadding);

          bool omatch = (ddata_cp_map[i] == odata_cp_map[i]);
          // std::cout << "dmatch: " << i << std::endl;
          if (!omatch) {
            match_index.push_back(i);
            xpaddings.push_back(dxpadding);
            ypaddings.push_back(dypadding);
            zpaddings.push_back(dzpadding);
            interp_errors.push_back(interp_error);
            match_count++;
          }
        }
      }
    }
    detection_time += timer.stop("pattern match");
    std::cout << "match_count: " << match_count << std::endl;
    std::cout << "detection_time = " << detection_time << std::endl;
    SZ::writefile("match_index.dat", match_index.data(), match_index.size());
    SZ::writefile("xpaddings.dat", xpaddings.data(), xpaddings.size());
    SZ::writefile("ypaddings.dat", ypaddings.data(), ypaddings.size());
    SZ::writefile("zpaddings.dat", zpaddings.data(), zpaddings.size());
    SZ::writefile("cp_map.dat", ddata_cp_map.data(), num_elements);
    SZ::writefile("cp_map_orig.dat", odata_cp_map.data(), num_elements);
    SZ::writefile("interp_errors.dat", interp_errors.data(),
                  interp_errors.size());
  }

  return 0;
}