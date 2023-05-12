#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "critical_points.hpp"
#include "posterization.hpp"
#include <algorithm>
#include <array>
#include <iostream>
#include <vector>

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cout << "Usage: ./interp_detect <input_file> <N> <dim1 dim2 ...> "
              << std::endl;
    return 0;
  }
  int N = atoi(argv[2]);
  int num_elements = 1;
  std::vector<int> global_dimensions;
  for (int i = 0; i < N; i++) {
    global_dimensions.push_back(atoi(argv[3 + i]));
    num_elements *= global_dimensions[i];
  }

  // float *data = new float[num_elements];
  std::vector<float> ddata(num_elements);
  SZ::readfile<float>(argv[1], num_elements, ddata.data());
  // constrcut the critical point map;
  CriticalPointsCalculator cp_calculator(ddata.data(), N,
                                         global_dimensions.data());
  std::vector<int> critical_points_map =
      cp_calculator.get_critical_points_map();
  SZ::writefile<int>("critical_points_map.dat", critical_points_map.data(),
                     num_elements);
  std::vector<int> cp_indexes;
  std::vector<int> cp_xpaddings;
  std::vector<int> cp_ypaddings;
  std::vector<int> cp_zpaddings;
  int match_count = 0;
  if (N == 2) {
    cp_calculator.pattern_match_global(cp_indexes, cp_xpaddings, cp_ypaddings,
                                       match_count);
    SZ::writefile<int>("cp_indexes.dat", cp_indexes.data(), cp_indexes.size());
    SZ::writefile<int>("cp_xpaddings.dat", cp_xpaddings.data(),
                       cp_xpaddings.size());
    SZ::writefile<int>("cp_ypaddings.dat", cp_ypaddings.data(),
                       cp_ypaddings.size());
  } else if (N == 3) {
    cp_calculator.pattern_match_global(cp_indexes, cp_xpaddings, cp_ypaddings,
                                       cp_zpaddings, match_count);
    SZ::writefile<int>("cp_indexes.dat", cp_indexes.data(), cp_indexes.size());
    SZ::writefile<int>("cp_xpaddings.dat", cp_xpaddings.data(),
                       cp_xpaddings.size());
    SZ::writefile<int>("cp_ypaddings.dat", cp_ypaddings.data(),
                       cp_ypaddings.size());
    SZ::writefile<int>("cp_zpaddings.dat", cp_zpaddings.data(),
                       cp_zpaddings.size());
  } else {
    std::cout << "N = " << N << " is not supported!" << std::endl;
    return 0;
  }
  std::cout << "match_count = " << match_count << std::endl;

  // interpolation walk to evaluate the interpolation error arount the critical
  // points
  InterpolationWalk<float> interpolation_walk(ddata.data(), N,
                                              global_dimensions.data());

  if (N == 2) {
    std::vector<float> interpolation_error;
    for (int i = 0; i < cp_indexes.size(); i++) {
      int cp_index = cp_indexes[i];
      int cp_xpadding = cp_xpaddings[i];
      int cp_ypadding = cp_ypaddings[i];
      int idx = cp_index % global_dimensions[0];
      int idy = cp_index / global_dimensions[0];
      float error = 0;
      //   int count = 0;
      for (int j = 1; j <= cp_xpadding; j++) {
        for (int k = 1; k <= cp_ypadding; k++) {
          error += interpolation_walk.interp_walk(idx + j, idy + k);
          error += interpolation_walk.interp_walk(idx + j, idy - k);
          error += interpolation_walk.interp_walk(idx - j, idy + k);
          error += interpolation_walk.interp_walk(idx - j, idy - k);
          //   count += 4;
        }
      }
      //   if (i ==1)
      //   {
      //         std::cout << "count = " << count << std::endl;
      //         std::cout << "padding size" << 4*cp_xpadding*cp_ypadding -1 <<
      //         std::endl;
      //   }
      interpolation_error.push_back(error / (4 * cp_xpadding * cp_ypadding));
    }

    float max_error = *std::max_element(interpolation_error.begin(),
                                        interpolation_error.end());
    float min_error = *std::min_element(interpolation_error.begin(),
                                        interpolation_error.end());
    std::cout << "min_error = " << min_error << std::endl;
    std::cout << "max_error = " << max_error << std::endl;
    SZ::writefile<float>("interpolation_error.dat", interpolation_error.data(),
                         interpolation_error.size());
    std::sort(interpolation_error.begin(), interpolation_error.end());
    float median = interpolation_error[interpolation_error.size() / 2];
    std::cout << "median = " << median << std::endl;
  } else if (N == 3) {
    std::vector<float> interpolation_error;
    for (int i = 0; i < cp_indexes.size(); i++) {
      int cp_index = cp_indexes[i];
      int cp_xpadding = cp_xpaddings[i];
      int cp_ypadding = cp_ypaddings[i];
      int cp_zpadding = cp_zpaddings[i];
      int idx = cp_index % global_dimensions[0];
      int idy = (cp_index / global_dimensions[0]) % global_dimensions[1];
      int idz = cp_index / (global_dimensions[0] * global_dimensions[1]);
      float error = 0;
      // int count = 0;
      for (int j = 1; j <= cp_xpadding; j++) {
        for (int k = 1; k <= cp_ypadding; k++) {
          for (int l = 1; l <= cp_zpadding; l++) {
            error += interpolation_walk.interp_walk(idx + j, idy + k, idz + l);
            error += interpolation_walk.interp_walk(idx + j, idy + k, idz - l);
            error += interpolation_walk.interp_walk(idx + j, idy - k, idz + l);
            error += interpolation_walk.interp_walk(idx + j, idy - k, idz - l);
            error += interpolation_walk.interp_walk(idx - j, idy + k, idz + l);
            error += interpolation_walk.interp_walk(idx - j, idy + k, idz - l);
            error += interpolation_walk.interp_walk(idx - j, idy - k, idz + l);
            error += interpolation_walk.interp_walk(idx - j, idy - k, idz - l);
            // count += 8;
          }
        }
      }
      // if (i ==1)
      // {
      //     std::cout << "count = " << count << std::endl;
      //     std::cout << "padding size" <<
      //     8*cp_xpadding*cp_ypadding*cp_zpadding << std::endl;
      // }
      interpolation_error.push_back(
          error / (8 * cp_xpadding * cp_ypadding * cp_zpadding));
    }
    float min_error = *std::min_element(interpolation_error.begin(),
                                        interpolation_error.end());
    float max_error = *std::max_element(interpolation_error.begin(),
                                        interpolation_error.end());

    std::cout << "min_error = " << min_error << std::endl;
    std::cout << "max_error = " << max_error << std::endl;
    SZ::writefile<float>("interpolation_error.dat", interpolation_error.data(),
                         interpolation_error.size());
    std::sort(interpolation_error.begin(), interpolation_error.end());
    float median = interpolation_error[interpolation_error.size() / 2];
    std::cout << "median = " << median << std::endl;
  }

  return 0;
}