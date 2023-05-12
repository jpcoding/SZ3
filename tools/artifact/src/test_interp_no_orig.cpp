#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "critical_points.hpp"
#include "posterization.hpp"
#include <algorithm>
#include <array>
#include <data_preprocess.hpp>
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

  float interp_threshold = 1e-2;

  if (argc == (3 + N + 1)) {
    interp_threshold = atof(argv[argc - 1]);
  }

  if (argc == (3 + N + 2)) {
    interp_threshold = atof(argv[argc - 1]);
  }
  std::cout<<"input interp_threshold = "<<interp_threshold<<std::endl;

  // float *data = new float[num_elements];
  std::vector<float> ddata(num_elements);
  SZ::readfile<float>(argv[1], num_elements, ddata.data());
  // normalize the data;
  float global_max;
  float global_min;
  float original_range = normalization(ddata, global_max, global_min);
  if(original_range == -1){
    std::cout<<"original range is 0, please check the input data!"<<std::endl;
    return 0;
  }
  std::cout<<"original range = "<<original_range<<std::endl;
interp_threshold = interp_threshold / original_range;
std::cout<<"normalized interp_threshold = "<<interp_threshold<<std::endl;
  
  // constrcut the critical point map;
  CriticalPointsCalculator cp_calculator(ddata.data(), N,
                                         global_dimensions.data());
  std::vector<int> critical_points_map =
      cp_calculator.get_critical_points_map();
  cp_calculator.set_local_range(1e-4);
  SZ::writefile<int>("critical_points_map.dat", critical_points_map.data(),
                     num_elements);
  std::vector<int> cp_indexes;
  std::vector<int> cp_xpaddings;
  std::vector<int> cp_ypaddings;
  std::vector<int> cp_zpaddings;
  std::vector<float> interp_errors;
  int match_count = 0;
  if (N == 2) {
    cp_calculator.pattern_match_global(cp_indexes, cp_xpaddings, cp_ypaddings,
                                       interp_errors, match_count,
                                       interp_threshold);
    SZ::writefile<int>("cp_indexes.dat", cp_indexes.data(), cp_indexes.size());
    SZ::writefile<int>("cp_xpaddings.dat", cp_xpaddings.data(),
                       cp_xpaddings.size());
    SZ::writefile<int>("cp_ypaddings.dat", cp_ypaddings.data(),
                       cp_ypaddings.size());
  } else if (N == 3) {
    cp_calculator.pattern_match_global(cp_indexes, cp_xpaddings, cp_ypaddings,
                                       cp_zpaddings, interp_errors, match_count,
                                       interp_threshold);
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

  return 0;
}