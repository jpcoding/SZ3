#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include <cstdlib>
#include <iostream>
#include <posterization.hpp>
#include <unordered_map>
#include <vector>
#include <data_preprocess.hpp>
int main(int argc, char **argv) {
  if (argc < 2) {
    std::cout << "Usage: ./interp_detect <input_file> <N> <dim1 dim2 ...> "
              << std::endl;
    return 0;
  }
  int N = atoi(argv[2]);
  int num_elements = 1;
  std::vector<int> global_dimensions;
  int i;
  for (i = 0; i < N; i++) {
    global_dimensions.push_back(atoi(argv[3 + i]));
    num_elements *= global_dimensions[i];
  }
  // std::cout << "i: " << argc << std::endl;
  std::cout << "num_elements: " << num_elements << std::endl;
  float threhshold = 1e-5;

  if (argc == (3 + N + 1)) {
    threhshold = atof(argv[argc - 1]);
  }

  if (argc == (3 + N + 2)) {
    threhshold = atof(argv[argc - 1]);
  }

  std::cout << "input threhshold: " << threhshold << std::endl;

  std::vector<float> data(num_elements);
  SZ::readfile<float>(argv[1], num_elements, data.data());
  float data_min, dat_max;
  float range = normalization( data, data_min, dat_max);
  threhshold = threhshold / range;
  std::cout << "factoirzed threshold: " << threhshold << std::endl;
  SZ::Timer timer;
  // construct posterization analyzer
  Posterization<float> posterization_analyzer(data.data(), N,
                                              global_dimensions.data());
          
  timer.start();
  std::vector<int> segmentation_map =
      posterization_analyzer.get_segmentation_map(threhshold);
  timer.stop("Posterization map generation ");
  timer.start();
  posterization_analyzer.evaluate();
  timer.stop("Posterization evaluation ");

  // rearrange the label from 0 to the uniuqe label count-1
  timer.start();
  std::unordered_map<int, int> label_map;
  label_map.reserve(segmentation_map.size());
  int current_label = 0;
  // for (int i = 0; i < segmentation_map.size(); i++) {
  //   if (label_map.find(segmentation_map[i]) == label_map.end()) {
  //     label_map[segmentation_map[i]] = current_label;
  //     current_label++;
  //   }
  //   segmentation_map[i] = label_map[segmentation_map[i]];
  // }
 for (int i = 0; i < segmentation_map.size(); i++) {
  auto emplace_result = label_map.emplace(segmentation_map[i], current_label);
  if (emplace_result.second) {
    current_label++;
  }
  segmentation_map[i] = emplace_result.first->second;
}
  timer.stop("Rearrange label ");

  SZ::writefile("segmentation_map.dat", segmentation_map.data(),
                segmentation_map.size());
  return 0;
}