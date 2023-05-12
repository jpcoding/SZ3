#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "block_detection.hpp"
#include <array>
#include <iostream>
#include <map>
#include <vector>

void majority_vote(int detect_block_size_x, double detect_ratio_x,
                   int detect_block_size_y, double detect_ratio_y,
                   int detect_block_size_z, double detect_ratio_z) {
  int block_sizes[3] = {detect_block_size_x, detect_block_size_y,
                        detect_block_size_z};
  double detect_ratios[3] = {detect_ratio_x, detect_ratio_y, detect_ratio_z};
  int max_count = 1;
  int max_index = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = i; j < 3; j++) {
      max_count = 1;
      if (block_sizes[i] == block_sizes[j]) {
        max_count++;
        max_index = i;
      }
    }
  }
  double max_ratio = 0.0;
  for (int i = 0; i < 3; i++) {
    max_ratio = detect_ratios[i];
    for (int j = i; j < 3; j++) {
      if (detect_ratios[i] > detect_ratios[j]) {
        max_ratio = detect_ratios[i];
        max_index = i;
      }
    }
  }
  if (max_count > 1) {
    std::cout << "majority vote: " << block_sizes[max_index] << std::endl;
  } else {
    std::cout << "no majority vote, use the largest block size: "
              << block_sizes[0] << std::endl;
    std::cout << "majority vote: " << block_sizes[max_index] << std::endl;
  }
}

int main(int argc, char **argv) {
  SZ::Timer timer;

  if (argc < 2) {
    std::cout
        << "Usage: ./test_block_detection <input_file> <N> <dim1 dim2 ...> "
        << std::endl;
    return 0;
  }

  int N = atoi(argv[2]);
  int num_elements = 1;
  std::vector<int> global_dimensions;
  for (int i = 0; i < N; i++) {
    global_dimensions.push_back(atoi(argv[3 + i]));
    num_elements *= global_dimensions[i];
    std::cout << "dim " << i << " = " << global_dimensions[i] << std::endl;
  }

  std::vector<float> input_data(num_elements);
  SZ::readfile<float>(argv[1], num_elements, input_data.data());

  float data_min, dat_max;
  // float range = normalization(input_data, data_min, dat_max);
  // block detection
  float threhshold = 1e-5;
  if (argc == (3 + N + 1)) {
    threhshold = atof(argv[argc - 1]);
  }

  if (argc == (3 + N + 2)) {
    threhshold = atof(argv[argc - 1]);
  }
  std::cout << "input threhshold: " << threhshold << std::endl;
  // threhshold = threhshold / range;
  // std::cout << "factoirzed threshold: " << threhshold << std::endl;
  int detect_blocksize;
  int min_block_size = 4;
  int max_block_size = 16;
  if (N == 2) {
    // if input data is 2d, sample on two directions, sample 10%
    // 1. detect along y direction
    double sample_ratio = 0.2;
    int y_sample_stride = 3;
    int num_ysample = global_dimensions[1] /y_sample_stride;
    int y_sample_begin = 1;
    int y_sample_end = num_ysample;
    double detect_ratio_y = 0;
    int detect_block_size_y = 0;
    timer.start();

    for (detect_blocksize = min_block_size; detect_blocksize < max_block_size;
         detect_blocksize++) {
      int current_block_size = detect_blocksize;
      int current_effective_block_count = 0;
      int current_artifact_block_count = 0;
      double current_detect_ratio = 0;
      std::cout << "current block size: " << current_block_size << std::endl;
      for (int i = y_sample_begin; i < y_sample_end; i++) {
        int begin = i * y_sample_stride * global_dimensions[0];
        std::cout<<"begin = "<<begin<<std::endl;
        int sample_size = global_dimensions[0];
        int stride = 1;
        int i_effective_block_count = 0;
        int i_artifact_block_count = 0;
        auto block_detection = BlockDetection<float>(input_data.data(), begin,
                                                     sample_size, stride);
        block_detection.block_detection(threhshold, detect_blocksize,
                                        i_effective_block_count,
                                        i_artifact_block_count);                               
        std::cout<<"i_effective_block_count = "<<i_effective_block_count<<std::endl;
        std::cout<<"i_artifact_block_count = "<<i_artifact_block_count<<std::endl;
        current_effective_block_count += i_effective_block_count;
        current_artifact_block_count += i_artifact_block_count;

        std::cout<<"current_effective_block_count = "<<current_effective_block_count<<std::endl;
        std::cout<<"current_artifact_block_count = "<<current_artifact_block_count<<std::endl;
      }
      current_detect_ratio = (1.0*current_artifact_block_count)/(1.0*current_effective_block_count);
      std::cout<<"current_effective_block_count = "<<current_effective_block_count<<std::endl;
      std::cout<<"current_artifact_block_count = "<<current_artifact_block_count<<std::endl;
      std::cout << "block size = "<< current_block_size  << " current detect ratio: " << current_detect_ratio <<std::endl;
      if (current_detect_ratio > detect_ratio_y) {
        detect_ratio_y = current_detect_ratio;
        detect_block_size_y = current_block_size;
      }
    }
    timer.stop("block detection");
    std::cout << "block size = " << detect_block_size_y << std::endl;
    std::cout << "detect ratio = " << detect_ratio_y << std::endl;

  // 2. detect along x direction
    int x_sample_stride = 4;
    int num_xsample = global_dimensions[1] /x_sample_stride;
    int x_sample_begin = 1;
    int x_sample_end =num_xsample;
    double detect_ratio_x= 0;
    int detect_block_size_x = 0;

    for (detect_blocksize = min_block_size; detect_blocksize < max_block_size;
         detect_blocksize++) {
      int current_block_size = detect_blocksize;
      int current_effective_block_count = 0;
      int current_artifact_block_count = 0;
      double current_detect_ratio = 0;
      for (int i = x_sample_begin; i < x_sample_end; i++) {
        int begin = i * x_sample_stride;
        int sample_size = global_dimensions[1];
        int stride = global_dimensions[0];
        int i_effective_block_count = 0;
        int i_artifact_block_count = 0;
        auto block_detection = BlockDetection<float>(input_data.data(), begin,
                                                     sample_size, stride);
        block_detection.block_detection(threhshold, detect_blocksize,
                                        i_effective_block_count,
                                        i_artifact_block_count);
        current_effective_block_count += i_effective_block_count;
        current_artifact_block_count += i_artifact_block_count;
      }
      current_detect_ratio = (double)current_artifact_block_count /
                             (double)(current_effective_block_count);
      if (current_detect_ratio > detect_ratio_x) {
        detect_ratio_x = current_detect_ratio;
        detect_block_size_x = current_block_size;
      }
    }

    timer.stop("block detection");
    std::cout << "block size = " << detect_block_size_x << std::endl;
    std::cout << "detect ratio = " << detect_ratio_x << std::endl;

    std::cout << "block size = " << detect_block_size_y << std::endl;
    std::cout << "detect ratio = " << detect_ratio_y << std::endl;

    // decision on 2D data
    if (detect_ratio_x > detect_ratio_y) {
      std::cout << "block size = " << detect_block_size_x << std::endl;
      std::cout << "detect ratio = " << detect_ratio_x << std::endl;
    } else {
      std::cout << "block size = " << detect_block_size_y << std::endl;
      std::cout << "detect ratio = " << detect_ratio_y << std::endl;
    }
  } else if (N == 3) {
    // block detection
    // if input data is 3d, sample on three directions, sample 10%
    int num_xsample = global_dimensions[0] / 20;
    int num_ysample = global_dimensions[1] / 20;
    int num_zsample = global_dimensions[2] / 20;

    int x_stride = global_dimensions[0] / num_xsample;
    int y_stride = global_dimensions[1] / num_ysample;
    int z_stride = global_dimensions[2] / num_zsample;

    double block_detect_threshold = 1e-5;

    // 1. detect along z direction
    int detect_block_size_z = 0;
    double detect_ratio_z = 0;
    for (int i = 1; i < num_xsample - 1; i++) {
      for (int j = 1; j < num_ysample - 1; j++) {
        int begin = i * x_stride + j * y_stride * global_dimensions[0];
        int sample_size = global_dimensions[2];
        int stride = global_dimensions[0] * global_dimensions[1];
        auto block_detection = BlockDetection<float>(input_data.data(), begin,
                                                     sample_size, stride);
        block_detection.set_flush_threshold(block_detect_threshold);
        double i_detect_ratio = 0;
        auto iblock_size = block_detection.block_detection(
            block_detect_threshold, i_detect_ratio);
        if (i_detect_ratio > detect_ratio_z) {
          detect_ratio_z = i_detect_ratio;
          detect_block_size_z = iblock_size;
        }
      }
    }
    std::cout << "z block size = " << detect_block_size_z << std::endl;
    std::cout << "z detect ratio = " << detect_ratio_z << std::endl;

    // 2. detect along y direction
    int detect_block_size_y = 0;
    double detect_ratio_y = 0;
    for (int i = 1; i < num_xsample - 1; i++) {
      for (int j = 1; j < num_zsample - 1; j++) {
        int begin = i * x_stride +
                    j * z_stride * global_dimensions[0] * global_dimensions[1];
        int sample_size = global_dimensions[1];
        int stride = global_dimensions[0];
        auto block_detection = BlockDetection<float>(input_data.data(), begin,
                                                     sample_size, stride);
        block_detection.set_flush_threshold(block_detect_threshold);
        double i_detect_ratio;
        auto iblock_size = block_detection.block_detection(
            block_detect_threshold, i_detect_ratio);
        if (i_detect_ratio > detect_ratio_y) {
          detect_ratio_y = i_detect_ratio;
          detect_block_size_y = iblock_size;
        }
      }
    }
    std::cout << "y block size = " << detect_block_size_y << std::endl;
    std::cout << "y detect ratio = " << detect_ratio_y << std::endl;

    // 3. detect along x direction
    int detect_block_size_x = 0;
    double detect_ratio_x = 0;
    for (int i = 1; i < num_ysample - 1; i++) {
      for (int j = 1; j < num_zsample - 1; j++) {
        int begin = i * y_stride * global_dimensions[0] +
                    j * z_stride * global_dimensions[1] * global_dimensions[0];
        int sample_size = global_dimensions[0];
        int stride = 1;
        auto block_detection = BlockDetection<float>(input_data.data(), begin,
                                                     sample_size, stride);
        block_detection.set_flush_threshold(block_detect_threshold);
        double i_detect_ratio;
        auto iblock_size = block_detection.block_detection(
            block_detect_threshold, i_detect_ratio);
        if (i_detect_ratio > detect_ratio_x) {
          detect_ratio_x = i_detect_ratio;
          detect_block_size_x = iblock_size;
        }
      }
    }
    std::cout << "block size x= " << detect_block_size_x << std::endl;
    std::cout << "detect ratio x = " << detect_ratio_x << std::endl;

    std::cout << "block size y = " << detect_block_size_y << std::endl;
    std::cout << "detect ratio y = " << detect_ratio_y << std::endl;

    std::cout << "block size z = " << detect_block_size_z << std::endl;
    std::cout << "detect ratio z = " << detect_ratio_z << std::endl;

    majority_vote(detect_block_size_x, detect_ratio_x, detect_block_size_y,
                  detect_ratio_y, detect_block_size_z, detect_ratio_z);

    return 0;
  }
}