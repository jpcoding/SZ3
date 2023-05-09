#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "block_detection.hpp"
#include <array>
#include <iostream>
#include <vector>

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
  float range = normalization( input_data, data_min, dat_max);
  // block detection

  if (N == 2) {
    // if input data is 2d, sample on two directions, sample 10%
    int num_xsample = global_dimensions[0] / 10;
    int num_ysample = global_dimensions[1] / 10;

    // 1. detect along x direction
    int y_stride = global_dimensions[1] / num_ysample;
    double detect_ratio_y = 0;
    double detect_block_size_y = 0;
    timer.start();
    for (int i = 0; i < num_ysample; i++) {
      int begin = i * y_stride * global_dimensions[0];
      int sample_size = global_dimensions[0];
      int stride = 1;
      auto block_detection =
          BlockDetection<float>(input_data.data(), begin, sample_size, stride);

      double i_detect_ratio;
      auto iblock_size = block_detection.block_detection(1e-5, i_detect_ratio);
      if (i_detect_ratio > detect_ratio_y) {
        detect_ratio_y = i_detect_ratio;
        detect_block_size_y = iblock_size;
      }
    }
    timer.stop("block detection");
    std::cout << "block size = " << detect_block_size_y << std::endl;
    std::cout << "detect ratio = " << detect_ratio_y << std::endl;

    // 1. detect along x direction
    int x_stride = global_dimensions[0] / num_xsample;
    double detect_ratio_x = 0;
    double detect_block_size_x = 0;
    timer.start();
    for (int i = 0; i < num_xsample; i++) {
      int begin = i * x_stride;
      int sample_size = global_dimensions[1];
      int stride = global_dimensions[0];
      auto block_detection =
          BlockDetection<float>(input_data.data(), begin, sample_size, stride);

      double i_detect_ratio;
      auto iblock_size = block_detection.block_detection(1e-5, i_detect_ratio);
      if (i_detect_ratio > detect_ratio_x) {
        detect_ratio_x = i_detect_ratio;
        detect_block_size_x = iblock_size;
      }
    }
    timer.stop("block detection");
    std::cout << "block size = " << detect_block_size_x << std::endl;
    std::cout << "detect ratio = " << detect_ratio_x << std::endl;

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
    int num_xsample = global_dimensions[0] / 10;
    int num_ysample = global_dimensions[1] / 10;
    int num_zsample = global_dimensions[2] / 10;

    int x_stride = global_dimensions[0] / num_xsample;
    int y_stride = global_dimensions[1] / num_ysample;
    int z_stride = global_dimensions[2] / num_zsample;

    double block_detect_threshold = 1e-5;

    // 1. detect along z direction
    int detect_block_size_z = 0;
    double detect_ratio_z = 0;
    for (int i = 0; i < num_xsample; i++) {
      for (int j = 0; j < num_ysample; j++) {
        int begin = i * x_stride + j * y_stride * global_dimensions[0];
        int sample_size = global_dimensions[2];
        int stride = global_dimensions[0] * global_dimensions[1];
        auto block_detection = BlockDetection<float>(input_data.data(), begin,
                                                     sample_size, stride);
        double i_detect_ratio=0;
        auto iblock_size =
            block_detection.block_detection(block_detect_threshold, i_detect_ratio);
        if (i_detect_ratio > detect_ratio_z) {
          detect_ratio_z = i_detect_ratio;
          detect_block_size_z = iblock_size;
        }
      }
    }
    std::cout << "block size = " << detect_block_size_z << std::endl;
    std::cout << "detect ratio = " << detect_ratio_z << std::endl;

    // 2. detect along y direction
    int detect_block_size_y = 0;
    double detect_ratio_y = 0;
    for (int i = 0; i < num_xsample; i++) {
      for (int j = 0; j < num_zsample; j++) {
        int begin = i * x_stride +
                    j * z_stride * global_dimensions[0] * global_dimensions[1];
        int sample_size = global_dimensions[1];
        int stride = global_dimensions[0];
        auto block_detection = BlockDetection<float>(input_data.data(), begin,
                                                     sample_size, stride);
        double i_detect_ratio;
        auto iblock_size =
            block_detection.block_detection(block_detect_threshold, i_detect_ratio);
        if (i_detect_ratio > detect_ratio_y) {
          detect_ratio_y = i_detect_ratio;
          detect_block_size_y = iblock_size;
        }
      }
    }
    std::cout << "block size = " << detect_block_size_y << std::endl;
    std::cout << "detect ratio = " << detect_ratio_y << std::endl;

    // 3. detect along x direction
    int detect_block_size_x = 0;
    double detect_ratio_x = 0;
    for (int i = 0; i < num_ysample; i++) {
      for (int j = 0; j < num_zsample; j++) {
        int begin = i * y_stride * global_dimensions[0] +
                    j * z_stride * global_dimensions[1] * global_dimensions[0];
        int sample_size = global_dimensions[0];
        int stride = 1;
        auto block_detection = BlockDetection<float>(input_data.data(), begin,
                                                     sample_size, stride);
        double i_detect_ratio;
        auto iblock_size =
            block_detection.block_detection(block_detect_threshold, i_detect_ratio);
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


    // decision on 3D data
    if (detect_ratio_x > detect_ratio_y) {
      if (detect_ratio_x > detect_ratio_z) {
        std::cout << "block size = " << detect_block_size_x << std::endl;
        std::cout << "detect ratio = " << detect_ratio_x << std::endl;
      } else {
        std::cout << "block size = " << detect_block_size_z << std::endl;
        std::cout << "detect ratio = " << detect_ratio_z << std::endl;
      }
    } 
    else {
      if (detect_ratio_y > detect_ratio_z) {
        std::cout << "block size = " << detect_block_size_y << std::endl;
        std::cout << "detect ratio = " << detect_ratio_y << std::endl;
      } else {
        std::cout << "block size = " << detect_block_size_z << std::endl;
        std::cout << "detect ratio = " << detect_ratio_z << std::endl;
      }
    }
  }

  // int begin = 1600*global_dimensions[0];
  // int end = global_dimensions[0];
  // int stride = 1;
  // auto block_detection = BlockDetection<float>(input_data.data(), begin,
  // end, stride); timer.start(); auto block_size =
  // block_detection.block_detection(1e-5); timer.stop("block detection");
  // std::cout << "block size" << block_size << std::endl;

  // detect on 3d data

  return 0;
}