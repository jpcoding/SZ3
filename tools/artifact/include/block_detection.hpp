#include "SZ3/utils/FileUtil.hpp"
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>
#include <data_preprocess.hpp>

// only works on the sampled 1D data from ND data

template <class T> class BlockDetection {
public:
  BlockDetection(T *data, int begin, int sample_size, int stride) {
    this->input_data = data;
    this->sample_size = sample_size;
    this->begin = begin;
    this->input_stride = stride;
    // central derivative is used for block detection
    calculate_central_derivative();
    // for (int i = 0; i < sample_size; i++) {
    //   std::cout<<derivative[i] <<"  ";
    //   if (i % 10 == 0) {
    //     std::cout << std::endl;
    //   }
    // }
    // exit(0);
  }

  ~BlockDetection() {}

  // Called by the user
  // given a threshold, return the block detection result.
  // By default, the block size is from 4-16 (inclusive)
  // The block size is the number of elements in the block
  // The threshold is the threshold to determine if a block is an artifact
  // if the ratio of artifact blocks to effective blocks is larger than the
  // previous one, then the block size is updated

  int block_detection(T threshold, double &artifact_ratio ) {

    int effective_block_count = 0;
    int artifact_block_count = 0;
    int block_size = 0;
    int max_block_size = std::min(16, sample_size);
    artifact_ratio = 0.0;
    for (int i = 4; i <= max_block_size; i++) {
      double i_ratio = 0.0;
      int i_effective_block_count = 0;
      int i_artifact_block_count = 0;
      try_block_detection(derivative.data(), i, threshold,
                          i_effective_block_count, i_artifact_block_count);
      i_ratio =
          (double)i_artifact_block_count / (double)i_effective_block_count;

      if (i_ratio > artifact_ratio) {
        artifact_ratio = i_ratio;
        effective_block_count = i_effective_block_count;
        artifact_block_count = i_artifact_block_count;
        block_size = i;
        std::cout << "block size = " << i << "\nratio = " << i_ratio << std::endl;
      std::cout << "effective_block_count = " << i_effective_block_count
                << "\nartifact_block_count = " << i_artifact_block_count
                << std::endl;
      }
    }

    return block_size;
  }

  // int block_detection(int block_size, T threshold) {
  //   int effective_block_count = 0;
  //   int artifact_block_count = 0;
  //   try_block_detection(input_data, block_size, threshold,
  //                       effective_block_count, artifact_block_count);
  //   return artifact_block_count;
  // }

private:
  T *input_data;                 // memory address of global data
  int begin;                     // global begin index of the input 1D slice
  int sample_size;              // number of elements in the input 1D slice
  int input_stride;              // stride of the input 1D slice
  std::vector<T> derivative;     // derivative of the input data 1D slice
  double flush_threshold = 1e-7; // any value below this will be treated as 0

  // Input:
  // data: detivative of a 1D array of the original data.
  // block_size: the size of the block.
  // threshold: the threshold to determine if a block is an artifact.
  // Output:
  // effective_block_count: the number of blocks that are above the flush.
  // artifact_block_count: the number of blocks that are determined as artifact
  // blocks.
  //  Make sure you are use the correct input data!!! (derivative).
  void try_block_detection(T *data, int block_size, T threshold,
                           int &effective_block_count,
                           int &artifact_block_count) {
    int n_block = sample_size / block_size;
    effective_block_count = 0;
    artifact_block_count = 0;
    for (int i = 0; i < n_block; i++) {
      int block_begin = i * block_size;
      int block_end = block_begin + block_size;
      T block_max = block_abs_max(data, block_begin, block_end);
      if (block_max > flush_threshold) {
        effective_block_count++;
        // std::cout<< "block data " << std::abs(data[block_begin]) << "  " << std::abs(data[block_end-1])<<std::endl;
        if ((std::abs(data[block_begin]) > threshold) ||
            (std::abs(data[block_end-1]) > threshold)) {
          
          auto detect_indicator = block_range(data, block_begin+1, block_end-1);
          // std::cout << i <<"detect_indicator = " << detect_indicator << std::endl;
          // T block_std = block_std(data, block_begin, block_end);
          // std::cout << "indicator = " << detect_indicator << " - " <<  threshold << " = "  << (detect_indicator - threshold) << std::endl;
          if (detect_indicator < threshold) {
            artifact_block_count++;
          }
        }
      }
    }
  }

  // derivative calculation
  void calculate_central_derivative() {
    derivative.resize(sample_size);
    derivative[0]=0;
    derivative[sample_size-1]=0;
    for (int i = 1; i < sample_size - 1; i ++) {
      derivative[i] =
          (input_data[begin + (i+1) * input_stride] - input_data[begin + (i-1) * input_stride]) / 2;
    }
    // if(begin==0 && sample_size==100)
    // {
    //   SZ::writefile("derivative.dat", derivative.data(), derivative.size());
    //   for (int i = 0; i < sample_size; i ++) {
    //     std::cout << input_data[begin+i*input_stride]<<std::endl;
    //   }
    // }
  }

  void calculate_forward_derivate() {
    derivative.assign(sample_size, 0);
    for (int i = 0; i < sample_size - 1; i ++) {
      derivative[i] = (input_data[begin + (i + 1)*input_stride] - input_data[begin + i*input_stride]);
    }
  }

  void calculate_backward_derivative() {
    derivative.assign(sample_size, 0);
    for (int i = 1; i < sample_size; i ++) {
      derivative[i] = (input_data[begin + i*input_stride] - input_data[begin + (i - 1)*input_stride]);
    }
  }

  T block_abs_max(T *data, int block_begin, int block_end, int stride = 1) {
    T max = 0;
    for (int i = block_begin; i < block_end; i += stride) {
      max = std::max(max, std::abs(data[i]));
    }
    return max;
  }

  // candidates for block detrcer function
  T block_range(T *data, int block_begin, int block_end, int stride = 1) {
    T max = -std::numeric_limits<T>::max();
    T min = std::numeric_limits<T>::max(); 
    for (int i = block_begin; i < block_end; i += stride) {
      max = std::max(max, data[i]);
      min = std::min(min, data[i]);
    }
    // std::cout << "max = " << max << " min = " << min << std::endl;
    return (max - min);
  }

  double block_std(T *data, int block_begin, int block_end, int stride = 1) {
    double mean = 0;
    double std = 0;
    double count = 0;
    for (int i = block_begin; i < block_end; i += stride) {
      mean += data[i];
      count++;
    }
    mean /= count;
    for (int i = block_begin; i < block_end; i += stride) {
      std += (data[i] - mean) * (data[i] - mean);
    }
    std /= count;
    return std;
  }

  double block_mean(T *data, int block_begin, int block_end, int stride = 1) {
    double mean = 0;
    int count = 0;
    for (int i = block_begin; i < block_end; i += stride) {
      mean += data[i];
      count++;
    }
    mean /= count;
    return mean;
  }
};
