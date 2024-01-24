//
// Created by Kai Zhao on 4/20/20.
//
#include <algorithm>
#include <cstddef>
#include <vector>
#include <iostream>
#include <cmath>
namespace SZ {

template <class T, SZ::uint N> class ExtractionLevel {

public:
  ExtractionLevel(const size_t *original_dims, int extraction_stride,
                  const T *input_data) {

    this->stride = extraction_stride;
    this->data = input_data;
    this->original_dims.resize(N);
    this->extracted_dims.resize(N);
    for (int i = 0; i < N; i++) {
      this->original_dims[i] = original_dims[i]; // the last one is the fastest dimension
      this->extracted_dims[i] = original_dims[i] / extraction_stride;
      this->num_elements *= original_dims[i];
      this->extracted_size *= extracted_dims[i];
    }
  }

  ExtractionLevel() {}

  ~ExtractionLevel() {}

  inline std::vector<size_t> get_coordinates(size_t index) {
    std::vector<size_t> coordinates(N, 0);
    if (N == 1) {
      coordinates[0] = index ;
    } else if (N == 2) {
      coordinates[0] = (index % extracted_dims[1]);
      coordinates[1] = (index / extracted_dims[1]);
    } else if (N == 3) {
      coordinates[0] = (index % extracted_dims[2]);
      coordinates[1] =
          ((index / extracted_dims[2]) % extracted_dims[1]);
      coordinates[2] =
          (index / (extracted_dims[1] * extracted_dims[2]));

    } else {
      std::cout << "N is not supported in get_coordinates" << std::endl;
    }
    return coordinates;
  }

  inline size_t sampleIdx_to_originalIdx(size_t index) {
    std::vector<size_t> coordinates = get_coordinates(index);
    size_t original_index = 0;
    if (N == 1) {
      original_index = coordinates[0]*stride;
    } else if (N == 2) {
      original_index = coordinates[0] * original_dims[1] *stride + coordinates[1] * stride;
    } else if (N == 3) {
      original_index = coordinates[0]*stride * original_dims[1] * original_dims[2] +
                       coordinates[1]*stride * original_dims[2] + coordinates[2]*stride;
    } else {
      std::cout << "N is not supported in sampleIdx_to_originalIdx"
                << std::endl;
    }
    return original_index;
  }

  std::vector<T> get_extracted_data() {
    extracted_data.resize(extracted_size);
    for (size_t i = 0; i < extracted_size; i++) {
      extracted_data[i] = data[sampleIdx_to_originalIdx(i)];
    }
    return extracted_data;
  }

  std::vector<size_t> get_extracted_dims() { return extracted_dims; }

  int get_highest_level()
  {
    size_t max_dim = *std::max_element(extracted_dims.begin(), extracted_dims.end());
    return (int) ceil(log2(max_dim));
  }

private:
  std::vector<size_t> original_dims;
  std::vector<size_t> extracted_dims;
  size_t num_elements = 1;
  size_t extracted_size = 1;
  int stride;
  int block_size;
  const T *data;
  std::vector<T> extracted_data; // the extracted data
};
} // namespace SZ

