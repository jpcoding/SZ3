#include "SZ3/utils/Interpolators.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

template <class T> class InterpolationWalk {

public:
  InterpolationWalk(T *data, int N, int *global_dimensions) {
    this->input_data = data;
    this->N = N;
    this->num_elements = 1;
    this->global_dimensions.resize(N);
    for (int i = 0; i < N; i++) {
      this->num_elements *= global_dimensions[i];
      this->global_dimensions[i] = global_dimensions[i];
    }
  }

  // Deconstructor
  ~InterpolationWalk() {  }

  // Artifact detection walks

  T interp_walk(T *ptr) {
    if (N == 2) {
      return interp_walk(input_data, (ptr - input_data) % global_dimensions[0],
                         (ptr - input_data) / global_dimensions[0]);
    }
  }

  T interp_walk(int global_index) {
    if (N == 2) {
      return interp_walk(input_data, global_index % global_dimensions[0],
                         global_index / global_dimensions[0]);
    }
    else if (N==3)
    {
      int idx = global_index % global_dimensions[0];
      int idy = (global_index / global_dimensions[0]) % global_dimensions[1];
      int idz = global_index / (global_dimensions[0] * global_dimensions[1]);
      return interp_walk(input_data, idx, idy, idz);
    }
  }

  T interp_walk(int idx, int idy) { return interp_walk(input_data, idx, idy); }

private:
  int num_elements;
  std::vector<int> global_dimensions;
  T *input_data;
  int N;

  // check index is valid
  inline bool is_valid_index(int index) {
    if (index < 0 || index >= num_elements) {
      return false;
    }
    return true;
  }

  // detection help functions
  T get_value(int idx, int idy) {
    return input_data[idx + idy * global_dimensions[0]];
  }

  T get_value(int idx, int idy, int idz) {
    return input_data[idx + idy * global_dimensions[0] +
                      idz * global_dimensions[0] * global_dimensions[1]];
  }

  // get level for 1D and any other global index
  int get_level(const int index) {
    int level = 1;
    int stride = 1 << level;
    if (N == 1) {
      while (((index & (stride - 1)) == 0)) {
        level++;
        stride = 1 << level;
      }
      return level;
    }
    if (N == 2) {
      int x = index % global_dimensions[0];
      int y = index / global_dimensions[0];
      while (((x & (stride - 1)) == 0) && ((y & (stride - 1)) == 0)) {
        level++;
        stride = 1 << level;
      }
      return level;
    } else if (N == 3) {
      int x = index % global_dimensions[0];
      int y = (index / global_dimensions[0]) % global_dimensions[1];
      int z = index / (global_dimensions[0] * global_dimensions[1]);
      while (((x & (stride - 1)) == 0) && ((y & (stride - 1)) == 0) &&
             ((z & (stride - 1)) == 0)) {
        level++;
        stride = 1 << level;
      }
      return level;
    }
  }

  // get level for 2D
  int get_level(const int x, const int y) {
    int level = 1;
    int stride = 1 << level;
    while (((x & (stride - 1)) == 0) && ((y & (stride - 1)) == 0)) {
      level++;
      stride = 1 << level;
    }
    return level;
  }
  // get level for 3D
  int get_level(int x, int y, int z) {
    int level = 1;
    int stride = 1 << level;
    while (((x & (stride - 1)) == 0) && ((y & (stride - 1)) == 0) &&
           ((z & (stride - 1)) == 0)) {
      level++;
      stride = 1 << level;
    }
    return level;
  }

  int cart2global(int idx, int idy) { return idx + idy * global_dimensions[0]; }

  int cart2global(int idx, int idy, int idz) {
    return idx + idy * global_dimensions[0] +
           idz * global_dimensions[0] * global_dimensions[1];
  }

  // interp walk on 2D data
  // return the absolute minimum error of the interpolation after trying all the
  // possible interpolation methods
  T interp_walk(T *data, int idx, int idy) {
    int global_index = cart2global(idx, idy);
    T value = data[global_index];
    int level = get_level(idx, idy);
    int stride = 1 << (level - 1);
    T error = std::numeric_limits<T>::max();
    // construct the grind of the interpolation points
    std::array<int, 4> interp_dir1 = {
        cart2global(idx - 3 * stride, idy), cart2global(idx - stride, idy),
        cart2global(idx + stride, idy), cart2global(idx + 3 * stride, idy)};
    std::array<int, 4> interp_dir2 = {
        cart2global(idx, idy - 3 * stride), cart2global(idx, idy - stride),
        cart2global(idx, idy + stride), cart2global(idx, idy + 3 * stride)};

    // validate the interpolation points
    std::array<bool, 4> valid_interp_dir1 = {
        is_valid_index(interp_dir1[0]), is_valid_index(interp_dir1[1]),
        is_valid_index(interp_dir1[2]), is_valid_index(interp_dir1[3])};
    std::array<bool, 4> valid_interp_dir2 = {
        is_valid_index(interp_dir2[0]), is_valid_index(interp_dir2[1]),
        is_valid_index(interp_dir2[2]), is_valid_index(interp_dir2[3])};

    // Try all the possible interpolation methods
    // Deal with corner cases first
    // 1. right lorenzo
    // if (valid_interp_dir1[1]) {
    //   error = std::min(error, std::abs(value - data[interp_dir1[1]]));
    // }
    // if (valid_interp_dir2[1]) {
    //   error = std::min(error, std::abs(value - data[interp_dir2[1]]));
    // }
    // // 2. right linear
    // if (valid_interp_dir1[0] && valid_interp_dir1[1]) {
    //   error = std::min(
    //       error, std::abs(value - SZ::interp_linear1(data[interp_dir1[0]],
    //                                                  data[interp_dir1[1]])));
    // }
    // if (valid_interp_dir2[0] && valid_interp_dir2[1]) {
    //   error = std::min(
    //       error, std::abs(value - SZ::interp_linear1(data[interp_dir2[0]],
    //                                                  data[interp_dir2[1]])));
    // }
    // 3. regular linear interpolation
    if (valid_interp_dir1[1] && valid_interp_dir1[2]) {
      error = std::min(
          error,
          std::abs(value - (data[interp_dir1[1]] + data[interp_dir1[2]]) / 2));
    }
    if (valid_interp_dir2[1] && valid_interp_dir2[2]) {
      error = std::min(
          error,
          std::abs(value - (data[interp_dir2[1]] + data[interp_dir2[2]]) / 2));
    }
    // // 4. left quadratic interpolation
    // if (valid_interp_dir1[1] && valid_interp_dir1[2] && valid_interp_dir1[3]) {
    //   // left quadratic
    //   error = std::min(
    //       error, std::abs(value - SZ::interp_quad_1(data[interp_dir1[1]],
    //                                                 data[interp_dir1[2]],
    //                                                 data[interp_dir1[3]])));
    // }
    // if (valid_interp_dir2[1] && valid_interp_dir2[2] && valid_interp_dir2[3]) {
    //   // left quadratic
    //   error = std::min(
    //       error, std::abs(value - SZ::interp_quad_1(data[interp_dir2[1]],
    //                                                 data[interp_dir2[2]],
    //                                                 data[interp_dir2[3]])));
    // }
    // // 5. right quadratic interpolation1
    // if (valid_interp_dir1[0] && valid_interp_dir1[1] && valid_interp_dir1[2]) {
    //   // left quadratic
    //   error = std::min(
    //       error, std::abs(value - SZ::interp_quad_2(data[interp_dir1[0]],
    //                                                 data[interp_dir1[1]],
    //                                                 data[interp_dir1[2]])));
    // }
    // if (valid_interp_dir2[0] && valid_interp_dir2[1] && valid_interp_dir2[2]) {
    //   // left quadratic
    //   error = std::min(
    //       error, std::abs(value - SZ::interp_quad_2(data[interp_dir2[0]],
    //                                                 data[interp_dir2[1]],
    //                                                 data[interp_dir2[2]])));
    // }
    // 6. regular quadratic interpolation and cubic interpolation
    if (valid_interp_dir1[0] && valid_interp_dir1[1] && valid_interp_dir1[2] &&
        valid_interp_dir1[3]) {
      // cubic
      error = std::min(
          error, std::abs(value - SZ::interp_cubic(data[interp_dir1[0]],
                                                   data[interp_dir1[1]],
                                                   data[interp_dir1[2]],
                                                   data[interp_dir1[3]])));
    }
    if (valid_interp_dir2[0] && valid_interp_dir2[1] && valid_interp_dir2[2] &&
        valid_interp_dir2[3]) {
      // cubic
      error = std::min(
          error, std::abs(value - SZ::interp_cubic(data[interp_dir2[0]],
                                                   data[interp_dir2[1]],
                                                   data[interp_dir2[2]],
                                                   data[interp_dir2[3]])));
    }

    // // 7. right corner quodratic interpolation that uses -5xStride
    // int corner_idx = cart2global(idx - 5 * stride, idy);
    // int corner_idy = cart2global(idx, idy - 5 * stride);
    // if (is_valid_index(corner_idx) && valid_interp_dir1[0] &&
    //     valid_interp_dir1[1]) {
    //   error = std::min(
    //       error, std::abs(value - SZ::interp_quad_3(data[corner_idx],
    //                                                 data[interp_dir1[0]],
    //                                                 data[interp_dir1[1]])));
    // }
    // if (is_valid_index(corner_idy) && valid_interp_dir2[0] &&
    //     valid_interp_dir2[1]) {
    //   error = std::min(
    //       error, std::abs(value - SZ::interp_quad_3(data[corner_idy],
    //                                                 data[interp_dir2[0]],
    //                                                 data[interp_dir2[1]])));
    // }
    return error; // return the minimum possible error
  }

  T interp_walk(T* data, int idx, int idy, int idz)
  {
    int global_idx = cart2global(idx, idy, idz);
    T value = data[global_idx];
    T error = std::numeric_limits<T>::max();
    int level = get_level(idx, idy, idz);
    int stride = 1 << level;
    int stride3x = 3 * stride;

    // Try all the possible interpolation methods
    // Get all neighbors
    std::array<int, 4> interp_dir1 = {cart2global(idx -stride3x , idy, idz),
                                      cart2global(idx -stride , idy, idz),
                                      cart2global(idx +stride , idy, idz),
                                      cart2global(idx +stride3x , idy, idz)};
    std::array<int, 4> interp_dir2 = {cart2global(idx, idy -stride3x, idz),
                                      cart2global(idx, idy -stride, idz),
                                      cart2global(idx, idy +stride, idz),
                                      cart2global(idx, idy +stride3x, idz)};
    std::array<int, 4> interp_dir3 = {cart2global(idx, idy, idz -stride3x),
                                      cart2global(idx, idy, idz -stride),
                                      cart2global(idx, idy, idz +stride),
                                      cart2global(idx, idy, idz +stride3x)};      
    std::array<bool, 4> valid_interp_dir1 = {is_valid_index(interp_dir1[0]),
                                             is_valid_index(interp_dir1[1]),
                                             is_valid_index(interp_dir1[2]),
                                             is_valid_index(interp_dir1[3])};
    std::array<bool, 4> valid_interp_dir2 = {is_valid_index(interp_dir2[0]),
                                              is_valid_index(interp_dir2[1]),
                                              is_valid_index(interp_dir2[2]),
                                              is_valid_index(interp_dir2[3])};
    std::array<bool, 4> valid_interp_dir3 = {is_valid_index(interp_dir3[0]),
                                              is_valid_index(interp_dir3[1]),
                                              is_valid_index(interp_dir3[2]),
                                              is_valid_index(interp_dir3[3])};  
    // 1. check linear interpolation 
    if (valid_interp_dir1[1] && valid_interp_dir1[2]) {
      // left linear
      error = std::min(
          error, std::abs(value - SZ::interp_linear(data[interp_dir1[1]],
                                                     data[interp_dir1[2]])));
    }
    if (valid_interp_dir2[1] && valid_interp_dir2[2]) {
      // left linear
      error = std::min(
          error, std::abs(value - SZ::interp_linear(data[interp_dir2[1]],
                                                     data[interp_dir2[2]])));
    }
    if (valid_interp_dir3[1] && valid_interp_dir3[2]) {
      // left linear
      error = std::min(
          error, std::abs(value - SZ::interp_linear(data[interp_dir3[1]],
                                                     data[interp_dir3[2]])));
    }

    // 2. check cubic interpolation
    if (valid_interp_dir1[0] && valid_interp_dir1[1] && valid_interp_dir1[2] &&
        valid_interp_dir1[3]) {
      // cubic
      error = std::min(
          error, std::abs(value - SZ::interp_cubic(data[interp_dir1[0]],
                                                   data[interp_dir1[1]],
                                                   data[interp_dir1[2]],
                                                   data[interp_dir1[3]])));
    }
    if (valid_interp_dir2[0] && valid_interp_dir2[1] && valid_interp_dir2[2] &&
        valid_interp_dir2[3]) {
      // cubic
      error = std::min(
          error, std::abs(value - SZ::interp_cubic(data[interp_dir2[0]],
                                                   data[interp_dir2[1]],
                                                   data[interp_dir2[2]],
                                                   data[interp_dir2[3]])));
    }
    if (valid_interp_dir3[0] && valid_interp_dir3[1] && valid_interp_dir3[2] &&
        valid_interp_dir3[3]) {
      // cubic
      error = std::min(
          error, std::abs(value - SZ::interp_cubic(data[interp_dir3[0]],
                                                   data[interp_dir3[1]],
                                                   data[interp_dir3[2]],
                                                   data[interp_dir3[3]])));
    }
    return error;
  }
};