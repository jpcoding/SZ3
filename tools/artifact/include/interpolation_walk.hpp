#include "SZ3/api/sz.hpp"
#include "SZ3/api/sz_config.hpp"
#include <array>
#include <arrray>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

template <class T, uint N> class interpolation_walk {

public:
  interpolation_walk(std::array<size_t, N> global_dimensions, T *data) {
    this->global_dimensions = global_dimensions;
    this->input_data = data;
    num_elements = 1;
    for (int i = 0; i < N; i++) {
      num_elements *= global_dimensions[i];
    }
  }

  // Artifact detection walks
  double interp_walk(int idx, int idy)
  {
    return (input_data, idx, idy);
  }

private:
  size_t num_elements;
  std::array<size_t, N> global_dimensions;
  T *input_data;

  // check index is valid
  inline bool is_valid_index(int index) {
    if (index < 0 || index >= num_elements) {
      return false;
    }
    return true;
  }

  // detection help functions
  T get_value(int idx, int idy) {
    return data[idx + idy * global_dimensions[0]];
  }

  T get_value(int idx, int idy, int idz) {
    return data[idx + idy * global_dimensions[0] +
                idz * global_dimensions[0] * global_dimensions[1]];
  }

  // get level for 1D and any other global index
  int get_level(const size_t index) {
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
  double interp_walk(T *data, int idx, int idy) {
    int global_index = cart2global(idx, idy);
    T value = get_value(idx, idy);
    int level = get_level(idx, idy);
    int stride = 1 << (level - 1);
    double error = std::numeric_limits<double>::max();
    // construct the grind of the interpolation points
    std::array<int, 4> interp_dir1 = {
        cart2global(idx - 3 * stride, idy), cart2global(idx - stride, idy),
        cart2global(idx + stride, idy), cart2global(idx + 3 * stride, idy)};
    std::array<int, 4> interp_dir2 = {
        cart2global(idx, idy - 3 * stride), cart2global(idx, idy - stride),
        cart2global(idx, idy + stride), cart2global(idx, idy + 3 * stride)};
    // Try all the possible interpolation methods
    // Deal with corner cases first
    // 1. right lorenzo
    if (is_valid_index(interp_dir1[1])) {
      error = std::min(error, std::abs(value - data[interp_dir1[1]]));
    }
    if (is_valid_index(interp_dir2[1])) {
      error = std::min(error, std::abs(value - data[interp_dir2[1]]));
    }
    // 2. right linear
    if (is_valid_index(interp_dir1[0]), is_valid_index(interp_dir1[1])) {
      error = std::min(
          error, std::abs(value - SZ::interp_linear1(data[interp_dir1[0]],
                                                     data[interp_dir1[2]])));
    }
    if (is_valid_index(interp_dir2[0]), is_valid_index(interp_dir2[1])) {
      error = std::min(
          error, std::abs(value - SZ::interp_linear1(data[interp_dir2[0]],
                                                     data[interp_dir2[2]])));
    }
    // 3. normal linear interpolation
    if (is_valid_index(interp_dir1[1]) && is_valid_index(interp_dir1[2])) {
      error = std::min(
          error, std::abs(value -
                          (data[interp_dir1[1]] + data[interp_dir1[2]]) / 2.0));
    }
    if (is_valid_index(interp_dir2[1]) && is_valid_index(interp_dir2[2])) {
      error = std::min(
          error, std::abs(value -
                          (data[interp_dir2[1]] + data[interp_dir2[2]]) / 2.0));
    }
    // 4. left quadratic interpolation
    if (is_valid_index(interp_dir1[1]) && is_valid_index(interp_dir1[2]) &&
        is_valid_index(interp_dir1[3])) {
      // left quadratic
      error = std::min(
          error, std::abs(value - SZ::interp_quad_1(data[interp_dir1[1]],
                                                    data[interp_dir1[2]],
                                                    data[interp_dir1[3]])));
    }
    if (is_valid_index(interp_dir2[1]) && is_valid_index(interp_dir2[2]) &&
        is_valid_index(interp_dir2[3])) {
      // left quadratic
      error = std::min(
          error, std::abs(value - SZ::interp_quad_1(data[interp_dir2[1]],
                                                    data[interp_dir2[2]],
                                                    data[interp_dir2[3]])));
    }
    // 5. right quadratic interpolation1
    if (is_valid_index(interp_dir1[0]) && is_valid_index(interp_dir1[1]) &&
        is_valid_index(interp_dir1[2])) {
      // left quadratic
      error = std::min(
          error, std::abs(value - SZ::interp_quad_2(data[interp_dir1[0]],
                                                    data[interp_dir1[1]],
                                                    data[interp_dir1[2]])));
    }
    if (is_valid_index(interp_dir2[0]) && is_valid_index(interp_dir2[1]) &&
        is_valid_index(interp_dir2[2])) {
      // left quadratic
      error = std::min(
          error, std::abs(value - SZ::interp_quad_2(data[interp_dir2[0]],
                                                    data[interp_dir2[1]],
                                                    data[interp_dir2[2]])));
    }
    // 6. normal quadratic interpolation and cubic interpolation
    if (is_valid_index(interp_dir1[0]) && is_valid_index(interp_dir1[1]) &&
        is_valid_index(interp_dir1[2]) && is_valid_index(interp_dir1[3])) {
      // cubic
      error = std::min(
          error, std::abs(value - SZ::interp_cubic(data[interp_dir1[0]],
                                                   data[interp_dir1[1]],
                                                   data[interp_dir1[2]],
                                                   data[interp_dir1[3]])));
    }
    if (is_valid_index(interp_dir1[0]) && is_valid_index(interp_dir1[1]) &&
        is_valid_index(interp_dir1[2]) && is_valid_index(interp_dir1[3])) {
      // cubic
      error = std::min(
          error, std::abs(value - SZ::interp_cubic(data[interp_dir1[0]],
                                                   data[interp_dir1[1]],
                                                   data[interp_dir1[2]],
                                                   data[interp_dir1[3]])));
    }

    // 7. right corner quodratic interpolation that uses -5xStride
    int corner_idx = cart2global(idx - 5 * stride, idy);
    int corner_idy = cart2global(idx, idy - 5 * stride);
    if (is_valid_index(corner_idx) && is_valid_index(interp_dir1[0]) &&
        is_valid_index(interp_dir1[1])) {
      error = std::min(
          error, std::abs(value - SZ::interp_quad_3(data[corner_idx],
                                                    data[interp_dir1[0]],
                                                    data[interp_dir1[1]])));
    }
    if (is_valid_index(corner_idy) && is_valid_index(interp_dir2[0]) &&
        is_valid_index(interp_dir2[1])) {
      {
        error = std::min(
            error, std::abs(value - SZ::interp_quad_3(data[corner_idy],
                                                      data[interp_dir2[0]],
                                                      data[interp_dir2[1]])));
      }
      return error; // return the minimum possible error
    }
  }
};