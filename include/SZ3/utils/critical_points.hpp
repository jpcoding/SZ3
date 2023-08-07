#include "interpolation_walk.hpp"
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace SZ {
template <class T> class CriticalPointsCalculator {

public:
  CriticalPointsCalculator(T *inputdata, int N, unsigned long *global_dimensions) {
    this->inputdata = inputdata;
    this->N = N;
    this->num_elements = 1;
    this->global_dimensions.resize(N);

    for (int i = 0; i < N; i++) {
      this->global_dimensions[i] = global_dimensions[i];
      this->num_elements *= global_dimensions[i];
    }
    this->interps = new InterpolationWalk<T>(inputdata, N, global_dimensions);
  }

  // Deconstructor
  ~CriticalPointsCalculator() { delete interps; }

  // calculate critical points map and return it
  std::vector<int> get_critical_points_map() {
    calculate_critial_points_map();
    return critical_points_map;
  }

  // match pattern

  bool error_check(const T *odata, const T *ddata, const int index,
                   const int xpadding, const int ypadding) {
    T center_error = std::abs(odata[index] - ddata[index]);
    T current_error = 0;
    int idx = index % global_dimensions[0];
    int idy = index / global_dimensions[0];

    // int idx = index;
    for (int i = -xpadding; i <= xpadding; i++) {
      // idx = index + i;
      for (int j = -ypadding; j <= ypadding; j++) {
        if (i == 0 && j == 0) {
          continue;
        }
        // idx += j * global_dimensions[0];
        // int idx = index + i + j * global_dimensions[0];
        int current_index = (idx + i) + (idy + j) * global_dimensions[0];
        current_error = std::abs(odata[current_index] - ddata[current_index]);
        if (current_error >= center_error) {
          return false;
        }
      }
    }
    return true;
  }

  bool error_check2(const T *odata, const T *ddata, const int index,
                    const int xpadding, const int ypadding) {
    T center_error = std::abs(odata[index] - ddata[index]);
    T current_error = 0;
    // check x direction
    for (int i = -xpadding; i <= xpadding; i++) {
      int left = std::abs(odata[index + i] - ddata[index + i]);
      int right = std::abs(odata[index - i] - ddata[index - i]);
      if (left > center_error || right > center_error) {
        return false;
      }
    }
    // check y direction
    for (int i = -ypadding; i <= ypadding; i++) {
      int up = std::abs(odata[index + i * global_dimensions[0]] -
                        ddata[index + i * global_dimensions[0]]);
      int down = std::abs(odata[index - i * global_dimensions[0]] -
                          ddata[index - i * global_dimensions[0]]);
      if (up > center_error || down > center_error) {
        return false;
      }
    }
    return true;
  }

  bool error_check(const T *odata, const T *ddata, const int index,
                   const int xpadding, const int ypadding, const int zpadding) {
    T center_error = std::abs(odata[index] - ddata[index]);
    T current_error = 0;
    for (int i = -xpadding; i <= xpadding; i++) {
      for (int j = -ypadding; j <= ypadding; j++) {
        for (int k = -zpadding; k <= zpadding; k++) {
          if (i == 0 && j == 0 && k == 0) {
            continue;
          }
          int idx = index + i + j * global_dimensions[0] +
                    k * global_dimensions[0] * global_dimensions[1];
          current_error = std::abs(odata[idx] - ddata[idx]);
          if (current_error >= center_error) {
            return false;
          }
        }
      }
    }
    return true;
  }

  bool cp_check(const T *data, const int index, const int xpadding,
                const int ypadding) {
    T center_value = std::abs(data[index]);
    T current_value = 0;
    for (int i = -xpadding; i <= xpadding; i++) {
      for (int j = -ypadding; j <= ypadding; j++) {
        if (i == 0 && j == 0) {
          continue;
        }
        int idx = index + i + j * global_dimensions[0];
        current_value = std::abs(data[idx]);
        if (current_value >= center_value) {
          return false;
        }
      }
    }
    return true;
  }

  // This is used to check if the point is maximum or minimum on the given data
  // paddings are gurannteed to be valid (i.e. not out of bound)
  // 2d version
  bool padding_check(const T *data, const int index, const int sign,
                     const int xpadding, const int ypadding, T tol = 1e-5) {

    // if sign is positive, then we are looking for maximum
    if (sign > 0) {
      T center_value = data[index] - tol;
      // check x direction
      for (int i = -xpadding; i <= xpadding; i++) {
        int left = data[index + i];
        int right = data[index - i];
        if (left > center_value || right > center_value) {
          return false;
        }
      }
      // check y direction
      for (int i = -ypadding; i <= ypadding; i++) {
        int up = data[index + i * global_dimensions[0]];
        int down = data[index - i * global_dimensions[0]];
        if (up > center_value || down > center_value) {
          return false;
        }
      }
      return true;
    }
    // if sign is negative, then we are looking for minimum
    else if (sign < 0) {
      T center_value = data[index] + tol;
      // check x direction
      for (int i = -xpadding; i <= xpadding; i++) {
        int left = data[index + i];
        int right = data[index - i];
        if (left < center_value || right < center_value) {
          return false;
        }
      }
      // check y direction
      for (int i = -ypadding; i <= ypadding; i++) {
        int up = data[index + i * global_dimensions[0]];
        int down = data[index - i * global_dimensions[0]];
        if (up < center_value || down < center_value) {
          return false;
        }
      }
      return true;
    } else {
      return false;
    }
  }

  // 3d version
  // This is the same as 2d version, but with z direction
  bool padding_check(const T *data, const int index, const int sign,
                     const int xpadding, const int ypadding, const int zpadding,
                     T tol = 1e-5) {
    // T center_value = data[index];
    // if sign is positive, then we are looking for maximum
    if (sign > 0) {
      T center_value = data[index] - tol;
      ;
      // check x direction
      for (int i = -xpadding; i <= xpadding; i++) {
        int left = data[index + i];
        int right = data[index - i];
        if (left > center_value || right > center_value) {
          return false;
        }
      }
      // check y direction
      for (int i = -ypadding; i <= ypadding; i++) {
        int up = data[index + i * global_dimensions[0]];
        int down = data[index - i * global_dimensions[0]];
        if (up > center_value || down > center_value) {
          return false;
        }
      }
      // check z direction
      for (int i = -zpadding; i <= zpadding; i++) {
        int front =
            data[index + i * global_dimensions[0] * global_dimensions[1]] + tol;
        int back =
            data[index - i * global_dimensions[0] * global_dimensions[1]] + tol;
        if (front > center_value || back > center_value) {
          return false;
        }
      }
      return true;
    } else if (sign < 0) {
      T center_value = data[index] + tol;
      // check x direction
      for (int i = -xpadding; i <= xpadding; i++) {
        int left = data[index + i];
        int right = data[index - i];
        if (left < center_value || right < center_value) {
          return false;
        }
      }
      // check y direction
      for (int i = -ypadding; i <= ypadding; i++) {
        int up = data[index + i * global_dimensions[0]];
        int down = data[index - i * global_dimensions[0]];
        if (up < center_value || down < center_value) {
          return false;
        }
      }
      // check z direction
      for (int i = -zpadding; i <= zpadding; i++) {
        int front =
            data[index + i * global_dimensions[0] * global_dimensions[1]];
        int back =
            data[index - i * global_dimensions[0] * global_dimensions[1]];
        if (front < center_value || back < center_value) {
          return false;
        }
      }
      return true;
    } else {
      return false;
    }
  }

  int pattern_match_global(std::vector<int> &match_index,
                           std::vector<int> &xpaddings,
                           std::vector<int> &ypaddings, int &match_count) {

    match_index.resize(0);
    xpaddings.resize(0);
    ypaddings.resize(0);
    match_count = 0;
    for (int i = 0; i < num_elements; i++) {
      if (critical_points_map[i] == 4 || critical_points_map[i] == -4) {
        int xpadding = 0;
        int ypadding = 0;
        auto match_result =
            pattern_match2d(critical_points_map, i % global_dimensions[0],
                            i / global_dimensions[0], xpadding, ypadding);
        if (match_result == 1) {
          match_index.push_back(i);
          xpaddings.push_back(xpadding);
          ypaddings.push_back(ypadding);
          match_count++;
        }
      }
    }
    // std::cout << "match_count: " << match_count << std::endl;
    return match_count;
  }
  int pattern_match_global(std::vector<int> &match_index,
                           std::vector<int> &xpaddings,
                           std::vector<int> &ypaddings,
                           std::vector<int> &zpaddings, int &match_count) {

    match_index.resize(0);
    xpaddings.resize(0);
    ypaddings.resize(0);
    zpaddings.resize(0);
    match_count = 0;
    for (int i = 0; i < num_elements; i++) {
      if (critical_points_map[i] == 6 || critical_points_map[i] == -6) {
        int xpadding = 0;
        int ypadding = 0;
        int zpadding = 0;
        int idx = i % global_dimensions[0];
        int idy = (i / global_dimensions[0]) % global_dimensions[1];
        int idz = i / (global_dimensions[0] * global_dimensions[1]);
        auto match_result = pattern_match3d(critical_points_map, idx, idy, idz,
                                            xpadding, ypadding, zpadding);
        if (match_result == 1) {
          match_index.push_back(i);
          xpaddings.push_back(xpadding);
          ypaddings.push_back(ypadding);
          zpaddings.push_back(zpadding);
          match_count++;
        }
      }
    }
    // std::cout << "match_count: " << match_count << std::endl;
    return match_count;
  }

  int pattern_match_global(std::vector<int> &match_index,
                           std::vector<int> &xpaddings,
                           std::vector<int> &ypaddings,
                           std::vector<T> &interp_errors, int &match_count,
                           T interp_threshold) {

    match_index.reserve(num_elements / 4);
    xpaddings.reserve(num_elements / 4);
    ypaddings.reserve(num_elements / 4);
    interp_errors.reserve(num_elements / 4);
    match_count = 0;
    for (int i = 0; i < num_elements; i++) {
      if (critical_points_map[i] == 4 || critical_points_map[i] == -4) {
        int xpadding = 0;
        int ypadding = 0;
        T interp_error = 0;
        auto match_result = pattern_match2d_interp(
            critical_points_map, i % global_dimensions[0],
            i / global_dimensions[0], xpadding, ypadding, interp_error,
            interp_threshold);
        if (match_result == 1) {
          match_index.push_back(i);
          xpaddings.push_back(xpadding);
          ypaddings.push_back(ypadding);
          interp_errors.push_back(interp_error);
          match_count++;
        }
      }
    }
    // std::cout << "match_count: " << match_count << std::endl;
    return match_count;
  }

  int pattern_match_global(std::vector<int> &match_index,
                           std::vector<int> &xpaddings,
                           std::vector<int> &ypaddings,
                           std::vector<int> &zpaddings,
                           std::vector<T> &interp_errors, int &match_count,
                           T interp_threshold) {
    match_index.reserve(num_elements / 4);
    xpaddings.reserve(num_elements / 4);
    ypaddings.reserve(num_elements / 4);
    zpaddings.reserve(num_elements / 4);
    interp_errors.reserve(num_elements / 4);
    match_count = 0;
    for (int i = 0; i < num_elements; i++) {
      if (critical_points_map[i] == 6 || critical_points_map[i] == -6) {
        int xpadding = 0;
        int ypadding = 0;
        int zpadding = 0;
        T interp_error = 0;
        int idx = i % global_dimensions[0];
        int idy = (i / global_dimensions[0]) % global_dimensions[1];
        int idz = i / (global_dimensions[0] * global_dimensions[1]);
        auto match_result = pattern_match3d_interp(
            critical_points_map, idx, idy, idz, xpadding, ypadding, zpadding,
            interp_error, interp_threshold);
        if (match_result == 1) {
          match_index.push_back(i);
          xpaddings.push_back(xpadding);
          ypaddings.push_back(ypadding);
          zpaddings.push_back(zpadding);
          interp_errors.push_back(interp_error);
          match_count++;
        }
      }
    }
    // std::cout << "match_count: " << match_count << std::endl;
    return match_count;
  }

  bool try_match2d(std::vector<int> &cpmap, int index, int &xpadding,
                   int &ypadding, int max_padding = 5) {
    return pattern_match2d(cpmap, index % global_dimensions[0],
                           index / global_dimensions[0], xpadding, ypadding,
                           max_padding);
  }

  bool try_match3d(std::vector<int> &cpmap, int index, int &xpadding,
                   int &ypadding, int &zpadding, int max_padding = 7) {

    int idx = index % global_dimensions[0];
    int idy = (index / global_dimensions[0]) % global_dimensions[1];
    int idz = index / (global_dimensions[0] * global_dimensions[1]);
    return pattern_match3d(cpmap, idx, idy, idz, xpadding, ypadding, zpadding,
                           max_padding);
  }

  bool try_match2d(std::vector<int> &cpmap, int idx, int idy, int &xpadding,
                   int &ypadding, int max_padding = 5) {
    return pattern_match2d(cpmap, idx, idy, xpadding, ypadding, max_padding);
  }

  bool try_match3d(std::vector<int> &cpmap, int idx, int idy, int idz,
                   int &xpadding, int &ypadding, int &zpadding,
                   int max_padding = 5) {
    return pattern_match3d(cpmap, idx, idy, idz, xpadding, ypadding, zpadding,
                           max_padding);
  }

  bool try_match2d_interp(std::vector<int> &cpmap, int index, int &xpadding,
                          int &ypadding, T &interp_error, T interp_threshold,
                          int max_padding = 5) {
    return pattern_match2d_interp(
        cpmap, index % global_dimensions[0], index / global_dimensions[0],
        xpadding, ypadding, interp_error, interp_threshold, max_padding);
  }

  bool try_match3d_interp(std::vector<int> &cpmap, int index, int &xpadding,
                          int &ypadding, int &zpadding, T &interp_error,
                          T interp_threshold, int max_padding = 7) {
    int idx = index % global_dimensions[0];
    int idy = (index / global_dimensions[0]) % global_dimensions[1];
    int idz = index / (global_dimensions[0] * global_dimensions[1]);
    return pattern_match3d_interp(cpmap, idx, idy, idz, xpadding, ypadding,
                                  zpadding, interp_error, interp_threshold,
                                  max_padding);
  }

  void set_local_range(T range) { local_val_tol = range; }
  void set_cp_map_tol(T tol) { cp_map_tol = tol; }

  int get_sign_single_point(T *data, int index, double tol) {
    return calculate_single_point(data, index, tol);
  }

private:
  std::vector<int> global_dimensions;
  size_t num_elements;
  T *inputdata;
  int N;
  std::vector<int> critical_points_map;
  InterpolationWalk<T> *interps;
  double flush_threshold = 1e-7;
  double local_val_tol = 1e-3;
  double cp_map_tol = 1e-5;

  inline bool is_valid_index(int idx, int idy) {
    return idx >= 0 && idx < global_dimensions[0] && idy >= 0 &&
           idy < global_dimensions[1];
  }

  inline bool is_valid_index(int idx, int idy, int idz) {
    return idx >= 0 && idx < global_dimensions[0] && idy >= 0 &&
           idy < global_dimensions[1] && idz >= 0 && idz < global_dimensions[2];
  }

  void calculate_critial_points_map() {
    this->critical_points_map.resize(num_elements, 0);
    for (int i = 0; i < num_elements; i++) {
      this->critical_points_map[i] =
          calculate_single_point(this->inputdata, i, this->cp_map_tol);
    }
  }

  inline T *get_value(int idx, int idy) {
    return inputdata + idx + idy * global_dimensions[0];
  }

  inline T *get_value(int idx, int idy, int idz) {

    return inputdata + idx + idy * global_dimensions[0] +
           idz * global_dimensions[0] * global_dimensions[1];
  }

  // Calculate the gradient degree for a single point.
  // 4 if it is a local max or max.
  // the cornor cases is left for future work(the boundary points).
  // int calculate_single_point(int global_index) {

  //   if (N == 2) {
  //     T tol = cp_map_tol;

  //     int idx = global_index % global_dimensions[0];
  //     int idy = global_index / global_dimensions[0];

  //     if (idx == 0 || idy == 0 || idx == global_dimensions[0] - 1 ||
  //         idy == global_dimensions[1] - 1) {
  //       return 0;
  //     } else {
  //       T current_value = data[global_index];
  //       T value_minus_tol = current_value - tol;
  //       T value_plus_tol = current_value + tol;
  //       // T up = *get_value(idx, idy + 1);
  //       // T down = *get_value(idx, idy - 1);
  //       // T left = *get_value(idx - 1, idy);
  //       // T right = *get_value(idx + 1, idy);
  //       T up = data[global_index + global_dimensions[0]];
  //       T down = data[global_index - global_dimensions[0]];
  //       T left = data[global_index - 1];
  //       T right = data[global_index + 1];
  //       int sign1 = (value_minus_tol > up) - (value_plus_tol < up);
  //       int sign2 = (value_minus_tol > down) - (value_plus_tol < down);
  //       int sign3 = (value_minus_tol > left) - (value_plus_tol < left);
  //       int sign4 = (value_minus_tol > right) - (value_plus_tol < right);
  //       return sign1 + sign2 + sign3 + sign4;
  //     }
  //   } else if (N == 3) {
  //     T tol = cp_map_tol;

  //     int idx = global_index % global_dimensions[0];
  //     int idy = (global_index / global_dimensions[0]) % global_dimensions[1];
  //     int idz = global_index / (global_dimensions[0] * global_dimensions[1]);

  //     if (idx == 0 || idy == 0 || idz == 0 || idx == global_dimensions[0] - 1
  //     ||
  //         idy == global_dimensions[1] - 1 || idz == global_dimensions[2] - 1)
  //         {
  //       return 0;
  //     } else {
  //       T current_value = data[global_index];
  //       T value_minus_tol = current_value - tol;
  //       T value_plus_tol = current_value + tol;
  //       // T up = *get_value(idx, idy + 1, idz);
  //       // T down = *get_value(idx, idy - 1, idz);
  //       // T left = *get_value(idx - 1, idy, idz);
  //       // T right = *get_value(idx + 1, idy, idz);
  //       // T front = *get_value(idx, idy, idz + 1);
  //       // T back = *get_value(idx, idy, idz - 1);
  //       T up = data[global_index + global_dimensions[0]];
  //       T down = data[global_index - global_dimensions[0]];
  //       T left = data[global_index - 1];
  //       T right = data[global_index + 1];
  //       int dim0xdim1 = global_dimensions[0] * global_dimensions[1];
  //       T front = data[global_index + dim0xdim1];
  //       T back = data[global_index - dim0xdim1];
  //       int sign1 = (value_minus_tol > up) - (value_plus_tol < up);
  //       int sign2 = (value_minus_tol > down) - (value_plus_tol < down);
  //       int sign3 = (value_minus_tol > left) - (value_plus_tol < left);
  //       int sign4 = (value_minus_tol > right) - (value_plus_tol < right);
  //       int sign5 = (value_minus_tol > front) - (value_plus_tol < front);
  //       int sign6 = (value_minus_tol > back) - (value_plus_tol < back);
  //       return sign1 + sign2 + sign3 + sign4 + sign5 + sign6;
  //     }
  //   } else {
  //     std::cout << "CriticalPointsCalculator: N is not 2 or 3" << std::endl;
  //     return 0;
  //   }
  // }

  int calculate_single_point(T *data, int global_index, double tol) {
    if (N == 2) {
      // T tol = cp_map_tol;

      int idx = global_index % global_dimensions[0];
      int idy = global_index / global_dimensions[0];

      if (idx == 0 || idy == 0 || idx == global_dimensions[0] - 1 ||
          idy == global_dimensions[1] - 1) {
        return 0;
      } else {
        T current_value = data[global_index];
        T value_minus_tol = current_value - tol;
        T value_plus_tol = current_value + tol;
        // T up = *get_value(idx, idy + 1);
        // T down = *get_value(idx, idy - 1);
        // T left = *get_value(idx - 1, idy);
        // T right = *get_value(idx + 1, idy);
        T up = data[global_index + global_dimensions[0]];
        T down = data[global_index - global_dimensions[0]];
        T left = data[global_index - 1];
        T right = data[global_index + 1];
        int sign1 = (value_minus_tol > up) - (value_plus_tol < up);
        int sign2 = (value_minus_tol > down) - (value_plus_tol < down);
        int sign3 = (value_minus_tol > left) - (value_plus_tol < left);
        int sign4 = (value_minus_tol > right) - (value_plus_tol < right);
        return sign1 + sign2 + sign3 + sign4;
      }
    } else if (N == 3) {
      // T tol = cp_map_tol;

      int idx = global_index % global_dimensions[0];
      int idy = (global_index / global_dimensions[0]) % global_dimensions[1];
      int idz = global_index / (global_dimensions[0] * global_dimensions[1]);

      if (idx == 0 || idy == 0 || idz == 0 || idx == global_dimensions[0] - 1 ||
          idy == global_dimensions[1] - 1 || idz == global_dimensions[2] - 1) {
        return 0;
      } else {
        T current_value = data[global_index];
        T value_minus_tol = current_value - tol;
        T value_plus_tol = current_value + tol;
        // T up = *get_value(idx, idy + 1, idz);
        // T down = *get_value(idx, idy - 1, idz);
        // T left = *get_value(idx - 1, idy, idz);
        // T right = *get_value(idx + 1, idy, idz);
        // T front = *get_value(idx, idy, idz + 1);
        // T back = *get_value(idx, idy, idz - 1);
        T up = data[global_index + global_dimensions[0]];
        T down = data[global_index - global_dimensions[0]];
        T left = data[global_index - 1];
        T right = data[global_index + 1];
        int dim0xdim1 = global_dimensions[0] * global_dimensions[1];
        T front = data[global_index + dim0xdim1];
        T back = data[global_index - dim0xdim1];
        int sign1 = (value_minus_tol > up) - (value_plus_tol < up);
        int sign2 = (value_minus_tol > down) - (value_plus_tol < down);
        int sign3 = (value_minus_tol > left) - (value_plus_tol < left);
        int sign4 = (value_minus_tol > right) - (value_plus_tol < right);
        int sign5 = (value_minus_tol > front) - (value_plus_tol < front);
        int sign6 = (value_minus_tol > back) - (value_plus_tol < back);
        return sign1 + sign2 + sign3 + sign4 + sign5 + sign6;
      }
    } else {
      std::cout << "CriticalPointsCalculator: N is not 2 or 3" << std::endl;
      return 0;
    }
  }

  //
  //  pattern matching

  inline int map_value(int idx, int idy) {
    return critical_points_map[idx + idy * global_dimensions[0]];
  }

  inline int map_value(int idx, int idy, int idz) {
    return critical_points_map[idx + idy * global_dimensions[0] +
                               idz * global_dimensions[0] *
                                   global_dimensions[1]];
  }

  inline int map_value(std::vector<int> &cmap, int idx, int idy) {
    return cmap[idx + idy * global_dimensions[0]];
  }

  inline int map_value(std::vector<int> &cmap, int idx, int idy, int idz) {
    return cmap[idx + idy * global_dimensions[0] +
                idz * global_dimensions[0] * global_dimensions[1]];
  }
  //
  //  pattern matching
  // Input
  //  cpmap: critical points map.
  //  idx, idy: index of the critical point.
  //  xpadding, ypadding: padding of the critical point.
  //  max_padding: maximum padding to search.
  // Output
  //  return for match.
  bool pattern_match2d(std::vector<int> &cpmap, int idx, int idy, int &xpadding,
                       int &ypadding, int max_padding = 5) {
    xpadding = 0;
    ypadding = 0;
    int global_index = idx + idy * global_dimensions[0];
    int dimx = global_dimensions[0];
    int dimy = global_dimensions[1];
    int pad_value = 0;
    for (int i = 1; i <= max_padding; ++i) {
      // dir1 index
      if ((is_valid_index(idx - i, idy)) && is_valid_index(idx + i, idy)) {
        if (map_value(cpmap, idx - i, idy) == 2 &&
            map_value(cpmap, idx + i, idy) == 2 &&
            map_value(cpmap, idx, idy) == 4) {
          xpadding = i;
        } else if (map_value(cpmap, idx - i, idy) == -2 &&
                   map_value(cpmap, idx + i, idy) == -2 &&
                   map_value(cpmap, idx, idy) == -4) {
          xpadding = i;
        } else {
          break;
        }
      }
    }

    if (xpadding == 0)
      return false;

    for (int i = 1; i <= max_padding; ++i) {
      // dir2 index
      if ((is_valid_index(idx, idy - i)) && is_valid_index(idx, idy + i)) {
        if (map_value(cpmap, idx, idy - i) == 2 &&
            map_value(cpmap, idx, idy + i) == 2 &&
            map_value(cpmap, idx, idy) == 4) {
          ypadding = i;
        } else if (map_value(cpmap, idx, idy - i) == -2 &&
                   map_value(cpmap, idx, idy + i) == -2 &&
                   map_value(cpmap, idx, idy) == -4) {
          ypadding = i;
        } else {
          break;
        }
      }
    }

    if (ypadding == 0)
      return false;

    // check the signs of the corner points
    int center = map_value(cpmap, idx, idy);
    int new_xpadding = xpadding;
    int new_ypadding = ypadding;
    if (center > 0) {
      for (int i = -xpadding; i <= xpadding; ++i) {
        if (i == 0)
          continue;
        for (int j = -ypadding; j <= ypadding; ++j) {
          if (j == 0)
            continue;
          if (map_value(cpmap, idx + i, idy + j) != 0) {
            new_ypadding = std::min(std::abs(j) - 1, new_ypadding);
            new_xpadding = std::min(std::abs(i) - 1, new_xpadding);
          }
        }
      }
    } else if (center < 0) {
      for (int i = -xpadding; i <= xpadding; ++i) {
        if (i == 0)
          continue;
        for (int j = -ypadding; j <= ypadding; ++j) {
          if (j == 0)
            continue;
          if (map_value(cpmap, idx + i, idy + j) != 0) {
            new_ypadding = std::min(std::abs(j) - 1, new_ypadding);
            new_xpadding = std::min(std::abs(i) - 1, new_xpadding);
          }
        }
      }
    } else {
      return false;
    }

    xpadding = new_xpadding;
    ypadding = new_ypadding;

    T local_abs_max = 0;
    auto get_local_value_range = [this, idx, idy, xpadding,
                                  ypadding](T &local_abs_max) {
      T min = std::numeric_limits<T>::max();
      T max = -std::numeric_limits<T>::max();
      for (int i = idx - xpadding; i <= idx + xpadding; ++i) {
        for (int j = idy - ypadding; j <= idy + ypadding; ++j) {
          T value = *get_value(i, j);
          if (value < min)
            min = value;
          if (value > max)
            max = value;
        }
      }
      local_abs_max = std::max(std::abs(min), std::abs(max));
      return (T)(max - min);
    };

    if (xpadding != 0 && ypadding != 0) {
      // T center_diff = 0 ;
      // center_diff = std::max(center_diff, std::abs(*get_value(idx,idy) -
      // *get_value(idx+xpadding,idy))); center_diff = std::max(center_diff,
      // std::abs(*get_value(idx,idy) - *get_value(idx-xpadding,idy)));
      // center_diff = std::max(center_diff, std::abs(*get_value(idx,idy) -
      // *get_value(idx,idy+ypadding))); center_diff = std::max(center_diff,
      // std::abs(*get_value(idx,idy) - *get_value(idx,idy-ypadding)));

      T local_value_range = get_local_value_range(local_abs_max);
      if (local_value_range > local_val_tol)
        // if (center_diff > local_val_tol)
        return 1;
      else
        return 0;
    } else
      return 0;
  }

  bool pattern_match3d(std::vector<int> &cpmap, int idx, int idy, int idz,
                       int &xpadding, int &ypadding, int &zpadding,
                       int max_padding = 5) {
    xpadding = 0;
    ypadding = 0;
    zpadding = 0;
    int padding_value = 0;

    int global_index = idx + idy * global_dimensions[0] +
                       idz * global_dimensions[0] * global_dimensions[1];
    int count = 0;
    // std::cout << "x y z " << idx << " " << idy << " " << idz << std::endl;
    for (int i = 1; i <= max_padding; ++i) {
      if (is_valid_index(idx - i, idy, idz) &&
          is_valid_index(idx + i, idy, idz)) {
        if ((map_value(cpmap, idx - i, idy, idz) == -4) &&
            (map_value(cpmap, idx + i, idy, idz) == -4) &&
            (map_value(cpmap, idx, idy, idz) == -6)) {
          xpadding = i;
        } else if ((map_value(cpmap, idx - i, idy, idz) == 4) &&
                   (map_value(cpmap, idx + i, idy, idz) == 4) &&
                   (map_value(cpmap, idx, idy, idz) == 6)) {
          xpadding = i;
        }
        { break; }
      }
    }
    if (xpadding == 0)
      return false;

    for (int i = 1; i <= max_padding; ++i) {
      if (is_valid_index(idx, idy - i, idz) &&
          is_valid_index(idx, idy + i, idz)) {
        if ((map_value(cpmap, idx, idy - i, idz) == -4) &&
            (map_value(cpmap, idx, idy + i, idz) == -4) &&
            (map_value(cpmap, idx, idy, idz) == -6)) {
          ypadding = i;
        } else if ((map_value(cpmap, idx, idy - i, idz) == 4) &&
                   (map_value(cpmap, idx, idy + i, idz) == 4) &&
                   (map_value(cpmap, idx, idy, idz) == 6)) {
          ypadding = i;
        }
        { break; }
      }
    }
    if (ypadding == 0)
      return false;

    for (int i = 1; i <= max_padding; ++i) {
      if (is_valid_index(idx, idy, idz - i) &&
          is_valid_index(idx, idy, idz + i)) {
        if ((map_value(cpmap, idx, idy, idz - i) == -4) &&
            (map_value(cpmap, idx, idy, idz + i) == -4) &&
            (map_value(cpmap, idx, idy, idz) == -6)) {
          zpadding = i;
        } else if ((map_value(cpmap, idx, idy, idz - i) == 4) &&
                   (map_value(cpmap, idx, idy, idz + i) == 4) &&
                   (map_value(cpmap, idx, idy, idz) == 6)) {
          zpadding = i;
        }
        { break; }
      }
    }
    if (zpadding == 0)
      return false;

    int new_xpadding = xpadding;
    int new_ypadding = ypadding;
    int new_zpadding = zpadding;

    // check the signs of the corner points
    // they have to be the same sign or 0
    int center_sign = map_value(cpmap, idx, idy, idz);
    if (center_sign > 0) {
      for (int i = -xpadding; i <= xpadding; i++) {
        if (i == 0)
          continue;
        for (int j = -ypadding; j <= ypadding; j++) {
          if (j == 0)
            continue;
          for (int k = -zpadding; k <= zpadding; k++) {
            if (k == 0)
              continue;
            if (map_value(cpmap, idx + i, idy + j, idz + k) != 0) {
              new_ypadding = std::min(std::abs(j) - 1, new_ypadding);
              new_xpadding = std::min(std::abs(i) - 1, new_xpadding);
              new_zpadding = std::min(std::abs(k) - 1, new_zpadding);
            }
          }
        }
      }
    } else if (center_sign < 0) {
      for (int i = -xpadding; i <= xpadding; i++) {
        if (i == 0)
          continue;
        for (int j = -ypadding; j <= ypadding; j++) {
          if (j == 0)
            continue;
          for (int k = -zpadding; k <= zpadding; k++) {
            if (k == 0)
              continue;
            if (map_value(cpmap, idx + i, idy + j, idz + k) != 0) {
              new_ypadding = std::min(std::abs(j) - 1, new_ypadding);
              new_xpadding = std::min(std::abs(i) - 1, new_xpadding);
              new_zpadding = std::min(std::abs(k) - 1, new_zpadding);
            }
          }
        }
      }
    } else {
      return false;
    }

    xpadding = new_xpadding;
    ypadding = new_ypadding;
    zpadding = new_zpadding;

    // if (idx == 103 && idy == 211 && idz == 11) {
    //   std::cout << "global index = " << global_index << std::endl;
    //   std::cout << "map_value = " << map_value(idx, idy, idz) << std::endl;
    //   std::cout << "xpadding = " << xpadding << std::endl;
    //   std::cout << "ypadding = " << ypadding << std::endl;
    //   std::cout << "zpadding = " << zpadding << std::endl;
    // }
    T local_abs_max = 0;
    auto get_local_value_range = [this, idx, idy, idz, xpadding, ypadding,
                                  zpadding](T &local_abs_max) {
      T min = std::numeric_limits<T>::max();
      T max = -std::numeric_limits<T>::max();
      for (int i = idx - xpadding; i <= idx + xpadding; ++i) {
        for (int j = idy - ypadding; j <= idy + ypadding; ++j) {
          for (int k = idz - zpadding; k <= idz + zpadding; ++k) {
            T value = *get_value(i, j, k);
            if (value < min)
              min = value;
            if (value > max)
              max = value;
          }
        }
      }
      local_abs_max = std::max(std::abs(min), std::abs(max));
      return (T)(max - min);
    };

    if ((xpadding != 0 && ypadding != 0 && zpadding != 0)) {
      T local_value_range = get_local_value_range(local_abs_max);
      if (local_value_range > local_val_tol)
        return true;
      else
        return false;
    } else
      return false;
  }

  bool pattern_match2d_interp(std::vector<int> &cpmap, int idx, int idy,
                              int &xpadding, int &ypadding,
                              T &local_interp_error, T interp_error_threshold,
                              int max_padding = 5) {
    xpadding = 0;
    ypadding = 0;
    local_interp_error = 0;
    int padding_value = 0;

    int global_index = idx + idy * global_dimensions[0];

    for (int i = 1; i <= max_padding; ++i) {
      // dir1 index
      if ((is_valid_index(idx - i, idy)) && is_valid_index(idx + i, idy)) {
        if (map_value(cpmap, idx - i, idy) == 2 &&
            map_value(cpmap, idx + i, idy) == 2 &&
            map_value(cpmap, idx, idy) == 4) {
          xpadding = i;
          // local_interp_error += interps->interp_walk(idx - i, idy);
          // local_interp_error += interps->interp_walk(idx + i, idy);
          local_interp_error += std::pow(interps->interp_walk(idx - i, idy), 2);
          local_interp_error += std::pow(interps->interp_walk(idx + i, idy), 2);
        } else if (map_value(cpmap, idx - i, idy) == -2 &&
                   map_value(cpmap, idx + i, idy) == -2 &&
                   map_value(cpmap, idx, idy) == -4) {
          xpadding = i;
          // local_interp_error += interps->interp_walk(idx - i, idy);
          // local_interp_error += interps->interp_walk(idx + i, idy);
          local_interp_error += std::pow(interps->interp_walk(idx - i, idy), 2);
          local_interp_error += std::pow(interps->interp_walk(idx + i, idy), 2);
        } else {
          break;
        }
      }
    }

    if (xpadding == 0)
      return false;

    for (int i = 1; i <= max_padding; ++i) {
      // dir1 index
      if ((is_valid_index(idx, idy - i)) && is_valid_index(idx, idy + i)) {
        if (map_value(cpmap, idx, idy - i) == 2 &&
            map_value(cpmap, idx, idy + i) == 2 &&
            map_value(cpmap, idx, idy) == 4) {
          ypadding = i;
          // local_interp_error += interps->interp_walk(idx, idy - i);
          // local_interp_error += interps->interp_walk(idx, idy + i);
          local_interp_error += std::pow(interps->interp_walk(idx, idy - i), 2);
          local_interp_error += std::pow(interps->interp_walk(idx, idy + i), 2);
        } else if (map_value(cpmap, idx, idy - i) == -2 &&
                   map_value(cpmap, idx, idy + i) == -2 &&
                   map_value(cpmap, idx, idy) == -4) {
          ypadding = i;
          local_interp_error += std::pow(interps->interp_walk(idx, idy - i), 2);
          local_interp_error += std::pow(interps->interp_walk(idx, idy + i), 2);
        } else {
          break;
        }
      }
    }
    if (ypadding == 0)
      return false;

    local_interp_error += std::pow(interps->interp_walk(idx, idy), 2);

    // T local_abs_max = 0;
    // auto get_local_value_range = [this, idx, idy, xpadding,
    //                               ypadding](T &local_abs_max) {
    //   T min = std::numeric_limits<T>::max();
    //   T max = -std::numeric_limits<T>::max();
    //   for (int i = idx - xpadding; i <= idx + xpadding; ++i) {
    //     for (int j = idy - ypadding; j <= idy + ypadding; ++j) {
    //       T value = *get_value(i, j);
    //       if (value < min)
    //         min = value;
    //       if (value > max)
    //         max = value;
    //     }
    //   }
    //   local_abs_max = std::max(std::abs(min), std::abs(max));
    //   return (T)(max - min);
    // };

    T local_abs_max = 0;
    T local_std = 0;
    auto get_local_value_range = [this, idx, idy, xpadding,
                                  ypadding](T &local_abs_max, T &local_std) {
      T min = std::numeric_limits<T>::max();
      T max = -std::numeric_limits<T>::max();
      T sum = 0;
      int local_count = ((2 * xpadding + 1) * (2 * ypadding + 1));
      int interp_count = 0;
      for (int i = idx - xpadding; i <= idx + xpadding; ++i) {
        for (int j = idy - ypadding; j <= idy + ypadding; ++j) {
          T value = *get_value(i, j);
          if (value < min)
            min = value;
          if (value > max)
            max = value;
          sum += value;
        }
      }
      T mean = sum / local_count;
      for (int i = idx - xpadding; i <= idx + xpadding; ++i) {
        for (int j = idy - ypadding; j <= idy + ypadding; ++j) {
          T value = *get_value(i, j);
          local_std += (value - mean) * (value - mean);
        }
      }
      local_std = std::sqrt(local_std / local_count);
      local_abs_max = std::max(std::abs(min), std::abs(max));
      return (T)(max - min);
    };

    if (xpadding != 0 && ypadding != 0) {
      T local_value_range = get_local_value_range(local_abs_max, local_std);
      local_interp_error /= (T)(2 * xpadding + 2 * ypadding + 1);
      local_interp_error = std::sqrt(local_interp_error);
      if (local_value_range > local_val_tol &&
          local_interp_error < interp_error_threshold)
        // if (local_std > local_val_tol &&
        //     local_interp_error < interp_error_threshold)
        return true;
      else
        return false;
    } else
      return false;
  }

  bool pattern_match3d_interp(std::vector<int> &cpmap, int idx, int idy,
                              int idz, int &xpadding, int &ypadding,
                              int &zpadding, T &interp_error,
                              T interp_error_threshold, int max_padding = 5) {
    xpadding = 0;
    ypadding = 0;
    zpadding = 0;
    interp_error = 0;
    int padding_value = 0;
    int global_index = idx + idy * global_dimensions[0] +
                       idz * global_dimensions[0] * global_dimensions[1];
    for (int i = 1; i <= max_padding; ++i) {
      if (is_valid_index(idx - i, idy, idz) &&
          is_valid_index(idx + i, idy, idz)) {
        if ((map_value(cpmap, idx - i, idy, idz) == -4) &&
            (map_value(cpmap, idx + i, idy, idz) == -4) &&
            (map_value(cpmap, idx, idy, idz) == -6)) {
          xpadding = i;
          // interp_error += interps->interp_walk(idx - i, idy, idz);
          // interp_error += interps->interp_walk(idx + i, idy, idz);
          interp_error += std::pow(interps->interp_walk(idx - i, idy, idz), 2);
          interp_error += std::pow(interps->interp_walk(idx + i, idy, idz), 2);
        } else if ((map_value(cpmap, idx - i, idy, idz) == 4) &&
                   (map_value(cpmap, idx + i, idy, idz) == 4) &&
                   (map_value(cpmap, idx, idy, idz) == 6)) {
          xpadding = i;
          // interp_error += interps->interp_walk(idx - i, idy, idz);
          // interp_error += interps->interp_walk(idx + i, idy, idz);
          interp_error += std::pow(interps->interp_walk(idx - i, idy, idz), 2);
          interp_error += std::pow(interps->interp_walk(idx + i, idy, idz), 2);
        }
        { break; }
      }
    }
    // if (idx == 103 && idy == 211 && idz == 11) {
    //     std::cout << "global index = " << global_index << std::endl;
    //     std::cout << "map_value = " << map_value(idx, idy, idz) << std::endl;
    //     std::cout << "xpadding = " << xpadding << std::endl;
    //     // std::cout << "local_value_range = " << local_value_range <<
    //     std::endl;
    // }

    if (xpadding == 0)
      return false;

    for (int i = 1; i <= max_padding; ++i) {
      if (is_valid_index(idx, idy - i, idz) &&
          is_valid_index(idx, idy + i, idz)) {
        if ((map_value(cpmap, idx, idy - i, idz) == -4) &&
            (map_value(cpmap, idx, idy + i, idz) == -4) &&
            (map_value(cpmap, idx, idy, idz) == -6)) {
          ypadding = i;
          // interp_error += interps->interp_walk(idx, idy - i, idz);
          // interp_error += interps->interp_walk(idx, idy + i, idz);
          interp_error += std::pow(interps->interp_walk(idx, idy - i, idz), 2);
          interp_error += std::pow(interps->interp_walk(idx, idy + i, idz), 2);
        } else if ((map_value(cpmap, idx, idy - i, idz) == 4) &&
                   (map_value(cpmap, idx, idy + i, idz) == 4) &&
                   (map_value(cpmap, idx, idy, idz) == 6)) {
          ypadding = i;
          // interp_error += interps->interp_walk(idx, idy - i, idz);
          // interp_error += interps->interp_walk(idx, idy + i, idz);
          interp_error += std::pow(interps->interp_walk(idx, idy - i, idz), 2);
          interp_error += std::pow(interps->interp_walk(idx, idy + i, idz), 2);
        }
        { break; }
      }
    }
    //       if (idx == 103 && idy == 211 && idz == 11) {
    //   std::cout << "ypadding = " << ypadding << std::endl;
    //   // std::cout << "local_value_range = " << local_value_range <<
    //   std::endl;
    // }
    if (ypadding == 0)
      return false;

    for (int i = 1; i <= max_padding; ++i) {
      if (is_valid_index(idx, idy, idz - i) &&
          is_valid_index(idx, idy, idz + i)) {
        if ((map_value(cpmap, idx, idy, idz - i) == -4) &&
            (map_value(cpmap, idx, idy, idz + i) == -4) &&
            (map_value(cpmap, idx, idy, idz) == -6)) {
          zpadding = i;
          // interp_error += interps->interp_walk(idx, idy, idz - i);
          // interp_error += interps->interp_walk(idx, idy, idz + i);
          interp_error += std::pow(interps->interp_walk(idx, idy, idz - i), 2);
          interp_error += std::pow(interps->interp_walk(idx, idy, idz + i), 2);
        } else if ((map_value(cpmap, idx, idy, idz - i) == 4) &&
                   (map_value(cpmap, idx, idy, idz + i) == 4) &&
                   (map_value(cpmap, idx, idy, idz) == 6)) {
          zpadding = i;
          // interp_error += interps->interp_walk(idx, idy, idz - i);
          // interp_error += interps->interp_walk(idx, idy, idz + i);
          interp_error += std::pow(interps->interp_walk(idx, idy, idz - i), 2);
          interp_error += std::pow(interps->interp_walk(idx, idy, idz + i), 2);
        }
        { break; }
      }
    }
    if (zpadding == 0)
      return false;

    //           if (idx == 103 && idy == 211 && idz == 11) {
    //   std::cout << "zpadding = " << zpadding << std::endl;
    //   // std::cout << "local_value_range = " << local_value_range <<
    //   std::endl;
    // }

    // if (zpadding == 0)
    //   return false;

    //      if (idx == 103 && idy == 211 && idz == 11) {
    //     std::cout << "global index = " << global_index << std::endl;
    //     std::cout << "map_value = " << map_value(idx, idy, idz) << std::endl;
    //     std::cout << "xpadding = " << xpadding << std::endl;
    //     std::cout << "ypadding = " << ypadding << std::endl;
    //     std::cout << "zpadding = " << zpadding << std::endl;
    //     // std::cout << "local_value_range = " << local_value_range <<
    //     std::endl;
    //   }

    T local_abs_max = 0;
    T local_std = 0;
    // auto get_local_value_range = [this, idx, idy, idz, xpadding, ypadding,
    //                               zpadding](T &local_abs_max, T &local_std) {
    //   T min = std::numeric_limits<T>::max();
    //   T max = -std::numeric_limits<T>::max();
    //   double sum = 0;
    //   int local_count = ((2 * xpadding + 1) * (2 * ypadding + 1) *
    //                   (2 * zpadding + 1));
    //   int interp_count = 0;
    //   for (int i = idx - xpadding; i <= idx + xpadding; ++i) {
    //     for (int j = idy - ypadding; j <= idy + ypadding; ++j) {
    //       for (int k = idz - zpadding; k <= idz + zpadding; ++k) {
    //         T value = *get_value(i, j, k);
    //         sum += value;
    //         if (value < min)
    //           min = value;
    //         if (value > max)
    //           max = value;
    //       }
    //     }
    //   }
    //   double mean = sum / local_count;
    //   for (int i = idx - xpadding; i <= idx + xpadding; ++i) {
    //     for (int j = idy - ypadding; j <= idy + ypadding; ++j) {
    //       for (int k = idz - zpadding; k <= idz + zpadding; ++k) {
    //         T value = *get_value(i, j, k);
    //         local_std += (value - mean) * (value - mean);
    //       }
    //     }
    //   }
    //   local_std = std::sqrt(local_std / local_count);
    //   local_abs_max = std::max(std::abs(min), std::abs(max));
    //   return (T)(max - min);
    // };

    auto get_local_value_range = [this, idx, idy, idz, xpadding, ypadding,
                                  zpadding](T &local_abs_max, T &local_std) {
      T min = std::numeric_limits<T>::max();
      T max = -std::numeric_limits<T>::max();
      double sum = 0;
      int local_count =
          ((2 * xpadding + 1) * (2 * ypadding + 1) * (2 * zpadding + 1));
      int interp_count = 0;
      for (int i = idx - xpadding; i <= idx + xpadding; ++i) {
        for (int j = idy - ypadding; j <= idy + ypadding; ++j) {
          for (int k = idz - zpadding; k <= idz + zpadding; ++k) {
            T value = *get_value(i, j, k);
            if (value < min) {
              min = value;
            }
            if (value > max) {
              max = value;
            }
          }
        }
      }
      local_abs_max = std::max(std::abs(min), std::abs(max));
      return (T)(max - min);
    };

    // auto get_local_value_range = [this, idx, idy, idz, xpadding, ypadding,
    //                               zpadding](T &local_abs_max) {
    //   T min = std::numeric_limits<T>::max();
    //   T max = -std::numeric_limits<T>::max();
    //   T range = 0;
    //   T value = *get_value(idx, idy, idz);
    //   T left = std::abs(*get_value(idx - xpadding, idy, idz)-value);
    //   T right = std::abs(*get_value(idx + xpadding, idy, idz)-value);
    //   T up = std::abs(*get_value(idx, idy - ypadding, idz)-value);
    //   T down = std::abs(*get_value(idx, idy + ypadding, idz)-value);
    //   T front = std::abs(*get_value(idx, idy, idz - zpadding)-value);
    //   T back = std::abs(*get_value(idx, idy, idz + zpadding)-value);
    //   range = std::max(left, right);
    //   range = std::max(range, up);
    //   range = std::max(range, down);
    //   range = std::max(range, front);
    //   range = std::max(range, back);
    //   // local_abs_max = range;
    //   return range;
    // };

    if ((xpadding != 0 && ypadding != 0 && zpadding != 0)) {

      T local_value_range = get_local_value_range(local_abs_max, local_std);
      interp_error += std::pow(interps->interp_walk(idx, idy, idz), 2);
      interp_error = std::sqrt(interp_error /
                               (2.0 * (xpadding + ypadding + zpadding) + 1));

      if (local_value_range > local_val_tol &&
          interp_error < interp_error_threshold) {
        // if (local_value_range > local_val_tol &&
        //     interp_error < interp_error_threshold) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }
};
} // namespace SZ