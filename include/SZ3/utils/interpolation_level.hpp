#include "Interpolators.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

// the dimension is the same order with the dimension of SZ3 global_dimensions 
namespace SZ {
template <typename T> class InterpolationLevel {
public:
  InterpolationLevel(int N, size_t *global_dimensions) {
    this->N = N;
    this->num_elements = 1;
    this->global_dimensions.resize(N);
    for (int i = 0; i < N; i++) {
      this->num_elements *= global_dimensions[i];
      this->global_dimensions[i] = global_dimensions[i]; // the last one is the fastest dimension
      // std::cout << "global_dimensions[" << i << "] = " << global_dimensions[i]
      //           << std::endl;
    }
  }
  InterpolationLevel() {}

  ~InterpolationLevel() {}

  inline int get_level(const size_t index) {
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
      int x = index % global_dimensions[1];
      int y = index / global_dimensions[1];
      while (((x & (stride - 1)) == 0) && ((y & (stride - 1)) == 0)) {
        level++;
        stride = 1 << level;
      }
      return level;
    } else if (N == 3) {
      int x = index % global_dimensions[2];
      int y = (index / global_dimensions[2]) % global_dimensions[1];
      int z = index / (global_dimensions[2] * global_dimensions[1]);
      while (((x & (stride - 1)) == 0) && ((y & (stride - 1)) == 0) &&
             ((z & (stride - 1)) == 0)) {
        level++;
        stride = 1 << level;
      }
      return level;
    } else {
      return -1;
    }
  }

  inline std::vector<size_t> get_coordinates(size_t index) {
    std::vector<size_t> coordinates(N, 0);
    if (N == 1) {
      coordinates[0] = index;
    } else if (N == 2) {
      coordinates[0] = index % global_dimensions[1];
      coordinates[1] = index / global_dimensions[1];
    } else if (N == 3) {
      coordinates[0] = index % global_dimensions[2];
      coordinates[1] = (index / global_dimensions[2]) % global_dimensions[1];
      coordinates[2] = index / (global_dimensions[2] * global_dimensions[1]);
    }
    return coordinates;
  }

  template <typename T_id> inline size_t get_index(T_id x, T_id y, T_id z) {
    size_t index = 0;
    index = x + y * global_dimensions[2] +
            z * global_dimensions[2] * global_dimensions[1];
    return index;
  }

  template <typename T_id> inline size_t get_index(T_id x, T_id y) {
    size_t index = 0;
    index = x + y * global_dimensions[1];
    return index;
  }

  template <typename T_map>
  inline int
  label_neighbor_points(const size_t index, std::vector<T_map> &map_data,
                        const T_map fill_value, const int fill_radius) {
    if (N == 2) {
      int x = index % global_dimensions[1];
      int y = index / global_dimensions[1];
      int x_start = std::max(0, x - fill_radius);
      int x_end = std::min<int>(global_dimensions[1] - 1, x + fill_radius);
      int y_start = std::max(0, y - fill_radius);
      int y_end = std::min<int>(global_dimensions[0] - 1, y + fill_radius);
      for (int i = x_start; i <= x_end; i++) {
        for (int j = y_start; j <= y_end; j++) {
          size_t neighbor_index = get_index(i, j);
          map_data[neighbor_index] = fill_value;
        }
      }
      return 0;
    } else if (N == 3) {
      int x = index % global_dimensions[2];
      int y = (index / global_dimensions[2]) % global_dimensions[1];
      int z = index / (global_dimensions[2] * global_dimensions[1]);
      int x_start = std::max(0, x - fill_radius);
      int x_end = std::min<int>(global_dimensions[2] - 1, x + fill_radius);
      int y_start = std::max(0, y - fill_radius);
      int y_end = std::min<int>(global_dimensions[1] - 1, y + fill_radius);
      int z_start = std::max(0, z - fill_radius);
      int z_end = std::min<int>(global_dimensions[0] - 1, z + fill_radius);
      for (int i = x_start; i <= x_end; i++) {
        for (int j = y_start; j <= y_end; j++) {
          for (int k = z_start; k <= z_end; k++) {
            size_t neighbor_index = get_index(i, j, k);
            map_data[neighbor_index] = fill_value;
          }
        }
      }
      return 0;
    } else {
      std::cout << "Error: unsupported dimension!" << std::endl;
      exit(0);
    }
  }

  template <typename T_map>
  int inline extract_index(const std::vector<T_map> &map_vector,
                           const T_map extract_value,
                           std::vector<size_t> &extract_index_vector) {
    extract_index_vector.reserve(map_vector.size() / 2);
    for (size_t i = 0; i < map_vector.size(); i++) {
      if (map_vector[i] == extract_value) {
        extract_index_vector.push_back(i);
      }
    }
    return 0;
  }

  template <typename T_map>
  int inline map_to_control_blocks(
      const std::vector<T_map> &map_vector,
      const std::vector<size_t> &control_points, const T_map extract_value,
      const int radius, std::vector<size_t> &extract_index_vector) {
    extract_index_vector.reserve(map_vector.size() / 2);
    if (N == 2) {
      for (auto &index : control_points) {
        if (map_vector[index] == extract_value) {
          // push all the points in the block
          for (int i = -radius; i <= radius; i++) {
            for (int j = -radius; j <= radius; j++) {
              size_t neighbor_index = i + j * global_dimensions[1];
              extract_index_vector.push_back(map_vector[neighbor_index]);
            }
          }
        }
      }
    } else if (N == 3) {
      for (auto &index : control_points) {
        if (map_vector[index] == extract_value) {
          // push all the points in the block
          for (int i = -radius; i <= radius; i++) {
            for (int j = -radius; j <= radius; j++) {
              for (int k = -radius; k <= radius; k++) {
                size_t neighbor_index =
                    i + j * global_dimensions[2] +
                    k * global_dimensions[2] * global_dimensions[1];
                extract_index_vector.push_back(map_vector[neighbor_index]);
              }
            }
          }
        }
      }
    }
    else 
    {
      std::cout << "Error: unsupported dimension!" << std::endl;
      exit(0);
    }
    return 0;
  }

  template <typename T_map> 
  int inline control_blocks_to_map(
      const std::vector<size_t> &control_points, const T_map fill_value,
      const int radius, std::vector<T_map> &map_vector) {     
    if (N == 2) {
      for (auto &index : control_points) {
        int x = index % global_dimensions[1];
        int y = index / global_dimensions[1];
        // change the value of the points in the block
        for (int i = -radius; i <= radius; i++) {
          for (int j = -radius; j <= radius; j++) {
            size_t neighbor_index = i + j * global_dimensions[0];
            map_vector[neighbor_index] = fill_value;
          }
        }
      }
    } else if (N == 3) {
      for (auto &index : control_points) {
        int x = index % global_dimensions[2];
        int y = (index / global_dimensions[2]) % global_dimensions[1];
        int z = index / (global_dimensions[2] * global_dimensions[1]);
        // change the value of the points in the block
        for (int i = -radius; i <= radius; i++) {
          for (int j = -radius; j <= radius; j++) {
            for (int k = -radius; k <= radius; k++) {
              size_t neighbor_index =
                  i + j * global_dimensions[2] +
                  k * global_dimensions[2] * global_dimensions[1];
              map_vector[neighbor_index] = fill_value;
            }
          }
        }
      }
    }
    else 
    {
      std::cout << "Error: unsupported dimension!" << std::endl;
      exit(0);
    }
    return 0;
  }

  // inline int get_level(const size_t index) {
  //   int level = 1;
  //   int stride = 1 << level;
  //   switch (N) {
  //     case 1:
  //       while ((index & (stride - 1)) == 0) {
  //         level++;
  //         stride <<= 1;
  //       }
  //       return level;
  //     case 2: {
  //       int x = index & (global_dimensions[0] - 1);
  //       int y = index >> __builtin_ctz(global_dimensions[0]);
  //       while ((x & (stride - 1)) == 0 && (y & (stride - 1)) == 0) {
  //         level++;
  //         stride <<= 1;
  //       }
  //       return level;
  //     }
  //     case 3: {
  //       int x = index & (global_dimensions[0] - 1);
  //       int y = (index >> __builtin_ctz(global_dimensions[0])) &
  //               (global_dimensions[1] - 1);
  //       int z = index >> (__builtin_ctz(global_dimensions[0]) +
  //                         __builtin_ctz(global_dimensions[1]));
  //       while ((x & (stride - 1)) == 0 && (y & (stride - 1)) == 0 &&
  //              (z & (stride - 1)) == 0) {
  //         level++;
  //         stride <<= 1;
  //       }
  //       return level;
  //     }
  //     default:
  //       return -1;
  //   }
  // }

private:
  int N;
  size_t num_elements;
  std::vector<size_t> global_dimensions;
};

} // namespace SZ