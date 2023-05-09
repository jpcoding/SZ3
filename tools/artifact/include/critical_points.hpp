#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

template <class T> class CriticalPointsCalculator {

public:
  CriticalPointsCalculator(T *inputdata, int N, int *global_dimensions) {
    this->data = inputdata;
    this->N = N;
    this->num_elements = 1;
    this->global_dimensions.resize(N);

    for (int i = 0; i < N; i++) {
      this->global_dimensions[i] = global_dimensions[i];
      this->num_elements *= global_dimensions[i];
    }
  }

  // Deconstructor
  ~CriticalPointsCalculator() {}

  // calculate critical points map and return it
  std::vector<int> get_critical_points_map() {
    calculate_critial_points_map();
    return critical_points_map;
  }

  // match pattern

  int pattern_match_global(std::vector<int> &match_index,
                           std::vector<int> &xpaddings,
                           std::vector<int> &ypaddings, int &match_count) {

    match_index.resize(0);
    xpaddings.resize(0);
    ypaddings.resize(0);
    match_count = 0;
    for (int i = 0; i < num_elements; i++) {
      if (critical_points_map[i] == 4) {
        int xpadding = 0;
        int ypadding = 0;
        auto match_result =
            pattern_match2d(critical_points_map, i % global_dimensions[0],
                            i / global_dimensions[0], xpadding, ypadding);
        // if (i == 74160) {
        //   std::cout << "i: " << i << std::endl;
        //   std::cout << i % global_dimensions[0] << std::endl;
        //   std::cout << i / global_dimensions[0] << std::endl;
        //   std::cout << "critical_points_map[i]: " << critical_points_map[i]
        //             << std::endl;
        //   std::cout << "map value:"
        //             << map_value(i % global_dimensions[0],
        //                          i / global_dimensions[0])
        //             << std::endl;
        //   std::cout << "get_value(i % global_dimensions[0], i / "
        //                "global_dimensions[0]): "
        //             << get_value(i % global_dimensions[0],
        //                          i / global_dimensions[0])
        //             << std::endl;

        //   pattern_match2d(critical_points_map, i % global_dimensions[0],
        //                          i / global_dimensions[0], xpadding,
        //                          ypadding);
        //   std::cout << "xpadding: " << xpadding << std::endl;
        //   std::cout << "ypadding: " << ypadding << std::endl;
        //   std::cout << "match_result: " << match_result << std::endl;
        // }

        if (match_result == 1) {
          match_index.push_back(i);
          xpaddings.push_back(xpadding);
          ypaddings.push_back(ypadding);
          match_count++;
        }
      }
    }
    std::cout << "match_count: " << match_count << std::endl;
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
      if (critical_points_map[i] == 6) {
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
    std::cout << "match_count: " << match_count << std::endl;
    return match_count;
  }

private:
  std::vector<int> global_dimensions;
  size_t num_elements;
  T *data;
  int N;
  std::vector<int> critical_points_map;

  inline bool is_valid_index(int idx, int idy) {
    return idx >= 0 && idx < global_dimensions[0] && idy >= 0 &&
           idy < global_dimensions[1];
  }

  inline bool is_valid_index(int idx, int idy, int idz) {
    return idx >= 0 && idx < global_dimensions[0] && idy >= 0 &&
           idy < global_dimensions[1] && idz >= 0 && idz < global_dimensions[2];
  }

  void calculate_critial_points_map() {
    critical_points_map.resize(num_elements, 0);
    for (int i = 0; i < num_elements; i++) {
      critical_points_map[i] = calculate_single_point(i);
    }
  }

  inline T get_value(int idx, int idy) {
    return data[idx + idy * global_dimensions[0]];
  }

  inline T get_value(int idx, int idy, int idz) {
    return data[idx + idy * global_dimensions[0] +
                idz * global_dimensions[0] * global_dimensions[1]];
  }

  // Calculate the gradient degree for a single point.
  // 4 if it is a local max or max.
  // the cornor cases is left for future work(the boundary points).
  int calculate_single_point(int global_index) {
    if (N == 2) {
      int idx = global_index % global_dimensions[0];
      int idy = global_index / global_dimensions[0];
      if (idx == 0 || idy == 0 || idx == global_dimensions[0] - 1 ||
          idy == global_dimensions[1] - 1) {
        // int sign2 = (data[global_index]> get_value(idx+1, idy))-
        // (data[global_index]< get_value(idx+1, idy)); int sign3 =
        // (data[global_index]> get_value(idx, idy-1))- (data[global_index]<
        // get_value(idx, idy-1)); int sign4 = (data[global_index]>
        // get_value(idx, idy+1))- (data[global_index]< get_value(idx, idy+1));
        // return std::abs(sign2+sign3+sign4);
        return 0;
      }
      // else if (idy==0)
      // {
      //     int sign1 = (data[global_index]> get_value(idx-1, idy)) -
      //     (data[global_index]< get_value(idx-1, idy)); int sign2 =
      //     (data[global_index]> get_value(idx+1, idy))- (data[global_index]<
      //     get_value(idx+1, idy)); int sign4 = (data[global_index]>
      //     get_value(idx, idy+1))- (data[global_index]< get_value(idx,
      //     idy+1)); return std::abs(sign1+sign2+sign4);
      // }
      // else if (idx == global_dimensions[0]-1 )
      // {
      //     int sign1 = (data[global_index]> get_value(idx-1, idy)) -
      //     (data[global_index]< get_value(idx-1, idy)); int sign3 =
      //     (data[global_index]> get_value(idx, idy-1))- (data[global_index]<
      //     get_value(idx, idy-1)); int sign4 = (data[global_index]>
      //     get_value(idx, idy+1))- (data[global_index]< get_value(idx,
      //     idy+1)); return std::abs(sign1+sign3+sign4);
      // }
      // else if (idy == global_dimensions[1]-1 )
      // {
      //     int sign1 = (data[global_index]> get_value(idx-1, idy)) -
      //     (data[global_index]< get_value(idx-1, idy)); int sign2 =
      //     (data[global_index]> get_value(idx+1, idy))- (data[global_index]<
      //     get_value(idx+1, idy)); int sign3 = (data[global_index]>
      //     get_value(idx, idy-1))- (data[global_index]< get_value(idx,
      //     idy-1)); return std::abs(sign1+sign2+sign3);
      // }
      else {
        int sign1 = (data[global_index] > get_value(idx - 1, idy)) -
                    (data[global_index] < get_value(idx - 1, idy));
        int sign2 = (data[global_index] > get_value(idx + 1, idy)) -
                    (data[global_index] < get_value(idx + 1, idy));
        int sign3 = (data[global_index] > get_value(idx, idy - 1)) -
                    (data[global_index] < get_value(idx, idy - 1));
        int sign4 = (data[global_index] > get_value(idx, idy + 1)) -
                    (data[global_index] < get_value(idx, idy + 1));
        return std::abs(sign1 + sign2 + sign3 + sign4);
      }
    } else if (N == 3) {
      int idx = global_index % global_dimensions[0];
      int idy = (global_index / global_dimensions[0]) % global_dimensions[1];
      int idz = global_index / (global_dimensions[0] * global_dimensions[1]);
      if (idx == 0 || idy == 0 || idz == 0 || idx == global_dimensions[0] - 1 ||
          idy == global_dimensions[1] - 1 || idz == global_dimensions[2] - 1) {
        return 0;
      } else {
        int sign1 = (data[global_index] > get_value(idx - 1, idy, idz)) -
                    (data[global_index] < get_value(idx - 1, idy, idz));
        int sign2 = (data[global_index] > get_value(idx + 1, idy, idz)) -
                    (data[global_index] < get_value(idx + 1, idy, idz));
        int sign3 = (data[global_index] > get_value(idx, idy - 1, idz)) -
                    (data[global_index] < get_value(idx, idy - 1, idz));
        int sign4 = (data[global_index] > get_value(idx, idy + 1, idz)) -
                    (data[global_index] < get_value(idx, idy + 1, idz));
        int sign5 = (data[global_index] > get_value(idx, idy, idz - 1)) -
                    (data[global_index] < get_value(idx, idy, idz - 1));
        int sign6 = (data[global_index] > get_value(idx, idy, idz + 1)) -
                    (data[global_index] < get_value(idx, idy, idz + 1));
        return std::abs(sign1 + sign2 + sign3 + sign4 + sign5 + sign6);
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

    for (int i = 1; i <= max_padding; ++i) {
      // dir1 index
      if (is_valid_index(idx - i, idy) && is_valid_index(idx + i, idy) &&
          map_value(idx - i, idy) == 2 && map_value(idx + i, idy) == 2) {
        xpadding = i;
      } else {
        { break; }
      }
    }

    for (int i = 1; i <= max_padding; ++i) {
      // dir1 index
      if (is_valid_index(idx, idy - i) && is_valid_index(idx, idy + i) &&
          map_value(idx, idy - i) == 2 && map_value(idx, idy + i) == 2) {
        ypadding = i;
      } else {
        break;
      }
    }
    T local_abs_max = 0;
    auto get_local_value_range = [this, idx, idy, xpadding, ypadding](T &local_abs_max ) {
      T min = std::numeric_limits<T>::max();
      T max = -std::numeric_limits<T>::max();
      for (int i = idx - xpadding; i < idx + xpadding; ++i) {
        for (int j = idy - ypadding; j < idy + ypadding; ++j) {
          T value = get_value(i, j);
          if (value < min)
            min = value;
          if (value > max)
            max = value;
        }
      }
      local_abs_max = std::max(std::abs(min), std::abs(max));
      return (T)(max - min);
    };

    if (xpadding != 0 && ypadding != 0)
    { T local_value_range = get_local_value_range(local_abs_max);
      if (local_value_range > (T)(1e-5) && local_abs_max > (T)1e-5) return 1;
      else return 0;
    }
    else
      return 0;
  }

  bool pattern_match3d(std::vector<int> &cpmap, int idx, int idy, int idz,
                       int &xpadding, int &ypadding, int &zpadding,
                       int max_padding = 7) {
    xpadding = 0;
    ypadding = 0;
    zpadding = 0;
    int global_index = idx + idy * global_dimensions[0] +
                       idz * global_dimensions[0] * global_dimensions[1];

    for (int i = 1; i <= max_padding; ++i) {
      if (is_valid_index(idx - i, idy, idz) &&
          is_valid_index(idx + i, idy, idz) &&
          map_value(idx - i, idy, idz) >= 2 &&
          map_value(idx + i, idy, idz) >= 2) {
        xpadding = i;
      } else {
        break;
      }
    }

    for (int i = 1; i <= max_padding; ++i) {
      if (is_valid_index(idx, idy - i, idz) &&
          is_valid_index(idx, idy + i, idz) &&
          map_value(idx, idy - i, idz) >= 2 &&
          map_value(idx, idy + i, idz) >= 2) {
        ypadding = i;
      } else {
        break;
      }
    }

    for (int i = 1; i <= max_padding; ++i) {
      if (is_valid_index(idx, idy, idz - i) &&
          is_valid_index(idx, idy, idz + i) &&
          map_value(idx, idy, idz - i) >= 2 &&
          map_value(idx, idy, idz + i) >= 2) {
        zpadding = i;
      } else {
        break;
      }
    }

    if (idx == 148 && idy == 88 && idz == 64) {
      std::cout << "global index = " << global_index << std::endl;
      std::cout << "map_value = " << map_value(idx, idy, idz) << std::endl;
      std::cout << "xpadding = " << xpadding << std::endl;
      std::cout << "ypadding = " << ypadding << std::endl;
      std::cout << "zpadding = " << zpadding << std::endl;
    }
    T local_abs_max = 0;
    auto get_local_value_range = [this, idx, idy, idz, xpadding, ypadding,
                               zpadding](T &local_abs_max) {
      T min = std::numeric_limits<T>::max();
      T max = -std::numeric_limits<T>::max();
      for (int i = idx - xpadding; i < idx + xpadding; ++i) {
        for (int j = idy - ypadding; j < idy + ypadding; ++j) {
          for (int k = idz - zpadding; k < idz + zpadding; ++k) {
            T value = get_value(i, j, k);
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

    if ((xpadding != 0 && ypadding != 0) || (xpadding != 0 && zpadding != 0) ||
        (ypadding != 0 && zpadding != 0)) {
      T local_value_range = get_local_value_range(local_abs_max);
      if (local_value_range > (T)(1e-5) && local_abs_max > (T)1e-5)
        return true;
      else
        return false;
    } else
      return false;
  }

  // TODO: strict pattern matching
  bool pattern_match2d_strict(std::vector<int> &cpmap, int idx, int idy,
                              int &xpadding, int &ypadding,
                              int max_padding = 5) {
    xpadding = 0;
    ypadding = 0;

    for (int i = 1; i <= max_padding; ++i) {
      if (is_valid_index(idx - i, idy) && is_valid_index(idx + i, idy) &&
          (map_value(idx - i, idy) == 2) && (map_value(idx + i, idy) == 2)) {
        xpadding = i;
      } else {
        break;
      }
    }

    for (int i = 1; i <= max_padding; ++i) {
      if (is_valid_index(idx, idy - i) && is_valid_index(idx, idy + i) &&
          (map_value(idx, idy - i) == 2) && (map_value(idx, idy + i) == 2)) {
        ypadding = i;
      } else {
        break;
      }
    }

    if (xpadding != 0 && ypadding != 0) {
      ypadding = std::min(xpadding, ypadding);
      xpadding = ypadding;
    } else {
      return 0;
    }
    // check other places
    int global_index = idx + idy * global_dimensions[0];

    if (global_index == 46244) {
      for (int i = 1; i <= xpadding; i++) {
        for (int j = 1; j <= xpadding; j++) {
          std::cout << idx - i << " " << idy - j << " " << idx + i << " "
                    << idy + j << std::endl;
          std::cout << map_value(idx - i, idy - j) << " "
                    << map_value(idx - i, idy + j) << " "
                    << map_value(idx + i, idy - j) << " "
                    << map_value(idx + i, idy + j) << std::endl;
        }
      }
    }

    if (xpadding == 0 && ypadding == 0)
      return 0;
    else
      return 1;
  }
};
