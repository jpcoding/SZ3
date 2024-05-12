#ifndef SZ_COMPENSATE_2D_HPP
#define SZ_COMPENSATE_2D_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <iostream>
#include <numeric>
#include <vector>

namespace SZ {
template <typename T>
class Compensate2D {
 public:
  Compensate2D(
      T* data, int* quant_inds, const size_t* begin, const size_t* end,
      const std::array<int, 3>& dims,
      const std::array<size_t, 3>& dimension_offsets,
      int data_interp_direction, int data_interp_stride, 
      int data_plane_dir0, int data_plane_dir1, 
      int data_plane_dir0_stride, int data_plane_dir1_stride,
       const int quantizer_radius) :
      quant_inds(quant_inds),
      begin(begin),
      end(end),
      dims(dims),
      dimension_offsets(dimension_offsets),
      data_interp_direction(data_interp_direction),
      data_interp_stride(data_interp_stride),
      data_plane_dir0(data_plane_dir0),
      data_plane_dir1(data_plane_dir1),
      data_dir0_stride(data_plane_dir0_stride),
      data_dir1_stride(data_plane_dir1_stride),
      quantizer_radius(quantizer_radius)
  {
    boundary_map = build_boundary_map(
        quant_inds, begin, end, dims, data_plane_dir0, data_plane_dir1,
        data_plane_dir0_stride, data_plane_dir1_stride, quantizer_radius);
  }

  void process(T max_compensation, int max_extend = 5)
  {
    compensation3d_data_2d_plane(max_compensation, max_extend);
  }

  std::vector<char> build_boundary_map(
      int* quant_inds, const size_t* begin, const size_t* end,
      const std::array<int, 3>& dims, int plane_dir0, int plane_dir1,
      int plane_dir0_stride, int plane_dir1_stride, const int quantizer_radius)
  {
    size_t plane_dim0 =
        (end[dims[plane_dir0]] - begin[dims[plane_dir0]]) / plane_dir0_stride +
        1;  // stride = plane dim1
    size_t plane_dim1 =
        (end[dims[plane_dir1]] - begin[dims[plane_dir1]]) / plane_dir1_stride +
        1;  // stride = 1
    std::vector<char> boundary_map(plane_dim0 * plane_dim1, 0);
    boundary_x_min = 0;
    boundary_x_max = plane_dim0 - 1;
    boundary_y_min = 0;
    boundary_y_max = plane_dim1 - 1;
    for (size_t i = 1; i < plane_dim0 - 1; i++) {
      size_t idx_i = begin[dims[plane_dir0]] + i * plane_dir0_stride;

      for (size_t j = 1; j < plane_dim1 - 1; j++) {
        size_t idx_j = idx_i + j * plane_dir1_stride;
        /*
        check the neighbors of X
        # A #
        B X C
        # D #
        */
        int X = quant_inds[idx_j] == 0 ? quantizer_radius : quant_inds[idx_j];
        int A = quant_inds[idx_j - plane_dir1_stride] != 0
                    ? quant_inds[idx_j - plane_dir1_stride]
                    : quantizer_radius;
        int B = quant_inds[idx_j - plane_dir0_stride] != 0
                    ? quant_inds[idx_j - plane_dir0_stride]
                    : quantizer_radius;
        int C = quant_inds[idx_j + plane_dir1_stride] != 0
                    ? quant_inds[idx_j + plane_dir1_stride]
                    : quantizer_radius;
        int D = quant_inds[idx_j + plane_dir0_stride] != 0
                    ? quant_inds[idx_j + plane_dir0_stride]
                    : quantizer_radius;
        if ((X != A) || (X != B) || (X != C) || (X != D)) {
          boundary_map[i * plane_dim1 + j] = 1;
          // update boundary
          boundary_x_min = std::min(boundary_x_min, (int)i);
          boundary_x_max = std::max(boundary_x_max, (int)i);
          boundary_y_min = std::min(boundary_y_min, (int)j);
          boundary_y_max = std::max(boundary_y_max, (int)j);
        }
      }
    }
    for (size_t i = 1; i < plane_dim0 - 1; i++) {
      for (size_t j = 1; j < plane_dim1 - 1; j++) {
        if (quant_inds[begin[dims[plane_dir0]] + i * plane_dir0_stride +
                       j * plane_dir1_stride] != 0) {
          int left = quant_inds[begin[dims[plane_dir0]] + (i - 1) * plane_dir0_stride +
                                j * plane_dir1_stride];
          int right = quant_inds[begin[dims[plane_dir0]] + (i + 1) * plane_dir0_stride +
                                 j * plane_dir1_stride];
          int up = quant_inds[begin[dims[plane_dir0]] + i * plane_dir0_stride +
                              (j - 1) * plane_dir1_stride];
          int down = quant_inds[begin[dims[plane_dir0]] + i * plane_dir0_stride +
                                (j + 1) * plane_dir1_stride];
          if (left == 0 && right == 0 && up == 0 && down == 0) {
            boundary_map[i * plane_dim1 + j] = 0;
            boundary_map[(i - 1) * plane_dim1 + j] = 0;
            boundary_map[(i + 1) * plane_dim1 + j] = 0;
            boundary_map[i * plane_dim1 + j - 1] = 0;
            boundary_map[i * plane_dim1 + j + 1] = 0;
          }
        }
      }
    }
    return boundary_map;
  }

  std::array<int, 4> get_distance_to_boundary(
      const int x, const int y, int max_extend = 5)
  {
    int dim0 = dimx;
    int dim1 = dimy;
    auto boundary = boundary_map.data();
    // ignore 4 borders
    if (x == 0 || x == dim0 - 1 || y == 0 || y == dim1 - 1) {
      return {max_extend, max_extend, max_extend, max_extend};
    }
    // left right up down
    if (boundary[x * dim1 + y] == 0) { return {0, 0, 0, 0}; }
    std::array<int, 4> distance = {0, 0, 0, 0};  // left right up down
    char current_status = boundary[x * dim1 + y];
    int dx, dy;
    // left
    dy = y;
    dx = x - 1;
    while (boundary[dx * dim1 + dy] == current_status) {
      dx--;
      distance[0]++;
      if (distance[0] == max_extend) { break; }
    }
    // right
    dy = y;
    dx = x + 1;
    while (boundary[dx * dim1 + dy] == current_status) {
      dx++;
      distance[1]++;
      if (distance[1] == max_extend) { break; }
    }
    // up
    dx = x;
    dy = y - 1;
    while (boundary[dx * dim1 + dy] == current_status) {
      dy--;
      distance[2]++;
      if (distance[2] == max_extend) { break; }
    }
    // down
    dx = x;
    dy = y + 1;
    while (boundary[dx * dim1 + dy] == current_status) {
      dy++;
      distance[3]++;
      if (distance[3] == max_extend) { break; }
    }
    return distance;
  }

  std::array<char, 4> get_compensation_direction(int x, int y)
  {
    // this will seach the quantization plane,
    // so, use the original data offset to search
    // size_t offset_x = begin[plane_dir0];
    // size_t offset_y = begin[plane_dir1];
    int stride_x = data_dir0_stride;
    int stride_y = data_dir1_stride;
    size_t offset_x = dimension_offsets[data_plane_dir0];
    size_t offset_y = dimension_offsets[data_plane_dir1];

    std::array<char, 4> compensation_direction = {
        0, 0, 0, 0};  // left right up down
    int current_quant =
        quant_inds[x * stride_x * offset_x + y * stride_y * offset_y];
    // move to left to check the compensation direction
    int tx = x;
    int ty = y - 1;
    int* quant_ptr =
        quant_inds + x * stride_x * offset_x + y * stride_y * offset_y;
    while (ty > 0) {
      if (*quant_ptr != current_quant) {
        compensation_direction[0] = get_sign(current_quant - *quant_ptr);
        break;
      }
      ty--;
      quant_ptr -= stride_y * offset_y;
    }
    // move to right to check the compensation direction
    tx = x;
    ty = y + 1;
    quant_ptr = quant_inds + x * stride_x * offset_x + y * stride_y * offset_y;
    while (ty < dimy) {
      if (*quant_ptr != current_quant) {
        compensation_direction[1] =
            get_sign(-current_quant + *quant_ptr);  // always look backward
        break;
      }
      ty++;
      quant_ptr += stride_y * offset_y;
    }
    // move to up to check the compensation direction
    tx = x - 1;
    ty = y;
    quant_ptr = quant_inds + x * stride_x * offset_x + y * stride_y * offset_y;
    while (tx > 0) {
      if (*quant_ptr != current_quant) {
        compensation_direction[2] = get_sign(current_quant - *quant_ptr);
        break;
      }
      tx--;
      quant_ptr -= stride_x * offset_x;
    }
    // move to down to check the compensation direction
    tx = x + 1;
    ty = y;
    quant_ptr = quant_inds + x * stride_x * offset_x + y * stride_y * offset_y;
    while (tx < dimx) {
      if (*quant_ptr != current_quant) {
        compensation_direction[3] =
            get_sign(-current_quant + *quant_ptr);  // always look backward
        break;
      }
      tx++;
      quant_ptr += stride_x * offset_x;
    }
    return compensation_direction;
  }

  std::array<int, 4> get_boundary_info(int x, int y, int max_extend = 5)
  {
    // this will search the built boundary map
    // so, use the offset of the plane
    int stride_x = dimy;
    int stride_y = 1;
    std::array<int, 4> boundary_info = {0, 0, 0, 0};  // left right up down
    char current_value = boundary_map[x * stride_x + y * stride_y];
    int tx, ty;
    tx = x;
    ty = y - 1;
    char* boundary_ptr = boundary_map.data() + x * stride_x + y * stride_y;
    while (*boundary_ptr == current_value) {
      ty--;
      boundary_ptr -= stride_y;
      if (ty < 0) break;
    }
    boundary_info[0] = std::min((y - ty - 1) / 2, max_extend);
    tx = x;
    ty = y + 1;
    boundary_ptr = boundary_map.data() + x * stride_x + y * stride_y;
    while (*boundary_ptr == current_value) {
      ty++;
      boundary_ptr += stride_y;
      if (ty >= dimy) break;
    }
    boundary_info[1] = std::min((ty - y - 1) / 2, max_extend);
    tx = x - 1;
    ty = y;
    boundary_ptr = boundary_map.data() + x * stride_x + y * stride_y;
    while (*boundary_ptr == current_value) {
      tx--;
      boundary_ptr -= stride_x;
      if (tx < 0) break;
    }
    boundary_info[2] = std::min((x - tx - 1) / 2, max_extend);
    tx = x + 1;
    ty = y;
    boundary_ptr = boundary_map.data() + x * stride_x + y * stride_y;
    while (*boundary_ptr == current_value) {
      tx++;
      boundary_ptr += stride_x;
      if (tx >= dimx) break;
    }
    boundary_info[3] = std::min((tx - x - 1) / 2, max_extend);
    return boundary_info;
  }


  std::array<int, 4> get_quant_change_distance(int x, int y, int max_extend = 5)
  { 
    std::array<int, 4> distance = {0, 0, 0, 0};  // left right up down
    int stride_x = data_dir0_stride;
    int stride_y = data_dir1_stride;
    size_t offset_x = dimension_offsets[data_plane_dir0];
    size_t offset_y = dimension_offsets[data_plane_dir1];
    int current_quant = quant_inds[x * stride_x * offset_x + y * stride_y * offset_y];
    // move to left to check the compensation direction
    int tx = x;
    int ty = y - 1;
    int* quant_ptr = quant_inds + x * stride_x * offset_x + y * stride_y * offset_y;
    while (ty > 0) {
      if (*quant_ptr != current_quant) {
        distance[0]++;
        if (distance[0] == max_extend) { break; }
      }
      ty--;
      quant_ptr -= stride_y * offset_y;
    }
    // move to right to check the compensation direction
    tx = x;
    ty = y + 1;
    quant_ptr = quant_inds + x * stride_x * offset_x + y * stride_y * offset_y;
    while (ty < dimy) {
      if (*quant_ptr != current_quant) {
        distance[1]++;
        if (distance[1] == max_extend) { break; }
      }
      ty++;
      quant_ptr += stride_y * offset_y;
    }
    // move to up to check the compensation direction
    tx = x - 1;
    ty = y;
    quant_ptr = quant_inds + x * stride_x * offset_x + y * stride_y * offset_y;
    while (tx > 0) {
      if (*quant_ptr != current_quant) {
        distance[2]++;
        if (distance[2] == max_extend) { break; }
      }
      tx--;
      quant_ptr -= stride_x * offset_x;
    }
    // move to down to check the compensation direction
    tx = x + 1;
    ty = y;
    quant_ptr = quant_inds + x * stride_x * offset_x + y * stride_y * offset_y;
    while (tx < dimx) {
      if (*quant_ptr != current_quant) {
        distance[3]++;
        if (distance[3] == max_extend) { break; }
      }
      tx++;
      quant_ptr += stride_x * offset_x;
    }
    return distance;
  }

  template <typename T_array>
  char argmin(std::array<T_array, 4> arr)
  {
    T_array min_val = arr[0];
    char min_idx = 0;
    for (int i = 1; i < 4; i++) {
      if (arr[i] < min_val) {
        min_val = arr[i];
        min_idx = i;
      }
    }
    return min_idx;
  }

  template <typename T_array>
  char argmax(std::array<T_array, 4> arr)
  {
    T_array max_val = arr[0];
    char max_idx = 0;
    for (int i = 1; i < 4; i++) {
      if (arr[i] > max_val) {
        max_val = arr[i];
        max_idx = i;
      }
    }
    return max_idx;
  }

  int compensation3d_data_2d_plane(double max_compensation, int max_extend = 5)
  {
    // this fucntion will change the data directly
    // construct a 2d plane for the compensation
    int distance_to_boundary = max_extend;
    std::vector<T> compensation_plane(
        dimx * dimy, 0);  // same with the input data type

    int num_planes = dims[data_interp_direction] / data_interp_stride;

    for (size_t k =
             (begin[dims[data_interp_direction]]
                  ? begin[dims[data_interp_direction]] + data_interp_stride * 2
                  : data_interp_stride);
         k <= end[dims[data_interp_direction]]; k += data_interp_stride * 2) {
      {
        size_t begin_offset =
            (begin[dims[data_plane_dir0]]
                 ? begin[dims[data_plane_dir0]] + data_dir0_stride
                 : 0) *
                dimension_offsets[dims[data_plane_dir0]] +
            (begin[dims[data_plane_dir1]]
                 ? begin[dims[data_plane_dir1]] + data_dir1_stride
                 : 0) *
                dimension_offsets[dims[data_plane_dir1]] +
            k * dimension_offsets[dims[data_interp_direction]];

        // clean the compensation plane
        std::fill(compensation_plane.begin(), compensation_plane.end(), 0);

        for (int i = max_extend; i < dimx - 1; i++) {
          for (int j = max_extend; j < dimy - 1; j++) {
            // get the distance to the boundary
            auto distance_to_boundary =
                get_distance_to_boundary(i, j, max_extend);
            if (*std::min_element(distance_to_boundary.begin(), distance_to_boundary.end()) >= max_extend) {
              // do nothing
              continue;
            }
            else if (boundary_map[i * dimy + j] == 1) {
              // get the compensation direction
              auto compensation_direction = get_compensation_direction(i, j);
              // get the boundary info
              auto boundary_info = get_boundary_info(i, j, max_extend);
              auto quant_change_distance = get_quant_change_distance(i, j, max_extend);
              // get the compensation value
              T compensation_value = 0;

              if (*std::max_element(boundary_info.begin(),boundary_info.end()) == 0) {
                auto direction = argmin(quant_change_distance);
                float sign = pow(-1, direction + 1)* compensation_direction[direction]; 
                T reduction = 0.9;
                compensation_plane[i * dimy + j] =
                    sign * reduction * max_compensation;
              }
              else {
                T reduction =
                    (max_extend - *std::min_element(boundary_info.begin(), boundary_info.end())) / max_extend;
                char dir_idx = argmin(quant_change_distance);
                char compensation_dir = compensation_direction[dir_idx];
                char opposite_dir = opposite_direction[dir_idx];
                float sign = pow(-1, dir_idx + 1) * compensation_direction[dir_idx];
                compensation_plane[i * dimy + j] = sign * reduction * max_compensation;
              }
            }
            else {
              // these points are not on the boundary
              // only compensate the cloest direction
              auto compensation_direction = get_compensation_direction(i, j);
              char dir_idx = argmin(distance_to_boundary);
              int distance = distance_to_boundary[dir_idx];
              char opposite_dir = opposite_direction[dir_idx];
              int opposite_distance = distance_to_boundary[opposite_dir];
              int bound_width_max = 4;


              int current_width = distance+opposite_distance;
              int width =
                  std::min(current_width, bound_width_max);
              float sign =
                  pow(-1, dir_idx + 1) * 1.0 * compensation_direction[dir_idx];
              T magnitude;
              if(compensation_direction[dir_idx] == compensation_direction[opposite_dir]){
                // if the compensation direction is the same with the opposite direction
                // then we need to reduce the compensation
                T mid = (current_width +2.0)/2.0;
                T magnitude = 1/mid/mid*(mid-distance-1)*(mid-distance-1);
              }
              else{
                current_width = width;
                T magnitude = 1.0 / (2.0+current_width) / (2.0+current_width)/(current_width+2 - distance -1)*(current_width+2 - distance -1); 
              }
              compensation_plane[i * dimy + j] = sign * magnitude * max_compensation;
            }
          }
        }
        // now update the data
        for (int i = max_extend; i < dimx - 1; i++) {
          for (int j = max_extend; j < dimy - 1; j++) {
            data
                [begin_offset +
                 i * data_dir0_stride * dimension_offsets[data_plane_dir0] +
                 j * data_dir1_stride * dimension_offsets[data_plane_dir1]] +=
                compensation_plane[i * dimy + j];
          }
        }
        // clean the compensation plane 
        std::fill(compensation_plane.begin(), compensation_plane.end(), 0);
      }
    }
    return 0;
  }

    template <typename T_val>
    char get_sign(T_val val)
    {
      return (val > 0) - (val < 0);
    }

   private:
    T* data;
    int* quant_inds;
    const size_t* begin;
    const size_t* end;
    const std::array<int, 3>& dims;
    const std::array<size_t, 3>& dimension_offsets;      // offset of each dimension
    int data_interp_direction;  // x, y, z
    int data_interp_stride;     // x, y, z stride
    int data_plane_dir0;        // x
    int data_plane_dir1;        // y
    int data_dir0_stride;       // x stride
    int data_dir1_stride;       // y stride

    // 2d plane information
    int dimx, dimy;              // x, y dimension
    int stride_x, stride_y = 1;  // x, y stride

    const int quantizer_radius;
    std::vector<char> boundary_map;
    int defualt_max_extend = 5;

    int boundary_x_min = 0;
    int boundary_x_max = 0;
    int boundary_y_min = 0;
    int boundary_y_max = 0;

    // direction utilities
    // index 0: left, 1: right, 2: up, 3: down
    std::array<char, 4> opposite_direction = {1, 0, 3, 2};  //
  };
}  // namespace SZ

#endif  // SZ_COMPENSATE_2D_HPP