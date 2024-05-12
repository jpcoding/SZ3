#ifndef SZ_COMPENSATE_ZERO_QUANT_HPP
#define SZ_COMPENSATE_ZERO_QUANT_HPP

#include <stddef.h> 
#include <array>
#include <cmath>
#include <vector> 

namespace SZ {

template<typename T>
double block_compensation_1d(T *data, int* ordered_quant, size_t begin, size_t end, size_t stride)
{
    // if a index is zero, check the four interpolators 
    size_t n = (end - begin) / stride + 1;
    if (n <= 1) { return 0; }
    double predict_error = 0;
    int quant_compensation = 0;
    size_t stride3x = 3 * stride;
    size_t stride5x = 5 * stride;
    T *d;
    size_t i;
    int radius = 1<<15;
    d = data + begin + stride;
    for (i = 3; i + 3 < n; i += 2) {
        d = data + begin + i * stride;
        int quant = ordered_quant[d-data];
        if (quant != radius) {
            continue;
        }
        T left =  (*(d + stride) - 2* *(d - stride)+ *(d - stride3x));
        T right = (*(d - stride) - 2* *(d + stride)+ *(d + stride3x));
        if(left * right < 0) {
            *d = (*(d - stride) + *(d + stride)) / 2;
        }
    }
    return 0.0;
}

template <typename T>
std::vector<char> compensation_zero_quant_by_plane(const T* original_data,
                T*data, const double tol, const double compensation,int scan_block_size,
                int* quant_inds, 
                const size_t* begin, 
                const size_t* end, 
                const std::array<int,3> & dims,
                std::array<size_t,3> & dimension_offsets, 
                int interp_direction,
                int interpolation_stride,   
                int plane_dir0, int plane_dir1, 
                int plane_dir0_stride, int plane_dir1_stride,
                const double compensation_max, 
                const int zero_quant_index
                )
{
    // for each slice, search concentrates region and compensate, tag that region 
          size_t plane_dim0 = (end[dims[plane_dir0]] - begin[dims[plane_dir0]]) / plane_dir0_stride + 1; // stride = plane dim1
      size_t plane_dim1 = (end[dims[plane_dir1]] - begin[dims[plane_dir1]]) / plane_dir1_stride + 1; // stride = 1

      int num_blocks_dir0 = (plane_dim0 + scan_block_size - 1) / scan_block_size;
      int num_blocks_dir1 = (plane_dim1 + scan_block_size - 1) / scan_block_size;

      std::vector<char> compensation_tags;

      
      int block_size = scan_block_size;

      int stride2x = interpolation_stride * 2;
      int stride3x = interpolation_stride * 3;
      for (size_t k =
               (begin[dims[interp_direction]] ? begin[dims[interp_direction]] + interpolation_stride + interpolation_stride : interpolation_stride);
           k <= end[dims[interp_direction]]; k += stride2x) {
        // slice begin 

        size_t begin_offset = (begin[dims[plane_dir0]] ? begin[dims[plane_dir0]] + plane_dir0_stride : 0) *
                                  dimension_offsets[dims[plane_dir0]] +
                              (begin[dims[plane_dir1]] ? begin[dims[plane_dir1]] + plane_dir1_stride : 0) *
                                  dimension_offsets[dims[plane_dir1]] +
                              k * dimension_offsets[dims[interp_direction]];
        size_t idx, idy, idxx, idyy;

        // on each slice, use stride*block size to scan
        for (int i = 0; i < plane_dim0; i += block_size) {
            idx = begin_offset + i * plane_dir0_stride * dimension_offsets[dims[plane_dir0]]; // block begin i 
            for(int j = 0; j < plane_dim1; j += block_size) {
                idy = idx + j * plane_dir1_stride * dimension_offsets[dims[plane_dir1]]; // block begin j
                // scan block
                bool compensate_zero_quant_block = true;
                int count_zero_quant_block = 0;
                bool error_direction = false; // 0: dir0, 1: dir1
                std::vector<bool> quant_zero_block(block_size*block_size, false);
                std::vector<char> error_directions(block_size*block_size, 0);
                // scan for quant info first 
                for(int ii = 0; ii < block_size; ii++) {
                    if(ii + i >= plane_dim0) {
                        break;
                    }
                    idxx = idy + ii * plane_dir0_stride * dimension_offsets[dims[plane_dir0]]; // block begin ii
                    int current_quant = quant_inds[idxx]; 
                    for(int jj = 0; jj < block_size; jj++) {
                        if(jj + j >= plane_dim1) {
                            break;
                        }
                        idyy = idxx + jj * plane_dir1_stride * dimension_offsets[dims[plane_dir1]]; // block begin jj
                        // check if the index is zero
                        if(quant_inds[idyy] == current_quant) {
                            count_zero_quant_block++;
                        }
                    }
                }

                if(count_zero_quant_block < block_size * block_size) {
                    compensate_zero_quant_block = false;
                }
                // check error direction
                int count_positive = 0;
                int count_negative = 0;
                if(compensate_zero_quant_block) {
                    // compensate the block
                    for(int ii = 0; ii < block_size; ii++) {
                        if(ii + i >= plane_dim0) {
                            break;
                        }
                        idxx = idy + ii * plane_dir0_stride * dimension_offsets[dims[plane_dir0]]; // block begin ii
                        for(int jj = 0; jj < block_size; jj++) {
                            if(jj + j >= plane_dim1) {
                                break;
                            }
                            idyy = idxx + jj * plane_dir1_stride * dimension_offsets[dims[plane_dir1]]; // block begin jj
                            // check error direction 
                            if(original_data[idyy]  - data[idyy] >=tol) {
                                count_positive++;
                            } else if(original_data[idyy] -  data[idyy] <= -tol) {
                                count_negative++;
                            }
                        }
                    }
                    if(count_positive == block_size * block_size || count_negative == block_size * block_size) {
                        char compensation_direction = count_positive == block_size * block_size ? 1 : -1;
                        compensation_tags.push_back(compensation_direction);
                        for(int ii = 0; ii < block_size; ii++) {
                            if(ii + i >= plane_dim0) {
                                break;
                            }
                            idxx = idy + ii * plane_dir0_stride * dimension_offsets[dims[plane_dir0]]; // block begin ii
                            for(int jj = 0; jj < block_size; jj++) {
                                if(jj + j >= plane_dim1) {
                                    break;
                                }
                                idyy = idxx + jj * plane_dir1_stride * dimension_offsets[dims[plane_dir1]]; // block begin jj
                                // compensate the block
                                data[idyy] += compensation_direction * compensation;
                            }
                        }
                    } 
                }
            }
        }
                 
}
return compensation_tags; 
}




}

#endif
