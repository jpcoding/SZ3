#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include "SZ3/quantizer/Quantizer.hpp"

namespace SZ {

// function to compensate a line of data
// data: the data to be compensated
// quant_inds: the quantization indices, line start index
// data_stride: the stride of the data
// quant_stride: the stride of the quantization indices

template <typename T>
inline void boundary_compensation(
    T* data, int start, int end, int sign, int data_stride, double compensation_max)
{
    int size = end - start;
    if(size<5)
    {
        for(int i = 0; i < size; i++)
        {
            data[(i+start)*data_stride] = 1.0*(i+1)*(i+1)/(1.0*size*size)*sign*compensation_max;
        }
    }
    else {
        for(int i = 0; i < 5; i++)
        {
            data[(end-5+i)*data_stride] =(i+1)*(i+1)/(25.0)*sign*compensation_max;
        }
    }
}

template <typename T>
inline void reverse_array(T* data, int start, int end, int stride)
{
    for(int i = 0; i < (end-start)/2; i++)
    {
        T temp = data[(start+i)*stride];
        data[(start+i)*stride] = data[(end-i-1)*stride];
        data[(end-i-1)*stride] = temp;
    }
}




template <typename T>
inline void regular_compensation_linear(
    T* data, int start, int end, int sign, int stride, double compensation_max)
{
    int top = 0;
    int size = end - start;
    if(size%2==1)
    {
        int half_size = size/2;
        double slope = (top+1.0)/half_size;
        for (int i = 0; i < half_size; i++)
        {
            data[(i+start)*stride] = (-1+i*slope)*sign*compensation_max*-1;
            data[(start+size-i-1)*stride] = (-1+i*slope)*sign*compensation_max*-1;
        }
        data[(half_size+start)*stride] = top*sign*compensation_max*-1;
    }
    else
    {
        int half_size = size/2; 
        double slope = (top+1.0)/(half_size-1);
        for (int i = 0; i < half_size; i++)
        {
            data[(i+start)*stride] = (-1+i*slope)*sign*compensation_max*-1;
            data[(start+size-i-1)*stride] = (-1+i*slope)*sign*compensation_max*-1;
        }
    }
}


template <typename T>
void regular_compensation_cubic(
    T* data, int start, int end, int sign, int stride, double compensation_max)
{
    int size = end - start;
    if(size<10)
    {
        int top = 0;
        double mid = (size-1)*1.0/2;
        for (int i =0; i<size; i++)
        {
            data[(i+start)*stride] = ((-1-top)/(mid*mid)*(i-mid)*(i-mid)+top)*sign*compensation_max*-1;
        }
        if(size%2==0)
        {
            data[(start+size/2)*stride] = top*sign*compensation_max*-1;
        }
    }
    else {
        double mid = 5;
        for (int i =0; i<5; i++)
        {
            data[(i+start)*stride] = ((-1)/(mid*mid)*(i-mid)*(i-mid)+0)*sign*compensation_max*-1;
            data[(start+size-i-1)*stride] = ((-1)/(mid*mid)*(i-mid)*(i-mid)+0)*sign*compensation_max*-1;
        }
    }

}


// quant_inds: the quantization indices, the beginning of the line
template <typename T>
int compensate_line(
    T* compensation, int* quant_inds, int compensation_stride, int quant_stride,int length,
    double compensation_max)
{
  int defualt_extend_length = 5;
  int curr_sign, prev_sign, next_sign;
  int quant_base = 1 << 15;
  int curr_idx = 0, prev_idx, next_idx;  // line index, increased by 1
  int curr_quant_idx = 0, prev_quant_idx,
      next_quant_idx;  // quantization index, increased by quant_stride
  while (curr_idx < length) {
    prev_idx = curr_idx;
    curr_idx += 1;
    while (curr_idx < length &&
           quant_inds[curr_idx * quant_stride] == quant_inds[(curr_idx-1) * quant_stride]) {
      curr_idx += 1;
    }
    next_idx = curr_idx;
    // left boundary cases
    if (prev_idx == 0) {
      if (next_idx == length) { return 0; }
      prev_sign = 0;
      int extend_size = std::min(defualt_extend_length, next_idx);
      next_sign = (quant_inds[next_idx * quant_stride] -
                   quant_inds[(next_idx - 1) * quant_stride]);
      next_sign = (next_sign == 0) ? 0 : ((next_sign > 0) ? 1 : -1);
    
      boundary_compensation(
          compensation, next_idx - extend_size, curr_idx, next_sign, compensation_stride,
          compensation_max);
      continue; // go to next line 
    }
    else{
        prev_sign = (quant_inds[prev_idx * quant_stride] - quant_inds[(prev_idx - 1) * quant_stride]);
        prev_sign = (prev_sign == 0) ? 0 : ((prev_sign > 0) ? 1 : -1);
        
    }
    // normal cases 
    if (abs(prev_sign )>1)
    {
        prev_sign =0;
    }
    // right boundary cases
    if (next_idx == length)
    {
        next_sign = 0;
        int extend_size = std::min(defualt_extend_length, length-prev_idx);
        boundary_compensation(
          compensation, prev_idx, prev_idx + extend_size, prev_sign*-1, compensation_stride,
          compensation_max);
        // reverse order here with the stride 

        // std::reverse((compensation + prev_idx)*compensation_stride, compensation + (prev_idx + extend_size)*compensation_stride);

        reverse_array(compensation, prev_idx, prev_idx + extend_size, compensation_stride);

        continue; // next interval
    }
    else
    {   
        next_sign = (quant_inds[next_idx * quant_stride] - quant_inds[(next_idx - 1) * quant_stride]);
        next_sign = (next_sign == 0) ? 0 : ((next_sign > 0) ? 1 : -1);
    }
    // regular cases 
    if(abs(next_sign)>1)
    {
        next_sign = 0;
    }
    int interval_size = next_idx - prev_idx; 
    if(prev_sign == next_sign && prev_sign != 0 && next_sign != 0)
    {
        if(interval_size >1)
        {
            for( int i = 0; i < interval_size; i++)
            {
                compensation[(prev_idx + i)*compensation_stride] = (i*1.0/(interval_size -1)*2*compensation_max -compensation_max)*prev_sign;
            }
        }
        else
        {
            compensation[prev_idx*compensation_stride] = (-compensation_max)*prev_sign;
            
        }
    }
    else 
    {
        if (interval_size <5)
        {
            for( int i = 0; i < interval_size; i++)
            {
                compensation[(prev_idx + i)*compensation_stride] = (compensation_max)*next_sign;
            }
        }
        else
        {
        //     regular_compensation_linear(
        //         compensation, prev_idx, next_idx, next_sign,compensation_stride , compensation_max);
            regular_compensation_cubic(
                compensation, prev_idx, next_idx, next_sign,compensation_stride , compensation_max);
            
        }
    }
  }
  return 0;
}


// function to compensate a plane of data
// data: the data to be compensated
// quant_inds: the quantization indices, organized same with data
// dims: the dimensions of the plane, last dimension stride is 1
// length: the length of the plane
template <typename T>
int compensation_3d_(
    T* data, int* quant_inds, const size_t *dims, int last_interp_direction,
    double compensation_max)
{
    // 
    size_t plane_dim1, plane_dim2, interp_dim;
    size_t plane_stride1, plane_stride2, interp_stride;
    if(last_interp_direction ==0)
    {
        plane_dim1 = dims[1];
        plane_dim2 = dims[2];
        plane_stride1 = dims[2];
        plane_stride2 = 1;
        interp_stride = dims[1]*dims[2];

    }
    else if(last_interp_direction ==1)
    {
        plane_dim1 = dims[0];
        plane_dim2 = dims[2];
        plane_stride1 = dims[2]*dims[1];
        plane_stride2 = 1;
        interp_stride = dims[2];
    }
    else if(last_interp_direction ==2)
    {
        plane_dim1 = dims[0];
        plane_dim2 = dims[1];
        plane_stride1 = dims[2]*dims[1];
        plane_stride2 = dims[2];
        interp_stride = 1;
    }
    else {
        std::cout << "Error: last_interp_direction should be 0, 1, or 2" << std::endl;
        exit(0);
        return -1;
    }
    // number of slices 
    size_t num_slices = dims[last_interp_direction]/2;
    interp_dim = dims[last_interp_direction];
    
    int interp_slice_stride = 2;
    
    std::vector<T> compensation_plane_dir1(plane_dim1*plane_dim2, 0);
    std::vector<T> compensation_plane_dir2(plane_dim1*plane_dim2, 0);
    std::vector<int> quant_line1(plane_dim2, 0);
    std::vector<int> quant_line2(plane_dim1, 0);

    std::vector<T> compensation_line2(plane_dim1, 0);

    std::cout << "plane_dim1: " << plane_dim1 << " plane_dim2: " << plane_dim2 << std::endl;
    std::cout << "plane_stride1: " << plane_stride1 << " plane_stride2: " << plane_stride2 << std::endl;
    std::cout << "interp_stride: " << interp_stride << std::endl;
    std::cout << "interp_dim: " << interp_dim << std::endl;

    for(int i = 1; i <=dims[last_interp_direction]; i+=interp_slice_stride) // slices along the last interp dim
    {
        // data plane begin 
        T* data_plane_begin = data + i*interp_stride; 
        // get direction 1 compensation 
        for(int j = 0; j < plane_dim1; j++) // dir1 
        {
            std::fill(quant_line1.begin(), quant_line1.end(), 0);
           // beginning of the quant line
           int count_jump = 0;
           int quant_max = 0;
            for(int k = 0; k < plane_dim2; k++) // dir2 
            {
                int this_quant = quant_inds[i*interp_stride + j*plane_stride1 + k*plane_stride2];
                if(this_quant ==0)
                {
                    quant_line1[k] =0;
                }
                else
                {
                    quant_line1[k] = 32768-this_quant;
                }
                if(quant_line1[k] != 0)
                {
                    count_jump++;
                }
                quant_max = std::max(quant_max, std::abs(quant_line1[k]));
            }
            if(quant_max<=2 || count_jump<=2)
            {
                continue;
            }
            T* compensation_line_begin = compensation_plane_dir1.data() + j*plane_stride1;

                compensate_line(
                compensation_line_begin, quant_line1.data(), 
                1, 1,plane_dim2, compensation_max);
        }
        // get direction 2 compensation
        for(int j = 0; j < plane_dim2; j++) // dir2 
        {
            // beginning of the quant line
            // int* quant_line_begin = quant_inds + i*interp_stride + j*plane_stride2;
            int quant_max  = 0;
            int quant_jump = 0;
            for(int k = 0; k < plane_dim1; k++) // dir2 
            {
                size_t this_quant = quant_inds[i*interp_stride + j*plane_stride2+ k*plane_stride1];
                if(this_quant ==0)
                {
                    quant_line2[k] =0;
                }
                else
                {
                    quant_line2[k] = 32768-this_quant;
                }
                quant_max = std::max(quant_max, std::abs(quant_line2[k]));
                if(quant_line2[k] != 0)
                {
                    quant_jump++;
                }
            
                compensation_line2[k] = compensation_plane_dir2[k*plane_stride1 + j];
            }
            if(quant_max<=2 || quant_jump<=2)
            {
                continue;
            }

            // beginning of the compensarion line 
            compensate_line(
                compensation_line2.data(), quant_line2.data(), 
                1, 1,plane_dim1, compensation_max);
            // copy back 
            for(int k = 0; k < plane_dim1; k++) // dir2 
            {
                compensation_plane_dir2[k*plane_stride1 + j] = compensation_line2[k];
            }
            // compensation_line2.assign(plane_dim1, 0);
            std::fill(quant_line2.begin(), quant_line2.end(), 0);   

        }
        // // get the avg of the two compensation planes
        for (int i = 0; i< plane_dim1*plane_dim2; i++)
        {
            compensation_plane_dir1[i] = (compensation_plane_dir1[i] + compensation_plane_dir2[i])/2*0.8;
        }

        // apply the compensation to the data plane
        for(int j = 0; j < plane_dim1; j++) // dir1 
        {
            for(int k = 0; k < plane_dim2; k++) // dir2 
            {
                data_plane_begin[j*plane_stride1 + k*plane_stride2] -= compensation_plane_dir1[j*plane_dim2 + k];
            }
        }
        // clear the compensation planes
        std::fill(compensation_plane_dir1.begin(), compensation_plane_dir1.end(), 0);
        std::fill(compensation_plane_dir2.begin(), compensation_plane_dir2.end(), 0);
    }
    return 0;

}

template <typename T>
int compensation_3d(T*data, int* quant_inds, 
                const size_t* begin, 
                const size_t* end, 
                const int* dims, 
                const size_t* dimension_offsets, 
                int interp_direction,
                int interpolation_stride,   
                int plane_dir0, int plane_dir1, 
                int plane_dir0_stride, int plane_dir1_stride,
                 const double &compensation_max, 
                const int & quantizer_radius)
{

      size_t plane_dim0 = (end[dims[plane_dir0]] - begin[dims[plane_dir0]]) / plane_dir0_stride + 1; // stride = plane dim1
      size_t plane_dim1 = (end[dims[plane_dir1]] - begin[dims[plane_dir1]]) / plane_dir1_stride + 1; // stride = 1

      std::vector<T> plane_data(plane_dim1 * plane_dim0, 0);
      std::vector<T> plane_data_dirction2(plane_dim1 * plane_dim0, 0);

      std::vector<T> compensation_line_dir0(plane_dim1, 0);
      std::vector<T> compensation_line_dir1(plane_dim0, 0);

      std::vector<int> quant_line_dir0(plane_dim1, 0);
        std::vector<int> quant_line_dir1(plane_dim0, 0);


      int stride2x = interpolation_stride * 2;
      int *my_quant_inds = quant_inds;
      for (size_t i =
               (begin[dims[interp_direction]] ? begin[dims[interp_direction]] + interpolation_stride + interpolation_stride : interpolation_stride);
           i <= end[dims[interp_direction]]; i += stride2x) {
        // compensate plane size

        size_t begin_offset = (begin[dims[plane_dir0]] ? begin[dims[plane_dir0]] + plane_dir0_stride : 0) *
                                  dimension_offsets[dims[plane_dim0]] +
                              (begin[dims[plane_dir1]] ? begin[dims[plane_dir1]] + plane_dir1_stride : 0) *
                                  dimension_offsets[dims[plane_dir1]] +
                              i * dimension_offsets[dims[interp_direction]];
        // compensate the first direction
        for (int j = 0; j < plane_dim0; j++) {
          // clean quant inds
          int *current_quant_inds = my_quant_inds + begin_offset +
                                    j * dimension_offsets[dims[plane_dir0]] * plane_dir0_stride;
          int quant_max = 0;
          int quant_jump = 0;
          int jump_lower_bound = (int)(plane_dim1*0.005);
          for (int k = 0; k < plane_dim1; k++) {
            if (*current_quant_inds == quantizer_radius) {
              *current_quant_inds = 0;
            }
            else {
              *current_quant_inds =
                  *current_quant_inds - quantizer_radius;
            }
            quant_max = std::max(quant_max, std::abs(*current_quant_inds));
            if (*current_quant_inds != 0) {
              quant_jump++;
            }
            current_quant_inds += dimension_offsets[dims[plane_dir1]] * plane_dir1_stride;
          }
          if(quant_max<=1 )
          {
              continue;
          }
          compensate_line(
              plane_data.data() + j * plane_dim1,
              my_quant_inds + begin_offset +
                  j * dimension_offsets[dims[plane_dir0]] * plane_dir0_stride,
              1, dimension_offsets[dims[plane_dir1]] * plane_dir1_stride,
              plane_dim1, compensation_max);
            // compensate_line(
            //     compensation_line_dir0.data(), quant_line_dir0.data(), 
            //     1, 1,plane_dim1, compensation_max);
            // for (int k = 0; k < plane_dim1; k++) {
            //     plane_data[j * plane_dim1 + k] = compensation_line_dir0[k];
            // }
            // std::fill(compensation_line_dir0.begin(), compensation_line_dir0.end(), 0);
        }

        for (size_t j = 0; j < plane_dim1; j++) {

            int *current_quant_inds = my_quant_inds + begin_offset +
                                        j * dimension_offsets[dims[plane_dir1]] * plane_dir1_stride;
            int quant_max = 0;
            int quant_jump = 0;
            int jump_lower_bound = (int)(plane_dim0*0.005);
            for (int k = 0; k < plane_dim0; k++) {
              if (*current_quant_inds != 0) {
                quant_jump++;
              }
              quant_max = std::max(quant_max, std::abs(*current_quant_inds));
              current_quant_inds += dimension_offsets[dims[plane_dir0]] * plane_dir0_stride;
            }
            if(quant_max<=1 )
            {
                continue;
            }
          compensate_line(
              plane_data_dirction2.data() + j,
              my_quant_inds + begin_offset +
                  j * dimension_offsets[dims[plane_dir1]] * plane_dir1_stride,
              plane_dim1, dimension_offsets[dims[plane_dir0]] * plane_dir0_stride,
              plane_dim0, compensation_max);
            // compensate_line(
            //     compensation_line_dir1.data(), quant_line_dir1.data(), 
            //     1, 1,plane_dim0, compensation_max);
            // for (int k = 0; k < plane_dim0; k++) {
            //     plane_data_dirction2[k * plane_dim1 + j] = compensation_line_dir1[k];
            // }
            // std::fill(compensation_line_dir1.begin(), compensation_line_dir1.end(), 0);
        }
        // // avg 
        for (int j = 0; j < plane_dim1*plane_dim0; j++)
        {
          plane_data[j] = (plane_data[j] + plane_data_dirction2[j])/2*0.8;
        }
        // add plane compensation to the data
        for (size_t j = 0; j < plane_dim0; j++) {
          for (size_t k = 0; k < plane_dim1; k++) {
            data[begin_offset + j * dimension_offsets[dims[plane_dir0]] * plane_dir0_stride +
                 k * dimension_offsets[dims[plane_dir1]] * plane_dir1_stride] +=
                plane_data[j * plane_dim1 + k];
          }
        }
        std::fill(plane_data.begin(), plane_data.end(), 0);
        std::fill(
            plane_data_dirction2.begin(), plane_data_dirction2.end(), 0);
      }

    return 0; 

}


}  // namespace SZ3