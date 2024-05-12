#ifndef SZ3_POSTPROCESS
#define SZ3_POSTPROCESS
#include "SZ3/def.hpp"
#include<vector>
#include<array>
#include<iostream> 
#include "SZ3/quantization/Quantizer.hpp"



namespace SZ {



template <class T, uint N>
class Postprocessor {
  template <uint NN = N>
  typename std::enable_if<NN == 1, double>::type block_post_process(
      T* data, std::array<size_t, N> begin, std::array<size_t, N> end,
       const std::string& interp_func,
      const int direction, size_t stride = 1)
  {
    return 0;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 2, double>::type block_post_process(
      T* data, std::array<size_t, N> begin, std::array<size_t, N> end,
       const std::string& interp_func,
      const int direction, size_t stride = 1)
  {
    return 0;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 3, double>::type block_post_process(
      T* data, std::array<size_t, N> begin, std::array<size_t, N> end,
       const std::string& interp_func,
      const int direction, size_t stride = 1)
  {
    size_t stride2x = stride * 2;
    auto default_eb = quantizer.get_eb();
    const std::array<int, N> dims = dimension_sequences[direction];

    compensation_3d(
        data, aux_quant_inds.data(), begin.data(), end.data(), dims,
        dimension_offsets, 0, stride, 1, 2, stride2x, stride2x,
        quantizer.get_eb(), quantizer.get_radius());

    compensation_3d(
        data, aux_quant_inds.data(), begin.data(), end.data(), dims,
        dimension_offsets, 1, stride, 0, 2, stride, stride2x,
        quantizer.get_eb(), quantizer.get_radius());

    compensation_3d(
        data, aux_quant_inds.data(), begin.data(), end.data(), dims,
        dimension_offsets, 2, stride, 0, 1, stride, stride, quantizer.get_eb(),
        quantizer.get_radius());

    // std::cout<<"quantizer radius = " << quantizer.get_radius() << std::endl;
    // writefile("first_compressed.dat", data, num_elements);

    return 0;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 4, double>::type block_post_process(
      T* data, std::array<size_t, N> begin, std::array<size_t, N> end,
       const std::string& interp_func,
      const int direction, size_t stride = 1)
  {
    return 0;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 1, uchar*>::type error_sample_compress(
      const T* orig_data_ptr, T* data, size_t& compressed_error_size,
      int direction_sequence_id,
      const std::vector<std::array<int, 1>>& dimension_sequences,
      std::array<size_t, 1>& global_dimensions,
      std::array<size_t, 1>& dimension_offsets, uchar*& buffer_pos)
  {
    return nullptr;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 2, uchar*>::type error_sample_compress(
      const T* orig_data_ptr, T* data, size_t& compressed_error_size,
      int direction_sequence_id,
      const std::vector<std::array<int, 2>>& dimension_sequences,
      std::array<size_t, 2>& global_dimensions,
      std::array<size_t, 2>& dimension_offsets, uchar*& buffer_pos)
  {
    return nullptr;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 3, uchar*>::type error_sample_compress(
      const T* orig_data_ptr, T* data, size_t& compressed_error_size,
      int direction_sequence_id,
      const std::vector<std::array<int, 3>>& dimension_sequences,
      std::array<size_t, 3>& global_dimensions,
      std::array<size_t, 3>& dimension_offsets, uchar*& buffer_pos)
  {
    std::cout << "compile for 3D cases\n";
    // write the downsampled and decompressed data for the last level
    auto error_compressor = SZ::ErrorCompressor<T, 3>(
        num_elements, global_dimensions.data(), dimension_offsets.data(),
        direction_sequence_id);
    compressed_error_size = 0;
    auto dims = dimension_sequences[direction_sequence_id];
    int interp_direction = dims[2];
    int interp_dir_stride = 2;
    int plane_dir1 = dims[0];
    int plane_dir2 = dims[1];
    int plane_dir1_stride = 1;
    int plane_dir2_stride = 1;
    int plane_sample_stride = 4;

    uchar* error_compressed_data = error_compressor.error_compensation_3d(
        orig_data_ptr, data, compressed_error_size, interp_direction,
        interp_dir_stride, plane_dir1, plane_dir1_stride, plane_dir2,
        plane_dir2_stride, plane_sample_stride);
    // std::cout << error_compressed_data << std::endl;
    // append the compressed error data to the buffer
    // std::cout<<"complete writing error data\n";
    // write(compressed_error_size, buffer_pos);
    // write(error_compressed_data, compressed_error_size, buffer_pos);
    // delete[] error_compressed_data;
    // return nullptr;
    return error_compressed_data;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 4, uchar*>::type error_sample_compress(
      const T* orig_data_ptr, T* data, size_t& compressed_error_size,
      int direction_sequence_id,
      const std::vector<std::array<int, 4>>& dimension_sequences,
      std::array<size_t, 4>& global_dimensions,
      std::array<size_t, 4>& dimension_offsets, uchar*& buffer_pos)
  {
    return nullptr;
  }
};

}  // namespace SZ

#endif  // SZ3_POSTPROCESS