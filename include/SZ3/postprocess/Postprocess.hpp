#ifndef SZ3_POSTPROCESS
#define SZ3_POSTPROCESS

#include <array>
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/postprocess/error_compression.hpp"

namespace SZ {
template <class T, uint N>
class SZPostprocessor {
 public:
  SZPostprocessor(){}

  void post_process(T* processed_data, Config conf) {
    std::cout << "abs eb = " << conf.absErrorBound << std::endl;
    std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
    blocksize = conf.interpBlockSize;
    interpolator_id = conf.interpAlgo;
    std::cout << "interpolator_id = " << interpolator_id << std::endl;
    direction_sequence_id = conf.interpDirection;
    std::cout << "direction sequence id = " << direction_sequence_id << std::endl;
    aux_quant_inds = conf.PASS_DATA.aux_quant_inds_ptr;
    direction_sequence_id = conf.interpDirection;
    quantizer = SZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2);
    quantizer.set_eb(conf.absErrorBound);
    std::cout << "at the begin quantizer.get_eb() = " << quantizer.get_eb() << std::endl;

    int interpolator_id = conf.interpAlgo;
    // double eb_reduction_factor;

    init();

    int start_level = 1; 
    int end_level = 1; 


    for (uint level = start_level;
         level > 0 && level <= end_level; level--) {

        size_t stride = 1U << (level - 1);
        quantizer.set_eb(level_ebs[level - 1]);
        std::cout << "current quantizer.get_eb() = " << quantizer.get_eb() << std::endl;
        // change blocksize to the whole dataset size
        int blocksize_copy = blocksize;
        blocksize = *std::max_element(global_dimensions.begin(), global_dimensions.end())+1;
        std::cout << "blocksize = " << blocksize << std::endl;
        auto inter_block_range_post =
          std::make_shared<SZ::multi_dimensional_range<T, N>>(
              processed_data, std::begin(global_dimensions), std::end(global_dimensions),
              blocksize * stride, 0);
        auto inter_begin_post = inter_block_range_post->begin();
        auto inter_end_post = inter_block_range_post->end();
        for (auto block = inter_begin_post; block != inter_end_post; ++block) {
          auto end_idx = block.get_global_index();
          for (int i = 0; i < N; i++) {
            end_idx[i] += blocksize * stride;
            if (end_idx[i] > global_dimensions[i] - 1) {
              end_idx[i] = global_dimensions[i] - 1;
            }
          }
          block_smoothing(
              processed_data, block.get_global_index(), end_idx,
              interpolators[interpolator_id], direction_sequence_id, stride);
        }
        blocksize = blocksize_copy;
      }   


  }



  private:
  void init()
  {
    assert(
        blocksize % 2 == 0 &&
        "Interpolation block size should be even numbers");
    num_elements = 1;
    interpolation_level = -1;
    for (int i = 0; i < N; i++) {
      if (interpolation_level < ceil(log2(global_dimensions[i]))) {
        interpolation_level = (uint)ceil(log2(global_dimensions[i]));
      }
      num_elements *= global_dimensions[i];
    }

    dimension_offsets[N - 1] = 1;
    for (int i = N - 2; i >= 0; i--) {
      dimension_offsets[i] =
          dimension_offsets[i + 1] * global_dimensions[i + 1];
      // std::cout << "dimension_offsets[" << i << "] = " << dimension_offsets[i]
      //           << std::endl;
    }

    dimension_sequences = std::vector<std::array<int, N>>();
    auto sequence = std::array<int, N>();
    for (int i = 0; i < N; i++) { sequence[i] = i; }
    do {
      dimension_sequences.push_back(sequence);
    } while (std::next_permutation(sequence.begin(), sequence.end()));

    level_ebs.resize(interpolation_level);
    level_ebs[0] = quantizer.get_eb(); 
    std::cout << "quantizer.get_eb() = " << quantizer.get_eb() << std::endl; 
    for (int i = 1; i < interpolation_level; i++) {
      level_ebs[i] = level_ebs[i - 1] / eb_factors[1];
      std::cout<<"level_ebs["<<i+1<<"] = "<<level_ebs[i]<<std::endl;
    }
  }



    template <uint NN = N>
    typename std::enable_if<NN == 1, double>::type block_smoothing(
        T* data, std::array<size_t, N> begin, std::array<size_t, N> end,
        const std::string& interp_func, const int direction, size_t stride = 1)
    {
      return 0;
    }

    template <uint NN = N>
    typename std::enable_if<NN == 2, double>::type block_smoothing(
        T* data, std::array<size_t, N> begin, std::array<size_t, N> end,
        const std::string& interp_func, const int direction, size_t stride = 1)
    {
      return 0;
    }

    template <uint NN = N>
    typename std::enable_if<NN == 3, double>::type block_smoothing(
        T* data, std::array<size_t, N> begin, std::array<size_t, N> end,
        const std::string& interp_func, const int direction, size_t stride = 1)
    {
      size_t stride2x = stride * 2;
      auto default_eb = quantizer.get_eb();
      const std::array<int, N> dims = dimension_sequences[direction];

      compensation_3d(
          data, (*aux_quant_inds).data(), begin.data(), end.data(), dims,
          dimension_offsets, 0, stride, 1, 2, stride2x, stride2x,
          quantizer.get_eb(), quantizer.get_radius());

      compensation_3d(
          data, (*aux_quant_inds).data(), begin.data(), end.data(), dims,
          dimension_offsets, 1, stride, 0, 2, stride, stride2x,
          quantizer.get_eb(), quantizer.get_radius());

      compensation_3d(
          data, (*aux_quant_inds).data(), begin.data(), end.data(), dims,
          dimension_offsets, 2, stride, 0, 1, stride, stride,
          quantizer.get_eb(), quantizer.get_radius());

      return 0;
    }

    template <uint NN = N>
    typename std::enable_if<NN == 4, double>::type block_smoothing(
        T* data, std::array<size_t, N> begin, std::array<size_t, N> end,
        const std::string& interp_func, const int direction, size_t stride = 1)
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

 private:
  int interpolation_level = -1;
  uint blocksize;
  int interpolator_id;
  double eb_ratio = 0.5;
  std::vector<std::string> interpolators = {"linear", "cubic"};
  std::array<double, N> eb_factors = {pow(sqrt(1.5),N), pow(1.2808688457449497,N)};
  std::vector<int> quant_inds;
  size_t quant_index = 0;  // for decompress
  double max_error;
  SZ::LinearQuantizer<T> quantizer;
  size_t num_elements;
  std::array<size_t, N> global_dimensions;
  std::array<size_t, N> dimension_offsets;
  std::vector<std::array<int, N>> dimension_sequences;
  int direction_sequence_id;


  double linear_interp_eb_factor = sqrt(1.5);
  double cubic_interp_eb_factor = 1.2808688457449497;
  std::vector<double> level_ebs; // error bound for each level
  std::shared_ptr<std::vector<int>> aux_quant_inds = nullptr;
};

}  // namespace SZ

#endif  // SZ3_POSTPROCESS