#ifndef _SZ_TINY_COMPRESSSOR_HPP
#define _SZ_TINY_COMPRESSSOR_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <future>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "SZ3/compressor/SZInterpCompressorHelp.hpp"
#include "SZ3/utils/Statistic.hpp"
// #include "SZ3/compressor/SZInterpolation_postprocess.hpp"
#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/Accumulator.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"

namespace SZ {
template <class T, uint N, class Quantizer, class Encoder, class Lossless>
class SZTinyInterpolationCompressor {
 public:
  SZTinyInterpolationCompressor(
      Quantizer quantizer, Encoder encoder, Lossless lossless) :
      quantizer(quantizer), encoder(encoder), lossless(lossless)
  {
    static_assert(
        std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
        "must implement the quatizer interface");
    static_assert(
        std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
        "must implement the encoder interface");
    static_assert(
        std::is_base_of<concepts::LosslessInterface, Lossless>::value,
        "must implement the lossless interface");
  }

  T *decompress(uchar const *cmpData, const size_t &cmpSize, size_t num)
  {
    T *dec_data = new T[num];
    return decompress(cmpData, cmpSize, dec_data);
  }

  T *decompress(uchar const *cmpData, const size_t &cmpSize, T *decData)
  {
    return nullptr;
  }

  void init(const Config &conf, T *data )
  {
    std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
    blocksize = conf.interpBlockSize;
    interpolator_id = conf.interpAlgo;
    direction_sequence_id = conf.interpDirection;
    working_data = data;
    init();
  }

  // add the following vars to conf 
  // for each level, add direction, interp method, error bound factor
  
  void try_compress( int start_level, int end_level, const Config &conf, int sapling_stride, double &error, double &ratio)
  {
    // update the level configuration using the updated conf
    // level interp direction, interp method, error bound, etc
    quant_inds.resize(0); // reset the quant_inds
    quant_inds.reserve(num_elements);


      for (uint level = start_level;
         level >=end_level && level <= interpolation_level; level--) {
      current_level = level;
      quantizer.set_eb(level_abs_ebs[level - 1]);
      size_t stride = 1U << (level - 1);
      {
        auto inter_block_range =
            std::make_shared<SZ::multi_dimensional_range<T, N>>(
                working_data, std::begin(global_dimensions),
                std::end(global_dimensions), blocksize * stride, 0);

        auto inter_begin = inter_block_range->begin();
        auto inter_end = inter_block_range->end();

        // get the block size and sample the block at each level 

        // block grind calculation

        for (auto block = inter_begin; block != inter_end; ++block) {
          auto end_idx = block.get_global_index();
          for (int i = 0; i < N; i++) {
            end_idx[i] += blocksize * stride;
            if (end_idx[i] > global_dimensions[i] - 1) {
              end_idx[i] = global_dimensions[i] - 1;
            }
          }
          error += block_interpolation(
              working_data, block.get_global_index(), end_idx, PB_predict_overwrite,
              interpolators[interpolator_id], direction_sequence_id, stride);
        }
      }

      // get the size that is pushed into the quant_inds 
      size_t original_size = quant_inds.size(); 
      size_t compressed_size =0; 
      std::vector<uchar> buffer; 
      buffer.resize(original_size * sizeof(int));
      uchar *buffer_pos = buffer.data();
      encoder.preprocess_encode(quant_inds, 0);
      encoder.save(buffer_pos);
      encoder.encode(quant_inds, buffer_pos);
      encoder.postprocess_encode();

      lossless.compress(buffer.data(), buffer.size(), compressed_size);
      // calculate the error and ratio
      ratio = compressed_size*1.0 / (original_size * sizeof(T));

      
    }







    // try compress the given levels

    // calculate the ratio and error(sum of error) 


  }

  uchar *compress(
      const Config &conf, T *data, size_t &compressed_size,
      bool tuning = false)
  {

    quant_pred_start_level = conf.quantization_prediction_start_level;

    eb_factors[interpolator_id] = conf.eb_factor;

    init();
    // original_variance = SZ::data_variance(data, num_elements);
    // std::cout << "original variance = " << original_variance << std::endl;
    // For data range check.
    auto orig_min_max = std::minmax_element(data, data + num_elements);
    original_min = *orig_min_max.first;
    original_max = *orig_min_max.second;
    std::cout << "original max " << original_max << std::endl;
    std::cout << "original min " << original_min << std::endl;

    original_range = original_max - original_min;

    Timer timer;

    quant_inds.reserve(num_elements);
    size_t interp_compressed_size = 0;

    double eb_input = quantizer.get_eb();
    double eb_final;  // eb for the highest level
    // double eb_reduction_factor;
    if (interpolator_id == 0) {
      eb_final = eb_input /
                 pow(linear_interp_eb_factor, (interpolation_level - 1) * N);
      eb_reduction_factor = pow(linear_interp_eb_factor, N);
    }
    else {
      eb_final = eb_input /
                 pow(cubic_interp_eb_factor, (interpolation_level - 1) * N);
      eb_reduction_factor = pow(cubic_interp_eb_factor, N);
    }

    current_level = interpolation_level;
    quantize(0, *data, 0);

    // Timer timer;
    timer.start();

    for (int i = 0; i < level_abs_ebs.size(); i++) {
      level_abs_ebs[i] *= quantizer.get_eb();
    }

    for (uint level = interpolation_level;
         level > 0 && level <= interpolation_level; level--) {
      current_level = level;
      quantizer.set_eb(level_abs_ebs[level - 1]);
      size_t stride = 1U << (level - 1);
      {
        auto inter_block_range =
            std::make_shared<SZ::multi_dimensional_range<T, N>>(
                data, std::begin(global_dimensions),
                std::end(global_dimensions), blocksize * stride, 0);

        auto inter_begin = inter_block_range->begin();
        auto inter_end = inter_block_range->end();

        // get the block size and sample the block at each level 

        for (auto block = inter_begin; block != inter_end; ++block) {
          auto end_idx = block.get_global_index();
          for (int i = 0; i < N; i++) {
            end_idx[i] += blocksize * stride;
            if (end_idx[i] > global_dimensions[i] - 1) {
              end_idx[i] = global_dimensions[i] - 1;
            }
          }
          block_interpolation(
              data, block.get_global_index(), end_idx, PB_predict_overwrite,
              interpolators[interpolator_id], direction_sequence_id, stride);
        }
      }
    }

    assert(quant_inds.size() <= num_elements);

    encoder.preprocess_encode(quant_inds, 0);
    size_t bufferSize = 1.5 * (quantizer.size_est() + encoder.size_est() +
                               sizeof(T) * quant_inds.size());

    // TODO: change to smart pointer here
    uchar *buffer = new uchar[bufferSize];
    uchar *buffer_pos = buffer;

    std::cout << "bufferSize = " << bufferSize << std::endl;

    write(global_dimensions.data(), N, buffer_pos);
    write(blocksize, buffer_pos);
    write(interpolator_id, buffer_pos);
    write(direction_sequence_id, buffer_pos);

    // add auxilliary array
    write(detection_block_size, buffer_pos);
    // num_detection_block = significant_block_id.size();
    // if(num_detection_block == 0) num_detection_block =
    // flushed_block_id.size();

    // add additional variable
    write(original_max, buffer_pos);
    write(original_min, buffer_pos);

    write(quant_pred_on, buffer_pos);
    write(quant_pred_start_level, buffer_pos);
    write(post_process_on, buffer_pos);

    quantizer.save(buffer_pos);
    quantizer.postcompress_data();

    timer.start();
    encoder.save(buffer_pos);
    encoder.encode(quant_inds, buffer_pos);
    encoder.postprocess_encode();
    //            timer.stop("Coding");
    assert(buffer_pos - buffer < bufferSize);

    // writefile("compressed.dat",data, num_elements);
    timer.start();
    uchar *lossless_data =
        lossless.compress(buffer, buffer_pos - buffer, compressed_size);

    lossless.postcompress_data(buffer);

    compressed_size += interp_compressed_size;
    return lossless_data;
  }

 private:
  enum PredictorBehavior { PB_predict_overwrite, PB_predict, PB_recover };

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

    // post process utils
    // aux_quant_inds.resize(num_elements,0);
    // aux_quant_inds_ptr =
    // std::make_shared<std::vector<int>>(aux_quant_inds.begin(),
    // aux_quant_inds.end());

    aux_quant_inds_ptr = std::make_shared<std::vector<int>>();
    aux_quant_inds_ptr->resize(num_elements, 0);

    level_abs_ebs.resize(interpolation_level);
    level_interp_ids.resize(interpolation_level);
    std::fill(
        level_interp_ids.begin(), level_interp_ids.end(), interpolator_id);
    level_abs_ebs[0] = 1.0;
    for (int i = 1; i < interpolation_level; i++) {
      level_abs_ebs[i] =
          level_abs_ebs[i - 1] / eb_factors[level_interp_ids[i - 1]];
    }
  }

  inline void range_check(T &d)
  {
    if (d > original_max) { d = original_max; }
    if (d < original_min) { d = original_min; }
  }

  inline int backward_compensate_pred(
      size_t idx, size_t offset1, size_t offset2, T &pred,
      const double compensation)
  {
    return 0;
  }

  inline void quantize(size_t idx, T &d, T pred)
  {
    quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
  }

  inline void recover(size_t idx, T &d, T pred)
  {
    d = quantizer.recover(pred, quant_inds[quant_index++]);
  };

  //   double block_interpolation_1d(
  //       T *data, size_t begin, size_t end, size_t stride,
  //       const std::string &interp_func, const PredictorBehavior pb,
  //       bool quant_pred = false, size_t offset1 = 0, size_t offset2 = 0,
  //       bool is_left_boundary = true,
  //       bool use_begin_cross= false,
  //       bool use_end_cubic = false)
  //   {
  //     quant_pred = (quant_pred_on == true) && (is_left_boundary == false) &&
  //     (current_level <= quant_pred_start_level); size_t n = (end - begin) /
  //     stride + 1; if (n <= 1) { return 0; } double predict_error = 0;

  //     int quant_compensation = 0;

  //     size_t stride3x = 3 * stride;
  //     size_t stride5x = 5 * stride;

  //     if (interp_func == "linear" || n < 5) {
  //       if (pb == PB_predict_overwrite) {
  //         for (size_t i = 1; i + 1 < n; i += 2) {
  //           T *d = data + begin + i * stride;
  //           T d_copy = *d;
  //           T pred = interp_linear(*(d - stride), *(d + stride));

  //           if (quant_pred == true) {
  //             double tol = quantizer.get_eb() / linear_interp_eb_factor;
  //             double compensation =
  //                 region_error_control_eb_compensation * quantizer.get_eb();
  //             quant_compensation = backward_compensate_pred(
  //                 d - data, offset1, offset2, pred, compensation);
  //           }
  //           quantize(d - data, *d, pred);
  //           if(quant_inds.back()==0)
  //           {
  //             (*aux_quant_inds_ptr)[d - data] = 0;
  //           }
  //           else {
  //             (*aux_quant_inds_ptr)[d - data] = quant_inds.back() +
  //             quant_compensation;
  //           }

  //           #ifdef SZ_ANALYSIS
  //           my_compensation_label[d - data] = quant_compensation;
  //           #endif

  //         }
  //         if (n % 2 == 0) {
  //           T *d = data + begin + (n - 1) * stride;
  //           T pred = interp_linear(*(d - stride), *(d - stride));
  //           #ifdef SZ_ANALYSIS
  //           my_pred[d - data] = pred;
  //           #endif
  //           if ( quant_pred == true) {
  //             double tol = quantizer.get_eb() / linear_interp_eb_factor;
  //             double compensation =
  //                 region_error_control_eb_compensation * quantizer.get_eb();
  //             quant_compensation = backward_compensate_pred(
  //                 d - data, offset1, offset2, pred, compensation);
  //             // backward_compensate(d - stride-data, d - stride-data, pred,
  //             // error_recorder.data(), compensation, tol);
  //           }
  //           quantize(d - data, *d, pred);
  //           if(quant_inds.back()==0)
  //           {
  //             (*aux_quant_inds_ptr)[d - data] = 0;
  //           }
  //           else {
  //             (*aux_quant_inds_ptr)[d - data] = quant_inds.back() +
  //             quant_compensation;
  //           }
  //         }
  //       }
  //       else {
  //         for (size_t i = 1; i + 1 < n; i += 2) {
  //           T *d = data + begin + i * stride;
  //           T pred = interp_linear(*(d - stride), *(d + stride));
  //           // pred prediction
  //           if(quant_pred == true)
  //           {
  //             double compensation = region_error_control_eb_compensation *
  //             quantizer.get_eb(); quant_compensation =
  //             backward_compensate_pred(d - data, offset1, offset2, pred,
  //             compensation);
  //           }

  //           recover(d - data, *d, pred);

  //           if(quant_inds[quant_index-1]==0)
  //           {
  //             (*aux_quant_inds_ptr)[d - data] = 0;
  //           }
  //           else {
  //             (*aux_quant_inds_ptr)[d - data] = quant_inds[quant_index-1] +
  //             quant_compensation;
  //           }

  //         }
  //         if (n % 2 == 0) {
  //           T *d = data + begin + (n - 1) * stride;
  //           T pred = *(d - stride);
  //           if(quant_pred == true)
  //           {
  //             double tol = quantizer.get_eb() / linear_interp_eb_factor;
  //             double compensation = region_error_control_eb_compensation *
  //             quantizer.get_eb(); quant_compensation =
  //             backward_compensate_pred(d - data, offset1, offset2, pred,
  //             compensation);
  //           }
  //           recover(d - data, *d, pred);
  //           if(quant_inds[quant_index-1]==0)
  //           {
  //             (*aux_quant_inds_ptr)[d - data] = 0;
  //           }
  //           else {
  //             (*aux_quant_inds_ptr)[d - data] = quant_inds[quant_index-1] +
  //             quant_compensation;
  //           }
  //         }
  //       }
  //     }
  //     else {
  //       if (pb == PB_predict_overwrite) {
  //         double dafault_eb = quantizer.get_eb();
  //         T *d;
  //         size_t i;
  //         quant_compensation =0;
  //         d = data + begin + stride;
  //         T pred;
  //         if (use_begin_cross ==false) {pred = interp_quad_1(*(d - stride),
  //         *(d + stride), *(d + stride3x));} else {pred = interp_cubic(*(d -
  //         stride3x), *(d - stride), *(d + stride), *(d + stride3x),
  //         use_natural_cubic);}

  //         if(quant_pred == true)
  //         {
  //           double compensation = region_error_control_eb_compensation *
  //           quantizer.get_eb(); quant_compensation =
  //           backward_compensate_pred(d - data, offset1, offset2, pred,
  //           compensation);
  //         }
  //         quantize(d - data, *d, pred);
  //         (*aux_quant_inds_ptr)[d - data] = quant_inds.back();
  //         if(quant_inds.back()==0)
  //         {
  //           (*aux_quant_inds_ptr)[d - data] = 0;
  //         }
  //         else {
  //           (*aux_quant_inds_ptr)[d - data] = quant_inds.back() +
  //           quant_compensation;
  //         }

  //         // d = data + begin + stride;
  //         for (i = 3; i + 3 < n; i += 2) {
  //           quant_compensation =0;
  //           d = data + begin + i * stride;
  //           T pred = interp_cubic(
  //               *(d - stride3x), *(d - stride), *(d + stride), *(d +
  //               stride3x),use_natural_cubic);

  //           if ( quant_pred == true ) {
  //             double compensation = region_error_control_eb_compensation *
  //             quantizer.get_eb(); quant_compensation =
  //             backward_compensate_pred(d - data, offset1, offset2, pred,
  //             compensation);
  //           }
  //           quantize(d - data, *d, pred);
  //           if(quant_inds.back()==0)
  //           {
  //             (*aux_quant_inds_ptr)[d - data] = 0;
  //           }
  //           else {
  //             (*aux_quant_inds_ptr)[d - data] = quant_inds.back() +
  //             quant_compensation;
  //           }

  //         }

  //         d = data + begin + i * stride;
  //         quant_compensation =0;

  //         if(use_end_cubic == false)
  //         { pred = interp_quad_2(*(d - stride3x), *(d - stride), *(d +
  //         stride));
  //         }
  //         else{
  //           pred = interp_cubic(*(d - stride3x), *(d - stride),
  //                               *(d + stride), *(d + stride3x),
  //                               use_natural_cubic);
  //           }

  //         if(quant_pred == true)
  //         {
  //           double compensation = region_error_control_eb_compensation *
  //           quantizer.get_eb(); quant_compensation =
  //           backward_compensate_pred(d - data, offset1, offset2, pred,
  //           compensation);
  //         }
  //         quantize(d - data, *d, pred);
  //         if(quant_inds.back()==0)
  //         {
  //           (*aux_quant_inds_ptr)[d - data] = 0;
  //         }
  //         else {
  //           (*aux_quant_inds_ptr)[d - data] = quant_inds.back() +
  //           quant_compensation;
  //         }

  //         if (n % 2 == 0) {
  //           d = data + begin + (n - 1) * stride;
  //           // quantize(d - data, *d, *(d - stride));
  //           T pred = *(d - stride);
  //           quant_compensation = 0;

  //           if ( quant_pred == true) {
  //             double tol = quantizer.get_eb() / linear_interp_eb_factor;
  //             double compensation =
  //                 region_error_control_eb_compensation * quantizer.get_eb();
  //             quant_compensation = backward_compensate_pred(
  //                 d - data, offset1, offset2, pred, compensation);
  //             // backward_compensate(d - stride-data, d - stride-data, pred,
  //             // error_recorder.data(), compensation, tol);
  //           }
  //           quantize(d - data, *d, pred);
  //           if(quant_inds.back()==0)
  //           {
  //             (*aux_quant_inds_ptr)[d - data] = 0;
  //           }
  //           else {
  //             (*aux_quant_inds_ptr)[d - data] = quant_inds.back() +
  //             quant_compensation;
  //           }
  //           #ifdef SZ_ANALYSIS
  //           my_compensation_label[d - data] = quant_compensation;
  //           #endif
  //           quantizer.set_eb(dafault_eb);
  //         }
  //       }
  //       else {
  //         T *d;
  //         size_t i;
  //         quant_compensation = 0;

  //         d = data + begin + stride;
  //         T pred;

  //         if (use_begin_cross ==false) {pred = interp_quad_1(*(d - stride),
  //         *(d + stride), *(d + stride3x));} else {pred = interp_cubic(*(d -
  //         stride3x), *(d - stride), *(d + stride), *(d +
  //         stride3x),use_natural_cubic);}

  //         if(quant_pred == true)
  //         {
  //           double compensation = region_error_control_eb_compensation *
  //           quantizer.get_eb(); quant_compensation =
  //           backward_compensate_pred(d - data, offset1, offset2, pred,
  //           compensation);
  //         }
  //         recover(d - data, *d, pred);
  //         if(quant_inds[quant_index-1]==0)
  //         {
  //           (*aux_quant_inds_ptr)[d - data] = 0;
  //         }
  //         else {
  //           (*aux_quant_inds_ptr)[d - data] = quant_inds[quant_index-1] +
  //           quant_compensation;
  //         }

  //         for (i = 3; i + 3 < n; i += 2) {
  //           d = data + begin + i * stride;
  //           quant_compensation = 0;
  //           T pred = interp_cubic(
  //               *(d - stride3x), *(d - stride), *(d + stride), *(d +
  //               stride3x),use_natural_cubic);
  //           if(quant_pred == true)
  //           {
  //             double compensation = region_error_control_eb_compensation *
  //             quantizer.get_eb(); quant_compensation =
  //             backward_compensate_pred(d - data, offset1, offset2, pred,
  //             compensation);
  //           }
  //           recover(d - data, *d, pred);
  //           if(quant_inds[quant_index-1]==0)
  //           {
  //             (*aux_quant_inds_ptr)[d - data] = 0;
  //           }
  //           else {
  //             (*aux_quant_inds_ptr)[d - data] = quant_inds[quant_index-1] +
  //             quant_compensation;
  //           }
  //         }

  //         d = data + begin + i * stride;
  //         quant_compensation = 0;

  //         if(use_end_cubic == false)
  //         {
  //           pred = interp_quad_2(*(d - stride3x), *(d - stride), *(d +
  //           stride));
  //         }
  //         else {
  //           pred = interp_cubic(*(d - stride3x), *(d - stride),
  //                                   *(d + stride), *(d +
  //                                   stride3x),use_natural_cubic);
  //         }

  //         if(quant_pred == true)
  //         {
  //           double compensation = region_error_control_eb_compensation *
  //           quantizer.get_eb(); quant_compensation =
  //           backward_compensate_pred(d - data, offset1, offset2, pred,
  //           compensation);
  //         }
  //         recover(d - data, *d, pred); // end
  //         if(quant_inds[quant_index-1]==0)
  //         {
  //           (*aux_quant_inds_ptr)[d - data] = 0;
  //         }
  //         else {
  //           (*aux_quant_inds_ptr)[d - data] = quant_inds[quant_index-1] +
  //           quant_compensation;
  //         }
  //         // recover(
  //         //     d - data, *d,
  //         //     interp_quad_2(*(d - stride3x), *(d - stride), *(d +
  //         stride))); // end

  //         if (n % 2 == 0) {
  //           d = data + begin + (n - 1) * stride;
  //           T pred = *(d - stride);
  //           quant_compensation = 0;
  //           if(quant_pred == true)
  //           {
  //             double compensation = region_error_control_eb_compensation *
  //             quantizer.get_eb(); quant_compensation =
  //             backward_compensate_pred(d - data, offset1, offset2, pred,
  //             compensation);
  //           }
  //           recover(d - data, *d, pred); // last element on the line sample
  //           if(quant_inds[quant_index-1]==0)
  //           {
  //             (*aux_quant_inds_ptr)[d - data] = 0;
  //           }
  //           else {
  //             (*aux_quant_inds_ptr)[d - data] = quant_inds[quant_index-1] +
  //             quant_compensation;
  //           }
  //         }
  //       }
  //     }

  //     return predict_error;
  //   }

  double block_interpolation_1d(
      T *data, size_t begin, size_t end, size_t stride,
      const std::string &interp_func, const PredictorBehavior pb,
      bool quant_pred = false, size_t offset1 = 0, size_t offset2 = 0,
      bool is_left_boundary = true, bool use_fft_interp = false)
  {
    size_t n = (end - begin) / stride + 1;
    if (n <= 1) { return 0; }
    double predict_error = 0;

    size_t stride3x = 3 * stride;
    size_t stride5x = 5 * stride;
    if (interp_func == "linear" || n < 5) {
      if (pb == PB_predict_overwrite) {
        for (size_t i = 1; i + 1 < n; i += 2) {
          T *d = data + begin + i * stride;
          quantize(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
        }
        if (n % 2 == 0) {
          T *d = data + begin + (n - 1) * stride;
          if (n < 4) { quantize(d - data, *d, *(d - stride)); }
          else {
            quantize(
                d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
          }
        }
      }
      else {
        for (size_t i = 1; i + 1 < n; i += 2) {
          T *d = data + begin + i * stride;
          recover(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
        }
        if (n % 2 == 0) {
          T *d = data + begin + (n - 1) * stride;
          if (n < 4) { recover(d - data, *d, *(d - stride)); }
          else {
            recover(
                d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
          }
        }
      }
    }
    else {
      if (pb == PB_predict_overwrite) {
        T *d;
        size_t i;
        for (i = 3; i + 3 < n; i += 2) {
          d = data + begin + i * stride;

          quantize(
              d - data, *d,
              interp_cubic(
                  *(d - stride3x), *(d - stride), *(d + stride),
                  *(d + stride3x), use_natural_cubic));
        }
        d = data + begin + stride;

        quantize(
            d - data, *d,
            interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

        d = data + begin + i * stride;

        quantize(
            d - data, *d,
            interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
        if (n % 2 == 0) {
          d = data + begin + (n - 1) * stride;

          quantize(
              d - data, *d,
              interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
        }
      }
      else {
        T *d;

        size_t i;
        for (i = 3; i + 3 < n; i += 2) {
          d = data + begin + i * stride;
          recover(
              d - data, *d,
              interp_cubic(
                  *(d - stride3x), *(d - stride), *(d + stride),
                  *(d + stride3x), use_natural_cubic));
        }
        d = data + begin + stride;

        recover(
            d - data, *d,
            interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

        d = data + begin + i * stride;
        recover(
            d - data, *d,
            interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));

        if (n % 2 == 0) {
          d = data + begin + (n - 1) * stride;
          recover(
              d - data, *d,
              interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
        }
      }
    }

    return predict_error;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 1, double>::type block_interpolation(
      T *data, std::array<size_t, N> begin, std::array<size_t, N> end,
      const PredictorBehavior pb, const std::string &interp_func,
      const int direction, size_t stride = 1)
  {
    return block_interpolation_1d(
        data, begin[0], end[0], stride, interp_func, pb);
  }

  template <uint NN = N>
  typename std::enable_if<NN == 2, double>::type block_interpolation(
      T *data, std::array<size_t, N> begin, std::array<size_t, N> end,
      const PredictorBehavior pb, const std::string &interp_func,
      const int direction, size_t stride = 1)
  {
    double predict_error = 0;
    size_t stride2x = stride * 2;
    // auto default_eb = quantizer.get_eb();

    // quantizer.set_eb(default_eb * c1);
    const std::array<int, N> dims = dimension_sequences[direction];
    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
         j <= end[dims[1]]; j += stride2x) {
      size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] +
                            j * dimension_offsets[dims[1]];
      predict_error += block_interpolation_1d(
          data, begin_offset,
          begin_offset +
              (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
          stride * dimension_offsets[dims[0]], interp_func, pb);
    }
    // quantizer.set_eb(default_eb);
    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
         i <= end[dims[0]]; i += stride) {
      size_t begin_offset = i * dimension_offsets[dims[0]] +
                            begin[dims[1]] * dimension_offsets[dims[1]];
      predict_error += block_interpolation_1d(
          data, begin_offset,
          begin_offset +
              (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
          stride * dimension_offsets[dims[1]], interp_func, pb);
    }
    return predict_error;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 3, double>::type block_interpolation(
      T *data, std::array<size_t, N> begin, std::array<size_t, N> end,
      const PredictorBehavior pb, const std::string &interp_func,
      const int direction, size_t stride = 1)
  {
    double predict_error = 0;
    size_t stride2x = stride * 2;

    auto default_eb = quantizer.get_eb();

    const std::array<int, N> dims = dimension_sequences[direction];

    bool use_begin_cross = (begin[dims[0]] != 0) && (use_cross_block_cubic);
    bool is_cubic_end_boundary =
        (end[dims[0]] + 1 == global_dimensions[dims[0]]);
    bool is_next_to_bound =
        (end[dims[0]] + stride2x + 1 > global_dimensions[dims[0]]);
    bool is_cubic_end_boundary_ = is_cubic_end_boundary || is_next_to_bound;
    bool use_end_cross = (!is_cubic_end_boundary_) && use_cross_block_cubic;

    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
         j <= end[dims[1]]; j += stride2x) {
      for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
           k <= end[dims[2]]; k += stride2x) {
        size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] +
                              j * dimension_offsets[dims[1]] +
                              k * dimension_offsets[dims[2]];
        // bool high_level_predict = (j%4==0) && (k%4==0);
        // high_level_predict = 0;
        bool is_left_boundary = (j == 0) || (k == 0);
        use_end_cross =
            (j % stride2x == 0) && (k % stride2x == 0) && use_end_cross;
        bool is_cubic_end_boundary =
            (end[dims[0]] + 1 == global_dimensions[dims[0]]) ||
            (end[dims[0]] + 1 != global_dimensions[dims[0]] &&
             end[dims[0]] + stride > global_dimensions[dims[0]] - 1) ||
            !(use_cross_block_cubic);
        predict_error += block_interpolation_1d(
            data, begin_offset,
            begin_offset +
                (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
            stride * dimension_offsets[dims[0]], interp_func, pb,
            quant_pred_on, dimension_offsets[dims[1]] * (stride2x),
            dimension_offsets[dims[2]] * (stride2x), is_left_boundary,
            use_begin_cross, use_end_cross);
      }
    }

    use_begin_cross = (begin[dims[1]] != 0) && (use_cross_block_cubic);
    is_cubic_end_boundary = (end[dims[1]] + 1 == global_dimensions[dims[1]]);
    // ||(!use_cross_block_cubic);
    is_next_to_bound =
        (end[dims[1]] + stride2x + 1 > global_dimensions[dims[1]]);
    is_cubic_end_boundary_ = is_cubic_end_boundary || is_next_to_bound;
    use_end_cross = (!is_cubic_end_boundary_) && use_cross_block_cubic;

    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
         i <= end[dims[0]]; i += stride) {
      for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
           k <= end[dims[2]]; k += stride2x) {
        size_t begin_offset = i * dimension_offsets[dims[0]] +
                              begin[dims[1]] * dimension_offsets[dims[1]] +
                              k * dimension_offsets[dims[2]];
        // if( i%stride2x ==0 || k%stride2x ==0)
        // {
        //   quantizer.set_eb(default_eb*0.8);
        // }
        bool is_left_boundary = (i == 0) || (k == 0);
        use_end_cross =
            (i % stride2x == 0) && (k % stride2x == 0) && use_end_cross;
        predict_error += block_interpolation_1d(
            data, begin_offset,
            begin_offset +
                (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
            stride * dimension_offsets[dims[1]], interp_func, pb,
            quant_pred_on, dimension_offsets[dims[0]] * (stride),
            dimension_offsets[dims[2]] * (stride2x), is_left_boundary,
            use_begin_cross, use_end_cross);
        // quantizer.set_eb(default_eb);
      }
    }

    use_begin_cross = (begin[dims[2]] != 0) && (use_cross_block_cubic);
    is_cubic_end_boundary = (end[dims[2]] + 1 == global_dimensions[dims[2]]);
    is_next_to_bound =
        (end[dims[2]] + stride2x + 1 > global_dimensions[dims[2]]);
    is_cubic_end_boundary_ = is_cubic_end_boundary || is_next_to_bound;
    use_end_cross = (!is_cubic_end_boundary_) && use_cross_block_cubic;

    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
         i <= end[dims[0]]; i += stride) {
      for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0);
           j <= end[dims[1]]; j += stride) {
        // if (i==0 && current_level==1) std::cout <<"j = " << j << std::endl;
        size_t begin_offset = i * dimension_offsets[dims[0]] +
                              j * dimension_offsets[dims[1]] +
                              begin[dims[2]] * dimension_offsets[dims[2]];
        bool is_left_boundary =
            (i == 0) || (j == 0);  // this is for quant prediction.
        use_end_cross =
            (i % stride2x == 0) && (j % stride2x == 0) && use_end_cross;
        predict_error += block_interpolation_1d(
            data, begin_offset,
            begin_offset +
                (end[dims[2]] - begin[dims[2]]) * dimension_offsets[dims[2]],
            stride * dimension_offsets[dims[2]], interp_func, pb,
            quant_pred_on, dimension_offsets[dims[0]] * stride,
            dimension_offsets[dims[1]] * stride, is_left_boundary,
            use_begin_cross, use_end_cross);
      }
    }
    quantizer.set_eb(default_eb);
    return predict_error;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 4, double>::type block_interpolation(
      T *data, std::array<size_t, N> begin, std::array<size_t, N> end,
      const PredictorBehavior pb, const std::string &interp_func,
      const int direction, size_t stride = 1)
  {
    double predict_error = 0;
    size_t stride2x = stride * 2;
    max_error = 0;
    const std::array<int, N> dims = dimension_sequences[direction];
    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
         j <= end[dims[1]]; j += stride2x) {
      for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
           k <= end[dims[2]]; k += stride2x) {
        for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
             t <= end[dims[3]]; t += stride2x) {
          size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] +
                                j * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
          predict_error += block_interpolation_1d(
              data, begin_offset,
              begin_offset +
                  (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
              stride * dimension_offsets[dims[0]], interp_func, pb);
        }
      }
    }
    max_error = 0;
    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
         i <= end[dims[0]]; i += stride) {
      for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
           k <= end[dims[2]]; k += stride2x) {
        for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
             t <= end[dims[3]]; t += stride2x) {
          size_t begin_offset = i * dimension_offsets[dims[0]] +
                                begin[dims[1]] * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
          predict_error += block_interpolation_1d(
              data, begin_offset,
              begin_offset +
                  (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
              stride * dimension_offsets[dims[1]], interp_func, pb);
        }
      }
    }
    max_error = 0;
    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
         i <= end[dims[0]]; i += stride) {
      for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0);
           j <= end[dims[1]]; j += stride) {
        for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
             t <= end[dims[3]]; t += stride2x) {
          size_t begin_offset = i * dimension_offsets[dims[0]] +
                                j * dimension_offsets[dims[1]] +
                                begin[dims[2]] * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
          predict_error += block_interpolation_1d(
              data, begin_offset,
              begin_offset +
                  (end[dims[2]] - begin[dims[2]]) * dimension_offsets[dims[2]],
              stride * dimension_offsets[dims[2]], interp_func, pb);
        }
      }
    }

    max_error = 0;
    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
         i <= end[dims[0]]; i += stride) {
      for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0);
           j <= end[dims[1]]; j += stride) {
        for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride : 0);
             k <= end[dims[2]]; k += stride) {
          size_t begin_offset = i * dimension_offsets[dims[0]] +
                                j * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                begin[dims[3]] * dimension_offsets[dims[3]];
          predict_error += block_interpolation_1d(
              data, begin_offset,
              begin_offset +
                  (end[dims[3]] - begin[dims[3]]) * dimension_offsets[dims[3]],
              stride * dimension_offsets[dims[3]], interp_func, pb);
        }
      }
    }
    return predict_error;
  }

  int interpolation_level = -1;
  uint blocksize;
  int interpolator_id;
  double eb_ratio = 0.5;
  std::vector<std::string> interpolators = {"linear", "cubic"};
  double linear_interp_eb_factor = sqrt(1.5);
  double cubic_interp_eb_factor = 1.2808688457449497;

  // double cubic_interp_eb_factor = std::cbrt(2);
  double cunic_interp_eb_factor_natural = 1.2932517156377563;
  std::array<double, 2> eb_factors = {
      pow(linear_interp_eb_factor, N), pow(cunic_interp_eb_factor_natural, N)};
  std::vector<int> quant_inds;
  size_t quant_index = 0;  // for decompress
  double max_error;
  Quantizer quantizer;
  Encoder encoder;
  Lossless lossless;
  size_t num_elements;
  std::array<size_t, N> global_dimensions;
  std::array<size_t, N> dimension_offsets;
  std::vector<std::array<int, N>> dimension_sequences;
  int direction_sequence_id;

  int current_level = 0;

  std::vector<double> level_abs_ebs;
  std::vector<int> level_interp_ids;

  int detection_block_size = 4;
  double detection_threshold = 0.9;
  double detection_eb_rate;
  double noise_rate = 0;
  double eb_reduction_factor = 1.0;

  T original_max;    // for data range check
  T original_min;    // for data range check
  T original_range;  // for data range check

  // original data copy
  const T *orig_data_ptr;
  double original_variance;
  T *working_data;

  // std::vector<int> aux_quant_inds;
  std::shared_ptr<std::vector<int>> aux_quant_inds_ptr;
  double region_error_control_eb_compensation = 2.0;  // for quant prediction
  bool quant_pred_on = 0;
  int quant_pred_start_level = 3;
  bool post_process_on = true;

  bool use_natural_cubic = false;
  bool use_cross_block_cubic = false;

  // Analysis utils;
  // This is for conditional compilation not comments
  // #ifndef SZ_ANALYSIS
  // #define SZ_ANALYSIS
};
};  // namespace SZ

#endif
