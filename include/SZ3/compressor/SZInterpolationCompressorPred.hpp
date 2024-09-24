#ifndef _SZ_INTERPOLATION_COMPRESSSORPRED_HPP
#define _SZ_INTERPOLATION_COMPRESSSORPRED_HPP

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <future>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <vector>
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"


namespace SZ3 {
template <class T, uint N, class Quantizer, class Encoder, class Lossless>
class SZInterpolationCompressorPred{
 public:
  SZInterpolationCompressorPred(
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
    size_t remaining_length = cmpSize;
    uchar *buffer = lossless.decompress(cmpData, remaining_length);
    uchar const *buffer_pos = buffer;

    read(global_dimensions.data(), N, buffer_pos, remaining_length);
    read(blocksize, buffer_pos, remaining_length);
    read(interpolator_id, buffer_pos, remaining_length);
    read(direction_sequence_id, buffer_pos, remaining_length);

    Timer timer;
    timer.start();

    // read smoothing and pred quant
    read(quant_pred_on, buffer_pos);
    read(quant_pred_start_level, buffer_pos);

    std::cout << "quant_pred_on = " << quant_pred_on << std::endl;
    std::cout << "quant_pred_start_level = " << quant_pred_start_level << std::endl;

    init();


    quantizer.load(buffer_pos, remaining_length);
    encoder.load(buffer_pos, remaining_length);
    quant_inds =
        encoder.decode(buffer_pos, num_elements);

    encoder.postprocess_decode();

    lossless.postdecompress_data(buffer);
    
    double eb_input = quantizer.get_eb();
    double eb_final;

    std::cout << "start decompression\n";
    current_level = interpolation_level;
    // *decData = quantizer.recover(0, quant_inds[quant_index++]);
    recover(0, *decData, 0);


    for (uint level = interpolation_level;
         level > 0 && level <= interpolation_level; level--) {
      if (level >= 3) {
          quantizer.set_eb(eb_input * eb_ratio);
      } else {
          quantizer.set_eb(eb_input);
      }
      current_level = level;
    
    { 
      size_t stride = 1U << (level - 1);
      auto inter_block_range =
          std::make_shared<multi_dimensional_range<T, N>>(
              decData, std::begin(global_dimensions),
              std::end(global_dimensions), stride * blocksize, 0);
      auto inter_begin = inter_block_range->begin();
      auto inter_end = inter_block_range->end();
      for (auto block = inter_begin; block != inter_end; ++block) {
        auto end_idx = block.get_global_index();
        for (int i = 0; i < N; i++) {
          end_idx[i] += stride * blocksize;
          if (end_idx[i] > global_dimensions[i] - 1) {
            end_idx[i] = global_dimensions[i] - 1;
          }
        }

        block_interpolation(
            decData, block.get_global_index(), end_idx, PB_recover,
            interpolators[interpolator_id], direction_sequence_id, stride);
      }
    }

    }
    quantizer.postdecompress_data();
    // writefile("quant_de.dat", aux_quant_inds_ptr->data(), num_elements);


    return decData;
  }


  uchar *compress(const Config &conf, T *data, size_t &compressed_size, bool tuning = false)
  {
    std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
    blocksize = conf.interpBlockSize;
    interpolator_id = conf.interpAlgo;
    direction_sequence_id = conf.interpDirection;
    quant_pred_on = conf.quantization_prediction_on;
    quant_pred_start_level = conf.quantization_prediction_start_level;

    // interpolators tuning 
    use_cross_block_cubic = conf.corss_block_cubic;
    use_natural_cubic = conf.use_natural_cubic;

    init();

    Timer timer;

    quant_inds.reserve(num_elements);
    size_t interp_compressed_size = 0;


    double eb_input = quantizer.get_eb();
    double eb_final;  // eb for the highest level


    // quant_inds.push_back(quantizer.quantize_and_over*data, 0));
    current_level = interpolation_level;
    quantize(0, *data, 0);

    // Timer timer;
    timer.start();


    for (uint level = interpolation_level;
         level > 0 && level <= interpolation_level; level--) {
      if (level >= 3) {
          quantizer.set_eb(eb_input * eb_ratio);
      } else {
          quantizer.set_eb(eb_input);
      }
      current_level = level;
      size_t stride = 1U << (level - 1);

      {
      auto inter_block_range =
          std::make_shared<multi_dimensional_range<T, N>>(
              data, std::begin(global_dimensions), std::end(global_dimensions),
              blocksize * stride, 0);

      auto inter_begin = inter_block_range->begin();
      auto inter_end = inter_block_range->end();

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

    std::cout << "compression loop = " << timer.stop() << std::endl;

    assert(quant_inds.size() <= num_elements);

    encoder.preprocess_encode(quant_inds, 0);
    size_t bufferSize = 1.5 * (quantizer.size_est() + encoder.size_est() +
                               sizeof(T) * quant_inds.size());

    // TODO: change to smart pointer here
    uchar *buffer = new uchar[bufferSize];
    uchar *buffer_pos = buffer;

    // std::cout << "bufferSize = " << bufferSize << std::endl; 

    write(global_dimensions.data(), N, buffer_pos);
    write(blocksize, buffer_pos);
    write(interpolator_id, buffer_pos);
    write(direction_sequence_id, buffer_pos);


    // write smoothing and pred quant
    write(quant_pred_on, buffer_pos);
    write(quant_pred_start_level, buffer_pos);

    // write interps 
    // write(use_cross_block_cubic, buffer_pos);
    // write(use_natural_cubic, buffer_pos);

    if(tuning == false) quantizer.print();
    quantizer.save(buffer_pos);
    quantizer.postcompress_data();
    timer.start();
    encoder.save(buffer_pos);
    encoder.encode(quant_inds, buffer_pos);
    encoder.postprocess_encode();
    //            timer.stop("Coding");
    assert(buffer_pos - buffer < bufferSize);
    
    timer.start();
    uchar *lossless_data =
        lossless.compress(buffer, buffer_pos - buffer, compressed_size);
    lossless.postcompress_data(buffer);


// writefile("decompressed.dat", data, num_elements);
#ifdef SZ_ANALYSIS
    // writefile("pred.dat", my_pred.data(), num_elements);
    writefile("quant.dat", aux_quant_inds_ptr->data(), num_elements);
    writefile("quant_processed.dat", my_quant_inds_copy.data(), num_elements);
    writefile("decompressed.dat", data, num_elements);
    // writefile("level.dat", my_level.data(), num_elements);
    // writefile(
    //     "interp_direction.dat", my_interp_direction.data(), num_elements);
    // writefile(
    //     "compensation_label.int32", my_compensation_label.data(),
    //     my_compensation_label.size());
    // std::cout << "[ANALYSIS COMPILATION MODE]" << std::endl;

    // Try huffman encoding on inorder quantization and level-wise quantization
    // integers
    // SZ::HuffmanEncoder<int> huff_coding = SZ::HuffmanEncoder<int>();

    // uchar *test_buffer = new uchar[bufferSize];
    // uchar *test_buffer_pos = test_buffer;
    // huff_coding.preprocess_encode(aux_quant_inds, 0);
    // huff_coding.save(test_buffer_pos);
    // huff_coding.encode(aux_quant_inds, test_buffer_pos);
    // huff_coding.postprocess_encode();
    // size_t comsize = 0;
    // uchar *lossless_data2 =
    //     lossless.compress(test_buffer, test_buffer_pos - test_buffer,
    //     comsize);
    // // lossless.postcompress_data(test_buffer);

    // std::cout << "[inorder]comsize = " << comsize << std::endl;
    // free(test_buffer);

    // // try it on level-wise quantization integers
    // SZ::HuffmanEncoder<int> huff_coding2 = SZ::HuffmanEncoder<int>();

    // uchar *test_buffer2 = new uchar[bufferSize];
    // uchar *test_buffer_pos2 = test_buffer2;
    // huff_coding2.preprocess_encode(quant_inds, 0);
    // huff_coding2.save(test_buffer_pos2);
    // huff_coding2.encode(quant_inds, test_buffer_pos2);
    // huff_coding2.postprocess_encode();

    // size_t comsize2 = 0;
    // uchar *lossless_data3 = lossless.compress(
    //     test_buffer2, test_buffer_pos2 - test_buffer2, comsize2);
    // // lossless.postcompress_data(test_buffer2);

    // std::cout << "[level-wise]comsize = " << comsize2 << std::endl;

    // std::cout << "[level-wise] huffman encoding size = " << test_buffer_pos2
    // - test_buffer2 << std::endl;
    // free(test_buffer2);

    // writefile("quant_level.dat", quant_inds.data(), quant_inds.size());

#endif

    compressed_size += interp_compressed_size;
    return lossless_data;
  }

  std::shared_ptr<std::vector<int>> get_aux_quant_inds_ptr()
  {
    // used to pass the auxilliary quantization indices to the post process
    // or other modules
    return aux_quant_inds_ptr;
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
    }

    dimension_sequences = std::vector<std::array<int, N>>();
    auto sequence = std::array<int, N>();
    for (int i = 0; i < N; i++) { sequence[i] = i; }
    do {
      dimension_sequences.push_back(sequence);
    } while (std::next_permutation(sequence.begin(), sequence.end()));

    // post process utils 
    // aux_quant_inds.resize(num_elements,0);
    // aux_quant_inds_ptr = std::make_shared<std::vector<int>>(aux_quant_inds.begin(), aux_quant_inds.end());

    aux_quant_inds_ptr = std::make_shared<std::vector<int>>();
    if(quant_pred_on == true) aux_quant_inds_ptr->resize(num_elements, 0);

    

#ifdef SZ_ANALYSIS
    my_level.resize(num_elements);
    my_quant_inds_copy.resize(num_elements);
    my_pred.resize(num_elements);
    my_pred[0] = 0;
    my_level[0] = interpolation_level;
    my_pred_noise.resize(num_elements, 0);
    my_interp_direction.resize(num_elements, 0);
    my_compensation_label.resize(num_elements, 0);
#endif

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
    // return 0;
    // compensation = ? * eb; 0.5, 1.1, 1.5, 2
    //    0   1  2
    // 0 *A  *B *C
    // 1 *D  *E *F
    // 2 *G  *H  X
    // one layer is enough
    int A = idx - offset1 * 2 - offset2 * 2;
    int B = idx - offset1 * 2 - offset2;
    int C = idx - offset1 * 2;
    int D = idx - offset1 - offset2 * 2;
    int E = idx - offset1 - offset2;
    int F = idx - offset1;
    int G = idx - offset2 * 2;
    int H = idx - offset2;
    if ((*aux_quant_inds_ptr)[E] == INT_MIN || (*aux_quant_inds_ptr)[F] == INT_MIN ||
        (*aux_quant_inds_ptr)[H] == INT_MIN) {
        std::cout << "idx = " << idx << std::endl;
        std::cout << "current_level = " << current_level << std::endl;
        std::cout << "offset1 = " << offset1 << std::endl;
        std::cout << "offset2 = " << offset2 << std::endl;
        std::cout << "E = " << E << std::endl;
        exit(0);
      return 0;
    }

    if ((*aux_quant_inds_ptr)[E] == 0 || (*aux_quant_inds_ptr)[F] == 0 ||
        (*aux_quant_inds_ptr)[H] == 0) {
      return 0;
    }
    // int quant_A = (*aux_quant_inds_ptr)[A] - quantizer.get_radius();
    // int quant_B = (*aux_quant_inds_ptr)[B] - quantizer.get_radius();
    // int quant_C = (*aux_quant_inds_ptr)[C] - quantizer.get_radius();
    // int quant_D = (*aux_quant_inds_ptr)[D] - quantizer.get_radius(); // be cautious about out-of-bound access
    int quant_E = (*aux_quant_inds_ptr)[E] - quantizer.get_radius();
    int quant_F = (*aux_quant_inds_ptr)[F] - quantizer.get_radius();
    // int quant_G = (*aux_quant_inds_ptr)[G] - quantizer.get_radius();
    int quant_H = (*aux_quant_inds_ptr)[H] - quantizer.get_radius();

    int quant_compensate = 0;
    if (quant_H > 0 && quant_F > 0) {
      quant_compensate = (quant_H + quant_F - quant_E);
      pred += quant_compensate * compensation;
    }
    else if (quant_H < 0 && quant_F < 0) {
      quant_compensate = (quant_H + quant_F - quant_E);
      pred+= quant_compensate * compensation;

    }
    else {
      return 0;
    }
    return quant_compensate;
  }



  inline void quantize(size_t idx, T &d, T pred)
  {
    quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
    // (*aux_quant_inds_ptr)[idx] = quant_inds.back();
#ifdef SZ_ANALYSIS
    // my_level[idx] = current_level;
    my_quant_inds_copy[idx] = quant_inds.back(); // post quant index
    // my_pred[idx] = pred;
    // my_interp_direction[idx] = my_current_interp_direction;
#endif
  }


  inline void quant_pred_quantize(size_t idx, T &d, T &pred, size_t offset1, size_t offset2,  bool quant_record)
  {
    int quant_compensation = 0;
    if(quant_record == false) 
    {
    double compensation =
                region_error_control_eb_compensation * quantizer.get_eb();
    quant_compensation = backward_compensate_pred(
        idx, offset1, offset2, pred, compensation);
    }
    
    quantize(idx, d, pred);

    if(quant_inds.back()==0)
    {
      (*aux_quant_inds_ptr)[idx] = 0;
    }
    else {
      (*aux_quant_inds_ptr)[idx] = quant_inds.back() + quant_compensation;
    }
    #ifdef SZ_ANALYSIS
    // my_level[idx] = current_level;
    // my_quant_inds_copy[idx] = quant_inds.back();
    // my_pred[idx] = pred;
    // my_interp_direction[idx] = my_current_interp_direction;
    #endif
    
  }

  inline void recover(size_t idx, T &d, T pred)
  {
    d = quantizer.recover(pred, quant_inds[quant_index++]);
  };


  inline void quant_pred_recover(size_t idx, T &d, T pred,size_t offset1, size_t offset2,  bool quant_record)
  {
    int quant_compensation = 0;
    if(quant_record == false)
    {double  compensation =
                region_error_control_eb_compensation * quantizer.get_eb();
    quant_compensation = backward_compensate_pred(
        idx, offset1, offset2, pred, compensation); 
    }

    recover(idx, d, pred);

    if(quant_inds[quant_index-1]==0)
    {
      (*aux_quant_inds_ptr)[idx] = 0;
    }
    else {
      (*aux_quant_inds_ptr)[idx] = quant_inds[quant_index-1] + quant_compensation;
    }
  }

  double block_interpolation_1d(
      T *data, size_t begin, size_t end, size_t stride,
      const std::string &interp_func, const PredictorBehavior pb,
      bool quant_pred = false, size_t offset1 = 0, size_t offset2 = 0, 
      bool is_left_boundary = true,
      bool use_begin_cross= false, 
      bool use_end_cubic = false)
  {
    quant_pred = (quant_pred_on == true) && (current_level <= quant_pred_start_level);
    bool quant_record = (is_left_boundary == true) ;
    size_t n = (end - begin) / stride + 1;
    if (n <= 1) { return 0; }
    double predict_error = 0;

    size_t stride3x = 3 * stride;
    size_t stride5x = 5 * stride;

    if (interp_func == "linear" || n < 5) {
      if (pb == PB_predict_overwrite) {
        for (size_t i = 1; i + 1 < n; i += 2) {
          T *d = data + begin + i * stride;
          T pred = interp_linear(*(d - stride), *(d + stride));

          if (quant_pred == true) {
            quant_pred_quantize(d-data, *d, pred, offset1, offset2, quant_record);        
          }
          else{
            quantize(d - data, *d, pred);
          }


        }
        if (n % 2 == 0) {
          T *d = data + begin + (n - 1) * stride;
          T pred ; 
          if(n<4) {pred = *(d - stride);}
          else {pred = interp_linear1(*(d - stride3x), *(d - stride));}          
          
          if (quant_pred == true) {
            quant_pred_quantize(d-data, *d, pred, offset1, offset2, quant_record);        
          }
          else{
            quantize(d - data, *d, pred);
          }

        }
      }
      else {
        for (size_t i = 1; i + 1 < n; i += 2) {
          T *d = data + begin + i * stride;
          T pred = interp_linear(*(d - stride), *(d + stride));

          // pred prediction 
          if (quant_pred == true) {
            quant_pred_recover(d-data, *d, pred, offset1, offset2, quant_record);        
          }
          else{
            recover(d - data, *d, pred);
          }

        }
        if (n % 2 == 0) {
          T *d = data + begin + (n - 1) * stride;
          T pred ; 
          if(n<4) {pred = *(d - stride);}
          else {pred = interp_linear1(*(d - stride3x), *(d - stride));}
          if (quant_pred == true) {
            quant_pred_recover(d-data, *d, pred, offset1, offset2, quant_record);        
          }
          else{
            recover(d - data, *d, pred);
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
          T pred = interp_cubic(
              *(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x),use_natural_cubic);

          if (quant_pred == true) {
            quant_pred_quantize(d-data, *d, pred, offset1, offset2, quant_record);        
          }
          else{
            quantize(d - data, *d, pred);
          }
        }
        T pred;
        d = data + begin + stride;
        if (use_begin_cross ==false) {pred = interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x));}
        else {pred = interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x), use_natural_cubic);}

        if (quant_pred == true) {
          quant_pred_quantize(d-data, *d, pred, offset1, offset2, quant_record);        
        }
        else{
          quantize(d - data, *d, pred);
        }



        d = data + begin + i * stride;

        if(use_end_cubic == false)
        { pred = interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride));
        }
        else{ 
          use_cross_end_idx[d-data] = (begin + i * stride)/dimension_offsets[0];
          pred = interp_cubic(*(d - stride3x), *(d - stride), 
                              *(d + stride), *(d + stride3x),
                              use_natural_cubic);
        }         

        if (quant_pred == true) {
          quant_pred_quantize(d-data, *d, pred, offset1, offset2, quant_record);        
        }
        else{
          quantize(d - data, *d, pred);
        }


        if (n % 2 == 0) {
          d = data + begin + (n - 1) * stride;
          T pred = interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride));

          if (quant_pred == true) {
            quant_pred_quantize(d-data, *d, pred, offset1, offset2, quant_record);        
          }
          else{
            quantize(d - data, *d, pred);
          }

        }
      }
      else {
        T *d;
        size_t i;
        for (i = 3; i + 3 < n; i += 2) {
          d = data + begin + i * stride;
          T pred = interp_cubic(
              *(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x),use_natural_cubic);
          if (quant_pred == true) {
            quant_pred_recover(d-data, *d, pred, offset1, offset2, quant_record);        
          }
          else{
            recover(d - data, *d, pred);
          }
        }
        T pred;
        d = data + begin + stride;  
        if (use_begin_cross ==false) {pred = interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x));}
        else {pred = interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x),use_natural_cubic);}

        if (quant_pred == true) {
          quant_pred_recover(d-data, *d, pred, offset1, offset2, quant_record);        
        }
        else{
          recover(d - data, *d, pred);
        }   


        d = data + begin + i * stride;

        if(use_end_cubic == false) 
        {
          pred = interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride));
        }
        else {
          pred = interp_cubic(*(d - stride3x), *(d - stride), 
                                  *(d + stride), *(d + stride3x),use_natural_cubic);
        }

        if (quant_pred == true) {
          quant_pred_recover(d-data, *d, pred, offset1, offset2, quant_record);        
        }
        else{
          recover(d - data, *d, pred);
        }

        if (n % 2 == 0) {
          d = data + begin + (n - 1) * stride;
          T pred = interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride));
          if (quant_pred == true) {
            quant_pred_recover(d-data, *d, pred, offset1, offset2, quant_record);        
          }
          else{
            recover(d - data, *d, pred);
          }
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

#ifdef SZ_ANALYSIS
    my_current_interp_direction = 1;
#endif



    const std::array<int, N> dims = dimension_sequences[direction];

    bool use_begin_cross = (begin[dims[0]] != 0) && (use_cross_block_cubic); 
    bool is_cubic_end_boundary = (end[dims[0]]+1 == global_dimensions[dims[0]]);
    bool is_next_to_bound = (end[dims[0]] + stride2x+1 > global_dimensions[dims[0]]);
    bool is_cubic_end_boundary_ = is_cubic_end_boundary || is_next_to_bound;
    bool use_end_cross = (!is_cubic_end_boundary_) && use_cross_block_cubic;



    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
         j <= end[dims[1]]; j += stride2x) {
      for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
           k <= end[dims[2]]; k += stride2x) {
        size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] +
                              j * dimension_offsets[dims[1]] +
                              k * dimension_offsets[dims[2]];
        bool is_left_boundary = (j==0) || (k==0);
        use_end_cross =  (j%stride2x==0) && (k%stride2x==0) && use_end_cross;
        bool is_cubic_end_boundary = (end[dims[0]] +1 == global_dimensions[dims[0]]) 
                                      ||(end[dims[0]] +1 != global_dimensions[dims[0]] && end[dims[0]] +stride > global_dimensions[dims[0]]-1)
                                      || !(use_cross_block_cubic);
        predict_error += block_interpolation_1d(
            data, begin_offset,
            begin_offset +
                (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
            stride * dimension_offsets[dims[0]], interp_func, pb, quant_pred_on,
            dimension_offsets[dims[1]] * (stride2x),
            dimension_offsets[dims[2]] * (stride2x), is_left_boundary, use_begin_cross, use_end_cross);

      }
    }



#ifdef SZ_ANALYSIS
    my_current_interp_direction = 2;
#endif
    use_begin_cross = (begin[dims[1]] !=0) && (use_cross_block_cubic); 
    is_cubic_end_boundary = (end[dims[1]]+1 == global_dimensions[dims[1]]);
                                    // ||(!use_cross_block_cubic);
    is_next_to_bound = (end[dims[1]] + stride2x+1 > global_dimensions[dims[1]]);
    is_cubic_end_boundary_ = is_cubic_end_boundary || is_next_to_bound;
    use_end_cross = (!is_cubic_end_boundary_) && use_cross_block_cubic;

    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
         i <= end[dims[0]]; i += stride) {
      for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
           k <= end[dims[2]]; k += stride2x) {
        size_t begin_offset = i * dimension_offsets[dims[0]] +
                              begin[dims[1]] * dimension_offsets[dims[1]] +
                              k * dimension_offsets[dims[2]];

        bool is_left_boundary = (i==0) || (k==0); 
        use_end_cross =  (i%stride2x==0) && (k%stride2x==0) && use_end_cross;
        predict_error += block_interpolation_1d(
            data, begin_offset,
            begin_offset +
                (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
            stride * dimension_offsets[dims[1]], interp_func, pb, quant_pred_on,
            dimension_offsets[dims[0]] * (stride),
            dimension_offsets[dims[2]] * (stride2x),is_left_boundary, use_begin_cross, use_end_cross);
      }
    }
  

#ifdef SZ_ANALYSIS
    my_current_interp_direction = 3;
#endif
    use_begin_cross = (begin[dims[2]] !=0) && (use_cross_block_cubic);
    is_cubic_end_boundary = (end[dims[2]]+1 == global_dimensions[dims[2]]);
    is_next_to_bound = (end[dims[2]] + stride2x+1 > global_dimensions[dims[2]]);
    is_cubic_end_boundary_ = is_cubic_end_boundary || is_next_to_bound;
    use_end_cross = (!is_cubic_end_boundary_) && use_cross_block_cubic;


    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
         i <= end[dims[0]]; i += stride) {
      for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0);
           j <= end[dims[1]]; j += stride) {
        size_t begin_offset = i * dimension_offsets[dims[0]] +
                              j * dimension_offsets[dims[1]] +
                              begin[dims[2]] * dimension_offsets[dims[2]];
        bool is_left_boundary = (i==0) || (j==0); // this is for quant prediction.
        use_end_cross =  (i%stride2x==0) && (j%stride2x==0) && use_end_cross;
        predict_error += block_interpolation_1d(
            data, begin_offset,
            begin_offset +
                (end[dims[2]] - begin[dims[2]]) * dimension_offsets[dims[2]],
            stride * dimension_offsets[dims[2]], interp_func, pb, quant_pred_on,
            dimension_offsets[dims[0]] * stride,
            dimension_offsets[dims[1]] * stride, is_left_boundary, use_begin_cross, use_end_cross);
      }
    }
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



  // at the lowest level, we construct a map to store the representitive's
  // quant index this is a bit map 1 for non-zero quant index, 0 otherwise only
  // for slices with

  int interpolation_level = -1;
  uint blocksize;
  int interpolator_id;
  double eb_ratio = 0.5;
  std::vector<std::string> interpolators = {"linear", "cubic"};
  double linear_interp_eb_factor = sqrt(1.5);
  double cubic_interp_eb_factor = 1.2808688457449497;

  // double cubic_interp_eb_factor = std::cbrt(2);
  double cunic_interp_eb_factor_natural= 1.2932517156377563;
  std::array<double,2> eb_factors = {pow(linear_interp_eb_factor,N), 
                                    pow(cunic_interp_eb_factor_natural,N)};
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
  // added for artifact mitigation
  int current_level = 0;



  std::vector<double> level_abs_ebs; 
  std::vector<int> level_interp_ids;


  double eb_reduction_factor = 1.0;

  double current_base_eb;


  T original_max;    // for data range check
  T original_min;    // for data range check
  T original_range;  // for data range check


  bool artifact_sift_on;

  // original data copy
  const T* orig_data_ptr;
  double original_variance;

  // to record post process
  std::vector<char> post_process_tags; // 0, -1 and 1;

  // post_process utils 

  // std::vector<int> aux_quant_inds;
  std::shared_ptr<std::vector<int>> aux_quant_inds_ptr;
  double region_error_control_eb_compensation = 2.0; // for quant prediction 
  bool quant_pred_on = 0;
  int quant_pred_start_level = 3;
  // bool post_process_on = true;  

  // interpolators
  bool use_cross_block_cubic = false;
  bool use_natural_cubic = false;


   

// Analysis utils;
// This is for conditional compilation not comments
// #ifndef SZ_ANALYSIS
// #define SZ_ANALYSIS

std::vector<int> use_cross_end_idx;
std::vector<int> order_idx; 
#ifdef SZ_ANALYSIS
  std::vector<int> my_level;
  std::vector<int> my_quant_inds_copy;
  std::vector<T> my_pred;
  // std::vector<T> my_pred_noise;
  std::vector<T> my_pred_noise;
  std::vector<int> my_interp_direction;
  int my_current_interp_direction = 0;
  std::vector<int> my_compensation_label;

#endif
};
};  // namespace SZ

#endif
