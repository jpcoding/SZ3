#ifndef _SZ_INTERPOLATION_COMPRESSSOR_HPP
#define _SZ_INTERPOLATION_COMPRESSSOR_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include "SZ3/utils/Statistic.hpp"

#include "SZ3/compressor/SZInterpCompressorHelp.hpp"
// #include "SZ3/compressor/SZInterpolation_postprocess.hpp"
#include "SZ3/postprocess/postprocess.hpp"
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
#include "SZ3/utils/critical_points.hpp"
#include "SZ3/utils/interpolation_level.hpp"

namespace SZ {
template <class T, uint N, class Quantizer, class Encoder, class Lossless>
class SZInterpolationCompressor {
 public:
  SZInterpolationCompressor(
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
    read(sift_mode, buffer_pos, remaining_length);
    read(block_flush_on, buffer_pos, remaining_length);
    read(block_sift_on, buffer_pos, remaining_length);

    read(use_stochastic_predict, buffer_pos);
    read(use_stochastic_quantize, buffer_pos);
    read(use_stochastic_decompress, buffer_pos);
    read(use_stochastic_eb, buffer_pos);
    read(normal_mean, buffer_pos);
    read(normal_std, buffer_pos);
    read(uniform_lower, buffer_pos);
    read(bernoulli_p, buffer_pos);
    read(random_seed, buffer_pos);

    read(region_error_control_on, buffer_pos);
    read(region_error_control_mode, buffer_pos);
    read(region_error_control_threshold, buffer_pos);
    read(error_starting_level, buffer_pos);
    read(error_map_stride, buffer_pos);
    read(region_error_control_eb_compensation, buffer_pos);
    read(region_error_control_eb_reduction, buffer_pos);

    if (region_error_control_on == true) {
      //   std::cout << "region_error_control_on = " << region_error_control_on
      //             << std::endl;
      //   std::cout << "region_error_control_mode = " <<
      //   REGION_ERROR_CONTROL_MODE_STR[region_error_control_mode]
      //             << std::endl;

      // size_t error_control_tags_size =0;
      // read(error_control_tags_size,buffer_pos);
      // region_error_tags.resize(error_control_tags_size);
      // convertByteArray2IntArray_fast_1b_sz(error_control_tags_size,
      // buffer_pos,
      //                                     (error_control_tags_size-1)/8+1,
      //                                     region_error_tags.data());
      // writefile("region_error_tags_decompress.dat",
      // region_error_tags.data(),region_error_tags.size()); size_t
      // huffman_size = 0; size_t zstd_size = 0; size_t region_error_tags_size
      // = 0; read(region_error_tags_size, buffer_pos); read(huffman_size,
      // buffer_pos); read(zstd_size, buffer_pos); std::cout <<
      // "region_error_tags_size = " << region_error_tags_size << std::endl;
      // std::cout << "huffman_size = " << huffman_size << std::endl;
      // std::cout << "zstd_size = " << zstd_size << std::endl;
      // std::vector<uchar> compressed_map(zstd_size);
      // read( compressed_map.data(), zstd_size,buffer_pos);
      // // lossless decompress
      // uchar *decompressed_map_data =
      //     lossless.decompress(compressed_map.data(), zstd_size);
      // // huffman decode
      // uchar const *decompressed_map_ptr = decompressed_map_data;

      // auto map_decoder = SZ::HuffmanEncoder<uchar>();
      // map_decoder.load(decompressed_map_ptr, zstd_size);
      // region_error_tags.resize(0);
      // region_error_tags = map_decoder.decode(decompressed_map_ptr,
      // region_error_tags_size); map_decoder.postprocess_decode();
      // lossless.postdecompress_data(decompressed_map_data);
      size_t tag_map_size = 0;
      read(tag_map_size, buffer_pos);
      std::cout << "ptr position = " << remaining_length << std::endl;
      auto map_decoder = SZ::HuffmanEncoder<uchar>();
      map_decoder.load(buffer_pos, remaining_length);
      region_error_tags.resize(0);
      region_error_tags = map_decoder.decode(buffer_pos, tag_map_size);
      std::cout << "compressed size = " << cmpSize << "\n";
      std::cout << "ptr position = " << remaining_length << std::endl;
      map_decoder.postprocess_decode();
    }
    // read auxilliary data
    read(detection_block_size, buffer_pos);
    read(num_detection_block, buffer_pos);
    // std::cout<<"num_detection_block = " << num_detection_block <<std::endl;
    Timer timer;
    timer.start();

    if (block_flush_on == 1) {
      flushed_block_id = std::vector<uchar>(num_detection_block);
      convertByteArray2IntArray_fast_1b_sz(
          num_detection_block, buffer_pos, (num_detection_block - 1) / 8 + 1,
          flushed_block_id.data());
    }
    // flushed_block_id = std::vector<uchar>(num_detection_block);
    // convertByteArray2IntArray_fast_1b_sz(num_detection_block, buffer_pos,
    // (num_detection_block - 1) / 8 + 1, flushed_block_id.data());
    if (block_sift_on == 1) {
      significant_block_id = std::vector<uchar>(num_detection_block);
      convertByteArray2IntArray_fast_1b_sz(
          num_detection_block, buffer_pos, (num_detection_block - 1) / 8 + 1,
          significant_block_id.data());
    }

    // read additional variable
    read(noise_rate, buffer_pos);
    read(detection_eb_rate, buffer_pos);
    read(original_max, buffer_pos);
    read(original_min, buffer_pos);
    // std::cout << "detection_eb_rate = " << detection_eb_rate << std::endl;
    // std::cout << "noise_rate = " << noise_rate << std::endl;

    pred_noise_thresh =
        std::max(std::abs(original_max), std::abs(original_min)) *
        pred_noise_thresh_ratio;
    original_range = original_max - original_min;

    srand(3333);
    init();

    bool post_block_noise_rng_on = false;
    read(post_block_noise_rng_on, buffer_pos);
    std::vector<uchar> post_block_error_tags;
    if (post_block_noise_rng_on == true) {
      // restore the blockinfo here
      std::cout << "post_block_noise_rng_on = " << post_block_noise_rng_on
                << std::endl;
      size_t post_block_error_tags_size = 0;
      read(post_block_error_tags_size, buffer_pos);
      post_block_error_tags.resize(post_block_error_tags_size);
      convertByteArray2IntArray_fast_1b_sz(
          post_block_error_tags_size, buffer_pos,
          (post_block_error_tags_size - 1) / 8 + 1,
          post_block_error_tags.data());
    }

    // timer.stop("read auxilliary data");

    size_t num_flushed_elements = 0;
    if (block_flush_on == 1 || block_sift_on == 1) {
      // num_flushed_elements = compute_auxilliary_data_decompress(decData);
      num_flushed_elements = compute_auxilliary_data_decompress(
          decData, N, global_dimensions.data(), detection_block_size,
          num_elements, block_sift_on, block_flush_on, 0,
          flushed_block_id, flushed_block, significant_block_id,
          significant_block);
    }

    artifact_sift_on = block_sift_on;

    quantizer.load(buffer_pos, remaining_length);
    encoder.load(buffer_pos, remaining_length);
    quant_inds =
        encoder.decode(buffer_pos, num_elements - num_flushed_elements);

    encoder.postprocess_decode();

    lossless.postdecompress_data(buffer);
    double eb_input = quantizer.get_eb();
    double eb_final;
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

    std::cout << "start decompression\n";
    // *decData = quantizer.recover(0, quant_inds[quant_index++]);
    recover(0, *decData, 0);

    std::uniform_int_distribution<int> direction_choice(
        0, dimension_sequences.size() - 1);

    for (uint level = interpolation_level;
         level > 0 && level <= interpolation_level; level--) {
      // if (level >= 3) {
      //     quantizer.set_eb(eb * eb_ratio);
      // } else {
      //     quantizer.set_eb(eb);
      // }

      // direction_sequence_id = direction_choice(mt);
      // std::cout << "direction_sequence_id = " << direction_sequence_id <<
      // std::endl; for (int i = 0; i < N; i++)
      // {
      //     std::cout << dimension_sequences[direction_sequence_id][i] << " ";
      // }
      // std::cout << std::endl;

      // current_level = level;
      // current_base_eb = eb_final;
      // quantizer.set_eb(eb_final);
      // eb_final *= eb_reduction_factor;



      // if (level <= 2) {
      //   quantizer.set_eb(eb_input);
      // }

      size_t stride = 1U << (level - 1);
      auto inter_block_range =
          std::make_shared<SZ::multi_dimensional_range<T, N>>(
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
    quantizer.postdecompress_data();

    // postprocess to remove artifacts
    if (post_block_noise_rng_on == true) {
      // create a blockinfo vector and read the data
      std::cout << "post_block_noise_rng_on = " << post_block_noise_rng_on
                << std::endl;
      double noise_base = noise_rate * eb_input;
      // copy the dimensions to conf.dims.dat
      decompress_opt_reginal_post_process(
          decData, post_block_error_tags, N, global_dimensions.data(),
          detection_block_size, noise_base, mt);
    }
    //            timer.stop("Interpolation Decompress");
    // std::cout << "counter = " << counter << std::endl;
    writefile(
        "quant_inds.decompress.dat", quant_inds.data(), quant_inds.size());

    // Debug use
    writefile(
        "flushed_block_id.decompress.dat", flushed_block_id.data(),
        num_detection_block);
    writefile(
        "significant_block_id.decompress.dat", significant_block_id.data(),
        num_detection_block);
    writefile(
        "flushed_block.decompress.dat", flushed_block.data(), num_elements);
    writefile(
        "significant_block.decompress.dat", significant_block.data(),
        num_elements);

    return decData;
  }

  // compress given the error bound
  uchar *compress(const Config &conf, T *data, size_t &compressed_size)
  {
    operating_data_ptr = data;
    
    std::cout << "original variance = " << original_variance << std::endl;

    std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
    blocksize = conf.interpBlockSize;
    interpolator_id = conf.interpAlgo;
    direction_sequence_id = conf.interpDirection;
    // assign additional variable
    detection_block_size = conf.detection_block_size;
    detection_threshold = conf.detection_threshold;
    detection_eb_rate = conf.detection_eb_rate;
    noise_rate = conf.noise_rate;
    sift_mode = conf.block_sift_mode;
    block_flush_on = conf.block_flush_on;
    block_sift_on = conf.block_sift_on;

    use_stochastic_quantize = conf.use_stochastic_quantize;
    use_stochastic_predict = conf.use_stochastic_predict;
    use_stochastic_decompress = conf.use_stochastic_decompress;
    use_stochastic_eb = conf.use_stochastic_eb;
    normal_std = conf.normal_std;
    normal_mean = conf.normal_mean;
    uniform_lower = conf.uniform_lower;
    bernoulli_p = conf.bernoulli_p;
    random_seed = conf.random_seed;

    region_error_control_on = conf.region_error_on;
    region_error_control_mode = conf.region_error_control_mode;
    region_error_control_threshold = conf.region_error_control_threshold;
    error_starting_level = conf.region_error_control_start_level;
    error_map_stride = conf.regional_error_block_size;
    region_error_control_eb_compensation =
        conf.region_error_control_eb_compensation;
    region_error_control_eb_reduction = conf.region_error_control_eb_reduction;

    if (block_sift_on || block_flush_on)
      std::cout << "detection_block_size = " << detection_block_size
                << std::endl;
    if (block_sift_on)
      std::cout << "detection_threshold = " << detection_threshold
                << std::endl;
    if (block_sift_on)
      std::cout << "detection_eb_rate = " << detection_eb_rate << std::endl;
    if (noise_rate != 0)
      std::cout << "noise_rate = " << noise_rate << std::endl;

    orig_data_copy = (const T *) conf.PASS_DATA.original_data_prt;



    init();

    original_variance = SZ::data_variance(data, num_elements);


    if (region_error_control_on == true) {
      std::cout << "region_error_control_on = " << region_error_control_on
                << std::endl;
      std::cout << "region_error_control_mode = "
                << REGION_ERROR_CONTROL_MODE_STR[region_error_control_mode]
                << std::endl;
      region_error_tags_inorder.resize(num_elements, 0);
    }


    // For data range check.
    auto orig_min_max = std::minmax_element(data, data + num_elements);
    original_min = *orig_min_max.first;
    original_max = *orig_min_max.second;
    std::cout << "original max " << original_max << std::endl;
    std::cout << "original min " << original_min << std::endl;
    pred_noise_thresh =
        std::max(std::abs(original_max), std::abs(original_min)) *
        pred_noise_thresh_ratio;
    original_range = original_max - original_min;

    Timer timer;
    if (block_flush_on == 1 || block_sift_on == 1 ) {
      timer.start();
      compute_auxilliary_data(
          conf, data, detection_block_size, num_elements, num_detection_block,
          flushed_block_id, flushed_block, significant_block_id,
          significant_block);
      timer.stop("Auxilliary Data Compress");
    }

    artifact_sift_on = block_sift_on;
    quant_inds.reserve(num_elements);
    size_t interp_compressed_size = 0;

    rand_collector.resize(num_elements);
    error_recorder.resize(num_elements, 0);

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

    // quant_inds.push_back(quantizer.quantize_and_over*data, 0));
    quantize(0, *data, 0);

    // Timer timer;
    timer.start();

    // double reduction_factor;
    // double real_eb_ratio;
    // if( interpolators[interpolator_id] == "linear")
    // {
    //     reduction_factor = sqrt(27/8);
    // }
    // else
    // {
    //     reduction_factor = sqrt(4.462681);
    // }
    // real_eb_ratio = pow(1/reduction_factor, interpolation_level-1);

    std::uniform_int_distribution<int> direction_choice(
        0, dimension_sequences.size() - 1);

    for (uint level = interpolation_level;
         level > 0 && level <= interpolation_level; level--) {
      // direction_sequence_id = direction_choice(mt);
      // std::cout << "direction_sequence_id = " << direction_sequence_id <<
      // std::endl; for (int i = 0; i < N; i++)
      // {
      //     std::cout << dimension_sequences[direction_sequence_id][i] << " ";
      // }
      // std::cout << std::endl;
      // if (level >= 3) {
      //     quantizer.set_eb(eb_input * eb_ratio);
      // } else {
      //     quantizer.set_eb(eb_input);
      // }
     
      
      current_level = level;
      current_base_eb = eb_final;
      quantizer.set_eb(eb_final);
      eb_final *= eb_reduction_factor;

            // linear increase 
      // current_level = level; 
      // double current_eb = eb_input + (1e-32 - eb_input)/(interpolation_level -1)*(level-1);
      // quantizer.set_eb(current_eb);

      // exponent 
      // current_level = level;
      // double current_eb = eb_input * std::exp(1) * std::exp(-1.0*level);
      // quantizer.set_eb(eb_input/(level*level*level));

      // if(level >1)
      // {
      //   quantizer.set_eb(1e-32*eb_input);
      // }
      // else {
      //   quantizer.set_eb(eb_input);
      // }


      // if (level <= 2) {
      //   quantizer.set_eb(eb_input);
      // }

      size_t stride = 1U << (level - 1);

      // printf("blocksize = %d\n", blocksize);

      auto inter_block_range =
          std::make_shared<SZ::multi_dimensional_range<T, N>>(
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

        // define pred error collector for blocks

        block_interpolation(
            data, block.get_global_index(), end_idx, PB_predict_overwrite,
            interpolators[interpolator_id], direction_sequence_id, stride);
      }
    }

    std::vector<size_t> error_control_map_index;

    assert(quant_inds.size() <= num_elements);

    encoder.preprocess_encode(quant_inds, 0);
    size_t bufferSize = 1.5 * (quantizer.size_est() + encoder.size_est() +
                               sizeof(T) * quant_inds.size());

    uchar *buffer = new uchar[bufferSize];
    uchar *buffer_pos = buffer;

    std::cout << "bufferSize = " << bufferSize << std::endl; 

    write(global_dimensions.data(), N, buffer_pos);
    write(blocksize, buffer_pos);
    write(interpolator_id, buffer_pos);
    write(direction_sequence_id, buffer_pos);
    write(sift_mode, buffer_pos);
    write(block_flush_on, buffer_pos);
    write(block_sift_on, buffer_pos);

    write(use_stochastic_predict, buffer_pos);
    write(use_stochastic_quantize, buffer_pos);
    write(use_stochastic_decompress, buffer_pos);
    write(use_stochastic_eb, buffer_pos);
    write(normal_mean, buffer_pos);
    write(normal_std, buffer_pos);
    write(uniform_lower, buffer_pos);
    write(bernoulli_p, buffer_pos);
    write(random_seed, buffer_pos);

    write(region_error_control_on, buffer_pos);
    write(region_error_control_mode, buffer_pos);
    write(region_error_control_threshold, buffer_pos);
    write(error_starting_level, buffer_pos);
    write(error_map_stride, buffer_pos);
    write(region_error_control_eb_compensation, buffer_pos);
    write(region_error_control_eb_reduction, buffer_pos);

    if (region_error_control_on == true) {
      // create another huffman encoder for the error control map
      auto map_encoder = SZ::HuffmanEncoder<uchar>();
      std::cout << "region_error_control_on = " << region_error_control_on
                << std::endl;
      std::cout << "region_error_control_mode = "
                << REGION_ERROR_CONTROL_MODE_STR[region_error_control_mode]
                << std::endl;
      std::cout << "region_error_tags size" << region_error_tags.size()
                << "\n";
      writefile(
          "region_error_tags_compress.dat", region_error_tags.data(),
          region_error_tags.size());
      write(region_error_tags.size(), buffer_pos);
      // creat a new buffer
      std::vector<uchar> region_error_tags_bitmap;
      region_error_tags_bitmap = region_error_tags;
      // save as 2-bit error map
      // Bitmap_to_Bytes(region_error_tags, region_error_tags_bitmap);
      // std::cout << "region_error_tags_bitmap size" <<
      // region_error_tags_bitmap.size() << "\n";

      // write the map size to the buffer
      // write(region_error_tags_bitmap.size(),buffer_pos);
      // huffman encoding for buffer_error_map
      map_encoder.preprocess_encode(region_error_tags_bitmap, 0);
      map_encoder.save(buffer_pos);
      map_encoder.encode(region_error_tags_bitmap, buffer_pos);
      map_encoder.postprocess_encode();

      //
      //   // huffman encoding for buffer_error_map
      //   uchar *buffer_error_map = new
      //   uchar[region_error_tags_bitmap.size()]; uchar *buffer_error_map_pos
      //   = buffer_error_map;
      //   map_encoder.preprocess_encode(region_error_tags_bitmap, 0);
      //   map_encoder.save(buffer_error_map_pos);
      //   map_encoder.encode(region_error_tags_bitmap, buffer_error_map_pos);
      //   map_encoder.postprocess_encode();
      //   write(region_error_tags.size(), buffer_pos);
      //   write((size_t) (buffer_error_map_pos - buffer_error_map),
      //   buffer_pos);
      //   // write(buffer_error_map, buffer_error_map_pos - buffer_error_map,
      //   buffer_pos);
      //   //compress with zstd on the buffer_error_map
      //   size_t compressed_size = 0;

      //   uchar *lossless_map_data =
      //     lossless.compress(buffer_error_map, buffer_error_map_pos -
      //     buffer_error_map, compressed_size);
      //   std::cout << "compresssed size = " << compressed_size << "\n";

      //   write(compressed_size, buffer_pos);
      //   write(lossless_map_data, compressed_size, buffer_pos);
      //   lossless.postcompress_data(lossless_map_data);
      //   // std::cout << "huffman encoding size = " << buffer_error_map_pos -
      //   buffer_error_map << "\n";
      //   // // use zstd on the buffer_pos_map
      //   // size_t compressed_size = 0;
      //   // uchar *lossless_map_data =
      //   //   lossless.compress(buffer_error_map_pos, buffer_error_map_pos -
      //   buffer_error_map, compressed_size);
      //   // lossless.postcompress_data(buffer_error_map);
      //   // convertIntArray2ByteArray_fast_1b_to_result_sz(
      //   //     region_error_tags.data(), region_error_tags.size(),
      //   buffer_pos);
      //  free(buffer_error_map);
    }

    // add auxilliary array
    write(detection_block_size, buffer_pos);
    // num_detection_block = significant_block_id.size();
    // if(num_detection_block == 0) num_detection_block =
    // flushed_block_id.size();
    write(num_detection_block, buffer_pos);
    if (block_flush_on) {
      std::cout << "writing flushed_block_id" << std::endl;
      convertIntArray2ByteArray_fast_1b_to_result_sz(
          flushed_block_id.data(), flushed_block_id.size(), buffer_pos);
    }
    if (block_sift_on) {
      std::cout << "writing block_sift_id" << std::endl;
      convertIntArray2ByteArray_fast_1b_to_result_sz(
          significant_block_id.data(), significant_block_id.size(),
          buffer_pos);
    }
    // add additional variable
    write(noise_rate, buffer_pos);
    write(detection_eb_rate, buffer_pos);
    write(original_max, buffer_pos);
    write(original_min, buffer_pos);


    // Postprocess to remove artifact
    write(conf.post_block_noise_rng_on, buffer_pos);
    if (conf.post_block_noise_rng_on == true) {
      std::cout << "post_block_noise_rng_on = " << conf.post_block_noise_rng_on
                << std::endl;
      // create a new noise uniform distribution
      double noise_base = noise_rate * eb_input;
      post_block_error_tags = compress_opt_reginal_post_process(
          data, error_recorder.data(), conf, detection_block_size, noise_base,
          mt);
      writefile(
          "post_block_error_tags.dat", post_block_error_tags.data(),
          post_block_error_tags.size());
      size_t post_block_error_tags_size = post_block_error_tags.size();
      write(post_block_error_tags_size, buffer_pos);
      SZ::convertIntArray2ByteArray_fast_1b_to_result_sz(
          post_block_error_tags.data(), post_block_error_tags.size(),
          buffer_pos);
    }

    // store post process post_process_tags
    // create a huffman encoder to encode the post_process_tags 
    if(post_process_tags.size() > 0)
    {
    SZ::HuffmanEncoder<char> post_process_encoder = SZ::HuffmanEncoder<char>();
    post_process_encoder.preprocess_encode(post_process_tags, 0);
    post_process_encoder.save(buffer_pos);
    post_process_encoder.encode(post_process_tags, buffer_pos);
    post_process_encoder.postprocess_encode();
    }
    std::cout << "post_process_tags size = " << post_process_tags.size() << "\n"; 


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
    //            timer.stop("Lossless");

    // post process
    // writefile("compressed.dat", data, num_elements);
    // if(N==3)
    // {
    //   std::cout << "3D post process" << std::endl;
    //   compensation_3d2(
    //   data, my_quant_inds.data(), conf.dims.data(), 0,quantizer.get_eb());
    // }

    writefile("post_compressed.dat", data, num_elements);

    // writefile("rand.dat", rand_collector.data(), num_elements);
    // writefile("pred_noise.dat", my_pred_noise.data(), num_elements);
    // writefile("error.dat", error_recorder.data(), num_elements);
    writefile(
        "error_map.dat", region_error_tags_inorder.data(),
        region_error_tags_inorder.size());
    // writefile("quant_inds.compress.dat", quant_inds.data(),
    // quant_inds.size());

    //   writefile("flushed_block_id.compress.dat", flushed_block_id.data(),
    //           num_detection_block);
    // writefile("significant_block_id.compress.dat",
    //           significant_block_id.data(), num_detection_block);
    // writefile("flushed_block.compress.dat", flushed_block.data(),
    //           num_elements);
    // writefile("significant_block.compress.dat", significant_block.data(),
    // num_elements);

    // std::cout << "rand counter = " << counter << std::endl;

#ifdef SZ_ANALYSIS
    writefile("pred.dat", my_pred.data(), num_elements);
    writefile("quant.dat", my_quant_inds.data(), num_elements);
    writefile("quant_processed.dat", my_quant_inds_copy.data(), num_elements);
    writefile("decompressed.dat", data, num_elements);
    writefile("level.dat", my_level.data(), num_elements);
    writefile(
        "interp_direction.dat", my_interp_direction.data(), num_elements);
    writefile(
        "compensation_label.int32", my_compensation_label.data(),
        my_compensation_label.size());
    std::cout << "[ANALYSIS COMPILATION MODE]" << std::endl;

    // Try huffman encoding on inorder quantization and level-wise quantization
    // integers
    // SZ::HuffmanEncoder<int> huff_coding = SZ::HuffmanEncoder<int>();

    // uchar *test_buffer = new uchar[bufferSize];
    // uchar *test_buffer_pos = test_buffer;
    // huff_coding.preprocess_encode(my_quant_inds, 0);
    // huff_coding.save(test_buffer_pos);
    // huff_coding.encode(my_quant_inds, test_buffer_pos);
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

  // const T* get_decompressed_data_from_compression()
  // {
  //   return operating_data_ptr;
  // }

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
      std::cout << "dimension_offsets[" << i << "] = " << dimension_offsets[i]
                << std::endl;
    }

    dimension_sequences = std::vector<std::array<int, N>>();
    auto sequence = std::array<int, N>();
    for (int i = 0; i < N; i++) { sequence[i] = i; }
    do {
      dimension_sequences.push_back(sequence);
    } while (std::next_permutation(sequence.begin(), sequence.end()));

    // error map construction

    size_t error_map_size =
        (size_t)std::ceil(1.0 * global_dimensions[0] / error_map_stride) *
        std::ceil(1.0 * global_dimensions[1] / error_map_stride) *
        std::ceil(1.0 * global_dimensions[2] / error_map_stride);
    error_map.resize(error_map_size);

    mt = std::mt19937(random_seed);
    mt_uni_int_dist = std::mt19937(10000);
    uniform_dist = std::uniform_real_distribution<T>(0, 1);
    // uniform_dist = std::uniform_real_distribution<T>(0.45, 0.55);

    normal_dist = std::normal_distribution<T>(normal_mean, normal_std);

    bernoulli_dist = std::bernoulli_distribution(bernoulli_p);

    uniform_dist_int32 = std::uniform_int_distribution<int>(-1, 1);
    std::vector<double> randome_predit_weights{0.3, 0.4, 0.3};
    std::discrete_distribution<> discrete_dist_predict(
        randome_predit_weights.begin(), randome_predit_weights.end());

#ifdef SZ_ANALYSIS
    my_level.resize(num_elements);
    my_quant_inds.resize(num_elements);
    my_quant_inds_copy.resize(num_elements);
    my_pred.resize(num_elements);
    my_pred[0] = 0;
    my_level[0] = interpolation_level;
    my_pred_noise.resize(num_elements, 0);
    my_interp_direction.resize(num_elements, 0);
    my_compensation_label.resize(num_elements, 0);
#endif

    interp_level_calculator =
        InterpolationLevel<T>(N, global_dimensions.data());
    if (region_error_control_on == true) {
      error_control_map.resize(num_elements, 0);
    }
    if (N == 3) region_error_tags.reserve(63 / 64 * num_elements);
    if (N == 2) region_error_tags.reserve(15 / 16 * num_elements);
  }

  inline void range_check(T &d)
  {
    if (d > original_max) { d = original_max; }
    if (d < original_min) { d = original_min; }
  }

  inline int backward_compensate_pred(
      size_t idx, size_t offset1, size_t offset2, T &pred, T *data,
      const double compensation)
  {
    // return 0;
    // compensation = ? * eb; 0.5, 1.1, 1.5, 2
    //    0   1  2
    // 0 *A  *B *C
    // 1 *D  *E *F
    // 2 *G  *H  X
    int A = idx - offset1 * 2 - offset2 * 2;
    int B = idx - offset1 * 2 - offset2;
    int C = idx - offset1 * 2;
    int D = idx - offset1 - offset2 * 2;
    int E = idx - offset1 - offset2;
    int F = idx - offset1;
    int G = idx - offset2 * 2;
    int H = idx - offset2;

    // if (my_quant_inds[A] == 0 || my_quant_inds[B] == 0 ||
    //     my_quant_inds[C] == 0 || my_quant_inds[D] == 0 ||
    //     my_quant_inds[E] == 0 || my_quant_inds[F] == 0 ||
    //     my_quant_inds[G] == 0 || my_quant_inds[H] == 0) {
    //   return 0;
    // }

    if (my_quant_inds[E] == 0 || my_quant_inds[F] == 0 ||
        my_quant_inds[H] == 0) {
      return 0;
    }
    int quant_A = my_quant_inds[A] - (1 << 15);
    int quant_B = my_quant_inds[B] - (1 << 15);
    int quant_C = my_quant_inds[C] - (1 << 15);
    int quant_D = my_quant_inds[D] - (1 << 15);
    int quant_E = my_quant_inds[E] - (1 << 15);
    int quant_F = my_quant_inds[F] - (1 << 15);
    int quant_G = my_quant_inds[G] - (1 << 15);
    int quant_H = my_quant_inds[H] - (1 << 15);

    int quant_compensate = 0;
    if (quant_H > 0 && quant_F > 0) {
      quant_compensate = (quant_H + quant_F - quant_E);
      pred += quant_compensate * compensation;
      
    }
    else if (quant_H < 0 && quant_F < 0) {
      quant_compensate = (quant_H + quant_F - quant_E);
      pred += quant_compensate * compensation;
    }
    else {
      return 0;
    }
    return quant_compensate;
  }

  inline int backward_compensate_decompressed(const size_t &idx, const size_t &offset1, const size_t &offset2, 
          const T &current_copy,T *data, 
          const double &tol  , const double &compensation)
  {
    // position
    //H D E 
    //F A B 
    //G C X 
    // X is current idx
    T error = current_copy - data[idx]; 
    if(std::abs(error)<tol) return 0;
    size_t A = idx - offset1 - offset2; // index
    size_t B = idx - offset1; // index 
    size_t C = idx - offset2; // index
    size_t D = idx - offset1 - offset2*2; // index
    size_t E = idx - offset1*2; // index
    size_t F = idx - offset1*2 - offset2; // index
    size_t G = idx - offset2*2; // index
    size_t H = idx - offset2*2 - offset1*2; // index

    int sign = error > 0 ? 1 : -1; 
    data[A] = data[A] + sign * compensation;
    data[B] = data[B] + sign * compensation;
    data[C] = data[C] + sign * compensation;
    data[idx] = data[idx] + sign * compensation;
    data[D] = data[D] + sign * compensation;
    data[E] = data[E] + sign * compensation;
    data[F] = data[F] + sign * compensation;
    data[G] = data[G] + sign * compensation;
    data[H] = data[H] + sign * compensation;
    return sign;
  }

  inline char backward_compensate_decompressed_3d(const size_t &idx, const size_t &offset1, const size_t &offset2, const size_t &offset3,
          const T &current_copy,T *data, 
          const double &tol  , const double &compensation)
  {
    // position
    // suppose offset1 is the along the interpolation direction 
    // offset2 and offset3 are the other two directions 
    // current plane 
    // A B 
    // C X 
    // X is current idx 
    // the previous plane
    // D E
    // F G

    T error = current_copy - data[idx];
    int count_pos = 0;
    int count_neg = 0;
    int thresh = 5;

    size_t A = idx - offset2 - offset3; // index
    size_t B = idx - offset2; // index
    size_t C = idx - offset3; // index

    size_t D = idx - offset1 - offset2 - offset3; // index
    size_t E = idx - offset1 - offset2; // index
    size_t F = idx - offset1 - offset3; // index
    size_t G = idx - offset1; // index

    size_t indexs[8] = {idx, A, B, C, D, E, F, G};

    for(int i = 0; i < 8; i++)
    {
      error = orig_data_copy[indexs[i]] - data[indexs[i]];
      if(std::abs(error) < tol) continue;
      if(error > 0) count_pos++;
      else count_neg++;
    }

    if(count_pos < thresh && count_neg < thresh) return 0;

    if(count_pos >= thresh ) 
    {
      for(int i = 0; i < 8; i++)
      {
        data[indexs[i]] = data[indexs[i]] + compensation;
      }
      return 1;
    }
    if(count_neg >= thresh ) 
    {
      for(int i = 0; i < 8; i++)
      {
        data[indexs[i]] = data[indexs[i]] - compensation;
      }
      return -1;
    }
    return 0;

    // char sign = error > 0 ? 1 : -1;
    // data[A] = data[A] + sign * compensation;
    // data[B] = data[B] + sign * compensation;
    // data[C] = data[C] + sign * compensation;
    // data[D] = data[D] + sign * compensation;
    // data[E] = data[E] + sign * compensation;
    // data[F] = data[F] + sign * compensation;
    // data[G] = data[G] + sign * compensation;
    // data[idx] = data[idx] + sign * compensation;
    // return sign; 
  }
  

  inline void quantize(size_t idx, T &d, T pred)
  {
    T d_copy = d;
    if (block_flush_on == 1 && flushed_block[idx]) { d = 0; }
    else {
      auto default_eb = quantizer.get_eb();
      if (artifact_sift_on == true && (current_level == 1) &&
          significant_block[idx] == 1) {
        quantizer.set_eb(default_eb * detection_eb_rate);
      }
      else {        
        quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
      }
      quantizer.set_eb(default_eb);
    }
    // range_check(d);

#ifdef SZ_ANALYSIS
    my_level[idx] = current_level;
    my_quant_inds[idx] = quant_inds.back();
    my_quant_inds_copy[idx] = quant_inds.back();
    // my_pred[idx] = pred;
    my_interp_direction[idx] = my_current_interp_direction;
#endif

    // recorder error
    error_recorder[idx] = d_copy - d;
  }

  inline void recover(size_t idx, T &d, T pred)
  {
    if (block_flush_on == 1 && flushed_block[idx]) { d = 0; }
    else {
      auto default_eb = quantizer.get_eb();
      if (artifact_sift_on == true && (current_level == 1) &&
          significant_block[idx] == 1) {
        quantizer.set_eb(default_eb * detection_eb_rate);
      }
      if (use_stochastic_quantize == 1 && current_level >= 1) {
        d = quantizer.recover(pred, quant_inds[quant_index++]);
      }
      else {
        d = quantizer.recover(pred, quant_inds[quant_index++]);
      }
      quantizer.set_eb(default_eb);
    }
    // range_check(d);
  };

  double block_interpolation_1d(
      T *data, size_t begin, size_t end, size_t stride,
      const std::string &interp_func, const PredictorBehavior pb,
      bool error_tune = false, size_t offset1 = 0, size_t offset2 = 0, bool high_level_predict = false,
      bool use_fft_interp = false)
  {
    size_t n = (end - begin) / stride + 1;
    if (n <= 1) { return 0; }
    double predict_error = 0;

    int quant_compensation = 0;

    size_t stride3x = 3 * stride;
    size_t stride5x = 5 * stride;
    region_error_counter = 0;

    // if(use_fft_interp = true && current_level == 1 && n > 5)
    // {
    //   T *d = data + begin; 
    //   size_t interplator_offset = stride; 
    //   interp_fft_fftw( d, stride, 1, n+1, quantizer, quant_inds);
    //   if (n % 2 == 0) {
    //       T *d = data + begin + (n - 1) * stride;
    //       T pred = interp_linear(*(d - stride), *(d - stride));
    //       quantize(d - data, *d, pred);
    //   }

    //   return 0;

    // }
    if (interp_func == "linear" || n < 5) {
      if (pb == PB_predict_overwrite) {
        for (size_t i = 1; i + 1 < n; i += 2) {
          T *d = data + begin + i * stride;
          T d_copy = *d;
          T pred = interp_linear(*(d - stride), *(d + stride));
          my_pred[d - data] = pred;
          if (error_tune == 1 && current_level <= 3) {
            double tol = quantizer.get_eb() / linear_interp_eb_factor;
            double compensation =
                region_error_control_eb_compensation * quantizer.get_eb();
            quant_compensation = backward_compensate_pred(
                d - data, offset1, offset2, pred, data, compensation);
            // backward_compensate(d - stride-data, d - stride-data, pred,
            // error_recorder.data(), compensation, tol);
            
          }
          quantize(d - data, *d, pred);
          my_quant_inds[d - data] = quant_inds.back() + quant_compensation;
          my_compensation_label[d - data] = quant_compensation;

          // quant_inds.back() -= x;

          // quantize(d - data, *d, interp_linear(*(d - stride), *(d +
          // stride)));
        }
        if (n % 2 == 0) {
          T *d = data + begin + (n - 1) * stride;
          // quantize(d - data, *d, *(d - stride));
          T pred = interp_linear(*(d - stride), *(d - stride));
          my_pred[d - data] = pred;
          if (error_tune == 1 && current_level <= 3) {
            double tol = quantizer.get_eb() / linear_interp_eb_factor;
            double compensation =
                region_error_control_eb_compensation * quantizer.get_eb();
            quant_compensation = backward_compensate_pred(
                d - data, offset1, offset2, pred, data, compensation);
            // backward_compensate(d - stride-data, d - stride-data, pred,
            // error_recorder.data(), compensation, tol);
          }
          quantize(d - data, *d, pred);
          my_quant_inds[d - data] = quant_inds.back() + quant_compensation;
          my_compensation_label[d - data] = quant_compensation;
          // quantize(d - data, *d, *(d - stride));

          // if (n < 4) {
          //     quantize(d - data, *d, *(d - stride));
          // } else {
          //     quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d -
          //     stride)));
          // }
        }
      }
      else {
        for (size_t i = 1; i + 1 < n; i += 2) {
          T *d = data + begin + i * stride;
          recover(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
        }
        if (n % 2 == 0) {
          T *d = data + begin + (n - 1) * stride;
          recover(d - data, *d, *(d - stride));
          // if (n < 4) {
          //     recover(d - data, *d, *(d - stride));
          // } else {
          //     recover(d - data, *d, interp_linear1(*(d - stride3x), *(d -
          //     stride)));
          // }
        }
      }
    }
    else {
      if (pb == PB_predict_overwrite) {
        double dafault_eb = quantizer.get_eb();
        T *d;
        size_t i;
        d = data + begin + stride;
        if( current_level ==1 && 0 )
        {// check range of this line 
          // check the variance of this line 

            T max = *d;
            T min = *d;
            // T second_max = *d;
            // T second_min = *d;
            // // find the max,sencond_max, min, second_min
            // for (i = 0; i < n; i+=2) {
            //   d = data + begin + i * stride;
            //   if(*d > max) max = *d;
            //   if(*d < min) min = *d;
            //   if(*d > second_max && *d < max) second_max = *d;
            //   if(*d < second_min && *d > min) second_min = *d;
            // }

            // T range = max - min;
            // T rnage2 = second_max - second_min;

            // std::cout << "range1 = " << range << "\nrange2 = " << rnage2 << std::endl;
            // if(range > 1.5*rnage2)
            // {
            //   quantizer.set_eb(dafault_eb*0.5);
            // }


            double mean = 0;
            int count = 0;
            for (i = 0; i < n; i+=2) {
              d = data + begin + i * stride;
              if(*d > max) max = *d;
              if(*d < min) min = *d;
              mean += *d;
              count++;
            }
            mean = mean/count;


            T range = max - min; 
            double  variance1 = 0;
            for (i = 0; i < n; i+=2) {
              d = data + begin + i * stride;
              variance1 += (*d - mean)*(*d - mean);
            }
            variance1 = variance1/(count-1);
            if(variance1 < dafault_eb*dafault_eb/original_range/original_range)
            {

            }

            else {
            // remove max and min and calculate another variance
            mean = (mean*count- max - min)/(count-2);
            double variance2 = 0;
            for (i = 0; i < n; i+=2) {
              d = data + begin + i * stride;
              variance2 += (*d - mean)*(*d - mean);
            }
            variance2 -= (max - mean)*(max - mean);
            variance2 -= (min - mean)*(min - mean);
            variance2 = variance2/(count-3);
            // std::cout << "variance1 = " << variance1 << "\nvariance2 = "<< variance2 << std::endl;

            if( variance1 > 1.1*variance2)
            {
              quantizer.set_eb(dafault_eb*0.6);
            }
            }
        }



        

        for (i = 3; i + 3 < n; i += 2) {
          d = data + begin + i * stride;
          // quantize(
          //     d - data, *d,
          //     interp_cubic(
          //         *(d - stride3x), *(d - stride), *(d + stride),
          //         *(d + stride3x)));
          T pred = interp_cubic(
              *(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x));
           
          my_pred[d - data] = pred;
          if (error_tune == 1 && current_level <= 3) {
            double tol = quantizer.get_eb() / linear_interp_eb_factor;
            double compensation =
                region_error_control_eb_compensation * quantizer.get_eb();
            quant_compensation = backward_compensate_pred(
                d - data, offset1, offset2, pred, data, compensation);

            // quant_compensation = quant_inds.back();
          }
          T d_copy = *d;  
          quantize(d - data, *d, pred);
          my_quant_inds[d - data] = quant_inds.back() + quant_compensation;
          my_compensation_label[d - data] = quant_compensation;
          
          
          // if(current_level ==1 && high_level_predict == true)
          // {
          //   bool quant_condition =  my_quant_inds[d - data] == quantizer.get_radius()
          //   && my_quant_inds[d - data-offset1 ] == quantizer.get_radius()
          //   && my_quant_inds[d - data-offset2] == quantizer.get_radius()
          //   && my_quant_inds[d - data - offset1-offset2] == quantizer.get_radius();

          //   quant_condition = quant_condition
          //                   && my_quant_inds[d - data - 2* offset1] == quantizer.get_radius()
          //                   && my_quant_inds[d - data -2* offset2] == quantizer.get_radius()
          //                   && my_quant_inds[d - data - 2* offset1 - offset2] == quantizer.get_radius()
          //                   && my_quant_inds[d - data - 2* offset2 - offset1] == quantizer.get_radius()
          //                   && my_quant_inds[d - data - 2* offset1 - 2*offset2] == quantizer.get_radius();
          
          //   if(quant_condition == true)
          //   {char tmp_tag = (char) backward_compensate_decompressed(d-data, offset1,  offset2, 
          //             d_copy, data, quantizer.get_eb()*0.5, 0.5*quantizer.get_eb());
          //   post_process_tags.push_back(tmp_tag);
          //   }
          // }
          // if((i+1)%4 ==0 && high_level_predict == true && current_level ==1)
          // {
          //   size_t stride2x = 2*stride; 
          //   char tmp_tag = (char) backward_compensate_decompressed_3d(d-data, stride2x, offset1,  offset2, 
          //             d_copy, data, quantizer.get_eb()*0.5, 0.5*quantizer.get_eb());
          //   post_process_tags.push_back(tmp_tag);
          // }
        }
        d = data + begin + stride;
        quantize(
            d - data, *d,
            interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
        my_pred[d - data] =
            interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x));

        d = data + begin + i * stride;
        quantize(
            d - data, *d,
            interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
        my_pred[d - data] =
            interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride));
        if (n % 2 == 0) {
          d = data + begin + (n - 1) * stride;
          // quantize(d - data, *d, *(d - stride));
          T pred = *(d - stride);
          my_pred[d - data] = pred;
          if (error_tune == 1 && current_level <= 3) {
            double tol = quantizer.get_eb() / linear_interp_eb_factor;
            double compensation =
                region_error_control_eb_compensation * quantizer.get_eb();
            quant_compensation = backward_compensate_pred(
                d - data, offset1, offset2, pred, data, compensation);
            // backward_compensate(d - stride-data, d - stride-data, pred,
            // error_recorder.data(), compensation, tol);
          }
          quantize(d - data, *d, pred);
          my_quant_inds[d - data] = quant_inds.back() + quant_compensation;
          my_compensation_label[d - data] = quant_compensation;

          quantizer.set_eb(dafault_eb);

          // quantize(d - data, *d, *(d - stride));
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
                  *(d + stride3x)));
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
          recover(d - data, *d, *(d - stride));
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
    // quantizer.set_eb(default_eb * c2);

#ifdef SZ_ANALYSIS
    my_current_interp_direction = 1;
#endif

    const std::array<int, N> dims = dimension_sequences[direction];
    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
         j <= end[dims[1]]; j += stride2x) {
      for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
           k <= end[dims[2]]; k += stride2x) {
        size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] +
                              j * dimension_offsets[dims[1]] +
                              k * dimension_offsets[dims[2]];
        // if( k%stride2x ==0 || j%stride2x ==0)
        // {
        //   quantizer.set_eb(default_eb*0.8);
        // }
        bool high_level_predict = (j%4==0) && (k%4==0);
        high_level_predict = 0;

        predict_error += block_interpolation_1d(
            data, begin_offset,
            begin_offset +
                (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
            stride * dimension_offsets[dims[0]], interp_func, pb, 1,
            dimension_offsets[dims[1]] * (stride2x),
            dimension_offsets[dims[2]] * (stride2x), high_level_predict);
        quantizer.set_eb(default_eb);
      }
    }



    // quantizer.set_eb(default_eb * c1);
    // if(current_level ==1) quantizer.set_eb(default_eb*0.5);
#ifdef SZ_ANALYSIS
    my_current_interp_direction = 2;
#endif
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

        bool high_level_predict = (i%2==0) && (k%4==0);
        high_level_predict = 0;
        predict_error += block_interpolation_1d(
            data, begin_offset,
            begin_offset +
                (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
            stride * dimension_offsets[dims[1]], interp_func, pb, 1,
            dimension_offsets[dims[0]] * (stride),
            dimension_offsets[dims[2]] * (stride2x),high_level_predict);
        quantizer.set_eb(default_eb);
      }
    }
  

    // if(current_level ==1) quantizer.set_eb(0.5*default_eb);
    // if(current_level ==1)
    // {
    //   compensation_3d(data,my_quant_inds.data(),
    //         begin.data(), end.data(), dims.data(),
    //         dimension_offsets.data(),
    //         1,stride, 
    //         0, 2, stride, stride2x, quantizer.get_eb(),
    //         quantizer.get_radius()*0.5);
    // } 

#ifdef SZ_ANALYSIS
    my_current_interp_direction = 3;
#endif
    // if(current_level ==1) quantizer.set_eb(default_eb*0.5);
    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
         i <= end[dims[0]]; i += stride) {
      // ) std::cout <<"i = " << i << std::endl;
      for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0);
           j <= end[dims[1]]; j += stride) {
        // if (i==0 && current_level==1) std::cout <<"j = " << j << std::endl;
        size_t begin_offset = i * dimension_offsets[dims[0]] +
                              j * dimension_offsets[dims[1]] +
                              begin[dims[2]] * dimension_offsets[dims[2]];
        // std::vector<size_t> cords =
        // interp_level_calculator.get_coordinates(begin_offset);
        // std::cout<<"cords = "<< cords[0] << " " << cords[1] << " " <<
        // cords[2] << std::endl;
        // if( i%stride2x ==1 && j%stride2x ==1)
        // {
        //   quantizer.set_eb(default_eb*0.65);
        // }
        bool high_level_predict = (i%2==0) && (j%2==0);
        high_level_predict = 0;
        bool use_fft_interp = false ; 
        predict_error += block_interpolation_1d(
            data, begin_offset,
            begin_offset +
                (end[dims[2]] - begin[dims[2]]) * dimension_offsets[dims[2]],
            stride * dimension_offsets[dims[2]], interp_func, pb, 1,
            dimension_offsets[dims[0]] * stride,
            dimension_offsets[dims[1]] * stride, high_level_predict, use_fft_interp);
        quantizer.set_eb(default_eb);
        // if this is decompress and the level 1, we need to special treatment
      }
    }
    quantizer.set_eb(default_eb);


    if(current_level ==1)
    {writefile("compressed.dat", data, num_elements);
     writefile("quant_inds_pre.dat", my_quant_inds.data(), num_elements);
      }



    if(current_level ==1)
    {
      compensation_3d(data,my_quant_inds.data(),
            begin.data(), end.data(), dims,
            dimension_offsets,
            0, stride, 
            1, 2, stride2x, stride2x, quantizer.get_eb(),
            quantizer.get_radius());
    }


    if(current_level ==1)
    {

      compensation_3d(data,my_quant_inds.data(),
            begin.data(), end.data(), dims,
            dimension_offsets,
            1,stride, 
            0, 2, stride, stride2x, quantizer.get_eb(),
            quantizer.get_radius());
    } 
    if(current_level ==1)
    {
      compensation_3d(data,my_quant_inds.data(),
            begin.data(), end.data(), dims,
            dimension_offsets,
            2,stride, 
            0, 1, stride, stride, quantizer.get_eb(),
            quantizer.get_radius());
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

  std::vector<uchar> int2uchar(int *int_data, const size_t datasize)
  {
    std::vector<uchar> uchar_data(datasize, 0);
    for (size_t i = 0; i < datasize; i++) {
      uchar_data[i] =
          (int_data[i] == 4 || int_data[i] == -4 || int_data[i] == -6 ||
           int_data[i] == 6);
    }
    return uchar_data;
  }

  // at the lowest level, we construct a map to store the representitive's
  // quant index this is a bit map 1 for non-zero quant index, 0 otherwise only
  // for slices with

  int interpolation_level = -1;
  uint blocksize;
  int interpolator_id;
  double eb_ratio = 0.5;
  std::vector<std::string> interpolators = {"linear", "cubic"};
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
  // double c = sqrt(4.4159889); // deprecated
  // double c1 = 1.0 / sqrt(1.640625); //deprecated
  // double c2 = 1.0 / 1.640625; //deprecated
  // double c3 = 1.0 / sqrt(4.4159889); //deprecated

  double linear_interp_eb_factor = sqrt(1.5);
  double cubic_interp_eb_factor = 1.2808688457449497;

  int detection_block_size = 4;
  double detection_threshold = 0.9;
  double detection_eb_rate;
  double noise_rate = 0;
  double eb_reduction_factor = 1.0;

  std::vector<uchar> flushed_block;
  std::vector<uchar> flushed_block_id;
  std::vector<uchar> significant_block;     // per datapoint
  std::vector<uchar> significant_block_id;  // per block
  uint8_t sift_mode = SZ::BLOCK_SIFT_MODE::RANGE;
  double current_base_eb;
  bool block_flush_on;
  bool block_sift_on;
  size_t num_detection_block = 0;

  T original_max;    // for data range check
  T original_min;    // for data range check
  T original_range;  // for data range check
  T pred_noise_thresh;
  double pred_noise_thresh_ratio = 1e-4;

  bool artifact_sift_on;

  bool use_stochastic_quantize = false;
  bool use_stochastic_predict = false;
  bool use_stochastic_decompress = false;
  bool use_stochastic_eb = false;

  std::mt19937 mt;

  std::mt19937 mt_uni_int_dist;
  int random_seed = 2333;
  std::uniform_real_distribution<T> uniform_dist;
  std::normal_distribution<T> normal_dist;
  std::bernoulli_distribution bernoulli_dist;
  std::uniform_real_distribution<T> uniform_dist2;
  std::uniform_int_distribution<> uniform_dist_int32;
  std::discrete_distribution<> discrete_dist_predict;

  float normal_std = 1.0;
  float normal_mean = 0.0;
  float uniform_lower = 0.0;
  float bernoulli_p = 0;

  std::vector<T> rand_collector;
  int counter = 0;

  // use to choose the index
  std::uniform_int_distribution<int> uni_choose =
      std::uniform_int_distribution<int>(0, 1);

  // Postprocess varables
  std::vector<uchar> post_block_error_tags;
  std::vector<T> error_recorder;

  // region error tunning variables

  bool region_error_control_on = false;
  uint8_t region_error_control_mode = SZ::REGION_ERROR_CONTROL_MODE::REDUCE_EB;
  int error_starting_level = 3;
  std::vector<uchar> region_error_tags;
  std::vector<uchar> region_error_tags_inorder;
  size_t region_error_tags_index = 0;
  int region_error_counter = 0;
  double region_error_control_threshold = 0.5;
  double region_error_control_eb_compensation = 1.0;
  double region_error_control_eb_reduction = 0.5;
  int regional_error_block_size = 2;
  std::vector<int> error_map;
  int error_map_stride = 4;

  std::vector<size_t> error_control_index;
  std::vector<uchar> error_control_stride;
  std::vector<uchar> error_control_map;  // 0 for no control, 1 for control
  // std::vector<uchar> slice_error_tags;

  InterpolationLevel<T> interp_level_calculator;
  T *operating_data_ptr = nullptr;

  // original data copy
  const T* orig_data_copy;
  double original_variance;

  // to record post process
  std::vector<char> post_process_tags; // 0, -1 and 1;

// Analysis utils;
// This is for conditional compilation not comments
// #ifndef SZ_ANALYSIS
// #define SZ_ANALYSIS
#ifdef SZ_ANALYSIS
  std::vector<int> my_level;
  std::vector<int> my_quant_inds;
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
