#ifndef SZ3_SZINTERP_HPP
#define SZ3_SZINTERP_HPP

#include "SZ3/compressor/SZInterpolationCompressor.hpp"
#include "SZ3/compressor/deprecated/SZBlockInterpolationCompressor.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/api/impl/SZLorenzoReg.hpp"
#include "SZ3/utils/Extraction_level.hpp"
#include "SZ3/utils/Metrics.hpp"
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <memory>


template<class T, SZ::uint N>
char *SZ_compress_Interp(SZ::Config &conf, T *data, size_t &outSize) {


    assert(N == conf.N);
    assert(conf.cmprAlgo == SZ::ALGO_INTERP);
    SZ::calAbsErrorBound(conf, data);

    auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            SZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
    char *cmpData = (char *) sz.compress(conf, data, outSize);
    conf.PASS_DATA.aux_quant_inds_ptr = sz.get_aux_quant_inds_ptr();
    // std::cout <<"ptr address = " << sz.get_aux_quant_inds_ptr() << std::endl;
    // std::cout <<"ptr address = " << conf.PASS_DATA.aux_quant_inds_ptr << std::endl;
    return cmpData;
}


template<class T, SZ::uint N>
void SZ_decompress_Interp(SZ::Config &conf, char *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == SZ::ALGO_INTERP);
    SZ::uchar const *cmpDataPos = (SZ::uchar *) cmpData;
    auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            SZ::LinearQuantizer<T>(),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
    sz.decompress(cmpDataPos, cmpSize, decData);
    conf.PASS_DATA.aux_quant_inds_ptr = sz.get_aux_quant_inds_ptr();
}


template<class T, SZ::uint N>
double do_not_use_this_interp_compress_block_test(T *data, std::vector<size_t> dims, size_t num,
                                                  double eb, int interp_op, int direction_op, int block_size) {

    std::vector<T> data1(data, data + num);
    size_t outSize = 0;

    SZ::Config conf;
    conf.absErrorBound = eb;
    conf.setDims(dims.begin(), dims.end());
    conf.blockSize = block_size;
    conf.interpAlgo = interp_op;
    conf.interpDirection = direction_op;
    auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            SZ::LinearQuantizer<T>(eb),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
    char *cmpData = (char *) sz.compress(conf, data1.data(), outSize, true);
    delete[]cmpData;
    auto compression_ratio = num * sizeof(T) * 1.0 / outSize;
    return compression_ratio;
}

template<class T, SZ::uint N>
char *SZ_compress_Interp_lorenzo(SZ::Config &conf, T *data, size_t &outSize) {
    assert(conf.cmprAlgo == SZ::ALGO_INTERP_LORENZO);

    SZ::Timer timer(true);

    SZ::calAbsErrorBound(conf, data);

    size_t sampling_num, sampling_block;
    std::vector<size_t> sample_dims(N);
    std::vector<T> sampling_data = SZ::sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);

    double best_lorenzo_ratio = 0, best_interp_ratio = 0, ratio;
    size_t sampleOutSize;
    char *cmprData;
    SZ::Config lorenzo_config = conf;
    {
        //test lorenzo
        lorenzo_config.cmprAlgo = SZ::ALGO_LORENZO_REG;
        lorenzo_config.setDims(sample_dims.begin(), sample_dims.end());
        lorenzo_config.lorenzo = true;
        lorenzo_config.lorenzo2 = true;
        lorenzo_config.regression = false;
        lorenzo_config.regression2 = false;
        lorenzo_config.openmp = false;
        lorenzo_config.blockSize = 5;
//        lorenzo_config.quantbinCnt = 65536 * 2;
        std::vector<T> data1(sampling_data);
        cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, data1.data(), sampleOutSize);
        delete[]cmprData;
//    printf("Lorenzo ratio = %.2f\n", ratio);
        best_lorenzo_ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
    }

    {
        //tune interp which first? 
        for (auto &interp_op: {SZ::INTERP_ALGO_LINEAR, SZ::INTERP_ALGO_CUBIC}) {
            ratio = do_not_use_this_interp_compress_block_test<T, N>(sampling_data.data(), sample_dims, sampling_num, conf.absErrorBound,
                                                                     interp_op, conf.interpDirection, sampling_block);
            if (ratio > best_interp_ratio) {
                best_interp_ratio = ratio;
                conf.interpAlgo = interp_op;
            }
        }

        int direction_op = SZ::factorial((int)N) - 1;
        ratio = do_not_use_this_interp_compress_block_test<T, N>(sampling_data.data(), sample_dims, sampling_num, conf.absErrorBound,
                                                                     conf.interpAlgo, direction_op, sampling_block);
        if (ratio > best_interp_ratio * 1.02) {
                best_interp_ratio = ratio;
                conf.interpDirection = direction_op;
        }

        std::cout << "choose direction = "<< direction_op << std::endl;
        std::cout << "choose interp = "<< SZ::INTERP_ALGO_STR[conf.interpAlgo] << std::endl;

    }

    bool useInterp = !(best_lorenzo_ratio > best_interp_ratio && best_lorenzo_ratio < 80 && best_interp_ratio < 80);

    if (useInterp) {
        conf.cmprAlgo = SZ::ALGO_INTERP;
        double tuning_time = timer.stop();
        return SZ_compress_Interp<T, N>(conf, data, outSize);
    } else {
        //further tune lorenzo
        if (N == 3) {
            float pred_freq, mean_freq;
            T mean_guess;
            lorenzo_config.quantbinCnt = SZ::optimize_quant_invl_3d<T>(data, conf.dims[0], conf.dims[1], conf.dims[2],
                                                                       conf.absErrorBound, pred_freq, mean_freq, mean_guess);
            lorenzo_config.pred_dim = 2;
            cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), sampleOutSize);
            delete[]cmprData;
            ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
            if (ratio > best_lorenzo_ratio * 1.02) {
                best_lorenzo_ratio = ratio;
            } else {
                lorenzo_config.pred_dim = 3;
            }
        }

        if (conf.relErrorBound < 1.01e-6 && best_lorenzo_ratio > 5 && lorenzo_config.quantbinCnt != 16384) {
            auto quant_num = lorenzo_config.quantbinCnt;
            lorenzo_config.quantbinCnt = 16384;
            cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), sampleOutSize);
            delete[]cmprData;
            ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
            if (ratio > best_lorenzo_ratio * 1.02) {
                best_lorenzo_ratio = ratio;
            } else {
                lorenzo_config.quantbinCnt = quant_num;
            }
        }
        lorenzo_config.setDims(conf.dims.begin(), conf.dims.end());
        conf = lorenzo_config;
        double tuning_time = timer.stop();
        return SZ_compress_LorenzoReg<T, N>(conf, data, outSize);
    }


}




template<class T, SZ::uint N>
char *SZ_compress_Interp___(SZ::Config &conf, T *data, size_t &outSize, bool opt = 1)
{
    assert(conf.cmprAlgo == SZ::ALGO_INTERP_LORENZO);

    SZ::Timer timer(true);

    SZ::calAbsErrorBound(conf, data);
    
    // how to sample data ?
    // sample data with stride
    // std::vector<T> sampling_data = SZ::sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);
    size_t sampleOutSize =0;
    char *cmprData;

    SZ::Config conf_copy = conf; // copy config file
    int tunning_level = conf.region_error_control_start_level;
    int tunning_stride = (int)1 << (tunning_level-1);

    // extract sample data
    auto level_extractor = SZ::ExtractionLevel<T, N>(conf.dims.data(), tunning_stride, data);
    std::vector<T> sample_data = level_extractor.get_extracted_data();
    std::vector<size_t> sample_dims_ = level_extractor.get_extracted_dims();
    size_t sample_num = sample_data.size();
    conf_copy.blockSize = conf.blockSize / tunning_stride;
    conf_copy.setDims(sample_dims_.begin(), sample_dims_.end());
    conf_copy.region_error_control_start_level = 2;
    // make sure the highest level is larger than 2
    if(level_extractor.get_highest_level()<=2)
    {
        std::cout << "highest level is smaller than 2, not enough data\n";
        exit(1);
    }


    int best_pred_error_thresh_idx = 0;
    double best_pred_error_thresh_ratio = 0;
    double best_ratio = 0;
    double best_indicator = 0;
    if(conf_copy.region_error_control_mode == SZ::REDUCE_EB)
    {
        std::vector<double> pred_error_thresh = {0.1, 0.2, 0.3, 0.4, 0.5,0.6,0.7};
        std::vector<double> reduction_ratio = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
        for(auto & thresh: pred_error_thresh)
        {
            conf_copy.region_error_control_threshold = thresh;
            std::vector<T> sample_data_copy(sample_data);
            char *sample_cmprData = SZ_compress_Interp<T, N>(conf_copy, sample_data_copy.data(), sampleOutSize);
            delete[]sample_cmprData;
            double current_ratio = sample_num * 1.0 * sizeof(T) / sampleOutSize;
            if(current_ratio> best_ratio)
            {
                best_ratio = current_ratio;
                best_pred_error_thresh_ratio = thresh;
            }

        }
        conf.region_error_control_threshold = best_pred_error_thresh_ratio;
    }
    else if (conf_copy.region_error_control_mode == SZ::COMPENSATE_EB)
    {
        // TODO: compensate error bound
        std::vector<double> pred_error_thresh = {0.1, 0.2, 0.3, 0.4, 0.5,0.6,0.7};
        std::vector<double> compensate_ratio = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}; //??
        double best_compensate_ratio = 0;
        // for(auto & thresh: pred_error_thresh) 
        for (auto & compensate : compensate_ratio)
        {
            conf_copy.region_error_control_eb_compensation = compensate;
            // conf_copy.region_error_control_threshold = thresh;
            std::vector<T> sample_data_copy(sample_data);
            char *sample_cmprData = SZ_compress_Interp<T, N>(conf_copy, sample_data_copy.data(), sampleOutSize);
            delete[]sample_cmprData;
            double current_psnr = SZ::METRICS::CALCULATE_PSNR<T>(sample_data.data(), sample_data_copy.data(), sample_num);
            double current_ratio = sample_num * 1.0 * sizeof(T) / sampleOutSize;
            if(current_psnr> best_indicator)
            {
                best_indicator = current_psnr;
                best_compensate_ratio = compensate;
            }
        }
        conf.region_error_control_eb_compensation = best_compensate_ratio;
        std::cout << "best compensate ratio: " << best_compensate_ratio << std::endl;

    }


    return SZ_compress_Interp<T, N>(conf, data, outSize);
}

#endif
