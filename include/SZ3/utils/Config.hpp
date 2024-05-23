//
// Created by Kai Zhao on 4/28/20.
//

#ifndef SZ_Config_HPP
#define SZ_Config_HPP

#include <iostream>
#include <memory>
#include <vector>
#include <numeric>
#include "SZ3/def.hpp"
#include "MemoryUtil.hpp"
#include "SZ3/utils/inih/INIReader.h"

namespace SZ {

    enum EB {
        EB_ABS, EB_REL, EB_PSNR, EB_L2NORM, EB_ABS_AND_REL, EB_ABS_OR_REL
    };
    constexpr const char *EB_STR[] = {"ABS", "REL", "PSNR", "NORM", "ABS_AND_REL", "ABS_OR_REL"};
    constexpr EB EB_OPTIONS[] = {EB_ABS, EB_REL, EB_PSNR, EB_L2NORM, EB_ABS_AND_REL, EB_ABS_OR_REL};

    enum ALGO {
        ALGO_LORENZO_REG, ALGO_INTERP_LORENZO, ALGO_INTERP
    };
    constexpr const char *ALGO_STR[] = {"ALGO_LORENZO_REG", "ALGO_INTERP_LORENZO", "ALGO_INTERP"};
    constexpr const ALGO ALGO_OPTIONS[] = {ALGO_LORENZO_REG, ALGO_INTERP_LORENZO, ALGO_INTERP};

    enum INTERP_ALGO {
        INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC
    };
    constexpr const char *INTERP_ALGO_STR[] = {"INTERP_ALGO_LINEAR", "INTERP_ALGO_CUBIC"};
    constexpr INTERP_ALGO INTERP_ALGO_OPTIONS[] = {INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC};

    enum BLOCK_SIFT_MODE{
        VARIANCE, RANGE, BLOCK_MAX,ISOVALUE
    };

    enum REGION_ERROR_CONTROL_MODE{
        REDUCE_EB, COMPENSATE_EB
    };
    constexpr const char *REGION_ERROR_CONTROL_MODE_STR[] = {"REDUCE_EB", "COMPENSATE_EB"};
    constexpr REGION_ERROR_CONTROL_MODE REGION_ERROR_CONTROL_MODE_OPTIONS[] = {REDUCE_EB, COMPENSATE_EB};

    constexpr const char *BLOCK_SIFT_MODE_STR[] = {"VARIANCE", "RANGE", "BLOCK_MAX", "ISOVALUE"};
    constexpr BLOCK_SIFT_MODE BLOCK_SIFT_MODE_OPTIONS[] = {VARIANCE, RANGE, BLOCK_MAX, ISOVALUE};



    template<class T>
    const char *enum2Str(T e) {
        if (std::is_same<T, ALGO>::value) {
            return ALGO_STR[e];
        } else if (std::is_same<T, INTERP_ALGO>::value) {
            return INTERP_ALGO_STR[e];
        } else if (std::is_same<T, EB>::value) {
            return EB_STR[e];
        } else {
            printf("invalid enum type for enum2Str()\n ");
            exit(0);
        }
    }

    class Config {
    public:
        template<class ... Dims>
        Config(Dims ... args) {
            dims = std::vector<size_t>{static_cast<size_t>(std::forward<Dims>(args))...};
            N = dims.size();
            num = std::accumulate(dims.begin(), dims.end(), (size_t) 1, std::multiplies<size_t>());
            blockSize = (N == 1 ? 128 : (N == 2 ? 16 : 6));
            pred_dim = N;
            stride = blockSize;
        }

        template<class Iter>
        size_t setDims(Iter begin, Iter end) {
            dims = std::vector<size_t>(begin, end);
            N = dims.size();
            num = std::accumulate(dims.begin(), dims.end(), (size_t) 1, std::multiplies<size_t>());
            return num;
        }

        void loadcfg(const std::string &cfgpath) {
            INIReader cfg(cfgpath);

            if (cfg.ParseError() != 0) {
                std::cout << "Can't load cfg file " << cfgpath << std::endl;
                exit(0);
            }

            auto cmprAlgoStr = cfg.Get("GlobalSettings", "CmprAlgo", "");
            if (cmprAlgoStr == ALGO_STR[ALGO_LORENZO_REG]) {
                cmprAlgo = ALGO_LORENZO_REG;
            } else if (cmprAlgoStr == ALGO_STR[ALGO_INTERP_LORENZO]) {
                cmprAlgo = ALGO_INTERP_LORENZO;
            } else if (cmprAlgoStr == ALGO_STR[ALGO_INTERP]) {
                cmprAlgo = ALGO_INTERP;
            }
            auto ebModeStr = cfg.Get("GlobalSettings", "ErrorBoundMode", "");
            if (ebModeStr == EB_STR[EB_ABS]) {
                errorBoundMode = EB_ABS;
            } else if (ebModeStr == EB_STR[EB_REL]) {
                errorBoundMode = EB_REL;
            } else if (ebModeStr == EB_STR[EB_PSNR]) {
                errorBoundMode = EB_PSNR;
            } else if (ebModeStr == EB_STR[EB_L2NORM]) {
                errorBoundMode = EB_L2NORM;
            } else if (ebModeStr == EB_STR[EB_ABS_AND_REL]) {
                errorBoundMode = EB_ABS_AND_REL;
            } else if (ebModeStr == EB_STR[EB_ABS_OR_REL]) {
                errorBoundMode = EB_ABS_OR_REL;
            }
            absErrorBound = cfg.GetReal("GlobalSettings", "AbsErrorBound", absErrorBound);
            relErrorBound = cfg.GetReal("GlobalSettings", "RelErrorBound", relErrorBound);
            psnrErrorBound = cfg.GetReal("GlobalSettings", "PSNRErrorBound", psnrErrorBound);
            l2normErrorBound = cfg.GetReal("GlobalSettings", "L2NormErrorBound", l2normErrorBound);

            openmp = cfg.GetBoolean("GlobalSettings", "OpenMP", openmp);
            lorenzo = cfg.GetBoolean("AlgoSettings", "Lorenzo", lorenzo);
            lorenzo2 = cfg.GetBoolean("AlgoSettings", "Lorenzo2ndOrder", lorenzo2);
            regression = cfg.GetBoolean("AlgoSettings", "Regression", regression);
            regression2 = cfg.GetBoolean("AlgoSettings", "Regression2ndOrder", regression2);

            auto interpAlgoStr = cfg.Get("AlgoSettings", "InterpolationAlgo", "");
            if (interpAlgoStr == INTERP_ALGO_STR[INTERP_ALGO_LINEAR]) {
                interpAlgo = INTERP_ALGO_LINEAR;
            } else if (interpAlgoStr == INTERP_ALGO_STR[INTERP_ALGO_CUBIC]) {
                interpAlgo = INTERP_ALGO_CUBIC;
            }
            interpDirection = cfg.GetInteger("AlgoSettings", "InterpolationDirection", interpDirection);
            interpBlockSize = cfg.GetInteger("AlgoSettings", "InterpolationBlockSize", interpBlockSize);
            blockSize = cfg.GetInteger("AlgoSettings", "BlockSize", blockSize);
            quantbinCnt = cfg.GetInteger("AlgoSettings", "QuantizationBinTotal", quantbinCnt);
            // additional variable
            detection_block_size = cfg.GetInteger("ArtifactSettings", "DetectionBlockSize", detection_block_size);
            detection_threshold = cfg.GetReal("ArtifactSettings", "DetectionThreshold", detection_threshold);
            detection_eb_rate = cfg.GetReal("ArtifactSettings", "DetectionEBRate", detection_eb_rate);
            noise_rate = cfg.GetReal("ArtifactSettings", "NoiseRate", noise_rate);
            auto block_sift_mode_str = cfg.Get("ArtifactSettings", "block_sift_mode", "");
            if(block_sift_mode_str == BLOCK_SIFT_MODE_STR[VARIANCE]){
                block_sift_mode = VARIANCE;
            }
            else if(block_sift_mode_str == BLOCK_SIFT_MODE_STR[RANGE])
            {
                block_sift_mode = RANGE;
            }
            else if(block_sift_mode_str == BLOCK_SIFT_MODE_STR[BLOCK_MAX])
            {
                block_sift_mode = BLOCK_MAX;
            }
            else if (block_sift_mode_str == BLOCK_SIFT_MODE_STR[ISOVALUE])
            {
                block_sift_mode = ISOVALUE;
            }
            block_flush_on = cfg.GetBoolean("ArtifactSettings", "block_flush_on", block_flush_on);
            block_sift_on = cfg.GetBoolean("ArtifactSettings", "block_sift_on", block_sift_on);
            block_iso_on = cfg.GetBoolean("ArtifactSettings", "block_iso_on", block_iso_on);
            block_isovalue =  cfg.GetReal("ArtifactSettings", "block_isovalue", block_isovalue);
            
            
            // pattern_check_on = cfg.GetBoolean("ArtifactSettings", "pattern_check_on", pattern_check_on);
            // block_diff_on = cfg.GetBoolean("ArtifactSettings", "block_diff_on", block_diff_on);
            // diff_thresh = cfg.GetReal("ArtifactSettings", "block_diff_thresh", diff_thresh);
            // use_stochastic_quantize = cfg.GetBoolean("ArtifactSettings", "use_stochastic_quantize", use_stochastic_quantize);
            // use_stochastic_decompress = cfg.GetBoolean("ArtifactSettings", "use_stochastic_decompress", use_stochastic_decompress);
            // use_stochastic_predict = cfg.GetBoolean("ArtifactSettings", "use_stochastic_predict", use_stochastic_predict);
            // use_stochastic_eb = cfg.GetBoolean("ArtifactSettings", "use_stochastic_eb", use_stochastic_eb);
            // use_random_predictor = cfg.GetBoolean("ArtifactSettings", "use_random_predictor", use_random_predictor);
            // normal_std = cfg.GetReal("ArtifactSettings", "normal_std", normal_std);
            // uniform_lower = cfg.GetReal("ArtifactSettings", "uniform_lower", uniform_lower);
            // normal_mean = cfg.GetReal("ArtifactSettings", "normal_mean", normal_mean);
            // bernoulli_p = cfg.GetReal("ArtifactSettings", "bernoulli_p", bernoulli_p);
            // random_seed = cfg.GetInteger("ArtifactSettings", "random_seed", random_seed);

            // block_noise_rng_threshold = cfg.GetReal("ArtifactSettings", "block_noise_rng_threshold", block_noise_rng_threshold);
            // post_block_noise_rng_on = cfg.GetBoolean("ArtifactSettings", "post_block_noise_rng_on", post_block_noise_rng_on);

            // region_error_on = cfg.GetBoolean("ArtifactSettings", "region_error_on", region_error_on);
            // region_error_control_start_level = cfg.GetInteger("ArtifactSettings", "region_error_control_start_level", region_error_control_start_level);
            // regional_error_block_size = cfg.GetInteger("ArtifactSettings", "regional_error_block_size", regional_error_block_size);
            // region_error_control_threshold = cfg.GetReal("ArtifactSettings", "region_error_control_threshold", region_error_control_threshold);
            auto region_error_control_mode_str = cfg.Get("ArtifactSettings", "region_error_control_mode", "");
            if(region_error_control_mode_str == REGION_ERROR_CONTROL_MODE_STR[REDUCE_EB]){
                region_error_control_mode = REDUCE_EB;
            }
            else if(region_error_control_mode_str == REGION_ERROR_CONTROL_MODE_STR[COMPENSATE_EB])
            {
                region_error_control_mode = COMPENSATE_EB;
            }
            region_error_control_eb_reduction = cfg.GetReal("ArtifactSettings", "region_error_control_eb_reduction", region_error_control_eb_reduction);
            region_error_control_eb_compensation = cfg.GetReal("ArtifactSettings", "region_error_control_eb_compensation", region_error_control_eb_compensation);
            quantization_prediction_on = cfg.GetBoolean("ArtifactSettings", "quantization_prediction_on", quantization_prediction_on);
            quantization_prediction_start_level = cfg.GetInteger("ArtifactSettings", "quantization_prediction_start_level", quantization_prediction_start_level);
            error_smoothing = cfg.GetBoolean("ArtifactSettings", "error_smoothing", error_smoothing);
            compress_error = cfg.GetBoolean("ArtifactSettings", "compress_error", compress_error);
            
            post_process_on = cfg.GetBoolean("ArtifactSettings", "post_process_on", post_process_on);
            corss_block_cubic = cfg.GetBoolean("ArtifactSettings", "corss_block_cubic", corss_block_cubic);
            use_natural_cubic = cfg.GetBoolean("ArtifactSettings", "use_natural_cubic", use_natural_cubic);
            eb_factor = cfg.GetReal("ArtifactSettings", "eb_factor", eb_factor);
        }


        void save(unsigned char *&c) {
            write(N, c);
            write(dims.data(), dims.size(), c);
            write(num, c);
            write(cmprAlgo, c);
            write(errorBoundMode, c);
            write(absErrorBound, c);
            write(relErrorBound, c);
            write(lorenzo, c);
            write(lorenzo2, c);
            write(regression, c);
            write(regression2, c);
            write(interpAlgo, c);
            write(interpDirection, c);
            write(interpBlockSize, c);
            write(lossless, c);
            write(encoder, c);
            write(quantbinCnt, c);
            write(blockSize, c);
            write(stride, c);
            write(pred_dim, c);
            write(openmp, c);
            // add additional variable
            // write(detection_block_size, c);
            // write(detection_threshold, c);
            // write(detection_eb_rate, c);
            // write(noise_rate, c);
            // write(block_sift_mode, c);
            // write(block_sift_on, c);
            // write(block_flush_on, c);
            // write(use_stochastic_decompress, c);
            // write(use_stochastic_quantize, c);
            // write(use_stochastic_predict, c);
            // write(use_stochastic_eb, c);
            // write(random_seed, c);
            // write(post_block_noise_rng_on, c);
            // write(block_noise_rng_threshold, c);


            

        };

        void load(const unsigned char *&c) {
            read(N, c);
            dims.resize(N);
            read(dims.data(), N, c);
            read(num, c);
            read(cmprAlgo, c);
            read(errorBoundMode, c);
            read(absErrorBound, c);
            read(relErrorBound, c);
            read(lorenzo, c);
            read(lorenzo2, c);
            read(regression, c);
            read(regression2, c);
            read(interpAlgo, c);
            read(interpDirection, c);
            read(interpBlockSize, c);
            read(lossless, c);
            read(encoder, c);
            read(quantbinCnt, c);
            read(blockSize, c);
            read(stride, c);
            read(pred_dim, c);
            read(openmp, c);
            // add additional variable
            // read(detection_block_size, c);
            // read(detection_threshold, c);
            // read(detection_eb_rate, c);
            // read(noise_rate, c);
            // read(block_sift_mode, c);
            // read(block_sift_on, c);
            // read(block_flush_on, c);
            // read(use_stochastic_decompress, c);
            // read(use_stochastic_quantize, c);
            // read(use_stochastic_predict, c);
            // read(use_stochastic_eb, c);
            // read(random_seed, c);
            // read(post_block_noise_rng_on, c);
            // read(block_noise_rng_threshold, c);
        }

        void print() {
            printf("CmprAlgo = %s\n", enum2Str((ALGO) cmprAlgo));
        }

        static size_t size_est() {
            return sizeof(size_t) * 5 + sizeof(double) * 4 + sizeof(bool) * 5 + sizeof(uint8_t) * 6 + sizeof(int) * 5 + 50; //50 is for redundancy
        }

        char N;
        std::vector<size_t> dims;
        size_t num;
        uint8_t cmprAlgo = ALGO_INTERP_LORENZO;
        uint8_t errorBoundMode = EB_ABS;
        double absErrorBound;
        double relErrorBound;
        double psnrErrorBound;
        double l2normErrorBound;
        bool lorenzo = true;
        bool lorenzo2 = false;
        bool regression = true;
        bool regression2 = false;
        bool openmp = false;
        uint8_t lossless = 1; // 0-> skip lossless(use lossless_bypass); 1-> zstd
        uint8_t encoder = 1;// 0-> skip encoder; 1->HuffmanEncoder; 2->ArithmeticEncoder
        uint8_t interpAlgo = INTERP_ALGO_CUBIC;
        uint8_t interpDirection = 0;
        int interpBlockSize = 32;
        int quantbinCnt = 65536;
        int blockSize;
        int stride; //not used now
        int pred_dim; // not used now
        // for artifact mitagation
        int detection_block_size = 4;
        double detection_threshold = 0.9;
        double detection_eb_rate = 1.0 / sqrt(4.4159889);
        double noise_rate = 0;
        uint8_t block_sift_mode = RANGE;
        bool block_flush_on = 0;
        bool block_sift_on =0;
        bool block_iso_on = 0;
        bool pattern_check_on = 0;
        double block_isovalue=0;
        double pattern_eb_rate = 1e-23;

        bool block_diff_on = 0;
        double diff_thresh = 0.1;

        bool use_stochastic_quantize=0;
        bool use_stochastic_predict=0;
        bool use_stochastic_decompress=0;
        bool use_stochastic_eb = 0;
        bool use_random_predictor = 0;
        float normal_mean = 0.0;
        float normal_std = 1.0;
        float uniform_lower = 0.0;
        float bernoulli_p = 0.5;
        int random_seed = 2333;

        double block_noise_rng_threshold = 0.0;
        bool post_block_noise_rng_on = 0;

        // for region error control
        bool region_error_on = 1;
        int region_error_control_start_level = 3;
        int regional_error_block_size = 4; 
        double region_error_control_threshold = 0.9;
        uint8_t region_error_control_mode = REDUCE_EB;
        double region_error_control_eb_reduction = 0.5;
        double region_error_control_eb_compensation = 0.5;

        bool quantization_prediction_on = 1;
        int quantization_prediction_start_level = 3; 

        bool error_smoothing = 1;

        bool compress_error = false;

        bool post_process_on = false;
        bool corss_block_cubic = false; 
        bool use_natural_cubic = false;

        double eb_factor = 2.0;

        std::vector<int> level_dirctions; 
        std::vector<double> level_eb_factors; 


        // for passing data through the workflow that has a config object 
        struct PASS_DATA{
        // std::unique_ptr<void> processed_data_prt;
        std::shared_ptr<void> processed_data_prt;
        const void* original_data_prt;
        std::shared_ptr<std::vector<int>> aux_quant_inds_ptr;
        } PASS_DATA;
    };


}

#endif //SZ_CONFIG_HPP