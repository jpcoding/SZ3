#include "SZ3/compressor/SZInterpolationCompressor.hpp"
#include "SZ3/compressor/deprecated/SZBlockInterpolationCompressor.hpp"
#include "SZ3/compressor/SZGeneralCompressor.hpp"
#include "SZ3/frontend/SZMetaFrontend.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/predictor/ComposedPredictor.hpp"
#include "SZ3/predictor/SimplePredictor.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Verification.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/Config.hpp"
#include <cstdio>
#include <iostream>
#include <cmath>
#include <memory>
#include <type_traits>
#include <sstream>


template<class T, uint N, class ... Dims>
SZ::CompressStats
interp_compress_decompress(char *path, T *data, size_t num, double eb, int interp_op,
                           int direction_op, int block_size, int interp_block_size, Dims ... args) {
    SZ::CompressStats compressStats;
    std::string compressed_file_name;
    {

        std::cout << "****************** Interp Compression ****************" << std::endl;
        std::cout << "Interp Op          = " << interp_op << std::endl
                  << "Direction          = " << direction_op << std::endl
                  << "SZ block size      = " << block_size << std::endl
                  << "Interp block size  = " << interp_block_size << std::endl;

        SZ::Timer timer(true);
        size_t compressed_size = 0;
        std::unique_ptr<unsigned char[]> compressed;


        auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};
        auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer<T>(eb),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd(),
                dims,
                interp_block_size,
                interp_op,
                direction_op
        );
        compressed.reset(sz.compress(data, compressed_size));


        double compress_time = timer.stop();
        auto compression_ratio = num * sizeof(T) * 1.0 / compressed_size;
        compressStats.compress_time = compress_time;
        compressed_file_name = SZ::compressed_path(path, false);

        std::cout << "Compression time = " << compress_time << "s" << std::endl;
        std::cout << "Compression size = " << compressed_size << std::endl;
        std::cout << "Compression ratio = " << compression_ratio << std::endl;
        std::cout << "Compression file = " << compressed_file_name << std::endl;

        SZ::writefile(compressed_file_name.c_str(), compressed.get(), compressed_size);
    }
    {
        std::cout << "***************** Interp Decompression ****************" << std::endl;
        size_t compressed_size = 0;
        auto compressed = SZ::readfile<SZ::uchar>(compressed_file_name.c_str(), compressed_size);
//        remove(compressed_file_name.c_str());

        SZ::Timer timer(true);
        std::unique_ptr<T[]> dec_data;

        auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};
        auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer<T>(eb),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd(),
                dims,
                interp_block_size,
                interp_op,
                direction_op
        );
        dec_data.reset(sz.decompress(compressed.get(), compressed_size));


        compressStats.decompress_time = timer.stop();
        std::string decompressed_file_name = SZ::decompressed_path(path, false);

        std::cout << "Decompression time = " << compressStats.decompress_time << "s" << std::endl;
        std::cout << "Decompression file = " << decompressed_file_name << std::endl;
        SZ::writefile(decompressed_file_name.c_str(), dec_data.get(), num);

        size_t num1 = 0;
        auto ori_data = SZ::readfile<T>(path, num1);
        assert(num1 == num);
        double psnr, nrmse;
        SZ::verify<T>(ori_data.get(), dec_data.get(), num, psnr, nrmse);

        auto compression_ratio = num * sizeof(T) * 1.0 / compressed_size;
        printf("PSNR = %f, NRMSE = %.10G, Compression Ratio = %.2f\n", psnr, nrmse, compression_ratio);

        compressStats.psnr = psnr;
        compressStats.nrmse = nrmse;
        compressStats.ratio = compression_ratio;
        return compressStats;

    }

}

template<class T, uint N>
double interp_compress_block_version(T *data, std::array<size_t, N> dims, size_t num, double eb, int interp_level,
                                     int interp_op, int direction_op,
                                     int block_size, int interp_block_size) {

    std::cout << "****************** Interp Compression ****************" << std::endl;
    std::cout << "Interp Op          = " << interp_op << std::endl
              << "Direction          = " << direction_op << std::endl
              << "SZ block size      = " << block_size << std::endl
              << "Interp block size  = " << interp_block_size << std::endl;
//              << "SZ_mode            = " << sz_op << std::endl;

    SZ::Timer timer(true);

    std::vector<T> data1(data, data + num);
    size_t compressed_size = 0;

    SZ::Config conf(eb, dims);
    conf.block_size = block_size;
    conf.stride = conf.block_size;
    auto sz = SZ::make_sz_fast_block_interpolation_compressor(
            conf,
            SZ::SimplePredictor<T, N>(eb),
            SZ::LinearQuantizer<T>(eb),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd(),
            interp_op,
            direction_op,
            interp_level
    );

    sz.compress(data1.data(), compressed_size);

    double compression_time = timer.stop();

    auto compression_ratio = num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compressed size = " << compressed_size << std::endl;
    std::cout << "Compression ratio = " << compression_ratio << std::endl;
    std::cout << "Interp compression time = " << compression_time
              << " Ratio = " << compression_ratio
              << " Params = " << interp_level
              << " " << interp_op
              << " " << direction_op
              << " " << block_size
              << " " << interp_block_size
              //              << " " << sz_op
              << std::endl;


    std::cout << "****************** Interp end ****************" << std::endl;
    return compression_ratio;
}

template<typename T, uint N>
SZ::CompressStats
lorenzo_compress_decompress_3d(char *path, T *data, size_t num_elements, SZ::Config conf, bool decompress) {
    SZ::Timer timer(true);
    SZ::CompressStats compressStats;
    std::cout << "***************** Lorenzo Compression ****************" << std::endl;

    auto quantizer = SZ::LinearQuantizer<T>(conf.absErrorBound, conf.quant_state_num / 2);
    auto sz = make_sz_general_compressor(conf, make_sz_meta_frontend(conf, quantizer), SZ::HuffmanEncoder<int>(),
                                         SZ::Lossless_zstd());

    size_t compressed_size = 0;
    SZ::uchar *compressed = sz.compress(data, compressed_size);

    compressStats.compress_time = timer.stop();
    compressStats.ori_bytes = num_elements * sizeof(T);
    compressStats.compress_bytes = compressed_size;
    compressStats.ratio = compressStats.ori_bytes * 1.0 / compressStats.compress_bytes;

    std::string compressed_file_name = SZ::compressed_path(path, false);

    std::cout << "Compression time = " << compressStats.compress_time << "s" << std::endl;
    std::cout << "Compression size = " << compressed_size << std::endl;
    std::cout << "Compression ratio = " << compressStats.ratio << std::endl;


    delete[]compressed;
    if (decompress) {
        std::cout << "Compression file = " << compressed_file_name << std::endl;
        SZ::writefile(compressed_file_name.c_str(), compressed, compressed_size);

        std::cout << "***************** Lorenzo Decompression ****************" << std::endl;
        auto compressed = SZ::readfile<SZ::uchar>(compressed_file_name.c_str(), compressed_size);
//        remove(compressed_file_name.c_str());

        timer.start();
        quantizer = SZ::LinearQuantizer<T>(conf.absErrorBound, conf.quant_state_num / 2);
        sz = make_sz_general_compressor(make_sz_meta_frontend(conf, quantizer), SZ::HuffmanEncoder<int>(),
                                        SZ::Lossless_zstd());

        T *dec_data = sz->decompress(compressed.get(), compressed_size);

        compressStats.decompress_time = timer.stop();
        std::string decompressed_file_name = SZ::decompressed_path(path, false);

        std::cout << "Decompression time = " << compressStats.decompress_time << "s" << std::endl;
        std::cout << "Decompression file = " << decompressed_file_name << std::endl;
        SZ::writefile(decompressed_file_name.c_str(), dec_data, num_elements);

        size_t num1 = 0;
        auto ori_data = SZ::readfile<T>(path, num1);
        assert(num1 == num_elements);
        SZ::verify<T>(ori_data.get(), dec_data, num_elements, compressStats.psnr, compressStats.nrmse);

        printf("PSNR = %f, NRMSE = %.10G, Compression Ratio = %.2f\n",
               compressStats.psnr, compressStats.nrmse, compressStats.ratio);

    } else {
        std::cout << "***************** Lorenzo End ****************" << std::endl;
    }

    return compressStats;
}


template<class T, uint N, class ... Dims>
void interp_lorenzo_tuning(char *path, double reb, bool enable_lorenzo, Dims ... args) {
    std::cout << "================================ BEGIN TUNING ================================" << std::endl;

    size_t num = 0;
    auto data = SZ::readfile<T>(path, num);
    double eb = reb * SZ::data_range(data.get(), num);
    auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};

    SZ::Timer timer(true);

    size_t sampling_num, sampling_block;
    std::array<size_t, N> sample_dims;
    std::vector<T> sampling_data = SZ::sampling<T, N>(data.get(), dims, sampling_num, sample_dims, sampling_block);
//    printf("%lu %lu %lu %lu %lu\n", sampling_data.size(), sampling_num, sample_dims[0], sample_dims[1], sample_dims[2]);

    SZ::CompressStats lorenzo_stats;
    SZ::Config lorenzo_config(eb, sample_dims);
    if (enable_lorenzo) {
        if (N != 3) {
            printf("Lorenzo can only be enabled in 3D mode");
            exit(0);
        }
        lorenzo_config.enable_2ndlorenzo = true;
        lorenzo_config.enable_regression = false;
        lorenzo_config.block_size = 5;
        lorenzo_config.quant_state_num = 65536 * 2;

        lorenzo_stats = lorenzo_compress_decompress_3d(path, sampling_data.data(), sampling_num, lorenzo_config, false);
        printf("Lorenzo ratio = %.2f, compress_time:%.3f\n",
               lorenzo_stats.ratio, lorenzo_stats.compress_time);
    }
    double best_lorenzo_ratio = lorenzo_stats.ratio;
    double best_interp_ratio = 0, ratio;

    int interp_level = -1, interp_op, direction_op = 0, block_size = sampling_block, interp_block_size = sampling_block;
    for (int i = 0; i < 2; i++) {
        ratio = interp_compress_block_version<T, N>(sampling_data.data(), sample_dims, sampling_num, eb, interp_level,
                                                    i, direction_op,
                                                    block_size, interp_block_size);
        if (ratio > best_interp_ratio) {
            best_interp_ratio = ratio;
            interp_op = i;
        }
    }
    std::cout << "interp best interp_op = " << interp_op << " , best ratio = " << best_interp_ratio << std::endl;

    ratio = interp_compress_block_version<T, N>(sampling_data.data(), sample_dims, sampling_num, eb, interp_level,
                                                interp_op, 5,
                                                block_size, interp_block_size);
    if (ratio > best_interp_ratio * 1.02) {
        best_interp_ratio = ratio;
        direction_op = 5;
    }
    std::cout << "interp best direction_op = " << direction_op << " , best ratio = " << best_interp_ratio
              << std::endl;

    bool lorenzo = enable_lorenzo && lorenzo_stats.ratio > best_interp_ratio && lorenzo_stats.ratio < 80 &&
                   best_interp_ratio < 80;
    printf("\nLorenzo compression ratio = %.2f\n", lorenzo_stats.ratio);
    printf("Interp compression ratio = %.2f\n", best_interp_ratio);
    printf("choose %s\n", lorenzo ? "lorenzo" : "interp");

    if (lorenzo) {
        lorenzo_config.quant_state_num = SZ::optimize_quant_invl_3d<T>(data.get(), dims[0], dims[1], dims[2], eb);
        lorenzo_config.pred_dim = 2;
        lorenzo_stats = lorenzo_compress_decompress_3d(path, sampling_data.data(), sampling_num, lorenzo_config,
                                                       false);
        printf("Lorenzo, pred_dim=2, ratio = %.2f, compress_time:%.3f\n",
               lorenzo_stats.ratio, lorenzo_stats.compress_time);
        if (lorenzo_stats.ratio > best_lorenzo_ratio * 1.02) {
            best_lorenzo_ratio = lorenzo_stats.ratio;
        } else {
            lorenzo_config.pred_dim = 3;
        }

        int quant_num = lorenzo_config.quant_state_num;
        if (reb < 1.01e-6 && best_lorenzo_ratio > 5) {
            lorenzo_config.quant_state_num = 16384;
            lorenzo_stats = lorenzo_compress_decompress_3d(path, sampling_data.data(), sampling_num, lorenzo_config,
                                                           false);
            printf("Lorenzo, quant_bin=8192, ratio = %.2f, compress_time:%.3f\n",
                   lorenzo_stats.ratio, lorenzo_stats.compress_time);
            if (lorenzo_stats.ratio > best_lorenzo_ratio * 1.02) {
                best_lorenzo_ratio = lorenzo_stats.ratio;
            } else {
                lorenzo_config.quant_state_num = quant_num;
            }
        }

        double tuning_time = timer.stop();
        std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
        std::cout << "================================ END TUNING ================================" << std::endl
                  << std::endl;
        std::cout << "================================ BEGIN SZ-Lorenzo ================================" << std::endl;

        lorenzo_config.dims = dims;
        lorenzo_config.num = num;
        auto result = lorenzo_compress_decompress_3d(path, data.get(), num, lorenzo_config, true);
        std::cout << "\nTotal compress time (include tuning) = " << tuning_time + result.compress_time << std::endl;
        std::cout << "Total decompress time = " << result.decompress_time << std::endl;
        std::cout << "==================================== END SZ-Lorenzo ==================================="
                  << std::endl;

    } else {
        double tuning_time = timer.stop();
        std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
        std::cout << "====================================== END TUNING ======================================"
                  << std::endl << std::endl;
        std::cout << "==================================== BEGIN SZ-Interp ==================================="
                  << std::endl;

        block_size = 6;
        interp_block_size = 32;
        auto result = interp_compress_decompress<T, N>(path, data.get(), num, eb, interp_op, direction_op,
                                                       block_size, interp_block_size, args...);
        std::cout << "\nTotal compress time (include tuning) = " << tuning_time + result.compress_time << std::endl;
        std::cout << "Total decompress time = " << result.decompress_time << std::endl;
        std::cout << "==================================== END SZ-Interp ==================================="
                  << std::endl;
    }
}