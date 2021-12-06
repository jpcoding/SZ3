#include "frontend/SZ3QoIFrontend.hpp"
#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "predictor/RegressionPredictor.hpp"
#include "predictor/ComposedPredictor.hpp"
#include "quantizer/QoIIntegerQuantizer.hpp"
#include "qoi/XSquare.hpp"
#include "utils/FileUtil.hpp"
#include "utils/Config.hpp"
#include "def.hpp"
#include <cstdio>
#include <iostream>
#include <memory>

#include "sz_backend.hpp"


template<typename T, uint N>
float SZ_compress_build_frontend(std::unique_ptr<T[]> const &data, const SZ::Config<T, N> &conf) {
    auto quantizer = SZ::VariableEBLinearQuantizer<T, T>(conf.quant_state_num / 2);
    // auto quantizer = SZ::LinearQuantizer<T>(conf.eb, conf.quant_state_num / 2);
    auto quantizer_eb = SZ::EBLogQuantizer<T>(conf.qoi_eb_base, conf.qoi_log_base, conf.qoi_quant_state_num / 2);
    auto qoi = SZ::QoI_X_Square<T>(conf.qoi_eb, conf.num, conf.eb);
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N>>> predictors;

    int use_single_predictor =
            (conf.enable_lorenzo + conf.enable_2ndlorenzo + conf.enable_regression) == 1;
    if (conf.enable_lorenzo) {
        if (use_single_predictor) {
            return SZ_compress_build_backend<T>(data, conf,
                                                make_sz3_qoi_frontend(conf, SZ::LorenzoPredictor<T, N, 1>(conf.eb),
                                                                  quantizer, quantizer_eb, qoi));
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.eb));
        }
    }
    if (conf.enable_2ndlorenzo) {
        if (use_single_predictor) {
            return SZ_compress_build_backend<T>(data, conf,
                                                make_sz3_qoi_frontend(conf, SZ::LorenzoPredictor<T, N, 2>(conf.eb),
                                                                  quantizer, quantizer_eb, qoi));
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 2>>(conf.eb));
        }
    }
    if (conf.enable_regression) {
        if (use_single_predictor) {
            return SZ_compress_build_backend<T>(data, conf,
                                                make_sz3_qoi_frontend(conf, SZ::RegressionPredictor<T, N>(conf.block_size,
                                                                                                      conf.eb),
                                                                  quantizer, quantizer_eb, qoi));
        } else {
            predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N>>(conf.block_size, conf.eb));
        }
    }

    return SZ_compress_build_backend<T>(data, conf,
                                        make_sz3_qoi_frontend(conf, SZ::ComposedPredictor<T, N>(predictors), quantizer,
                                                                  quantizer_eb, qoi));
}


template<class T, uint N>
float SZ_compress_parse_args(int argc, char **argv, int argp, std::unique_ptr<T[]> &data, float qoi_eb, float global_eb,
                             std::array<size_t, N> dims) {
    SZ::Config<float, N> conf(global_eb, dims);

    float max = data[0];
    float min = data[0];
    for (int i = 1; i < conf.num; i++) {
        if (max < data[i]) max = data[i];
        if (min > data[i]) min = data[i];
    }
    float eb = global_eb * (max - min);

    conf.encoder_op = 4;
    // qoi = x^2
    auto max_val = std::max(fabs(max), fabs(min));
    conf.qoi_eb = max_val * max_val * qoi_eb;
    conf.eb = eb;

    if (argp < argc) {
        int block_size = atoi(argv[argp++]);
        conf.block_size = block_size;
        conf.stride = block_size;
    }

    int lorenzo_op = 1;
    if (argp < argc) {
        lorenzo_op = atoi(argv[argp++]);
        conf.enable_lorenzo = lorenzo_op == 1 || lorenzo_op == 3;
        conf.enable_2ndlorenzo = lorenzo_op == 2 || lorenzo_op == 3;
    }

    int regression_op = 1;
    if (argp < argc) {
        regression_op = atoi(argv[argp++]);
        conf.enable_regression = regression_op == 1;
    }

    if (argp < argc) {
        conf.encoder_op = atoi(argv[argp++]);
        if (conf.encoder_op == 0) {
            conf.quant_state_num = 250;
        } else if (conf.encoder_op == 2) {
            conf.quant_state_num = 1024;
        }
    }

    if (argp < argc) {
        conf.lossless_op = atoi(argv[argp++]);
    }

    if (argp < argc) {
        conf.quant_state_num = atoi(argv[argp++]);
    }

    auto ratio = SZ_compress_build_frontend(data, conf);
    return ratio;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "SZ v" << SZ_versionString() << std::endl;
        std::cout << "usage: " << argv[0] <<
                  " data_file -num_dim dim0 .. dimn relative_eb [blocksize lorenzo_op regression_op encoder_op lossless_op quant_state_num]"
                  << std::endl;
        std::cout << "example: " << argv[0] <<
                  " qmcpack.dat -3 33120 69 69 1e-3 [6 1 1 1 1 32768]" << std::endl;
        return 0;
    }

    size_t num = 0;
    auto data = SZ::readfile<float>(argv[1], num);
    src_file_name = argv[1];
    std::cout << "Read " << num << " elements\n";

    int dim = atoi(argv[2] + 1);
    assert(1 <= dim && dim <= 4);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }

    float qoi_eb = atof(argv[argp++]);
    float global_eb = atof(argv[argp++]);

    if (dim == 1) {
        SZ_compress_parse_args<float, 1>(argc, argv, argp, data, qoi_eb, global_eb, std::array<size_t, 1>{dims[0]});
    } else if (dim == 2) {
        SZ_compress_parse_args<float, 2>(argc, argv, argp, data, qoi_eb, global_eb, std::array<size_t, 2>{dims[0], dims[1]});
    } else if (dim == 3) {
        SZ_compress_parse_args<float, 3>(argc, argv, argp, data, qoi_eb, global_eb,
                                         std::array<size_t, 3>{dims[0], dims[1], dims[2]});
    } else if (dim == 4) {
        SZ_compress_parse_args<float, 4>(argc, argv, argp, data, qoi_eb, global_eb,
                                         std::array<size_t, 4>{dims[0], dims[1], dims[2], dims[3]});
    }

    return 0;
}