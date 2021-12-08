#ifndef SZ3_FRONTEND
#define SZ3_FRONTEND
/**
 * This module contains SZ3 Predictor and Quantizer for QOI preservation
 */

#include "Frontend.hpp"
#include "def.hpp"
#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "quantizer/Quantizer.hpp"
#include "utils/Iterator.hpp"
#include "utils/Config.hpp"
#include "utils/MemoryUtil.hpp"

namespace SZ {


    template<class T, uint N, class Predictor, class Quantizer, class Quantizer_EB, class QoI>
    class SZ3QoIFrontend : public concepts::FrontendInterface<T, N> {
    public:
        SZ3QoIFrontend(const Config<T, N> &conf, Predictor predictor, Quantizer quantizer, Quantizer_EB quantizer_eb, QoI qoi) :
                fallback_predictor(LorenzoPredictor<T, N, 1>(conf.eb)),
                predictor(predictor),
                quantizer(quantizer),
                block_size(conf.block_size),
                stride(conf.stride),
                global_dimensions(conf.dims),
                num_elements(conf.num),
                quantizer_eb(quantizer_eb),
                qoi(qoi) {
        }

        std::vector<int> compress(T *data) {
            // 0 ~ num_elements - 1: quant_inds_eb
            // num_elements ~ 2*num_elements - 1: quant_inds_data
            std::vector<int> quant_inds(num_elements * 2);

            auto block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                    data, std::begin(global_dimensions), std::end(global_dimensions), stride, 0);

            auto element_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                    data, std::begin(global_dimensions), std::end(global_dimensions), 1, 0);

            predictor.precompress_data(block_range->begin());
            quantizer.precompress_data();
            size_t quant_count = 0;
            for (auto block = block_range->begin(); block != block_range->end(); ++block) {

                element_range->update_block_range(block, block_size);

                concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
                if (!predictor.precompress_block(element_range)) {
                    predictor_withfallback = &fallback_predictor;
                }
                predictor_withfallback->precompress_block_commit();

                for (auto element = element_range->begin(); element != element_range->end(); ++element) {
                    auto ori_data = *element;
                    // interpret the error bound for current data based on qoi
                    auto eb = qoi.interpret_eb(ori_data);
                    quant_inds[quant_count] = quantizer_eb.quantize_and_overwrite(eb);
                    quant_inds[num_elements + quant_count] = quantizer.quantize_and_overwrite(
                            *element, predictor_withfallback->predict(element), eb);
                    quant_count ++;
                    // update cumulative tolerance if needed 
                    qoi.update_tolerance(ori_data, *element);
                }
                qoi.postcompress_block();
            }

            predictor.postcompress_data(block_range->begin());
            quantizer.postcompress_data();
            return quant_inds;
        }

        T *decompress(std::vector<int> &quant_inds) {

            int const *quant_inds_eb_pos = (int const *) quant_inds.data();
            int const *quant_inds_pos = quant_inds_eb_pos + num_elements;
            std::array<size_t, N> intra_block_dims;
            auto dec_data = new T[num_elements];
            auto block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                    dec_data, std::begin(global_dimensions), std::end(global_dimensions), stride, 0);

            auto element_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                    dec_data, std::begin(global_dimensions), std::end(global_dimensions), 1, 0);

            predictor.predecompress_data(block_range->begin());
            quantizer.predecompress_data();

            for (auto block = block_range->begin(); block != block_range->end(); ++block) {

                element_range->update_block_range(block, block_size);

                concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
                if (!predictor.predecompress_block(element_range)) {
                    predictor_withfallback = &fallback_predictor;
                }
                for (auto element = element_range->begin(); element != element_range->end(); ++element) {
                    auto eb = quantizer_eb.recover(*(quant_inds_eb_pos++));
                    *element = quantizer.recover(predictor_withfallback->predict(element), *(quant_inds_pos++), eb);
                }
            }
            predictor.postdecompress_data(block_range->begin());
            quantizer.postdecompress_data();
            return dec_data;
        }

        void save(uchar *&c) {
            write(global_dimensions.data(), N, c);
            write(block_size, c);

            predictor.save(c);
            quantizer_eb.save(c);
            quantizer.save(c);
        }

        void load(const uchar *&c, size_t &remaining_length) {
            read(global_dimensions.data(), N, c, remaining_length);
            num_elements = 1;
            for (const auto &d: global_dimensions) {
                num_elements *= d;
                std::cout << d << " ";
            }
            std::cout << std::endl;
            read(block_size, c, remaining_length);
            stride = block_size;
            predictor.load(c, remaining_length);
            quantizer_eb.load(c, remaining_length);
            quantizer.load(c, remaining_length);
        }

        void print() {}

        void clear() {
            predictor.clear();
            fallback_predictor.clear();
            quantizer.clear();
            quantizer_eb.clear();
        }

        int get_radius() const { 
            return quantizer.get_radius(); 
        }

        size_t get_num_elements() const { 
            return num_elements; 
        }

    private:
        Predictor predictor;
        LorenzoPredictor<T, N, 1> fallback_predictor;
        Quantizer_EB quantizer_eb;
        Quantizer quantizer;
        QoI qoi;
        uint block_size;
        uint stride;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;
    };

    template<class T, uint N, class Predictor, class Quantizer, class Quantizer_EB, class QoI>
    SZ3QoIFrontend<T, N, Predictor, Quantizer, Quantizer_EB, QoI>
    make_sz3_qoi_frontend(const Config<T, N> &conf, Predictor predictor, Quantizer quantizer, Quantizer_EB quantizer_eb, QoI qoi) {
        return SZ3QoIFrontend<T, N, Predictor, Quantizer, Quantizer_EB, QoI>(conf, predictor, quantizer, quantizer_eb, qoi);
    }
}

#endif
