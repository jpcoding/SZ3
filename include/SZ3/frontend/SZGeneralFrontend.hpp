#ifndef SZ3_FRONTEND
#define SZ3_FRONTEND
/**
 * This module is the implementation of general frontend in SZ3
 */

#include "Frontend.hpp"
#include "SZ3/def.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

#include "SZ3/postprocess/postprocess.hpp"

namespace SZ {


    template<class T, uint N, class Predictor, class Quantizer>
    class SZGeneralFrontend : public concepts::FrontendInterface<T, N> {
    public:

        SZGeneralFrontend(const Config &conf, Predictor predictor, Quantizer quantizer) :
                fallback_predictor(LorenzoPredictor<T, N, 1>(conf.absErrorBound)),
                predictor(predictor),
                quantizer(quantizer),
                block_size(conf.blockSize),
                num_elements(conf.num) {
            std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
        }

        ~SZGeneralFrontend() = default;

        std::vector<int> compress(T *data) {
            std::vector<int> quant_inds(num_elements);
            
            //ANALYSIS COMPILATION
            #ifdef SZ_ANALYSIS
            std::cout << "SZ_ANALYSIS: SZGeneralFrontend::compress" << std::endl;
            std::vector<T> my_pred(num_elements);
            std::vector<int> my_quant_inds(num_elements);
            #endif 

            auto block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                    data, std::begin(global_dimensions), std::end(global_dimensions), block_size, 0);

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
                    auto pred = predictor_withfallback->predict(element);
                    quant_inds[quant_count++] = quantizer.quantize_and_overwrite(
                            *element, pred);
                    #ifdef SZ_ANALYSIS
                    my_pred[element.get_offset()] = pred;
                    my_quant_inds[element.get_offset()] = quant_inds[quant_count-1];
                    #endif
                }
            }
            predictor.postcompress_data(block_range->begin());
            quantizer.postcompress_data();
            #ifdef SZ_ANALYSIS
            std::cout<<"SZ_ANALYSIS: SZGeneralFrontend::compress: writing files"<<std::endl;
            writefile("pred.dat", my_pred.data(), num_elements);
            writefile("quant.dat", my_quant_inds.data(), num_elements);
            writefile("decompressed.dat", data, num_elements);
            #endif
        if(N==3)
        {
        std::cout << "3D post process" << std::endl;
        compensation_3d_(
        data, my_quant_inds.data(), global_dimensions.data(), 0,
        quantizer.get_eb()*0.5);
        writefile("postdecompressed.dat", data, num_elements);
        }
        


            return quant_inds;
        }

        T *decompress(std::vector<int> &quant_inds, T *dec_data) {

            int const *quant_inds_pos = (int const *) quant_inds.data();
            std::array<size_t, N> intra_block_dims;
//            auto dec_data = new T[num_elements];
            auto block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                    dec_data, std::begin(global_dimensions), std::end(global_dimensions), block_size, 0);

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
                    *element = quantizer.recover(predictor_withfallback->predict(element), *(quant_inds_pos++));
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
            quantizer.save(c);
        }

        void load(const uchar *&c, size_t &remaining_length) {
            read(global_dimensions.data(), N, c, remaining_length);
            num_elements = 1;
            for (const auto &d: global_dimensions) {
                num_elements *= d;
            }
            read(block_size, c, remaining_length);
            predictor.load(c, remaining_length);
            quantizer.load(c, remaining_length);
        }

        size_t size_est() {
            return quantizer.size_est();
        }

        void print() {
//            predictor.print();
//            quantizer.print();
        }

        void clear() {
            predictor.clear();
            fallback_predictor.clear();
            quantizer.clear();
        }

        int get_radius() const { return quantizer.get_radius(); }

        size_t get_num_elements() const { return num_elements; };

    private:
        Predictor predictor;
        LorenzoPredictor<T, N, 1> fallback_predictor;
        Quantizer quantizer;
        uint block_size;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;

    };

    template<class T, uint N, class Predictor, class Quantizer>
    SZGeneralFrontend<T, N, Predictor, Quantizer>
    make_sz_general_frontend(const Config &conf, Predictor predictor, Quantizer quantizer) {
        return SZGeneralFrontend<T, N, Predictor, Quantizer>(conf, predictor, quantizer);
    }
}

#endif