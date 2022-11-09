#ifndef _SZ_INTERPOLATION_COMPRESSSOR_HPP
#define _SZ_INTERPOLATION_COMPRESSSOR_HPP

#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"
#include <cstring>
#include <cmath>
#define BIT_PER_BYTE 8

namespace SZ {
    template<class T, uint N, class Quantizer, class Encoder, class Lossless>
    class SZInterpolationCompressor {
    public:


        SZInterpolationCompressor(Quantizer quantizer, Encoder encoder, Lossless lossless) :
                quantizer(quantizer), encoder(encoder), lossless(lossless) {

            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");
        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, size_t num) {
            T *dec_data = new T[num];
            return decompress(cmpData, cmpSize, dec_data);
        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, T *decData) {
            size_t remaining_length = cmpSize;
            uchar *buffer = lossless.decompress(cmpData, remaining_length);
            uchar const *buffer_pos = buffer;

            read(global_dimensions.data(), N, buffer_pos, remaining_length);
            read(blocksize, buffer_pos, remaining_length);
            read(interpolator_id, buffer_pos, remaining_length);
            read(direction_sequence_id, buffer_pos, remaining_length);
            init();

            if (N==3)
            {
                std::vector<unsigned char> sift_index_bytes;
                my_sift_index.reserve(num_elements);
                size_t num_sift_index_bytes = num_elements/sizeof(unsigned char)/8+(num_elements%(sizeof(unsigned char)*BIT_PER_BYTE)!=0);
                sift_index_bytes.reserve(num_sift_index_bytes);
                for (size_t i =0; i< num_sift_index_bytes; i++)
                {
                    sift_index_bytes.push_back(0);
                }
                read(sift_index_bytes.data(), num_sift_index_bytes, buffer_pos,remaining_length);
                size_t byte_index=0;
                size_t data_index;
                for(size_t i =0; i<num_elements; i+=sizeof(unsigned char)*BIT_PER_BYTE)
                {
                    for(int j =0; j<sizeof(unsigned char)*BIT_PER_BYTE; j++)
                    {
                         data_index = i+j;
                         byte_index = data_index/8/sizeof(unsigned char);
                        if(data_index==num_elements)
                        {
                            break;
                        }
                        unsigned char tmp = sift_index_bytes[byte_index]<<(j);
                        tmp = tmp>>(sizeof(unsigned char)*BIT_PER_BYTE-1);
                        my_sift_index[data_index] = (tmp==1);
                    }
                }
                writefile("sift_index_decompressed.dat", sift_index_bytes.data(), num_sift_index_bytes);
            }

            quantizer.load(buffer_pos, remaining_length);
            encoder.load(buffer_pos, remaining_length);
            quant_inds = encoder.decode(buffer_pos, num_elements);

            encoder.postprocess_decode();

            lossless.postdecompress_data(buffer);
            double eb = quantizer.get_eb();

            *decData = quantizer.recover(0, quant_inds[quant_index++]);


            if( interpolators[interpolator_id] == "linear")
            {
                reduction_factor = sqrt(19/8);
            }
            else 
            {
                reduction_factor = sqrt(4.4159889);
            }
            real_eb_ratio = pow(1/reduction_factor, interpolation_level-1);

            

            for (uint level = interpolation_level; level > 0 && level <= interpolation_level; level--) {
                // if (level >= 3) {
                //     quantizer.set_eb(eb * eb_ratio);
                // } else {
                //     quantizer.set_eb(eb);
                // }
                quantizer.set_eb(eb * real_eb_ratio);
                real_eb_ratio *= reduction_factor;
                current_level = level;

                size_t stride = 1U << (level - 1);
                auto inter_block_range = std::make_shared<
                        SZ::multi_dimensional_range<T, N>>(decData,
                                                           std::begin(global_dimensions), std::end(global_dimensions),
                                                           stride * blocksize, 0);
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
                    block_interpolation(decData, block.get_global_index(), end_idx, PB_recover,
                                        interpolators[interpolator_id], direction_sequence_id, stride);

                }
            }
            quantizer.postdecompress_data();
//            timer.stop("Interpolation Decompress");

            return decData;
        }

        // compress given the error bound
        uchar *compress(const Config &conf, T *data, size_t &compressed_size) {
            std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
            blocksize = conf.interpBlockSize;
            interpolator_id = conf.interpAlgo;
            direction_sequence_id = conf.interpDirection;
            var_percentage = conf.var_percentage;

            init();

            for (int i =0; i< global_dimensions.size(); i++)
                printf("dims %ld \n", global_dimensions[i]);

            std::vector<size_t> my_test_index;
            if (N==3)
            {
                size_t xdim = global_dimensions[0];
                size_t ydim = global_dimensions[1];
                size_t zdim = global_dimensions[2];
                size_t x_block =16;
                size_t y_block =16;
                size_t z_block =16;
                size_t i, j, k;
                int ii, jj, kk;
                size_t x , y , z;
                printf("dims 0 %ld \n", xdim);
                printf("dims 1 %ld \n", ydim);
                printf("dims 2 %ld \n", zdim);
                size_t zy = zdim*ydim;
                size_t current_index;
                // printf("calculating\n");
                size_t ioffset, joffset;
                // printf("zy = %ld\nzdim = %ld",zy,zdim);
                double block_avg=0;
                double block_var=0;
                double block_max=data[0];
                double block_min=data[0];
                std::vector<double> block_vars;
                std::vector<double> block_range;

                for ( i = 0; i<xdim; i+=x_block)
                {
                    // ioffset = i*zy;
                    for ( j = 0; j <ydim; j+=y_block )
                    {
                        // joffset = j*zdim;
                        for ( k=0; k<zdim; k+=z_block)
                        {
                            for(ii=0; ii<x_block; ii++)
                            {
                                x = (i+ii<xdim)?(i+ii):xdim-1; 
                                ioffset = x*zy;
                                for (jj=0; jj<y_block; jj++)
                                {
                                    y = (j+jj<ydim)?(j+jj):ydim-1;
                                    joffset = y*zdim;
                                    for(kk=0 ; kk<z_block; kk++)
                                    {
                                        z=(k+kk<zdim)?(k+kk):zdim-1;
                                        current_index = ioffset+joffset+z;
                                        block_avg+=data[current_index];
                                        if(data[current_index]>block_max)
                                            block_max=data[current_index];
                                        if(data[current_index]<block_min)
                                            block_min=data[current_index];
                                    }
                                }
                            }
                            block_range.push_back(block_max-block_min);
                            block_avg= block_avg/(x_block*y_block*z_block);

                            for(ii=0; ii<x_block; ii++)
                            {
                                x = (i+ii<xdim)?(i+ii):xdim-1; 
                                ioffset = x*zy;
                                for (jj=0; jj<y_block; jj++)
                                {
                                    y = (j+jj<ydim)?(j+jj):ydim-1;
                                    joffset = y*zdim;
                                    for(kk=0 ; kk<z_block; kk++)
                                    {
                                        z=(k+kk<zdim)?(k+kk):zdim-1;
                                        current_index = ioffset+joffset+z;
                                        block_var+=(data[current_index]-block_avg)*(data[current_index]-block_avg);
                                    }
                                }
                            }
                            block_var=block_var/(x_block*y_block*z_block);
                            block_vars.push_back(block_var);
                            
                        }
                    }
                }
                
                double percentile = var_percentage;
                //  percentile = percentile*block_vars.size();
                size_t threshold_index = size_t (percentile*block_vars.size());
                std::vector<double> block_vars_copy;
                for (i= 0; i<block_vars.size(); i++)
                {
                    block_vars_copy.push_back(block_vars[i]);
                }

                std::sort(block_vars_copy.begin(), block_vars_copy.end());
                double threshold = block_vars_copy[threshold_index];
                std::cout<< "block_vars.size() = "<<block_vars.size()<<std::endl;

                printf("variance %2.f %% threshold = %f \n",percentile*100, threshold);

                my_sift_index.reserve(num_elements);
                sift_index=0;
               size_t ptr_block_vars=0;
               for ( i = 0; i<xdim; i+=x_block)
                {
                    // ioffset = i*zy;
                    for ( j = 0; j <ydim; j+=y_block )
                    {
                        // joffset = j*zdim;
                        for ( k=0; k<zdim; k+=z_block)
                        {
                            block_var=block_vars[ptr_block_vars];
                            for(ii=0; ii<x_block; ii++)
                            {
                                x = (i+ii<xdim)?(i+ii):xdim-1; 
                                ioffset = x*zy;
                                for (jj=0; jj<y_block; jj++)
                                {
                                    y = (j+jj<ydim)?(j+jj):ydim-1;
                                    joffset = y*zdim;
                                    for(kk=0 ; kk<z_block; kk++)
                                    {
                                        z=(k+kk<zdim)?(k+kk):zdim-1;
                                        current_index = ioffset+joffset+z;
                                        if (block_var>threshold)
                                        {   my_sift_index[current_index]=1;
                                            sift_index+=1;
                                        }
                                        else{
                                            my_sift_index[current_index]=0;
                                        }
                                    }
                                }
                            }
                            ptr_block_vars++;
                        }
                    }
                }

                // printf("done calculating\n");
                writefile("block_vars.dat",block_vars.data(), block_vars.size());
                


            }
            printf("sift_index %ld = \n", sift_index);
            
            quant_inds.reserve(num_elements);
            // my_level.reserve(num_elements);
            // my_index.reserve(num_elements);
            size_t interp_compressed_size = 0;

            double eb = quantizer.get_eb();

            quant_inds.push_back(quantizer.quantize_and_overwrite(*data, 0));
            // my_level[0] = 0;

            Timer timer;
            timer.start();


            if( interpolators[interpolator_id] == "linear")
            {
                reduction_factor = sqrt(19/8);
            }
            else 
            {
                reduction_factor = sqrt(4.4159889);
            }
            real_eb_ratio = pow(1/reduction_factor, interpolation_level-1);

            for (uint level = interpolation_level; level > 0 && level <= interpolation_level; level--) {
                // if (level >= 3) {
                //     quantizer.set_eb(eb * eb_ratio);
                // } else {
                //     quantizer.set_eb(eb);
                // }
                current_level = level;
                quantizer.set_eb(eb * real_eb_ratio);
                real_eb_ratio *= reduction_factor;
                size_t stride = 1U << (level - 1);

                auto inter_block_range = std::make_shared<
                        SZ::multi_dimensional_range<T, N>>(data, std::begin(global_dimensions),
                                                           std::end(global_dimensions),
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

                    block_interpolation(data, block.get_global_index(), end_idx, PB_predict_overwrite,
                                        interpolators[interpolator_id], direction_sequence_id, stride);
                }
            }
            assert(quant_inds.size() == num_elements);
//            timer.stop("Prediction & Quantization");

//            writefile("pred.dat", preds.data(), num_elements);
//            writefile("quant.dat", quant_inds.data(), num_elements);
            encoder.preprocess_encode(quant_inds, 0);
            size_t bufferSize = 1.2 * (quantizer.size_est() + encoder.size_est() + sizeof(T) * quant_inds.size());

            uchar *buffer = new uchar[bufferSize];
            uchar *buffer_pos = buffer;

            write(global_dimensions.data(), N, buffer_pos);
            write(blocksize, buffer_pos);
            write(interpolator_id, buffer_pos);
            write(direction_sequence_id, buffer_pos);

            if (N==3)
            {
                std::vector<unsigned char> sift_index_bytes;
                size_t num_sift_index_bytes = num_elements/sizeof(unsigned char)/BIT_PER_BYTE+(num_elements%(sizeof(unsigned char)*BIT_PER_BYTE)!=0);
                std::cout<<"num_sift_index_bytes = "<< num_sift_index_bytes<<std::endl;
                sift_index_bytes.reserve(num_sift_index_bytes);
                for (size_t i =0; i< num_sift_index_bytes; i++)
                {
                    sift_index_bytes.push_back(0);
                }
                size_t byte_index=0;
                size_t data_index;
                for(size_t i =0; i<num_elements; i+=sizeof(unsigned char)*BIT_PER_BYTE)
                {
                    for(int j =0; j<sizeof(unsigned char)*BIT_PER_BYTE; j++)
                    {
                        data_index = i+j;
                        byte_index = data_index/sizeof(unsigned char)/8;
                        if(data_index==num_elements)
                        {
                            std::cout<<"complted convertion"<<std::endl;
                            break;
                        }
                        unsigned char tmp = (1==my_sift_index[data_index]);
                        tmp = tmp<<(sizeof(unsigned char)*BIT_PER_BYTE-1-j);
                        sift_index_bytes.at(byte_index) = sift_index_bytes.at(byte_index) | tmp;
                    }
                }
                write(sift_index_bytes.data(), num_sift_index_bytes, buffer_pos);
                writefile("sift_index.dat", sift_index_bytes.data(), num_sift_index_bytes);
                std::cout<<"byte_index = "<<byte_index<<std::endl;
                std::cout<<"data_index = "<<data_index<<std::endl;
                std::cout<<"sift_index_bytes.size() "<<sift_index_bytes.size()<<std::endl;
            }

            quantizer.save(buffer_pos);
            quantizer.postcompress_data();


            timer.start();
            encoder.save(buffer_pos);
            encoder.encode(quant_inds, buffer_pos);
            encoder.postprocess_encode();
//            timer.stop("Coding");
            assert(buffer_pos - buffer < bufferSize);

            timer.start();
            uchar *lossless_data = lossless.compress(buffer,
                                                     buffer_pos - buffer,
                                                     compressed_size);
            lossless.postcompress_data(buffer);
//            timer.stop("Lossless");

            compressed_size += interp_compressed_size;
            writefile("compressed.dat",data,num_elements);
            // writefile("unpred_index.dat",quantizer.get_unpred_idx().data(),quantizer.get_unpred_idx().size());
            // writefile("level.dat", my_level.data(), num_elements);
            // writefile("index.dat", my_index.data(), num_elements);
            // printf("# of unpredict %ld \n", quantizer.get_unpred_idx().size());
            return lossless_data;
        }

    private:

        enum PredictorBehavior {
            PB_predict_overwrite, PB_predict, PB_recover
        };

        void init() {
            assert(blocksize % 2 == 0 && "Interpolation block size should be even numbers");
            num_elements = 1;
            interpolation_level = -1;
            for (int i = 0; i < N; i++) {
                if (interpolation_level < ceil(log2(global_dimensions[i]))) {
                    interpolation_level = (uint) ceil(log2(global_dimensions[i]));
                }
                num_elements *= global_dimensions[i];
            }

            dimension_offsets[N - 1] = 1;
            for (int i = N - 2; i >= 0; i--) {
                dimension_offsets[i] = dimension_offsets[i + 1] * global_dimensions[i + 1];
            }

            dimension_sequences = std::vector<std::array<int, N>>();
            auto sequence = std::array<int, N>();
            for (int i = 0; i < N; i++) {
                sequence[i] = i;
            }
            do {
                dimension_sequences.push_back(sequence);
            } while (std::next_permutation(sequence.begin(), sequence.end()));
        }

        inline void quantize(size_t idx, T &d, T pred) {
            double default_eb = quantizer.get_eb();
            if(my_sift_index[idx] && current_level ==1 )
            {
                quantizer.set_eb(default_eb/reduction_factor/2);
                quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
                quantizer.set_eb(default_eb);
                
            }
            // else if( my_sift_index[idx] )
            // {
            //     quantizer.set_eb(default_eb/reduction_factor/5);
            //     quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
            //     quantizer.set_eb(default_eb);
            // }
            else{
                quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
            }
            // my_level[idx] = current_level;
            
        }

        inline void recover(size_t idx, T &d, T pred) {
            double default_eb = quantizer.get_eb(); 
            if(my_sift_index[idx] && current_level ==1 )
            {
                quantizer.set_eb(default_eb/reduction_factor/2);
                d = quantizer.recover(pred, quant_inds[quant_index++]);
                quantizer.set_eb(default_eb);
            }
            else
            {
                d = quantizer.recover(pred, quant_inds[quant_index++]);
            }
        };


        double block_interpolation_1d(T *data, size_t begin, size_t end, size_t stride,
                                      const std::string &interp_func,
                                      const PredictorBehavior pb) {
            size_t n = (end - begin) / stride + 1;
            if (n <= 1) {
                return 0;
            }
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
                            quantize(d - data, *d, *(d - stride));
                    }
                } else {
                    for (size_t i = 1; i + 1 < n; i += 2) {
                        T *d = data + begin + i * stride;
                        recover(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                    }
                    if (n % 2 == 0) {
                        T *d = data + begin + (n - 1) * stride;
                            recover(d - data, *d, *(d - stride));
                    }
                }
            } else {
                if (pb == PB_predict_overwrite) {

                    T *d;
                    size_t i;
                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;

                            quantize(d - data, *d,
                                 interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));

                    }
                    d = data + begin + stride;
                    quantize(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                    d = data + begin + i * stride;
                    quantize(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));

                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        quantize(d - data, *d,  *(d - stride));
                    }

                } else {
                    T *d;

                    size_t i;
                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        recover(d - data, *d, interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                    }
                    d = data + begin + stride;
                    recover(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                    d = data + begin + i * stride;
                    recover(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        recover(d - data, *d,  *(d - stride));
                    }
                }
            }

            return predict_error;
        }

        template<uint NN = N>
        typename std::enable_if<NN == 1, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, size_t stride = 1) {
            return block_interpolation_1d(data, begin[0], end[0], stride, interp_func, pb);
        }

        template<uint NN = N>
        typename std::enable_if<NN == 2, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, size_t stride = 1) {
            double predict_error = 0;
            size_t stride2x = stride * 2;
            const std::array<int, N> dims = dimension_sequences[direction];
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]];
                predict_error += block_interpolation_1d(data, begin_offset,
                                                        begin_offset +
                                                        (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
                                                        stride * dimension_offsets[dims[0]], interp_func, pb);
            }
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]];
                predict_error += block_interpolation_1d(data, begin_offset,
                                                        begin_offset +
                                                        (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
                                                        stride * dimension_offsets[dims[1]], interp_func, pb);
            }
            return predict_error;
        }

        template<uint NN = N>
        typename std::enable_if<NN == 3, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, size_t stride = 1) {
            double predict_error = 0;
            size_t stride2x = stride * 2;
            const std::array<int, N> dims = dimension_sequences[direction];
            double tmp = quantizer.get_eb();
            // double reduction_factor = sqrt(4.4159889);
            double C = sqrt(4.4159889);
            double C1 = sqrt(1.640625);
            double C2 = C1*C1;

            quantizer.set_eb(tmp/C2);
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                           k * dimension_offsets[dims[2]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[0]] - begin[dims[0]]) *
                                                            dimension_offsets[dims[0]],
                                                            stride * dimension_offsets[dims[0]], interp_func, pb);
                }
            }

            quantizer.set_eb(tmp/C1);
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                          k * dimension_offsets[dims[2]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[1]] - begin[dims[1]]) *
                                                            dimension_offsets[dims[1]],
                                                            stride * dimension_offsets[dims[1]], interp_func, pb);
                }
            }
            quantizer.set_eb(tmp);
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                    size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                          begin[dims[2]] * dimension_offsets[dims[2]];

                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[2]] - begin[dims[2]]) *
                                                            dimension_offsets[dims[2]],
                                                            stride * dimension_offsets[dims[2]], interp_func, pb);
                }
            }
            return predict_error;
        }


        template<uint NN = N>
        typename std::enable_if<NN == 4, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, size_t stride = 1) {
            double predict_error = 0;
            size_t stride2x = stride * 2;
            max_error = 0;
            const std::array<int, N> dims = dimension_sequences[direction];
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                         t <= end[dims[3]]; t += stride2x) {
                        size_t begin_offset =
                                begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[0]] - begin[dims[0]]) *
                                                                dimension_offsets[dims[0]],
                                                                stride * dimension_offsets[dims[0]], interp_func, pb);
                    }
                }
            }
            max_error = 0;
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                         t <= end[dims[3]]; t += stride2x) {
                        size_t begin_offset =
                                i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[1]] - begin[dims[1]]) *
                                                                dimension_offsets[dims[1]],
                                                                stride * dimension_offsets[dims[1]], interp_func, pb);
                    }
                }
            }
            max_error = 0;
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                    for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                         t <= end[dims[3]]; t += stride2x) {
                        size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                              begin[dims[2]] * dimension_offsets[dims[2]] +
                                              t * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[2]] - begin[dims[2]]) *
                                                                dimension_offsets[dims[2]],
                                                                stride * dimension_offsets[dims[2]], interp_func, pb);
                    }
                }
            }

            max_error = 0;
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                    for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride : 0); k <= end[dims[2]]; k += stride) {
                        size_t begin_offset =
                                i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                begin[dims[3]] * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[3]] - begin[dims[3]]) *
                                                                dimension_offsets[dims[3]],
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
        std::vector<int> quant_inds;
        size_t quant_index = 0; // for decompress
        double max_error;
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;
        std::array<size_t, N> dimension_offsets;
        std::vector<std::array<int, N>> dimension_sequences;
        int direction_sequence_id;
        int current_level; 
        std::vector<size_t> my_level;
        std::vector<int> my_index; 
        std::vector<bool> my_sift_index;
        size_t sift_index=0;
        double reduction_factor;
        double real_eb_ratio;
        double var_percentage;
    };


};


#endif

