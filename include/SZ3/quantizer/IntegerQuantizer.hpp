#ifndef _SZ_INTEGER_QUANTIZER_HPP
#define _SZ_INTEGER_QUANTIZER_HPP

#include <cmath>
#include <cstddef>
#include <cstring>
#include <cassert>
#include <iostream>
#include <math.h>
#include <vector>
#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/Timer.hpp"
#include <random>

namespace SZ {

    template<class T>
    class LinearQuantizer : public concepts::QuantizerInterface<T> {
    public:
        LinearQuantizer() : error_bound(1), error_bound_reciprocal(1), radius(32768) {}

        LinearQuantizer(double eb, int r = 32768) : error_bound(eb),
                                                    error_bound_reciprocal(1.0 / eb),
                                                    radius(r) {
            assert(eb != 0);
        }

        int get_radius() const { return radius; }

        double get_eb() const { return error_bound; }

        void set_eb(double eb) {
            error_bound = eb;
            error_bound_reciprocal = 1.0 / eb;
        }

        // quantize the data with a prediction value, and returns the quantization index
        int quantize(T data, T pred) {
            T diff = data - pred;
            int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
            if (quant_index < this->radius * 2) {
                quant_index >>= 1;
                int half_index = quant_index;
                quant_index <<= 1;
                int quant_index_shifted;
                if (diff < 0) {
                    quant_index = -quant_index;
                    quant_index_shifted = this->radius - half_index;
                } else {
                    quant_index_shifted = this->radius + half_index;
                }
                T decompressed_data = pred + quant_index * this->error_bound;
                if (fabs(decompressed_data - data) > this->error_bound) {
                    return 0;
                } else {
                    return quant_index_shifted;
                }
            } else {
                return 0;
            }
        }

        // quantize the data with a prediction value, and returns the quantization index and the decompressed data
        // int quantize(T data, T pred, T& dec_data);
        int quantize_and_overwrite(T &data, T pred,bool save_unpred=true) {
            T diff = data - pred;
            int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
            if (quant_index < this->radius * 2) {
                quant_index >>= 1; 
                int half_index = quant_index;
                quant_index <<= 1;  
                int quant_index_shifted;
                if (diff < 0) {
                    quant_index = -quant_index;
                    quant_index_shifted = this->radius - half_index; 
                } else {
                    quant_index_shifted = this->radius + half_index; // 1+ radius
                }
                T decompressed_data = pred + quant_index * this->error_bound ;
                if (fabs(decompressed_data - data) > this->error_bound) {
                    // printf("data: %.7f, pred: %.7f, diff: %.7f, quant_index: %d, decompressed_data: %.7f, eb: %.7f, error: %.7f \n", data, pred,
                        //    diff, quant_index, decompressed_data, this->error_bound, fabs(decompressed_data - data));
                    if(save_unpred)
                        unpred.push_back(data);
                    return 0;
                } else {
                    data = decompressed_data;
                    return quant_index_shifted;
                }
                return quant_index_shifted;
            } else {
                if(save_unpred)
                    unpred.push_back(data);
                return 0;
            }
        }


        int quantize_and_overwrite_prob(T &data, T pred, T probability, bool save_unpred=true) {
            T diff = data - pred;
            int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
            // 1 or 2 
            if (quant_index < this->radius * 2) {
                quant_index >>= 1; 
                int half_index = quant_index; 
                // introcude random quantization here
                T distance = fabs(diff) - half_index*2 * this->error_bound;
                int direction = distance>0? 1:-1; // direction of quantization, if orig is 2, 1 means toward 4, -1 means toward 0 
                int rand_quant = fabs(distance)/(2.0*this->error_bound) > probability? 1:0; // magnitude of quantization 0 or 1,1 means quantization to a farther integer
                quant_index  += direction*rand_quant;
                half_index = quant_index;
                quant_index <<= 1;

                int quant_index_shifted;
                if (diff < 0) {
                    quant_index = -quant_index;
                    quant_index_shifted = this->radius - half_index; 
                } else {
                    quant_index_shifted = this->radius + half_index; // 1+ radius
                }
                T decompressed_data = pred + quant_index * this->error_bound ;
                if (fabs(decompressed_data - data) > 2*this->error_bound) {
                    // printf("data: %.7f, pred: %.7f, diff: %.7f, quant_index: %d, decompressed_data: %.7f, eb: %.7f, error: %.7f \n", data, pred,
                        //    diff, quant_index, decompressed_data, this->error_bound, fabs(decompressed_data - data));
                    if(save_unpred)
                        unpred.push_back(data);
                    return 0;
                } else {
                    data = decompressed_data;
                    return quant_index_shifted;
                }
                return quant_index_shifted;
            } else {
                if(save_unpred)
                    unpred.push_back(data);
                return 0;
            }
        }


        int quantize_and_overwrite_with_noise(T &data, T pred, T  rand_noise, bool save_unpred=true) {
            T diff = data - pred;
            int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
            if (quant_index < this->radius * 2) {
                quant_index >>= 1;
                int half_index = quant_index;
                quant_index <<= 1;
                int quant_index_shifted;
                if (diff < 0) {
                    quant_index = -quant_index;
                    quant_index_shifted = this->radius - half_index;
                } else {
                    quant_index_shifted = this->radius + half_index;
                }
                T decompressed_data = pred + quant_index * this->error_bound+ rand_noise; 
                T updated_bound = this->error_bound + std::fabs(rand_noise);
                if (fabs(decompressed_data - data) > updated_bound) {
                    if(save_unpred)
                        unpred.push_back(data);
                    return 0;
                } else {
                    data = decompressed_data;
                    return quant_index_shifted;
                }
                data = decompressed_data;
                return quant_index_shifted;
            } else {
                if(save_unpred)
                    unpred.push_back(data);
                return 0;
            }
        }


        int stochastic_quantize_and_overwrite(T &data, T pred, int rand_index , bool save_unpred=true)
        {
            T diff = data - pred; 
            // SZ::Timer timer;
            // timer.start();
            // srand(index);
            // int quantization_rand = rand() & 1;
            int quantization_rand = rand_index; // rand_index = 1 or 0 
            // timer.stop("random generator");
            int quant_index = (int) std::floor(std::fabs(diff)/this->error_bound);
            int shifted_quant_index;
            // int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal); 
            if (quant_index < this->radius )
            {       
                if (diff < 0)
                {
                    quant_index = -quant_index-quantization_rand;
                    shifted_quant_index = quant_index + this->radius;
                }
                else
                {
                    quant_index = quant_index+quantization_rand;
                    shifted_quant_index = this->radius + quant_index;
                }
                T decompressed_data = pred + quant_index * this->error_bound;
                if (fabs(decompressed_data - data) > this->error_bound) {
                    if(save_unpred)
                        unpred.push_back(data);
                    std::cout<<"unpred"<<std::endl;
                    return 0;
                } else {
                    data = decompressed_data;
                    return shifted_quant_index;
                }

            }else {
                if(save_unpred)
                    unpred.push_back(data);
                return 0;
            }
        }

        /*
        Do not use this
        */
        T stochastic_recover(T pred, int shifted_quant_index, size_t data_index)
        {
            if (shifted_quant_index !=0) {
                // srand(data_index);
                // int quantization_rand = rand() & 1;
                int quant_index = shifted_quant_index - this->radius;
                // if (quant_index < 0)
                //     quant_index = quant_index - quantization_rand;
                // else
                //     quant_index = quant_index + quantization_rand;
                T decompressed_data = pred + quant_index * this->error_bound;
                return decompressed_data;
            } else {
                return recover_unpred();
            }
        }


        int quantize_and_overwrite(T ori, T pred, T &dest,bool save_unpred=true) {
            T diff = ori - pred;
            int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
            if (quant_index < this->radius * 2) {
                quant_index >>= 1;
                int half_index = quant_index;
                quant_index <<= 1;
                int quant_index_shifted;
                if (diff < 0) {
                    quant_index = -quant_index;
                    quant_index_shifted = this->radius - half_index;
                } else {
                    quant_index_shifted = this->radius + half_index;
                }
                T decompressed_data = pred + quant_index * this->error_bound;
                if (fabs(decompressed_data - ori) > this->error_bound) {
                    if(save_unpred)
                        unpred.push_back(ori);
                    dest = ori;
                    return 0;
                } else {
                    dest = decompressed_data;
                    return quant_index_shifted;
                }
            } else {
                if(save_unpred)
                    unpred.push_back(ori);
                dest = ori;
                return 0;
            }
        }


        void insert_unpred(T ori){
            unpred.push_back(ori);
        }

        void insert_unpred_idx(size_t ori){
            unpred_idx.push_back(ori);
        }

        void print_unpred(){
            for(auto x:unpred)
                std::cout<<x<<std::endl;
        }

        // recover the data using the quantization index
        T recover(T pred, int quant_index) {
            if (quant_index) {
                return recover_pred(pred, quant_index);
            } else {
                return recover_unpred();
            }
        }

        T recover_with_noise(T pred, int quant_index, T rand_noise) {
            if (quant_index) {
                return recover_pred(pred, quant_index)+ rand_noise;
            } else {
                return recover_unpred();
            }
        }



        T stochastic_recover_pred(T pred, int quant_index) {
            // return pred + quant_index * this->error_bound;
            // int index = ; 
            return pred + (quant_index + (rand() & 1)) * this->error_bound;
        }

        T recover_pred(T pred, int quant_index) {
            return pred + 2 * (quant_index - this->radius) * this->error_bound;
        }

        T recover_unpred() {
            return unpred[index++];
        }

        size_t recover_unpred_idx() {
            return unpred_idx[idx_index++];
        }

        size_t current_unpred_idx() {
            return unpred_idx[idx_index];
        }

        std::vector<size_t> get_unpred_idx()
        {
            return unpred_idx;
        }


        size_t size_est() {
            return unpred.size() * sizeof(T);
        }

        void save(unsigned char *&c) const {
            // std::string serialized(sizeof(uint8_t) + sizeof(T) + sizeof(int),0);
            c[0] = 0b00000010;
            c += 1;
            // std::cout << "saving eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            *reinterpret_cast<double *>(c) = this->error_bound;
            c += sizeof(double);
            *reinterpret_cast<int *>(c) = this->radius;
            c += sizeof(int);
            *reinterpret_cast<size_t *>(c) = unpred.size();
            c += sizeof(size_t);
            memcpy(c, unpred.data(), unpred.size() * sizeof(T));
            c += unpred.size() * sizeof(T);
            *reinterpret_cast<size_t *>(c) = unpred_idx.size();
            c += sizeof(size_t);
            memcpy(c, unpred_idx.data(), unpred_idx.size() * sizeof(size_t));
            c += unpred_idx.size() * sizeof(size_t);

            std::cout << "unpred size" << unpred.size() << std::endl;


        };

        void load(const unsigned char *&c, size_t &remaining_length) {
            assert(remaining_length > (sizeof(uint8_t) + sizeof(T) + sizeof(int)));
            c += sizeof(uint8_t);
            remaining_length -= sizeof(uint8_t);
            this->error_bound = *reinterpret_cast<const double *>(c);
            this->error_bound_reciprocal = 1.0 / this->error_bound;
            c += sizeof(double);
            this->radius = *reinterpret_cast<const int *>(c);
            c += sizeof(int);
            size_t unpred_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            this->unpred = std::vector<T>(reinterpret_cast<const T *>(c), reinterpret_cast<const T *>(c) + unpred_size);
            c += unpred_size * sizeof(T);

            size_t unpred_idx_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            this->unpred_idx = std::vector<size_t>(reinterpret_cast<const size_t *>(c), reinterpret_cast<const size_t *>(c) + unpred_idx_size);
            c += unpred_idx_size * sizeof(size_t);
            // std::cout << "loading: eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            // reset index
            index = 0;
            idx_index =0;
        }

        void print() {
            printf("[IntegerQuantizer] error_bound = %.8G, radius = %d, unpred = %lu\n", error_bound, radius, unpred.size());
        }

        void clear() {
            unpred.clear();
            index = 0;
        }


        virtual void postcompress_data() {
        }

        virtual void postdecompress_data() {
        }

        virtual void precompress_data() {};

        virtual void predecompress_data() {};


    private:
        std::vector<T> unpred;
        std::vector<size_t> unpred_idx;
        size_t index = 0; // used in decompression only
        size_t idx_index = 0; 

        double error_bound;
        double error_bound_reciprocal;
        int radius; // quantization interval radius
    };

}
#endif
