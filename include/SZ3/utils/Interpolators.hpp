//
// Created by Kai Zhao on 9/1/20.
//

#ifndef SZ_INTERPOLATORS_HPP
#define SZ_INTERPOLATORS_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <complex> 
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "fftw3.h"
namespace SZ {
    template<class T>
    inline T interp_linear(T a, T b) {
        return (a + b) / 2;
    }

    template<class T>
    inline T interp_linear1(T a, T b) {
        return -0.5 * a + 1.5 * b;
    }

    template<class T>
    inline T interp_quad_1(T a, T b, T c) {
        return (3 * a + 6 * b - c) / 8;
    }

    template<class T>
    inline T interp_quad_2(T a, T b, T c) {
        return (-a + 6 * b + 3 * c) / 8;
    }

    template<class T>
    inline T interp_quad_3(T a, T b, T c) {
        return (3 * a - 10 * b + 15 * c) / 8;
    }

    template<class T>
    inline T interp_cubic(T a, T b, T c, T d) {
        return (-a + 9 * b + 9 * c - d) / 16;
    }

    template<class T>
    inline T interp_cubic_front(T a, T b, T c, T d) {
        return (5 * a + 15 * b - 5 * c + d) / 16;
    }

    template<class T>
    inline T interp_cubic_front_2(T a, T b, T c, T d) {
        return ( a + 6 * b - 4 * c + d) / 4;
    }

    template<class T>
    inline T interp_cubic_back_1(T a, T b, T c, T d) {
        return (a - 5 * b + 15 * c + 5 * d) / 16;
    }

    template<class T>
    inline T interp_cubic_back_2(T a, T b, T c, T d) {
        return (-5 * a + 21 * b - 35 * c + 35 * d) / 16;
    }

    template<class T>
    inline T interp_cubic2(T a, T b, T c, T d) {
        return (-3 * a + 23 * b + 23 * c - 3 * d) / 40;
    }

    template<class T>
    inline T interp_akima(T a, T b, T c, T d) {
        T t0 = 2 * b - a - c;
        T t1 = 2 * c - b - d;
        T abt0 = fabs(t0);
        T abt1 = fabs(t1);
        if (fabs(abt0 + abt1) > 1e-9) {
            return (b + c) / 2 + (t0 * abt1 + t1 * abt0) / 8 / (abt0 + abt1);
        } else {
            return (b + c) / 2;
        }
    }

    template<class T>
    inline T interp_pchip(T a, T b, T c, T d) {
        T pchip = (b + c) / 2;
        if ((b - a < 0) == (c - b < 0) && fabs(c - a) > 1e-9) {
            pchip += 1 / 4 * (b - a) * (c - b) / (c - a);
        }
        if ((c - b < 0) == (d - c < 0) && fabs(d - b) > 1e-9) {
            pchip -= 1 / 4 * (c - b) * (d - c) / (d - b);
        }
        return pchip;
    }
    

    template <class T>
    inline T trigonometric_interp(T a, T b, T c, T d) 
    {
        T x0 = (a* std::cos(0*0* 2*M_PI/4) + b* std::cos(0*1*2*M_PI/4) + c* std::cos(0*2* 2*M_PI/4) + d* std::cos(0*3*2*M_PI/4))/4.0;
        T x1 = (a* std::cos(1*0* 2*M_PI/4) + b* std::cos(1*1*2*M_PI/4) + c* std::cos(1*2* 2*M_PI/4) + d* std::cos(1*3*2*M_PI/4))/4.0;
        T x2 = (a* std::cos(2*0* 2*M_PI/4) + b* std::cos(2*1*2*M_PI/4) + c* std::cos(2*2* 2*M_PI/4) + d* std::cos(2*3*2*M_PI/4))/4.0;
        T x3 = (a* std::cos(3*0* 2*M_PI/4) + b* std::cos(3*1*2*M_PI/4) + c* std::cos(3*2* 2*M_PI/4) + d* std::cos(3*3*2*M_PI/4))/4.0;
        T result = x0+x1*cos(3*M_PI/4)+x2*cos(2*3*M_PI/4)+x3*cos(3*3*M_PI/4);
        return result;
    }

    // template<class T> 
    // inline void fft_interp(T *data_start, size_t data_offset, int interp_interval, int num_interpolators, SZ::LinearQuantizer<T> quantizer, 
    //             std::vector<int> & quant_inds )
    // {
    //     // current implementation only support 1D interpolation to expand by a factor of 2 
    //     T *data = data_start; 
    //     int interpolator_offset = 2*data_offset;
    //     int full_length = num_interpolators * 2;
    //     // create a new array to store result of fft coefficients 
    //     std::vector<std::complex<T>> X(num_interpolators,0);
    //     // perform fft one the strided data
    //     for (int i = 0 ; i < num_interpolators; i++)
    //     {
    //         std::complex<double> sum = 0;
    //         for (int j = 0; j < num_interpolators; j++)
    //         {
    //             sum += data[j*interpolator_offset] * std::exp(std::complex<double>(0,-2*M_PI*i*j/num_interpolators));
    //         }
    //         X[i] = sum;
    //     }
    //     // padding the X array to do inverse fft 
    //     // X[0:num_interpolators/2] = X[0:num_interpolators/2]
    //     // X[full_length/2:full_length] = 0
    //     // not going to copy the data, just overwrite the data
    //     int nyquist = num_interpolators/2; 
    //     T *current_data = data_start+data_offset; 
    //     for (int i = 1; i < full_length; i+=2)
    //     {
    //         double pred = 0;
    //         for (int j = 0; j <= nyquist; j++)
    //         {
    //             // the first half of X
    //             pred += (X[j] * std::exp(std::complex<double>(0,2*M_PI*i*j/full_length))).real();
    //             // the second half of X
    //             int second_half_index = j+num_interpolators;
    //             pred += (X[second_half_index] * std::exp(std::complex<double>(0,2*M_PI*(i)*(second_half_index)/full_length))).real();
    //         }
    //         pred = pred/num_interpolators;
    //         quant_inds.push_back(quantizer.quantize_and_overwrite(current_data, pred));
    //         current_data += data_offset*2;
    //     }
    // }

template <typename T>
int interp_fft_fftw(T *data_start, size_t data_offset, int interp_interval,
                 int num_interpolators, SZ::LinearQuantizer<T> quantizer, 
                std::vector<int> & quant_inds)
{
    // current implementation only support 1D interpolation to expand by a factor of 2
    T *data = data_start;
    int interpolator_offset = 2*data_offset;
    int full_length = num_interpolators * 2; // only interpolate the odd index inbetween
    fftw_complex *in, *out, *ifft_in, *ifft_out;
    fftw_plan p;
    // auto timer = SZ::Timer(); 
    // timer.start();
    in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_interpolators);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_interpolators);
    ifft_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * full_length);
    ifft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * full_length);
    for (size_t i = 0; i < num_interpolators; i++) {
      in[i][0] = data[i*interpolator_offset];
      in[i][1] = 0;
    }
    p = fftw_plan_dft_1d(num_interpolators, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    // std::cout << "fftw time = " << timer.stop() << std::endl;


    // timer.start();
    int nyquist_idx = num_interpolators/2; 
    for(int i = 0; i <= nyquist_idx; i++) {
        ifft_in[i][0] = out[i][0];
        ifft_in[i][1] = out[i][1];
    }
    for(int i = nyquist_idx+1; i < num_interpolators; i++) {
        ifft_in[i+num_interpolators][0] = out[i][0];
        ifft_in[i+num_interpolators][1] = out[i][1];
    }
    if(num_interpolators % 2 == 0) {
        ifft_in[nyquist_idx][0] = out[nyquist_idx][0]/2;
        ifft_in[nyquist_idx][1] = out[nyquist_idx][1]/2;
        ifft_in[nyquist_idx+num_interpolators][0] = ifft_in[nyquist_idx][0];
        ifft_in[nyquist_idx+num_interpolators][1] = ifft_in[nyquist_idx][1];
    }
    p = fftw_plan_dft_1d(full_length, ifft_in, ifft_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    // std::cout << "fftw ifft time = " << timer.stop() << std::endl;
    // // write the fft result to a file
    // std::ofstream fft_file("fftw_fft_in.txt");
    // for (size_t i = 0; i < full_length; i++) {
    //     fft_file << ifft_in[i][0] << "  " << ifft_in[i][1]<<  std::endl;
    // }
    // fft_file.close();

    // std::ofstream fft_file2("fftw_fft_out.txt");
    // for (size_t i = 0; i < full_length; i++) {
    //     fft_file2 << ifft_out[i][0]/full_length << "  " << ifft_out[i][1]/full_length<<  std::endl;
    // }
    // fft_file2.close();

    for (int i = 1; i < full_length-1; i+=2) {
        T pred = ifft_out[i][0]/full_length*2;
        quant_inds.push_back(quantizer.quantize_and_overwrite(data_start[i*data_offset], pred));
    }

    // free fftw
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    fftw_free(ifft_in);
    fftw_free(ifft_out);

    return 0;
}
}

#endif //SZ_INTERPOLATORS_HPP
