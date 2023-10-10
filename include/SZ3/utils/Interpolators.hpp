//
// Created by Kai Zhao on 9/1/20.
//

#ifndef SZ_INTERPOLATORS_HPP
#define SZ_INTERPOLATORS_HPP

#include <algorithm>
#include <array>
#include <cmath>
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
}

#endif //SZ_INTERPOLATORS_HPP
