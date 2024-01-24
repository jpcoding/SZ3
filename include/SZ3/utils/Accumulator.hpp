//
// Created by Pu Jiao on 1/22/2023.
//

#ifndef SZ_Accumulator_HPP
#define SZ_Accumulator_HPP

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace SZ {


template <typename T>
class Accumulator {

    Accumulator() {};
    Accumulator(size_t possible_size) {
        collected_data.reserve(possible_size);
    }

    ~Accumulator() {};

    void collect_data(T data) {
        collected_data.push_back(data);
    }

    T get_mean() {
        double sum = 0;
        for (size_t i = 0; i < collected_data.size(); i++) {
            sum += collected_data[i];
        }
        return (T) (sum / (1.0*collected_data.size()));
    }

    


private:
    size_t count;
    std::vector<T> collected_data;


};
}
#endif //SZ_Accumulator_HPP



