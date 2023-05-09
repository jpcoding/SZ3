#include <algorithm>
#include <iostream>
#include <vector>

template <class T>
inline T  normalization(std::vector<T> &data, T &min, T  &max) {
    max = *std::max_element(data.begin(), data.end());
    min = *std::min_element(data.begin(), data.end());
    float range = max - min;
    for (auto &d : data) {
      d = (d - min) / range;
    }
    return range;
}

template <class T>
inline void normalization(std::vector<T> &data, T &min, T  &max, T &range) {
    max = *std::max_element(data.begin(), data.end());
    min = *std::min_element(data.begin(), data.end());
    range = max - min;
    for (auto &d : data) {
      d = (d - min) / range;
    }
}