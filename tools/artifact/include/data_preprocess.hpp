#include <algorithm>
#include <iostream>
#include <vector>

template <class T>
inline T  normalization(std::vector<T> &data, T &min, T  &max) {
    max = data[0];
    min = data[0];
    for (int i; i < data.size(); i++) {
      if(data[i] > max) max = data[i];
      if(data[i] < min) min = data[i];
    }
    T range = max - min;
    if (range == 0) {
      return -1;
    }
    for (int i; i < data.size(); i++) {
      data[i] = (data[i] - min) / range;
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