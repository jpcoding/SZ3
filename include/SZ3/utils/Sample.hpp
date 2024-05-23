
#ifndef SZ_SAMPLE_HPP
#define SZ_SAMPLE_HPP

#include <vector>

template <typename T>
int level_sampleing(
    T *data, T *sample_data, size_t num_elements, size_t sample_rate)
{
  size_t sample_num = num_elements / sample_rate;
  for (size_t i = 0; i < sample_num; i++) {
    sample_data[i] = data[i * sample_rate];
  }
  return sample_num;
}

#endif
