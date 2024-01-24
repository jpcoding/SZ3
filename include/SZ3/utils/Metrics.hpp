#include "Config.hpp"
#include "SZ3/def.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace SZ {
namespace METRICS {

template <typename T>
double CALCULATE_PSNR_MSE(const T *odata, const T *ddata, size_t num_elements) {
  double sum = 0;
  for (size_t i = 0; i < num_elements; i++) {
    double diff = odata[i] - ddata[i];
    sum += diff * diff;
  }
  return sum / num_elements;
}

template <typename T>
double CALCULATE_PSNR(const T *odata, const T *ddata, size_t num_elements) {
  T max = *std::max_element(odata, odata + num_elements);
  T min = *std::min_element(odata, odata + num_elements);
  double rng = max - min;
  double mse = CALCULATE_PSNR_MSE(odata, ddata, num_elements);
  return 20 * log10(rng) - 10 * log10(mse);
}

} // namespace METRICS
} // namespace SZ
