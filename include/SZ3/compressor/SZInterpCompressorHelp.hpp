#ifndef  SZ_INTERPOLATION_COMPRESSOR_HELP_HPP
#define  SZ_INTERPOLATION_COMPRESSOR_HELP_HPP
#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/predictor/ComposedPredictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <exception>
#include <random>
#include <sys/types.h>
#include <vector>
/*
This file is only meant for the helper functions for
SZInterpolationCompressor.hpp Do not use the functions in this files in
eleswhere.
*/

namespace SZ {

template <typename T>
void compute_auxilliary_data(const SZ::Config conf, T *data,
                             const int detection_block_size,
                             const size_t num_elements,
                             size_t &num_detection_block,
                             std::vector<uchar> &flushed_block_id,
                             std::vector<uchar> &flushed_block,
                             std::vector<uchar> &significant_block_id,
                             std::vector<uchar> &significant_block) {
  double detection_threshold = conf.detection_threshold;
  bool block_sift_on = conf.block_sift_on;
  bool block_flush_on = conf.block_flush_on;
  bool block_iso_on = conf.block_iso_on;
  double eb_input = conf.absErrorBound;
  double isovalue = conf.block_isovalue;
  uint sift_mode = conf.block_sift_mode;
  uint N = conf.N;
  if (detection_block_size == 1) {
    num_detection_block = num_elements;
    if (block_flush_on == 1 && block_sift_on == 1) {
      // Timer timer;
      // timer.start();
      flushed_block_id = std::vector<uchar>(num_elements, 0);
      significant_block_id = std::vector<uchar>(num_elements, 0);
      double threshold;
      // timer.stop("SIFT: define variables ");
      { // destroy the copy after use
        std::vector<T> block_significance_tmp(data, data + num_elements);
        // timer.start();
        std::nth_element(block_significance_tmp.begin(),
                         block_significance_tmp.begin() +
                             (size_t)(detection_threshold * num_elements),
                         block_significance_tmp.end());
        // timer.stop("SIFT: select nth ");
        threshold = block_significance_tmp[(size_t)(detection_threshold *
                                                    num_elements)];
      }
      // timer.start();
      for (int i = 0; i < num_elements; i++) {
        if (data[i] > threshold) {
          significant_block_id[i] = 1;
        }
        if (fabs(data[i]) < eb_input * 0.1) {
          flushed_block_id[i] = 1;
        }
      }
      // timer.stop("SIFT: sift");
      significant_block = significant_block_id;
      flushed_block = flushed_block_id;

    } else if (block_flush_on == 1 && block_sift_on == 0) {
      flushed_block_id = std::vector<uchar>(num_elements, 0);
      for (int i = 0; i < num_elements; i++) {
        T *data_pos = data + i;
        if (fabs(*data_pos) < eb_input * 0.1) {
          flushed_block_id[i] = 1;
        }
      }
      flushed_block = flushed_block_id;

    } else if (block_flush_on == 0 && block_sift_on == 1) {
      significant_block_id = std::vector<uchar>(num_elements, 0);

      double threshold;
      Timer timer;
      { // destroy the copy after use

        timer.start();
        std::vector<T> block_significance_tmp(data, data + num_elements);
        // std::sort(block_significance_tmp.begin(),
        // block_significance_tmp.end());
        std::nth_element(block_significance_tmp.begin(),
                         block_significance_tmp.begin() +
                             (size_t)(detection_threshold * num_elements),
                         block_significance_tmp.end());
        threshold = block_significance_tmp[(size_t)(detection_threshold *
                                                    num_elements)];
        timer.stop("SIFT:compute threshold ");
      }
      timer.start();
      for (int i = 0; i < num_elements; i++) {
        T *data_pos = data + i;
        if (*data_pos > threshold) {
          significant_block_id[i] = 1;
        }
      }
      timer.stop("SIFT: lebel ");
      significant_block = significant_block_id;
    }
  } else if (N == 2) {
    int block_size = detection_block_size;
    const auto &dims = conf.dims;
    int nx = (int)ceil(dims[0] * 1.0 / block_size);
    int ny = (int)ceil(dims[1] * 1.0 / block_size);
    flushed_block = std::vector<uchar>(num_elements, 0);
    flushed_block_id = std::vector<uchar>(nx * ny, 0);
    num_detection_block = nx * ny;
    Timer timer;
    timer.start();
    // compute statistics
    int block_id = 0;
    int flushed_count = 0;
    // using variance or value range for significant blocks
    std::vector<double> block_significance;
    T *x_data_pos = data;
    for (int i = 0; i < nx; i++) {
      int block_size_i = (i + 1) * block_size > dims[0]
                             ? dims[0] - i * block_size
                             : block_size;
      T *y_data_pos = x_data_pos;
      for (int j = 0; j < ny; j++) {
        int block_size_j = (j + 1) * block_size > dims[1]
                               ? dims[1] - j * block_size
                               : block_size;
        T block_max = -std::numeric_limits<T>::max();
        T block_min = std::numeric_limits<T>::max();
        T block_distance_to_isovalue = std::numeric_limits<T>::max();
        // compute average
        double sum = 0;
        T *xx_data_pos = y_data_pos;
        for (int ii = 0; ii < block_size_i; ii++) {
          T *yy_data_pos = xx_data_pos;
          for (int jj = 0; jj < block_size_j; jj++) {
            sum += *yy_data_pos;
            if (*yy_data_pos > block_max)
              block_max = *yy_data_pos;
            if (*yy_data_pos < block_min)
              block_min = *yy_data_pos;
            yy_data_pos++;
          }
          xx_data_pos += dims[1];
        }
        T max_abs_val = std::max(fabs(block_max), fabs(block_min));
        T block_range = block_max - block_min;
        if (block_flush_on && max_abs_val < eb_input * 0.1) {
          T *xx_data_pos = y_data_pos;
          for (int ii = 0; ii < block_size_i; ii++) {
            T *yy_data_pos = xx_data_pos;
            for (int jj = 0; jj < block_size_j; jj++) {
              *yy_data_pos = 0;
              flushed_block[yy_data_pos - data] = 1;
              yy_data_pos++;
            }
            xx_data_pos += dims[1];
          }
          flushed_block_id[block_id] = 1;
          flushed_count++;
          // skip the variance
          if (block_sift_on)
            block_significance.push_back(0);
        } else if (block_sift_on) {
          // compute block significance
          if (sift_mode == SZ::BLOCK_SIFT_MODE::BLOCK_MAX) {
            block_significance.push_back(max_abs_val);
          } else if (sift_mode == SZ::BLOCK_SIFT_MODE::RANGE) {
            block_significance.push_back(block_range);
          } else if (sift_mode == SZ::BLOCK_SIFT_MODE::VARIANCE) {
            double average = sum / block_size_i / block_size_j;
            double variance = 0;
            T *xx_data_pos = y_data_pos;
            for (int ii = 0; ii < block_size_i; ii++) {
              T *yy_data_pos = xx_data_pos;
              for (int jj = 0; jj < block_size_j; jj++) {
                variance += (*yy_data_pos - average) * (*yy_data_pos - average);
                yy_data_pos++;
              }
              xx_data_pos += dims[1];
            }
            variance /= (block_size_i * block_size_j);
            block_significance.push_back(variance);
          }
        }
        block_id++;
        y_data_pos += block_size;
      }
      x_data_pos += block_size * dims[1];
    }
    std::cout << "flushed_count = " << flushed_count
              << ", percent = " << flushed_count * 1.0 / (nx * ny) << std::endl;
    if (block_sift_on == true) {
      double threshold, percent;
      percent = detection_threshold;
      {
        auto block_significance_tmp(block_significance);
        // std::sort(block_significance_tmp.begin(),
        // block_significance_tmp.end());
        std::nth_element(block_significance_tmp.begin(),
                         block_significance_tmp.begin() +
                             (size_t)(detection_threshold * num_elements),
                         block_significance_tmp.end());
        threshold = block_significance_tmp[(size_t)(percent *
                                                    block_significance.size())];
      }
      std::cout << percent * 100 << "% threshold = " << threshold << std::endl;
      significant_block = std::vector<uchar>(num_elements, 0);
      significant_block_id = std::vector<uchar>(nx * ny, 0);
      block_id = 0;
      x_data_pos = data;
      for (int i = 0; i < nx; i++) {
        int block_size_i = (i + 1) * block_size > dims[0]
                               ? dims[0] - i * block_size
                               : block_size;
        T *y_data_pos = x_data_pos;
        for (int j = 0; j < ny; j++) {
          int block_size_j = (j + 1) * block_size > dims[1]
                                 ? dims[1] - j * block_size
                                 : block_size;
          if (block_significance[block_id] > threshold) {
            T *xx_data_pos = y_data_pos;
            for (int ii = 0; ii < block_size_i; ii++) {
              T *yy_data_pos = xx_data_pos;
              for (int jj = 0; jj < block_size_j; jj++) {
                significant_block[yy_data_pos - data] = 1;
                yy_data_pos++;
              }
              xx_data_pos += dims[1];
            }
            significant_block_id[block_id] = 1;
          }
          block_id++;
          y_data_pos += block_size;
        }
        x_data_pos += block_size * dims[1];
      }
    }
    timer.stop("detect");
  } else if (N == 3) {
    int block_size = detection_block_size;
    const auto &dims = conf.dims;
    int nx = (int)ceil(dims[0] * 1.0 / block_size);
    int ny = (int)ceil(dims[1] * 1.0 / block_size);
    int nz = (int)ceil(dims[2] * 1.0 / block_size);
    flushed_block = std::vector<uchar>(num_elements, 0);
    flushed_block_id = std::vector<uchar>(nx * ny * nz, 0);
    num_detection_block = nx * ny * nz;
    // int actual_block_size;
    Timer timer;
    timer.start();
    // compute statistics
    int block_id = 0;
    int flushed_count = 0;
    // using variance for significant blocks
    std::vector<double> block_significance;
    block_significance.reserve(nx * ny * nz);
    T *x_data_pos = data;
    for (int i = 0; i < nx; i++) {
      int block_size_i = (i + 1) * block_size > dims[0]
                             ? dims[0] - i * block_size
                             : block_size;
      T *y_data_pos = x_data_pos;
      for (int j = 0; j < ny; j++) {
        int block_size_j = (j + 1) * block_size > dims[1]
                               ? dims[1] - j * block_size
                               : block_size;
        T *z_data_pos = y_data_pos;
        for (int k = 0; k < nz; k++) {
          int block_size_k = (k + 1) * block_size > dims[2]
                                 ? dims[2] - k * block_size
                                 : block_size;
          T block_max = -std::numeric_limits<T>::max();
          T block_min = std::numeric_limits<T>::max();
          T block_isovalue_distance = std::numeric_limits<T>::max();
          // actual_block_size = 0;
          // compute average
          double sum = 0;
          T *xx_data_pos = z_data_pos;
          for (int ii = 0; ii < block_size_i; ii++) {
            T *yy_data_pos = xx_data_pos;
            for (int jj = 0; jj < block_size_j; jj++) {
              T *zz_data_pos = yy_data_pos;
              for (int kk = 0; kk < block_size_k; kk++) {
                sum += *zz_data_pos;
                if (*zz_data_pos > block_max)
                  block_max = *zz_data_pos;
                if (*zz_data_pos < block_min)
                  block_min = *zz_data_pos;
                if (block_iso_on == true) {
                  block_isovalue_distance =
                      std::min<T>(fabs(block_isovalue_distance),
                                  fabs(isovalue - *zz_data_pos));
                }
                zz_data_pos++;
              }
              yy_data_pos += dims[2];
            }
            xx_data_pos += dims[1] * dims[2];
          }
          T max_abs_val = std::max(fabs(block_max), fabs(block_min));

          T block_range = block_max - block_min;
          if (block_flush_on == true && max_abs_val < eb_input * 0.1) {
            T *xx_data_pos = z_data_pos;
            for (int ii = 0; ii < block_size_i; ii++) {
              T *yy_data_pos = xx_data_pos;
              for (int jj = 0; jj < block_size_j; jj++) {
                T *zz_data_pos = yy_data_pos;
                for (int kk = 0; kk < block_size_k; kk++) {
                  *zz_data_pos = 0;
                  flushed_block[zz_data_pos - data] = 1;
                  zz_data_pos++;
                }
                yy_data_pos += dims[2];
              }
              xx_data_pos += dims[1] * dims[2];
            }
            flushed_block_id[block_id] = 1;
            flushed_count++;
            // skip the variance
            if (block_sift_on == true) {
              block_significance.push_back(0);
            }
          } else if (block_sift_on == true) {
            // compute variance
            if (sift_mode == SZ::BLOCK_SIFT_MODE::BLOCK_MAX) {
              block_significance.push_back(max_abs_val);
            } else if (sift_mode == SZ::BLOCK_SIFT_MODE::RANGE) {
              block_significance.push_back(block_range);
            } else if (sift_mode == SZ::BLOCK_SIFT_MODE::ISOVALUE) {
              block_significance.push_back(-block_isovalue_distance);
            } else if (sift_mode == SZ::BLOCK_SIFT_MODE::VARIANCE) {
              double average = sum / block_size_i / block_size_j / block_size_k;
              double variance = 0;
              T *xx_data_pos = z_data_pos;
              for (int ii = 0; ii < block_size_i; ii++) {
                T *yy_data_pos = xx_data_pos;
                for (int jj = 0; jj < block_size_k; jj++) {
                  T *zz_data_pos = yy_data_pos;
                  for (int kk = 0; kk < block_size_k; kk++) {
                    variance +=
                        (*zz_data_pos - average) * (*zz_data_pos - average);
                    zz_data_pos++;
                  }
                  yy_data_pos += dims[2];
                }
                xx_data_pos += dims[1] * dims[2];
              }
              variance /= (block_size_i * block_size_j * block_size_k);
              block_significance.push_back(variance);
            }
          }
          block_id++;
          z_data_pos += block_size;
        }
        y_data_pos += block_size * dims[2];
      }
      x_data_pos += block_size * dims[1] * dims[2];
    }
    std::cout << "flushed_count = " << flushed_count
              << ", percent = " << flushed_count * 1.0 / (nx * ny * nz)
              << std::endl;
    if (block_sift_on == true) {
      double threshold;
      double percent = detection_threshold;

      { // destroy the copy after use
        // Timer timer;
        // timer.start();
        auto block_significance_tmp(block_significance);
        // timer.stop("SIFT: copy data");
        // timer.start();
        // std::sort(block_significance_tmp.begin(),
        // block_significance_tmp.end());
        std::nth_element(block_significance_tmp.begin(),
                         block_significance_tmp.begin() +
                             (size_t)(detection_threshold * num_elements),
                         block_significance_tmp.end());
        // timer.stop("SIFT: sort time");
        threshold = block_significance_tmp[(size_t)(percent *
                                                    block_significance.size())];
      }
      std::cout << percent * 100 << "% " << SZ::BLOCK_SIFT_MODE_STR[sift_mode]
                << " threshold = " << threshold << std::endl;
      significant_block = std::vector<uchar>(num_elements, 0); // per datapoint
      significant_block_id = std::vector<uchar>(nx * ny * nz, 0); // per block
      block_id = 0;
      x_data_pos = data;
      for (int i = 0; i < nx; i++) {
        int block_size_i = (i + 1) * block_size > dims[0]
                               ? dims[0] - i * block_size
                               : block_size;
        T *y_data_pos = x_data_pos;
        for (int j = 0; j < ny; j++) {
          int block_size_j = (j + 1) * block_size > dims[1]
                                 ? dims[1] - j * block_size
                                 : block_size;
          T *z_data_pos = y_data_pos;
          for (int k = 0; k < nz; k++) {
            int block_size_k = (k + 1) * block_size > dims[2]
                                   ? dims[2] - k * block_size
                                   : block_size;
            if (block_significance[block_id] > threshold) {
              T *xx_data_pos = z_data_pos;
              for (int ii = 0; ii < block_size_i; ii++) {
                T *yy_data_pos = xx_data_pos;
                for (int jj = 0; jj < block_size_j; jj++) {
                  T *zz_data_pos = yy_data_pos;
                  for (int kk = 0; kk < block_size_k; kk++) {
                    significant_block[zz_data_pos - data] = 1;
                    zz_data_pos++;
                  }
                  yy_data_pos += dims[2];
                }
                xx_data_pos += dims[1] * dims[2];
              }
              significant_block_id[block_id] = 1;
            }
            block_id++;
            z_data_pos += block_size_k;
          }
          y_data_pos += block_size_j * dims[2];
        }
        x_data_pos += block_size_i * dims[1] * dims[2];
      }
    }
    timer.stop("detect");
  }
}

template <typename T>
size_t compute_auxilliary_data_decompress(
    T* data, uint N, size_t *global_dims,
    const int detection_block_size, 
    const size_t num_elements,
    const bool block_sift_on, 
    const bool block_flush_on,
    const bool block_iso_on,
    std::vector<uchar> &flushed_block_id, 
    std::vector<uchar> &flushed_block,
    std::vector<uchar> &significant_block_id,
    std::vector<uchar> &significant_block) {
  size_t num_flushed_elements;
  std::cout << "compute aux " << std::endl;

  std::vector<size_t> global_dimensions(N);
  for (int i = 0; i < N; i++) {
    global_dimensions[i] = global_dims[i];
  }

  // special case when blocksize = 1
  if (detection_block_size == 1) {
    if (block_sift_on == 1) {
      significant_block = significant_block_id;
      num_flushed_elements = 0;
    }
    if (block_flush_on == 1) {
      num_flushed_elements = 0;
      flushed_block = flushed_block_id;
      for (int i = 0; i < num_elements; i++) {
        if (flushed_block[i])
          num_flushed_elements++;
      }
    }
    std::cout << "num_flushed_elements " << num_flushed_elements << std::endl;
    return num_flushed_elements;
  }

  if (N == 2) {
    int block_size = detection_block_size;
    const auto &dims = global_dimensions;
    int nx = (int)ceil(dims[0] * 1.0 / block_size);
    int ny = (int)ceil(dims[1] * 1.0 / block_size);
    flushed_block = std::vector<uchar>(num_elements, 0);
    significant_block = std::vector<uchar>(num_elements, 0);
    int block_id = 0;
    num_flushed_elements = 0;
    const T *x_data_pos = data;
    for (int i = 0; i < nx; i++) {
      int block_size_i = (i + 1) * block_size > dims[0]
                             ? dims[0] - i * block_size
                             : block_size;
      const T *y_data_pos = x_data_pos;
      for (int j = 0; j < ny; j++) {
        int block_size_j = (j + 1) * block_size > dims[1]
                               ? dims[1] - j * block_size
                               : block_size;
        if (block_flush_on && block_sift_on && flushed_block_id[block_id] &&
            significant_block_id[block_id]) {
          int actual_block_size = 0;
          const T *xx_data_pos = y_data_pos;
          for (int ii = 0; ii < block_size_i; ii++) {
            const T *yy_data_pos = xx_data_pos;
            for (int jj = 0; jj < block_size_j; jj++) {
              flushed_block[yy_data_pos - data] = 1;
              actual_block_size++;
              significant_block[yy_data_pos - data] = 1;
              yy_data_pos += 1;
            }
            xx_data_pos += dims[1];
          }
          num_flushed_elements += actual_block_size;
        } else if (block_flush_on && flushed_block_id[block_id]) {
          int actual_block_size = 0;
          const T *xx_data_pos = y_data_pos;
          for (int ii = 0; ii < block_size_i; ii++) {
            const T *yy_data_pos = xx_data_pos;
            for (int jj = 0; jj < block_size_j; jj++) {
              flushed_block[yy_data_pos - data] = 1;
              yy_data_pos += 1;
            }
            xx_data_pos += dims[1];
          }
          num_flushed_elements += actual_block_size;
        } else if (block_sift_on && significant_block_id[block_id]) {
          const T *xx_data_pos = y_data_pos;
          for (int ii = 0; ii < block_size_i; ii++) {
            const T *yy_data_pos = xx_data_pos;
            for (int jj = 0; jj < block_size_j; jj++) {
              significant_block[yy_data_pos - data] = 1;
              yy_data_pos += 1;
            }
            xx_data_pos += dims[1];
          }
        }
        block_id++;
        y_data_pos += block_size_j;
      }
      x_data_pos += block_size_i * dims[1];
    }
  } else if (N == 3) {
    int block_size = detection_block_size;
    const auto &dims = global_dimensions;
    int nx = (int)ceil(1.0 * dims[0] / block_size);
    int ny = (int)ceil(1.0 * dims[1] / block_size);
    int nz = (int)ceil(1.0 * dims[2] / block_size);
    flushed_block = std::vector<uchar>(num_elements, 0);
    significant_block = std::vector<uchar>(num_elements, 0);
    int block_id = 0;
    num_flushed_elements = 0;
    const T *x_data_pos = data;
    for (int i = 0; i < nx; i++) {
      int block_size_i = (i + 1) * block_size > dims[0]
                             ? dims[0] - i * block_size
                             : block_size;
      const T *y_data_pos = x_data_pos;
      for (int j = 0; j < ny; j++) {
        int block_size_j = (j + 1) * block_size > dims[1]
                               ? dims[1] - j * block_size
                               : block_size;
        const T *z_data_pos = y_data_pos;
        for (int k = 0; k < nz; k++) {
          int block_size_k = (k + 1) * block_size > dims[2]
                                 ? dims[2] - k * block_size
                                 : block_size;
          if (block_flush_on && block_sift_on && flushed_block_id[block_id] &&
              significant_block_id[block_id]) {
            int actual_block_size = 0;
            const T *xx_data_pos = z_data_pos;
            for (int ii = 0; ii < block_size_i; ii++) {
              const T *yy_data_pos = xx_data_pos;
              for (int jj = 0; jj < block_size_j; jj++) {
                const T *zz_data_pos = yy_data_pos;
                for (int kk = 0; kk < block_size_k; kk++) {
                  flushed_block[zz_data_pos - data] = 1;
                  significant_block[zz_data_pos - data] = 1;
                  zz_data_pos++;
                  actual_block_size++;
                }
                yy_data_pos += dims[2];
              }
              xx_data_pos += dims[1] * dims[2];
            }
            num_flushed_elements += actual_block_size;
          } else if (block_flush_on && flushed_block_id[block_id]) {
            int actual_block_size = 0;
            const T *xx_data_pos = z_data_pos;
            for (int ii = 0; ii < block_size_i; ii++) {
              const T *yy_data_pos = xx_data_pos;
              for (int jj = 0; jj < block_size_j; jj++) {
                const T *zz_data_pos = yy_data_pos;
                for (int kk = 0; kk < block_size_k; kk++) {
                  flushed_block[zz_data_pos - data] = 1;
                  zz_data_pos++;
                  actual_block_size++;
                }
                yy_data_pos += dims[2];
              }
              xx_data_pos += dims[1] * dims[2];
            }
            num_flushed_elements += actual_block_size;
          } else if (block_sift_on && significant_block_id[block_id]) {
            const T *xx_data_pos = z_data_pos;
            for (int ii = 0; ii < block_size_i; ii++) {
              const T *yy_data_pos = xx_data_pos;
              for (int jj = 0; jj < block_size_j; jj++) {
                const T *zz_data_pos = yy_data_pos;
                for (int kk = 0; kk < block_size_k; kk++) {
                  significant_block[zz_data_pos - data] = 1;
                  zz_data_pos++;
                }
                yy_data_pos += dims[2];
              }
              xx_data_pos += dims[1] * dims[2];
            }
          }
          block_id++;
          z_data_pos += block_size_k;
        }
        y_data_pos += block_size_j * dims[2];
      }
      x_data_pos += block_size_i * dims[1] * dims[2];
    }
  }
  return num_flushed_elements;
}

template <typename T>
size_t compute_aux_diff_compress(const SZ::Config conf, const T *data,
                                 const int detection_blocksize,
                                 std::vector<uchar> &significant_block,
                                 std::vector<uchar> &significant_index) {
  // 2D AND 3D ONLY
  if (conf.N == 2) {
    size_t significant_point_counter = 0;
    int block_size = detection_blocksize;
    size_t num_elements = conf.num;
    auto dims = conf.dims;
    // dims[-1] is the fastest changing dimension
    int nx = (int)ceil(1.0 * dims[0] / block_size);
    int ny = (int)ceil(1.0 * dims[1] / block_size);
    // resize the two vectors
    significant_block.resize(nx * ny, 0);
    significant_index.resize(num_elements, 0);
    // compute the root of sum of squares of differences of each block
    // and mark the significant blocks accoring to the threshold from config
    const T *x_data_pos = data;
    for (int i = 0; i < nx; i++) {
      int block_size_i = (i + 1) * block_size > dims[0]
                             ? dims[0] - i * block_size
                             : block_size;
      const T *y_data_pos = x_data_pos;
      for (int j = 0; j < ny; j++) {
        int block_size_j = (j + 1) * block_size > dims[1]
                               ? dims[1] - j * block_size
                               : block_size;
        // compute the first derivative of x direction
        double diff_x = 0;
        double diff_y = 0;
        const T *xx_data_pos = y_data_pos; // block begin index
        for (int ii = 1; ii < block_size_i - 1; ii++) {
          const T *yy_data_pos = xx_data_pos; // block begin index
          for (int jj = 1; jj < block_size_j - 1; jj++) {
            // along y direction
            diff_y += pow(fabs(*(yy_data_pos + 1) - *(yy_data_pos - 1)), 2);
            // along x direction
            diff_x += pow(
                fabs(*(yy_data_pos + dims[1]) - *(yy_data_pos - dims[1])), 2);
            yy_data_pos += 1;
          }
          xx_data_pos += dims[1];
        }
        // check if this block is significantly unsmooth
        // if the central diff is larger than the threshold, then it is
        // significant
        diff_x = sqrt(diff_x) / ((block_size_i - 2) * (block_size_j - 2));
        diff_y = sqrt(diff_y) / ((block_size_i - 2) * (block_size_j - 2));
        if (diff_x > conf.diff_thresh || diff_y > conf.diff_thresh) {
          significant_block[i * ny + j] = 1;
          // mark the significant index
          for (int ii = 0; ii < block_size_i; ii++) {
            for (int jj = 0; jj < block_size_j; jj++) {
              significant_index[(i * block_size_i + ii) * dims[1] +
                                j * block_size_j + jj] = 1;
              significant_point_counter++;
            }
          }
        }
        y_data_pos += block_size_j;
      }
      x_data_pos += block_size_i * dims[1];
    }
    // return the number of significant blocks
    return significant_point_counter;
  }

  if (conf.N == 3) {
    size_t significant_point_counter = 0;
    int block_size = detection_blocksize;
    size_t num_elements = conf.num;
    auto dims = conf.dims;
    // dims[-1] is the fastest changing dimension
    int nx = (int)ceil(1.0 * dims[0] / block_size);
    int ny = (int)ceil(1.0 * dims[1] / block_size);
    int nz = (int)ceil(1.0 * dims[2] / block_size);
    // resize the two vectors
    significant_block.resize(nx * ny * nz, 0);
    significant_index.resize(num_elements, 0);
    // compute the root of sum of squares of differences of each block
    // and mark the significant blocks accoring to the threshold from config
    const T *x_data_pos = data;
    for (int i = 0; i < nx; i++) {
      int block_size_i = (i + 1) * block_size > dims[0]
                             ? dims[0] - i * block_size
                             : block_size;
      const T *y_data_pos = x_data_pos;
      for (int j = 0; j < ny; j++) {
        int block_size_j = (j + 1) * block_size > dims[1]
                               ? dims[1] - j * block_size
                               : block_size;
        const T *z_data_pos = y_data_pos;
        for (int k = 0; k < nz; k++) {
          int block_size_k = (k + 1) * block_size > dims[2]
                                 ? dims[2] - k * block_size
                                 : block_size;
          // compute the first derivative of x direction
          double diff_x = 0;
          double diff_y = 0;
          double diff_z = 0;
          const T *xx_data_pos = z_data_pos; // block begin index
          for (int ii = 1; ii < block_size_i - 1; ii++) {
            const T *yy_data_pos = xx_data_pos; // block begin index
            for (int jj = 1; jj < block_size_j - 1; jj++) {
              const T *zz_data_pos = yy_data_pos; // block begin index
              for (int kk = 1; kk < block_size_k - 1; kk++) {
                // along y direction
                diff_x += pow((*(zz_data_pos + 1) - *(zz_data_pos - 1)) / 2, 2);
                // along x direction
                diff_y += pow(
                    (*(zz_data_pos + dims[2]) - *(zz_data_pos - dims[2])) / 2,
                    2);
                // along z direction
                diff_z += pow((*(zz_data_pos + dims[1] * dims[2]) -
                               *(zz_data_pos - dims[1] * dims[2])) /
                                  2,
                              2);
                zz_data_pos += 1;
              }
              yy_data_pos += dims[2];
            }
            xx_data_pos += dims[1] * dims[2];
          }
          // check if this block is significantly unsmooth
          // if the central diff is larger than the threshold, then it is
          // significant
          int num_sum =
              (block_size_i - 2) * (block_size_j - 2) * (block_size_k - 2);
          if (num_sum > 0) {
            diff_x = sqrt(diff_x / num_sum);
            diff_y = sqrt(diff_y / num_sum);
            diff_z = sqrt(diff_z / num_sum);
          } else {
            diff_x = 0;
            diff_y = 0;
            diff_z = 0;
          }

          if (diff_x > conf.diff_thresh || diff_y > conf.diff_thresh ||
              diff_z > conf.diff_thresh) {
            significant_block[i * ny * nz + j * nz + k] = 1;
            // mark the significant index
            for (int ii = 0; ii < block_size_i; ii++) {
              for (int jj = 0; jj < block_size_j; jj++) {
                for (int kk = 0; kk < block_size_k; kk++) {
                  significant_index[(i * block_size_i + ii) * dims[1] *
                                        dims[2] +
                                    j * block_size_j * dims[2] +
                                    k * block_size_k + kk] = 1;
                  significant_point_counter++;
                }
              }
            }
          }
          z_data_pos += block_size_k;
        }
        y_data_pos += block_size_j * dims[2];
      }
      x_data_pos += block_size_i * dims[1] * dims[2];
    }
    // return the number of significant blocks
    std::cout << "significant_point_counter: " << significant_point_counter
              << std::endl;
    return significant_point_counter;
  } else {
    std::cout << "Error: only support 2D and 3D data" << std::endl;
    exit(0);
  }
};

/*
This function recover the significant_index using significant_block with given
blcoksize
*/
size_t compute_aux_diff_decompress(size_t *global_dims, int N,
                                   const int detection_blocksize,
                                   std::vector<uchar> &significant_block,
                                   std::vector<uchar> &significant_index) {

  size_t num_elements = 1;
  for (int i = 0; i < N; i++) {
    num_elements *= global_dims[i];
  }
  if (N == 2) {
    size_t significant_point_counter = 0;
    int block_size = detection_blocksize;
    auto dims = global_dims;
    // dims[-1] is the fastest changing dimension
    int nx = (int)ceil(1.0 * dims[0] / block_size);
    int ny = (int)ceil(1.0 * dims[1] / block_size);
    // resize the two vectors
    significant_index.resize(num_elements, 0);
    // compute the root of sum of squares of differences of each block
    // and mark the significant blocks accoring to the threshold from config
    size_t block_id = 0;
    size_t index_id = 0;
    for (int i = 0; i < nx; i++) {
      int block_size_i = std::min(block_size, (int)dims[0] - i * block_size);
      for (int j = 0; j < ny; j++) {
        int block_size_j = std::min(block_size, (int)dims[1] - j * block_size);
        if (significant_block[block_id] == 1) {
          // mark the significant index
          index_id = (i * block_size_i) * dims[1] + j * block_size_j;
          for (int ii = 0; ii < block_size_i; ii++) {
            for (int jj = 0; jj < block_size_j; jj++) {
              significant_index[index_id] = 1;
              significant_point_counter++;
              index_id++;
            }
            index_id += dims[1];
          }
        }
        block_id++;
      }
      block_id += ny;
    }
    // return the number of signific·ant blocks
    return significant_point_counter;
  } else if (N == 3) {
    size_t significant_point_counter = 0;
    int block_size = detection_blocksize;
    auto dims = global_dims;
    ;
    // dims[-1] is the fastest changing dimension
    int nx = (int)ceil(1.0 * dims[0] / block_size);
    int ny = (int)ceil(1.0 * dims[1] / block_size);
    int nz = (int)ceil(1.0 * dims[2] / block_size);
    // resize the two vectors
    significant_index.resize(num_elements, 0);
    // compute the root of sum of squares of differences of each block
    // and mark the significant blocks accoring to the threshold from config
    for (int i = 0; i < nx; i++) {
      int block_size_i = (i + 1) * block_size > dims[0]
                             ? dims[0] - i * block_size
                             : block_size;
      for (int j = 0; j < ny; j++) {
        int block_size_j = (j + 1) * block_size > dims[1]
                               ? dims[1] - j * block_size
                               : block_size;
        for (int k = 0; k < nz; k++) {
          int block_size_k = (k + 1) * block_size > dims[2]
                                 ? dims[2] - k * block_size
                                 : block_size;
          if (significant_block[i * ny * nz + j * nz + k] == 1) {
            // mark the significant index
            for (int ii = 0; ii < block_size_i; ii++) {
              for (int jj = 0; jj < block_size_j; jj++) {
                for (int kk = 0; kk < block_size_k; kk++) {
                  significant_index[(i * block_size_i + ii) * dims[1] *
                                        dims[2] +
                                    j * block_size_j * dims[2] +
                                    k * block_size_k + kk] = 1;
                  significant_point_counter++;
                }
              }
            }
          }
        }
      }
    }
    // return the number of significant blocks
    return significant_point_counter;
  } else {
    std::cout << "Error: only support 2D and 3D data" << std::endl;
    return 0;
    exit(0);
  }
};

template <typename T>
std::vector<uchar> compress_opt_reginal_post_process(
    T *data, T *error, const SZ::Config conf, const int detection_blocksize,
    const double noise, std::mt19937 &mt, bool overwrite = true) {
  // 2D AND 3D ONLY
  std::vector<uchar> block_info;
  double block_noise_rng_threshold = conf.block_noise_rng_threshold;

  std::uniform_real_distribution<double> post_noise(-noise, noise);
  // double real_noise = dist(mt);

  if (conf.N == 2) {
    int block_size = detection_blocksize;
    auto dims = conf.dims;
    // dims[-1] is the fastest changing dimension
    int nx = (int)ceil(1.0 * dims[0] / block_size);
    int ny = (int)ceil(1.0 * dims[1] / block_size);
    block_info.resize(nx * ny, 0);
    // compute the variance of error in each block
    const T *x_error_pos, *y_error_pos, *z_error_pos;
    T *x_data_pos, *y_data_pos, *z_data_pos;
    const T *xx_error_pos, *yy_error_pos, *zz_error_pos;
    T *xx_data_pos, *yy_data_pos, *zz_data_pos;
    T local_error_index = 0;
    // int blcok_size_i = 0;
    // int blcok_size_j = 0;
    for (int i = 0; i < nx; i++) {
      int block_size_i = (i + 1) * block_size > dims[0]
                             ? dims[0] - i * block_size
                             : block_size;
      y_error_pos = x_error_pos;
      y_data_pos = x_data_pos;
      for (int j = 0; j < ny; j++) {
        int block_size_j = (j + 1) * block_size > dims[1]
                               ? dims[1] - j * block_size
                               : block_size;
        // compute the value range of the error of this block
        T max_error = 0;
        T min_error = 0;
        xx_error_pos = y_error_pos; // block begin index
        for (int ii = 0; ii < block_size_i; ii++) {
          yy_error_pos = xx_error_pos; // block begin index
          for (int jj = 0; jj < block_size_j; jj++) {
            if (max_error < *yy_error_pos) {
              max_error = *yy_error_pos;
            }
            if (min_error > *yy_error_pos) {
              min_error = *yy_error_pos;
            }
            yy_error_pos++;
          }
          xx_error_pos += dims[1];
        }
        T error_range = max_error - min_error;
        /*
        if the block error range exceeds the threshold,
        then we need to do the post processing on the decompressed data
        otherwise, we do nothing.
        Post-processing: (1) tag this blcok with a soomthing factor
                         (2) smooth the data in this block on demand, overwrite
        flag is used to control this
        */
        if (error_range > block_noise_rng_threshold) {
          // mark the block
          block_info[i * ny + j] = 1;
          // post processing the data in this block on demand
          if (overwrite == true) {
            xx_data_pos = y_data_pos; // block begin index
            for (int ii = 0; ii < block_size_i; ii++) {
              yy_data_pos = xx_data_pos; // block begin index
              for (int jj = 0; jj < block_size_j; jj++) {
                *yy_data_pos = *yy_data_pos + post_noise(mt) * error_range;
                yy_data_pos++;
              }
              xx_data_pos += dims[1];
            }
          }
        }
        y_error_pos += block_size_j;
        y_data_pos += block_size_j;
      }
      x_error_pos += block_size_i * dims[1];
      x_data_pos += block_size_i * dims[1];
    }
  } else if (conf.N == 3) {
    int block_size = detection_blocksize;
    auto dims = conf.dims;
    std::cout << "block_noise_rng_threshold: " << block_noise_rng_threshold
              << std::endl;
    // dims[-1] is the fastest changing dimension
    int nx = (int)ceil(1.0 * dims[0] / block_size);
    int ny = (int)ceil(1.0 * dims[1] / block_size);
    int nz = (int)ceil(1.0 * dims[2] / block_size);
    block_info.resize(nx * ny * nz, 0);
    // compute the variance of error in each block
    T *x_data_pos, *y_data_pos, *z_data_pos;
    const T *x_error_pos, *y_error_pos, *z_error_pos;
    T *xx_data_pos, *yy_data_pos, *zz_data_pos;
    const T *xx_error_pos, *yy_error_pos, *zz_error_pos;
    x_data_pos = data;
    x_error_pos = error;
    T local_error_index = 0;
    for (int i = 0; i < nx; i++) {
      int block_size_i = (i + 1) * block_size > dims[0]
                             ? dims[0] - i * block_size
                             : block_size;
      y_data_pos = x_data_pos;
      y_error_pos = x_error_pos;
      for (int j = 0; j < ny; j++) {
        int block_size_j = (j + 1) * block_size > dims[1]
                               ? dims[1] - j * block_size
                               : block_size;
        z_data_pos = y_data_pos;
        z_error_pos = y_error_pos;
        for (int k = 0; k < nz; k++) {
          int block_size_k = (k + 1) * block_size > dims[2]
                                 ? dims[2] - k * block_size
                                 : block_size;
          // compute the value range of the error of this block
          T max_error = *z_error_pos;
          T min_error = *z_error_pos;
          xx_error_pos = z_error_pos; // block begin index
          for (int ii = 0; ii < block_size_i; ii++) {
            yy_error_pos = xx_error_pos; // block begin index
            for (int jj = 0; jj < block_size_j; jj++) {
              zz_error_pos = yy_error_pos; // block begin index
              for (int kk = 0; kk < block_size_k; kk++) {
                if (max_error < *zz_error_pos) {
                  max_error = *zz_error_pos;
                }
                if (min_error > *zz_error_pos) {
                  min_error = *zz_error_pos;
                }
                zz_error_pos++;
              }
              yy_error_pos += dims[2];
            }
            xx_error_pos += dims[1] * dims[2];
          }
          T error_range = max_error - min_error;
          /*
          if the block error range exceeds the threshold,
          then we need to do the post processing on the decompressed data
          otherwise, we do nothing.
          Post-processing: (1) tag this blcok with a soomthing factor
                           (2) smooth the data in this block on demand,
          overwrite flag is used to control this
          */
          if (error_range > block_noise_rng_threshold) {
            // std::cout << "error_range: " << error_range << std::endl;
            // mark the block
            block_info[i * ny * nz + j * nz + k] = 1;
            // post processing the data in this block on demand
            if (overwrite == true) {
              xx_data_pos = z_data_pos; // block begin index
              for (int ii = 0; ii < block_size_i; ii++) {
                yy_data_pos = xx_data_pos; // block begin index
                for (int jj = 0; jj < block_size_j; jj++) {
                  zz_data_pos = yy_data_pos; // block begin index
                  for (int kk = 0; kk < block_size_k; kk++) {
                    *zz_data_pos = *zz_data_pos + post_noise(mt) * error_range;
                    zz_data_pos++;
                  }
                  yy_data_pos += dims[2];
                }
                xx_data_pos += dims[1] * dims[2];
              }
            }
          }
          z_data_pos += block_size_k;
          z_error_pos += block_size_k;
        }
        y_data_pos += block_size_j * dims[2];
        y_error_pos += block_size_j * dims[2];
      }
      x_data_pos += block_size_i * dims[1] * dims[2];
      x_error_pos += block_size_i * dims[1] * dims[2];
    }

  } else {
    std::cout << "Error: only support 2D and 3D data" << std::endl;
    exit(0);
  }

  return block_info;
}

template <typename T>
void decompress_opt_reginal_post_process(T *data,
                                         std::vector<uchar> &block_info, uint N,
                                         const size_t *dims,
                                         const int detection_blocksize,
                                         const double noise, std::mt19937 mt) {
  // 2D AND 3D ONLY
  // Overwrite the decompressed data
  // This is the overwrite part of compress_opt_reginal_post_process
  // double block_noise_rng_threshold = conf.block_noise_rng_threshold;
  std::uniform_real_distribution<double> post_noise(-noise, noise);
  int block_size = detection_blocksize;
  if (N == 2) {
    // auto dims = conf.dims;
    // dims[-1] is the fastest changing dimension
    int nx = (int)ceil(1.0 * dims[0] / block_size);
    int ny = (int)ceil(1.0 * dims[1] / block_size);
    T *xx_data_pos, *yy_data_pos;
    // add noise to tagged blocks
    for (int i = 0; i < nx; i++) {
      int block_size_i = (i + 1) * block_size > dims[0]
                             ? dims[0] - i * block_size
                             : block_size;
      for (auto j = 0; j < ny; j++) {
        int block_size_j = (j + 1) * block_size > dims[1]
                               ? dims[1] - j * block_size
                               : block_size;
        if (block_info[i * ny + j] == 1) {
          // post processing the data in this block on demand
          xx_data_pos = data + (i * block_size_i) * dims[1] + j * block_size_j;
          for (int ii = 0; ii < block_size_i; ii++) {
            yy_data_pos = xx_data_pos; // block begin index
            for (int jj = 0; jj < block_size_j; jj++) {
              *yy_data_pos = *yy_data_pos + post_noise(mt);
              yy_data_pos++;
            }
            xx_data_pos += dims[1];
          }
        }
      }
    }
  } else if (N == 3) {
    // same procedure with 2d version
    // auto dims = conf.dims;
    // dims[-1] is the fastest changing dimension
    int nx = (int)ceil(1.0 * dims[0] / block_size);
    int ny = (int)ceil(1.0 * dims[1] / block_size);
    int nz = (int)ceil(1.0 * dims[2] / block_size);
    // add noise to tagged blocks
    T *xx_data_pos, *yy_data_pos, *zz_data_pos;
    for (int i = 0; i < nx; i++) {
      int block_size_i = (i + 1) * block_size > dims[0]
                             ? dims[0] - i * block_size
                             : block_size;
      for (auto j = 0; j < ny; j++) {
        int block_size_j = (j + 1) * block_size > dims[1]
                               ? dims[1] - j * block_size
                               : block_size;
        for (auto k = 0; k < nz; k++) {
          int block_size_k = (k + 1) * block_size > dims[2]
                                 ? dims[2] - k * block_size
                                 : block_size;
          if (block_info[i * ny * nz + j * nz + k] == 1) {
            // post processing the data in this block on demand
            xx_data_pos = data + (i * block_size_i) * dims[1] * dims[2] +
                          j * block_size_j * dims[2] + k * block_size_k;
            for (int ii = 0; ii < block_size_i; ii++) {
              yy_data_pos = xx_data_pos; // block begin index
              for (int jj = 0; jj < block_size_j; jj++) {
                zz_data_pos = yy_data_pos; // block begin index
                for (int kk = 0; kk < block_size_k; kk++) {
                  *zz_data_pos = *zz_data_pos + post_noise(mt);
                  zz_data_pos++;
                }
                yy_data_pos += dims[2];
              }
              xx_data_pos += dims[1] * dims[2];
            }
          }
        }
      }
    }
  } else {
    std::cout << "Error: only support 2D and 3D data" << std::endl;
    exit(0);
  }
}

} // namespace SZ
#endif
