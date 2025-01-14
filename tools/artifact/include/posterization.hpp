#include "DisjointSet.hpp"
#include "SZ3/utils/Timer.hpp"
#include <iostream>
#include <vector>

template <class T> class Posterization {

public:
  Posterization(T *data, int dim, int *dims) {
    this->input_data = data;
    this->num_elements = 1;
    this->N = dim;
    this->global_dimensions.resize(N);
    for (int i = 0; i < N; i++) {
      this->num_elements *= dims[i];
      this->global_dimensions[i] = dims[i];
    }

    posterization_dsu = new DisjointSet(num_elements);  }

  ~Posterization() {
    delete posterization_dsu;
  }

  std::vector<int> get_global_dimensions() { return global_dimensions; }

  std::vector<int> get_segmentation_map(T threshold = 1e-5) {
    posterization_threshold = threshold;

    if (N == 2) {
      Segmentation2D(input_data, threshold);
    } else if (N == 3) {
      Segmentation3D(input_data, threshold);
    }

    return posterization_dsu->get_map();
  }

  void set_flush_threshold(T threshold) { posterization_threshold = threshold; }

  void evaluate() {
    int background_size = 0;
    int label_num = 0;
    std::vector<int> label_set(num_elements);
    std::vector<int> result_segmentation_map = posterization_dsu->get_map();
    for (int i = 0; i < num_elements; i++) {
      int root = result_segmentation_map[i];
      label_set[root] += 1;
    }

    for (int i = 0; i < num_elements; i++) {
      if (label_set[i] > 0) {
        label_num++;
      }
      if (label_set[i] > background_size) {
        background_size = label_set[i];
      }
    }
    std::cout << "background size = " << background_size << std::endl;
    std::cout << "label count =  " << label_num << std::endl;
    std::cout << "dataset size =  " << num_elements << std::endl;
    double lable2size_ratio = (double)(1.0 * (label_num - 1)) /
                              (double)(1.0 * num_elements - background_size);
    std::cout << "without background: lable / size =  " << lable2size_ratio
              << std::endl;
    std::cout << "with background: lable / size =  "
              << (double)((label_num) / (1.0 * num_elements)) << std::endl;
  }

private:
  T *input_data;
  int num_elements;
  std::vector<int> global_dimensions;

  T posterization_threshold;
  int N;
  DisjointSet *posterization_dsu;
  // std::unordered_map<int, int> root_size; // root label and size of the tree
  // std::map<int, int> root_size; // root label and size of the tree

  /*
  def segmentation(data, threshold):
    h, w = data.shape[:2]
    num_pixels = h * w
    dsu = DisjointSet(num_pixels)

    for i in range(h):
        for j in range(w):
            if j < w - 1 and abs(data[i, j] - data[i, j+1]) <= threshold:
                dsu.union(i*w+j, i*w+(j+1))
            if i < h - 1 and abs(data[i, j] - data[i+1, j]) <= threshold:
                dsu.union(i*w+j, (i+1)*w+j)

    labels = {}
    current_label = 0
    for i in range(h):
        for j in range(w):
            root = dsu.find(i*w+j)
            if root not in labels:
                labels[root] = current_label
                current_label += 1
            data[i, j] = labels[root]

    return data
  */
  void Segmentation2D(T *data, T threshold) {
    // python image view
    int h = global_dimensions[1];
    int w = global_dimensions[0];

    for (int i = 0; i < h; i++) {
      for (int j = 0; j < w; j++) {
        if (j < w - 1 &&
            std::abs(data[i * w + j] -
                     data[posterization_dsu->find(i * w + (j + 1))]) <=
                threshold) {
          posterization_dsu->union_(i * w + j, i * w + (j + 1));
        }
        if (i < h - 1 &&
            std::abs(data[i * w + j] -
                     data[posterization_dsu->find((i + 1) * w + j)]) <=
                threshold) {
          posterization_dsu->union_(i * w + j, (i + 1) * w + j);
        }
      }
    }

  }

  void Segmentation3D(T *data, T threshold) {
    int d = global_dimensions[2];
    int h = global_dimensions[1];
    int w = global_dimensions[0];
    SZ::Timer timer;
    timer.start();
    for (int k = 0; k < d; k++) {
      for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
          if (j < w - 1 &&
              std::abs(
                  data[k * h * w + i * w + j] -
                  data[posterization_dsu->find(k * h * w + i * w + (j + 1))]) <=
                  threshold) {
            posterization_dsu->union_(k * h * w + i * w + j,
                                      k * h * w + i * w + (j + 1));
          }
          if (i < h - 1 &&
              std::abs(
                  data[k * h * w + i * w + j] -
                  data[posterization_dsu->find(k * h * w + (i + 1) * w + j)]) <=
                  threshold) {
            posterization_dsu->union_(k * h * w + i * w + j,
                                      k * h * w + (i + 1) * w + j);
          }
          if (k < d - 1 &&
              std::abs(
                  data[k * h * w + i * w + j] -
                  data[posterization_dsu->find((k + 1) * h * w + i * w + j)]) <=
                  threshold) {
            posterization_dsu->union_(k * h * w + i * w + j,
                                      (k + 1) * h * w + i * w + j);
          }
        }
      }
    }

    timer.stop("nested loop for segmentation");

  }
};
