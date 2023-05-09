#include "DisjointSet.hpp"
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <typeinfo>
#include <unordered_map>
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

    // root_label.reserve(num_elements);
    root_size.reserve( num_elements);
  }

  ~Posterization() { delete posterization_dsu; }

  std::vector<int> get_global_dimensions() { return global_dimensions; }

  std::vector<int> get_segmentation_map(T threshold = 1e-5) {
    // std::vector<T> datacopy(input_data, input_data + num_elements);
    posterization_threshold = threshold;
    result_segmentation_map.resize(num_elements);
    // std::copy(this->input_data, this->input_data + num_elements,
    //           result_segmentation_map.begin());
    // std::vector<int> segmentation_map(num_elements);
    if (N == 2) {
      Segmentation2D(input_data, threshold);
    } else if (N == 3) {
      Segmentation3D(input_data, threshold);
    }

    // std::transform(datacopy.begin(), datacopy.end(),
    // segmentation_map.begin(),
    //                [](T i) { return (int)i; });

    return result_segmentation_map;
  }

  void evaluation(const std::vector<T> &segmentation_map) {
    T label_max =
        *std::max_element(segmentation_map.begin(), segmentation_map.end());

    std::cout << "label count =  " << (int)(label_max + 1) << std::endl;
    double lable2size_ratio = (double)(label_max + 1) / (1.0 * num_elements);
    std::cout << "lable / size =  " << lable2size_ratio << std::endl;
  }

  void evaluate(bool remove_background = true) {
    if (remove_background)
      evaluate_no_background();
    else
      evaluate_with_background();
  }

  void evaluate_no_background() {
    int background_count = 0;
    for (const auto &pair : root_size) {
      if (pair.second > background_count) {
        background_count = pair.second;
      }
    }
    // background_count = posterization_dsu->findLargestTreeSize();
    std::cout << "background label count = " << background_count
              << std::endl;
    std::cout << "label count =  " << root_size.size() << std::endl;
    std::cout << "num_elements =  " << num_elements << std::endl;
    double lable2size_ratio =
        (double)(1.0 * (root_size.size() - 1)) /
        (double)(1.0 * num_elements - background_count);
    std::cout << "lable / size =  " << lable2size_ratio << std::endl;
  }

  void evaluate_with_background() {

    std::cout << "label count =  " << root_size.size() << std::endl;
    double lable2size_ratio =
        (double)(root_size.size()) / (1.0 * num_elements);
    std::cout << "lable / size =  " << lable2size_ratio << std::endl;
  }

private:
  T *input_data;
  int num_elements;
  std::vector<int> global_dimensions;
  std::vector<int> result_segmentation_map;
  T posterization_threshold;
  int N;
  DisjointSet *posterization_dsu;
  // std::unordered_map<int, int> root_label;
  // std::unordered_map<int, std::array<int,2>> root_label;
  std::unordered_map<int, int> root_size; // root label and size of the tree

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

    posterization_dsu = new DisjointSet(num_elements);
    for (int i = 0; i < h; i++) {
      for (int j = 0; j < w; j++) {
        if (j < w - 1 &&
            std::abs(data[i * w + j] - data[i * w + (j + 1)]) <= threshold) {
          posterization_dsu->union_(i * w + j, i * w + (j + 1));
        }
        if (i < h - 1 &&
            std::abs(data[i * w + j] - data[(i + 1) * w + j]) <= threshold) {
          posterization_dsu->union_(i * w + j, (i + 1) * w + j);
        }
      }
    }

       for (int i=0; i< num_elements; i++)
    {
      int root = posterization_dsu->find(i);
        auto it = root_size.find(root);
        if (it == root_size.end()) {
          root_size.insert({root, 1});
          // current_label += 1;
          // root_size.insert({root, 1});
        }
        result_segmentation_map[i] = root;
        root_size[root] += 1;
        // root_size[root] += 1;
    }


  // two has table version, not efficient 
    // int current_label = 0;
    // for (int i=0; i< num_elements; i++)
    // {
    //   int root = posterization_dsu->find(i);
    //     auto it = root_label.find(root);
    //     if (it == root_label.end()) {
    //       root_label.insert({root, {current_label, 1}});
    //       current_label += 1;
    //       // root_size.insert({root, 1});
    //     }
    //     data[i] = root_label[root][0];
    //     root_label[root][1] += 1;
    //     // root_size[root] += 1;
    // }
  }

  void Segmentation3D(T *data, T threshold) {
    int d = global_dimensions[2];
    int h = global_dimensions[1];
    int w = global_dimensions[0];
    posterization_dsu = new DisjointSet(num_elements);
    for (int k = 0; k < d; k++) {
      for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
          if (j < w - 1 &&
              std::abs(data[k * h * w + i * w + j] -
                       data[k * h * w + i * w + (j + 1)]) <= threshold) {
            posterization_dsu->union_(k * h * w + i * w + j,
                                      k * h * w + i * w + (j + 1));
          }
          if (i < h - 1 &&
              std::abs(data[k * h * w + i * w + j] -
                       data[k * h * w + (i + 1) * w + j]) <= threshold) {
            posterization_dsu->union_(k * h * w + i * w + j,
                                      k * h * w + (i + 1) * w + j);
          }
          if (k < d - 1 &&
              std::abs(data[k * h * w + i * w + j] -
                       data[(k + 1) * h * w + i * w + j]) <= threshold) {
            posterization_dsu->union_(k * h * w + i * w + j,
                                      (k + 1) * h * w + i * w + j);
          }
        }
      }
    }
    int current_label = 0;
    // for (int i=0; i< num_elements; i++)
    // {
    //   int root = posterization_dsu->find(i);
    //     auto it = root_label.find(root);
    //     if (it == root_label.end()) {
    //       root_label.insert({root, current_label});
    //       current_label += 1;
    //       // root_size.insert({root, 1});
    //     }
    //     data[i] = root_label[root];
    //     // root_size[root] += 1;
    // }

           for (int i=0; i< num_elements; i++)
    {
      int root = posterization_dsu->find(i);
        auto it = root_size.find(root);
        if (it == root_size.end()) {
          root_size.insert({root, 1});
          // current_label += 1;
          // root_size.insert({root, 1});
        }
        result_segmentation_map[i] = root;
        root_size[root] += 1;
        // root_size[root] += 1;
    }

    //     for (int i=0; i< num_elements; i++)
    // {
    //   int root = posterization_dsu->find(i);
    //     auto it = root_label.find(root);
    //     if (it == root_label.end()) {
    //       root_label.insert({root, {current_label, 1}});
    //       current_label += 1;
    //       // root_size.insert({root, 1});
    //     }
    //     data[i] = root_label[root][0];
    //     root_label[root][1] += 1;
    //     // root_size[root] += 1;
    // }
  }
};
