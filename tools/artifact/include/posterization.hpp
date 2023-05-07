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
  }

  ~Posterization() {}

  std::vector<int> get_global_dimensions() { return global_dimensions; }

  std::vector<T> get_segmentation_map(T threshold=1e-5) {
    // std::vector<T> datacopy(input_data, input_data + num_elements);
    posterization_threshold = threshold;
    result_segmentation_map.resize(num_elements);
    std::copy(this->input_data, this->input_data + num_elements, result_segmentation_map.begin());
    // std::vector<int> segmentation_map(num_elements);
    if (N == 2) {
      Segmentation2D(result_segmentation_map.data(), threshold);
    } else if (N == 3) {
      Segmentation3D(result_segmentation_map.data(), threshold);
    }

    // std::transform(datacopy.begin(), datacopy.end(),
    // segmentation_map.begin(),
    //                [](T i) { return (int)i; });

    return result_segmentation_map;
  }

  void evaluation(const std::vector<T> &segmentation_map )
  {
    T label_max = *std::max_element(segmentation_map.begin(), segmentation_map.end());

    std::cout<<"label count =  "<< (int)(label_max+1) <<std::endl;
    double lable2size_ratio = (double)(label_max+1) /(1.0*num_elements) ;
    std::cout<<"lable / size =  "<< lable2size_ratio <<std::endl;

  }

  void evaluate()
  {
    evaluate(result_segmentation_map,posterization_threshold);
  }

  void evaluate(std::vector<T> &segmantation_map,T &threshold)
  {
    // find out the most frequent label, treat it as the background
    std::unordered_map<T, int> label_count;
    for (int i = 0; i < num_elements; i++)
    {
      label_count[segmantation_map[i]]++;
    }
    T background_label = 0;
    int max_count = 0;
    for (auto it = label_count.begin(); it != label_count.end(); it++)
    {
      if (it->second > max_count)
      {
        max_count = it->second;
        background_label = it->first;
      }
    }
    std::cout << "background label = " << background_label << std::endl;
    std::cout << "background label count = " << max_count << std::endl;
    // remove the background label

    std::cout<<"label count =  "<< label_count.size()-1 <<std::endl;
    double lable2size_ratio = (double)(1.0*(label_count.size()-1))/(double)(1.0*num_elements) ;
    std::cout<<"lable / size =  "<< lable2size_ratio <<std::endl;
  }

private:
  T *input_data;
  int num_elements;
  std::vector<int> global_dimensions;
  std::vector<T> result_segmentation_map;
  T posterization_threshold;
  int N;

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

    auto dsu = DisjointSet(num_elements);
    for (int i = 0; i < h; i++) {
      for (int j = 0; j < w; j++) {
        if (j < w - 1 &&
            std::abs(data[i * w + j] - data[i * w + (j + 1)]) <= threshold) {
          dsu.union_(i * w + j, i * w + (j + 1));
        }
        if (i < h - 1 &&
            std::abs(data[i * w + j] - data[(i + 1) * w + j]) <= threshold) {
          dsu.union_(i * w + j, (i + 1) * w + j);
        }
      }
    }
    std::unordered_map<int, int> labels;
    int current_label = 0;
    for (int i = 0; i < h; i++) {
      for (int j = 0; j < w; j++) {
        int root = dsu.find(i * w + j);
        auto it = labels.find(root);
        if (it == labels.end()) {
          labels.insert({root, current_label});
          current_label += 1;
        }
        data[i * w + j] = labels[root];
      }
    }
  }

  void Segmentation3D(T *data, T threshold) {
    int d = global_dimensions[2];
    int h = global_dimensions[1];
    int w = global_dimensions[0];
    auto dsu = DisjointSet(num_elements);
    for (int k = 0; k < d; k++) {
      for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
          if (j < w - 1 &&
              std::abs(data[k * h * w + i * w + j] -
                       data[k * h * w + i * w + (j + 1)]) <= threshold) {
            dsu.union_(k * h * w + i * w + j, k * h * w + i * w + (j + 1));
          }
          if (i < h - 1 &&
              std::abs(data[k * h * w + i * w + j] -
                       data[k * h * w + (i + 1) * w + j]) <= threshold) {
            dsu.union_(k * h * w + i * w + j, k * h * w + (i + 1) * w + j);
          }
          if (k < d - 1 &&
              std::abs(data[k * h * w + i * w + j] -
                       data[(k + 1) * h * w + i * w + j]) <= threshold) {
            dsu.union_(k * h * w + i * w + j, (k + 1) * h * w + i * w + j);
          }
        }
      }
    }
    std::unordered_map<int, int> labels;
    int current_label = 0;
    for (int k = 0; k < d; k++) {
      for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
          int root = dsu.find(k * h * w + i * w + j);
          auto it = labels.find(root);
          if (it == labels.end()) {
            labels.insert({root, current_label});
            current_label += 1;
          }
          data[k * h * w + i * w + j] = labels[root];
        }
      }
    }
  }


  //
};
