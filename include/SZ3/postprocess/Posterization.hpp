#ifndef SZ_POSTERIZATION_HPP
#define SZ_POSTERIZATION_HPP

#include "DisjointSet.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include <iostream>
#include <vector>



namespace SZ {
template <class T> 
class Posterization {

public:
  Posterization(T *data, int dim, int *dims) {
    this->input_data = data;
    this->num_elements = 1;
    this->N = dim;
    this->global_dimensions.resize(N);
    for (int i = 0; i < N; i++) {
      this->num_elements *= dims[i];
      this->global_dimensions[i] = dims[i]; // the last dimension is the fastest changing dimension 
    }
    posterization_dsu = new DisjointSet(num_elements);  }

  ~Posterization() {
    delete posterization_dsu;
  }

  std::vector<int> get_global_dimensions() { return global_dimensions; }

  void perform_posterization(T threshold = 1) {
    posterization_threshold = threshold;
    if (N == 2) {
      Segmentation2D(input_data, threshold);
    } else if (N == 3) {
      Segmentation3D(input_data, threshold);
    }
  }

  std::vector<int> get_segmentation_map(T threshold = 1) {
    posterization_threshold = threshold;
    perform_posterization(threshold);
    return posterization_dsu->get_map();
  }

  std::vector<int> get_segmentation_map(bool test = false) {
    // posterization_dsu->writefiles();
    return posterization_dsu->get_map();
  }

  void set_flush_threshold(T threshold) { posterization_threshold = threshold; }

  int get_second_largest_set_size( T threshold = 1)
  {
    perform_posterization(threshold);
    // posterization_dsu->findSecondLargestTreeSize()
    return posterization_dsu->findSecondLargestTreeSize();
  }

  int get_num_of_components()
  {
    return posterization_dsu->get_num_of_components();
  }

  int get_largest_set_size( T threshold = 1)
  {
    perform_posterization(threshold);
    return posterization_dsu->findLargestTreeSize();
  }
  
  int get_background_area()
  {
    return posterization_dsu->get_background_area();
  }

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


  // void Segmentation2D(T *data, T threshold) {
  //   // python image view
  //   // int h = global_dimensions[1]; // dim1 stride is w 
  //   // int w = global_dimensions[0]; //dim0 stride is 1
  //   int h = global_dimensions[1];
  //   int w = global_dimensions[0]; // the last dimension is the fastest changing dimension 

  //   int 

  //   for (int i = 0; i < h; i++) {
  //     for (int j = 0; j < w; j++) {
  //       if (j < w - 1 &&
  //           std::abs(data[i * w + j] -
  //                    data[posterization_dsu->find(i * w + (j + 1))]) <=
  //               threshold) {
  //         posterization_dsu->union_(i * w + j, i * w + (j + 1));
  //       }
  //       if (i < h - 1 &&
  //           std::abs(data[i * w + j] -
  //                    data[posterization_dsu->find((i + 1) * w + j)]) <=
  //               threshold) {
  //         posterization_dsu->union_(i * w + j, (i + 1) * w + j);
  //       }
  //     }
  //   }

  // }

  void Segmentation2D(T *data, T threshold)
  {

    int dim0 = global_dimensions[0];
    int dim1 = global_dimensions[1]; // fastest changing dimension is the last dimension


    for (int i = 0; i< dim0; i++)
    {
      for (int j = 0; j < dim1; j++)
      {
        if (j < dim1 - 1 &&
            std::abs(data[i * dim1 + j] -
                     data[posterization_dsu->find(i * dim1 + (j + 1))]) <
                threshold) {
          posterization_dsu->union_(i * dim1 + j, i * dim1 + (j + 1));
        }
        if (i < dim0 - 1 &&
            std::abs(data[i * dim1 + j] -
                     data[posterization_dsu->find((i + 1) * dim1 + j)]) <
                threshold) {
          posterization_dsu->union_(i * dim1 + j, (i + 1) * dim1 + j);
        }
      }
    }
  }


  // void Segmentation2D(T *data, T threshold)
  // {

    
  //   int dim0 = global_dimensions[0];
  //   int dim1 = global_dimensions[1]; // fastest changing dimension is the last dimension

  //   for (int i = 0; i < dim0; i++)
  //   {
  //     for (int j = 0; j< dim1 ; j++)
  //     {
  //       if(j+1<dim1)
  //       {
  //         int tmp_j = posterization_dsu->find(i * dim1 + j +1);
  //         T abs_diff = (T) std::abs(data[i * dim1 + j] - data[tmp_j]);
  //         if(abs_diff < threshold)
  //         {
  //           posterization_dsu->union_(i * dim1 + j, tmp_j);
  //         }
  //       }
  //       if(i+1<dim0)
  //       {
  //         int tmp_i = posterization_dsu->find((i+1) * dim1 + j);
  //         T abs_diff = (T) std::abs(data[i * dim1 + j] - data[tmp_i]);
  //         if(abs_diff < threshold)
  //         {
  //           posterization_dsu->union_(i * dim1 + j, tmp_i);
  //         }
  //       }
  //     }
  //   }

  // }

  void Segmentation3D(T *data, T threshold) {
    // int d = global_dimensions[2];
    // int h = global_dimensions[1];
    // int w = global_dimensions[0];

    int d = global_dimensions[0]; 
    int h = global_dimensions[1];
    int w = global_dimensions[2]; // the last dimension is the fastest changing dimension


    SZ::Timer timer;
    timer.start();
    for (int k = 0; k < d; k++) {
      for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
          if (j < w - 1 &&
              std::abs(
                  data[k * h * w + i * w + j] -
                  data[posterization_dsu->find(k * h * w + i * w + (j + 1))]) <
                  threshold) {
            posterization_dsu->union_(k * h * w + i * w + j,
                                      k * h * w + i * w + (j + 1));
          }
          if (i < h - 1 &&
              std::abs(
                  data[k * h * w + i * w + j] -
                  data[posterization_dsu->find(k * h * w + (i + 1) * w + j)]) <
                  threshold) {
            posterization_dsu->union_(k * h * w + i * w + j,
                                      k * h * w + (i + 1) * w + j);
          }
          if (k < d - 1 &&
              std::abs(
                  data[k * h * w + i * w + j] -
                  data[posterization_dsu->find((k + 1) * h * w + i * w + j)]) <
                  threshold) {
            posterization_dsu->union_(k * h * w + i * w + j,
                                      (k + 1) * h * w + i * w + j);
          }
        }
      }
    }
  }
};
}
#endif

