#ifndef SZ_DSJ_HPP
#define SZ_DSJ_HPP


#include <SZ3/utils/FileUtil.hpp>
#include <algorithm>
#include <cstddef>
#include <unordered_map>
#include <vector>

class DisjointSet {
public:
  DisjointSet(int num_elements) {
    this->num_elements = num_elements;
    this->parent.resize(num_elements);
    this->rank.resize(num_elements);
    make_set();
  }

  ~DisjointSet() {}

  int find(int x) {
    if (parent[x] != x) {
      parent[x] = find(parent[x]);
    }
    return parent[x];
  }

  void union_(int x, int y)
  {
    int xroot = find(x);
    int yroot = find(y);
    if (xroot != yroot)
    {
      if (rank[xroot]> rank[yroot])
      {
        parent[yroot] = xroot;
      }
      else
      {
        parent[xroot] = yroot; // path compression is used here. 
        if (rank[xroot] == rank[yroot])
        {
          rank[yroot] += 1;
        }
      }
    }
  }

  void writefiles() {
    SZ::writefile("parent.dat", parent.data(), parent.size());
    SZ::writefile("rank.dat", rank.data(), rank.size());
  }

  std::vector<int> get_map() { return parent; }

  int findLargestTreeSize() {
    int largestTreeSize = 0;
    std::unordered_map<int, int> treeSizes;

    // Iterate through all elements
    for (int i = 0; i < num_elements; i++) {
      int root = this->find(i);
      treeSizes[root]++; // Increment size of tree

      // Update largest tree size if necessary
      if (treeSizes[root] > largestTreeSize) {
        largestTreeSize = treeSizes[root];
      }
    }

    return largestTreeSize;
  }

  void set_num_of_components( int n) {
    this-> num_of_components = n;
  }

  int get_num_of_components() {
    return this->num_of_components;
  }

  void set_background_area( int n) {
    this->background_size = n;
  }

  int get_background_area() {
    return this->background_size;
  }

  size_t findSecondLargestTreeSize() {
      size_t largestTreeSize = 0;
      size_t secondLargestTreeSize = 0;
      std::unordered_map<int, int> treeSizes;
      std::vector<int> counts; 
      // Iterate through all elements
      for (int i = 0; i < num_elements; i++) {
        int root = this->find(i);
        treeSizes[root]++; // Increment size of tree
        // // Update largest and second largest tree size if necessary
        // if (treeSizes[root] > largestTreeSize) {
        //   secondLargestTreeSize = largestTreeSize;
        //   largestTreeSize = treeSizes[root];
        // } else if (treeSizes[root] != largestTreeSize && treeSizes[root] > secondLargestTreeSize) {
        //   secondLargestTreeSize = treeSizes[root];
        // }
      }
      counts.reserve(treeSizes.size());
      for (auto it = treeSizes.begin(); it != treeSizes.end(); it++) {
        counts.push_back(it->second);
      }
      std::sort(counts.begin(), counts.end(), std::greater<int>());
      if (counts.size() > 1) {
        secondLargestTreeSize = counts[1];
      }
      set_num_of_components(treeSizes.size());
      set_background_area(counts[0]);

      return secondLargestTreeSize;
    }

private:
  int num_elements;
  std::vector<int> parent;
  std::vector<int> rank;
  int background_size; 
  size_t label_num;
  int num_of_components; 

  void make_set() {
    for (int i = 0; i < this->num_elements; i++) {
      parent[i] = i;
      rank[i] = 0;
    }
  }
};



#endif
