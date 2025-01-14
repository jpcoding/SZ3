
#include <SZ3/utils/FileUtil.hpp>
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
        parent[xroot] = yroot;
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

private:
  int num_elements;
  std::vector<int> parent;
  std::vector<int> rank;

  void make_set() {
    for (int i = 0; i < this->num_elements; i++) {
      parent[i] = i;
      rank[i] = 0;
    }
  }
};
