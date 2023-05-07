#include<vector>

class DisjointSet
{
    public:
    DisjointSet(int num_elements)
    {
        this->num_elements = num_elements;
        this->parent.resize(num_elements);
        this->rank.resize(num_elements);
        make_set();
    }

    ~DisjointSet()
    {
        
    }

    int find(int x)
    {
        if (parent[x] != x)
        {
            parent[x] = find(parent[x]);
        }
        return parent[x];
    }

    void union_ (int x, int y)
    {
        int xroot = find(x);
        int yroot = find(y);

        if (xroot == yroot)
        {
            return;
        }
        if (this->rank[xroot] < this->rank[yroot])
        {
            parent[xroot] = yroot;
        }
        else if (this->rank[xroot] > this->rank[yroot])
        {
            parent[yroot] = xroot;
        }
        else
        {
            parent[yroot] = xroot;
            this->rank[xroot] += 1;
        }
    }


    private: 
    int num_elements;
    std::vector<int> parent;
    std::vector<int> rank;

    void make_set()
    {
        for (int i = 0; i < this->num_elements; i++)
        {
            parent[i] = i;
            rank[i] = 0;
        }
    }
};

