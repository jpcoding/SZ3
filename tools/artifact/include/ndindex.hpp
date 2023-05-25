#include<vector>
#include<array>

template <class T> 
class NDIndex {
    public:
    // For global view use 
    NDIndex(T*inputdata, int dim, int *shape) {
        size = 1;
        for (int i = 0; i < dim; i++) {
            size *= shape[i];
            this->shape.push_back(shape[i]);
        }
        data = inputdata;

    }

    // operations
    NDIndex<T> operator+(NDIndex<T> &other)
    {
        std::vector<T> result_data(size);
        NDIndex<T> result(result_data.data(),N, shape.data());
        for (int i = 0; i < size; i++)
        {
            result.data[i] = data[i] + other.data[i];
        }
        return result;
    }

    NDIndex<T> operator-(NDIndex<T> &other)
    {
        std::vector<T> result_data(size);
        NDIndex<T> result(result_data.data(),N, shape.data());
        for (int i = 0; i < size; i++)
        {
            result.data[i] = data[i] - other.data[i];
        }
        return result;
    }

    NDIndex<T> operator*(NDIndex<T> &other)
    {
        std::vector<T> result_data(size);
        NDIndex<T> result(result_data.data(),N, shape.data());
        for (int i = 0; i < size; i++)
        {
            result.data[i] = data[i] * other.data[i];
        }
        return result;
    }

    NDIndex<T> operator/(NDIndex<T> &other)
    {
        std::vector<T> result_data(size);
        NDIndex<T> result(result_data.data(),N, shape.data());
        for (int i = 0; i < size; i++)
        {
            result.data[i] = data[i] / other.data[i];
        }
        return result;
    }

    NDIndex<T> operator()(int startx, int starty, int endx, int endy, int stridex=1, int stridey=1)
    {
        std::vector<T> result_data(size);
        NDIndex<T> result(result_data.data(),N, shape.data());
        for (int i = 0; i < size; i++)
        {
            result.data[i] = data[i] / other.data[i];
        }
        return result;
    }








    private:
    int N;
    size_t size;
    std::vector<int> shape;
    T *data;
};
