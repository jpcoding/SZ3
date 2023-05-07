#include<vector>
#include<array>

template <class T, int N> 
class NDView {
    public:
    // For global view use 
    NDView(T*inputdata, int *shape) {
        num_elements = 1;
        for (int i = 0; i < N; i++) {
            num_elements *= shape[i];
            global_shape.push_back(shape[i]);
            local_shape.push_back(shape[i]);
        }
        data = inputdata;

    }
    // For local view use 
    NDView(T*inputdata, int *shape, int startx, int starty, int endx, int endy, int stridex=1, int stridey=1) {
        num_elements = 1;
        local_strides.resize(2);
        local_strides[0] = stridex;
        local_strides[1] = stridey;
        for (int i = 0; i < N; i++) {
            num_elements *= shape[i];
            global_shape.push_back(shape[i]);
        }
        data = inputdata;

    }

    // operations
    NDView<T,N> operator+(NDView<T,N> &other)
    {
        std::vector<T> result_data(num_elements);
        NDView<T,N> result(result_data.data(), shape);
        for (int i = 0; i < num_elements; i++)
        {
            result.data[i] = data[i] + other.data[i];
        }
        return result;
    }

    NDView<T,N> operator-(NDView<T,N> &other)
    {
        std::vector<T> result_data(num_elements);
        NDView<T,N> result(result_data.data(), shape);
        for (int i = 0; i < num_elements; i++)
        {
            result.data[i] = data[i] - other.data[i];
        }
        return result;
    }

    NDView<T,N> operator*(NDView<T,N> &other)
    {
        std::vector<T> result_data(num_elements);
        NDView<T,N> result(result_data.data(), shape);
        for (int i = 0; i < num_elements; i++)
        {
            result.data[i] = data[i] + other.data[i];
        }
        return result;
    }




    private:
    int num_elements;
    int *shape;
    T *data;
    std::vector<int> global_shape;
    // local view
    std::vector<int> local_shape;
    int local_offset;
    std::vector<int> local_strides;


};
