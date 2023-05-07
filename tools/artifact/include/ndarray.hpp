#include<vector>
#include<iostream>

template <class T>
class NDArray
{
    public:
    NDArray(int N, int *global_dimensions)
    {
        this->N = N;
        this->num_elements =1;
        this->global_dimensions.resize(N);
        for (int i = 0; i < N; i++)
        {
            this->global_dimensions[i] = global_dimensions[i];
            this->num_elements *= global_dimensions[i];
        }
        this->data = new T[num_elements];
    }
    NDArray(int N, int *global_dimensions, T *data)
    {
        this->N = N;
        this->num_elements =1;
        this->global_dimensions.resize(N);
        for (int i = 0; i < N; i++)
        {
            this->global_dimensions[i] = global_dimensions[i];
            this->num_elements *= global_dimensions[i];
        }
        this->data = data;
    }

    // Getters and Setters
    T get_value(int global_index)
    {
        return data[global_index];
    }
    T* operator()(int x)
    {
        return data+(x);
    }

    T* operator()(int x, int y)
    {
        return (data+x*global_dimensions[1] + y);
    }

    T* operator()(int x, int y, int z)
    {
        return data+(x*global_dimensions[1]*global_dimensions[2] + y*global_dimensions[2] + z);
    }
    
    // Setters
    void set_value(int global_index, T value)
    {
        data[global_index] = value;
    }

    std::vector<int> get_global_dimensions();
    size_t get_num_elements();
    T* get_data();
    // Deconstructor
    ~NDArray()
    {
        delete[] data;
    }

    // Slice functions
    NDArray<T> slice(int startx, int endx, int starty, int endy)
    {
        int new_global_dimensions[2] = {endx-startx, endy-starty};
        NDArray<T> new_array(2, new_global_dimensions);
        for(int i = startx; i < endx; i++)
        {
            for(int j = starty; j < endy; j++)
            {
                int new_global_index[2] = {i-startx, j-starty};
                int global_index[2] = {i, j};
                new_array.set_value(new_global_index, get_value(global_index));
            }
        }
        return new_array;
    }



    NDArray<T> slice_s(int startx, int endx, int stridex, int starty, int endy, int stridey);

    NDArray<T> slice(int startx, int endx, int starty, int endy, int startz, int endz);
    NDArray<T> slice(int startx, int endx, int starty, int endy, int startz, int endz, int startw, int endw);


    // basic arthmetic operations
    NDArray<T> operator+(NDArray<T> &other);
    NDArray<T> operator-(NDArray<T> &other);
    NDArray<T> operator*(NDArray<T> &other);
    NDArray<T> operator/(NDArray<T> &other);



    private:
    std::vector<int> global_dimensions;
    size_t num_elements;
    T *data;
    int N;
};



