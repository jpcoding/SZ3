#include <cstddef>
#include<iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <array>
#include "SZ3/api/sz.hpp"
#include "SZ3/api/sz_config.hpp"


template <class T>
class CriticalPointsCalculator
{

    public:
    CriticalPointsCalculator(T *data, int N, int *global_dimensions)
    {
        this->global_dimensions = global_dimensions;
        this->data = data;
        this->num_elements =1;
        this->global_dimensions.resize(N);
        for (int i = 0; i < N; i++)
        {
            this->global_dimensions[i] = global_dimensions[i];
            this->num_elements *= global_dimensions[i];
        }
    }

    std::vector<int>  get_critical_points_map()
    {
        std::vector<int> critical_points_map(num_elements, 0);
        for (int i = 0; i < num_elements; i++)
        {
            critical_points_map[i] = calculate_single_point(i);
        }
        return critical_points_map;
    }


    private:

    std::vector<int> global_dimensions;
    size_t num_elements;
    T *data;
    int N;


    size_t get_global_dimensions(std::array<size_t, N> global_dimensions)
    {
        this->global_dimensions = global_dimensions;
    }

    int calculate_single_point (int global_index)
    {
        if(N==2)
        {
            int idx = global_index % global_dimensions[0];
            int idy = global_index / global_dimensions[0];
            if(idx==0 || idy==0 ||idx==global_dimensions[0]-1 || idy==global_dimensions[1]-1)
            {
                return 0;
            }
            else
            {
                int sign1 = (data[global_index]> get_value(idx-1, idy)) - (data[global_index]< get_value(idx-1, idy));
                int sign2 = (data[global_index]> get_value(idx+1, idy))- (data[global_index]< get_value(idx+1, idy));
                int sign3 = (data[global_index]> get_value(idx, idy-1))- (data[global_index]< get_value(idx, idy-1));
                int sign4 = (data[global_index]> get_value(idx, idy+1))- (data[global_index]< get_value(idx, idy+1));
                return std::abs(sign1+sign2+sign3+sign4);
           }
        }
        else if(N==3)
        {
            int idx = global_index % global_dimensions[0];
            int idy = (global_index / global_dimensions[0]) % global_dimensions[1];
            int idz = global_index / (global_dimensions[0] * global_dimensions[1]);
            if(idx==0 || idy==0 || idz==0 || idx==global_dimensions[0]-1 
                || idy==global_dimensions[1]-1 || idz==global_dimensions[2]-1)
            {
                return 0;
            }
            else
            {
                int sign1 = (data[global_index]> get_value(idx-1, idy, idz)) - (data[global_index]< get_value(idx-1, idy, idz));
                int sign2 = (data[global_index]> get_value(idx+1, idy, idz))- (data[global_index]< get_value(idx+1, idy, idz));
                int sign3 = (data[global_index]> get_value(idx, idy-1, idz))- (data[global_index]< get_value(idx, idy-1, idz));
                int sign4 = (data[global_index]> get_value(idx, idy+1, idz))- (data[global_index]< get_value(idx, idy+1, idz));
                int sign5 = (data[global_index]> get_value(idx, idy, idz-1))- (data[global_index]< get_value(idx, idy, idz-1));
                int sign6 = (data[global_index]> get_value(idx, idy, idz+1))- (data[global_index]< get_value(idx, idy, idz+1));
                return std::abs(sign1+sign2+sign3+sign4+sign5+sign6);
            }
        }
        else
        {
            std::cout << "CriticalPointsCalculator: N is not 2 or 3" << std::endl;
            return 0;
        }
    }

    T get_value (int idx, int idy)
    {
        return data[idx + idy * global_dimensions[0]];

    }

    T get_value (int idx, int idy, int idz)
    {
        return data[idx + idy * global_dimensions[0] + idz * global_dimensions[0] * global_dimensions[1]];

    }


};

