#include <bits/types/cookie_io_functions_t.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstddef>
#include <array>
#include "SZ3/api/sz.hpp"
#include "critical_points.hpp"
#include "interpolation_walk.hpp"




int main (int argc, char **argv)
{
    if (argc != 2)
    {
        std::cout << "Usage: ./interp_detect <input_file> <N> <dim1 dim2 ...> " << std::endl;
        return 0;
    }
    int N = atoi(argv[2]);
    int num_elements = 1;
    std::vector<int> global_dimensions;
    for (int i = 0; i < N; i++)
    {
        global_dimensions.push_back(atoi(argv[3+i]));
        num_elements *= global_dimensions[i];
    }

    float *data = new float[num_elements];
    SZ::readfile<float>(argv[1], num_elements, data);
    // constrcut the critical point map;
    CriticalPointsCalculator cp_calculator(data, N, global_dimensions.data());
    std::vector<int> critical_points_map = cp_calculator.get_critical_points_map();


}