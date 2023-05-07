#include <iostream>
#include <vector>
#include <array>
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "critical_points.hpp"
#include "interpolation_walk.hpp"
#include "posterization.hpp"




int main (int argc, char **argv)
{
    if (argc < 2)
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

    // float *data = new float[num_elements];
    std::vector<float> ddata(num_elements);
    SZ::readfile<float>(argv[1], num_elements, ddata.data());
    // constrcut the critical point map;
    CriticalPointsCalculator cp_calculator(ddata.data(), N, global_dimensions.data());
    std::vector<int> critical_points_map = cp_calculator.get_critical_points_map();

    SZ::writefile<int>("critical_points_map.dat",  critical_points_map.data(),num_elements);
    std::vector<int> cp_index;
    int match_count = 0;
    if(N==2){
    std::vector<int> cp_xpaddings;
    std::vector<int> cp_ypaddings;
    cp_calculator.pattern_match_global(cp_index, cp_xpaddings, cp_ypaddings,match_count);
    SZ::writefile<int>("cp_index.dat",  cp_index.data(),cp_index.size());
    SZ::writefile<int>("cp_xpaddings.dat",  cp_xpaddings.data(),cp_xpaddings.size());
    SZ::writefile<int>("cp_ypaddings.dat",  cp_ypaddings.data(),cp_ypaddings.size());
    }
    else if(N==3){
    std::vector<int> cp_xpaddings;
    std::vector<int> cp_ypaddings;
    std::vector<int> cp_zpaddings;
    cp_calculator.pattern_match_global(cp_index, cp_xpaddings, cp_ypaddings, cp_zpaddings,match_count);
    SZ::writefile<int>("cp_index.dat",  cp_index.data(),cp_index.size());
    SZ::writefile<int>("cp_xpaddings.dat",  cp_xpaddings.data(),cp_xpaddings.size());
    SZ::writefile<int>("cp_ypaddings.dat",  cp_ypaddings.data(),cp_ypaddings.size());
    SZ::writefile<int>("cp_zpaddings.dat",  cp_zpaddings.data(),cp_zpaddings.size());
    }
    else{
        std::cout<<"N = "<<N<<" is not supported!"<<std::endl;
        return 0;
    }


    std::cout<<"match_count = "<<match_count<<std::endl;

    return 0;

}