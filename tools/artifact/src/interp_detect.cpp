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
    int N = atoi(argv[3]);
    int num_elements = 1;
    std::vector<int> global_dimensions;
    for (int i = 0; i < N; i++)
    {
        global_dimensions.push_back(atoi(argv[4+i]));
        num_elements *= global_dimensions[i];
    }

    // float *data = new float[num_elements];
    std::vector<float> odata(num_elements);
    std::vector<float> ddata(num_elements);
    SZ::readfile<float>(argv[1], num_elements, odata.data());
    SZ::readfile<float>(argv[2], num_elements, ddata.data());
    // constrcut the critical point map;
    CriticalPointsCalculator cp_calculator(ddata.data(), N, global_dimensions.data());
    std::vector<int> critical_points_map = cp_calculator.get_critical_points_map();

    SZ::writefile<int>("critical_points_map.dat",  critical_points_map.data(),num_elements);

    // construct the interpolation walk
    InterpolationWalk<float> interpolation_walk(ddata.data(), N, global_dimensions.data());

    InterpolationWalk<float> interpolation_walk_original(odata.data(), N, global_dimensions.data());

    std::vector<float> interp_errors_orig(num_elements,0);
    std::vector<float> interp_errors_decompressed(num_elements,0);
    std::vector<int> cp_indeice;

    for (int i =0; i < num_elements; i++)
    {
        if (critical_points_map[i] == 4)
        {
            interp_errors_decompressed[i]=(interpolation_walk.interp_walk(i));
            interp_errors_orig[i] = (interpolation_walk_original.interp_walk(i));
            cp_indeice.push_back(i);

        }
    }
    SZ::writefile("oerror.dat", interp_errors_orig.data(), interp_errors_orig.size());
    SZ::writefile("derror.dat", interp_errors_decompressed.data(), interp_errors_decompressed.size());
    SZ::writefile("cp_indeice.dat", cp_indeice.data(), cp_indeice.size());

    // construct the posterization
    Posterization<float> posterization(ddata.data(), N, global_dimensions.data());
    SZ::Timer timer;
    timer.start();
    auto segmentation_map =posterization.get_segmentation_map( 0.0001);
    
    timer.stop("posterization");
timer.start();
    posterization.evaluate();
    timer.stop("posterization evaluation");
    
    // posterization.evaluation( segmentation_map);

    // SZ::writefile("segmentation_map.dat", segmentation_map.data(), segmentation_map.size());

    return 0;

}