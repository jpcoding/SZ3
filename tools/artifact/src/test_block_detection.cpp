#include <iostream>
#include <vector>
#include <array>
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "block_detection.hpp"

int main(int argc, char** argv)
{
    SZ::Timer timer;

    if (argc < 2)
    {
        std::cout << "Usage: ./test_block_detection <input_file> <N> <dim1 dim2 ...> " << std::endl;
        return 0;
    }

    int N = atoi(argv[2]);
    int num_elements = 1;
    std::vector<int> global_dimensions;
    for (int i = 0; i < N; i++)
    {
        global_dimensions.push_back(atoi(argv[3 + i]));
        num_elements *= global_dimensions[i];
    }

    std::vector<float> input_data(num_elements);
    SZ::readfile<float>(argv[1], num_elements, input_data.data());

    // block detection
    int begin = 1600*global_dimensions[0];
    int end = global_dimensions[0];
    int stride = 1;
    auto block_detection = BlockDetection<float>(input_data.data(), begin, end, stride);
    timer.start();
    auto block_size = block_detection.block_detection(1e-5);
    timer.stop("block detection");
    std::cout << "block size" << block_size << std::endl;

    return 0;


}