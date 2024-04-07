#include <cstdio>
#include <filesystem>
#include <cstdlib>
#include <cmath>
#include <memory>
#include <string>
#include "SZ3/api/sz.hpp"
#include "SZ3/utils/FileUtil.hpp"



int main(int argc, char** argv)
{

    if (argc < 2) {
        printf("Usage: %s <input_file> <relative_error_bound>\n", argv[0]);
        printf("Example: %s testfloat_8_8_128.dat 32768 0 0\n", argv[0]);
        return 0;
    }

    std::filesystem::path p{argv[1]} ;
    if (!std::filesystem::exists(p)) {
        printf("File %s does not exist\n", argv[1]);
        return 0;
    }

    size_t file_size = std::filesystem::file_size(p)/sizeof(float);


    std::vector<int> quant_inds(file_size, 0);
    std::vector<float> input_data(file_size, 0);

    SZ::readfile(argv[1],  file_size, input_data.data());
    float max = *std::max_element(input_data.begin(), input_data.end());
    float min = *std::min_element(input_data.begin(), input_data.end());
    printf("max: %f, min: %f\n", max, min);
    double eb = atof(argv[2])*(max - min);
    printf("relative eb: %.6f\n", atof(argv[2]));
    printf("absolute eb: %.6f\n", eb);

    // create a linear quantizer
    auto quantizer = SZ::LinearQuantizer<float>();
    quantizer.set_eb(eb);

    // iterate the input data and quantize it
    for (size_t i = 0; i < file_size; i++) {
        quant_inds[i] = quantizer.quantize(input_data[i],0);
    }

    // write the quantized data to a file
    std::string output_file =  p.filename().string()+".quant.i32";
    printf("Writing quantized data to %s\n", output_file.c_str());
    SZ::writefile(output_file.c_str(), quant_inds.data(),file_size);

    return 0;
}


