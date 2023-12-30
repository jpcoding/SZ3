#include <cstddef>
#include<vector>
#include<iostream>
#include<random>
#include "SZ3/api/sz.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "omp.h"

auto global_rng = std::mt19937(1);

template <typename T> 
int create_initial_model(std::vector<T> &model, int Nt, int Ny, int Nx)
{
    int N = Nt*Nx*Ny;
    for (int i = 0; i < N; i++)
    {
        model[i] = std::uniform_real_distribution<>(10,30)(global_rng);
    }
    // set boundary conditions 
    // top -100.o C 
    // bottom 100.0 C
    // left 50.0 C
    // right -20.0 C
    for (int i = 0; i < Nt; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            model[i*Nx*Ny + j] = -100.0;
            model[(i+1)*Nx*Ny - Nx + j] = 100.0;
        }
    }
    for (int i = 0; i < Nt; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            model[i*Nx*Ny + j*Nx] = 50.0;
            model[i*Nx*Ny + j*Nx + Nx - 1] = -20.0;
        }
    }
    return 0;
}

template <typename T> 
int one_step_forward(std::vector<T> &model, int step,int Ny, int Nx)
{
    int current_start = step*Nx*Ny;
    int i ,j;
    // #pragma omp parallel for private (j) 
    for( i = 1; i<Ny-1; i++) 
    {
        for( j = 1; j<Nx-1; j++)
        {
            model[(step+1)*Ny*Nx + i*Nx + j] = 0.25*(model[current_start + i*Nx + j + 1] + model[current_start + i*Nx + j - 1] 
            + model[current_start + (i+1)*Nx + j] + model[current_start + (i-1)*Nx + j]);
        }
        
    }
    return 0;
}

// define the compressor
template <typename T>
void compress(T *data, int step,  int Ny, int Nx, double abs_eb, size_t &outSize , char** config_file ) {
    // get the copy of step data 
    std::vector<T> model(Ny*Nx, 0.0);
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            model[i*Nx + j] = data[step*Ny*Nx + i*Nx + j];
        }
    }

    auto conf = SZ::Config(Ny, Nx);
    conf.loadcfg(*config_file);
    // conf.errorBoundMode = SZ::EB_REL;
    conf.errorBoundMode = SZ::EB_ABS;
    conf.absErrorBound = abs_eb;
    // conf.relErrorBound = abs_eb;


    char *bytes = SZ_compress<T>(conf, model.data(), outSize);

    // T* cmpData = (T*) malloc(Ny*Nx*sizeof(T));
    T* cmpData = SZ_decompress<T>(conf, bytes, outSize);
 
    // copy data back to the original array]
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            data[step*Ny*Nx + i*Nx + j] = cmpData[i*Nx + j];
        }
    }
    delete[] bytes;
    delete[] cmpData;

}


int main(int argc, char **argv)
{
    int Nt = atoi(argv[3]);
    int Ny = atoi(argv[2]);
    int Nx = atoi(argv[1]);
    int N = Nt*Ny*Nx;
    double abs_eb = (double) atof(argv[4]);
    int outSize = 0;
    char *config_file = argv[5];
    std::vector<float> model(N, 0.0);
    create_initial_model<float>(model, Nt, Ny, Nx);
    // make a copy of the model
    std::vector<float> model_copy(N, 0.0);
    for (int i = 0; i < N; i++)
    {
        model_copy[i] = model[i];
    }

    std::cout << "[1]value at 0 20 20 " << model[0*Ny*Nx + 20*Nx + 20] << "\n";

    SZ::writefile("standard_model1.dat", model.data(), model.size());
    std::vector<size_t> compressed_sizes(Nt, 0);

    // get a standard model first
    std::cout << "model size is " << model.size() << "\n";
    for (int i = 0; i < Nt-1; i++)
    {
        one_step_forward<float>(model, i, Ny, Nx);
    }

    SZ::writefile("standard_model.dat", model.data(), model.size());

    // create_initial_model<float>(model, Nt, Ny, Nx);

    std::cout << "[2]value at 0 20 20 " << model_copy[0*Ny*Nx + 20*Nx + 20] << "\n";
    model = model_copy;
    for (int i = 0; i < Nt-1; i++)
    {
        compress<float>(model.data(), i, Ny, Nx, abs_eb, compressed_sizes[i], &config_file);
        one_step_forward<float>(model, i, Ny, Nx);
    }

    // compress the last slice 
    compress<float>(model.data(), Nt-1, Ny, Nx, abs_eb, compressed_sizes[Nt-1], &config_file);

    // write the compressed sizes to a file

    SZ::writefile("decompressed_data.dat", model.data(), model.size());

    // get the compression ratio 
    double  compression_ratio = 0.0;
    for (int i = 0; i < Nt; i++)
    {
        compression_ratio += (float)compressed_sizes[i];
    }
    compression_ratio = (Nt*Nx*Ny*sizeof(float))/compression_ratio;

    std::cout << "compression ratio is " << compression_ratio << "\n";


    
    return 0;
}



