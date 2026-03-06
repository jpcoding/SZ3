#include <cstddef>
#include<iostream> 
#include<vector> 
#include <string> 
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "mpi.h" 
#include <algorithm>
#include <filesystem>
#include <cmath>
#include "SZ3/api/sz.hpp"
namespace fs = std::filesystem;



template<class T>
void compress(const char *inPath, const char *cmpPath, SZ3::Config conf) {
    T *data = new T[conf.num];
    SZ3::readfile<T>(inPath, conf.num, data);

    size_t outSize;
    char *bytes = SZ_compress<T>(conf, data, outSize);

    char outputFilePath[1024];
    if (cmpPath == nullptr) {
        snprintf(outputFilePath, 1024, "%s.sz", inPath);
    } else {
        strcpy(outputFilePath, cmpPath);
    }
    SZ3::writefile(outputFilePath, bytes, outSize);

    // printf("compression ratio = %.2f \n", conf.num * 1.0 * sizeof(T) / outSize);
    // printf("compression time = %f\n", compress_time);
    // printf("compressed data file = %s\n", outputFilePath);

    delete[]data;
    delete[]bytes;
}


template<class T>
void decompress(const char *inPath, const char *cmpPath, const char *decPath,
                SZ3::Config conf,
                int binaryOutput = 1 , int printCmpResults= 0) {
    size_t cmpSize;
    auto cmpData = SZ3::readfile<char>(cmpPath, cmpSize);
    T *decData = SZ_decompress<T>(conf, cmpData.get(), cmpSize);
    char outputFilePath[1024];
    if (decPath == nullptr) {
        snprintf(outputFilePath, 1024, "%s.out", cmpPath);
    } else {
        strcpy(outputFilePath, decPath);
    }
    if (binaryOutput == 1) {
        SZ3::writefile<T>(outputFilePath, decData, conf.num);
    } else {
        SZ3::writeTextFile<T>(outputFilePath, decData, conf.num);
    }
    delete[]decData;
}

template<typename T>
double get_psnr(T* odata, T* ddata, size_t size)
{
	T omax = *std::max_element(odata, odata+size);
	T omin = *std::min_element(ddata, ddata+size);
	double range = omax - omin;
	double mse = 0;
	for(size_t i = 0; i< size; i++)
	{
		mse+=(odata[i]-ddata[i])*(odata[i]-ddata[i]);
	}
	mse = mse/size;
	double psnr = 20*std::log10(range) - 10*std::log10(mse);
	return psnr;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);  
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm comm = MPI_COMM_WORLD;

    SZ3::Config conf;
    size_t dim0, dim1, dim2; 
    std::string data_dir = argv[1];

    std::string output_dir = argv[2];

    dim0 = atoi(argv[3]); // fastest changing dimension 
    dim1 = atoi(argv[4]);
    dim2 = atoi(argv[5]); 

    // construct the config for all the data 
    conf = SZ3::Config(dim2, dim1, dim0);
    conf.loadcfg(argv[6]);
    double rel_eb = atof(argv[7]);

    std::string data_ext = ".f32";
    // create a global array to store the file paths
    std::vector<std::string> global_file_names;
    // read the file names from the directory

    for (const auto & entry : fs::directory_iterator(data_dir))
    {
        if (entry.path().extension() == data_ext)
        {
            // std::cout << entry.path().string() << std::endl;
            global_file_names.push_back(entry.path().string());
        }
    }

    
    // fallten the vector of strings to a char array 
    double local_min = std::numeric_limits<double>::max();
    double local_max = std::numeric_limits<double>::lowest();

    int num_files = global_file_names.size();  // Total number of files
    int files_per_rank = num_files / size;     // Basic number of files per rank
    int remainder = num_files % size;          // Remainder files to be distributed

    // Calculate the start and end file indices for each rank
    int local_file_index_start, local_file_index_end;

    if (rank < remainder) {
        // Ranks with extra files (remainder) get one extra file
        local_file_index_start = rank * (files_per_rank + 1);
        local_file_index_end = local_file_index_start + (files_per_rank + 1) - 1;
    } else {
        // Ranks without extra files
        local_file_index_start = rank * files_per_rank + remainder;
        local_file_index_end = local_file_index_start + files_per_rank - 1;
    }

    // std::cout << "Rank: " << rank << " Start: " << local_file_index_start << " End: " << local_file_index_end << std::endl;

    MPI_Barrier(comm);

    auto timer = SZ3::Timer();
    std::vector<float> data(conf.num); 

    // for (int i = local_file_index_start; i <= local_file_index_end; i++)
    // {
    //     std::string current_file = global_file_names[i];
    //     double file_min, file_max;
    //     // if(rank == 0)
    //     // {
    //     //     std::cout << "Rank: " << rank << " File: " << current_file << std::endl;
    //     // }
    //     SZ3::readfile<float>(current_file.c_str(), conf.num, data.data()); // should be swiched to posix io? 
    //     file_min = *std::min_element(data.begin(), data.end());
    //     file_max = *std::max_element(data.begin(), data.end());
    //     local_min = std::min(local_min, file_min);
    //     local_max = std::max(local_max, file_max);
    //     // std::cout << "Rank: " << rank << " File: " << current_file << " Min: " << file_min << " Max: " << file_max << std::endl;
    // }

    // double global_min, global_max;
    // MPI_Barrier(comm);
    // MPI_Reduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
    // MPI_Reduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    // MPI_Bcast(&global_min, 1, MPI_DOUBLE, 0, comm);
    // MPI_Bcast(&global_max, 1, MPI_DOUBLE, 0, comm);
    // MPI_Barrier(comm);
    // if(rank == 0)
    // {
    //     std::cout << "Global min: " << global_min << std::endl;
    //     std::cout << "Global max: " << global_max << std::endl;
    // }
    // get the global absolute eb 
    // double global_eb = (global_max - global_min) * rel_eb;
    // update config 
    // conf.absErrorBound = rel_eb;
    conf.errorBoundMode = SZ3::EB_REL; 
    conf.relErrorBound  = rel_eb;

    // compress the data 
    size_t local_outSize_sum = 0;
    size_t local_outSize = 0; 
    size_t global_outSize = 0; 
    size_t data_set_szie = conf.num * sizeof(float) *global_file_names.size() ;

    MPI_Barrier(comm);
    double time_start, time_end, compression_time;
    double read_orig_time = 0;
    double write_compressed_time = 0;
    double read_compressed_time = 0;
    double write_decompressed_time = 0;
    double sz_time_compress = 0;

    compression_time = 0;

    for (int i = local_file_index_start; i <= local_file_index_end; i++)
    {
        local_outSize = 0;
        std::string current_file = global_file_names[i];
        std::string outputFilePath = output_dir + "/" + fs::path(current_file).filename().string() + ".sz";
        time_start = MPI_Wtime();
        SZ3::readfile<float>(current_file.c_str(), conf.num, data.data());
        time_end = MPI_Wtime();
        read_orig_time += time_end - time_start;
        time_start = MPI_Wtime();

        timer.start();
        char *bytes = SZ_compress<float>(conf, data.data(), local_outSize);
        sz_time_compress += timer.stop();

        time_end = MPI_Wtime();
        compression_time += time_end - time_start;
        local_outSize_sum += local_outSize;
        time_start = MPI_Wtime();
        SZ3::writefile(outputFilePath.c_str(), bytes, local_outSize);
        time_end = MPI_Wtime();
        write_compressed_time += time_end - time_start;
        free(bytes);
    }

    MPI_Barrier(comm);
    double global_compression_time = 0;
    double global_sz_time_compress = 0;
    MPI_Reduce(&sz_time_compress, &global_sz_time_compress, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&local_outSize_sum, &global_outSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, comm);
    MPI_Reduce(&compression_time, &global_compression_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Barrier(comm);

        
    // decompress the data
    MPI_Barrier(comm);
    double decompression_time = 0;

    size_t cmpSize = 0;
   
    for (int i = local_file_index_start; i <= local_file_index_end; i++)
    {
        std::string current_file = global_file_names[i];
        std::string outputFilePath = output_dir + "/" + fs::path(current_file).filename().string() + ".sz";
        std::string decompressed_file = output_dir + "/" + fs::path(current_file).filename().string() + ".sz.out";
        time_start = MPI_Wtime();
        auto cmpData = SZ3::readfile<char>(outputFilePath.c_str(), cmpSize);
        time_end = MPI_Wtime();
        read_compressed_time += time_end - time_start;
        time_start = MPI_Wtime();
        float *decData = SZ_decompress<float>(conf, cmpData.get(), cmpSize);
        time_end = MPI_Wtime();
        decompression_time += time_end - time_start;
        time_start = MPI_Wtime();
        SZ3::writefile<float>(decompressed_file.c_str(), decData, conf.num);
        time_end = MPI_Wtime();
        write_decompressed_time += time_end - time_start;
        free(decData);
    }

    MPI_Barrier(comm);
    double global_decompression_time = 0;
    double global_read_orig_time = 0;
    double global_write_compressed_time = 0;
    double global_read_compressed_time = 0;
    double global_write_decompressed_time = 0;

    MPI_Reduce(&decompression_time, &global_decompression_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&read_orig_time, &global_read_orig_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&write_compressed_time, &global_write_compressed_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&read_compressed_time, &global_read_compressed_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&write_decompressed_time, &global_write_decompressed_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Barrier(comm);
    // calculate psnr
    //
    double local_psnr = 0;
   double current_psnr = 0; 
    std::vector<float> ddata(conf.num);
    for (int i = local_file_index_start; i <= local_file_index_end; i++)
    {
        std::string current_file = global_file_names[i];
        std::string outputFilePath = output_dir + "/" + fs::path(current_file).filename().string() + ".sz";
        std::string decompressed_file = output_dir + "/" + fs::path(current_file).filename().string() + ".sz.out";
	
	SZ3::readfile<float>(current_file.c_str(), conf.num, data.data());
	
	SZ3::readfile<float>(decompressed_file.c_str(), conf.num, ddata.data());

	current_psnr = get_psnr(data.data(), ddata.data(), conf.num);
	std::cout << "current psnr = " << current_psnr << std::endl;
	local_psnr += std::pow(10, -0.1*current_psnr);
    }

   double global_psnr = 0;
   MPI_Barrier(comm);
   MPI_Reduce(&local_psnr, &global_psnr, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
   global_psnr = 10*std::log10(num_files/global_psnr);

    if(rank == 0)
    {
        std::cout << "Compression time =  " << global_compression_time << std::endl;
        std::cout << "Decompression time = " << global_decompression_time << std::endl;
        std::cout << "Global compressed size: " << global_outSize << std::endl;
        std::cout << "Compression ratio: " << data_set_szie*1.0/ global_outSize << std::endl;
        std::cout << "Read original time: " << global_read_orig_time << std::endl;
        std::cout << "Write compressed time: " << global_write_compressed_time << std::endl;
        std::cout << "Read compressed time: " << global_read_compressed_time << std::endl;
        std::cout << "Write decompressed time: " << global_write_decompressed_time << std::endl;
        std::cout << "SZ compress time: " << global_sz_time_compress << std::endl;
	std::cout << "PSNR = " << global_psnr << std::endl;
    }


    MPI_Finalize();
    return 0;


}
