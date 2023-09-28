#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include <cstddef>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <exception>
#include <vector>
/*
This file is only meant for the helper functions for SZInterpolationCompressor.hpp
Do not use the functions in this files in eleswhere.
*/

namespace SZ
{
    template<typename T>
    size_t compute_aux_diff_compress(const SZ::Config conf, const T *data, const int detection_blocksize, 
    std::vector<uchar> &significant_block,  std::vector<uchar> &significant_index )
    {
        // 2D AND 3D ONLY
        if(conf.N ==2)
        {
            size_t significant_point_counter = 0; 
            int block_size = detection_blocksize;
            size_t num_elements = conf.num; 
            auto dims = conf.dims; 
            // dims[-1] is the fastest changing dimension 
            int nx = (int)ceil(1.0 * dims[0] / block_size);
            int ny = (int)ceil(1.0 * dims[1] / block_size);
            //resize the two vectors
            significant_block.resize(nx * ny, 0);
            significant_index.resize(num_elements, 0);  
            //compute the root of sum of squares of differences of each block 
            //and mark the significant blocks accoring to the threshold from config 
            const T *x_data_pos = data;
            for (int i = 0; i< nx; i++)
            {
                int block_size_i = (i + 1) * block_size > dims[0] ? dims[0] - i * block_size : block_size;
                const T *y_data_pos = x_data_pos;
                for (int j = 0; j< ny; j++)
                {
                    int block_size_j = (j + 1) * block_size > dims[1] ? dims[1] - j * block_size : block_size;
                    // compute the first derivative of x direction
                    double diff_x = 0;
                    double diff_y = 0; 
                    const T *xx_data_pos = y_data_pos; // block begin index 
                    for (int ii = 1; ii < block_size_i-1; ii++)
                    {
                        const T *yy_data_pos = xx_data_pos; // block begin index 
                        for (int jj = 1; jj < block_size_j-1; jj++)
                        {
                            //along y direction
                            diff_y += pow(fabs(*(yy_data_pos+1) - *(yy_data_pos-1)),2);
                            // along x direction 
                            diff_x += pow(fabs(*(yy_data_pos+dims[1]) - *(yy_data_pos-dims[1])),2);
                            yy_data_pos += 1;
                        }
                        xx_data_pos += dims[1];
                    }
                    // check if this block is significantly unsmooth 
                    // if the central diff is larger than the threshold, then it is significant 
                    diff_x = sqrt(diff_x)/((block_size_i-2)*(block_size_j-2));
                    diff_y = sqrt(diff_y)/((block_size_i-2)*(block_size_j-2));
                    if(diff_x > conf.diff_thresh || diff_y > conf.diff_thresh)
                    {
                        significant_block[i * ny + j] = 1;
                        // mark the significant index 
                        for (int ii = 0; ii < block_size_i; ii++)
                        {
                            for (int jj = 0; jj < block_size_j; jj++)
                            {
                                significant_index[(i * block_size_i + ii) * dims[1] + j * block_size_j + jj] = 1;   
                                significant_point_counter++;
                            }
                        }
                    }
                    y_data_pos += block_size_j;
                }
                x_data_pos += block_size_i * dims[1];
            }
            // return the number of significant blocks
            return significant_point_counter;
        }



        if (conf.N == 3)
        {
            size_t significant_point_counter = 0; 
            int block_size = detection_blocksize;
            size_t num_elements = conf.num; 
            auto dims = conf.dims; 
            // dims[-1] is the fastest changing dimension 
            int nx = (int)ceil(1.0 * dims[0] / block_size);
            int ny = (int)ceil(1.0 * dims[1] / block_size);
            int nz = (int)ceil(1.0 * dims[2] / block_size);
            //resize the two vectors
            significant_block.resize(nx * ny * nz, 0);
            significant_index.resize(num_elements, 0);  
            //compute the root of sum of squares of differences of each block 
            //and mark the significant blocks accoring to the threshold from config 
            const T *x_data_pos = data;
            for (int i = 0; i< nx; i++)
            {
                int block_size_i = (i + 1) * block_size > dims[0] ? dims[0] - i * block_size : block_size;
                const T *y_data_pos = x_data_pos;
                for (int j = 0; j< ny; j++)
                {
                    int block_size_j = (j + 1) * block_size > dims[1] ? dims[1] - j * block_size : block_size;
                    const T *z_data_pos = y_data_pos;
                    for (int k = 0; k< nz; k++)
                    {
                        int block_size_k = (k + 1) * block_size > dims[2] ? dims[2] - k * block_size : block_size;
                        // compute the first derivative of x direction
                        double diff_x = 0;
                        double diff_y = 0; 
                        double diff_z = 0; 
                        const T *xx_data_pos = z_data_pos; // block begin index 
                        for (int ii = 1; ii < block_size_i-1; ii++)
                        {
                            const T *yy_data_pos = xx_data_pos; // block begin index 
                            for (int jj = 1; jj < block_size_j-1; jj++)
                            {
                                const T *zz_data_pos = yy_data_pos; // block begin index 
                                for (int kk = 1; kk < block_size_k-1; kk++)
                                {
                                    //along y direction
                                    diff_x += pow((*(zz_data_pos+1) - *(zz_data_pos-1))/2,2);
                                    // along x direction 
                                    diff_y += pow((*(zz_data_pos+dims[2]) - *(zz_data_pos-dims[2]))/2,2);
                                    // along z direction 
                                    diff_z += pow((*(zz_data_pos+dims[1]*dims[2]) - *(zz_data_pos-dims[1]*dims[2]))/2,2);
                                    zz_data_pos += 1;
                                }
                                yy_data_pos += dims[2];
                            }
                            xx_data_pos += dims[1]*dims[2];
                        }
                        // check if this block is significantly unsmooth
                        // if the central diff is larger than the threshold, then it is significant
                        int num_sum = (block_size_i-2)*(block_size_j-2)*(block_size_k-2);
                        if(num_sum > 0)
                        {
                            diff_x = sqrt(diff_x/num_sum);
                            diff_y = sqrt(diff_y/num_sum);
                            diff_z = sqrt(diff_z/num_sum);
                        }
                        else
                        {
                            diff_x = 0;
                            diff_y = 0;
                            diff_z = 0;
                        }

                        if(diff_x > conf.diff_thresh || diff_y > conf.diff_thresh || diff_z > conf.diff_thresh)
                        {
                            significant_block[i * ny * nz + j * nz + k] = 1;
                            // mark the significant index 
                            for (int ii = 0; ii < block_size_i; ii++)
                            {
                                for (int jj = 0; jj < block_size_j; jj++)
                                {
                                    for (int kk = 0; kk < block_size_k; kk++)
                                    {
                                        significant_index[(i * block_size_i + ii) * dims[1]*dims[2] + j * block_size_j*dims[2] + k * block_size_k + kk] = 1;   
                                        significant_point_counter++;
                                    }
                                }
                            }
                        }
                        z_data_pos += block_size_k;
                    }
                    y_data_pos += block_size_j*dims[2];
                }
                x_data_pos += block_size_i * dims[1]*dims[2];
            }
            // return the number of significant blocks
            std::cout<< "significant_point_counter: " << significant_point_counter << std::endl;
            return significant_point_counter;
        }
        else {
            std::cout << "Error: only support 2D and 3D data" << std::endl;
            exit(0);
        }

    };

    /*
    This function recover the significant_index using significant_block with given blcoksize
    */
    size_t compute_aux_diff_decompress(size_t *global_dims, int N, const int detection_blocksize, 
    std::vector<uchar> &significant_block,  std::vector<uchar> &significant_index)
    {
        
        size_t num_elements = 1;
        for (int i = 0; i < N; i++)
        {
            num_elements *= global_dims[i];
        }
        if(N ==2)
        {
            size_t significant_point_counter = 0; 
            int block_size = detection_blocksize;
            auto dims = global_dims; 
            // dims[-1] is the fastest changing dimension 
            int nx = (int)ceil(1.0 * dims[0] / block_size);
            int ny = (int)ceil(1.0 * dims[1] / block_size);
            //resize the two vectors
            significant_index.resize(num_elements, 0);  
            //compute the root of sum of squares of differences of each block 
            //and mark the significant blocks accoring to the threshold from config 
            size_t block_id = 0;
            size_t index_id = 0; 
            for (int i = 0; i< nx; i++)
            {
                int block_size_i = std::min(block_size, (int)dims[0] - i * block_size);
                for (int j = 0; j< ny; j++)
                {
                    int block_size_j = std::min(block_size, (int)dims[1] - j * block_size);
                    if(significant_block[block_id] == 1)
                    {
                        // mark the significant index 
                        index_id = (i * block_size_i) * dims[1] + j * block_size_j;
                        for (int ii = 0; ii < block_size_i; ii++)
                        {
                            for (int jj = 0; jj < block_size_j; jj++)
                            {
                                significant_index[index_id] = 1;   
                                significant_point_counter++;
                                index_id++;
                            }
                            index_id += dims[1];
                        }
                    }
                    block_id++;
                }
                block_id += ny;
            }
            // return the number of signific·ant blocks
            return significant_point_counter;
        }
        else if (N ==3)
        {
            size_t significant_point_counter = 0; 
            int block_size = detection_blocksize;
            auto dims = global_dims;; 
            // dims[-1] is the fastest changing dimension 
            int nx = (int)ceil(1.0 * dims[0] / block_size);
            int ny = (int)ceil(1.0 * dims[1] / block_size);
            int nz = (int)ceil(1.0 * dims[2] / block_size);
            //resize the two vectors
            significant_index.resize(num_elements, 0);  
            //compute the root of sum of squares of differences of each block 
            //and mark the significant blocks accoring to the threshold from config
            for (int i = 0; i< nx; i++)
            {
                int block_size_i = (i + 1) * block_size > dims[0] ? dims[0] - i * block_size : block_size;
                for (int j = 0; j< ny; j++)
                {
                    int block_size_j = (j + 1) * block_size > dims[1] ? dims[1] - j * block_size : block_size;
                    for (int k = 0; k< nz; k++)
                    {
                        int block_size_k = (k + 1) * block_size > dims[2] ? dims[2] - k * block_size : block_size;
                        if(significant_block[i * ny * nz + j * nz + k] == 1)
                        {
                            // mark the significant index 
                            for (int ii = 0; ii < block_size_i; ii++)
                            {
                                for (int jj = 0; jj < block_size_j; jj++)
                                {
                                    for (int kk = 0; kk < block_size_k; kk++)
                                    {
                                        significant_index[(i * block_size_i + ii) * dims[1]*dims[2] + j * block_size_j*dims[2] + k * block_size_k + kk] = 1;   
                                        significant_point_counter++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // return the number of significant blocks
            return significant_point_counter;
        }
        else {
            std::cout << "Error: only support 2D and 3D data" << std::endl;
            return 0;
            exit(0);
        }

    };
} // namespace SZ




