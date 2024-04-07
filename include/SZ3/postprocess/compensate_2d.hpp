
class 2DCompensationusing 
{

};


std::vector<char> build_boundary_map(
                int* quant_inds,
                const size_t* begin, 
                const size_t* end, 
                const std::array<int,3> & dims,
                int plane_dir0, int plane_dir1, 
                int plane_dir0_stride, int plane_dir1_stride,
                const int quantizer_radius)
{
    size_t plane_dim0 = (end[dims[plane_dir0]] - begin[dims[plane_dir0]]) / plane_dir0_stride + 1; // stride = plane dim1
    size_t plane_dim1 = (end[dims[plane_dir1]] - begin[dims[plane_dir1]]) / plane_dir1_stride + 1; // stride = 1
    std::vector<char> boundary_map(plane_dim0*plane_dim1, 0);
    for(size_t i = 1; i < plane_dim0-1; i++)
    {
        size_t idx_i = begin[dims[plane_dir0]] + i*plane_dir0_stride;
        for(size_t j = 1; j < plane_dim1-1; j++)
        {
            size_t idx_j = idx_i + j*plane_dir1_stride;
            /*
            check the neighbors of X 
            # A #
            B X C
            # D #
            */
            int X = quant_inds[idx_j]==0? quantizer_radius: quant_inds[idx_j];
            int A = quant_inds[idx_j - plane_dir1_stride]!=0? quant_inds[idx_j - plane_dir1_stride]: quantizer_radius;
            int B = quant_inds[idx_j - plane_dir0_stride]!=0? quant_inds[idx_j - plane_dir0_stride]: quantizer_radius;
            int C = quant_inds[idx_j + plane_dir1_stride]!=0? quant_inds[idx_j + plane_dir1_stride]: quantizer_radius;
            int D = quant_inds[idx_j + plane_dir0_stride]!=0? quant_inds[idx_j + plane_dir0_stride]: quantizer_radius;
            if((X!=A) || (X!=B) || (X!=C) || (X!=D))
            {
                boundary_map[i*plane_dim1 +j] = 1;
            }
        }
    }
    return boundary_map;
}

std::array<char,4> get_compensation_direction(
                        int* quant_inds,
                        int x, int y, 
                        const size_t* begin, 
                        const size_t* end, 
                        const std::array<int,3> & dims,
                        int plane_dir0, int plane_dir1, 
                        int plane_dir0_stride, int plane_dir1_stride,
                        const int quantizer_radius)
{
    std::array<char,4> compensation_direction = {0, 0, 0, 0}; // left right up down



}

std::array<int, 4> get_distance_to_boundary(std::vector<char> & boundary, const int x, const int y, int dim0, int dim1 , int max_extend )
{
    // ignore 4 borders 
    if(x == 0 || x == dim0-1 || y == 0 || y == dim1-1)
    {
        return {max_extend, max_extend, max_extend, max_extend};
    }
    // left right up down
    if(boundary[x*dim1 + y] == 0)
    {
        return {0, 0, 0, 0};
    }
    std::array<int, 4> distance = {0, 0, 0, 0}; // left right up down 
    char current_status = boundary[x*dim1 + y];
    int dx, dy;
    // left
    dy = y;
    dx = x -1;
    while(boundary[dx*dim1 + dy] == current_status)
    {
        dx--;
        distance[0]++;
        if(distance[0] == max_extend)
        {
            break;
        }
    }
    // right
    dy = y;
    dx = x + 1;
    while(boundary[dx*dim1 + dy] == current_status)
    {
        dx++;
        distance[1]++;
        if(distance[1] == max_extend)
        {
            break;
        }
    }
    // up
    dx = x; 
    dy = y - 1;
    while(boundary[dx*dim1 + dy] == current_status)
    {
        dy--;
        distance[2]++;
        if(distance[2] == max_extend)
        {
            break;
        }
    }
    // down
    dx = x;
    dy = y + 1;
    while(boundary[dx*dim1 + dy] == current_status)
    {
        dy++;
        distance[3]++;
        if(distance[3] == max_extend)
        {
            break;
        }
    }
    return distance;
}












template <typename T>
int compensation3d_data_2d_plane(T*data, int* quant_inds, 
                const size_t* begin, 
                const size_t* end, 
                const std::array<int,3> & dims,
                std::array<size_t,3> & dimension_offsets, 
                int interp_direction,
                int interpolation_stride,   
                int plane_dir0, int plane_dir1, 
                int plane_dir0_stride, int plane_dir1_stride,
                const double compensation_max, 
                const int quantizer_radius)

{
    return 0;

}