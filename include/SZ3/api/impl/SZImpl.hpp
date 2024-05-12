#ifndef SZ3_IMPL_SZ_HPP
#define SZ3_IMPL_SZ_HPP

#include "SZ3/def.hpp"
#include "SZ3/api/impl/SZDispatcher.hpp"
#include "SZ3/api/impl/SZImplOMP.hpp"
#include "SZ3/utils/Timer.hpp"
#include <cmath>
#include <memory>

template<class T, SZ::uint N>
char *SZ_compress_impl(SZ::Config &conf, const T *data, size_t &outSize) {
#ifndef _OPENMP
    conf.openmp=false;
#endif
    if (conf.openmp) {
        //dataCopy for openMP is handled by each thread
        return SZ_compress_OMP<T, N>(conf, data, outSize);
    } else {
        // std::vector<T> dataCopy(data, data + conf.num);
        // std::unique_ptr<std::vector<T>> data_copy = std::make_unique<std::vector<T>>(conf.num);

        std::shared_ptr<std::vector<T>> data_copy = std::make_shared<std::vector<T>>(conf.num);

        std::copy(data, data + conf.num, data_copy->begin());
        conf.PASS_DATA.original_data_prt = (const void *)data;
        char *ret = SZ_compress_dispatcher<T, N>(conf, data_copy->data(), outSize);
        conf.PASS_DATA.processed_data_prt = std::static_pointer_cast<void>(data_copy); 
        return ret;
    }
}


template<class T, SZ::uint N>
void SZ_decompress_impl(SZ::Config &conf, char *cmpData, size_t cmpSize, T *decData) {


#ifndef _OPENMP
    conf.openmp=false;
#endif
    if (conf.openmp) {
        SZ_decompress_OMP<T, N>(conf, cmpData, cmpSize, decData);
    } else {
        SZ_decompress_dispatcher<T, N>(conf, cmpData, cmpSize, decData);
    }
}

#endif