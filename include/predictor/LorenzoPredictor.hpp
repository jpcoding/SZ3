#ifndef _SZ_LORENZO_PREDICTOR_HPP
#define _SZ_LORENZO_PREDICTOR_HPP

#include <cassert>
#include "def.hpp"
#include "utils/Iterator.hpp"

namespace SZ{

  // N-d lorenzo predictor
  template <class T, uint N>
  class LorenzoPredictor;

  namespace {
    template <class T, uint N>
    class LorenzoBase {
    public:
      static const uint8_t predictor_id = 0b00000001;
      using Range = multi_dimensional_range<T, N>;
      using iterator = typename multi_dimensional_range<T, N>::iterator;
      void precompress_data(const iterator&) {}
      void postcompress_data(const iterator&) {}
      void predecompress_data(const iterator&) {}
      void postdecompress_data(const iterator&) {}
      void precompress_block(const std::shared_ptr<Range>&) {}
      void predecompress_block(const std::shared_ptr<Range>&) {}

      /*
       * save doesn't need to store anything except the id
       */
      // std::string save() const {
      //   return std::string(1, predictor_id);
      // }
      void save(uchar*& c) const{
        c[0] = predictor_id;
        c += sizeof(uint8_t);
        *reinterpret_cast<T*>(c) = noise;
        c += sizeof(T);
      }

      /*
       * just verifies the ID, increments
       */
      // static LorenzoPredictor<T,N> load(const unsigned char*& c, size_t& remaining_length) {
      //   assert(remaining_length > sizeof(uint8_t));
      //   c += 1;
      //   remaining_length -= sizeof(uint8_t);
      //   return LorenzoPredictor<T,N>{};
      // }
      void load(const uchar*& c, size_t& remaining_length){
        c += sizeof(uint8_t);
        remaining_length -= sizeof(uint8_t);
        noise = *reinterpret_cast<const T*>(c);
        c += sizeof(T);
      }
      void print() const{
        std::cout << "Lorenzo predictor\n";
      }
    protected:
      T noise = 0;
    };
  }


  template <class T>
  class LorenzoPredictor<T, 1> : public LorenzoBase<T,1> {
  public:
    LorenzoPredictor(){
      this->noise = 0;
    }
    LorenzoPredictor(T eb){
      this->noise = 0.5*eb;
    }
    using iterator = typename LorenzoBase<T,1>::iterator;
    inline T estimate_error(const iterator& iter) const noexcept {
      return ABS(*iter - predict(iter)) + this->noise;
    }
    inline T predict(const iterator& iter) const noexcept { return iter.prev(0); };
  };

  template <class T>
  class LorenzoPredictor<T, 2>: public LorenzoBase<T,2> {
  public:
    LorenzoPredictor(){
      this->noise = 0;
    }
    LorenzoPredictor(T eb){
      this->noise = 0.81*eb;
    }
    using iterator = typename LorenzoBase<T,2>::iterator;
    inline T estimate_error(const iterator& iter) const noexcept {
      return ABS(*iter - predict(iter)) + this->noise;
    }
    inline T predict(const iterator& iter) const noexcept{
      return iter.prev(0, 1) + iter.prev(1, 0) - iter.prev(1, 1);
    };
  };

  template <class T>
  class LorenzoPredictor<T, 3>: public LorenzoBase<T,3> {
  public:
    LorenzoPredictor(){
      this->noise = 0;
    }
    LorenzoPredictor(T eb){
      this->noise = 1.22*eb;
    }
    using iterator = typename LorenzoBase<T,3>::iterator;
    inline T estimate_error(const iterator& iter) const noexcept {
      return ABS(*iter - predict(iter)) + this->noise;
    }
    inline T predict(const iterator& iter) const noexcept{
      return iter.prev(0, 0, 1) + iter.prev(0, 1, 0) + iter.prev(1, 0, 0) 
          - iter.prev(0, 1, 1) - iter.prev(1, 0, 1) - iter.prev(1, 1, 0)
          + iter.prev(1, 1, 1);
    };
  };

}
#endif