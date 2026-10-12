#ifndef SZ3_BIN_STREAM_HPP
#define SZ3_BIN_STREAM_HPP

#include <algorithm>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace SZ3 {

/**
 * Bins read in order, decoded a chunk at a time, so a decomposition that consumes them in order never holds them
 * all (port of sz3-rust's BinStream).
 */
template <class T>
class BinStream {
    static constexpr size_t kChunk = size_t(1) << 14;

   public:
    using Fill = std::function<void(T *, size_t)>;

    BinStream(size_t total, Fill fill) : total_(total), fill_(std::move(fill)) {}

    size_t size() const { return total_; }

    inline T next() {
        if (at_ < end_) return buf_[at_++];
        return refill_next();
    }

    /// The next chunk as [first, second), all of it consumed; throws at the end of the stream. For loops that keep
    /// their position in locals rather than calling next().
    __attribute__((noinline)) std::pair<const T *, const T *> next_chunk() {
        if (at_ < end_) {
            const size_t a = at_;
            at_ = end_;
            return {buf_.data() + a, buf_.data() + end_};
        }
        if (decoded_ == total_) throw std::out_of_range("SZ3: fewer bins than the decomposition consumes");
        refill();
        at_ = end_;
        return {buf_.data(), buf_.data() + end_};
    }

    /// Decodes what is left, so a damaged tail is still found.
    void finish() {
        while (decoded_ < total_) refill();
    }

   private:
    void refill() {
        const size_t m = std::min(total_ - decoded_, kChunk);
        buf_.resize(m);
        fill_(buf_.data(), m);
        decoded_ += m;
        at_ = 0;
        end_ = m;
    }

    __attribute__((noinline)) T refill_next() {
        if (decoded_ == total_) throw std::out_of_range("SZ3: fewer bins than the decomposition consumes");
        refill();
        return buf_[at_++];
    }

    std::vector<T> buf_;
    size_t at_ = 0, end_ = 0;
    size_t total_, decoded_ = 0;
    Fill fill_;
};

}  // namespace SZ3

#endif
