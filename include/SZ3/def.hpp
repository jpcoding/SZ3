#ifndef SZ3_DEF_HPP
#define SZ3_DEF_HPP

#include <cfloat>

namespace SZ3 {

typedef unsigned int uint;
typedef unsigned char uchar;
#define SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH \
    "The buffer for compressed data is not large enough."

#ifdef _MSC_VER
#define ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define ALWAYS_INLINE inline
#endif

// Inlines the calls of the statement it precedes. For the per-element calls of SZ3's iterators: Clang otherwise
// leaves a large lambda out of line and calls it per element, which GCC inlines.
#if defined(__clang__) && __clang_major__ >= 15
#define SZ3_INLINE_CALLS [[clang::always_inline]]
#else
#define SZ3_INLINE_CALLS
#endif

// In nofma(a * b) + c the product is rounded before the add, so the compiler cannot merge them into one
// FMA instruction, which rounds only once; otherwise a build with FMA and one without decompress different values.
template <class T>
ALWAYS_INLINE T nofma(T x) {
    // Not __builtin_assoc_barrier: GCC 13 and 14 drop it when they vectorize the loop, and fuse the product anyway.
    // An empty asm keeps the value in a register the compiler cannot see through, which a volatile would do only
    // through a store and a load per call.
#if (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || (defined(__i386__) && defined(__SSE2__)))
#if defined(__FMA__) || defined(__FMA4__)
    __asm__("" : "+x"(x));
#endif
    // without FMA instructions there is nothing to contract into, and the value stays free to vectorize
    return x;
#elif (defined(__GNUC__) || defined(__clang__)) && defined(__aarch64__)
    __asm__("" : "+w"(x));
    return x;
#else
    volatile T r = x;
    return r;
#endif
}

// Refuses x87 extended precision and GCC's fast math, which cannot be turned off for part of a program; otherwise
// builds decompress different values.
#if FLT_EVAL_METHOD > 0
#error "SZ3 does not support x87 floating point; on 32-bit x86 build with -msse2 -mfpmath=sse"
#endif
#if defined(__FAST_MATH__) && !defined(__clang__)
#error "SZ3 does not support GCC's -ffast-math or -Ofast"
#endif
}  // namespace SZ3

#endif
