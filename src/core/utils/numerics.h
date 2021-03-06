#ifndef SRC_NUMERICS_H_
#define SRC_NUMERICS_H_

#include <cstdlib>
#include <fstream>

#include <cstddef>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <cmath>

#include "Exception.h"

#ifdef __SSE3__
#define SSE_ENABLED
#define SIMD_ENABLED
#endif

#ifdef __AVX__
#undef SSE_ENABLED
#define AVX_ENABLED
#define SIMD_ENABLED
#endif

#ifdef SINGLE_PRECISION
    typedef float RealNumber;
#else
    typedef double RealNumber;
#endif
    
#ifdef SSE_ENABLED

#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include "sse_mathfun.h"

#ifdef SINGLE_PRECISION
    typedef __m128 SIMDRegister;
#define       SIMD_LOAD1(A)    _mm_load1_ps(A)
#define       SIMD_SET1(A)     _mm_set1_ps(A)
#define       SIMD_MUL(A,B)    _mm_mul_ps(A,B)
#define       SIMD_ADD(A,B)    _mm_add_ps(A,B)
#define       SIMD_MAX(A,B)    _mm_max_ps(A,B)
#define       SIMD_DIV(A,B)    _mm_div_ps(A,B)
#define       SIMD_LOG(A)      log_ps(A)
#define       SIMD_CMPNEQ(A,B) _mm_cmpneq_ps(A,B)
#define       SIMD_AND(A,B)    _mm_and_ps(A,B)
#else
    typedef __m128d SIMDRegister;
#define       SIMD_LOAD1(A)    _mm_load1_pd(A)
#define       SIMD_SET1(A)     _mm_set1_pd(A)
#define       SIMD_MUL(A,B)    _mm_mul_pd(A,B)
#define       SIMD_ADD(A,B)    _mm_add_pd(A,B)
#define       SIMD_MAX(A,B)    _mm_max_pd(A,B)
#define       SIMD_DIV(A,B)    _mm_div_pd(A,B)
#define       SIMD_CMPNEQ(A,B) _mm_cmpneq_pd(A,B)
#define       SIMD_AND(A,B)    _mm_and_pd(A,B)
#endif

#elif defined AVX_ENABLED
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <immintrin.h>
#include "avx_mathfun.h"

#ifdef SINGLE_PRECISION
    typedef __m256 SIMDRegister;
#define     SIMD_LOAD1(A)      _mm256_broadcast_ss(A)
#define     SIMD_SET1(A)       _mm256_set1_ps(A)
#define     SIMD_MUL(A,B)      _mm256_mul_ps(A,B)
#define     SIMD_ADD(A,B)      _mm256_add_ps(A,B)
#define     SIMD_MAX(A,B)      _mm256_max_ps(A,B)
#define     SIMD_DIV(A,B)      _mm256_div_ps(A,B)
#define     SIMD_LOG(A)        log256_ps(A)
#define     SIMD_CMPNEQ(A,B,C) _mm256_cmp_ps(A,B,0x14) // 0x14 neq
#define     SIMD_AND(A,B)      _mm256_and_ps(A,B)
#else
    typedef __m256d SIMDRegister;
#define     SIMD_LOAD1(A)      _mm256_broadcast_sd(A)
#define     SIMD_SET1(A)       _mm256_set1_pd(A)
#define     SIMD_MUL(A,B)      _mm256_mul_pd(A,B)
#define     SIMD_ADD(A,B)      _mm256_add_pd(A,B)
#define     SIMD_MAX(A,B)      _mm256_max_pd(A,B)
#define     SIMD_DIV(A,B)      _mm256_div_pd(A,B)
#define     SIMD_CMPNEQ(A,B,C) _mm256_cmp_pd(A,B,0x14)
#define     SIMD_AND(A,B)      _mm256_and_pd(A,B)
#endif

#endif
    
#ifdef SIMD_ENABLED
    
const size_t REALS_PER_SIMD_REGISTER = sizeof(SIMDRegister)/sizeof(RealNumber);

#ifndef SINGLE_PRECISION
inline SIMDRegister SIMD_LOG(SIMDRegister& inreg)
{
    SIMDRegister outreg;
    RealNumber* in  = (RealNumber*)&inreg;
    RealNumber* out = (RealNumber*)&outreg;
    for(size_t i = 0; i < REALS_PER_SIMD_REGISTER; i++)
    {
        out[i] = log(in[i]);
    }

    return outreg;
}
#endif
    
/**
 * Allocator for aligned data.
 *
 * Modified from the Mallocator from Stephan T. Lavavej.
 * <http://blogs.msdn.com/b/vcblog/archive/2008/08/28/the-mallocator.aspx>
 */
template <typename T, std::size_t Alignment>
class aligned_allocator
{
    public:
 
        // The following will be the same for virtually all allocators.
        typedef T * pointer;
        typedef const T * const_pointer;
        typedef T& reference;
        typedef const T& const_reference;
        typedef T value_type;
        typedef std::size_t size_type;
        typedef ptrdiff_t difference_type;
 
        T * address(T& r) const
        {
            return &r;
        }
 
        const T * address(const T& s) const
        {
            return &s;
        }
 
        std::size_t max_size() const
        {
            // The following has been carefully written to be independent of
            // the definition of size_t and to avoid signed/unsigned warnings.
            return (static_cast<std::size_t>(0) - static_cast<std::size_t>(1)) / sizeof(T);
        }
 
 
        // The following must be the same for all allocators.
        template <typename U>
        struct rebind
        {
            typedef aligned_allocator<U, Alignment> other;
        } ;
 
        bool operator!=(const aligned_allocator& other) const
        {
            return !(*this == other);
        }
 
        void construct(T * const p, const T& t) const
        {
            void * const pv = static_cast<void *>(p);
 
            new (pv) T(t);
        }
 
        void destroy(T * const p) const
        {
            p->~T();
        }
 
        // Returns true if and only if storage allocated from *this
        // can be deallocated from other, and vice versa.
        // Always returns true for stateless allocators.
        bool operator==(const aligned_allocator& other) const
        {
            return true;
        }
 
 
        // Default constructor, copy constructor, rebinding constructor, and destructor.
        // Empty for stateless allocators.
        aligned_allocator() { }
 
        aligned_allocator(const aligned_allocator&) { }
 
        template <typename U> aligned_allocator(const aligned_allocator<U, Alignment>&) { }
 
        ~aligned_allocator() { }
 
 
        // The following will be different for each allocator.
        T * allocate(const std::size_t n) const
        {
            // The return value of allocate(0) is unspecified.
            // Mallocator returns NULL in order to avoid depending
            // on malloc(0)'s implementation-defined behavior
            // (the implementation can define malloc(0) to return NULL,
            // in which case the bad_alloc check below would fire).
            // All allocators can return NULL in this case.
            if (n == 0) {
                return NULL;
            }
 
            // All allocators should contain an integer overflow check.
            // The Standardization Committee recommends that std::length_error
            // be thrown in the case of integer overflow.
            if (n > max_size())
            {
                throw std::length_error("aligned_allocator<T>::allocate() - Integer overflow.");
            }
 
            // Mallocator wraps malloc().
            void * const pv = _mm_malloc(n * sizeof(T), Alignment);
 
            // Allocators should throw std::bad_alloc in the case of memory allocation failure.
            if (pv == NULL)
            {
                throw std::bad_alloc();
            }
 
            return static_cast<T *>(pv);
        }
 
        void deallocate(T * const p, const std::size_t n) const
        {
            _mm_free(p);
        }
 
 
        // The following will be the same for all allocators that ignore hints.
        template <typename U>
        T * allocate(const std::size_t n, const U * /* const hint */) const
        {
            return allocate(n);
        }
 
 
        // Allocators are not required to be assignable, so
        // all allocators should have a private unimplemented
        // assignment operator. Note that this will trigger the
        // off-by-default (enabled under /Wall) warning C4626
        // "assignment operator could not be generated because a
        // base class assignment operator is inaccessible" within
        // the STL headers, but that warning is useless.
    private:
        aligned_allocator& operator=(const aligned_allocator&);
};

    typedef std::vector<RealNumber, aligned_allocator<RealNumber, sizeof(SIMDRegister)> > RealVector;
#else
    typedef std::vector<RealNumber> RealVector;
#endif

#endif /* SRC_NUMERICS_H_ */
