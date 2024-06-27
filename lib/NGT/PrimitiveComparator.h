//
// Copyright (C) 2015 Yahoo Japan Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#pragma once

#include	"NGT/defines.h"

#if defined(NGT_NO_AVX)
#warning "*** SIMD is *NOT* available! ***"
#else
#include	<immintrin.h>
#endif

namespace NGT {

  class MemoryCache {
  public:
    inline static void prefetch(unsigned char *ptr, const size_t byteSizeOfObject) {
#if !defined(NGT_NO_AVX)
      switch((byteSizeOfObject - 1) >> 6) {
      default:
      case 28: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 27: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 26: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 25: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 24: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 23: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 22: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 21: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 20: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 19: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 18: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 17: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 16: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 15: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 14: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 13: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 12: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 11: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 10: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 9: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 8: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 7: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 6: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 5: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 4: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 3: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 2: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 1: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
      case 0: _mm_prefetch(ptr, _MM_HINT_T0); ptr += 64;
	break;
      }
#endif
    }
    inline static void *alignedAlloc(const size_t allocSize) {
#ifdef NGT_NO_AVX
      return new uint8_t[allocSize];
#else
#if defined(NGT_AVX512)
      size_t alignment = 64;
      uint64_t mask = 0xFFFFFFFFFFFFFFC0;
#elif defined(NGT_AVX2)
      size_t alignment = 32;
      uint64_t mask = 0xFFFFFFFFFFFFFFE0;
#else
      size_t alignment = 16;
      uint64_t mask = 0xFFFFFFFFFFFFFFF0;
#endif
      uint8_t *p = new uint8_t[allocSize + alignment];
      uint8_t *ptr = p + alignment;
      ptr = reinterpret_cast<uint8_t*>((reinterpret_cast<uint64_t>(ptr) & mask));
      *p++ = 0xAB;
      while (p != ptr) *p++ = 0xCD;
      return ptr;
#endif
    }
    inline static void alignedFree(void *ptr) {
#ifdef NGT_NO_AVX
      delete[] static_cast<uint8_t*>(ptr);
#else
      uint8_t *p = static_cast<uint8_t*>(ptr);
      p--;
      while (*p == 0xCD) p--;
      if (*p != 0xAB) {
	NGTThrowException("MemoryCache::alignedFree: Fatal Error! Cannot find allocated address.");
      }
      delete[] p;
#endif
    }
  };

  class PrimitiveComparator {
  public:

    static double absolute(double v) { return fabs(v); }
    static int absolute(int v) { return abs(v); }

#if defined(NGT_NO_AVX)
    template <typename OBJECT_TYPE, typename COMPARE_TYPE>
    inline static double compareL2(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      const OBJECT_TYPE *last = a + size;
      const OBJECT_TYPE *lastgroup = last - 3;
      COMPARE_TYPE diff0, diff1, diff2, diff3;
      double d = 0.0;
      while (a < lastgroup) {
	diff0 = static_cast<COMPARE_TYPE>(a[0] - b[0]);
	diff1 = static_cast<COMPARE_TYPE>(a[1] - b[1]);
	diff2 = static_cast<COMPARE_TYPE>(a[2] - b[2]);
	diff3 = static_cast<COMPARE_TYPE>(a[3] - b[3]);
	d += diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
	a += 4;
	b += 4;
      }
      while (a < last) {
	diff0 = static_cast<COMPARE_TYPE>(*a++ - *b++);
	d += diff0 * diff0;
      }
      return sqrt(static_cast<double>(d));
    }

    inline static double compareL2(const uint8_t *a, const uint8_t *b, size_t size) {
      return compareL2<uint8_t, int>(a, b, size);
    }

    inline static double compareL2(const float *a, const float *b, size_t size) {
      return compareL2<float, double>(a, b, size);
    }
#ifdef NGT_HALF_FLOAT
    inline static double compareL2(const float16 *a, const float16 *b, size_t size) {
      return compareL2<float16, double>(a, b, size);
    }
#endif
#ifdef NGT_BFLOAT
    inline static double compareL2(const bfloat16 *a, const bfloat16 *b, size_t size) {
      return compareL2<bfloat16, float>(a, b, size);
    }
#endif
#else
    inline static double compareL2(const float *a, const float *b, size_t size) {
      const float *last = a + size;
#if defined(NGT_AVX512)
      __m512 sum512 = _mm512_setzero_ps();
      while (a < last) {
	__m512 v = _mm512_sub_ps(_mm512_loadu_ps(a), _mm512_loadu_ps(b));
	sum512 = _mm512_add_ps(sum512, _mm512_mul_ps(v, v));
	a += 16;
	b += 16;
      }

      __m256 sum256 = _mm256_add_ps(_mm512_extractf32x8_ps(sum512, 0), _mm512_extractf32x8_ps(sum512, 1));
      __m128 sum128 = _mm_add_ps(_mm256_extractf128_ps(sum256, 0), _mm256_extractf128_ps(sum256, 1));
#elif defined(NGT_AVX2)
      __m256 sum256 = _mm256_setzero_ps();
      __m256 v;
      while (a < last) {
	v = _mm256_sub_ps(_mm256_loadu_ps(a), _mm256_loadu_ps(b));
	sum256 = _mm256_add_ps(sum256, _mm256_mul_ps(v, v));
	a += 8;
	b += 8;
	v = _mm256_sub_ps(_mm256_loadu_ps(a), _mm256_loadu_ps(b));
	sum256 = _mm256_add_ps(sum256, _mm256_mul_ps(v, v));
	a += 8;
	b += 8;
      }
      __m128 sum128 = _mm_add_ps(_mm256_extractf128_ps(sum256, 0), _mm256_extractf128_ps(sum256, 1));
#else
      __m128 sum128 = _mm_setzero_ps();
      __m128 v;
      while (a < last) {
	v = _mm_sub_ps(_mm_loadu_ps(a), _mm_loadu_ps(b));
	sum128 = _mm_add_ps(sum128, _mm_mul_ps(v, v));
        a += 4;
        b += 4;
	v = _mm_sub_ps(_mm_loadu_ps(a), _mm_loadu_ps(b));
	sum128 = _mm_add_ps(sum128, _mm_mul_ps(v, v));
        a += 4;
        b += 4;
	v = _mm_sub_ps(_mm_loadu_ps(a), _mm_loadu_ps(b));
	sum128 = _mm_add_ps(sum128, _mm_mul_ps(v, v));
        a += 4;
        b += 4;
	v = _mm_sub_ps(_mm_loadu_ps(a), _mm_loadu_ps(b));
	sum128 = _mm_add_ps(sum128, _mm_mul_ps(v, v));
        a += 4;
        b += 4;
      }
#endif

      __attribute__((aligned(32))) float f[4];
      _mm_store_ps(f, sum128);

      double s = f[0] + f[1] + f[2] + f[3];
      return sqrt(s);
    }

#ifdef NGT_HALF_FLOAT
    inline static double compareL2(const float16 *a, const float16 *b, size_t size) {
      const float16 *last = a + size;
#if defined(NGT_AVX512)
      __m512 sum512 = _mm512_setzero_ps();
      while (a < last) {
	__m512 v = _mm512_sub_ps(_mm512_cvtph_ps(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(a))),
				 _mm512_cvtph_ps(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(b))));
	sum512 = _mm512_add_ps(sum512, _mm512_mul_ps(v, v));
	a += 16;
	b += 16;
      }

      __m256 sum256 = _mm256_add_ps(_mm512_extractf32x8_ps(sum512, 0), _mm512_extractf32x8_ps(sum512, 1));
      __m128 sum128 = _mm_add_ps(_mm256_extractf128_ps(sum256, 0), _mm256_extractf128_ps(sum256, 1));
#elif defined(NGT_AVX2)
      __m256 sum256 = _mm256_setzero_ps();
      __m256 v;
      while (a < last) {
	v = _mm256_sub_ps(_mm256_cvtph_ps(_mm_loadu_si128(reinterpret_cast<const __m128i*>(a))),
			  _mm256_cvtph_ps(_mm_loadu_si128(reinterpret_cast<const __m128i*>(b))));
	sum256 = _mm256_add_ps(sum256, _mm256_mul_ps(v, v));
	a += 8;
	b += 8;
	v = _mm256_sub_ps(_mm256_cvtph_ps(_mm_loadu_si128(reinterpret_cast<const __m128i*>(a))),
			  _mm256_cvtph_ps(_mm_loadu_si128(reinterpret_cast<const __m128i*>(b))));
	sum256 = _mm256_add_ps(sum256, _mm256_mul_ps(v, v));
	a += 8;
	b += 8;
      }
      __m128 sum128 = _mm_add_ps(_mm256_extractf128_ps(sum256, 0), _mm256_extractf128_ps(sum256, 1));
#else
      __m128 sum128 = _mm_setzero_ps();
      __m128 v;
      while (a < last) {
	__m128i va = _mm_load_si128(reinterpret_cast<const __m128i*>(a));
	__m128i vb = _mm_load_si128(reinterpret_cast<const __m128i*>(b));
	v = _mm_sub_ps(_mm_cvtph_ps(va), _mm_cvtph_ps(vb));
	sum128 = _mm_add_ps(sum128, _mm_mul_ps(v, v));
	va = _mm_srli_si128(va, 8);
	vb = _mm_srli_si128(vb, 8);
	v = _mm_sub_ps(_mm_cvtph_ps(va), _mm_cvtph_ps(vb));
	sum128 = _mm_add_ps(sum128, _mm_mul_ps(v, v));
        a += 8;
        b += 8;
      }
#endif
      __m128 tmp = _mm_hadd_ps(sum128, _mm_set1_ps(0));
      double s = _mm_cvtss_f32(_mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(0, 0, 0, 0))) +
	         _mm_cvtss_f32(_mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(0, 0, 0, 1)));
      return sqrt(s);
      //return s;
    }
#endif

#ifdef NGT_BFLOAT
    inline static double compareL2(const bfloat16 *a, const bfloat16 *b, size_t size) {
      const bfloat16 *last = a + size;
      __m512 sum512 = _mm512_setzero_ps();
      while (a < last) {
	__m512i av = _mm512_slli_epi32(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(a))), 16);
	__m512i bv = _mm512_slli_epi32(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(b))), 16);
	__m512 sub = _mm512_sub_ps(reinterpret_cast<__m512>(av), reinterpret_cast<__m512>(bv));
	sum512 = _mm512_fmadd_ps(sub, sub, sum512);
	a += 16;
	b += 16;
      }
      __m256 sum256 = _mm256_add_ps(_mm512_extractf32x8_ps(sum512, 0), _mm512_extractf32x8_ps(sum512, 1));
      __m128 sum128 = _mm_add_ps(_mm256_extractf128_ps(sum256, 0), _mm256_extractf128_ps(sum256, 1));
      __m128 tmp = _mm_hadd_ps(sum128, _mm_set1_ps(0));
      double d = _mm_cvtss_f32(_mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(0, 0, 0, 0))) +
	         _mm_cvtss_f32(_mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(0, 0, 0, 1)));
      //return sqrt(d);
      return d;
    }
#endif

    inline static double compareL2(const unsigned char *a, const unsigned char *b, size_t size) {
      __m128 sum = _mm_setzero_ps();
      const unsigned char *last = a + size;
      const unsigned char *lastgroup = last - 7;
      const __m128i zero = _mm_setzero_si128();
      while (a < lastgroup) {
	//__m128i x1 = _mm_cvtepu8_epi16(*reinterpret_cast<__m128i const*>(a));
	__m128i x1 = _mm_cvtepu8_epi16(_mm_loadu_si128((__m128i const*)a));
	//__m128i x2 = _mm_cvtepu8_epi16(*reinterpret_cast<__m128i const*>(b));
	__m128i x2 = _mm_cvtepu8_epi16(_mm_loadu_si128((__m128i const*)b));
	x1 = _mm_subs_epi16(x1, x2);
	__m128i v = _mm_mullo_epi16(x1, x1);
	sum = _mm_add_ps(sum, _mm_cvtepi32_ps(_mm_unpacklo_epi16(v, zero)));
	sum = _mm_add_ps(sum, _mm_cvtepi32_ps(_mm_unpackhi_epi16(v, zero)));
	a += 8;
	b += 8;
      }
      __attribute__((aligned(32))) float f[4];
      _mm_store_ps(f, sum);
      double s = f[0] + f[1] + f[2] + f[3];
      while (a < last) {
	int d = (int)*a++ - (int)*b++;
	s += d * d;
      }
      return sqrt(s);
    }


#endif

    template <typename OBJECT_TYPE>
    inline static double compareNormalizedL2(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      double v = 2.0 - 2.0 * compareDotProduct(a, b, size);
      if (v < 0.0) {
	return 0.0;
      } else {
	return sqrt(v);
      }
    }


#if defined(NGT_NO_AVX)
    template <typename OBJECT_TYPE, typename COMPARE_TYPE>
    static double compareL1(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      const OBJECT_TYPE *last = a + size;
      const OBJECT_TYPE *lastgroup = last - 3;
      COMPARE_TYPE diff0, diff1, diff2, diff3;
      double d = 0.0;
      while (a < lastgroup) {
	diff0 = (COMPARE_TYPE)(a[0] - b[0]);
	diff1 = (COMPARE_TYPE)(a[1] - b[1]);
	diff2 = (COMPARE_TYPE)(a[2] - b[2]);
	diff3 = (COMPARE_TYPE)(a[3] - b[3]);
	d += absolute(diff0) + absolute(diff1) + absolute(diff2) + absolute(diff3);
	a += 4;
	b += 4;
      }
      while (a < last) {
	diff0 = (COMPARE_TYPE)*a++ - (COMPARE_TYPE)*b++;
	d += absolute(diff0);
      }
      return d;
    }

    inline static double compareL1(const uint8_t *a, const uint8_t *b, size_t size) {
      return compareL1<uint8_t, int>(a, b, size);
    }

    inline static double compareL1(const float *a, const float *b, size_t size) {
      return compareL1<float, double>(a, b, size);
    }
#ifdef NGT_HALF_FLOAT
    inline static double compareL1(const float16 *a, const float16 *b, size_t size) {
      return compareL1<float16, double>(a, b, size);
    }
#endif
#ifdef NGT_BFLOAT
    inline static double compareL1(const bfloat16 *a, const bfloat16 *b, size_t size) {
      return compareL1<bfloat16, float>(a, b, size);
    }
#endif
#else
    inline static double compareL1(const float *a, const float *b, size_t size) {
      __m256 sum = _mm256_setzero_ps();
      const float *last = a + size;
      const float *lastgroup = last - 7;
      while (a < lastgroup) {
	__m256 x1 = _mm256_sub_ps(_mm256_loadu_ps(a), _mm256_loadu_ps(b));
	const __m256 mask = _mm256_set1_ps(-0.0f);
	__m256 v = _mm256_andnot_ps(mask, x1);
	sum = _mm256_add_ps(sum, v);
	a += 8;
	b += 8;
      }
      __attribute__((aligned(32))) float f[8];
      _mm256_store_ps(f, sum);
      double s = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7];
      while (a < last) {
	double d = fabs(*a++ - *b++);
	s += d;
      }
      return s;
    }
#ifdef NGT_HALF_FLOAT
    inline static double compareL1(const float16 *a, const float16 *b, size_t size) {
      __m256 sum = _mm256_setzero_ps();
      const float16 *last = a + size;
      const float16 *lastgroup = last - 7;
      while (a < lastgroup) {
	__m256 x1 = _mm256_sub_ps(_mm256_cvtph_ps(_mm_loadu_si128(reinterpret_cast<const __m128i*>(a))),
				  _mm256_cvtph_ps(_mm_loadu_si128(reinterpret_cast<const __m128i*>(b))));
	const __m256 mask = _mm256_set1_ps(-0.0f);
	__m256 v = _mm256_andnot_ps(mask, x1);
	sum = _mm256_add_ps(sum, v);
	a += 8;
	b += 8;
      }
      __attribute__((aligned(32))) float f[8];
      _mm256_store_ps(f, sum);
      double s = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7];
      while (a < last) {
	double d = fabs(*a++ - *b++);
	s += d;
      }
      return s;
    }
#endif
#ifdef NGT_BFLOAT
    inline static double compareL1(const bfloat16 *a, const bfloat16 *b, size_t size) {
      abort();
    }
#endif
    inline static double compareL1(const unsigned char *a, const unsigned char *b, size_t size) {
      __m128 sum = _mm_setzero_ps();
      const unsigned char *last = a + size;
      const unsigned char *lastgroup = last - 7;
      const __m128i zero = _mm_setzero_si128();
      while (a < lastgroup) {
	__m128i x1 = _mm_cvtepu8_epi16(*reinterpret_cast<__m128i const*>(a));
	__m128i x2 = _mm_cvtepu8_epi16(*reinterpret_cast<__m128i const*>(b));
	x1 = _mm_subs_epi16(x1, x2);
	x1 = _mm_sign_epi16(x1, x1);
	sum = _mm_add_ps(sum, _mm_cvtepi32_ps(_mm_unpacklo_epi16(x1, zero)));
	sum = _mm_add_ps(sum, _mm_cvtepi32_ps(_mm_unpackhi_epi16(x1, zero)));
	a += 8;
	b += 8;
      }
      __attribute__((aligned(32))) float f[4];
      _mm_store_ps(f, sum);
      double s = f[0] + f[1] + f[2] + f[3];
      while (a < last) {
	double d = fabs(static_cast<double>(*a++) - static_cast<double>(*b++));
	s += d;
      }
      return s;
    }
#endif

#if defined(NGT_NO_AVX) || !defined(__POPCNT__)
    inline static double popCount(uint32_t x) {
      x = (x & 0x55555555) + (x >> 1 & 0x55555555);
      x = (x & 0x33333333) + (x >> 2 & 0x33333333);
      x = (x & 0x0F0F0F0F) + (x >> 4 & 0x0F0F0F0F);
      x = (x & 0x00FF00FF) + (x >> 8 & 0x00FF00FF);
      x = (x & 0x0000FFFF) + (x >> 16 & 0x0000FFFF);
      return x;
    }

    template <typename OBJECT_TYPE>
    inline static double compareHammingDistance(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      const uint32_t *last = reinterpret_cast<const uint32_t*>(a + size);
      const uint32_t *uinta = reinterpret_cast<const uint32_t*>(a);
      const uint32_t *uintb = reinterpret_cast<const uint32_t*>(b);
      size_t count = 0;
      while( uinta < last ){
	count += popCount(*uinta++ ^ *uintb++);
      }

      return static_cast<double>(count);
    }
#else
template <typename OBJECT_TYPE>
inline static double compareHammingDistance(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size)
{
    size_t count = 0;
    const OBJECT_TYPE *last = a + size;
#if defined(__AVX512F__)
    while (a + 64 <= last)
    {
        __m512i vxor = _mm512_xor_si512(_mm512_loadu_si512(reinterpret_cast<const __m512i *>(a)), _mm512_loadu_si512(reinterpret_cast<const __m512i *>(b)));
        count += _mm_popcnt_u64(_mm512_extract_epi64(vxor, 0));
        count += _mm_popcnt_u64(_mm512_extract_epi64(vxor, 1));
        count += _mm_popcnt_u64(_mm512_extract_epi64(vxor, 2));
        count += _mm_popcnt_u64(_mm512_extract_epi64(vxor, 3));
        count += _mm_popcnt_u64(_mm512_extract_epi64(vxor, 4));
        count += _mm_popcnt_u64(_mm512_extract_epi64(vxor, 5));
        count += _mm_popcnt_u64(_mm512_extract_epi64(vxor, 6));
        count += _mm_popcnt_u64(_mm512_extract_epi64(vxor, 7));
        a += 64;
        b += 64;
    }
#endif

#if defined(__AVX512F__) || defined(__AVX2__)
    while (a + 32 <= last)
    {
        __m256i vxor = _mm256_xor_si256(_mm256_loadu_si256(reinterpret_cast<const __m256i *>(a)), _mm256_loadu_si256(reinterpret_cast<const __m256i *>(b)));
        count += _mm_popcnt_u64(_mm256_extract_epi64(vxor, 0));
        count += _mm_popcnt_u64(_mm256_extract_epi64(vxor, 1));
        count += _mm_popcnt_u64(_mm256_extract_epi64(vxor, 2));
        count += _mm_popcnt_u64(_mm256_extract_epi64(vxor, 3));
        a += 32;
        b += 32;
    }
#endif
    while (a < last)
    {
        __m128i vxor = _mm_xor_si128(_mm_loadu_si128(reinterpret_cast<const __m128i *>(a)), _mm_loadu_si128(reinterpret_cast<const __m128i *>(b)));
        count += _mm_popcnt_u32(_mm_extract_epi32(vxor, 0));
        count += _mm_popcnt_u32(_mm_extract_epi32(vxor, 1));
        count += _mm_popcnt_u32(_mm_extract_epi32(vxor, 2));
        count += _mm_popcnt_u32(_mm_extract_epi32(vxor, 3));
        a += 16;
        b += 16;
    }
    return static_cast<double>(count);
}

#endif

#if defined(NGT_NO_AVX) || !defined(__POPCNT__)
    template <typename OBJECT_TYPE>
      inline static double compareJaccardDistance(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      const uint32_t *last = reinterpret_cast<const uint32_t*>(a + size);

      const uint32_t *uinta = reinterpret_cast<const uint32_t*>(a);
      const uint32_t *uintb = reinterpret_cast<const uint32_t*>(b);
      size_t count = 0;
      size_t countDe = 0;
      while( uinta < last ){
	count   += popCount(*uinta   & *uintb);
	countDe += popCount(*uinta++ | *uintb++);
	count   += popCount(*uinta   & *uintb);
	countDe += popCount(*uinta++ | *uintb++);
      }

      return 1.0 - static_cast<double>(count) / static_cast<double>(countDe);
    }
#else
    template <typename OBJECT_TYPE>
      inline static double compareJaccardDistance(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      const uint64_t *last = reinterpret_cast<const uint64_t*>(a + size);

      const uint64_t *uinta = reinterpret_cast<const uint64_t*>(a);
      const uint64_t *uintb = reinterpret_cast<const uint64_t*>(b);
      size_t count = 0;
      size_t countDe = 0;
      while( uinta < last ){
	count   += _mm_popcnt_u64(*uinta   & *uintb);
	countDe += _mm_popcnt_u64(*uinta++ | *uintb++);
	count   += _mm_popcnt_u64(*uinta   & *uintb);
	countDe += _mm_popcnt_u64(*uinta++ | *uintb++);
      }

      return 1.0 - static_cast<double>(count) / static_cast<double>(countDe);
    }
#endif

    inline static double compareSparseJaccardDistance(const unsigned char *a, const unsigned char *b, size_t size) {
      std::cerr << "compareSparseJaccardDistance: Not implemented." << std::endl;
      abort();
    }

#ifdef NGT_HALF_FLOAT
    inline static double compareSparseJaccardDistance(const float16 *a, const float16 *b, size_t size) {
      std::cerr << "compareSparseJaccardDistance: Not implemented." << std::endl;
      abort();
    }
#endif

#ifdef NGT_BFLOAT
    inline static double compareSparseJaccardDistance(const bfloat16 *a, const bfloat16 *b, size_t size) {
      std::cerr << "compareSparseJaccardDistance: Not implemented." << std::endl;
      abort();
    }
#endif

    inline static double compareSparseJaccardDistance(const float *a, const float *b, size_t size) {
      size_t loca = 0;
      size_t locb = 0;
      const uint32_t *ai = reinterpret_cast<const uint32_t*>(a);
      const uint32_t *bi = reinterpret_cast<const uint32_t*>(b);
      size_t count = 0;
      while (locb < size && ai[loca] != 0 && bi[loca] != 0) {
	int64_t sub = static_cast<int64_t>(ai[loca]) - static_cast<int64_t>(bi[locb]);
	count += sub == 0;
	loca += sub <= 0;
	locb += sub >= 0;
      }
      while (ai[loca] != 0) {
	loca++;
      }
      while (locb < size && bi[locb] != 0) {
	locb++;
      }
      return 1.0 - static_cast<double>(count) / static_cast<double>(loca + locb - count);
    }

#if defined(NGT_NO_AVX)
   template <typename OBJECT_TYPE>
    inline static double compareDotProduct(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      double sum = 0.0;
      for (size_t loc = 0; loc < size; loc++) {
	sum += static_cast<float>(a[loc]) * static_cast<float>(b[loc]);
      }
      return sum;
    }

    template <typename OBJECT_TYPE>
    inline static double compareCosine(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      double normA = 0.0;
      double normB = 0.0;
      double sum = 0.0;
      for (size_t loc = 0; loc < size; loc++) {
	normA += static_cast<double>(a[loc]) * static_cast<double>(a[loc]);
	normB += static_cast<double>(b[loc]) * static_cast<double>(b[loc]);
	sum += static_cast<double>(a[loc]) * static_cast<double>(b[loc]);
      }

      double cosine = sum / sqrt(normA * normB);

      return cosine;
    }
#else
    inline static double compareDotProduct(const float *a, const float *b, size_t size) {
      const float *last = a + size;
#if defined(NGT_AVX512)
      __m512 sum512 = _mm512_setzero_ps();
      while (a < last) {
	sum512 = _mm512_add_ps(sum512, _mm512_mul_ps(_mm512_loadu_ps(a), _mm512_loadu_ps(b)));
	a += 16;
	b += 16;
      }
      __m256 sum256 = _mm256_add_ps(_mm512_extractf32x8_ps(sum512, 0), _mm512_extractf32x8_ps(sum512, 1));
      __m128 sum128 = _mm_add_ps(_mm256_extractf128_ps(sum256, 0), _mm256_extractf128_ps(sum256, 1));
#elif defined(NGT_AVX2)
      __m256 sum256 = _mm256_setzero_ps();
      while (a < last) {
	sum256 = _mm256_add_ps(sum256, _mm256_mul_ps(_mm256_loadu_ps(a), _mm256_loadu_ps(b)));
	a += 8;
	b += 8;
      }
      __m128 sum128 = _mm_add_ps(_mm256_extractf128_ps(sum256, 0), _mm256_extractf128_ps(sum256, 1));
#else
      __m128 sum128 = _mm_setzero_ps();
      while (a < last) {
	sum128 = _mm_add_ps(sum128, _mm_mul_ps(_mm_loadu_ps(a), _mm_loadu_ps(b)));
	a += 4;
	b += 4;
      }
#endif
      __attribute__((aligned(32))) float f[4];
      _mm_store_ps(f, sum128);
      double s = static_cast<double>(f[0]) + static_cast<double>(f[1]) + static_cast<double>(f[2]) + static_cast<double>(f[3]);
      return s;
    }

#ifdef NGT_HALF_FLOAT
    inline static double compareDotProduct(const float16 *a, const float16 *b, size_t size) {
      const float16 *last = a + size;
#if defined(NGT_AVX512)
      __m512 sum512 = _mm512_setzero_ps();
      while (a < last) {
	sum512 = _mm512_add_ps(sum512, _mm512_mul_ps(_mm512_cvtph_ps(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(a))),
						     _mm512_cvtph_ps(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(b)))));

	a += 16;
	b += 16;
      }
      __m256 sum256 = _mm256_add_ps(_mm512_extractf32x8_ps(sum512, 0), _mm512_extractf32x8_ps(sum512, 1));
      __m128 sum128 = _mm_add_ps(_mm256_extractf128_ps(sum256, 0), _mm256_extractf128_ps(sum256, 1));
#elif defined(NGT_AVX2)
      __m256 sum256 = _mm256_setzero_ps();
      while (a < last) {
	sum256 = _mm256_add_ps(sum256, _mm256_mul_ps(_mm256_cvtph_ps(_mm_loadu_si128(reinterpret_cast<const __m128i*>(a))),
						     _mm256_cvtph_ps(_mm_loadu_si128(reinterpret_cast<const __m128i*>(b)))));
	a += 8;
	b += 8;
      }
      __m128 sum128 = _mm_add_ps(_mm256_extractf128_ps(sum256, 0), _mm256_extractf128_ps(sum256, 1));
#else
      __m128 sum128 = _mm_setzero_ps();
      while (a < last) {
	__m128i va = _mm_load_si128(reinterpret_cast<const __m128i*>(a));
	__m128i vb = _mm_load_si128(reinterpret_cast<const __m128i*>(b));
	sum128 = _mm_add_ps(sum128, _mm_mul_ps(_mm_cvtph_ps(va), _mm_cvtph_ps(vb)));
	va = _mm_srli_si128(va, 8);
	vb = _mm_srli_si128(vb, 8);
	sum128 = _mm_add_ps(sum128, _mm_mul_ps(_mm_cvtph_ps(va), _mm_cvtph_ps(vb)));
	a += 8;
	b += 8;
      }
#endif
      __attribute__((aligned(32))) float f[4];
      _mm_store_ps(f, sum128);
      double s = static_cast<double>(f[0]) + static_cast<double>(f[1]) + static_cast<double>(f[2]) + static_cast<double>(f[3]);
      return s;
    }
#endif

#ifdef NGT_BFLOAT
    inline static double compareDotProduct(const bfloat16 *a, const bfloat16 *b, size_t size) {
      abort();
    }
#endif

    inline static double compareDotProduct(const unsigned char *a, const unsigned char *b, size_t size) {
      double sum = 0.0;
      for (size_t loc = 0; loc < size; loc++) {
	sum += static_cast<double>(a[loc]) * static_cast<double>(b[loc]);
      }
      return sum;
    }
    inline static double compareCosine(const float *a, const float *b, size_t size) {

      const float *last = a + size;
#if defined(NGT_AVX512)
      __m512 normA = _mm512_setzero_ps();
      __m512 normB = _mm512_setzero_ps();
      __m512 sum = _mm512_setzero_ps();
      while (a < last) {
	__m512 am = _mm512_loadu_ps(a);
	__m512 bm = _mm512_loadu_ps(b);
	normA = _mm512_add_ps(normA, _mm512_mul_ps(am, am));
	normB = _mm512_add_ps(normB, _mm512_mul_ps(bm, bm));
	sum = _mm512_add_ps(sum, _mm512_mul_ps(am, bm));
	a += 16;
	b += 16;
      }
      __m256 am256 = _mm256_add_ps(_mm512_extractf32x8_ps(normA, 0), _mm512_extractf32x8_ps(normA, 1));
      __m256 bm256 = _mm256_add_ps(_mm512_extractf32x8_ps(normB, 0), _mm512_extractf32x8_ps(normB, 1));
      __m256 s256 = _mm256_add_ps(_mm512_extractf32x8_ps(sum, 0), _mm512_extractf32x8_ps(sum, 1));
      __m128 am128 = _mm_add_ps(_mm256_extractf128_ps(am256, 0), _mm256_extractf128_ps(am256, 1));
      __m128 bm128 = _mm_add_ps(_mm256_extractf128_ps(bm256, 0), _mm256_extractf128_ps(bm256, 1));
      __m128 s128 = _mm_add_ps(_mm256_extractf128_ps(s256, 0), _mm256_extractf128_ps(s256, 1));
#elif defined(NGT_AVX2)
      __m256 normA = _mm256_setzero_ps();
      __m256 normB = _mm256_setzero_ps();
      __m256 sum = _mm256_setzero_ps();
      __m256 am, bm;
      while (a < last) {
	am = _mm256_loadu_ps(a);
	bm = _mm256_loadu_ps(b);
	normA = _mm256_add_ps(normA, _mm256_mul_ps(am, am));
	normB = _mm256_add_ps(normB, _mm256_mul_ps(bm, bm));
	sum = _mm256_add_ps(sum, _mm256_mul_ps(am, bm));
	a += 8;
	b += 8;
      }
      __m128 am128 = _mm_add_ps(_mm256_extractf128_ps(normA, 0), _mm256_extractf128_ps(normA, 1));
      __m128 bm128 = _mm_add_ps(_mm256_extractf128_ps(normB, 0), _mm256_extractf128_ps(normB, 1));
      __m128 s128 = _mm_add_ps(_mm256_extractf128_ps(sum, 0), _mm256_extractf128_ps(sum, 1));
#else
      __m128 am128 = _mm_setzero_ps();
      __m128 bm128 = _mm_setzero_ps();
      __m128 s128 = _mm_setzero_ps();
      __m128 am, bm;
      while (a < last) {
	am = _mm_loadu_ps(a);
	bm = _mm_loadu_ps(b);
	am128 = _mm_add_ps(am128, _mm_mul_ps(am, am));
	bm128 = _mm_add_ps(bm128, _mm_mul_ps(bm, bm));
	s128 = _mm_add_ps(s128, _mm_mul_ps(am, bm));
	a += 4;
	b += 4;
      }

#endif

      __attribute__((aligned(32))) float f[4];
      _mm_store_ps(f, am128);
      double na = f[0] + f[1] + f[2] + f[3];
      _mm_store_ps(f, bm128);
      double nb = f[0] + f[1] + f[2] + f[3];
      _mm_store_ps(f, s128);
      double s = f[0] + f[1] + f[2] + f[3];

      double cosine = s / sqrt(na * nb);
      return cosine;
    }

#ifdef NGT_HALF_FLOAT
    inline static double compareCosine(const float16 *a, const float16 *b, size_t size) {

      const float16 *last = a + size;
#if defined(NGT_AVX512)
      __m512 normA = _mm512_setzero_ps();
      __m512 normB = _mm512_setzero_ps();
      __m512 sum = _mm512_setzero_ps();
      while (a < last) {
	__m512 am = _mm512_cvtph_ps(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(a)));
	__m512 bm = _mm512_cvtph_ps(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(b)));
	normA = _mm512_add_ps(normA, _mm512_mul_ps(am, am));
	normB = _mm512_add_ps(normB, _mm512_mul_ps(bm, bm));
	sum = _mm512_add_ps(sum, _mm512_mul_ps(am, bm));
	a += 16;
	b += 16;
      }
      __m256 am256 = _mm256_add_ps(_mm512_extractf32x8_ps(normA, 0), _mm512_extractf32x8_ps(normA, 1));
      __m256 bm256 = _mm256_add_ps(_mm512_extractf32x8_ps(normB, 0), _mm512_extractf32x8_ps(normB, 1));
      __m256 s256 = _mm256_add_ps(_mm512_extractf32x8_ps(sum, 0), _mm512_extractf32x8_ps(sum, 1));
      __m128 am128 = _mm_add_ps(_mm256_extractf128_ps(am256, 0), _mm256_extractf128_ps(am256, 1));
      __m128 bm128 = _mm_add_ps(_mm256_extractf128_ps(bm256, 0), _mm256_extractf128_ps(bm256, 1));
      __m128 s128 = _mm_add_ps(_mm256_extractf128_ps(s256, 0), _mm256_extractf128_ps(s256, 1));
#elif defined(NGT_AVX2)
      __m256 normA = _mm256_setzero_ps();
      __m256 normB = _mm256_setzero_ps();
      __m256 sum = _mm256_setzero_ps();
      __m256 am, bm;
      while (a < last) {
	am = _mm256_cvtph_ps(_mm_loadu_si128(reinterpret_cast<const __m128i*>(a)));
	bm = _mm256_cvtph_ps(_mm_loadu_si128(reinterpret_cast<const __m128i*>(b)));
	normA = _mm256_add_ps(normA, _mm256_mul_ps(am, am));
	normB = _mm256_add_ps(normB, _mm256_mul_ps(bm, bm));
	sum = _mm256_add_ps(sum, _mm256_mul_ps(am, bm));
	a += 8;
	b += 8;
      }
      __m128 am128 = _mm_add_ps(_mm256_extractf128_ps(normA, 0), _mm256_extractf128_ps(normA, 1));
      __m128 bm128 = _mm_add_ps(_mm256_extractf128_ps(normB, 0), _mm256_extractf128_ps(normB, 1));
      __m128 s128 = _mm_add_ps(_mm256_extractf128_ps(sum, 0), _mm256_extractf128_ps(sum, 1));
#else
      __m128 am128 = _mm_setzero_ps();
      __m128 bm128 = _mm_setzero_ps();
      __m128 s128 = _mm_setzero_ps();
      __m128 am, bm;
      while (a < last) {
	__m128i va = _mm_load_si128(reinterpret_cast<const __m128i*>(a));
	__m128i vb = _mm_load_si128(reinterpret_cast<const __m128i*>(b));
	am = _mm_cvtph_ps(va);
	bm = _mm_cvtph_ps(vb);
	am128 = _mm_add_ps(am128, _mm_mul_ps(am, am));
	bm128 = _mm_add_ps(bm128, _mm_mul_ps(bm, bm));
	s128 = _mm_add_ps(s128, _mm_mul_ps(am, bm));
	va = _mm_srli_si128(va, 8);
	vb = _mm_srli_si128(vb, 8);
	am = _mm_cvtph_ps(va);
	bm = _mm_cvtph_ps(vb);
	am128 = _mm_add_ps(am128, _mm_mul_ps(am, am));
	bm128 = _mm_add_ps(bm128, _mm_mul_ps(bm, bm));
	s128 = _mm_add_ps(s128, _mm_mul_ps(am, bm));
	a += 8;
	b += 8;
      }

#endif

      __attribute__((aligned(32))) float f[4];
      _mm_store_ps(f, am128);
      double na = f[0] + f[1] + f[2] + f[3];
      _mm_store_ps(f, bm128);
      double nb = f[0] + f[1] + f[2] + f[3];
      _mm_store_ps(f, s128);
      double s = f[0] + f[1] + f[2] + f[3];

      double cosine = s / sqrt(na * nb);
      return cosine;
    }
#endif

#ifdef NGT_BFLOAT
    inline static double compareCosine(const bfloat16 *a, const bfloat16 *b, size_t size) {
      abort();
    }
#endif

    inline static double compareCosine(const unsigned char *a, const unsigned char *b, size_t size) {
      double normA = 0.0;
      double normB = 0.0;
      double sum = 0.0;
      for (size_t loc = 0; loc < size; loc++) {
	normA += static_cast<double>(a[loc]) * static_cast<double>(a[loc]);
	normB += static_cast<double>(b[loc]) * static_cast<double>(b[loc]);
	sum += static_cast<double>(a[loc]) * static_cast<double>(b[loc]);
      }

      double cosine = sum / sqrt(normA * normB);

      return cosine;
    }


#endif    // #if defined(NGT_NO_AVX)

    template <typename OBJECT_TYPE>
    inline static double compareAngleDistance(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      double cosine = compareCosine(a, b, size);
      if (cosine >= 1.0) {
	return 0.0;
      } else if (cosine <= -1.0) {
	return acos(-1.0);
      } else {
	return acos(cosine);
      }
    }

    template <typename OBJECT_TYPE>
    inline static double compareNormalizedAngleDistance(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      double cosine = compareDotProduct(a, b, size);
      if (cosine >= 1.0) {
	return 0.0;
      } else if (cosine <= -1.0) {
	return acos(-1.0);
      } else {
	return acos(cosine);
      }
    }

    // added by Nyapicom
    template <typename OBJECT_TYPE>
    inline static double comparePoincareDistance(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      // Unlike the other distance functions, this is not optimized...
      double a2 = 0.0;
      double b2 = 0.0;
      double c2 = compareL2(a, b, size);
      for(size_t i = 0; i < size; i++){
	a2 += static_cast<double>(a[i]) * static_cast<double>(a[i]);
	b2 += static_cast<double>(b[i]) * static_cast<double>(b[i]);
      }
      return std::acosh(1 + 2.0 * c2*c2 / (1.0 - a2) / (1.0 - b2));
    }

    // added by Nyapicom
    template <typename OBJECT_TYPE>
    inline static double compareLorentzDistance(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      // Unlike the other distance functions, this is not optimized...
      double sum = static_cast<double>(a[0]) * static_cast<double>(b[0]);
      for(size_t i = 1; i < size; i++){
	sum -= static_cast<double>(a[i]) * static_cast<double>(b[i]);
      }
      return std::acosh(sum);
    }

    template <typename OBJECT_TYPE>
    inline static double compareCosineSimilarity(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      auto v = 1.0 - compareCosine(a, b, size);
      return v < 0.0 ? -v : v;
    }

    template <typename OBJECT_TYPE>
    inline static double compareNormalizedCosineSimilarity(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      auto v = 1.0 - compareDotProduct(a, b, size);
      return v < 0.0 ? -v : v;
    }

    class L1Uint8 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareL1((const uint8_t*)a, (const uint8_t*)b, size);
      }
    };

    class L2Uint8 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareL2((const uint8_t*)a, (const uint8_t*)b, size);
      }
    };

    class HammingUint8 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareHammingDistance((const uint8_t*)a, (const uint8_t*)b, size);
      }
    };

    class JaccardUint8 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareJaccardDistance((const uint8_t*)a, (const uint8_t*)b, size);
      }
    };

    class SparseJaccardFloat {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareSparseJaccardDistance((const float*)a, (const float*)b, size);
      }
    };

    class L2Float {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareL2((const float*)a, (const float*)b, size);
      }
    };

    class NormalizedL2Float {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareNormalizedL2((const float*)a, (const float*)b, size);
      }
    };

    class L1Float {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareL1((const float*)a, (const float*)b, size);
      }
    };

    class CosineSimilarityFloat {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareCosineSimilarity((const float*)a, (const float*)b, size);
      }
    };

    class NormalizedCosineSimilarityFloat {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareNormalizedCosineSimilarity((const float*)a, (const float*)b, size);
      }
    };

    class AngleFloat {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareAngleDistance((const float*)a, (const float*)b, size);
      }
    };

    class NormalizedAngleFloat {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareNormalizedAngleDistance((const float*)a, (const float*)b, size);
      }
    };

    // added by Nyapicom
    class PoincareFloat {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::comparePoincareDistance((const float*)a, (const float*)b, size);
      }
    };

    // added by Nyapicom
    class LorentzFloat {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareLorentzDistance((const float*)a, (const float*)b, size);
      }
    };

#ifdef NGT_HALF_FLOAT
    class SparseJaccardFloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareSparseJaccardDistance((const float16*)a, (const float16*)b, size);
      }
    };

    class L2Float16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareL2((const float16*)a, (const float16*)b, size);
      }
    };

    class NormalizedL2Float16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareNormalizedL2((const float16*)a, (const float16*)b, size);
      }
    };

    class L1Float16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareL1((const float16*)a, (const float16*)b, size);
      }
    };

    class CosineSimilarityFloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareCosineSimilarity((const float16*)a, (const float16*)b, size);
      }
    };

    class NormalizedCosineSimilarityFloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareNormalizedCosineSimilarity((const float16*)a, (const float16*)b, size);
      }
    };

    class AngleFloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareAngleDistance((const float16*)a, (const float16*)b, size);
      }
    };

    class NormalizedAngleFloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareNormalizedAngleDistance((const float16*)a, (const float16*)b, size);
      }
    };

    // added by Nyapicom
    class PoincareFloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::comparePoincareDistance((const float16*)a, (const float16*)b, size);
      }
    };

    // added by Nyapicom
    class LorentzFloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareLorentzDistance((const float16*)a, (const float16*)b, size);
      }
    };
#endif
#ifdef NGT_BFLOAT
    class SparseJaccardBfloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	NGTThrowException("Not supported.");
      }
    };

    class L2Bfloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareL2((const bfloat16*)a, (const bfloat16*)b, size);
      }
    };

    class NormalizedL2Bfloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareNormalizedL2((const bfloat16*)a, (const bfloat16*)b, size);
      }
    };

    class L1Bfloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	NGTThrowException("Not supported.");
      }
    };

    class CosineSimilarityBfloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	NGTThrowException("Not supported.");
      }
    };

    class NormalizedCosineSimilarityBfloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	NGTThrowException("Not supported.");
      }
    };

    class AngleBfloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	NGTThrowException("Not supported.");
      }
    };

    class NormalizedAngleBfloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	NGTThrowException("Not supported.");
      }
    };

    // added by Nyapicom
    class PoincareBfloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	NGTThrowException("Not supported.");
      }
    };

    // added by Nyapicom
    class LorentzBfloat16 {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	NGTThrowException("Not supported.");
      }
    };
#endif


};


} // namespace NGT

