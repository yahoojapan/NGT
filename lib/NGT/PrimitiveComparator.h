//
// Copyright (C) 2015-2018 Yahoo Japan Corporation
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

#if !defined(NGT_AVX_DISABLED) && defined(__AVX__)
#include	<immintrin.h>
#else
#warning "*** AVX is *NOT* available! ***"
#endif


namespace NGT {

  class MemoryCache {
  public:
    inline const static uint32_t getPrefetchPos(size_t s) { return s < 2 ? s : 2; }
    inline static void prefetch(unsigned char *ptr, const size_t byteSizeOfObject) {
      _mm_prefetch(ptr, _MM_HINT_T0);
      ptr += 64;
      _mm_prefetch(ptr, _MM_HINT_T0);
    }
  };

  class PrimitiveComparator {
  public:

    static double absolute(double v) { return fabs(v); }
    static int absolute(int v) { return abs(v); }

#if defined(NGT_AVX_DISABLED) || !defined(__AVX__)
    template <typename OBJECT_TYPE, typename COMPARE_TYPE> 
    inline static double compareL2(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      const OBJECT_TYPE *last = a + size;
      const OBJECT_TYPE *lastgroup = last - 3;
      COMPARE_TYPE diff0, diff1, diff2, diff3;
      double d = 0.0;
      while (a < lastgroup) {
	diff0 = (COMPARE_TYPE)(a[0] - b[0]);
	diff1 = (COMPARE_TYPE)(a[1] - b[1]);
	diff2 = (COMPARE_TYPE)(a[2] - b[2]);
	diff3 = (COMPARE_TYPE)(a[3] - b[3]);
	d += diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
	a += 4;
	b += 4;
      }
      while (a < last) {
	diff0 = (COMPARE_TYPE)(*a++ - *b++);
	d += diff0 * diff0;
      }
      return sqrt((double)d);
    }

    inline static double compareL2(const uint8_t *a, const uint8_t *b, size_t size) {
      return compareL2<uint8_t, int>(a, b, size);
    }

    inline static double compareL2(const float *a, const float *b, size_t size) {
      return compareL2<float, double>(a, b, size);
    }

#else
    inline static double compareL2(const float *a, const float *b, size_t size) {
      __m256 sum = _mm256_setzero_ps();
      const float *last = a + size;
      const float *lastgroup = last - 31;
      __m256 v;
      while (a < lastgroup) {
	v = _mm256_sub_ps(_mm256_loadu_ps(a), _mm256_loadu_ps(b));
	sum = _mm256_add_ps(sum, _mm256_mul_ps(v, v));
	a += 8;
	b += 8;
	v = _mm256_sub_ps(_mm256_loadu_ps(a), _mm256_loadu_ps(b));
	sum = _mm256_add_ps(sum, _mm256_mul_ps(v, v));
	a += 8;
	b += 8;
	v = _mm256_sub_ps(_mm256_loadu_ps(a), _mm256_loadu_ps(b));
	sum = _mm256_add_ps(sum, _mm256_mul_ps(v, v));
	a += 8;
	b += 8;
	v = _mm256_sub_ps(_mm256_loadu_ps(a), _mm256_loadu_ps(b));
	sum = _mm256_add_ps(sum, _mm256_mul_ps(v, v));
	a += 8;
	b += 8;
      }
      __m128 v2;
      __m128 sum2 = _mm_add_ps(_mm256_extractf128_ps(sum, 0), _mm256_extractf128_ps(sum, 1));

      while (a < last) {
	v2 = _mm_sub_ps(_mm_loadu_ps(a), _mm_loadu_ps(b));
	sum2 = _mm_add_ps(sum2, _mm_mul_ps(v2, v2));
        a += 4;
        b += 4;
      }

      float f[4];
      _mm_store_ps(f, sum2);

      double s = f[0] + f[1] + f[2] + f[3];
      return sqrt(s);
    }

    inline static double compareL2(const unsigned char *a, const unsigned char *b, size_t size) {
      __m128 sum = _mm_setzero_ps();
      const unsigned char *last = a + size;
      const unsigned char *lastgroup = last - 7;
      const __m128i zero = _mm_setzero_si128();
      while (a < lastgroup) {
	__m128i x1 = _mm_cvtepu8_epi16(*(__m128i const*)a);
	__m128i x2 = _mm_cvtepu8_epi16(*(__m128i const*)b);
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
#if defined(NGT_AVX_DISABLED) || !defined(__AVX__)
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
    inline static double compareL1(const unsigned char *a, const unsigned char *b, size_t size) {
      __m128 sum = _mm_setzero_ps();
      const unsigned char *last = a + size;
      const unsigned char *lastgroup = last - 7;
      const __m128i zero = _mm_setzero_si128();
      while (a < lastgroup) {
	__m128i x1 = _mm_cvtepu8_epi16(*(__m128i const*)a);
	__m128i x2 = _mm_cvtepu8_epi16(*(__m128i const*)b);
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
	double d = fabs((double)*a++ - (double)*b++);
	s += d;
      }
      return s;
    }
#endif
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
      size_t byteSize = sizeof(OBJECT_TYPE) * size;
      size_t n = byteSize >> 2;
      uint32_t *uinta = (uint32_t*)a;
      uint32_t *uintb = (uint32_t*)b;
      size_t count = 0;
      for (size_t i = 0; i < n; i++) {
	count += popCount(*uinta++ ^ *uintb++);
      }
      return (double)count;
    }

#if defined(NGT_AVX_DISABLED) || !defined(__AVX__)
   template <typename OBJECT_TYPE> 
    inline static double compareDotProduct(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      double sum = 0.0F;
      for (size_t loc = 0; loc < size; loc++) {
	sum += (double)a[loc] * (double)b[loc];
      }
      return sum;
    }

    template <typename OBJECT_TYPE> 
    inline static double compareCosine(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      double normA = 0.0F;
      double normB = 0.0F;
      double sum = 0.0F;
      for (size_t loc = 0; loc < size; loc++) {
	normA += (double)a[loc] * (double)a[loc];
	normB += (double)b[loc] * (double)b[loc];
	sum += (double)a[loc] * (double)b[loc];
      }

      double cosine = sum / sqrt(normA * normB);

      return cosine;
    }
#else
    inline static double compareDotProduct(const float *a, const float *b, size_t size) {
      __m256 sum = _mm256_setzero_ps();
      const float *last = a + size;
      __m256 am, bm;
      const float *lastgroup = last - 31;
      while (a < lastgroup) {
	am = _mm256_loadu_ps(a);
	bm = _mm256_loadu_ps(b);
	sum = _mm256_add_ps(sum, _mm256_mul_ps(am, bm));
	a += 8;
	b += 8;

	am = _mm256_loadu_ps(a);
	bm = _mm256_loadu_ps(b);
	sum = _mm256_add_ps(sum, _mm256_mul_ps(am, bm));
	a += 8;
	b += 8;

	am = _mm256_loadu_ps(a);
	bm = _mm256_loadu_ps(b);
	sum = _mm256_add_ps(sum, _mm256_mul_ps(am, bm));
	a += 8;
	b += 8;

	am = _mm256_loadu_ps(a);
	bm = _mm256_loadu_ps(b);
	sum = _mm256_add_ps(sum, _mm256_mul_ps(am, bm));
	a += 8;
	b += 8;
      }
      __m128 am2, bm2;
      __m128 sum2 = _mm_add_ps(_mm256_extractf128_ps(sum, 0), _mm256_extractf128_ps(sum, 1));

      while (a < last) {
	am2 = _mm_loadu_ps(a);
	bm2 = _mm_loadu_ps(b);
	sum2 = _mm_add_ps(sum2, _mm_mul_ps(am2, bm2));
	a += 4;
	b += 4;
      }
      float f[4];
      _mm_store_ps(f, sum2);

      double s = f[0] + f[1] + f[2] + f[3];
      return s;
    }

    inline static double compareDotProduct(const unsigned char *a, const unsigned char *b, size_t size) {
      double sum = 0.0F;
      for (size_t loc = 0; loc < size; loc++) {
	sum += (double)a[loc] * (double)b[loc];
      }
      return sum;
    }

    inline static double compareCosine(const float *a, const float *b, size_t size) {

      __m256 normA = _mm256_setzero_ps();
      __m256 normB = _mm256_setzero_ps();
      __m256 sum = _mm256_setzero_ps();
      const float *last = a + size;
      const float *lastgroup = last - 7;
      while (a < lastgroup) {
	__m256 am = _mm256_loadu_ps(a);
	__m256 bm = _mm256_loadu_ps(b);
	normA = _mm256_add_ps(normA, _mm256_mul_ps(am, am));
	normB = _mm256_add_ps(normB, _mm256_mul_ps(bm, bm));
	sum = _mm256_add_ps(sum, _mm256_mul_ps(am, bm));
	a += 8;
	b += 8;
      }

      __attribute__((aligned(32))) float f[8];

      _mm256_store_ps(f, normA);
      double na = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7];
      _mm256_store_ps(f, normB);
      double nb = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7];
      _mm256_store_ps(f, sum);
      double s = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7];
      while (a < last) {
	double av = *a;
	double bv = *b;
	na += av * av;
	nb += bv * bv;
	s += av * bv;
	a++;
	b++;
      }

      double cosine = s / sqrt(na * nb);
      return cosine;
    }

    inline static double compareCosine(const unsigned char *a, const unsigned char *b, size_t size) {
      double normA = 0.0F;
      double normB = 0.0F;
      double sum = 0.0F;
      for (size_t loc = 0; loc < size; loc++) {
	normA += (double)a[loc] * (double)a[loc];
	normB += (double)b[loc] * (double)b[loc];
	sum += (double)a[loc] * (double)b[loc];
      }

      double cosine = sum / sqrt(normA * normB);

      return cosine;
    }
#endif    // #if defined(NGT_AVX_DISABLED) || !defined(__AVX__)

    template <typename OBJECT_TYPE> 
    inline static double compareAngleDistance(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      double cosine = compareCosine(a, b, size);
      if (cosine >= 1.0F) {
	return 0.0F;
      } else if (cosine <= -1.0F) {
	return acos(-1.0F);
      } else {
	return acos(cosine);
      }
    }

    template <typename OBJECT_TYPE> 
    inline static double compareNormalizedAngleDistance(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      double cosine = compareDotProduct(a, b, size);
      if (cosine >= 1.0F) {
	return 0.0F;
      } else if (cosine <= -1.0F) {
	return acos(-1.0F);
      } else {
	return acos(cosine);
      }
    }

    template <typename OBJECT_TYPE> 
    inline static double compareCosineSimilarity(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      return 1.0 - compareCosine(a, b, size);
    }

    template <typename OBJECT_TYPE> 
    inline static double compareNormalizedCosineSimilarity(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
      double v = 1.0 - compareDotProduct(a, b, size);
      return v < 0.0 ? 0.0 : v;
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

    class L2Float {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
#if defined(NGT_AVX_DISABLED) || !defined(__AVX__)
	return PrimitiveComparator::compareL2<float, double>((const float*)a, (const float*)b, size);
#else
	return PrimitiveComparator::compareL2((const float*)a, (const float*)b, size);
#endif
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

    class AngleFloat {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareAngleDistance((const float*)a, (const float*)b, size);
      }
    };

    class NormalizedCosineSimilarityFloat {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareNormalizedCosineSimilarity((const float*)a, (const float*)b, size);
      }
    };

    class NormalizedAngleFloat {
    public:
      inline static double compare(const void *a, const void *b, size_t size) {
	return PrimitiveComparator::compareNormalizedAngleDistance((const float*)a, (const float*)b, size);
      }
    };

};


} // namespace NGT

