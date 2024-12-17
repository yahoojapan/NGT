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

#include "NGT/defines.h"

namespace NGT {

class MemoryCache {
 public:
  inline static void prefetch(unsigned char *ptr, const size_t byteSizeOfObject) {}
  inline static void *alignedAlloc(const size_t allocSize) { return new uint8_t[allocSize]; }
  inline static void alignedFree(void *ptr) { delete[] static_cast<uint8_t *>(ptr); }
};

class PrimitiveComparator {
 public:
  static double absolute(double v) { return fabs(v); }
  static int absolute(int v) { return abs(v); }

  template <typename OBJECT_TYPE, typename COMPARE_TYPE>
  inline static double compareL2(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
    auto *last      = a + size;
    auto *lastgroup = last - 3;
    COMPARE_TYPE diff0, diff1, diff2, diff3;
    double d = 0.0;
    while (a < lastgroup) {
      diff0 = static_cast<COMPARE_TYPE>(a[0]) - b[0];
      diff1 = static_cast<COMPARE_TYPE>(a[1]) - b[1];
      diff2 = static_cast<COMPARE_TYPE>(a[2]) - b[2];
      diff3 = static_cast<COMPARE_TYPE>(a[3]) - b[3];
      d += diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
      a += 4;
      b += 4;
    }
    while (a < last) {
      diff0 = static_cast<COMPARE_TYPE>(*a++) - static_cast<COMPARE_TYPE>(*b++);
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
  inline static double compareL2(const quint8 *a, const quint8 *b, size_t size) {
    return compareL2<quint8, float>(a, b, size);
  }

  inline static double compareL2(const qsint8 *a, const qsint8 *b, size_t size) {
    auto *i8a = reinterpret_cast<const int8_t *>(a);
    auto *i8b = reinterpret_cast<const int8_t *>(b);
    double sum = 0.0;
    for (size_t loc = 0; loc < size; loc++) {
      auto sub = static_cast<int>(*i8a) - static_cast<int>(*i8b);
      sum += sub * sub;
      i8a++;
      i8b++;
    }
    return sqrt(sum);
  }

  inline static double compareL2(const qsint8 *a, const quint8 *b, size_t size) {
    auto *i8a  = reinterpret_cast<const int8_t *>(a);
    auto *i8b  = reinterpret_cast<const uint8_t *>(b);
    double sum = 0.0;
    for (size_t loc = 0; loc < size; loc++) {
      auto sub = static_cast<int>(*i8a) - static_cast<int>(*i8b);
      sum += sub * sub;
      i8a++;
      i8b++;
    }
    return sqrt(sum);
  }

  template <typename OBJECT_TYPE>
  inline static double compareNormalizedL2(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
    double v = 2.0 - 2.0 * compareDotProduct(a, b, size);
    if (v < 0.0) {
      return 0.0;
    } else {
      return sqrt(v);
    }
  }

  template <typename OBJECT_TYPE, typename COMPARE_TYPE>
  static double compareL1(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
    auto *last      = a + size;
    auto *lastgroup = last - 3;
    COMPARE_TYPE diff0, diff1, diff2, diff3;
    double d = 0.0;
    while (a < lastgroup) {
      diff0 = (COMPARE_TYPE)(a[0]) - b[0];
      diff1 = (COMPARE_TYPE)(a[1]) - b[1];
      diff2 = (COMPARE_TYPE)(a[2]) - b[2];
      diff3 = (COMPARE_TYPE)(a[3]) - b[3];
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

  inline static double compareL1(const int8_t *a, const int8_t *b, size_t size) {
    return compareL1<int8_t, int>(a, b, size);
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
  inline static double compareL1(const quint8 *a, const quint8 *b, size_t size) {
    return compareL1<quint8, float>(a, b, size);
  }
  inline static double compareL1(const qsint8 *a, const qsint8 *b, size_t size) {
    return compareL1<qsint8, float>(a, b, size);
  }

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
    const uint32_t *last  = reinterpret_cast<const uint32_t *>(a + size);
    const uint32_t *uinta = reinterpret_cast<const uint32_t *>(a);
    const uint32_t *uintb = reinterpret_cast<const uint32_t *>(b);
    size_t count          = 0;
    while (uinta < last) {
      count += popCount(*uinta++ ^ *uintb++);
    }
    return static_cast<double>(count);
  }

  template <typename OBJECT_TYPE>
  inline static double compareJaccardDistance(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
    const uint32_t *last = reinterpret_cast<const uint32_t *>(a + size);

    const uint32_t *uinta = reinterpret_cast<const uint32_t *>(a);
    const uint32_t *uintb = reinterpret_cast<const uint32_t *>(b);
    size_t count          = 0;
    size_t countDe        = 0;
    while (uinta < last) {
      count += popCount(*uinta & *uintb);
      countDe += popCount(*uinta++ | *uintb++);
      count += popCount(*uinta & *uintb);
      countDe += popCount(*uinta++ | *uintb++);
    }

    return 1.0 - static_cast<double>(count) / static_cast<double>(countDe);
  }

  inline static double compareSparseJaccardDistance(const unsigned char *a, const unsigned char *b,
                                                    size_t size) {
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

  inline static double compareSparseJaccardDistance(const qsint8 *a, const qsint8 *b, size_t size) {
    NGTThrowException("Not supported.");
  }
  inline static double compareSparseJaccardDistance(const float *a, const float *b, size_t size) {
    size_t loca        = 0;
    size_t locb        = 0;
    const uint32_t *ai = reinterpret_cast<const uint32_t *>(a);
    const uint32_t *bi = reinterpret_cast<const uint32_t *>(b);
    size_t count       = 0;
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

  template <typename OBJECT_TYPE>
  inline static double compareDotProduct(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
    double sum = 0.0;
    for (size_t loc = 0; loc < size; loc++) {
      sum += static_cast<float>(a[loc]) * static_cast<float>(b[loc]);
    }
    return sum;
  }

  inline static double compareDotProduct(const qsint8 *a, const quint8 *b, size_t size) {
    double sum = 0.0;
    for (size_t loc = 0; loc < size; loc++) {
      sum += static_cast<int32_t>(a[loc]) * static_cast<int32_t>(b[loc]);
    }
    return sum;
  }

  template <typename OBJECT_TYPE>
  inline static double compareCosine(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
    double normA = 0.0;
    double normB = 0.0;
    double sum   = 0.0;
    for (size_t loc = 0; loc < size; loc++) {
      normA += static_cast<double>(a[loc]) * static_cast<double>(a[loc]);
      normB += static_cast<double>(b[loc]) * static_cast<double>(b[loc]);
      sum += static_cast<double>(a[loc]) * static_cast<double>(b[loc]);
    }

    double cosine = sum / sqrt(normA * normB);

    return cosine;
  }

  inline static double compareNormalizedCosineSimilarity(const float *a, const float *b, size_t size) {
    auto v = 1.0 - compareDotProduct(a, b, size);
    return v < 0.0 ? -v : v;
  }
  inline static double compareNormalizedCosineSimilarity(const float16 *a, const float16 *b, size_t size) {
    auto v = 1.0 - compareDotProduct(a, b, size);
    return v < 0.0 ? -v : v;
  }
#ifdef NGT_BFLOAT
  inline static double compareNormalizedCosineSimilarity(const bfloat16 *a, const bfloat16 *b, size_t size) {
    auto v = 1.0 - compareDotProduct(a, b, size);
    return v < 0.0 ? -v : v;
  }
#endif
  inline static double compareNormalizedCosineSimilarity(const uint8_t *a, const uint8_t *b, size_t size) {
    auto v = 1.0 - compareDotProduct(a, b, size);
    return v < 0.0 ? -v : v;
  }
  inline static double compareNormalizedCosineSimilarity(const qsint8 *a, const qsint8 *b, size_t size) {
    float max = 127.0 * 255.0 * size;
    auto v    = max - compareDotProduct(a, b, size);
    return v;
  }
  inline static double compareNormalizedCosineSimilarity(const quint8 *a, const quint8 *b, size_t size) {
    float max = 255.0 * 255.0 * size;
    auto v    = max - compareDotProduct(a, b, size);
    return v;
  }

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
  inline static double compareNormalizedAngleDistance(const OBJECT_TYPE *a, const OBJECT_TYPE *b,
                                                      size_t size) {
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
    for (size_t i = 0; i < size; i++) {
      a2 += static_cast<double>(a[i]) * static_cast<double>(a[i]);
      b2 += static_cast<double>(b[i]) * static_cast<double>(b[i]);
    }
    return std::acosh(1 + 2.0 * c2 * c2 / (1.0 - a2) / (1.0 - b2));
  }

  // added by Nyapicom
  template <typename OBJECT_TYPE>
  inline static double compareLorentzDistance(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
    // Unlike the other distance functions, this is not optimized...
    double sum = static_cast<double>(a[0]) * static_cast<double>(b[0]);
    for (size_t i = 1; i < size; i++) {
      sum -= static_cast<double>(a[i]) * static_cast<double>(b[i]);
    }
    return std::acosh(sum);
  }

  template <typename OBJECT_TYPE>
  inline static double compareCosineSimilarity(const OBJECT_TYPE *a, const OBJECT_TYPE *b, size_t size) {
    auto v = 1.0 - compareCosine(a, b, size);
    return v < 0.0 ? -v : v;
  }

  class L1Uint8;
  class L2Uint8;
  class HammingUint8;
  class JaccardUint8;
  class SparseJaccardFloat;
  class L2Float;
  class NormalizedL2Float;
  class L1Float;
  class CosineSimilarityFloat;
  class NormalizedCosineSimilarityFloat;
  class AngleFloat;
  class NormalizedAngleFloat;
  class PoincareFloat;
  class LorentzFloat;
#ifdef NGT_HALF_FLOAT
  class SparseJaccardFloat16;
  class L2Float16;
  class NormalizedL2Float16;
  class L1Float16;
  class CosineSimilarityFloat16;
  class NormalizedCosineSimilarityFloat16;
  class AngleFloat16;
  class NormalizedAngleFloat16;
  class PoincareFloat16;
  class LorentzFloat16;
#endif
#ifdef NGT_BFLOAT
  class SparseJaccardBfloat16;
  class L2Bfloat16;
  class NormalizedL2Bfloat16;
  class L1Bfloat16;
  class CosineSimilarityBfloat16;
  class NormalizedCosineSimilarityBfloat16;
  class AngleBfloat16;
  class NormalizedAngleBfloat16;
  class PoincareBfloat16;
  class LorentzBfloat16;
#endif
  class SparseJaccardQsint8;
  class L2Qsint8;
  class NormalizedL2Qsint8;
  class L1Qsint8;
  class CosineSimilarityQsint8;
  class AngleQsint8;
  class NormalizedAngleQsint8;
  class PoincareQsint8;
  class LorentzQsint8;
  class InnerProductQsint8;
  class NormalizedCosineSimilarityQuint8;
  class NormalizedCosineSimilarityQsint8;
};

} // namespace NGT
