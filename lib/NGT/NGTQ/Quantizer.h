//
// Copyright (C) 2016 Yahoo Japan Corporation
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

#include	"NGT/Index.h"
#include	"NGT/ArrayFile.h"
#include	"NGT/Clustering.h"


#define NGTQ_SIMD_BLOCK_SIZE	16	
#define NGTQ_BATCH_SIZE		2	
#define NGTQG_PREFETCH
#if defined(NGT_AVX512)
#define NGTQG_AVX512	
#warning "AVX512 is available for NGTQG"
#elif defined(NGT_AVX2)
#define NGTQG_AVX2
#warning "AVX2 is available for NGTQG"
#else
#undef NGTQG_PREFETCH
#warning "SIMD is not available for NGTQG. NGTQG might not work well."
#endif


#ifdef NGT_SHARED_MEMORY_ALLOCATOR
#define		NGTQ_SHARED_INVERTED_INDEX	
#endif

namespace NGTQ {

template <typename T>
class InvertedIndexObject {
public:
  InvertedIndexObject(): id(0) { localID[0] = 0; }
  InvertedIndexObject<T>& operator=(const InvertedIndexObject<T> &iio) {
    std::cerr << "InvertedIndexObject: operator=() should not be used, because it cannot be implemented. Instead of this, use copy()." << std::endl;
    abort();
  }
  InvertedIndexObject<T>& operator=(InvertedIndexObject<T> &iio) {
    std::cerr << "InvertedIndexObject: operator=() should not be used, because it cannot be implemented. Instead of this, use copy()." << std::endl;
    abort();
  }
  void setID(uint32_t i) { id = i; }
  static size_t headerSize() {
    InvertedIndexObject<T> dummy;
    return (size_t)&dummy.localID[0] - (size_t)&dummy;
  }
  void clear(size_t n) {
    id = 0;
    for (size_t i = 0; i < n; i++) {
      localID[i] = 0;
    }
  }
  uint32_t	id;		
  T		localID[1];	
};
 
template <typename T>
class InvertedIndexEntry : public NGT::DynamicLengthVector<InvertedIndexObject<T> > {
public:
  typedef NGT::DynamicLengthVector<InvertedIndexObject<T>> PARENT;
#ifdef NGTQ_SHARED_INVERTED_INDEX
  InvertedIndexEntry(size_t n, NGT::ObjectSpace *os = 0):numOfLocalIDs(n) {
    PARENT::elementSize = getSizeOfElement();
  }
  InvertedIndexEntry(size_t n, SharedMemoryAllocator &allocator, NGT::ObjectSpace *os = 0):numOfLocalIDs(n) {
    PARENT::elementSize = getSizeOfElement();
  }
  void pushBack(SharedMemoryAllocator &allocator) { 
    PARENT::push_back(InvertedIndexObject<T>(), allocator);
    PARENT::back(allocator).clear(numOfLocalIDs);    
  }
  void pushBack(size_t id, SharedMemoryAllocator &allocator) { 
    pushBack(allocator);
    PARENT::back(allocator).setID(id);
  }
#else // NGTQ_SHARED_INVERTED_INDEX
  InvertedIndexEntry(NGT::ObjectSpace *os = 0):numOfLocalIDs(0) {}
  InvertedIndexEntry(size_t n, NGT::ObjectSpace *os = 0):numOfLocalIDs(n) {
    PARENT::elementSize = getSizeOfElement();
  }
  void pushBack() { 
    PARENT::push_back(InvertedIndexObject<T>());
    PARENT::back().clear(numOfLocalIDs);
  }
  void pushBack(size_t id) { 
    pushBack();
    PARENT::back().setID(id);
  }

  void serialize(std::ofstream &os, NGT::ObjectSpace *objspace = 0) {
    uint32_t sz = PARENT::size();
    NGT::Serializer::write(os, sz);
    if (numOfLocalIDs > 0xFFFF) {
      std::cerr << "num of local IDs is too large. " << numOfLocalIDs << std::endl;
      abort();
    }
    uint16_t nids = numOfLocalIDs;
    NGT::Serializer::write(os, nids);
    os.write(reinterpret_cast<char*>(PARENT::vector), PARENT::size() * PARENT::elementSize);
  }

  void deserialize(std::ifstream &is, NGT::ObjectSpace *objectspace = 0) {
    uint32_t sz;
    uint16_t nids;
    try {
      NGT::Serializer::read(is, sz);
      NGT::Serializer::read(is, nids);
    } catch(NGT::Exception &err) {
      std::stringstream msg;
      msg << "InvertedIndexEntry::deserialize: It might be caused by inconsistency of the valuable type of the inverted index size. " << err.what();
      NGTThrowException(msg);
    }
    numOfLocalIDs = nids;
    PARENT::elementSize = getSizeOfElement();
    PARENT::resize(sz);
    is.read(reinterpret_cast<char*>(PARENT::vector), sz * PARENT::elementSize);
  }

#endif 

  size_t getSizeOfElement() {
    InvertedIndexObject<T> dummy;
    size_t dsize = ((numOfLocalIDs * sizeof(T) - 1) / 4 + 1) * 4;
    return InvertedIndexObject<T>::headerSize() + dsize;
  }

  size_t numOfLocalIDs;
};

class LocalDatam {
public:
  LocalDatam(){};
  LocalDatam(size_t iii, size_t iil) : iiIdx(iii), iiLocalIdx(iil) {}
  size_t iiIdx;	
  size_t iiLocalIdx; 
};

template <typename TYPE, int SIZE>
class SerializableObject : public NGT::Object {
public:
  static size_t getSerializedDataSize() { return SIZE; }
};

 enum DataType {
   DataTypeUint8 = 0,
   DataTypeFloat = 1
 };

 enum DistanceType {
   DistanceTypeNone		= 0,
   DistanceTypeL1		= 1,
   DistanceTypeL2		= 2,
   DistanceTypeHamming		= 3,
   DistanceTypeAngle		= 4,
   DistanceTypeNormalizedCosine	= 5,
   DistanceTypeNormalizedL2	= 6,
   DistanceTypeCosine		= 7	
 };

 enum CentroidCreationMode {
   CentroidCreationModeDynamic		= 0,
   CentroidCreationModeStatic		= 1,
   CentroidCreationModeDynamicKmeans	= 2,
   CentroidCreationModeNone		= 9	
 };

 enum AggregationMode {
   AggregationModeApproximateDistance				= 0,
   AggregationModeApproximateDistanceWithLookupTable		= 1,
   AggregationModeApproximateDistanceWithCache			= 2,
   AggregationModeExactDistanceThroughApproximateDistance	= 3,
   AggregationModeExactDistance					= 4
 };

 class Property {
 public:
  Property() {
    // default values
    threadSize		= 32;
    globalRange		= 200;
    localRange		= 50;
    globalCentroidLimit	= 10000000;
    localCentroidLimit	= 1000000;
    dimension		= 0;
    dataSize		= 0;
    dataType		= DataTypeFloat;
    distanceType	= DistanceTypeNone;
    singleLocalCodebook = false;
    localDivisionNo	= 8;
    batchSize		= 1000;
    centroidCreationMode = CentroidCreationModeDynamic;
    localCentroidCreationMode = CentroidCreationModeDynamic;
    localIDByteSize	= 0;		// finally decided by localCentroidLimit
    localCodebookState	= false;	// not completed
    localClusteringSampleCoefficient = 10;	
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    invertedIndexSharedMemorySize = 512; // MB
#endif
  }

  void save(const string &path) {
    NGT::PropertySet prop;
    prop.set("ThreadSize", 	(long)threadSize);
    prop.set("GlobalRange", 	globalRange);
    prop.set("LocalRange", 	localRange);
    prop.set("GlobalCentroidLimit", (long)globalCentroidLimit);
    prop.set("LocalCentroidLimit", (long)localCentroidLimit);
    prop.set("Dimension", 	(long)dimension);
    prop.set("DataSize", 	(long)dataSize);
    prop.set("DataType", 	(long)dataType);
    prop.set("DistanceType",	(long)distanceType);
    prop.set("SingleLocalCodebook", (long)singleLocalCodebook);
    prop.set("LocalDivisionNo", (long)localDivisionNo);
    prop.set("BatchSize", 	(long)batchSize);
    prop.set("CentroidCreationMode", (long)centroidCreationMode);
    prop.set("LocalCentroidCreationMode", (long)localCentroidCreationMode);
    prop.set("LocalIDByteSize",	(long)localIDByteSize);	
    prop.set("LocalCodebookState", (long)localCodebookState);
    prop.set("LocalSampleCoefficient", (long)localClusteringSampleCoefficient);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    prop.set("InvertedIndexSharedMemorySize", 	(long)invertedIndexSharedMemorySize);
#endif
    prop.save(path + "/prf");
  }

  void setupLocalIDByteSize() {
    if (localCentroidLimit > 0xffff - 1) { 
      if (localIDByteSize == 2) {
	NGTThrowException("NGTQ::Property: The localIDByteSize is illegal for the localCentroidLimit.");
      }
      localIDByteSize = 4;
    } else {
      if (localIDByteSize == INT_MAX) {
	localIDByteSize = 4;
      } else if (localIDByteSize == 0) {
	localIDByteSize = 2;
      } else {
      }
    }
    if (localIDByteSize != 2 && localIDByteSize != 4) {
      NGTThrowException("NGTQ::Property: Fatal internal error! localIDByteSize should be 2 or 4.");
    }
  }

  void load(const string &path) {
    NGT::PropertySet prop;
    prop.load(path + "/prf");
    threadSize 		= prop.getl("ThreadSize", threadSize);
    globalRange 	= prop.getf("GlobalRange", globalRange);
    localRange		= prop.getf("LocalRange", localRange);
    globalCentroidLimit	= prop.getl("GlobalCentroidLimit", globalCentroidLimit);
    localCentroidLimit	= prop.getl("LocalCentroidLimit", localCentroidLimit);
    dimension		= prop.getl("Dimension", dimension);
    dataSize		= prop.getl("DataSize", dataSize);
    dataType		= (DataType)prop.getl("DataType", dataType);
    distanceType	= (DistanceType)prop.getl("DistanceType", distanceType);
    singleLocalCodebook	= prop.getl("SingleLocalCodebook", singleLocalCodebook);
    localDivisionNo	= prop.getl("LocalDivisionNo", localDivisionNo);
    batchSize		= prop.getl("BatchSize", batchSize);
    centroidCreationMode= (CentroidCreationMode)prop.getl("CentroidCreationMode", centroidCreationMode);
    localCentroidCreationMode = (CentroidCreationMode)prop.getl("LocalCentroidCreationMode", localCentroidCreationMode);
    localIDByteSize	= prop.getl("LocalIDByteSize", INT_MAX);
    localCodebookState	= prop.getl("LocalCodebookState", localCodebookState);
    localClusteringSampleCoefficient	= prop.getl("LocalSampleCoefficient", localClusteringSampleCoefficient);
    setupLocalIDByteSize();
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    invertedIndexSharedMemorySize
      = prop.getl("InvertedIndexSharedMemorySize", invertedIndexSharedMemorySize);
#endif
  }

  void setup(const Property &p) {
    threadSize		= p.threadSize;
    globalRange		= p.globalRange;
    localRange		= p.localRange;
    globalCentroidLimit	= p.globalCentroidLimit;
    localCentroidLimit	= p.localCentroidLimit;
    distanceType	= p.distanceType;
    singleLocalCodebook = p.singleLocalCodebook;
    localDivisionNo	= p.localDivisionNo;
    batchSize		= p.batchSize;
    centroidCreationMode = p.centroidCreationMode;
    localCentroidCreationMode = p.localCentroidCreationMode;
    localIDByteSize	= p.localIDByteSize;
    localCodebookState	= p.localCodebookState;
    localClusteringSampleCoefficient = p.localClusteringSampleCoefficient;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    invertedIndexSharedMemorySize = p.invertedIndexSharedMemorySize;
#endif
  }

  inline size_t getLocalCodebookNo() { return singleLocalCodebook ? 1 : localDivisionNo; }

  size_t	threadSize;
  float		globalRange;
  float		localRange;
  size_t	globalCentroidLimit;
  size_t	localCentroidLimit;
  size_t	dimension;
  size_t	dataSize;
  DataType	dataType;
  DistanceType	distanceType;
  bool		singleLocalCodebook;
  size_t	localDivisionNo;
  size_t	batchSize;
  CentroidCreationMode centroidCreationMode;
  CentroidCreationMode localCentroidCreationMode;
  size_t	localIDByteSize;
  bool		localCodebookState;
  size_t	localClusteringSampleCoefficient;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  size_t	invertedIndexSharedMemorySize;
#endif
};

#ifdef NGTQ_DISTANCE_ANGLE
 class LocalDistanceLookup {
 public:
 LocalDistanceLookup():a(0.0), b(0.0), sum(0.0){};
   void set(double pa, double pb, double psum) {a = pa; b = pb; sum = psum;}
   double a;
   double b;
   double sum;
 };
#endif

class QuantizedObjectDistance {
public:
  class DistanceLookupTable {
  public:
    DistanceLookupTable():localDistanceLookup(0) {}
    ~DistanceLookupTable() {
      if (localDistanceLookup != 0) {
	delete[] localDistanceLookup;
	localDistanceLookup = 0;
      }
    }
    bool isValid(size_t idx) { return flag[idx]; }
#ifndef NGTQ_DISTANCE_ANGLE
    void set(size_t idx, double d) { flag[idx] = true; localDistanceLookup[idx] = d; }
    double getDistance(size_t idx) { return localDistanceLookup[idx]; }
#endif
    void initialize(size_t s) {
      size = s;
#ifdef NGTQ_DISTANCE_ANGLE
      localDistanceLookup = new LocalDistanceLookup[size];
#else
      localDistanceLookup = new float[size];
#endif
      flag.resize(size, false);
    }
#ifdef NGTQ_DISTANCE_ANGLE
    LocalDistanceLookup	*localDistanceLookup;
#else
    float		*localDistanceLookup;
#endif
    size_t		size;
    vector<bool>	flag;	
  };

  class DistanceLookupTableUint8 {
  public:
    DistanceLookupTableUint8():localDistanceLookup(0) {}
    ~DistanceLookupTableUint8() {
      if (localDistanceLookup != 0) {
	delete[] localDistanceLookup;
	localDistanceLookup = 0;
	delete[] scales;
	delete[] offsets;
      }
    }
    void initialize(size_t numOfSubspaces, size_t localCodebookCentroidNo) {
      size_t alignedNumOfSubvectors = ((numOfSubspaces - 1) / NGTQ_BATCH_SIZE + 1) * NGTQ_BATCH_SIZE;
      size = alignedNumOfSubvectors * localCodebookCentroidNo;
      localDistanceLookup = new uint8_t[size];
      scales = new float[alignedNumOfSubvectors];
      offsets = new float[alignedNumOfSubvectors];
    }

    uint8_t		*localDistanceLookup;
    size_t		size;
    size_t		aslignedNumOfSubspaces;
    size_t		localCodebookCentroidNo;
    float		*scales;
    float		*offsets;
    float		totalOffset;
  };

  QuantizedObjectDistance(){}
  virtual ~QuantizedObjectDistance() {
    delete[] localCentroids;
  }

  virtual double operator()(NGT::Object &object, size_t objectID, void *localID) = 0;

  virtual double operator()(void *localID, DistanceLookupTable &distanceLUT) = 0;

  virtual void operator()(void *inv, float *distances, size_t size, DistanceLookupTableUint8 &distanceLUT) = 0;

  virtual double operator()(NGT::Object &object, size_t objectID, void *localID, DistanceLookupTable &distanceLUT) = 0;

  template <typename T>
  inline double getAngleDistanceUint8(NGT::Object &object, size_t objectID, T localID[]) {
    assert(globalCodebook != 0);
    NGT::PersistentObject &gcentroid = *globalCodebook->getObjectSpace().getRepository().get(objectID);
    size_t sizeOfObject = globalCodebook->getObjectSpace().getByteSizeOfObject();
    size_t localDataSize = sizeOfObject / localDivisionNo / sizeof(uint8_t);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    unsigned char *gcptr = &gcentroid.at(0, globalCodebook->getObjectSpace().getRepository().allocator);
#else
    unsigned char *gcptr = &gcentroid[0];
#endif
    unsigned char *optr = &((NGT::Object&)object)[0];
    double normA = 0.0F;
    double normB = 0.0F;
    double sum = 0.0F;
    for (size_t li = 0; li < localDivisionNo; li++) {
      size_t idx = localCodebookNo == 1 ? 0 : li;
      NGT::PersistentObject &lcentroid = *localCodebook[idx].getObjectSpace().getRepository().get(localID[li]);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      float *lcptr = (float*)&lcentroid.at(0, localCodebook[idx].getObjectSpace().getRepository().allocator);
#else
      float *lcptr = (float*)&lcentroid[0];
#endif
      float *lcendptr = lcptr + localDataSize;
      while (lcptr != lcendptr) {
	double a = *optr++;
	double b = *gcptr++ + *lcptr++;
	normA += a * a;
	normB += b * b;
	sum += a * b;
      }
    }
    double cosine = sum / (sqrt(normA) * sqrt(normB));
    if (cosine >= 1.0F) {
      return 0.0F;
    } else if (cosine <= -1.0F) {
      return acos(-1.0F);
    }
    return acos(cosine);
  }

#if defined(NGT_NO_AVX)
  template <typename T>
  inline double getL2DistanceUint8(NGT::Object &object, size_t objectID, T localID[]) {
    assert(globalCodebook != 0);
    NGT::PersistentObject &gcentroid = *globalCodebook->getObjectSpace().getRepository().get(objectID);
    size_t sizeOfObject = globalCodebook->getObjectSpace().getByteSizeOfObject();
    size_t localDataSize = sizeOfObject / localDivisionNo / sizeof(uint8_t);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    unsigned char *gcptr = &gcentroid.at(0, globalCodebook->getObjectSpace().getRepository().allocator);
#else
    unsigned char *gcptr = &gcentroid[0];
#endif
    unsigned char *optr = &((NGT::Object&)object)[0];
    double distance = 0.0;
    for (size_t li = 0; li < localDivisionNo; li++) {
      size_t idx = localCodebookNo == 1 ? 0 : li;
      NGT::PersistentObject &lcentroid = *localCodebook[idx].getObjectSpace().getRepository().get(localID[li]);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      float *lcptr = (float*)&lcentroid.at(0, localCodebook[idx].getObjectSpace().getRepository().allocator);
#else
      float *lcptr = (float*)&lcentroid[0];
#endif
      double d = 0.0;
      float *lcendptr = lcptr + localDataSize;
      while (lcptr != lcendptr) {
	double sub = ((int)*optr++ - (int)*gcptr++) - *lcptr++;
	d += sub * sub;
      }
      distance += d;
    }
    return sqrt(distance);
  }
#else 
  template <typename T>
  inline double getL2DistanceUint8(NGT::Object &object, size_t objectID, T localID[]) {
    assert(globalCodebook != 0);
    NGT::PersistentObject &gcentroid = *globalCodebook->getObjectSpace().getRepository().get(objectID);
    size_t sizeOfObject = globalCodebook->getObjectSpace().getByteSizeOfObject();
    size_t localDataSize = sizeOfObject / localDivisionNo / sizeof(uint8_t);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    unsigned char *gcptr = &gcentroid.at(0, globalCodebook->getObjectSpace().getRepository().allocator);
#else
    unsigned char *gcptr = &gcentroid[0];
#endif
    unsigned char *optr = &((NGT::Object&)object)[0];
    double distance = 0.0;
    for (size_t li = 0; li < localDivisionNo; li++) {
      size_t idx = localCodebookNo == 1 ? 0 : li;
      NGT::PersistentObject &lcentroid = *localCodebook[idx].getObjectSpace().getRepository().get(localID[li]);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      float *lcptr = (float*)&lcentroid.at(0, localCodebook[idx].getObjectSpace().getRepository().allocator);
#else
      float *lcptr = (float*)&lcentroid[0];
#endif

      float *lcendptr = lcptr + localDataSize - 3;
      __m128 sum = _mm_setzero_ps();
      while (lcptr < lcendptr) {
	__m128i x1 = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i const*)optr));
	__m128i x2 = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i const*)gcptr));
	x1 = _mm_sub_epi32(x1, x2);
	__m128 sub = _mm_sub_ps(_mm_cvtepi32_ps(x1), _mm_loadu_ps(lcptr));
	sum = _mm_add_ps(sum, _mm_mul_ps(sub, sub));
	optr += 4;
	gcptr += 4;
	lcptr += 4;
      }
      __attribute__((aligned(32))) float f[4];
      _mm_store_ps(f, sum);
      double d = f[0] + f[1] + f[2] + f[3];
      while (lcptr < lcendptr) {
	double sub = ((int)*optr++ - (int)*gcptr++) - *lcptr++;
	d += sub * sub;
      }
      distance += d;
    }
    distance = sqrt(distance);
    return distance;
  }
#endif 

  template <typename T>
  inline double getAngleDistanceFloat(NGT::Object &object, size_t objectID, T localID[]) {
    assert(globalCodebook != 0);
    NGT::PersistentObject &gcentroid = *globalCodebook->getObjectSpace().getRepository().get(objectID);
    size_t sizeOfObject = globalCodebook->getObjectSpace().getByteSizeOfObject();
    size_t localDataSize = sizeOfObject / localDivisionNo  / sizeof(float);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    float *gcptr = (float*)&gcentroid.at(0, globalCodebook->getObjectSpace().getRepository().allocator);
#else
    float *gcptr = (float*)&gcentroid[0];
#endif
    float *optr = (float*)&((NGT::Object&)object)[0];
    double normA = 0.0F;
    double normB = 0.0F;
    double sum = 0.0F;
    for (size_t li = 0; li < localDivisionNo; li++) {
      size_t idx = localCodebookNo == 1 ? 0 : li;
      NGT::PersistentObject &lcentroid = *localCodebook[idx].getObjectSpace().getRepository().get(localID[li]);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      float *lcptr = (float*)&lcentroid.at(0, localCodebook[idx].getObjectSpace().getRepository().allocator);
#else
      float *lcptr = (float*)&lcentroid[0];
#endif
      float *lcendptr = lcptr + localDataSize;
      while (lcptr != lcendptr) {
	double a = *optr++;
	double b = *gcptr++ + *lcptr++;
	normA += a * a;
	normB += b * b;
	sum += a * b;
      }
    }
    double cosine = sum / (sqrt(normA) * sqrt(normB));
    if (cosine >= 1.0F) {
      return 0.0F;
    } else if (cosine <= -1.0F) {
       return acos(-1.0F);
    }
    return acos(cosine);
  }

  template <typename T>
    inline double getL2DistanceFloat(NGT::Object &object, size_t objectID, T localID[]) {
    assert(globalCodebook != 0);
    NGT::PersistentObject &gcentroid = *globalCodebook->getObjectSpace().getRepository().get(objectID);
    size_t sizeOfObject = globalCodebook->getObjectSpace().getByteSizeOfObject();
    size_t localDataSize = sizeOfObject / localDivisionNo  / sizeof(float);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    float *gcptr = (float*)&gcentroid.at(0, globalCodebook->getObjectSpace().getRepository().allocator);
#else
    float *gcptr = (float*)&gcentroid[0];
#endif
    float *optr = (float*)&((NGT::Object&)object)[0];
    double distance = 0.0;
    for (size_t li = 0; li < localDivisionNo; li++) {
      size_t idx = localCodebookNo == 1 ? 0 : li;
      NGT::PersistentObject &lcentroid = *localCodebook[idx].getObjectSpace().getRepository().get(localID[li]);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      float *lcptr = (float*)&lcentroid.at(0, localCodebook[idx].getObjectSpace().getRepository().allocator);
#else
      float *lcptr = (float*)&lcentroid[0];
#endif
      float *lcendptr = lcptr + localDataSize;
      double d = 0.0;
      while (lcptr != lcendptr) {
	double sub = (*optr++ - *gcptr++) - *lcptr++;
	d += sub * sub;
      }
      distance += d;
    }
    distance = sqrt(distance);
    return distance;
  }

#ifdef NGTQ_DISTANCE_ANGLE
  inline void createDistanceLookup(NGT::Object &object, size_t objectID, DistanceLookupTable &distanceLUT) {
    assert(globalCodebook != 0);
    NGT::Object &gcentroid = (NGT::Object &)*globalCodebook->getObjectSpace().getRepository().get(objectID);
    size_t sizeOfObject = globalCodebook->getObjectSpace().getByteSizeOfObject();
    size_t localDataSize = sizeOfObject  / localDivisionNo / sizeof(float);
    float *optr = (float*)&((NGT::Object&)object)[0];
    float *gcptr = (float*)&gcentroid[0];
    LocalDistanceLookup *dlu = distanceLUT.localDistanceLookup;
    size_t oft = 0;
    for (size_t li = 0; li < localCodebookNo; li++, oft += localDataSize) {
      dlu++;  
      for (size_t k = 1; k < localCodebookCentroidNo; k++) {
	NGT::Object &lcentroid = (NGT::Object&)*localCodebook[li].getObjectSpace().getRepository().get(k);
	float *lcptr = (float*)&lcentroid[0];		
	float *lcendptr = lcptr + localDataSize;	
	float *toptr = optr + oft;
	float *tgcptr = gcptr + oft;
	double normA = 0.0F;
	double normB = 0.0F;
	double sum = 0.0F;
	while (lcptr != lcendptr) {
	  double a = *toptr++;
	  double b = *tgcptr++ + *lcptr++;
	  normA += a * a;
	  normB += b * b;
	  sum += a * b;
	}
	dlu->set(normA, normB, sum);
	dlu++;
      }
    }
  }
#else 
  inline void createDistanceLookup(NGT::Object &object, size_t objectID, DistanceLookupTable &distanceLUT) {
    assert(globalCodebook != 0);
    NGT::Object &gcentroid = (NGT::Object &)*globalCodebook->getObjectSpace().getRepository().get(objectID);
    size_t sizeOfObject = globalCodebook->getObjectSpace().getByteSizeOfObject();

    createFloatL2DistanceLookup(&((NGT::Object&)object)[0], sizeOfObject, &gcentroid[0], distanceLUT.localDistanceLookup); 

  }
#endif 

  inline void createFloatDotProductLookup(void *object, size_t sizeOfObject, void *globalCentroid, float *lut) {
    size_t localDataSize = sizeOfObject  / localDivisionNo / sizeof(float);

    float *optr = static_cast<float*>(object);
    float *gcptr = static_cast<float*>(globalCentroid);

    size_t oft = 0;
    float *lcptr = static_cast<float*>(&localCentroids[0]);
    for (size_t li = 0; li < localCodebookNo; li++, oft += localDataSize) {
      lut++;  
      lcptr += localDataSize;
      for (size_t k = 1; k < localCodebookCentroidNo; k++) {
	float *lcendptr = lcptr + localDataSize;	
	float *toptr = optr + oft;
	float *tgcptr = gcptr + oft;
	float d = 0.0;
	while (lcptr != lcendptr) {
	  float product = *toptr++ * (*lcptr++ + *tgcptr++);
	  d += product;
	}
	*lut++ = d;
      }
    }

  }



  inline void createFloatL2DistanceLookup(void *object, size_t sizeOfObject, void *globalCentroid, float *lut) {
    size_t localDataSize = sizeOfObject  / localDivisionNo / sizeof(float);

    float *optr = static_cast<float*>(object);
    float *gcptr = static_cast<float*>(globalCentroid);

    size_t oft = 0;
    float *lcptr = static_cast<float*>(&localCentroids[0]);
    for (size_t li = 0; li < localCodebookNo; li++, oft += localDataSize) {
      lut++;  
      lcptr += localDataSize;
      for (size_t k = 1; k < localCodebookCentroidNo; k++) {
	float *lcendptr = lcptr + localDataSize;	
	float *toptr = optr + oft;
	float *tgcptr = gcptr + oft;
	float d = 0.0;
	while (lcptr != lcendptr) {
	  float sub = *toptr++ - *tgcptr++ - *lcptr++;
	  d += sub * sub;
	}
	*lut++ = d;
      }
    }
  }


  inline void createDistanceLookup(NGT::Object &object, size_t objectID, DistanceLookupTableUint8 &distanceLUT) {
    assert(globalCodebook != 0);
    size_t sizeOfObject = dimension * sizeOfType;
    float dlutmp[distanceLUT.size];
#ifdef NGTQG_DOT_PRODUCT
    createFloatDotProductLookup(&((NGT::Object&)object)[0], sizeOfObject, &gcentroid[0], dlutmp);
#else
    createFloatL2DistanceLookup(&((NGT::Object&)object)[0], sizeOfObject, globalCentroid.data(), dlutmp);
#endif
    {
      uint8_t *cdlu = distanceLUT.localDistanceLookup;
      float *dlu = dlutmp;
      distanceLUT.totalOffset = 0.0;
      float min = FLT_MAX;
      float max = -FLT_MAX;
      for (size_t li = 0; li < localCodebookNo; li++) {    
	for (size_t k = 1; k < localCodebookCentroidNo; k++) {
	  if (dlu[k] > max) {
	    max = dlu[k];
	  }
	  if (dlu[k] < min) {
	    min = dlu[k];
	  }
	}
	dlu += localCodebookCentroidNo;
      }
      float offset = min;
      float scale = (max - min) / 255.0;
      dlu = dlutmp;
      for (size_t li = 0; li < localCodebookNo; li++) {
	for (size_t k = 1; k < localCodebookCentroidNo; k++) {
	  int32_t tmp = std::round((dlu[k] - offset) / scale);
	  assert(tmp >= 0 && tmp <= 255);
	  *cdlu++ = static_cast<uint8_t>(tmp);
	}
	dlu += localCodebookCentroidNo;
	distanceLUT.offsets[li] = offset;
	distanceLUT.scales[li] = scale;
	distanceLUT.totalOffset += offset;
      }
      if ((localCodebookNo & 0x1) != 0) {
	for (size_t k = 0; k < localCodebookCentroidNo; k++) {
	  *cdlu++ = 0;
        }
	distanceLUT.offsets[localCodebookNo] = 0.0;
	distanceLUT.scales[localCodebookNo] = 0.0;
      }

    }
  }

  void set(NGT::Index *gcb, NGT::Index lcb[], size_t dn, size_t lcn, size_t sizeoftype, size_t dim) {
    globalCodebook = gcb;
    localCodebook = lcb;
    localDivisionNo = dn;
    dimension = dim;
    assert(dimension % localDivisionNo == 0);
    localDataSize = dimension / localDivisionNo;
    sizeOfType = sizeoftype;
    set(lcb, lcn);
    if (globalCodebook->getObjectSpace().getRepository().size() == 2) {
      NGT::ObjectID id = 1;
      try {
	globalCodebook->getObjectSpace().getObject(id, globalCentroid);
      } catch (NGT::Exception &err) {
	std::cerr << "Cannot load the global centroid. id=" << id << std::endl;
      }
    }

    float *lc = new float[localCodebookNo * localCodebookCentroidNo * localDataSize];
    for (size_t li = 0; li < localCodebookNo; li++) {
      for (size_t k = 1; k < localCodebookCentroidNo; k++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	NGT::PersistentObject &lcentroid = *static_cast<NGT::PersistentObject*>(localCodebook[li].getObjectSpace().getRepository().get(k));
	memcpy(&lc[li * localCodebookCentroidNo * localDataSize + k * localDataSize],
	       &lcentroid.at(0, localCodebook[li].getObjectSpace().getRepository().allocator), localDataSize * sizeof(float));
#else
	NGT::Object &lcentroid = *static_cast<NGT::Object*>(localCodebook[li].getObjectSpace().getRepository().get(k));
	memcpy(&lc[li * localCodebookCentroidNo * localDataSize + k * localDataSize], &lcentroid[0], localDataSize * sizeof(float));
#endif
      }
    }

    localCentroids = lc;



  }

  void set(NGT::Index lcb[], size_t lcn) {
    localCodebookNo = lcn;
    localCodebookCentroidNo = lcb[0].getObjectRepositorySize();
  }

  void initialize(DistanceLookupTable &c) {
    c.initialize(localCodebookNo * localCodebookCentroidNo);
  }

  void initialize(DistanceLookupTableUint8 &c) {
    c.initialize(localCodebookNo, localCodebookCentroidNo);
  }

  NGT::Index	*globalCodebook;
  NGT::Index	*localCodebook;
  size_t	localDivisionNo;
  size_t	localCodebookNo;
  size_t	localCodebookCentroidNo;

  size_t	localDataSize;
  size_t	sizeOfType;
  size_t	dimension;
  vector<float>	globalCentroid;

  float		*localCentroids;	

  size_t	localCodebookCentroidNoSIMD;

};  

template <typename T>
class QuantizedObjectDistanceUint8 : public QuantizedObjectDistance {
public:

#ifdef NGTQ_DISTANCE_ANGLE
  inline double operator()(void *l, DistanceLookupTable &distanceLUT) {
    T *localID = static_cast<T*>(l);
    double normA = 0.0F;
    double normB = 0.0F;
    double sum = 0.0F;
    for (size_t li = 0; li < localDivisionNo; li++) {
      LocalDistanceLookup &ldl = *(distanceLUT.localDistanceLookup + li * localCodebookCentroidNo + localID[li]);
      normA += ldl.a;
      normB += ldl.b;
      sum += ldl.sum;
    }
    double cosine = sum / (sqrt(normA) * sqrt(normB));
    if (cosine >= 1.0F) {
      return 0.0F;
    } else if (cosine <= -1.0F) {
      return acos(-1.0F);
    }
    return acos(cosine);
  }
  inline double operator()(NGT::Object &object, size_t objectID, void *l) {
    return getAngleDistanceUint8(object, objectID, static_cast<T*>(l));
  }
  inline double operator()(NGT::Object &object, size_t objectID, void *l, DistanceLookupTable &distanceLUT) {
    cerr << "operator() is not implemented" << endl;
    abort();
    return 0.0;
  }
#else // NGTQ_DISTANCE_ANGLE
  inline double operator()(void *l, DistanceLookupTable &distanceLUT) {
    T *localID = static_cast<T*>(l);
    double distance = 0.0;
    for (size_t li = 0; li < localDivisionNo; li++) {
      distance += distanceLUT.getDistance(li * localCodebookCentroidNo + localID[li]);
    }
    return sqrt(distance);	
  }

  inline double operator()(NGT::Object &object, size_t objectID, void *l) {
    return getL2DistanceUint8(object, objectID, static_cast<T*>(l));
  }
  inline double operator()(NGT::Object &object, size_t objectID, void *l, DistanceLookupTable &distanceLUT) {
    T *localID = static_cast<T*>(l);
    NGT::PersistentObject &gcentroid = *globalCodebook->getObjectSpace().getRepository().get(objectID);
    size_t sizeOfObject = globalCodebook->getObjectSpace().getByteSizeOfObject();
    size_t localDataSize = sizeOfObject / localDivisionNo  / sizeof(uint8_t);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    unsigned char *gcptr = &gcentroid.at(0, globalCodebook->getObjectSpace().getRepository().allocator);
#else
    unsigned char *gcptr = &gcentroid[0];
#endif
    unsigned char *optr = &((NGT::Object&)object)[0];
    double distance = 0.0;
    for (size_t li = 0; li < localDivisionNo; li++) {
      if (distanceLUT.isValid(li * localCodebookCentroidNo + localID[li])) {
	distance += distanceLUT.getDistance(li * localCodebookCentroidNo + localID[li]);
	optr += localDataSize;
	gcptr += localDataSize;
      } else {
	size_t idx = localCodebookNo == 1 ? 0 : li;
	NGT::PersistentObject &lcentroid = *localCodebook[idx].getObjectSpace().getRepository().get(localID[li]);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	float *lcptr = (float*)&lcentroid.at(0, localCodebook[idx].getObjectSpace().getRepository().allocator);
#else
	float *lcptr = (float*)&lcentroid[0];
#endif
	double d = 0.0;
	float *lcendptr = lcptr + localDataSize;
	while (lcptr != lcendptr) {
	  double sub = ((int)*optr++ - (int)*gcptr++) - *lcptr++;
	  d += sub * sub;
	}
	distance += d;
	distanceLUT.set(li * localCodebookCentroidNo + localID[li], d);
      }
    }
    return sqrt(distance);	
  }
  inline void operator()(void *inv, float *distances, size_t size, DistanceLookupTableUint8 &distanceLUT) {
    cerr << "operator is not implemented" << endl;
    abort();
  }
#endif // NGTQ_DISTANCE_ANGLE

};

template <typename T>
class QuantizedObjectDistanceFloat : public QuantizedObjectDistance {
public:

#ifdef NGTQ_DISTANCE_ANGLE
  inline double operator()(void *l, DistanceLookupTable &distanceLUT) {
    T *localID = static_cast<T*>(l);
    double normA = 0.0F;
    double normB = 0.0F;
    double sum = 0.0F;
    for (size_t li = 0; li < localDivisionNo; li++) {
      LocalDistanceLookup &ldl = *(distanceLUT.localDistanceLookup + li * localCodebookCentroidNo + localID[li]);
      normA += ldl.a;
      normB += ldl.b;
      sum += ldl.sum;
    }
    double cosine = sum / (sqrt(normA) * sqrt(normB));
    if (cosine >= 1.0F) {
      return 0.0F;
    } else if (cosine <= -1.0F) {
      return acos(-1.0F);
    }
    return acos(cosine);
  }
  inline double operator()(NGT::Object &object, size_t objectID, void *l) {
    return getAngleDistanceFloat(object, objectID, static_cast<T*>(l));
  }
  inline double operator()(NGT::Object &object, size_t objectID, void *l, DistanceLookupTable &distanceLUT) {
    cerr << "operator() is not implemented." << endl;
    abort();
    return 0.0;
  }
#else 
  inline double operator()(void *l, DistanceLookupTable &distanceLUT) {
    T *localID = static_cast<T*>(l);
    float *lut = distanceLUT.localDistanceLookup;
    float *end = lut + localCodebookCentroidNo * localDivisionNo;
    double distance = 0.0;
    while (lut < end) {
      distance += *(lut + *localID);
      localID++;
      lut += localCodebookCentroidNo;
    }
    return sqrt(distance);	
  }


#if defined(NGTQG_AVX512) || defined(NGTQG_AVX2)
  inline void operator()(void *inv, float *distances, size_t size, DistanceLookupTableUint8 &distanceLUT) {
    uint8_t *localID = static_cast<uint8_t*>(inv);
    float *d = distances;
    
#if defined(NGTQG_AVX512)
    __m512i mask512x0F = _mm512_set1_epi16(0x000f);
    __m512i mask512xF0 = _mm512_set1_epi16(0x00f0);
    constexpr size_t step512 = 32;
    const size_t range512 = (localDivisionNo >> 2) * step512;
#endif 
    const __m256i mask256x0F = _mm256_set1_epi16(0x000f);
    const __m256i mask256xF0 = _mm256_set1_epi16(0x00f0);
    constexpr size_t step256 = 16;
    const size_t range256 = (((localDivisionNo - 1) >> 1) + 1) * step256; 
    auto *last = localID + range256 / NGTQ_SIMD_BLOCK_SIZE * size;
    while (localID < last) {
      uint8_t *lut = distanceLUT.localDistanceLookup;
      auto *lastgroup256 = localID + range256;
#if defined(NGTQG_AVX512)
      __m512i depu16 = _mm512_setzero_si512();
      auto *lastgroup512 = localID + range512;
      while (localID < lastgroup512) {
	__m512i lookupTable = _mm512_loadu_si512((__m512i const*)lut);
	_mm_prefetch(&localID[0] + 64 * 8, _MM_HINT_T0); 
	__m512i packedobj = _mm512_cvtepu8_epi16(_mm256_loadu_si256((__m256i const*)&localID[0]));  
	__m512i lo = _mm512_and_si512(packedobj, mask512x0F);  
	__m512i hi = _mm512_slli_epi16(_mm512_and_si512(packedobj, mask512xF0), 4);  
	__m512i obj = _mm512_or_si512(lo, hi);  
	__m512i vtmp = _mm512_shuffle_epi8(lookupTable, obj);  
        depu16 = _mm512_adds_epu16(depu16, _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(vtmp, 0)));
	depu16 = _mm512_adds_epu16(depu16, _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(vtmp, 1)));
	lut += (localCodebookCentroidNo - 1) * 4;
	localID += step512;
      } 
#else 
      __m256i depu16l = _mm256_setzero_si256();
      __m256i depu16h = _mm256_setzero_si256();
#endif 
      while (localID < lastgroup256) {
	__m256i lookupTable = _mm256_loadu_si256((__m256i const*)lut);
	_mm_prefetch(&localID[0] + 64 * 8, _MM_HINT_T0); 
	__m256i packedobj = _mm256_cvtepu8_epi16(_mm_loadu_si128((__m128i const*)&localID[0]));
	__m256i lo = _mm256_and_si256(packedobj, mask256x0F);
	__m256i hi = _mm256_slli_epi16(_mm256_and_si256(packedobj, mask256xF0), 4);
	__m256i obj = _mm256_or_si256(lo, hi);
	__m256i vtmp = _mm256_shuffle_epi8(lookupTable, obj);
	__attribute__((aligned(32))) uint8_t v[32];
	_mm256_storeu_si256((__m256i*)&v, vtmp);

#if defined(NGTQG_AVX512)
        depu16 = _mm512_adds_epu16(depu16, _mm512_cvtepu8_epi16(vtmp));
#else
	__attribute__((aligned(32))) uint16_t v16[16];
	_mm256_storeu_si256((__m256i*)&v16, depu16l);
	_mm256_storeu_si256((__m256i*)&v16, depu16h);
	depu16l = _mm256_adds_epu16(depu16l, _mm256_cvtepu8_epi16(_mm256_extractf128_si256(vtmp, 0)));
	depu16h = _mm256_adds_epu16(depu16h, _mm256_cvtepu8_epi16(_mm256_extractf128_si256(vtmp, 1)));
#endif
	lut += (localCodebookCentroidNo - 1) * 2;
	localID += 16;
      } 

#if defined(NGTQG_AVX512)
      __m512i lo = _mm512_cvtepu16_epi32(_mm512_extracti64x4_epi64(depu16, 0));
      __m512i hi = _mm512_cvtepu16_epi32(_mm512_extracti64x4_epi64(depu16, 1));
      __m512 distance = _mm512_cvtepi32_ps(_mm512_add_epi32(lo, hi));
      __m512 scale = _mm512_broadcastss_ps(*reinterpret_cast<__m128*>(&distanceLUT.scales[0]));
      distance = _mm512_mul_ps(distance, scale);
      distance = _mm512_add_ps(distance, _mm512_set1_ps(distanceLUT.totalOffset));
#if defined(NGTQG_DOT_PRODUCT)
      float one = 1.0;
      float two = 2.0;
      distance = _mm512_mul_ps(_mm512_sub_ps(_mm512_broadcastss_ps(*reinterpret_cast<__m128*>(&one)), distance), _mm512_broadcastss_ps(*reinterpret_cast<__m128*>(&two)));
#endif
      distance = _mm512_sqrt_ps(distance);  
      _mm512_storeu_ps(d, distance);
#else 
      __m256i lol = _mm256_cvtepu16_epi32(_mm256_extractf128_si256(depu16l, 0));
      __m256i loh = _mm256_cvtepu16_epi32(_mm256_extractf128_si256(depu16l, 1));
      __m256i hil = _mm256_cvtepu16_epi32(_mm256_extractf128_si256(depu16h, 0));
      __m256i hih = _mm256_cvtepu16_epi32(_mm256_extractf128_si256(depu16h, 1));
      __m256 distancel = _mm256_cvtepi32_ps(_mm256_add_epi32(lol, hil));
      __m256 distanceh = _mm256_cvtepi32_ps(_mm256_add_epi32(loh, hih));
      __attribute__((aligned(32))) float v32[8];
      _mm256_storeu_ps((float*)&v32, distancel);
      _mm256_storeu_ps((float*)&v32, distanceh);
      __m256 scalel = _mm256_broadcastss_ps(*reinterpret_cast<__m128*>(&distanceLUT.scales[0]));
      __m256 scaleh = _mm256_broadcastss_ps(*reinterpret_cast<__m128*>(&distanceLUT.scales[0]));
      distancel = _mm256_mul_ps(distancel, scalel);
      distancel = _mm256_add_ps(distancel, _mm256_set1_ps(distanceLUT.totalOffset));
      distanceh = _mm256_mul_ps(distanceh, scaleh);
      distanceh = _mm256_add_ps(distanceh, _mm256_set1_ps(distanceLUT.totalOffset));
#if defined(NGTQG_DOT_PRODUCT)
      float one = 1.0;
      float two = 2.0;
      distancel = _mm256_mul_ps(_mm256_sub_ps(_mm256_broadcastss_ps(*reinterpret_cast<__m128*>(&one)), distancel), _mm256_broadcastss_ps(*reinterpret_cast<__m128*>(&two)));
      distanceh = _mm256_mul_ps(_mm256_sub_ps(_mm256_broadcastss_ps(*reinterpret_cast<__m128*>(&one)), distanceh), _mm256_broadcastss_ps(*reinterpret_cast<__m128*>(&two)));
#endif
      distancel = _mm256_sqrt_ps(distancel);  
      distanceh = _mm256_sqrt_ps(distanceh);  
      _mm256_storeu_ps(d, distancel);
      _mm256_storeu_ps(d + 8, distanceh);
#endif 
      d += 16;
    } 
  }


#else 
  inline void operator()(void *inv, float *distances, size_t size, DistanceLookupTableUint8 &distanceLUT) {
    uint8_t *localID = static_cast<uint8_t*>(inv);
    size_t alignedNumOfSubvectors = ((localDivisionNo - 1) / NGTQ_BATCH_SIZE + 1) * NGTQ_BATCH_SIZE;
    size_t alignedSize = ((size - 1) / 2 + 1) * 2;
    uint32_t d[NGTQ_SIMD_BLOCK_SIZE];
    size_t didx = 0;
    size_t byteSize = alignedNumOfSubvectors * alignedSize / 2;
    auto *last = localID + byteSize;
    while (localID < last) {
      uint8_t *lut = distanceLUT.localDistanceLookup;
      memset(d, 0, sizeof(uint32_t) * NGTQ_SIMD_BLOCK_SIZE);
      for (size_t li = 0; li < alignedNumOfSubvectors; li++) {
	for (size_t i = 0; i < NGTQ_SIMD_BLOCK_SIZE; i++) {
	  uint8_t obj = *localID;
	  if (i % 2 == 0) {
	    obj &= 0x0f;
	  } else {
	    obj >>= 4;
	    localID++;
	  }
	  d[i] += *(lut + obj);
	}
	lut += localCodebookCentroidNo - 1;
      }
      for (size_t i = 0; i < NGTQ_SIMD_BLOCK_SIZE; i++) {
	distances[didx + i] = sqrt(static_cast<float>(d[i]) * distanceLUT.scales[0] + distanceLUT.totalOffset);
      }
      didx += NGTQ_SIMD_BLOCK_SIZE;
    }
  }
#endif 


  inline double operator()(NGT::Object &object, size_t objectID, void *l) {
    return getL2DistanceFloat(object, objectID, static_cast<T*>(l));
  }
  inline double operator()(NGT::Object &object, size_t objectID, void *l, DistanceLookupTable &distanceLUT) {
    T *localID = static_cast<T*>(l);
    NGT::PersistentObject &gcentroid = *globalCodebook->getObjectSpace().getRepository().get(objectID);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    float *gcptr = (float*)&gcentroid.at(0, globalCodebook->getObjectSpace().getRepository().allocator);
#else
    float *gcptr = (float*)&gcentroid[0];
#endif
    float *optr = (float*)&((NGT::Object&)object)[0];
    double distance = 0.0;
    for (size_t li = 0; li < localDivisionNo; li++) {
      size_t distanceLUTidx = li * localCodebookCentroidNo + localID[li];
      if (distanceLUT.isValid(distanceLUTidx)) {
	distance += distanceLUT.getDistance(distanceLUTidx);
	optr += localDataSize;
	gcptr += localDataSize;
      } else {
        size_t idx = li;
	NGT::PersistentObject &lcentroid = *localCodebook[idx].getObjectSpace().getRepository().get(localID[li]);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	float *lcptr = (float*)&lcentroid.at(0, localCodebook[idx].getObjectSpace().getRepository().allocator);
#else
	float *lcptr = (float*)&lcentroid[0];
#endif
#if defined(NGTQG_AVX512) || defined(NGTQG_AVX2)
	float *lcendptr = lcptr + localDataSize;
	__m256 sum256 = _mm256_setzero_ps();
	__m256 v;
	while (lcptr < lcendptr) {
	  v = _mm256_sub_ps(_mm256_sub_ps(_mm256_loadu_ps(optr), _mm256_loadu_ps(gcptr)), _mm256_loadu_ps(lcptr));
	  sum256 = _mm256_add_ps(sum256, _mm256_mul_ps(v, v));
	  optr += 8;
	  gcptr += 8;
	  lcptr += 8;
	}
	__m128 sum128 = _mm_add_ps(_mm256_extractf128_ps(sum256, 0), _mm256_extractf128_ps(sum256, 1));
	__attribute__((aligned(32))) float f[4];
	_mm_store_ps(f, sum128);
	double d = f[0] + f[1] + f[2] + f[3];
#else
	float *lcendptr = lcptr + localDataSize;
	double d = 0.0;
	while (lcptr != lcendptr) {
	  double sub = (*optr++ - *gcptr++) - *lcptr++;
	  d += sub * sub;
	}
#endif
	distance += d;
      }
    }
    return sqrt(distance);	
  }
#endif 

};


class Quantizer {
public:
  typedef ArrayFile<NGT::Object>	ObjectList;	

  Quantizer(DataType dt, size_t dim) {
    property.dimension = dim;
    property.dataType = dt;
    switch (property.dataType) {
    case DataTypeUint8:
      property.dataSize = sizeof(uint8_t) * property.dimension;
      break;
    case DataTypeFloat:
      property.dataSize = sizeof(float) * property.dimension;
      break;
    default:
      cerr << "Quantizer constructor: Inner error. Invalid data type." << endl;
      break;
    } 
  }

  virtual ~Quantizer() { }

  virtual void create(const string &index,
		      NGT::Property &globalPropertySet, 
		      NGT::Property &localPropertySet) = 0;
  virtual void insert(vector<pair<NGT::Object*, size_t> > &objects) = 0;
  virtual void insert(const string &line, vector<pair<NGT::Object*, size_t> > &objects, size_t id) = 0;
  virtual void insert(vector<float> &object, vector<pair<NGT::Object*, size_t> > &objects, size_t id) = 0;
  virtual void rebuildIndex() = 0;
  virtual void save() = 0;
  virtual void open(const string &index, NGT::Property &globalProperty) = 0;
  virtual void open(const string &index) = 0;
  virtual void close() = 0;
  virtual void closeCodebooks() = 0;
#ifdef NGTQ_SHARED_INVERTED_INDEX
  virtual void reconstructInvertedIndex(const string &indexFile) = 0;
#endif

  virtual void validate() = 0;

  virtual void extractInvertedIndexObject(InvertedIndexEntry<uint16_t> &invertedIndexObjects) = 0;
  virtual void extractInvertedIndex(std::vector<std::vector<uint32_t>> &invertedIndex) = 0;

  virtual NGT::Distance getApproximateDistance(NGT::Object &query, uint32_t globalID, uint16_t *localID, QuantizedObjectDistance::DistanceLookupTable &distanceLUT) { abort(); }

  virtual void search(NGT::Object *object, NGT::ObjectDistances &objs, size_t size,
		      size_t approximateSearchSize,
		      size_t codebookSearchSize, bool resultRefinement, bool lookUpTable,
		      double epsilon) = 0;

  virtual void search(NGT::Object *object, NGT::ObjectDistances &objs, size_t size,
		      size_t approximateSearchSize,
		      size_t codebookSearchSize, AggregationMode aggregationMode,
		      double epsilon) = 0;

  virtual void search(NGT::Object *object, NGT::ObjectDistances &objs, size_t size,
		      float expansion,
		      AggregationMode aggregationMode,
		      double epsilon) = 0;

  virtual void info(ostream &os, char mode) = 0;

  virtual NGT::Index & getLocalCodebook(size_t size) = 0;

  virtual void verify() = 0;

  virtual size_t getLocalCodebookSize(size_t size) = 0;

  virtual size_t getInstanceSharedMemorySize(ostream &os, SharedMemoryAllocator::GetMemorySizeType t = SharedMemoryAllocator::GetTotalMemorySize) = 0;

  NGT::Object *allocateObject(string &line, const string &sep) {
    return globalCodebook.allocateObject(line, " \t");
  }
  NGT::Object *allocateObject(vector<double> &obj) {
    return globalCodebook.allocateObject(obj);
  }
  void deleteObject(NGT::Object *object) { globalCodebook.deleteObject(object); }
  
  void setThreadSize(size_t size) { property.threadSize = size; }
  void setGlobalRange(float r) { property.globalRange = r; }
  void setLocalRange(float r) { property.localRange = r; }
  void setGlobalCentroidLimit(size_t s) { property.globalCentroidLimit = s; }
  void setLocalCentroidLimit(size_t s) { property.localCentroidLimit = s; }
  void setDimension(size_t s) { property.dimension = s; }
  void setDistanceType(DistanceType t) { property.distanceType = t; }

  string getRootDirectory() { return rootDirectory; }

  size_t getSharedMemorySize(ostream &os, SharedMemoryAllocator::GetMemorySizeType t = SharedMemoryAllocator::GetTotalMemorySize) {
    os << "Global centroid:" << endl;
    return globalCodebook.getSharedMemorySize(os, t) + getInstanceSharedMemorySize(os, t);
  }

  virtual QuantizedObjectDistance &getQuantizedObjectDistance() = 0;

  ObjectList	objectList;
  string	rootDirectory;

  Property	property;

  NGT::Index	globalCodebook;

  size_t	distanceComputationCount;

  size_t	localIDByteSize;
  NGT::ObjectSpace::ObjectType objectType;
  size_t	divisionNo;

};

class QuantizedObjectProcessingStream {
 public:
 QuantizedObjectProcessingStream(Quantizer &quantizer, size_t numOfObjects): stream(0) {
    initialize(quantizer.divisionNo);
    setStreamSize(numOfObjects);
    stream = new uint8_t[streamSize]();
  }

  QuantizedObjectProcessingStream(size_t numOfSubspaces): stream(0) {
    initialize(numOfSubspaces);
  }

  ~QuantizedObjectProcessingStream() {
    delete[] stream;
  }

  void initialize(size_t divisionNo) {
    alignedNumOfSubvectors = ((divisionNo - 1) / NGTQ_BATCH_SIZE + 1) * NGTQ_BATCH_SIZE;
    alignedBlockSize = NGTQ_SIMD_BLOCK_SIZE * alignedNumOfSubvectors;
  }

  void setStreamSize(size_t numOfObjects) {
    alignedNumOfObjects = (((numOfObjects - 1) / NGTQ_SIMD_BLOCK_SIZE + 1) * NGTQ_SIMD_BLOCK_SIZE);
    streamSize = alignedNumOfObjects * alignedNumOfSubvectors;
    return;
  }

  void arrangeQuantizedObject(size_t dataNo, size_t subvectorNo, uint8_t quantizedObject) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    abort();
#else
    size_t blkNo = dataNo / NGTQ_SIMD_BLOCK_SIZE;	
    size_t oft = dataNo - blkNo * NGTQ_SIMD_BLOCK_SIZE;	
    stream[blkNo * alignedBlockSize + NGTQ_SIMD_BLOCK_SIZE * subvectorNo + oft] = quantizedObject;
#endif
  }

  uint8_t* compressIntoUint4() {
    size_t idx = 0;
    size_t uint4StreamSize = streamSize / 2;
    uint8_t *uint4Objects = new uint8_t[uint4StreamSize]();
    while (idx < streamSize) {
      for (size_t lidx = 0; lidx < alignedNumOfSubvectors; lidx++) {
	for (size_t bidx = 0; bidx < NGTQ_SIMD_BLOCK_SIZE; bidx++) {
	  if (idx / 2 > uint4StreamSize) {
	    std::stringstream msg;
	    msg << "Quantizer::compressIntoUint4: Fatal inner error! " << (idx / 2) << ":" << uint4StreamSize;
	    NGTThrowException(msg);
	  }
	  if (idx % 2 == 0) {
	    uint4Objects[idx / 2] = stream[idx];
	  } else {
	    uint4Objects[idx / 2] |= (stream[idx] << 4);
	  }
	  idx++;
	}
      }
    }
    return uint4Objects;
  }
  
  uint8_t*	getStream() {
    auto s = stream;
    stream = 0;
    return s;
  }

  size_t getUint4StreamSize(size_t numOfObjects) {
    setStreamSize(numOfObjects);
    return streamSize / 2;
  }

  size_t getStreamSize(size_t numOfObjects) {
    setStreamSize(numOfObjects);
    return streamSize;
  }

  uint8_t	*stream;
  size_t	alignedNumOfSubvectors;
  size_t	alignedBlockSize;
  size_t	alignedNumOfObjects;
  size_t	streamSize;
};
 
class GenerateResidualObject {
public:
  virtual ~GenerateResidualObject() {}
  virtual void operator()(size_t objectID, size_t centroidID, 
			  vector<vector<pair<NGT::Object*, size_t> > > &localObjs) = 0;

  void set(NGT::Index &gc, NGT::Index lc[], size_t dn, size_t lcn,
	   Quantizer::ObjectList *ol) {
    globalCodebook = &(NGT::GraphAndTreeIndex&)gc.getIndex();
    divisionNo = dn;
    objectList = ol;
    set(lc, lcn);
  }
  void set(NGT::Index lc[], size_t lcn) {
    localCodebook.clear();
    localCodebookNo = lcn;
    for (size_t i = 0; i < localCodebookNo; ++i) {
      localCodebook.push_back(&(NGT::GraphAndTreeIndex&)lc[i].getIndex());
    }
  }

  NGT::GraphAndTreeIndex		*globalCodebook;
  vector<NGT::GraphAndTreeIndex*>	localCodebook;
  size_t				divisionNo;
  size_t				localCodebookNo;
  Quantizer::ObjectList			*objectList;
};

class GenerateResidualObjectUint8 : public GenerateResidualObject {
public:
  void operator()(size_t objectID, size_t centroidID, 
		  vector<vector<pair<NGT::Object*, size_t> > > &localObjs) {
    NGT::PersistentObject &globalCentroid = *globalCodebook->getObjectSpace().getRepository().get(centroidID);
    NGT::Object object(&globalCodebook->getObjectSpace());
    objectList->get(objectID, object, &globalCodebook->getObjectSpace());
    size_t sizeOfObject = globalCodebook->getObjectSpace().getByteSizeOfObject();
    size_t lsize = sizeOfObject / divisionNo;
    for (size_t di = 0; di < divisionNo; di++) {
      vector<double> subObject;
      subObject.resize(lsize);
      for (size_t d = 0; d < lsize; d++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	subObject[d] = (double)object[di * lsize + d] - 
	  (double)globalCentroid.at(di * lsize + d, globalCodebook->getObjectSpace().getRepository().allocator);
#else
	subObject[d] = (double)object[di * lsize + d] - (double)globalCentroid[di * lsize + d];
#endif
      }
      size_t idx = localCodebookNo == 1 ? 0 : di;
      NGT::Object *localObj = localCodebook[idx]->allocateObject(subObject);
      localObjs[idx].push_back(pair<NGT::Object*, size_t>(localObj, 0));
    }
  }
};

class GenerateResidualObjectFloat : public GenerateResidualObject {
public:
  void operator()(size_t objectID, size_t centroidID, 
		  vector<vector<pair<NGT::Object*, size_t> > > &localObjs) {
    NGT::PersistentObject &globalCentroid = *globalCodebook->getObjectSpace().getRepository().get(centroidID);
    NGT::Object object(&globalCodebook->getObjectSpace());
    objectList->get(objectID, object, &globalCodebook->getObjectSpace());
    size_t byteSizeOfObject = globalCodebook->getObjectSpace().getByteSizeOfObject();
    size_t localByteSize = byteSizeOfObject / divisionNo;
    size_t localDimension = localByteSize / sizeof(float);
    for (size_t di = 0; di < divisionNo; di++) {
      vector<double> subObject;
      subObject.resize(localDimension);
      float *subVector = static_cast<float*>(object.getPointer(di * localByteSize));
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      float *globalCentroidSubVector = static_cast<float*>(globalCentroid.getPointer(di * localByteSize, 
										     globalCodebook->getObjectSpace().getRepository().allocator));
#else
      float *globalCentroidSubVector = static_cast<float*>(globalCentroid.getPointer(di * localByteSize));
#endif
      for (size_t d = 0; d < localDimension; d++) {
	subObject[d] = (double)subVector[d] - (double)globalCentroidSubVector[d];
      }
      size_t idx = localCodebookNo == 1 ? 0 : di;
      NGT::Object *localObj = localCodebook[idx]->allocateObject(subObject);
      localObjs[idx].push_back(pair<NGT::Object*, size_t>(localObj, 0));
    }
  }
};

template <typename LOCAL_ID_TYPE>
class QuantizerInstance : public Quantizer {
public:

  typedef void (QuantizerInstance::*AggregateObjectsFunction)(NGT::ObjectDistance &, NGT::Object *, size_t size, NGT::ObjectSpace::ResultSet &, size_t);
  typedef InvertedIndexEntry<LOCAL_ID_TYPE>	IIEntry;

  QuantizerInstance(DataType dataType, size_t dimension, size_t nCodebooks):Quantizer(dataType, dimension) {
    if (nCodebooks < 1 || nCodebooks > 100000) {
      stringstream msg;
      msg << "Quantizer::Error. Invalid divion no. " << nCodebooks;
      NGTThrowException(msg);
    }
    property.localDivisionNo = nCodebooks;
    quantizedObjectDistance = 0;
    generateResidualObject = 0;
  }

  virtual ~QuantizerInstance() { close(); }

  void createEmptyIndex(const string &index,
			NGT::Property &globalProperty, 
			NGT::Property &localProperty) 
  {
    rootDirectory = index;
    NGT::Index::mkdir(rootDirectory);
    string global = rootDirectory + "/global";
    NGT::Index::mkdir(global);

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    NGT::GraphAndTreeIndex globalCodebook(global, globalProperty);
    globalCodebook.saveIndex(global);
    globalCodebook.close();
#else
    NGT::GraphAndTreeIndex globalCodebook(globalProperty);
    globalCodebook.saveIndex(global);
    globalCodebook.close();
#endif

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    size_t localCodebookNo = property.getLocalCodebookNo();
    for (size_t i = 0; i < localCodebookNo; ++i) {
      stringstream local;
      local << rootDirectory << "/local-" << i;
      NGT::Index::mkdir(local.str());
      NGT::GraphAndTreeIndex localCodebook(local.str(), localProperty);
      localCodebook.saveIndex(local.str());
    }
#else
    NGT::GraphAndTreeIndex localCodebook(localProperty);
    size_t localCodebookNo = property.getLocalCodebookNo();
    for (size_t i = 0; i < localCodebookNo; ++i) {
      stringstream local;
      local << rootDirectory << "/local-" << i;
      NGT::Index::mkdir(local.str());
      localCodebook.saveIndex(local.str());
    }
    localCodebook.close();
#endif
#ifdef NGTQ_SHARED_INVERTED_INDEX
    invertedIndex.open(index + "/ivt", property.invertedIndexSharedMemorySize);
#else
    ofstream of(rootDirectory + "/ivt");
    invertedIndex.serialize(of);
#endif
    string fname = rootDirectory + "/obj";
    if (property.dataSize == 0) {
      NGTThrowException("Quantizer: data size of the object list is 0.");
    }
    objectList.create(fname, property.dataSize);
    objectList.open(fname);
    objectList.close();

    property.save(rootDirectory);
  }

  void open(const string &index, NGT::Property &globalProperty) {
    open(index);
    globalCodebook.setProperty(globalProperty);
  }

  void open(const string &index) {
    rootDirectory = index;
    property.load(rootDirectory);
    string globalIndex = index + "/global";
    globalCodebook.open(globalIndex);
    size_t localCodebookNo = property.getLocalCodebookNo();

    localCodebook.resize(localCodebookNo);
    for (size_t i = 0; i < localCodebookNo; ++i) {
      stringstream localIndex;
      localIndex << index << "/local-" << i;
      localCodebook[i].open(localIndex.str());
    }
#ifdef NGTQ_SHARED_INVERTED_INDEX
    invertedIndex.open(index + "/ivt", 0);
#else
    ifstream ifs(index + "/ivt");
    if (!ifs) {
      cerr << "Cannot open " << index + "/ivt" << "." << endl;
      return;
    }
    invertedIndex.deserialize(ifs);
#endif
    objectList.open(index + "/obj");
    NGT::Property globalProperty;
    globalCodebook.getProperty(globalProperty);
    size_t sizeoftype = 0;
    if (globalProperty.objectType == NGT::Property::ObjectType::Float) {
      if (property.localIDByteSize == 4) {
	quantizedObjectDistance = new QuantizedObjectDistanceFloat<uint32_t>;
      } else if (property.localIDByteSize == 2) {
	quantizedObjectDistance = new QuantizedObjectDistanceFloat<uint16_t>;
      } else {
	abort();
      }
      generateResidualObject = new GenerateResidualObjectFloat;
      sizeoftype = sizeof(float);
    } else if (globalProperty.objectType == NGT::Property::ObjectType::Uint8) {
      if (property.localIDByteSize == 4) {
	quantizedObjectDistance = new QuantizedObjectDistanceUint8<uint32_t>;
      } else if (property.localIDByteSize == 2) {
	quantizedObjectDistance = new QuantizedObjectDistanceUint8<uint16_t>;
      } else {
	abort();
      }
      generateResidualObject = new GenerateResidualObjectUint8;
      sizeoftype = sizeof(uint8_t);
    } else {
      cerr << "NGTQ::open: Fatal Inner Error: invalid object type. " << globalProperty.objectType << endl;
      cerr << "   check NGT version consistency between the caller and the library." << endl;
      assert(0);
    }
    assert(quantizedObjectDistance != 0);
    quantizedObjectDistance->set(&globalCodebook, localCodebook.data(), property.localDivisionNo, property.getLocalCodebookNo(), sizeoftype, property.dimension);
    generateResidualObject->set(globalCodebook, localCodebook.data(), property.localDivisionNo, property.getLocalCodebookNo(), &objectList);
    localIDByteSize = property.localIDByteSize;
    objectType = globalProperty.objectType;
    divisionNo = property.localDivisionNo;
  }

  void save() {
#ifndef NGT_SHARED_MEMORY_ALLOCATOR
    string global = rootDirectory + "/global";
    globalCodebook.saveIndex(global);
    size_t localCodebookNo = property.getLocalCodebookNo();
    for (size_t i = 0; i < localCodebookNo; ++i) {
      stringstream local;
      local << rootDirectory << "/local-" << i;
      try {
	NGT::Index::mkdir(local.str());
      } catch (...) {}
      localCodebook[i].saveIndex(local.str());
    }
#endif // NGT_SHARED_MEMORY_ALLOCATOR
#ifndef NGTQ_SHARED_INVERTED_INDEX
    ofstream of(rootDirectory + "/ivt");
    invertedIndex.serialize(of);
#endif
    property.save(rootDirectory);
  }

  void closeCodebooks() {
    globalCodebook.close();
    for (size_t i = 0; i < localCodebook.size(); ++i) {
      localCodebook[i].close();
    }
  }

  void close() {
    objectList.close();
    closeCodebooks();
    if (quantizedObjectDistance != 0) {
      delete quantizedObjectDistance;
      quantizedObjectDistance = 0;
    }
    if (generateResidualObject != 0) {
      delete generateResidualObject;
      generateResidualObject = 0;
    }
#ifndef NGTQ_SHARED_INVERTED_INDEX
    invertedIndex.deleteAll();
#endif
  }

#ifdef NGTQ_SHARED_INVERTED_INDEX
  void reconstructInvertedIndex(const string &invertedFile) {
    size_t size = invertedIndex.size();
#ifdef NGTQ_RECONSTRUCTION_DISABLE
    cerr << "Reconstruction is disabled!!!!!" << endl;
    return;
#endif
    cerr << "reconstructing to reduce shared memory..." << endl;
    NGT::PersistentRepository<IIEntry>	tmpInvertedIndex;
    tmpInvertedIndex.open(invertedFile, 0);
    tmpInvertedIndex.reserve(size);
    for (size_t id = 0; id < size; ++id) {
      if (invertedIndex.isEmpty(id)) {
	continue;
      }
      if (id % 100000 == 0) {
	cerr << "Processed " << id << endl;
      }
      IIEntry *entry = new(tmpInvertedIndex.getAllocator()) InvertedIndexEntry<LOCAL_ID_TYPE>(localCodebook.size(), tmpInvertedIndex.getAllocator());
      size_t esize = (*invertedIndex.at(id)).size();
      (*entry).reserve(esize, tmpInvertedIndex.getAllocator());
      for (size_t i = 0; i < esize; ++i) {
	(*entry).pushBack(tmpInvertedIndex.getAllocator());
	entry->copy((*entry).at(i, tmpInvertedIndex.getAllocator()),
				  (*invertedIndex.at(id)).at(i, invertedIndex.getAllocator()));
      }
      tmpInvertedIndex.put(id, entry);
    }
    cerr << "verifying..." << endl;
    for (size_t id = 0; id < size; ++id) {
      if (invertedIndex.isEmpty(id)) {
	continue;
      }
      if (id % 100000 == 0) {
	cerr << "Processed " << id << endl;
      }
      IIEntry &sentry = *invertedIndex.at(id);
      IIEntry &dentry = *tmpInvertedIndex.at(id);
      size_t esize = sentry.size();
      if (esize != dentry.size()) {
	cerr << id << " : size is inconsistency" << endl;
      }
      for (size_t i = 0; i < esize; ++i) {
	InvertedIndexObject<LOCAL_ID_TYPE> &sobject = (*invertedIndex.at(id)).at(i, invertedIndex.getAllocator());
	InvertedIndexObject<LOCAL_ID_TYPE> &dobject = (*tmpInvertedIndex.at(id)).at(i, tmpInvertedIndex.getAllocator());
	if (sobject.id != dobject.id) {
	  cerr << id << "," << i << " : id is inconsistency" << endl;
	}
	for (size_t d = 0; d < property.localDivisionNo; ++d) {
	  if (sobject.localID[d] != dobject.localID[d]) {
	    cerr << id << "," << i << "," << d << " : local id is inconsistency" << endl;
	  }
	}
      }
    }

    tmpInvertedIndex.close();
  }
#endif

  void createIndex(NGT::GraphAndTreeIndex &codebook,
		   size_t centroidLimit,
		   const vector<pair<NGT::Object*, size_t> > &objects, 
		   vector<NGT::Index::InsertionResult> &ids, 
		   float &range)
  {
    if (centroidLimit > 0) {
      if (getNumberOfObjects(codebook) >= centroidLimit) {
	range = FLT_MAX;
	codebook.createIndex(objects, ids, range, property.threadSize);
      } else if (getNumberOfObjects(codebook) + objects.size() > centroidLimit) {
	auto start = objects.begin();
	do {
	  size_t s = centroidLimit - getNumberOfObjects(codebook);
	  auto end = start;
	  if (std::distance(objects.begin(), start) + s >= objects.size()) {
	    end = objects.end();
	  } else {
	    end += s;
	  }
	  vector<NGT::Index::InsertionResult> idstmp;
	  vector<pair<NGT::Object*, size_t> > objtmp;
	  std::copy(start, end, std::back_inserter(objtmp));
	  codebook.createIndex(objtmp, idstmp, range, property.threadSize);
	  assert(idstmp.size() == objtmp.size());
	  std::copy(idstmp.begin(), idstmp.end(), std::back_inserter(ids));
	  start = end;
	} while (start != objects.end() && centroidLimit - getNumberOfObjects(codebook) > 0);
	range = FLT_MAX;
	vector<NGT::Index::InsertionResult> idstmp;
	vector<pair<NGT::Object*, size_t> > objtmp;
	std::copy(start, objects.end(), std::back_inserter(objtmp));
	codebook.createIndex(objtmp, idstmp, range, property.threadSize);
	std::copy(idstmp.begin(), idstmp.end(), std::back_inserter(ids));
	assert(ids.size() == objects.size());
      } else {
	codebook.createIndex(objects, ids, range, property.threadSize);
      }
    } else {
      codebook.createIndex(objects, ids, range, property.threadSize);
    }
  }

  void setGlobalCodeToInvertedEntry(NGT::Index::InsertionResult &id, pair<NGT::Object*, size_t> &object, vector<LocalDatam> &localData) {
    size_t globalCentroidID = id.id;
    if (invertedIndex.isEmpty(globalCentroidID)) {
#ifdef NGTQ_SHARED_INVERTED_INDEX
      invertedIndex.put(globalCentroidID, new(invertedIndex.allocator) InvertedIndexEntry<LOCAL_ID_TYPE>(localCodebook.size(), invertedIndex.allocator));
#else
      invertedIndex.put(globalCentroidID, new InvertedIndexEntry<LOCAL_ID_TYPE>(localCodebook.size()));
#endif
    }
    assert(!invertedIndex.isEmpty(globalCentroidID));
    IIEntry &invertedIndexEntry = *invertedIndex.at(globalCentroidID);
    if (id.identical) {
      if (property.centroidCreationMode == CentroidCreationModeDynamic) {
	assert(invertedIndexEntry.size() != 0);
      }
#ifdef NGTQ_SHARED_INVERTED_INDEX
      invertedIndexEntry.pushBack(object.second, invertedIndex.allocator);
#else
      invertedIndexEntry.pushBack(object.second);
#endif
      if (property.centroidCreationMode == CentroidCreationModeStatic) {
	localData.push_back(LocalDatam(globalCentroidID, 
				       invertedIndexEntry.size() - 1)); 
      } else {
	if (id.distance != 0.0) {
	  localData.push_back(LocalDatam(globalCentroidID, 
					 invertedIndexEntry.size() - 1)); 
	}
      }
    } else {
      if (property.centroidCreationMode != CentroidCreationModeDynamic) {
	cerr << "Quantizer: Error! Although it is an original quantizer, object has been added to the global." << endl;
	cerr << "                  Specify the size limitation of the global." << endl;
	assert(id.identical);
	abort();
      }
      if (invertedIndexEntry.size() == 0) {
#ifdef NGTQ_SHARED_INVERTED_INDEX
	invertedIndexEntry.pushBack(object.second, invertedIndex.allocator);
#else
	invertedIndexEntry.pushBack(object.second);
#endif
      } else {
#ifdef NGTQ_SHARED_INVERTED_INDEX
	invertedIndexEntry.at(0, invertedIndex.allocator).setID(object.second);
#else
	invertedIndexEntry[0].setID(object.second);
#endif
      }
    }
  }

  void setSingleLocalCodeToInvertedIndexEntry(vector<NGT::GraphAndTreeIndex*> &lcodebook, vector<LocalDatam> &localData, vector<vector<pair<NGT::Object*, size_t> > > &localObjs) {
    float lr = property.localRange;
    size_t localCentroidLimit = property.localCentroidLimit;
    if (property.localCodebookState) {
      lr = FLT_MAX;	
      localCentroidLimit = 0;
    }
    vector<NGT::Index::InsertionResult> lids;
    createIndex(*lcodebook[0], localCentroidLimit, localObjs[0], lids, lr);
    size_t divisionNo = property.localDivisionNo;
    for (size_t i = 0; i < localData.size(); i++) {
      for (size_t di = 0; di < divisionNo; di++) {
	size_t id = lids[i * divisionNo + di].id;
	assert(!property.localCodebookState || id <= ((1UL << (sizeof(LOCAL_ID_TYPE) * 8)) - 1)); 
#ifdef NGTQ_SHARED_INVERTED_INDEX
	(*invertedIndex.at(localData[i].iiIdx)).at(localData[i].iiLocalIdx, invertedIndex.allocator).localID[di] = id;
#else
	(*invertedIndex.at(localData[i].iiIdx))[localData[i].iiLocalIdx].localID[di] = id;
#endif
      }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      localCodebook[0].deleteObject(localObjs[0][i].first);
#else
      if (lids[i].identical) {
	localCodebook[0].deleteObject(localObjs[0][i].first);
      }
#endif
    }
  }

  bool setMultipleLocalCodeToInvertedIndexEntry(vector<NGT::GraphAndTreeIndex*> &lcodebook, vector<LocalDatam> &localData, vector<vector<pair<NGT::Object*, size_t> > > &localObjs) {
    size_t localCodebookNo = property.getLocalCodebookNo();
    bool localCodebookFull = true;  
    for (size_t li = 0; li < localCodebookNo; ++li) {
      float lr = property.localRange;
      size_t localCentroidLimit = property.localCentroidLimit;
      if (property.localCentroidCreationMode == CentroidCreationModeDynamicKmeans) {
	localCentroidLimit *= property.localClusteringSampleCoefficient;
      }
      if (property.localCodebookState) {
	lr = FLT_MAX;	
	localCentroidLimit = 0;
      } else {
	if (property.localCentroidCreationMode == CentroidCreationModeDynamicKmeans) {
	  lr = -1.0;
	}
      }
      vector<NGT::Index::InsertionResult> lids;
      createIndex(*lcodebook[li], localCentroidLimit, localObjs[li], lids, lr);
      if (lr != FLT_MAX) { 
	localCodebookFull = false;
      }
      assert(localData.size() == lids.size());
      for (size_t i = 0; i < localData.size(); i++) {
	size_t id = lids[i].id;
	assert(!property.localCodebookState || id <= ((1UL << (sizeof(LOCAL_ID_TYPE) * 8)) - 1)); 
#ifdef NGTQ_SHARED_INVERTED_INDEX
	(*invertedIndex.at(localData[i].iiIdx)).at(localData[i].iiLocalIdx, invertedIndex.allocator).localID[li] = id;
#else
	(*invertedIndex.at(localData[i].iiIdx))[localData[i].iiLocalIdx].localID[li] = id;
#endif
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	localCodebook[li].deleteObject(localObjs[li][i].first);
#else
	if (lids[i].identical) {
	  localCodebook[li].deleteObject(localObjs[li][i].first);
	}
#endif
      } 
    } 
    return localCodebookFull;
  }
  
  void buildMultipleLocalCodebooks(NGT::Index *localCodebook, size_t localCodebookNo, size_t numberOfCentroids) {
    NGT::Clustering clustering;
    clustering.epsilonFrom = 0.10;
    clustering.epsilonTo = 0.50;
    clustering.epsilonStep = 0.05;
    clustering.maximumIteration = 20;
    for (size_t li = 0; li < localCodebookNo; ++li) {
      double diff = clustering.kmeansWithNGT(localCodebook[li], numberOfCentroids);
      if (diff > 0.0) {
	cerr << "Not converge. " << diff << endl;
      }
      cerr << "Clustering of the subvector is complete. " << localCodebook[li].getPath() << ":" << localCodebook[li].getObjectRepositorySize() << endl;
    }
  }

  void replaceInvertedIndexEntry(size_t localCodebookNo) {
    vector<LocalDatam> localData;
    for (size_t gidx = 1; gidx < invertedIndex.size(); gidx++) {
      if (invertedIndex.at(gidx) == 0) {
	std::cerr << "replaceInvertedIndexEntry: Warning. empty inverted index entry. " << gidx << ":" << invertedIndex.size() << std::endl;
	continue;
      }
      IIEntry &invertedIndexEntry = *invertedIndex.at(gidx);
      for (size_t oi = property.centroidCreationMode == CentroidCreationModeStatic ? 0 : 1;
	   oi < invertedIndexEntry.size(); oi++) {
	localData.push_back(LocalDatam(gidx, oi)); 
      }
    }
    vector<vector<pair<NGT::Object*, size_t> > > localObjs;
    localObjs.resize(localCodebookNo);	
    for (size_t i = 0; i < localData.size(); i++) {
      IIEntry &invertedIndexEntry = *invertedIndex.at(localData[i].iiIdx);
#ifdef NGTQ_SHARED_INVERTED_INDEX
      (*generateResidualObject)(invertedIndexEntry.at(localData[i].iiLocalIdx, invertedIndex.allocator).id,
				localData[i].iiIdx, // centroid:ID of global codebook
				localObjs);
#else
      (*generateResidualObject)(invertedIndexEntry[localData[i].iiLocalIdx].id,
				localData[i].iiIdx, // centroid:ID of global codebook
				localObjs);
#endif
    }
    vector<NGT::GraphAndTreeIndex*> lcodebook;
    for (size_t i = 0; i < localCodebookNo; i++) {
      lcodebook.push_back(&(NGT::GraphAndTreeIndex &)localCodebook[i].getIndex());
    }
    setMultipleLocalCodeToInvertedIndexEntry(lcodebook, localData, localObjs);
  }

  void insert(vector<pair<NGT::Object*, size_t> > &objects) {
    NGT::GraphAndTreeIndex &gcodebook = (NGT::GraphAndTreeIndex &)globalCodebook.getIndex();
    vector<NGT::GraphAndTreeIndex*> lcodebook;
    size_t localCodebookNo = property.getLocalCodebookNo();
    for (size_t i = 0; i < localCodebookNo; i++) {
      lcodebook.push_back(&(NGT::GraphAndTreeIndex &)localCodebook[i].getIndex());
    }
    float gr = property.globalRange;
    vector<NGT::Index::InsertionResult> ids;	
    createIndex(gcodebook, property.globalCentroidLimit, objects, ids, gr);
#ifdef NGTQ_SHARED_INVERTED_INDEX
    if (invertedIndex.getAllocatedSize() <= invertedIndex.size() + objects.size()) {
      invertedIndex.reserve(invertedIndex.getAllocatedSize() * 2);
    }
#else
    invertedIndex.reserve(invertedIndex.size() + objects.size());
#endif
    vector<LocalDatam> localData;
    for (size_t i = 0; i < ids.size(); i++) {
      setGlobalCodeToInvertedEntry(ids[i], objects[i], localData);
    } 
    vector<vector<pair<NGT::Object*, size_t> > > localObjs;
    localObjs.resize(property.getLocalCodebookNo());	
    for (size_t i = 0; i < localData.size(); i++) {
      IIEntry &invertedIndexEntry = *invertedIndex.at(localData[i].iiIdx);
#ifdef NGTQ_SHARED_INVERTED_INDEX
      (*generateResidualObject)(invertedIndexEntry.at(localData[i].iiLocalIdx, invertedIndex.allocator).id,
				localData[i].iiIdx, // centroid:ID of global codebook
				localObjs);
#else
      (*generateResidualObject)(invertedIndexEntry[localData[i].iiLocalIdx].id,
				localData[i].iiIdx, // centroid:ID of global codebook
				localObjs);
#endif
    }
    if (property.singleLocalCodebook) {
      // single local codebook
      setSingleLocalCodeToInvertedIndexEntry(lcodebook, localData, localObjs);
    } else { 
      // multiple local codebooks
      bool localCodebookFull = setMultipleLocalCodeToInvertedIndexEntry(lcodebook, localData, localObjs);
      if ((!property.localCodebookState) && localCodebookFull) { 
	if (property.localCentroidCreationMode == CentroidCreationModeDynamicKmeans) {
	  buildMultipleLocalCodebooks(localCodebook.data(), localCodebookNo, property.localCentroidLimit);
	  (*generateResidualObject).set(localCodebook.data(), localCodebookNo);
	  property.localCodebookState = true;	
	  localCodebookFull = false;		
	  replaceInvertedIndexEntry(localCodebookNo);
	} else {
	  property.localCodebookState = true;	
	  localCodebookFull = false;
	}
      }
    } 
    for (size_t i = 0; i < objects.size(); i++) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      globalCodebook.deleteObject(objects[i].first);
#else
      if (ids[i].identical == true) {
	globalCodebook.deleteObject(objects[i].first);
      }
#endif
    }
    objects.clear();
  }

  void insert(const string &line, vector<pair<NGT::Object*, size_t> > &objects, size_t count) {
    size_t id = count;
    if (count == 0) {
      id = objectList.size();
      id = id == 0 ? 1 : id;
    }
    NGT::Object *object = globalCodebook.allocateObject(line, " \t");
    objectList.put(id, *object, &globalCodebook.getObjectSpace());
    objects.push_back(pair<NGT::Object*, size_t>(object, id));
    if (objects.size() >= property.batchSize) {
      insert(objects);  // batch insert
    }
  }

  void insert(vector<float> &objvector, vector<pair<NGT::Object*, size_t> > &objects, size_t count) {
    size_t id = count;
    if (count == 0) {
      id = objectList.size();
      id = id == 0 ? 1 : id;
    }

    NGT::Object *object = globalCodebook.allocateObject(objvector);
    objectList.put(id, *object, &globalCodebook.getObjectSpace());

    objects.push_back(pair<NGT::Object*, size_t>(object, id));

    if (objects.size() >= property.batchSize) {
      insert(objects);   // batch insert
    }
  }

  void rebuildIndex() {
    vector<pair<NGT::Object*, size_t> > objects;
    size_t objectCount = objectList.size();
    size_t count = 0;
    for (size_t idx = 1; idx < objectCount; idx++) {
      count++;
      if (count % 100000 == 0) {
	  cerr << "Processed " << count;
	  cerr << endl;
      }
      NGT::Object *object = globalCodebook.getObjectSpace().allocateObject();
      objectList.get(idx, *object, &globalCodebook.getObjectSpace());
      objects.push_back(pair<NGT::Object*, size_t>(object, idx));
      if (objects.size() >= property.batchSize) {
	insert(objects);
      }
    }
    if (objects.size() >= 0) {
      insert(objects);
    }
  }

  void create(const string &index,
	      NGT::Property &globalProperty, 
	      NGT::Property &localProperty
	      ) {

    if (property.localCentroidLimit > ((1UL << (sizeof(LOCAL_ID_TYPE) * 8)) - 1)) {
      stringstream msg;
      msg << "Quantizer::Error. Local centroid limit is too large. " << property.localCentroidLimit << " It must be less than " << (1UL << (sizeof(LOCAL_ID_TYPE) * 8));
      NGTThrowException(msg);
    }

    NGT::Property gp;
    NGT::Property lp;

    gp.setDefault();
    lp.setDefault();

    gp.batchSizeForCreation = 500;
    gp.edgeSizeLimitForCreation = 0;
    gp.edgeSizeForCreation = 10;
    gp.graphType = NGT::Index::Property::GraphType::GraphTypeANNG;
    gp.insertionRadiusCoefficient = 1.1;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    gp.graphSharedMemorySize	= 512; // MB
    gp.treeSharedMemorySize	= 512; // MB
    gp.objectSharedMemorySize	= 512; // MB  512 is for up to 20M objects.
#endif

    lp.batchSizeForCreation = 500;
    lp.edgeSizeLimitForCreation = 0;
    lp.edgeSizeForCreation = 10;
    lp.graphType = NGT::Index::Property::GraphType::GraphTypeANNG;
    lp.insertionRadiusCoefficient = 1.1;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    lp.graphSharedMemorySize	= 128; // MB
    lp.treeSharedMemorySize	= 128; // MB
    lp.objectSharedMemorySize	= 128; // MB  128 is for up to 5M objects?
#endif

    gp.set(globalProperty);
    lp.set(localProperty);

    gp.edgeSizeForSearch = 40;	
    lp.edgeSizeForSearch = 40;	

    lp.objectType = NGT::Index::Property::ObjectType::Float;

    gp.dimension = property.dimension;
    if (gp.dimension == 0) {
      stringstream msg;
      msg << "NGTQ::Quantizer::create: specified dimension is zero!";
      NGTThrowException(msg);
    }
    if (property.localDivisionNo != 1 && property.dimension % property.localDivisionNo != 0) {
      stringstream msg;
      msg << "NGTQ::Quantizer::create: dimension and localDivisionNo are not proper. "
	  << property.dimension << ":" << property.localDivisionNo;
      NGTThrowException(msg);
    }
    lp.dimension = property.dimension / property.localDivisionNo;

    switch (property.dataType) {
    case DataTypeFloat:
      gp.objectType = NGT::Index::Property::ObjectType::Float;
      break;
    case DataTypeUint8:
      gp.objectType = NGT::Index::Property::ObjectType::Uint8;
      break;
    default:
      {
	stringstream msg;
	msg << "NGTQ::Quantizer::create: Inner error! Invalid data type.";
	NGTThrowException(msg);
      }
    }

    switch (property.distanceType) {
    case DistanceTypeL1:
      gp.distanceType = NGT::Index::Property::DistanceType::DistanceTypeL1;
      lp.distanceType = NGT::Index::Property::DistanceType::DistanceTypeL1;
      break;
    case DistanceTypeL2:
#ifdef NGTQ_DISTANCE_ANGLE
      {
	stringstream msg;
	msg << "NGTQ::Quantizer::create: L2 is unavailable!!! you have to rebuild.";
	NGTThrowException(msg);
      }
#endif
      gp.distanceType = NGT::Index::Property::DistanceType::DistanceTypeL2;
      lp.distanceType = NGT::Index::Property::DistanceType::DistanceTypeL2;
      break;
    case DistanceTypeHamming:
      gp.distanceType = NGT::Index::Property::DistanceType::DistanceTypeHamming;
      lp.distanceType = NGT::Index::Property::DistanceType::DistanceTypeHamming;
      break;
    case DistanceTypeAngle:
#ifndef NGTQ_DISTANCE_ANGLE
      {
	stringstream msg;
	msg << "NGTQ::Quantizer::create: Angle is unavailable!!! you have to rebuild.";
	NGTThrowException(msg);
      }
#endif
      gp.distanceType = NGT::Index::Property::DistanceType::DistanceTypeAngle;
      lp.distanceType = NGT::Index::Property::DistanceType::DistanceTypeAngle;
      break;
    case DistanceTypeNormalizedCosine:
      gp.distanceType = NGT::Index::Property::DistanceType::DistanceTypeNormalizedCosine;
      lp.distanceType = NGT::Index::Property::DistanceType::DistanceTypeL2;
      break;
    case DistanceTypeCosine:
      gp.distanceType = NGT::Index::Property::DistanceType::DistanceTypeCosine;
      lp.distanceType = NGT::Index::Property::DistanceType::DistanceTypeL2;
      break;
    case DistanceTypeNormalizedL2:
      gp.distanceType = NGT::Index::Property::DistanceType::DistanceTypeNormalizedL2;
      lp.distanceType = NGT::Index::Property::DistanceType::DistanceTypeL2;
      break;
    default:
      {
	stringstream msg;
	msg << "NGTQ::Quantizer::create Inner error! Invalid distance type.";
	NGTThrowException(msg);
      }
    }

    createEmptyIndex(index, gp, lp);
  }

  void validate() {
    size_t gcbSize = globalCodebook.getObjectRepositorySize();
    cerr << "global codebook size=" << gcbSize << endl;
    for (size_t gidx = 1; gidx < 4 && gidx < gcbSize; gidx++) {
      if (invertedIndex[gidx] == 0) {
	cerr << "something wrong" << endl;
	exit(1);
      }
      cerr << gidx << " inverted index size=" << (*invertedIndex[gidx]).size() << endl;
      if ((*invertedIndex[gidx]).size() == 0) {
	cerr << "something wrong" << endl;
	continue;
      }

      NGT::PersistentObject &gcentroid = *globalCodebook.getObjectSpace().getRepository().get(gidx);
      vector<double> gco;
      globalCodebook.getObjectSpace().getRepository().extractObject(&gcentroid, gco);
      cerr << "global centroid object(" << gco.size() << ")=";
      for (size_t i = 0; i < gco.size(); i++) {
	cerr << gco[i] << " ";
      }
      cerr << endl;

      {
#ifdef NGTQ_SHARED_INVERTED_INDEX
	InvertedIndexObject<LOCAL_ID_TYPE> &invertedIndexEntry = (*invertedIndex[gidx]).at(0, invertedIndex.allocator);
#else
	InvertedIndexObject<LOCAL_ID_TYPE> &invertedIndexEntry = (*invertedIndex[gidx])[0];
#endif
        if (invertedIndexEntry.id != gidx) {
	  cerr << "a global centroid id is wrong in the inverted index." << gidx << ":" << invertedIndexEntry.id << endl;
	  exit(1);
        }
      }
      NGT::Object *gcentroidFromList = globalCodebook.getObjectSpace().getRepository().allocateObject();
      objectList.get(gidx, *gcentroidFromList, &globalCodebook.getObjectSpace());
      vector<double> gcolist;
      globalCodebook.getObjectSpace().getRepository().extractObject(gcentroidFromList, gcolist);
      if (gco != gcolist) {
	cerr << "Fatal error! centroid in NGT is different from object list in NGTQ" << endl;
	exit(1);
      }
      vector<size_t> elements;
      for (size_t iidx = 0; iidx < (*invertedIndex[gidx]).size(); iidx++) {
#ifdef NGTQ_SHARED_INVERTED_INDEX
	InvertedIndexObject<LOCAL_ID_TYPE> &invertedIndexEntry = (*invertedIndex[gidx]).at(iidx, invertedIndex.allocator);
#else
	InvertedIndexObject<LOCAL_ID_TYPE> &invertedIndexEntry = (*invertedIndex[gidx])[iidx];
#endif
        elements.push_back(invertedIndexEntry.id);
	cerr << "  object ID=" << invertedIndexEntry.id;
	{
	  NGT::Object *o = globalCodebook.getObjectSpace().getRepository().allocateObject();
	  objectList.get(invertedIndexEntry.id, *o, &globalCodebook.getObjectSpace());
	  NGT::Distance distance = globalCodebook.getObjectSpace().getComparator()(*gcentroidFromList, *o);
	  cerr << ":distance=" << distance;
	}
	cerr << ":local codebook IDs=";
	for (size_t li = 0; li < property.localDivisionNo; li++) {
	  cerr << invertedIndexEntry.localID[li] << " ";
	}
	cerr << endl;
	for (size_t li = 0; li < property.localDivisionNo; li++) {
	  if (invertedIndexEntry.localID[li] == 0) {
	    if (property.centroidCreationMode != CentroidCreationModeStatic) {
	      if (iidx == 0) {
		break;
	      }
	    }
	    cerr << "local ID is unexpected zero." << endl;
	  }
	}
      }
      vector<size_t> ngid;
      {
	size_t resultSize = 30;
	size_t approximateSearchSize = 1000;
	size_t codebookSearchSize = 50;
	bool refine = true;
	bool lookuptable = false;
	double epsilon = 0.1;
	NGT::ObjectDistances objects;
	search(gcentroidFromList, objects, resultSize, approximateSearchSize, codebookSearchSize, 
	       refine, lookuptable, epsilon);
	for (size_t resulti = 0; resulti < objects.size(); resulti++) {
	  if (std::find(elements.begin(), elements.end(), objects[resulti].id) != elements.end()) {
	    cerr << "  ";
	  } else {
	    cerr << "x ";
	    ngid.push_back(objects[resulti].id);
	    NGT::ObjectDistances result;
	    NGT::Object *o = globalCodebook.getObjectSpace().getRepository().allocateObject();
	    objectList.get(ngid.back(), *o, &globalCodebook.getObjectSpace());
	    NGT::GraphAndTreeIndex &graphIndex = (NGT::GraphAndTreeIndex &)globalCodebook.getIndex();
	    graphIndex.searchForNNGInsertion(*o, result);
	    if (result[0].distance > objects[resulti].distance) {
	      cerr << " Strange! ";
	      cerr << result[0].distance << ":" << objects[resulti].distance << " ";
	    }
	    globalCodebook.getObjectSpace().getRepository().deleteObject(o);
	  }
	  cerr << "  search object " << resulti << " ID=" << objects[resulti].id << " distance=" << objects[resulti].distance << endl;
	}
      }
      globalCodebook.getObjectSpace().getRepository().deleteObject(gcentroidFromList);
    }
  }

  void searchGlobalCodebook(NGT::Object *query, size_t size, NGT::ObjectDistances &objects,
			    size_t &approximateSearchSize,
			    size_t codebookSearchSize, 
			    double epsilon) {

    NGT::SearchContainer sc(*query);
    sc.setResults(&objects);
    sc.size = codebookSearchSize;
    sc.radius = FLT_MAX;
    sc.explorationCoefficient = epsilon + 1.0;
    if (epsilon >= FLT_MAX) {
      globalCodebook.linearSearch(sc);
    } else {
      globalCodebook.search(sc);
    }

  }

  inline void aggregateObjectsWithExactDistance(NGT::ObjectDistance &globalCentroid, NGT::Object *query, size_t size, NGT::ObjectSpace::ResultSet &results, size_t approximateSearchSize) {
    NGT::ObjectSpace &objectSpace = globalCodebook.getObjectSpace();
    for (size_t j = 0; j < invertedIndex[globalCentroid.id]->size() && results.size() < approximateSearchSize; j++) {
#ifdef NGTQ_SHARED_INVERTED_INDEX
      InvertedIndexObject<LOCAL_ID_TYPE> &invertedIndexEntry = (*invertedIndex[globalCentroid.id]).at(j, invertedIndex.allocator);
#else
      InvertedIndexObject<LOCAL_ID_TYPE> &invertedIndexEntry = (*invertedIndex[globalCentroid.id])[j];
#endif
      double distance;
      if (invertedIndexEntry.localID[0] == 0) {
	distance = globalCentroid.distance;
      } else { 
	NGT::Object o(&objectSpace);
	objectList.get(invertedIndexEntry.id, (NGT::Object&)o, &objectSpace);
	distance = objectSpace.getComparator()(*query, (NGT::Object&)o);
      }  

      NGT::ObjectDistance obj;
      obj.id = invertedIndexEntry.id;
      obj.distance = distance;
      assert(obj.id > 0);
      results.push(obj);

    } 
  }

   inline void aggregateObjectsWithLookupTable(NGT::ObjectDistance &globalCentroid, NGT::Object *query, size_t size, NGT::ObjectSpace::ResultSet &results, size_t approximateSearchSize) {
     QuantizedObjectDistance::DistanceLookupTable distanceLUT;
     (*quantizedObjectDistance).initialize(distanceLUT);
     (*quantizedObjectDistance).createDistanceLookup(*query, globalCentroid.id, distanceLUT);
     for (size_t j = 0; j < invertedIndex[globalCentroid.id]->size() && results.size() < approximateSearchSize; j++) {
#ifdef NGTQ_SHARED_INVERTED_INDEX
       InvertedIndexObject<LOCAL_ID_TYPE> &invertedIndexEntry = (*invertedIndex[globalCentroid.id]).at(j, invertedIndex.allocator);
#else
       InvertedIndexObject<LOCAL_ID_TYPE> &invertedIndexEntry = (*invertedIndex[globalCentroid.id])[j];
#endif
       double distance;
       if (invertedIndexEntry.localID[0] == 0) {
	 distance = globalCentroid.distance;
       } else { 
	 distance = (*quantizedObjectDistance)(invertedIndexEntry.localID, distanceLUT);
       }  


       NGT::ObjectDistance obj;
       obj.id = invertedIndexEntry.id;
       obj.distance = distance;
       assert(obj.id > 0);
       results.push(obj);

     } 
  }

   void extractInvertedIndexObject(InvertedIndexEntry<uint16_t> &invertedIndexObjects) {
#ifdef NGTQ_SHARED_INVERTED_INDEX
     std::cerr << "not implemented" << std::endl;
     abort();
#else
    size_t lastID = 0;
    for (size_t gid = 0; gid < invertedIndex.size(); gid++) {
      if (invertedIndex[gid] == 0) {
	continue;
      }
      for (size_t idx = 0; idx < invertedIndex[gid]->size(); idx++) {
#ifdef NGTQ_SHARED_INVERTED_INDEX
        NGTQ::InvertedIndexObject<LOCAL_ID_TYPE> &entry = (*invertedIndex[gid]).at(idx, invertedIndex.allocator);
#else
        NGTQ::InvertedIndexObject<LOCAL_ID_TYPE> &entry = (*invertedIndex[gid])[idx];
#endif
	if (entry.id > lastID) {
	  lastID = entry.id;
	}
      }
    }
    invertedIndexObjects.resize(lastID + 1);
    for (size_t gid = 1; gid < invertedIndex.size(); gid++) {
      for (size_t idx = 0; idx < invertedIndex[gid]->size(); idx++) {
#ifdef NGTQ_SHARED_INVERTED_INDEX
        NGTQ::InvertedIndexObject<LOCAL_ID_TYPE> &entry = (*invertedIndex[gid]).at(idx, invertedIndex.allocator);
#else
        NGTQ::InvertedIndexObject<LOCAL_ID_TYPE> &entry = (*invertedIndex[gid])[idx];
#endif
	invertedIndexObjects[entry.id].id = entry.id;
	if (sizeof(entry.localID[0]) > sizeof(invertedIndexObjects[entry.id].localID[0])) {
	  std::cerr << "you should change the object ID type." << std::endl;
	  abort();
	}
	for (size_t i = 0; i < invertedIndexObjects.numOfLocalIDs; i++) {
	  invertedIndexObjects[entry.id].localID[i] = entry.localID[i];
	}
	assert(invertedIndexObjects[entry.id].localID[0] == entry.localID[0]);
      }
    } 
#endif 
  }

  void extractInvertedIndex(std::vector<std::vector<uint32_t>> &ii) { 
    ii.resize(invertedIndex.size());
    for (size_t gid = 1; gid < invertedIndex.size(); gid++) {
      if (invertedIndex[gid] == 0 || invertedIndex[gid]->size() == 0) {
	continue;
      }
      ii[gid].reserve(invertedIndex[gid]->size());
      for (size_t idx = 0; idx < invertedIndex[gid]->size(); idx++) {
#ifdef NGTQ_SHARED_INVERTED_INDEX
	NGTQ::InvertedIndexObject<LOCAL_ID_TYPE> &entry = (*invertedIndex[gid]).at(idx, invertedIndex.allocator);
#else
        NGTQ::InvertedIndexObject<LOCAL_ID_TYPE> &entry = (*invertedIndex[gid])[idx];
#endif
	ii[gid].push_back(entry.id);
      }
    }
  }

  inline NGT::Distance getApproximateDistance(NGT::Object &query, uint32_t globalID, LOCAL_ID_TYPE *localID, QuantizedObjectDistance::DistanceLookupTable &distanceLUT) {
    double distance;
      distance = (*quantizedObjectDistance)(query, globalID, localID, distanceLUT);

    return distance;
  }


   inline void aggregateObjectsWithCache(NGT::ObjectDistance &globalCentroid, NGT::Object *query, size_t size, NGT::ObjectSpace::ResultSet &results, size_t approximateSearchSize) {

     QuantizedObjectDistance::DistanceLookupTable distanceLUT;
     (*quantizedObjectDistance).initialize(distanceLUT);

     for (size_t j = 0; j < invertedIndex[globalCentroid.id]->size() && results.size() < approximateSearchSize; j++) {
#ifdef NGTQ_SHARED_INVERTED_INDEX
       InvertedIndexObject<LOCAL_ID_TYPE> &invertedIndexEntry = (*invertedIndex[globalCentroid.id]).at(j, invertedIndex.allocator);
#else
       InvertedIndexObject<LOCAL_ID_TYPE> &invertedIndexEntry = (*invertedIndex[globalCentroid.id])[j];
#endif
       double distance;
       if (invertedIndexEntry.localID[0] == 0) {
	 distance = globalCentroid.distance;
       } else { 
	 distance = getApproximateDistance(*query, globalCentroid.id, invertedIndexEntry.localID, distanceLUT);
       }  

       NGT::ObjectDistance obj;
       obj.id = invertedIndexEntry.id;
       obj.distance = distance;
       assert(obj.id > 0);
       results.push(obj);

     } 
  }


  inline void aggregateObjects(NGT::ObjectDistance &globalCentroid, NGT::Object *query, size_t size, NGT::ObjectSpace::ResultSet &results, size_t approximateSearchSize) {
    for (size_t j = 0; j < invertedIndex[globalCentroid.id]->size() && results.size() < approximateSearchSize; j++) {
#ifdef NGTQ_SHARED_INVERTED_INDEX
      InvertedIndexObject<LOCAL_ID_TYPE> &invertedIndexEntry = (*invertedIndex[globalCentroid.id]).at(j, invertedIndex.allocator);
#else
      InvertedIndexObject<LOCAL_ID_TYPE> &invertedIndexEntry = (*invertedIndex[globalCentroid.id])[j];
#endif
      double distance;
      if (invertedIndexEntry.localID[0] == 0) {
	distance = globalCentroid.distance;
      } else { 
	distance = (*quantizedObjectDistance)(*query, globalCentroid.id, invertedIndexEntry.localID);
      }  


      NGT::ObjectDistance obj;
      obj.id = invertedIndexEntry.id;
      obj.distance = distance;
      assert(obj.id > 0);
      results.push(obj);
      if (results.size() >= approximateSearchSize) {
	return;
      }

    } 
  }


  inline void aggregateObjects(NGT::Object *query, size_t size, NGT::ObjectDistances &objects, NGT::ObjectSpace::ResultSet &results, size_t approximateSearchSize, AggregateObjectsFunction aggregateObjectsFunction) {
    for (size_t i = 0; i < objects.size(); i++) {
      if (invertedIndex[objects[i].id] == 0) {
	if (property.centroidCreationMode == CentroidCreationModeDynamic) {
	  cerr << "Inverted index is empty. " << objects[i].id << endl;
	}
	continue;
      }
      ((*this).*aggregateObjectsFunction)(objects[i], query, size, results, results.size() == 0 ? INT_MAX : approximateSearchSize);
      if (results.size() >= approximateSearchSize) {
	return;
      }
    } 
  }

  void refineDistance(NGT::Object *query, NGT::ObjectDistances &results) {
     NGT::ObjectSpace &objectSpace = globalCodebook.getObjectSpace();
     for (auto i = results.begin(); i != results.end(); ++i) {
       NGT::ObjectDistance &result = *i;
       NGT::Object o(&objectSpace);
       objectList.get(result.id, (NGT::Object&)o, &objectSpace);
       double distance = objectSpace.getComparator()(*query, (NGT::Object&)o);
       result.distance = distance;
     }
     std::sort(results.begin(), results.end());
  }

  void search(NGT::Object *query, NGT::ObjectDistances &objs, 
	      size_t size, 
      	      float expansion,
	      AggregationMode aggregationMode,
	      double epsilon = FLT_MAX) {
    size_t approximateSearchSize = size * expansion;
    size_t codebookSearchSize = approximateSearchSize / (objectList.size() / globalCodebook.getObjectRepositorySize()) + 1;
    search(query, objs, size, approximateSearchSize, codebookSearchSize, aggregationMode, epsilon);
  }

  void search(NGT::Object *query, NGT::ObjectDistances &objs, 
	      size_t size, size_t approximateSearchSize,
	      size_t codebookSearchSize, bool resultRefinement,
	      bool lookUpTable = false,
	      double epsilon = FLT_MAX) {
    AggregationMode aggregationMode;
    if (resultRefinement) {
      aggregationMode = AggregationModeExactDistance;
    } else {
      if (lookUpTable) {
	aggregationMode = AggregationModeApproximateDistanceWithLookupTable;
      } else {
	aggregationMode = AggregationModeApproximateDistanceWithCache;
      }
    }
    search(query, objs, size, approximateSearchSize, codebookSearchSize, aggregationMode, epsilon);
  }

  void search(NGT::Object *query, NGT::ObjectDistances &objs, 
	      size_t size, size_t approximateSearchSize,
	      size_t codebookSearchSize, 
	      AggregationMode aggregationMode,
	      double epsilon = FLT_MAX) {
    if (aggregationMode == AggregationModeApproximateDistanceWithLookupTable) {
      if (property.dataType != DataTypeFloat) {
	NGTThrowException("NGTQ: Fatal inner error. the lookup table is only for dataType float!");
      }
    }
    NGT::ObjectDistances objects;
    searchGlobalCodebook(query, size, objects, approximateSearchSize, codebookSearchSize, epsilon);

    objs.clear();
    NGT::ObjectSpace::ResultSet results;
    distanceComputationCount = 0;

    AggregateObjectsFunction aggregateObjectsFunction = &QuantizerInstance::aggregateObjectsWithCache;
    switch(aggregationMode) {
    case AggregationModeExactDistance :
      aggregateObjectsFunction = &QuantizerInstance::aggregateObjectsWithExactDistance;
      break;
    case AggregationModeApproximateDistanceWithLookupTable :
      aggregateObjectsFunction = &QuantizerInstance::aggregateObjectsWithLookupTable;
      break;
    case AggregationModeExactDistanceThroughApproximateDistance :
    case AggregationModeApproximateDistanceWithCache :
      aggregateObjectsFunction = &QuantizerInstance::aggregateObjectsWithCache;
      break;
    case AggregationModeApproximateDistance :
      aggregateObjectsFunction = &QuantizerInstance::aggregateObjects;
      break;
    default:
      cerr << "NGTQ::Fatal Error. invalid aggregation mode. " << aggregationMode << endl;
      abort();
    }

    aggregateObjects(query, size, objects, results, approximateSearchSize, aggregateObjectsFunction);

    objs.resize(results.size());
    while (!results.empty()) {
      objs[results.size() - 1] = results.top();
      results.pop();
    }
    if (objs.size() > size) {
      objs.resize(size);
    }
    if (aggregationMode == AggregationModeExactDistanceThroughApproximateDistance) {
      refineDistance(query, objs);
    }
  }


  double calculateQuantizationError() {
    NGT::ObjectSpace &objectSpace = globalCodebook.getObjectSpace();
    double distance = 0.0;
    double globalDistance = 0.0;
    size_t count = 0;
    for (size_t gi = 0; gi < invertedIndex.size(); gi++) {
      if (invertedIndex[gi] != 0) {
	NGT::PersistentObject &gcentroid = *globalCodebook.getObjectSpace().getRepository().get(gi);
	for (size_t li = 0; li < invertedIndex[gi]->size(); li++) {
#ifdef NGTQ_SHARED_INVERTED_INDEX
          size_t id = invertedIndex[gi]->at(li, invertedIndex.allocator).id;
#else
	  size_t id = invertedIndex[gi]->at(li).id;
#endif
	  NGT::Object object(&objectSpace);
	  objectList.get(id, (NGT::Object&)object, &objectSpace);
#ifdef NGTQ_SHARED_INVERTED_INDEX
	  double d = (*quantizedObjectDistance)(object, gi, invertedIndex[gi]->at(li, invertedIndex.allocator).localID);
#else
	  double d = (*quantizedObjectDistance)(object, gi, invertedIndex[gi]->at(li).localID);
#endif
	  distance += d;
	  count++;
	  NGT::Distance gd = globalCodebook.getObjectSpace().getComparator()(object, gcentroid);
	  globalDistance += gd;
	}
      }
    }
    distance /= count;
    globalDistance /= count;
    std::cerr << distance << ":" << globalDistance << std::endl;
    return distance;
  }

  void info(ostream &os, char mode) {
    cerr << "info" << endl;
    os << "Inverted index size=" << invertedIndex.size() << endl;
    for (size_t i = 0; i < invertedIndex.size(); i++) {
      if (invertedIndex[i] != 0) {
	os << i << " " << invertedIndex[i]->size();
	if (mode == 'a' || mode == 'l') {
	  os << ": ";
	  for (size_t li = 0; li < invertedIndex[i]->size(); li++) {
#ifdef NGTQ_SHARED_INVERTED_INDEX
            os << invertedIndex[i]->at(li, invertedIndex.allocator).id << " ";
#else
	    os << invertedIndex[i]->at(li).id << " ";
#endif
	  }
	}
	os << endl;
      }
    }
    if (mode == 'a' || mode == 'e') {
      os << "Quantization Error=" << calculateQuantizationError() << endl;
    }
  }

  NGT::Index &getLocalCodebook(size_t idx) { return localCodebook[idx]; }

  void verify() {
    cerr << "sizeof(LOCAL_ID_TYPE)=" << sizeof(LOCAL_ID_TYPE) << endl;
    size_t objcount = objectList.size();
    cerr << "Object count=" << objcount << endl;
    size_t gcount = globalCodebook.getObjectRepositorySize();
    cerr << "Global codebook size=" << gcount << endl;
    size_t lcount = localCodebook[0].getObjectRepositorySize();
    cerr << "Local codebook size=" << lcount << endl;
    lcount *= 1.1;
    cerr << "Inverted index size=" << invertedIndex.size() << endl;

    cerr << "Started verifying global codebook..." << endl;
    vector<uint8_t> status;
    globalCodebook.verify(status);

    cerr << "Started verifing the inverted index." << endl;
    size_t errorCount = 0;
    for (size_t i = 1; i < invertedIndex.size(); i++) {
      if (i % 1000000 == 0) {
	cerr << "  verified " << i << " entries" << endl;
      }
      if (errorCount > 100) {
	cerr << "Too many errors. Stop..." << endl;
	return;
      }
      if (invertedIndex[i] == 0) {
	cerr << "Warning inverted index is zero. " << i << endl;
	continue;
      }
      for (size_t j = 1; j < invertedIndex[i]->size(); j++) {
#ifdef NGTQ_SHARED_INVERTED_INDEX
	InvertedIndexObject<LOCAL_ID_TYPE> &invertedIndexEntry = (*invertedIndex[i]).at(j, invertedIndex.allocator);
#else
	InvertedIndexObject<LOCAL_ID_TYPE> &invertedIndexEntry = (*invertedIndex[i])[j];
#endif
	if (invertedIndexEntry.id >= objcount) {
	  cerr << "The object ID of the inverted index entry is too big! " << invertedIndexEntry.id << ":" << objcount << endl;
	  cerr << "  No. of the wrong entry in the inverted index is " << i << endl;
	  errorCount++;
	}
	if (invertedIndexEntry.id == 0) {
	  cerr << "The object ID of the inverted index entry is zero! " << invertedIndexEntry.id << ":" << objcount << endl;
	  cerr << "  No. of the wrong entry in the inverted index is " << i << endl;
	  errorCount++;
	}
	for (size_t li = 0; li < property.localDivisionNo; li++) {
	  if (lcount != 0 && invertedIndexEntry.localID[li] >= lcount) {
	    cerr << "The local centroid ID of the inverted index entry is wrong. " << invertedIndexEntry.localID[li] << ":" << lcount << endl;
	    cerr << "  No. of the wrong entry in the inverted index is " << i << ". No. of local centroid db is " << li << endl;
	    errorCount++;
	  }
	  if (invertedIndexEntry.localID[li] == 0) {
	  }
	}
      }
    }
  }

  size_t getLocalCodebookSize(size_t size) { return localCodebook[size].getObjectRepositorySize(); }

  size_t getInstanceSharedMemorySize(ostream &os, SharedMemoryAllocator::GetMemorySizeType t = SharedMemoryAllocator::GetTotalMemorySize) {
#ifdef NGTQ_SHARED_INVERTED_INDEX
    size_t size = invertedIndex.getAllocator().getMemorySize(t);
#else
    size_t size = 0;
#endif
    os << "inverted=" << size << endl;
    os << "Local centroid:" << endl;
    for (size_t di = 0; di < localCodebook.size(); di++) {
      size += localCodebook[di].getSharedMemorySize(os, t);
    }
    return size;
  }

  size_t getNumberOfObjects(NGT::GraphAndTreeIndex &index) {
    return index.getObjectRepositorySize() == 0 ? 0 : static_cast<int>(index.getObjectRepositorySize()) - 1;
  }

  QuantizedObjectDistance &getQuantizedObjectDistance() { return *quantizedObjectDistance; }

#ifdef NGTQ_SHARED_INVERTED_INDEX
  NGT::PersistentRepository<IIEntry>	invertedIndex;
#else
  NGT::Repository<IIEntry>	invertedIndex;
#endif
  QuantizedObjectDistance	*quantizedObjectDistance;
  GenerateResidualObject	*generateResidualObject;
  std::vector<NGT::Index>	localCodebook;

};
  
class Quantization {
public:
  static Quantizer *generate(Property &property) {
    DataType dataType      = property.dataType;
    size_t dimension       = property.dimension;
    size_t divisionNo      = property.getLocalCodebookNo();
    size_t localIDByteSize = property.localIDByteSize;
    Quantizer *quantizer = 0;
    if (property.centroidCreationMode == CentroidCreationModeNone) {
      NGTThrowException("Centroid creation mode is not specified");
    } else {
      if (localIDByteSize == 4) {
	quantizer = new QuantizerInstance<uint32_t>(dataType, dimension, divisionNo);
      } else if (localIDByteSize == 2) {
	quantizer = new QuantizerInstance<uint16_t>(dataType, dimension, divisionNo);
      } else {
	std::stringstream msg;
	msg << "Not support the specified size of local ID. " << localIDByteSize;
	NGTThrowException(msg);
      }
    }

    return quantizer;
  }
};

 class Index {
 public:
   Index():quantizer(0) {}
   Index(const string& index, bool rdOnly = false):quantizer(0) { open(index, rdOnly); }
   ~Index() { close(); }



   static void create(const string &index, Property &property, 
		      NGT::Property &globalProperty,
		      NGT::Property &localProperty) {
     if (property.dimension == 0) {
       NGTThrowException("NGTQ::create: Error. The dimension is zero.");
     }
     property.setupLocalIDByteSize();
     NGTQ::Quantizer *quantizer = 
       NGTQ::Quantization::generate(property);
     try {
       quantizer->property.setup(property);
       quantizer->create(index, globalProperty, localProperty);
     } catch(NGT::Exception &err) {
       delete quantizer;
       throw err;
     }
     delete quantizer;
   }
#ifdef NGTQ_SHARED_INVERTED_INDEX
   static void compress(const string &indexFile) {
     Index index;
     index.open(indexFile);
     string tmpivt = indexFile + "/ivt-tmp";
     index.getQuantizer().reconstructInvertedIndex(tmpivt);
     index.close();
     string ivt = indexFile + "/ivt";
     unlink(ivt.c_str());
     rename(tmpivt.c_str(), ivt.c_str());
     string ivtc = ivt + "c";
     unlink(ivtc.c_str());
     string tmpivtc = tmpivt + "c";
     rename(tmpivtc.c_str(), ivtc.c_str());
   }
#endif

  static void append(const string &indexName,	// index file
		     const string &data,	// data file
		     size_t dataSize = 0	// data size
		     ) {
    NGTQ::Index index(indexName);
    istream *is;
    if (data == "-") {
      is = &cin;
    } else {
      ifstream *ifs = new ifstream;
      ifs->ifstream::open(data);
      if (!(*ifs)) {
	cerr << "Cannot open the specified file. " << data << endl;
	return;
      }
      is = ifs;
    }
    string line;
    vector<pair<NGT::Object*, size_t> > objects;
    size_t count = 0;
    // extract objects from the file and insert them to the object list.
    while(getline(*is, line)) {
      count++;
      index.insert(line, objects, 0);
      if (count % 10000 == 0) {
	  cerr << "Processed " << count;
	  cerr << endl;
      }
    }
    if (objects.size() > 0) {
      index.insert(objects);
    }
    cerr << "end of insertion. " << count << endl;
    if (data != "-") {
      delete is;
    }

    index.save();
    index.close();
  }

  static void rebuild(const string &indexName,		
		      const string &rebuiltIndexName	
		     ) {

    const string srcObjectList = indexName + "/obj";
    const string dstObjectList = rebuiltIndexName + "/obj";

    if (std::rename(srcObjectList.c_str(), dstObjectList.c_str()) != 0) {
      stringstream msg;
      msg << "Quantizer::rebuild: Cannot rename an object file. " << srcObjectList << "=>" << dstObjectList ;
      NGTThrowException(msg);
    }

    try {
      NGTQ::Index index(rebuiltIndexName);
    
      index.rebuildIndex();

      index.save();
      index.close();
    } catch(NGT::Exception &err) {
      std::rename(dstObjectList.c_str(), srcObjectList.c_str());
      throw err;
    }

  }

  void open(const string &index, bool readOnly = false) {
     close();
     NGT::Property globalProperty;
     globalProperty.clear();
     globalProperty.edgeSizeForSearch = 40;
     quantizer = getQuantizer(index, globalProperty);
     if (readOnly) {
       quantizer->closeCodebooks();
     }
   }

   void save() {
     getQuantizer().save();
   }

   void close() {
     if (quantizer != 0) {
       delete quantizer;
       quantizer = 0;
     }
   }
   void insert(string &line, vector<pair<NGT::Object*, size_t> > &objects, size_t id) {
     getQuantizer().insert(line, objects, id);
   }

   void insert(vector<pair<NGT::Object*, size_t> > &objects) {
     getQuantizer().insert(objects);
   }

   void rebuildIndex() {
     getQuantizer().rebuildIndex();
   }

   NGT::Object *allocateObject(string &line, const string &sep, size_t dimension) {
     return getQuantizer().allocateObject(line, sep);
   }

   NGT::Object *allocateObject(vector<double> &obj) {
     return getQuantizer().allocateObject(obj);
   }

   void deleteObject(NGT::Object *object) { getQuantizer().deleteObject(object); }

   void search(NGT::Object *object, NGT::ObjectDistances &objs, 
	       size_t size, size_t approximateSearchSize,
	       size_t codebookSearchSize, bool resultRefinement, 
	       bool lookUpTable, double epsilon) {
     getQuantizer().search(object, objs, size, approximateSearchSize, codebookSearchSize, 
			   resultRefinement, lookUpTable, epsilon);
   }

   void search(NGT::Object *object, NGT::ObjectDistances &objs, 
	       size_t size, float expansion,
	       AggregationMode aggregationMode,
	       double epsilon) {
     getQuantizer().search(object, objs, size, expansion, 
			   aggregationMode, epsilon);
   }

   void info(ostream &os, char mode) { getQuantizer().info(os, mode); }

   void verify() { getQuantizer().verify(); }

   NGTQ::Quantizer &getQuantizer() { 
     if (quantizer == 0) {
       NGTThrowException("NGTQ::Index: Not open.");
     }
     return *quantizer; 
   }

   size_t getGlobalCodebookSize() { return quantizer->globalCodebook.getObjectRepositorySize(); }
   size_t getLocalCodebookSize(size_t idx) { return quantizer->getLocalCodebookSize(idx); }

   size_t getSharedMemorySize(ostream &os, SharedMemoryAllocator::GetMemorySizeType t = SharedMemoryAllocator::GetTotalMemorySize) {
     return quantizer->getSharedMemorySize(os, t);
   }

 protected:

   static NGTQ::Quantizer *getQuantizer(const string &index) {
     NGT::Property globalProperty;
     globalProperty.clear();
     return getQuantizer(index, globalProperty);
   }

   static NGTQ::Quantizer *getQuantizer(const string &index, NGT::Property &globalProperty) {
     NGTQ::Property property;
     try {
       property.load(index);
     } catch (NGT::Exception &err) {
       stringstream msg;
       msg << "Quantizer::getQuantizer: Cannot load the property. " << index << " : " << err.what();
       NGTThrowException(msg);
     }
     NGTQ::Quantizer *quantizer = NGTQ::Quantization::generate(property);
     if (quantizer == 0) {
       NGTThrowException("NGTQ::Index: Cannot get quantizer.");
     }
     try {
       quantizer->open(index, globalProperty);
     } catch(NGT::Exception &err) {
       delete quantizer;
       throw err;
     }
     return quantizer;
   }

   NGTQ::Quantizer *quantizer;
 };

} // namespace NGTQ
