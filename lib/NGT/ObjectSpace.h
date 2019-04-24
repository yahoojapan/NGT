//
// Copyright (C) 2015-2019 Yahoo Japan Corporation
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

#include "PrimitiveComparator.h"

class ObjectSpace;

namespace NGT {

  class PersistentObjectDistances;
  class ObjectDistances : public vector<ObjectDistance> {
  public:
    ObjectDistances(NGT::ObjectSpace *os = 0) {}
    void serialize(ofstream &os, ObjectSpace *objspace = 0) { NGT::Serializer::write(os, (vector<ObjectDistance>&)*this);}
    void deserialize(ifstream &is, ObjectSpace *objspace = 0) { NGT::Serializer::read(is, (vector<ObjectDistance>&)*this);}

    void serializeAsText(ofstream &os, ObjectSpace *objspace = 0) { 
      NGT::Serializer::writeAsText(os, size());
      os << " ";
      for (size_t i = 0; i < size(); i++) {
	(*this)[i].serializeAsText(os);
	os << " ";
      }
    }
    void deserializeAsText(ifstream &is, ObjectSpace *objspace = 0) {
      size_t s;
      NGT::Serializer::readAsText(is, s);     
      resize(s);
      for (size_t i = 0; i < size(); i++) {
	(*this)[i].deserializeAsText(is);
      }
    }

    void moveFrom(priority_queue<ObjectDistance, vector<ObjectDistance>, less<ObjectDistance> > &pq) {
      this->clear();
      this->resize(pq.size());
      for (int i = pq.size() - 1; i >= 0; i--) {
	(*this)[i] = pq.top();
	pq.pop();
      }
      assert(pq.size() == 0);
    }

    void moveFrom(priority_queue<ObjectDistance, vector<ObjectDistance>, less<ObjectDistance> > &pq, double (&f)(double)) {
      this->clear();
      this->resize(pq.size());
      for (int i = pq.size() - 1; i >= 0; i--) {
	(*this)[i] = pq.top();
	(*this)[i].distance = f((*this)[i].distance);
	pq.pop();
      }
      assert(pq.size() == 0);
    }

    void moveFrom(priority_queue<ObjectDistance, vector<ObjectDistance>, less<ObjectDistance> > &pq, unsigned int id) {
      this->clear();
      if (pq.size() == 0) {
	return;
      }
      this->resize(id == 0 ? pq.size() : pq.size() - 1);
      int i = this->size() - 1;
      while (pq.size() != 0 && i >= 0) {
	if (pq.top().id != id) {
	  (*this)[i] = pq.top();
	  i--;
	}
	pq.pop();
      }
      if (pq.size() != 0 && pq.top().id != id) {
	cerr << "moveFrom: Fatal error: somethig wrong! " << pq.size() << ":" << this->size() << ":" << id << ":" << pq.top().id << endl;
	assert(pq.size() == 0 || pq.top().id == id);
      }
    }

    ObjectDistances &operator=(PersistentObjectDistances &objs);
  };

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  class PersistentObjectDistances : public Vector<ObjectDistance> {
  public:
    PersistentObjectDistances(SharedMemoryAllocator &allocator, NGT::ObjectSpace *os = 0) {}
    void serialize(ofstream &os, ObjectSpace *objectspace = 0) { NGT::Serializer::write(os, (Vector<ObjectDistance>&)*this); }
    void deserialize(ifstream &is, ObjectSpace *objectspace = 0) { NGT::Serializer::read(is, (Vector<ObjectDistance>&)*this); }
    void serializeAsText(ofstream &os, SharedMemoryAllocator &allocator, ObjectSpace *objspace = 0) { 
      NGT::Serializer::writeAsText(os, size());
      os << " ";
      for (size_t i = 0; i < size(); i++) {
	(*this).at(i, allocator).serializeAsText(os);
	os << " ";
      }
    }
    void deserializeAsText(ifstream &is, SharedMemoryAllocator &allocator, ObjectSpace *objspace = 0) {
      size_t s;
      is >> s;
      resize(s, allocator);
      for (size_t i = 0; i < size(); i++) {
	(*this).at(i, allocator).deserializeAsText(is);
      }
    }
    PersistentObjectDistances &copy(ObjectDistances &objs, SharedMemoryAllocator &allocator) {
      clear(allocator);
      reserve(objs.size(), allocator);
      for (ObjectDistances::iterator i = objs.begin(); i != objs.end(); i++) {
	push_back(*i, allocator);
      }
      return *this;
    }
  };
  typedef PersistentObjectDistances	GraphNode;

  inline ObjectDistances &ObjectDistances::operator=(PersistentObjectDistances &objs)
    {
      clear();
      reserve(objs.size());
      cerr << "not implemented" << endl;
      assert(0);
      return *this;
    }
#else // NGT_SHARED_MEMORY_ALLOCATOR
  typedef ObjectDistances	GraphNode;
#endif // NGT_SHARED_MEMORY_ALLOCATOR

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  class PersistentObject;
#else
  typedef Object	PersistentObject;
#endif

  class ObjectRepository;

  class ObjectSpace {
  public:
    class Comparator {
    public:
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    Comparator(size_t d, SharedMemoryAllocator &a) : dimension(d), allocator(a) {}
#else
    Comparator(size_t d) : dimension(d) {}
#endif
      virtual double operator()(Object &objecta, Object &objectb) = 0;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      virtual double operator()(Object &objecta, PersistentObject &objectb) = 0;
      virtual double operator()(PersistentObject &objecta, PersistentObject &objectb) = 0;
#endif
      size_t dimension;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      SharedMemoryAllocator &allocator;
#endif
      virtual ~Comparator(){}
    };
    enum DistanceType {
      DistanceTypeNone			= -1,
      DistanceTypeL1			= 0,
      DistanceTypeL2			= 1,
      DistanceTypeHamming		= 2,
      DistanceTypeAngle			= 3,
      DistanceTypeCosine		= 4,
      DistanceTypeNormalizedAngle	= 5,
      DistanceTypeNormalizedCosine	= 6
    };

    enum ObjectType {
      ObjectTypeNone	= 0,
      Uint8		= 1,
      Float		= 2
    };

    typedef priority_queue<ObjectDistance, vector<ObjectDistance>, less<ObjectDistance> > ResultSet;
    ObjectSpace(size_t d):dimension(d), paddedDimension(((d - 1) / 4 + 1) * 4), distanceType(DistanceTypeNone), comparator(0), normalization(false) {}
    virtual ~ObjectSpace() { if (comparator != 0) { delete comparator; } }
    
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    virtual void open(const string &f, size_t shareMemorySize) = 0;
    virtual Object *allocateObject(Object &o) = 0;
    virtual Object *allocateObject(PersistentObject &o) = 0;
    virtual PersistentObject *allocatePersistentObject(Object &obj) = 0;
    virtual void deleteObject(PersistentObject *) = 0;
    virtual void copy(PersistentObject &objecta, PersistentObject &objectb) = 0;
    virtual void show(ostream &os, PersistentObject &object) = 0;
    virtual size_t insert(PersistentObject *obj) = 0;
#else
    virtual size_t insert(Object *obj) = 0;
#endif

    Comparator &getComparator() { return *comparator; }

    virtual void serialize(const string &of) = 0;
    virtual void deserialize(const string &ifile) = 0;
    virtual void serializeAsText(const string &of) = 0;
    virtual void deserializeAsText(const string &of) = 0;
    virtual void readText(istream &is, size_t dataSize) = 0;
    virtual void appendText(istream &is, size_t dataSize) = 0;
    virtual void append(const float *data, size_t dataSize) = 0;
    virtual void append(const double *data, size_t dataSize) = 0;

    virtual void copy(Object &objecta, Object &objectb) = 0;

    virtual void linearSearch(Object &query, double radius, size_t size,  
			      ObjectSpace::ResultSet &results) = 0;

    virtual const std::type_info &getObjectType() = 0;
    virtual void show(ostream &os, Object &object) = 0;
    virtual size_t getSize() = 0;
    virtual size_t getSizeOfElement() = 0;
    virtual size_t getByteSizeOfObject() = 0;
    virtual Object *allocateNormalizedObject(const string &textLine, const string &sep) = 0;
    virtual Object *allocateNormalizedObject(const vector<double> &obj) = 0;
    virtual Object *allocateNormalizedObject(const vector<float> &obj) = 0;
    virtual Object *allocateNormalizedObject(const float *obj, size_t size) = 0;
    virtual PersistentObject *allocateNormalizedPersistentObject(const vector<double> &obj) = 0;
    virtual PersistentObject *allocateNormalizedPersistentObject(const vector<float> &obj) = 0;
    virtual void deleteObject(Object *po) = 0;
    virtual Object *allocateObject() = 0;
    virtual void remove(size_t id) = 0;

    virtual ObjectRepository &getRepository() = 0;

    virtual void setDistanceType(DistanceType t) = 0;

    virtual void *getObject(size_t idx) = 0;
    virtual void getObject(size_t idx, vector<float> &v) = 0;
    virtual void getObjects(const vector<size_t> &idxs, vector<vector<float>> &vs) = 0;

    size_t getDimension() { return dimension; }

    size_t getPaddedDimension() { return paddedDimension; }

    template <typename T>
    void normalize(T *data, size_t dim) {
      double sum = 0.0;
      for (size_t i = 0; i < dim; i++) {
	sum += (double)data[i] * (double)data[i];
      }
      if (sum == 0.0) {
	stringstream msg;
	msg << "ObjectSpace::normalize: Error! the object is an invalid zero vector for the cosine similarity or angle distance.";
	NGTThrowException(msg);
      }
      sum = sqrt(sum);
      for (size_t i = 0; i < dim; i++) {
	data[i] = (double)data[i] / sum;
      }
    }
    uint32_t getPrefetchOffset() { return prefetchOffset; }
    uint32_t setPrefetchOffset(size_t offset) {
      if (offset == 0) {
	prefetchOffset = floor(300.0 / (static_cast<float>(getPaddedDimension()) + 30.0) + 1.0);
      } else {
	prefetchOffset = offset;
      }
      return prefetchOffset;
    }
  protected:
    const size_t	dimension;
    const size_t	paddedDimension;
    DistanceType	distanceType;
    Comparator		*comparator;
    bool		normalization;
    uint32_t		prefetchOffset;
  };

  class BaseObject {
  public:
    virtual uint8_t &operator[](size_t idx) const = 0;
    void serialize(ostream &os, ObjectSpace *objectspace = 0) { 
      assert(objectspace != 0);
      size_t byteSize = objectspace->getByteSizeOfObject();
      NGT::Serializer::write(os, (uint8_t*)&(*this)[0], byteSize); 
    }
    void deserialize(istream &is, ObjectSpace *objectspace = 0) { 
      assert(objectspace != 0);
      size_t byteSize = objectspace->getByteSizeOfObject();
      assert(&(*this)[0] != 0);
      NGT::Serializer::read(is, (uint8_t*)&(*this)[0], byteSize); 
    }
    void serializeAsText(ostream &os, ObjectSpace *objectspace = 0) { 
      assert(objectspace != 0);
      const std::type_info &t = objectspace->getObjectType();
      size_t dimension = objectspace->getDimension();
      void *ref = (void*)&(*this)[0];
      if (t == typeid(uint8_t)) {
	NGT::Serializer::writeAsText(os, (uint8_t*)ref, dimension); 
      } else if (t == typeid(float)) {
	NGT::Serializer::writeAsText(os, (float*)ref, dimension); 
      } else if (t == typeid(double)) {
	NGT::Serializer::writeAsText(os, (double*)ref, dimension); 
      } else if (t == typeid(uint16_t)) {
	NGT::Serializer::writeAsText(os, (uint16_t*)ref, dimension); 
      } else if (t == typeid(uint32_t)) {
	NGT::Serializer::writeAsText(os, (uint32_t*)ref, dimension); 
      } else {
	cerr << "Object::serializeAsText: not supported data type. [" << t.name() << "]" << endl;
	assert(0);
      }
    }
    void deserializeAsText(ifstream &is, ObjectSpace *objectspace = 0) {
      assert(objectspace != 0);
      const std::type_info &t = objectspace->getObjectType();
      size_t dimension = objectspace->getDimension();
      void *ref = (void*)&(*this)[0];
      assert(ref != 0);
      if (t == typeid(uint8_t)) {
	NGT::Serializer::readAsText(is, (uint8_t*)ref, dimension); 
      } else if (t == typeid(float)) {
	NGT::Serializer::readAsText(is, (float*)ref, dimension); 
      } else if (t == typeid(double)) {
	NGT::Serializer::readAsText(is, (double*)ref, dimension); 
      } else if (t == typeid(uint16_t)) {
	NGT::Serializer::readAsText(is, (uint16_t*)ref, dimension); 
      } else if (t == typeid(uint32_t)) {
	NGT::Serializer::readAsText(is, (uint32_t*)ref, dimension); 
      } else {
	cerr << "Object::deserializeAsText: not supported data type. [" << t.name() << "]" << endl;
	assert(0);
      }
    }

  };

  class Object : public BaseObject {
  public:
    Object(NGT::ObjectSpace *os = 0):vector(0) {
      assert(os != 0);
      size_t s = os->getByteSizeOfObject();
      construct(s);
    }

    Object(size_t s):vector(0) {
      assert(s != 0);
      construct(s);
    }

    void copy(Object &o, size_t s) {
      assert(vector != 0);
      for (size_t i = 0; i < s; i++) {
	vector[i] = o[i];
      }
    }

    virtual ~Object() { clear(); }

    uint8_t &operator[](size_t idx) const { return vector[idx]; }

    void *getPointer(size_t idx = 0) const { return vector + idx; }

    static Object *allocate(ObjectSpace &objectspace) { return new Object(&objectspace); }
  private:
    void clear() {
      if (vector != 0) {
	MemoryCache::alignedFree(vector);
      }
      vector = 0;
    }

    void construct(size_t s) {
      assert(vector == 0);
      size_t allocsize = ((s - 1) / 16 + 1) * 16;	
      vector = static_cast<uint8_t*>(MemoryCache::alignedAlloc(allocsize));
      memset(vector, 0, allocsize);
    }

    uint8_t* vector;
  };


#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  class PersistentObject : public BaseObject {
  public:
    PersistentObject(SharedMemoryAllocator &allocator, NGT::ObjectSpace *os = 0):array(0) {
      assert(os != 0);
      size_t s = os->getByteSizeOfObject();
      construct(s, allocator);
    }
    PersistentObject(SharedMemoryAllocator &allocator, size_t s):array(0) {
      assert(s != 0);
      construct(s, allocator);
    }

    virtual ~PersistentObject() {}

    uint8_t &at(size_t idx, SharedMemoryAllocator &allocator) const { 
      uint8_t *a = (uint8_t *)allocator.getAddr(array);
      return a[idx];
    }
    uint8_t &operator[](size_t idx) const {
      cerr << "not implemented" << endl;
      assert(0);
      uint8_t *a = 0;
      return a[idx];
    }

    void *getPointer(size_t idx, SharedMemoryAllocator &allocator) {
      uint8_t *a = (uint8_t *)allocator.getAddr(array);
      return a + idx; 
    }

    // set v in objectspace to this object using allocator.
    void set(PersistentObject &po, ObjectSpace &objectspace);

    static off_t allocate(ObjectSpace &objectspace);

    void serializeAsText(ostream &os, SharedMemoryAllocator &allocator, 
			 ObjectSpace *objectspace = 0) { 
      serializeAsText(os, objectspace);
    }

    void serializeAsText(ostream &os, ObjectSpace *objectspace = 0);

    void deserializeAsText(ifstream &is, SharedMemoryAllocator &allocator, 
			   ObjectSpace *objectspace = 0) {
      deserializeAsText(is, objectspace);
    }

    void deserializeAsText(ifstream &is, ObjectSpace *objectspace = 0);

    void serialize(ostream &os, SharedMemoryAllocator &allocator, 
		   ObjectSpace *objectspace = 0) { 
      cerr << "serialize is not implemented" << endl;
      assert(0);
    }

  private:
    void construct(size_t s, SharedMemoryAllocator &allocator) {
      assert(array == 0);
      assert(s != 0);
      size_t allocsize = ((s - 1) / 16 + 1) * 16;
      array = allocator.getOffset(new(allocator) uint8_t[allocsize]);
    }
    off_t array;
  };
#endif // NGT_SHARED_MEMORY_ALLOCATOR

}

