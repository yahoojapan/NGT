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


#if !defined(NGT_AVX_DISABLED) && defined(__AVX__)
#include	<immintrin.h>
#endif

#include	"Common.h"

class ObjectSpace;

namespace NGT {
#pragma pack(2)
  class ObjectDistance {
  public:
    ObjectDistance():id(0), distance(0.0) {}
    ObjectDistance(unsigned int i, float d):id(i), distance(d) {}
    bool operator==(const ObjectDistance &o) const {
      return (distance == o.distance) && (id == o.id);
    }
    bool operator<(const ObjectDistance &o) const {
      if (distance == o.distance) {
	return id < o.id;
      } else {
	return distance < o.distance;
      }
    }
    bool operator>(const ObjectDistance &o) const {
      if (distance == o.distance) {
	return id > o.id;
      } else {
	return distance > o.distance;
      }
    }
    void serialize(ofstream &os) {
      NGT::Serializer::write(os, id);
      NGT::Serializer::write(os, distance);
    }
    void deserialize(ifstream &is) {
      NGT::Serializer::read(is, id);
      NGT::Serializer::read(is, distance);
    }

    void serializeAsText(ofstream &os) {
      os.unsetf(std::ios_base::floatfield);
      os << setprecision(8) << id << " " << distance;
    }

    void deserializeAsText(ifstream &is) {
      is >> id;
      is >> distance;
    }

    friend ostream &operator<<(ostream& os, const ObjectDistance &o) {
      os << o.id << " " << o.distance;
      return os;
    }
    friend istream &operator>>(istream& is, ObjectDistance &o) {
      is >> o.id;
      is >> o.distance;
      return is;
    }
    uint32_t		id;
    float		distance;
  };
#pragma pack()

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
      DistanceTypeNone		= -1,
      DistanceTypeL1		= 0,
      DistanceTypeL2		= 1,
      DistanceTypeHamming	= 2,
      DistanceTypeAngle		= 3,
      DistanceTypeCosine	= 4
    };
    typedef priority_queue<ObjectDistance, vector<ObjectDistance>, less<ObjectDistance> > ResultSet;
    ObjectSpace(size_t d):dimension(d), distanceType(DistanceTypeNone), comparator(0) {}
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
    virtual void appendText(ifstream &is, size_t dataSize) = 0;
    virtual void append(const float *data, size_t dataSize) = 0;
    virtual void copy(Object &objecta, Object &objectb) = 0;

    virtual void linearSearch(Object &query, double radius, size_t size,  
			      ObjectSpace::ResultSet &results) = 0;

    virtual const std::type_info &getObjectType() = 0;
    virtual void show(ostream &os, Object &object) = 0;
    virtual size_t getSize() = 0;
    virtual size_t getSizeOfElement() = 0;
    virtual size_t getByteSizeOfObject() = 0;
    virtual Object *allocateObject(const string &textLine, const string &sep) = 0;
    virtual Object *allocateObject(vector<double> &obj) = 0;
    virtual void deleteObject(Object *po) = 0;
    virtual Object *allocateObject() = 0;
    virtual void remove(size_t id) = 0;

    virtual ObjectRepository &getRepository() = 0;

    virtual void setDistanceType(DistanceType t) = 0;

    virtual void *getObject(size_t idx) = 0;

    size_t getDimension() { return dimension; }

  protected:
    const size_t	dimension;
    DistanceType	distanceType;
    Comparator		*comparator;

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
	delete[] vector;
      }
      vector = 0;
    }

    void construct(size_t s) {
      assert(vector == 0);
      vector = new uint8_t[s];
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

    ~PersistentObject() {}

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
      array = allocator.getOffset(new(allocator) uint8_t[s]);
    }
    off_t array;
  };
#else
  typedef Object	PersistentObject;
#endif // NGT_SHARED_MEMORY_ALLOCATOR

  /////////////////////////////////////////////
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  class ObjectRepository : 
  public PersistentRepository<PersistentObject> {
  public:
    typedef PersistentRepository<PersistentObject>	Parent;
    void open(const string &smfile, size_t sharedMemorySize) { 
      string file = smfile;
      file.append("po");
      Parent::open(file, sharedMemorySize);
    }
#else
  class ObjectRepository : public Repository<Object> {
  public:
    typedef Repository<Object>	Parent;
#endif
    ObjectRepository(size_t dim, const type_info &ot):dimension(dim), type(ot) { }

    void initialize() {
      deleteAll();
      Parent::push_back((PersistentObject*)0);
    }

    void serialize(const string &ofile, ObjectSpace *ospace) { 
      ofstream objs(ofile);
      if (!objs.is_open()) {
	stringstream msg;
	msg << "NGT::ObjectSpace: Cannot open the specified file " << ofile << ".";
	NGTThrowException(msg);
      }
      Parent::serialize(objs, ospace); 
    }

    void deserialize(const string &ifile, ObjectSpace *ospace) { 
      assert(ospace != 0);
      ifstream objs(ifile);
      if (!objs.is_open()) {
	stringstream msg;
	msg << "NGT::ObjectSpace: Cannot open the specified file " << ifile << ".";
	NGTThrowException(msg);
      }
      Parent::deserialize(objs, ospace);
    }

    void serializeAsText(const string &ofile, ObjectSpace *ospace) { 
      ofstream objs(ofile);
      if (!objs.is_open()) {
	stringstream msg;
	msg << "NGT::ObjectSpace: Cannot open the specified file " << ofile << ".";
	NGTThrowException(msg);
      }
      Parent::serializeAsText(objs, ospace); 
    }

    void deserializeAsText(const string &ifile, ObjectSpace *ospace) { 
      ifstream objs(ifile);
      if (!objs.is_open()) {
	stringstream msg;
	msg << "NGT::ObjectSpace: Cannot open the specified file " << ifile << ".";
	NGTThrowException(msg);
      }
      Parent::deserializeAsText(objs, ospace); 
    }

    void readText(istream &is, size_t dataSize = 0) {
      initialize();
      appendText(is, dataSize);
    }

    void appendText(istream &is, size_t dataSize = 0) {
      if (dimension == 0) {
	NGTThrowException("ObjectSpace::readText: Dimension is not specified.");
      }
      size_t prevDataSize = size();
      if (prevDataSize == 0) {
	// First entry should be always a dummy entry.
	// If it is empty, the dummy entry should be inserted.
	push_back((PersistentObject*)0);
      }
      if (dataSize > 0) {
	reserve(size() + dataSize);
      }
      string line;
      size_t lineNo = 0;
      while (getline(is, line)) {
	lineNo++;
	if (dataSize > 0 && (dataSize <= size() - prevDataSize)) {
	  cerr << "The size of data reached the specified size. The remaining data in the file are not inserted. " 
	       << dataSize << endl;
	  break;
	}
	vector<double> object;
	try {
	  extractObjectFromText(line, "\t ", object);
	  push_back((PersistentObject*)allocatePersistentObject(object));
	} catch (Exception &err) {
	  std::cerr << "ObjectSpace::readText: Warning! Invalid line. [" << line << "] Skip the line " << lineNo << " and continue." << std::endl;
	}
      }
    }

    template <typename T>
    void append(T *data, size_t objectCount) {
      if (dimension == 0) {
	NGTThrowException("ObjectSpace::readText: Dimension is not specified.");
      }
      if (size() == 0) {
	// First entry should be always a dummy entry.
	// If it is empty, the dummy entry should be inserted.
	push_back((PersistentObject*)0);
      }
      if (objectCount > 0) {
	reserve(size() + objectCount);
      }
      for (size_t idx = 0; idx < objectCount; idx++, data += dimension) {
	try {
	  push_back((PersistentObject*)allocatePersistentObject(data));
	} catch (Exception &err) {
	  std::cerr << "ObjectSpace::readText: Warning! Invalid data. Skip the data no. " << idx << " and continue." << std::endl;
	}
      }
    }

    Object *allocateObject() {
      return (Object*) new Object(getByteSize());
    }

    // This method is called during search to generate query.
    // Therefor the object is no persistent.
    Object *allocateObject(const string &textLine, const string &sep) {
      vector<double> object;
      extractObjectFromText(textLine, sep, object);
      Object *po = (Object*)allocateObject(object);
      return (Object*)po;
    }

    void extractObjectFromText(const string &textLine, const string &sep, vector<double> &object) {
      object.resize(dimension);
      vector<string> tokens;
      NGT::Common::tokenize(textLine, tokens, sep);
      if (dimension > tokens.size()) {
	stringstream msg;
	msg << "ObjectSpace::allocate: too few dimension. " << tokens.size() << ":" << dimension << ". " 
	    << textLine;
	NGTThrowException(msg);
      }
      size_t idx;
      for (idx = 0; idx < dimension; idx++) {
	if (tokens[idx].size() == 0) {
	  stringstream msg;
	  msg << "ObjectSpace::allocate: too few dimension. " << tokens.size() << ":" 
	      << dimension << ". "  << textLine;
	  NGTThrowException(msg);
        }
	char *e;
	object[idx] = strtod(tokens[idx].c_str(), &e);
	if (*e != 0) {
	  std::cerr << "ObjectSpace::readText: Warning! Not numerical value. [" << e << "]" << std::endl;
	  break;
	}
      }
    }

    template <typename T>
      Object *allocateObject(T *o, size_t size = 0) {
      Object *po = new Object(byteSize);
      if (size != 0 && dimension != size) {
	cerr << "ObjectSpace::allocateObject: Fatal error! dimension is invalid. The indexed objects=" 
	     << dimension << " The specified object=" << size << endl;
	assert(dimension == size);
      }
      void *object = (void*)(&(*po)[0]);
      if (type == typeid(uint8_t)) {
	uint8_t *obj = (uint8_t*)object;
	for (size_t i = 0; i < dimension; i++) {
	  obj[i] = (uint8_t)o[i];
	}
      } else if (type == typeid(float)) {
	float *obj = (float*)object;
	for (size_t i = 0; i < dimension; i++) {
	  obj[i] = (float)o[i];
	}
      } else {
	cerr << "ObjectSpace::allocate: Fatal error: unsupported type!" << endl;
	abort();
      }
      return (Object*)po;
    }

    template <typename T>
      Object *allocateObject(vector<T> &o) {
      return allocateObject(o.data(), o.size());
    }

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    PersistentObject *allocatePersistentObject(Object &o) {
      SharedMemoryAllocator &objectAllocator = getAllocator();
      size_t cpsize = dimension;
      if (type == typeid(uint8_t)) {
	cpsize *= sizeof(uint8_t);
      } else if (type == typeid(float)) {
	cpsize *= sizeof(float);
      } else {
	cerr << "ObjectSpace::allocate: Fatal error: unsupported type!" << endl;
	abort();
      }
      PersistentObject *po = new (objectAllocator) PersistentObject(objectAllocator, byteSize);
      void *dsto = &(*po).at(0, allocator);
      void *srco = &o[0];
      memcpy(dsto, srco, cpsize);
      return po;
    }

    template <typename T>
      PersistentObject *allocatePersistentObject(T *o, size_t size = 0) {
      SharedMemoryAllocator &objectAllocator = getAllocator();
      PersistentObject *po = new (objectAllocator) PersistentObject(objectAllocator, byteSize);
      if (size != 0 && dimension != size) {
	cerr << "ObjectSpace::allocateObject: Fatal error! dimension is invalid. The indexed objects=" 
	     << dimension << " The specified object=" << size << endl;
	assert(dimension == size);
      }
      void *object = (void*)(&(*po).at(0, allocator));
      if (type == typeid(uint8_t)) {
	uint8_t *obj = (uint8_t*)object;
	for (size_t i = 0; i < dimension; i++) {
	  obj[i] = (uint8_t)o[i];
	}
      } else if (type == typeid(float)) {
	float *obj = (float*)object;
	for (size_t i = 0; i < dimension; i++) {
	  obj[i] = (float)o[i];
	}
      } else {
	cerr << "ObjectSpace::allocate: Fatal error: unsupported type!" << endl;
	abort();
      }
      return po;
    }

    template <typename T>
    PersistentObject *allocatePersistentObject(vector<T> &o) {
      return allocatePersistentObject(o.data(), o.size());
    }

#else
    template <typename T>
    PersistentObject *allocatePersistentObject(vector<T> &o) {
      return allocateObject(o);
    }

    template <typename T>
      PersistentObject *allocatePersistentObject(T *o, size_t size = 0) {
      return allocateObject(o, size);
    }
#endif

    void deleteObject(Object *po) {
      delete po;
    }

    private:
    void extractObject(void *object, vector<double> &d) {
      if (type == typeid(uint8_t)) {
	uint8_t *obj = (uint8_t*)object;
	for (size_t i = 0; i < dimension; i++) {
	  d.push_back(obj[i]);
	}
      } else if (type == typeid(float)) {
	float *obj = (float*)object;
	for (size_t i = 0; i < dimension; i++) {
	  d.push_back(obj[i]);
	}
      } else {
	cerr << "ObjectSpace::allocate: Fatal error: unsupported type!" << endl;
	abort();
      }
    }
    public:
    void extractObject(Object *o, vector<double> &d) {
      void *object = (void*)(&(*o)[0]);
      extractObject(object, d);
    }

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    void extractObject(PersistentObject *o, vector<double> &d) {
      SharedMemoryAllocator &objectAllocator = getAllocator();
      void *object = (void*)(&(*o).at(0, objectAllocator));
      extractObject(object, d);
    }
#endif

    void setLength(size_t l) {
      byteSize = l;
    }

    size_t getByteSize() { return byteSize; }
    size_t insert(PersistentObject *obj) { return Parent::insert(obj); }
    size_t byteSize;		// the length of all of elements.
    const size_t dimension;
    const type_info &type;
  };

  template <typename OBJECT_TYPE, typename COMPARE_TYPE> 
    class ObjectSpaceT : public ObjectSpace, public ObjectRepository {
  public:

    class ComparatorL1 : public Comparator {
      public:
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
        ComparatorL1(size_t d, SharedMemoryAllocator &a) : Comparator(d, a) {}
	double operator()(Object &objecta, Object &objectb) {
	  return ObjectSpaceT::compareL1((OBJECT_TYPE*)&objecta[0], (OBJECT_TYPE*)&objectb[0], dimension);
	}
	double operator()(Object &objecta, PersistentObject &objectb) {
	  return ObjectSpaceT::compareL1((OBJECT_TYPE*)&objecta[0], (OBJECT_TYPE*)&objectb.at(0, allocator), dimension);
	}
	double operator()(PersistentObject &objecta, PersistentObject &objectb) {
	  return ObjectSpaceT::compareL1((OBJECT_TYPE*)&objecta.at(0, allocator), (OBJECT_TYPE*)&objectb.at(0, allocator), dimension);
	}
#else
        ComparatorL1(size_t d) : Comparator(d) {}
	double operator()(Object &objecta, Object &objectb) {
	  return ObjectSpaceT::compareL1((OBJECT_TYPE*)&objecta[0], (OBJECT_TYPE*)&objectb[0], dimension);
	}
#endif
    };

    class ComparatorL2 : public Comparator {
      public:
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
        ComparatorL2(size_t d, SharedMemoryAllocator &a) : Comparator(d, a) {}
	double operator()(Object &objecta, Object &objectb) {
	  return ObjectSpaceT::compareL2((OBJECT_TYPE*)&objecta[0], (OBJECT_TYPE*)&objectb[0], dimension);
	}
	double operator()(Object &objecta, PersistentObject &objectb) {
	  return ObjectSpaceT::compareL2((OBJECT_TYPE*)&objecta[0], (OBJECT_TYPE*)&objectb.at(0, allocator), dimension);
	}
	double operator()(PersistentObject &objecta, PersistentObject &objectb) {
	  return ObjectSpaceT::compareL2((OBJECT_TYPE*)&objecta.at(0, allocator), (OBJECT_TYPE*)&objectb.at(0, allocator), dimension);
	}
#else
        ComparatorL2(size_t d) : Comparator(d) {}
	double operator()(Object &objecta, Object &objectb) {
	  return ObjectSpaceT::compareL2((OBJECT_TYPE*)&objecta[0], (OBJECT_TYPE*)&objectb[0], dimension);
	}
#endif
    };

    class ComparatorHammingDistance : public Comparator {
      public:
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
        ComparatorHammingDistance(size_t d, SharedMemoryAllocator &a) : Comparator(d, a) {}
	double operator()(Object &objecta, Object &objectb) {
	  return ObjectSpaceT::compareHammingDistance((OBJECT_TYPE*)&objecta[0], (OBJECT_TYPE*)&objectb[0], dimension);
	}
	double operator()(Object &objecta, PersistentObject &objectb) {
	  return ObjectSpaceT::compareHammingDistance((OBJECT_TYPE*)&objecta[0], (OBJECT_TYPE*)&objectb.at(0, allocator), dimension);
	}
	double operator()(PersistentObject &objecta, PersistentObject &objectb) {
	  return ObjectSpaceT::compareHammingDistance((OBJECT_TYPE*)&objecta.at(0, allocator), (OBJECT_TYPE*)&objectb.at(0, allocator), dimension);
	}
#else
        ComparatorHammingDistance(size_t d) : Comparator(d) {}
	double operator()(Object &objecta, Object &objectb) {
	  return ObjectSpaceT::compareHammingDistance((OBJECT_TYPE*)&objecta[0], (OBJECT_TYPE*)&objectb[0], dimension);
	}
#endif
    };

    class ComparatorAngleDistance : public Comparator {
      public:
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
        ComparatorAngleDistance(size_t d, SharedMemoryAllocator &a) : Comparator(d, a) {}
	double operator()(Object &objecta, Object &objectb) {
	  return ObjectSpaceT::compareAngleDistance((OBJECT_TYPE*)&objecta[0], (OBJECT_TYPE*)&objectb[0], dimension);
	}
	double operator()(Object &objecta, PersistentObject &objectb) {
	  return ObjectSpaceT::compareAngleDistance((OBJECT_TYPE*)&objecta[0], (OBJECT_TYPE*)&objectb.at(0, allocator), dimension);
	}
	double operator()(PersistentObject &objecta, PersistentObject &objectb) {
	  return ObjectSpaceT::compareAngleDistance((OBJECT_TYPE*)&objecta.at(0, allocator), (OBJECT_TYPE*)&objectb.at(0, allocator), dimension);
	}
#else
        ComparatorAngleDistance(size_t d) : Comparator(d) {}
	double operator()(Object &objecta, Object &objectb) {
	  return ObjectSpaceT::compareAngleDistance((OBJECT_TYPE*)&objecta[0], (OBJECT_TYPE*)&objectb[0], dimension);
	}
#endif
    };

    class ComparatorCosineSimilarity : public Comparator {
      public:
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
        ComparatorCosineSimilarity(size_t d, SharedMemoryAllocator &a) : Comparator(d, a) {}
	double operator()(Object &objecta, Object &objectb) {
	  return ObjectSpaceT::compareCosineSimilarity((OBJECT_TYPE*)&objecta[0], (OBJECT_TYPE*)&objectb[0], dimension);
	}
	double operator()(Object &objecta, PersistentObject &objectb) {
	  return ObjectSpaceT::compareCosineSimilarity((OBJECT_TYPE*)&objecta[0], (OBJECT_TYPE*)&objectb.at(0, allocator), dimension);
	}
	double operator()(PersistentObject &objecta, PersistentObject &objectb) {
	  return ObjectSpaceT::compareCosineSimilarity((OBJECT_TYPE*)&objecta.at(0, allocator), (OBJECT_TYPE*)&objectb.at(0, allocator), dimension);
	}
#else
        ComparatorCosineSimilarity(size_t d) : Comparator(d) {}
	double operator()(Object &objecta, Object &objectb) {
	  return ObjectSpaceT::compareAngleDistance((OBJECT_TYPE*)&objecta[0], (OBJECT_TYPE*)&objectb[0], dimension);
	}
#endif
    };

    ObjectSpaceT(size_t d, const type_info &ot, DistanceType t) : ObjectSpace(d), ObjectRepository(d, ot) {
     size_t objectSize = 0;
     if (ot == typeid(uint8_t)) {
       objectSize = sizeof(uint8_t);
     } else if (ot == typeid(float)) {
       objectSize = sizeof(float);
     } else {
       stringstream msg;
       msg << "ObjectSpace::constructor: Not supported type. " << ot.name();
       NGTThrowException(msg);
     }
     setLength(objectSize * d);
     setDistanceType(t);
   }

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    void open(const string &f, size_t sharedMemorySize) { ObjectRepository::open(f, sharedMemorySize); }
    void copy(PersistentObject &objecta, PersistentObject &objectb) { objecta = objectb; }

    void show(ostream &os, PersistentObject &object) {
      const std::type_info &t = getObjectType();
      if (t == typeid(uint8_t)) {
	unsigned char *optr = static_cast<unsigned char*>(&object.at(0,allocator));
	for (size_t i = 0; i < getDimension(); i++) {
	  os << (int)optr[i] << " ";
	}
      } else if (t == typeid(float)) {
	float *optr = reinterpret_cast<float*>(&object.at(0,allocator));
	for (size_t i = 0; i < getDimension(); i++) {
	  os << optr[i] << " ";
	}
      } else {
	os << " not implement for the type.";
      }
    }

    Object *allocateObject(Object &o) { 
      Object *po = new Object(getByteSizeOfObject());
      for (size_t i = 0; i < getByteSizeOfObject(); i++) {
	(*po)[i] = o[i];
      }
      return po;
    }
    Object *allocateObject(PersistentObject &o) {
      PersistentObject &spo = (PersistentObject &)o;
      Object *po = new Object(getByteSizeOfObject());
      for (size_t i = 0; i < getByteSizeOfObject(); i++) {
	(*po)[i] = spo.at(i,ObjectRepository::allocator);
      }
      return (Object*)po;
    }
    void deleteObject(PersistentObject *po) {
      delete po;
    }
#endif // NGT_SHARED_MEMORY_ALLOCATOR

    void copy(Object &objecta, Object &objectb) {
      objecta.copy(objectb, getByteSizeOfObject());
    }

    void setDistanceType(DistanceType t) {
      if (comparator != 0) {
	delete comparator;
      }
      assert(ObjectSpace::dimension != 0);
      distanceType = t; 
      switch (distanceType) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      case DistanceTypeL1:
	comparator = new ObjectSpaceT::ComparatorL1(ObjectSpace::dimension, ObjectRepository::allocator);
	break;
      case DistanceTypeL2:
	comparator = new ObjectSpaceT::ComparatorL2(ObjectSpace::dimension, ObjectRepository::allocator);
	break;
      case DistanceTypeHamming:
	comparator = new ObjectSpaceT::ComparatorHammingDistance(ObjectSpace::dimension, ObjectRepository::allocator);
	break;
      case DistanceTypeAngle:
	comparator = new ObjectSpaceT::ComparatorAngleDistance(ObjectSpace::dimension, ObjectRepository::allocator);
	break;
      case DistanceTypeCosine:
	comparator = new ObjectSpaceT::ComparatorCosineSimilarity(ObjectSpace::dimension, ObjectRepository::allocator);
	break;
#else
      case DistanceTypeL1:
	comparator = new ObjectSpaceT::ComparatorL1(ObjectSpace::dimension);
	break;
      case DistanceTypeL2:
	comparator = new ObjectSpaceT::ComparatorL2(ObjectSpace::dimension);
	break;
      case DistanceTypeHamming:
	comparator = new ObjectSpaceT::ComparatorHammingDistance(ObjectSpace::dimension);
	break;
      case DistanceTypeAngle:
	comparator = new ObjectSpaceT::ComparatorAngleDistance(ObjectSpace::dimension);
	break;
      case DistanceTypeCosine:
	comparator = new ObjectSpaceT::ComparatorCosineSimilarity(ObjectSpace::dimension);
	break;
#endif
      default:
	cerr << "Distance type is not specified" << endl;
	assert(distanceType != DistanceTypeNone);
	abort();
      }
    }

    static COMPARE_TYPE absolute(int v) { return abs(v); }
    static COMPARE_TYPE absolute(double v) { return fabs(v); }
#if defined(NGT_AVX_DISABLED) || !defined(__AVX__)
    inline static double compareL2(OBJECT_TYPE *a, OBJECT_TYPE *b, size_t size) {
      assert(a != 0);
      assert(b != 0);
      OBJECT_TYPE *last = a + size;
      OBJECT_TYPE *lastgroup = last - 3;
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
#else
    inline static double compareL2(float *a, float *b, size_t size) {
      __m256 sum = _mm256_setzero_ps();
      float *last = a + size;
      float *lastgroup = last - 7;
      while (a < lastgroup) {
	__m256 v = _mm256_sub_ps(_mm256_loadu_ps(a), _mm256_loadu_ps(b));
	sum = _mm256_add_ps(sum, _mm256_mul_ps(v, v));
	a += 8;
	b += 8;
      }
      __attribute__((aligned(32))) float f[8];
      _mm256_store_ps(f, sum);
      double s = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7];
      while (a < last) {
	double d = *a++ - *b++;
	s += d * d;
      }
      return sqrt(s);
    }
    inline static double compareL2(unsigned char *a, unsigned char *b, size_t size) {
      __m128 sum = _mm_setzero_ps();
      unsigned char *last = a + size;
      OBJECT_TYPE *lastgroup = last - 7;
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
    static double compareL1(OBJECT_TYPE *a, OBJECT_TYPE *b, size_t size) {
      assert(a != 0);
      assert(b != 0);
      OBJECT_TYPE *last = a + size;
      OBJECT_TYPE *lastgroup = last - 3;
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
#else
    inline static double compareL1(float *a, float *b, size_t size) {
      __m256 sum = _mm256_setzero_ps();
      float *last = a + size;
      float *lastgroup = last - 7;
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
    inline static double compareL1(unsigned char *a, unsigned char *b, size_t size) {
      __m128 sum = _mm_setzero_ps();
      unsigned char *last = a + size;
      OBJECT_TYPE *lastgroup = last - 7;
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

    inline static double compareHammingDistance(OBJECT_TYPE *a, OBJECT_TYPE *b, size_t size) {
      assert(a != 0);
      assert(b != 0);
      size_t byteSize = sizeof(OBJECT_TYPE) * size;
      assert((0x3 & byteSize) == 0);
      size_t n = byteSize >> 2;
      uint32_t *uinta = (uint32_t*)a;
      uint32_t *uintb = (uint32_t*)b;
      size_t count = 0;
      for (size_t i = 0; i < n; i++) {
	count += popCount(*uinta++ ^ *uintb++);
      }
      return (double)count;
    }

    inline static double compareAngleDistance(OBJECT_TYPE *a, OBJECT_TYPE *b, size_t size) {
      // Calculate the norm of A and B (the supplied vector).
      double normA = 0.0F;
      double normB = 0.0F;
      double sum = 0.0F;
      for (size_t loc = 0; loc < size; loc++) {
	normA += (double)a[loc] * (double)a[loc];
	normB += (double)b[loc] * (double)b[loc];
	sum += (double)a[loc] * (double)b[loc];
      }

      assert(normA > 0.0F);
      assert(normB > 0.0F);

      // Compute the dot product of the two vectors. 
      double cosine = sum / (sqrt(normA) * sqrt(normB));
      // Compute the vector angle from the cosine value, and return.
      // Roundoff error could have put the cosine value out of range.
      // Handle these cases explicitly.
      if (cosine >= 1.0F) {
	return 0.0F;
      } else if (cosine <= -1.0F) {
	return acos (-1.0F);
      } else {
	return acos (cosine);
      }

    }

    inline static double compareCosineSimilarity(OBJECT_TYPE *a, OBJECT_TYPE *b, size_t size) {
      // Calculate the norm of A and B (the supplied vector).
      double normA = 0.0F;
      double normB = 0.0F;
      double sum = 0.0F;
      for (size_t loc = 0; loc < size; loc++) {
	normA += (double)a[loc] * (double)a[loc];
	normB += (double)b[loc] * (double)b[loc];
	sum += (double)a[loc] * (double)b[loc];
      }

      assert(normA > 0.0F);
      assert(normB > 0.0F);

      // Compute the dot product of the two vectors. 
      double cosine = sum / (sqrt(normA) * sqrt(normB));

      return 1.0 - cosine;
    }

    void serialize(const string &ofile) { ObjectRepository::serialize(ofile, this); }
    void deserialize(const string &ifile) { ObjectRepository::deserialize(ifile, this); }
    void serializeAsText(const string &ofile) { ObjectRepository::serializeAsText(ofile, this); }
    void deserializeAsText(const string &ifile) { ObjectRepository::deserializeAsText(ifile, this); }
    void readText(istream &is, size_t dataSize) { ObjectRepository::readText(is, dataSize); }
    void appendText(ifstream &is, size_t dataSize) { ObjectRepository::appendText(is, dataSize); }
    void append(const float *data, size_t dataSize) { ObjectRepository::append(data, dataSize); }
    


#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    PersistentObject *allocatePersistentObject(Object &obj) {
      return ObjectRepository::allocatePersistentObject(obj);
    }
    size_t insert(PersistentObject *obj) { return ObjectRepository::insert(obj); }
#else
    size_t insert(Object *obj) { return ObjectRepository::insert(obj); }
#endif

    void remove(size_t id) { ObjectRepository::remove(id); }

    void linearSearch(Object &query, double radius, size_t size,  ObjectSpace::ResultSet &results) {
      if (!results.empty()) {
	NGTThrowException("lenearSearch: results is not empty");
      }
      ObjectRepository &rep = *this;
      for (size_t idx = 0; idx < rep.size(); idx++) {
	if (rep[idx] == 0) {
	  continue;
	}
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	Distance d = (*comparator)((Object&)query, (PersistentObject&)*rep[idx]);
#else
	Distance d = (*comparator)((Object&)query, (Object&)*rep[idx]);
#endif
	if (radius < 0.0 || d <= radius) {
	  NGT::ObjectDistance obj(idx, d);
	  results.push(obj);
	  if (results.size() > size) {
	    results.pop();
	  }
	}
      }
      return;
    }

    void *getObject(size_t idx) {
      if (idx >= ObjectRepository::size() || idx == 0) {
	stringstream msg;
	msg << "NGT::ObjectSpaceT: The specified ID is out of the range. Object ID is greater than zero. " << idx << ":" << ObjectRepository::size() << ".";
	NGTThrowException(msg);
      }
      PersistentObject &obj = *(*this)[idx];
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      return reinterpret_cast<OBJECT_TYPE*>(&obj.at(0, allocator));
#else
      return reinterpret_cast<OBJECT_TYPE*>(&obj[0]);
#endif
    }

    Object *allocateObject() { return ObjectRepository::allocateObject(); }
    void deleteObject(Object *po) { ObjectRepository::deleteObject(po); }

    Object *allocateObject(const string &textLine, const string &sep) {
      return ObjectRepository::allocateObject(textLine, sep);
    }
    Object *allocateObject(vector<double> &obj) {
      return ObjectRepository::allocateObject(obj);
    }

    size_t getSize() { return ObjectRepository::size(); }
    size_t getSizeOfElement() { return sizeof(OBJECT_TYPE); }
    const std::type_info &getObjectType() { return typeid(OBJECT_TYPE); };
    size_t getByteSizeOfObject() { return getByteSize(); }

    ObjectRepository &getRepository() { return *this; };

    void show(ostream &os, Object &object) {
      const std::type_info &t = getObjectType();
      if (t == typeid(uint8_t)) {
	unsigned char *optr = static_cast<unsigned char*>(&object[0]);
	for (size_t i = 0; i < getDimension(); i++) {
	  os << (int)optr[i] << " ";
	}
      } else if (t == typeid(float)) {
	float *optr = reinterpret_cast<float*>(&object[0]);
	for (size_t i = 0; i < getDimension(); i++) {
	  os << optr[i] << " ";
	}
      } else {
	os << " not implement for the type.";
      }
    }
  };

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  // set v in objectspace to this object using allocator.
  inline void PersistentObject::set(PersistentObject &po, ObjectSpace &objectspace) {
    SharedMemoryAllocator &allocator = objectspace.getRepository().getAllocator();
    uint8_t *src = (uint8_t *)&po.at(0, allocator);
    uint8_t *dst = (uint8_t *)&(*this).at(0, allocator);
    memcpy(dst, src, objectspace.getByteSizeOfObject());
  }

  inline off_t PersistentObject::allocate(ObjectSpace &objectspace) {
    SharedMemoryAllocator &allocator = objectspace.getRepository().getAllocator();
    return allocator.getOffset(new(allocator) PersistentObject(allocator, &objectspace));
  }

  inline void PersistentObject::serializeAsText(ostream &os, ObjectSpace *objectspace) { 
    assert(objectspace != 0);
    SharedMemoryAllocator &allocator = objectspace->getRepository().getAllocator();
    const std::type_info &t = objectspace->getObjectType();
    void *ref = &(*this).at(0, allocator);
    size_t dimension = objectspace->getDimension();
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
      cerr << "ObjectT::serializeAsText: not supported data type. [" << t.name() << "]" << endl;
      assert(0);
    }
  }

  inline void PersistentObject::deserializeAsText(ifstream &is, ObjectSpace *objectspace) {
    assert(objectspace != 0);
    SharedMemoryAllocator &allocator = objectspace->getRepository().getAllocator();
    const std::type_info &t = objectspace->getObjectType();
    size_t dimension = objectspace->getDimension();
    void *ref = &(*this).at(0, allocator);
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

#endif
} // namespace NGT

