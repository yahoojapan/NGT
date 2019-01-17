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

namespace NGT {
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

    virtual PersistentObject *allocateNormalizedPersistentObject(const vector<double> &obj) {
      cerr << "ObjectRepository::allocateNormalizedPersistentObject: Fatal error! Something wrong!" << endl;
      abort();
    }

    virtual PersistentObject *allocateNormalizedPersistentObject(const vector<float> &obj) {
      cerr << "ObjectRepository::allocateNormalizedPersistentObject: Fatal error! Something wrong!" << endl;
      abort();
    }

    virtual PersistentObject *allocateNormalizedPersistentObject(const float *obj, size_t size) {
      cerr << "ObjectRepository::allocateNormalizedPersistentObject: Fatal error! Something wrong!" << endl;
      abort();
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
	  PersistentObject *obj = 0;
	  try {
	    obj = allocateNormalizedPersistentObject(object);
	  } catch (Exception &err) {
	    cerr << err.what() << " continue..." << endl;
	    obj = allocatePersistentObject(object);
	  }
	  push_back(obj);
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
	vector<double> object;
	object.reserve(dimension);
	for (size_t dataidx = 0; dataidx < dimension; dataidx++) {
	  object.push_back(data[dataidx]);
	}
	try {
	  PersistentObject *obj = 0;
	  try {
	    obj = allocateNormalizedPersistentObject(object);
	  } catch (Exception &err) {
	    cerr << err.what() << " continue..." << endl;
	    obj = allocatePersistentObject(object);
	  }
	  push_back(obj);

	} catch (Exception &err) {
	  std::cerr << "ObjectSpace::readText: Warning! Invalid data. Skip the data no. " << idx << " and continue." << std::endl;
	}
      }
    }

    Object *allocateObject() {
      return (Object*) new Object(byteSize);
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
      Object *allocateObject(const vector<T> &o) {
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
    PersistentObject *allocatePersistentObject(const vector<T> &o) {
      return allocatePersistentObject(o.data(), o.size());
    }

#else
    // ObjectRepository
    template <typename T>
    PersistentObject *allocatePersistentObject(const vector<T> &o) {
      return allocateObject(o);
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
    const size_t dimension;
    const type_info &type;
   protected:
    size_t byteSize;		// the length of all of elements.
  };

} // namespace NGT
