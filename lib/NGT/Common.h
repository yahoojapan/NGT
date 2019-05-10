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

#include	<string>
#include	<vector>
#include	<queue>
#include	<map>
#include	<iostream>
#include	<sstream>
#include	<fstream>
#include	<cassert>
#include	<cstdlib>
#include	<cmath>
#include	<cfloat>
#include	<climits>
#include	<iomanip>
#include	<algorithm>
#include	<typeinfo>

#include	<sys/time.h>

#include	"NGT/defines.h"
#include	"NGT/SharedMemoryAllocator.h"

#define ADVANCED_USE_REMOVED_LIST
#define	SHARED_REMOVED_LIST

using namespace std;

namespace NGT {
  typedef	unsigned int	ObjectID;
  typedef	float		Distance;

#define	NGTThrowException(MESSAGE)			throw NGT::Exception(__FILE__, (size_t)__LINE__, MESSAGE)
#define	NGTThrowSpecificException(MESSAGE, TYPE)	throw NGT::TYPE(__FILE__, (size_t)__LINE__, MESSAGE)

  class Exception : public std::exception {
  public:
    Exception():message("No message") {}
    Exception(const string &file, size_t line, stringstream &m) { set(file, line, m.str()); }
    Exception(const string &file, size_t line, const string &m) { set(file, line, m); }
    void set(const string &file, size_t line, const string &m) { 
      stringstream ss;
      ss << file << ":" << line << ": " << m;
      message = ss.str(); 
    }
    ~Exception() throw() {}
    Exception &operator=(const Exception &e) {
      message = e.message;
      return *this;
    }
    virtual const char *what() const throw() {
      return message.c_str();
    }
    string &getMessage() { return message; }
  protected:
    string message;
  };

  class Args : public map<string, string>
  {
   public:
    Args() {}
    Args(int argc, char **argv):
    option("a:b:c:d:e:f:g:hi:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:"
	   "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z:")
    {
      int opt;
      while ((opt = getopt(argc, argv, option)) != -1) {
	if ((char)opt == 'h') {
	  string str;
	  str.append(1, (char)opt);
	  insert(pair<string, string>(str, ""));
	  continue;
	}
	string str;
	str.append(1, (char)opt);
	insert(pair<string, string>(str, string(optarg)));
      }
      for (int i = 0; optind < argc; optind++, i++) {
	stringstream ss;
	ss << "#" << i;
	insert(pair<string, string>(ss.str(), string(argv[optind])));
      }
    }
    string &find(const char *s) { return get(s); }
    char getChar(const char *s, char v) {
      try {
	return get(s)[0];
      } catch (...) {
	return v;
      }
    }
    string getString(const char *s, const char *v) {
      try {
	return get(s);
      } catch (...) {
	return v;
      }
    }
    string &get(const char *s) {
      Args::iterator ai;
      ai = map<string, string>::find(string(s));
      if (ai == this->end()) {
	stringstream msg;
	msg << s << ": Not specified" << endl;
	NGTThrowException(msg.str());
      }
      return ai->second;
    }
    long getl(const char *s, long v) {
      char *e;
      long val;
      try {
	val = strtol(get(s).c_str(), &e, 10);
      } catch (...) {
	return v;
      }
      if (*e != 0) {
	stringstream msg;
	msg << "ARGS::getl: Illegal string. Option=-" << s << " Specified value=" << get(s) 
	    << " Illegal string=" << e << endl;
	NGTThrowException(msg.str());
      }
      return val;
    }
    float getf(const char *s, float v) {
      char *e;
      float val;
      try {
	val = strtof(get(s).c_str(), &e);
      } catch (...) {
	return v;
      }
      if (*e != 0) {
	stringstream msg;
	msg << "ARGS::getf: Illegal string. Option=-" << s << " Specified value=" << get(s) 
	    << " Illegal string=" << e << endl;
	NGTThrowException(msg.str());
      }
      return val;
    }
    const char *option;
  };

  class Common {
  public:
    static void tokenize(const string &str, vector<std::string> &token, const string seps) {
      string::size_type current = 0;
      string::size_type next;
      while ((next = str.find_first_of(seps, current)) != string::npos) {
	token.push_back(str.substr(current, next - current));
	current = next + 1;
      }
      string t = str.substr(current);
      token.push_back(t);
    }

    static double strtod(const string &str) {
      char *e;
      double val = std::strtod(str.c_str(), &e);
      if (*e != 0) {
	stringstream msg;
	msg << "Invalid string. " << e;
	NGTThrowException(msg);
      }
      return val;
    }

    static long strtol(const string &str, int base = 10) {
      char *e;
      long val = std::strtol(str.c_str(), &e, base);
      if (*e != 0) {
	stringstream msg;
	msg << "Invalid string. " << e;
	NGTThrowException(msg);
      }
      return val;
    }


    static string getProcessStatus(const string &stat) {
      pid_t pid = getpid();
      stringstream str;
      str << "/proc/" << pid << "/status";
      ifstream procStatus(str.str());
      if (!procStatus.fail()) {
	string line;
	while (getline(procStatus, line)) {
	  vector<string> tokens;
	  NGT::Common::tokenize(line, tokens, ": \t");
	  if (tokens[0] == stat) {
	    for (size_t i = 1; i < tokens.size(); i++) {
	      if (tokens[i].empty()) {
		continue;
	      }
	      return tokens[i];
	    }
	  }
	}
      }
      return "-1";
    }

    // size unit is kbyte
    static int getProcessVmSize() { return strtol(getProcessStatus("VmSize")); }
    static int getProcessVmPeak() { return strtol(getProcessStatus("VmPeak")); }
    static int getProcessVmRSS() { return strtol(getProcessStatus("VmRSS")); }
  };


  template <class TYPE> 
    class CompactVector {
  public:
    typedef TYPE *	iterator;

    CompactVector() : vector(0), vectorSize(0), allocatedSize(0){}
    virtual ~CompactVector() { clear(); }

    void clear() {
      if (vector != 0) {
	delete[] vector;
      }
      vector = 0;
      vectorSize = 0;
      allocatedSize = 0;
    }

    TYPE &front() { return vector[0]; }
    TYPE &back() { return vector[vectorSize - 1]; }
    bool empty() { return vector == 0; }
    iterator begin() { return &(vector[0]); }
    iterator end() { return begin() + vectorSize; }
    TYPE &operator[](size_t idx) const { return vector[idx]; }

    CompactVector &operator=(CompactVector<TYPE> &v) {
      assert((vectorSize == v.vectorSize) || (vectorSize == 0));
      if (vectorSize == v.vectorSize) {
	for (size_t i = 0; i < vectorSize; i++) {
	  vector[i] = v[i];
	}
	return *this;
      } else {
	reserve(v.vectorSize);
	assert(allocatedSize >= v.vectorSize);
	for (size_t i = 0; i < v.vectorSize; i++) {
	  push_back(v.at(i));
	}
	vectorSize = v.vectorSize;
	assert(vectorSize == v.vectorSize);
      }
      return *this;
    }

    TYPE &at(size_t idx) const {
      if (idx >= vectorSize) {
	stringstream msg;
	msg << "CompactVector: beyond the range. " << idx << ":" << vectorSize;
	NGTThrowException(msg);  
      }
      return vector[idx];
    }

    iterator erase(iterator b, iterator e) {
      iterator ret;
      e = end() < e ? end() : e;
      for (iterator i = b; i < e; i++) {
	ret = erase(i);
      }
      return ret;
    }

    iterator erase(iterator i) {
      iterator back = i;
      vectorSize--;
      iterator e = end();
      for (; i < e; i++) {
	*i = *(i + 1);
      }
      return ++back;
    }

    void pop_back() {
      if (vectorSize > 0) {
	vectorSize--;
      }
    }

    iterator insert(iterator &i, const TYPE &data) {
      if (size() == 0) {
	push_back(data);
	return end();
      }
      off_t oft = i - begin();
      extend();
      i = begin() + oft;
      iterator b = begin();
      for (iterator ci = end(); ci > i && ci != b; ci--) {
	*ci = *(ci - 1);
      }
      *i = data;
      vectorSize++;
      return i + 1;
    }

    void push_back(const TYPE &data) {
      extend();
      vector[vectorSize] = data;
      vectorSize++;
    }

    void reserve(size_t s) {
      if (s <= allocatedSize) {
	return;
      } else {
	TYPE *newptr = new TYPE[s];
	TYPE *dstptr = newptr;
	TYPE *srcptr = vector;
	TYPE *endptr = srcptr + vectorSize;
	while (srcptr < endptr) {
	  *dstptr++ = *srcptr;
	  (*srcptr).~TYPE();
	  srcptr++;
	}
	allocatedSize = s;
	if (vector != 0) {
	  delete[] vector;
	}
	vector = newptr;
      }
    }

    void resize(size_t s, TYPE v = TYPE()) {
      if (s > allocatedSize) {
	size_t asize = allocatedSize == 0 ? 1 : allocatedSize;
	while (asize < s) {
	  asize <<= 1;
	}
	reserve(asize);
	TYPE *base = vector;
	TYPE *dstptr = base + vectorSize;
	TYPE *endptr = base + s;
	for (; dstptr < endptr; dstptr++) {
	  *dstptr = v;
	}
      }
      vectorSize = s;
    }

    size_t size() const { return (size_t)vectorSize; }

    void extend() {
      if (vectorSize == allocatedSize) {
	if (vectorSize == 0) {
	  reserve(2);
	} else {
	  uint64_t size = vectorSize;
	  size <<= 1;
	  if (size > 0xffff) {
	    cerr << "CompactVector is too big. " << size << endl;
	    abort();
	  }
	  reserve(size);
	}
      }
    }

    TYPE *vector;
    uint16_t vectorSize;
    uint16_t allocatedSize;
  };

  class CompactString {
  public:
    CompactString():vector(0) { }

    CompactString(const CompactString &v):vector(0) { *this = v; }

    ~CompactString() { clear(); }

    void clear() {
      if (vector != 0) {
	delete[] vector;
      }
      vector = 0;
    }

    CompactString &operator=(const string &v) { return *this = v.c_str(); }

    CompactString &operator=(const CompactString &v) { return *this = v.vector; }

    CompactString &operator=(const char *str) {
      if (str == 0 || strlen(str) == 0) {
	clear();
	return *this;
      }
      if (size() != strlen(str)) {
	clear();
	vector = new char[strlen(str) + 1];
      }
      strcpy(vector, str);
      return *this;
    }

    char &at(size_t idx) const {
      if (idx >= size()) {
	NGTThrowException("CompactString: beyond the range");  
      }
      return vector[idx];
    }

    char *c_str() { return vector; }
    size_t size() const { 
      if (vector == 0) {
	return 0;
      } else {
	return (size_t)strlen(vector); 
      }
    }

    char *vector;
  };

  // BooleanSet has been already optimized.
  class BooleanSet {
  public:
    BooleanSet(size_t s) {
      size = (s >> 6) + 1; // 2^6=64
      bitvec.resize(size);
    }
    inline uint64_t getBitString(size_t i) { return (uint64_t)1 << (i & (64 - 1)); }
    inline uint64_t &getEntry(size_t i) { return bitvec[i >> 6]; }
    inline bool operator[](size_t i) {
      return (getEntry(i) & getBitString(i)) != 0;
    }
    inline void set(size_t i) {
      getEntry(i) |= getBitString(i);
    }
    inline void insert(size_t i) { set(i); }
    inline void reset(size_t i) {
      getEntry(i) &= getBitString(i);
    }
    vector<uint64_t>	bitvec;
    uint64_t		size;
  };


  class PropertySet : public map<string, string> {
  public:
    void set(const string &key, const string &value) {
      iterator it = find(key);
      if (it == end()) {
	insert(pair<string, string>(key, value));
      } else {
	(*it).second = value;
      }
    }
    template <class VALUE_TYPE> void set(const string &key, VALUE_TYPE value) {
      stringstream vstr;
      vstr << value;
      iterator it = find(key);
      if (it == end()) {
	insert(pair<string, string>(key, vstr.str()));
      } else {
	(*it).second = vstr.str();
      }
    }

    string get(const string &key) {
      iterator it = find(key);
      if (it != end()) {
	return it->second;
      }
      return "";
    }
    float getf(const string &key, float defvalue) {
      iterator it = find(key);
      if (it != end()) {
	char *e = 0;
	float val = strtof(it->second.c_str(), &e);
	if (*e != 0) {
	  cerr << "Warning: Illegal property. " << key << ":" << it->second << " (" << e << ")" << endl;
	  return defvalue;
	}
	return val;
      }
      return defvalue;
    }
    void updateAndInsert(PropertySet &prop) {
      for (map<string, string>::iterator i = prop.begin(); i != prop.end(); ++i) {
	set((*i).first, (*i).second);
      }
    }
    long getl(const string &key, long defvalue) {
      iterator it = find(key);
      if (it != end()) {
	char *e = 0;
	float val = strtol(it->second.c_str(), &e, 10);
	if (*e != 0) {
	  cerr << "Warning: Illegal property. " << key << ":" << it->second << " (" << e << ")" << endl;
	}
	return val;
      }
      return defvalue;
    }
    void load(const string &f) { 
      ifstream st(f); 
      if (!st) {
	stringstream msg;
	msg << "PropertySet::load: Cannot load the property file " << f << ".";
	NGTThrowException(msg);
      }
      load(st);
    }
    void save(const string &f) {    
      ofstream st(f); 
      if (!st) {
	stringstream msg;
	msg << "PropertySet::save: Cannot save. " << f << endl;
	NGTThrowException(msg);
      }
      save(st); 
    }    
    void save(ofstream &os) {
      for (map<string, string>::iterator i = this->begin(); i != this->end(); i++) {
	os << i->first << "\t" << i->second << std::endl;
      }
    }
    void load(ifstream &is) {
      string line;
      while (getline(is, line)) {
	vector<string> tokens;
	NGT::Common::tokenize(line, tokens, "\t");
	if (tokens.size() != 2) {
	  cerr << "Property file is illegal. " << line << endl;
	  continue;
	}
	set(tokens[0], tokens[1]);
      }
    }
  };

  namespace Serializer {
    static inline void read(istream &is, uint8_t *v, size_t s) {
      is.read((char*)v, s);
    }

    static inline void write(ostream &os, const uint8_t *v, size_t s) {
      os.write((const char*)v, s);
    }

    template <typename TYPE> void write(ostream &os, const TYPE v) {
      os.write((const char*)&v, sizeof(TYPE));
    }

    template <typename TYPE> void writeAsText(ostream &os, const TYPE v) {
      if (typeid(TYPE) == typeid(unsigned char)) {
	os << (int)v;
      } else {
	os << v;
      }
    }

    template <typename TYPE> void read(istream &is, TYPE &v) {
      is.read((char*)&v, sizeof(TYPE));
    }

    template <typename TYPE> void readAsText(istream &is, TYPE &v) {
      if (typeid(TYPE) == typeid(unsigned char)) {
	unsigned int tmp;
	is >> tmp;
	if (tmp > 255) {
	  cerr << "Error! Invalid. " << tmp << endl;
	}
	v = (TYPE)tmp;
      } else {
	is >> v;
      }
    }

    template <typename TYPE> void write(ostream &os, const vector<TYPE> &v) { 
      unsigned int s = v.size();
      write(os, s);
      for (unsigned int i = 0; i < s; i++) {
	write(os, v[i]);
      }
    }

    template <typename TYPE> void writeAsText(ostream &os, const vector<TYPE> &v) { 
      unsigned int s = v.size();
      os << s << " ";
      for (unsigned int i = 0; i < s; i++) {
	writeAsText(os, v[i]);
	os << " ";
      }
    }

    template <typename TYPE> void write(ostream &os, const CompactVector<TYPE> &v) { 
      unsigned int s = v.size();
      write(os, s);
      for (unsigned int i = 0; i < s; i++) {
	write(os, v[i]);
      }
    }

    template <typename TYPE> void writeAsText(ostream &os, const CompactVector<TYPE> &v) { 
      unsigned int s = v.size();
      for (unsigned int i = 0; i < s; i++) {
	writeAsText(os, v[i]);
	os << " ";
      }
    }

    template <typename TYPE> void writeAsText(ostream &os, TYPE *v, size_t s) { 
      os << s << " ";
      for (unsigned int i = 0; i < s; i++) {
	writeAsText(os, v[i]);
	os << " ";
      }
    }

    template <typename TYPE> void read(istream &is, vector<TYPE> &v) { 
      v.clear();
      unsigned int s;
      read(is, s);
      v.reserve(s);
      for (unsigned int i = 0; i < s; i++) {
	TYPE val;
	read(is, val);
	v.push_back(val);
      }
    }

    template <typename TYPE> void readAsText(istream &is, vector<TYPE> &v) { 
      v.clear();
      unsigned int s;
      is >> s;
      for (unsigned int i = 0; i < s; i++) {
	TYPE val;
	Serializer::readAsText(is, val);
	v.push_back(val);
      }
    }


    template <typename TYPE> void read(istream &is, CompactVector<TYPE> &v) { 
      v.clear();
      unsigned int s;
      read(is, s);
      v.reserve(s);
      for (unsigned int i = 0; i < s; i++) {
	TYPE val;
	read(is, val);
	v.push_back(val);
      }
    }

    template <typename TYPE> void readAsText(istream &is, CompactVector<TYPE> &v) { 
      v.clear();
      unsigned int s;
      is >> s;
      for (unsigned int i = 0; i < s; i++) {
	TYPE val;
	Serializer::readAsText(is, val);
	v.push_back(val);
      }
    }

    template <typename TYPE> void readAsText(istream &is,  TYPE *v, size_t s) { 
      unsigned int size;
      is >> size;
      if (s != size) {
	cerr << "readAsText: something wrong. " << size << ":" << s << endl;
	return;
      }
      for (unsigned int i = 0; i < s; i++) {
	TYPE val;
	Serializer::readAsText(is, val);
	v[i] = val;
      }
    }


  } // namespace Serialize


  class ObjectSpace;


#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  template <class TYPE> 
    class Vector {
  public:
    typedef TYPE *	iterator;

    Vector() : vector(0), vectorSize(0), allocatedSize(0) {}

    void clear(SharedMemoryAllocator &allocator) {
      if (vector != 0) {
	allocator.free(allocator.getAddr(vector));
      }
      vector = 0;
      vectorSize = 0;
      allocatedSize = 0;
    }

    TYPE &front(SharedMemoryAllocator &allocator) { return (*this).at(0, allocator); }
    TYPE &back(SharedMemoryAllocator &allocator) { return (*this).at(vectorSize - 1, allocator); }
    bool empty() { return vectorSize == 0; }
    iterator begin(SharedMemoryAllocator &allocator) { 
      return (TYPE*)allocator.getAddr((off_t)vector); 
    }
    iterator end(SharedMemoryAllocator &allocator) {
      return begin(allocator) + vectorSize;
    }

    Vector &operator=(Vector<TYPE> &v) {
      assert((vectorSize == v.vectorSize) || (vectorSize == 0));
      if (vectorSize == v.vectorSize) {
	for (size_t i = 0; i < vectorSize; i++) {
	  (*this)[i] = v[i];
	}
	return *this;
      } else {
	reserve(v.vectorSize);
	assert(allocatedSize >= v.vectorSize);
	for (size_t i = 0; i < v.vectorSize; i++) {
	  push_back(v.at(i));
	}
	vectorSize = v.vectorSize;
	assert(vectorSize == v.vectorSize);
      }
      return *this;
    }

    TYPE &at(size_t idx, SharedMemoryAllocator &allocator) {
      if (idx >= vectorSize) {
	stringstream msg;
	msg << "Vector: beyond the range. " << idx << ":" << vectorSize;
	NGTThrowException(msg);  
      }
      return *(begin(allocator) + idx);
    }

    iterator erase(iterator b, iterator e, SharedMemoryAllocator &allocator) {
      iterator ret;
      e = end(allocator) < e ? end(allocator) : e;
      for (iterator i = b; i < e; i++) {
	ret = erase(i, allocator);
      }
      return ret;
    }

    iterator erase(iterator i, SharedMemoryAllocator &allocator) {
      iterator back = i;
      vectorSize--;
      iterator e = end(allocator);
      for (; i < e; i++) {
	*i = *(i + 1);
      }
      return back;
    }

    void pop_back() {
      if (vectorSize > 0) {
	vectorSize--;
      }
    }
    iterator insert(iterator &i, const TYPE &data, SharedMemoryAllocator &allocator) {
      if (size() == 0) {
	push_back(data, allocator);
	return end(allocator);
      }
      off_t oft = i - begin(allocator);
      extend(allocator);
      i = begin(allocator) + oft;
      iterator b = begin(allocator);
      for (iterator ci = end(allocator); ci > i && ci != b; ci--) {
	*ci = *(ci - 1);
      }
      *i = data;
      vectorSize++;
      return i + 1;
    }

    void push_back(const TYPE &data, SharedMemoryAllocator &allocator) {
      extend(allocator);
      vectorSize++;
      (*this).at(vectorSize - 1, allocator) = data;
    }

    void reserve(size_t s, SharedMemoryAllocator &allocator) {
      if (s <= allocatedSize) {
	return;
      } else {
	TYPE *newptr = new(allocator) TYPE[s];
	TYPE *dstptr = newptr;
	TYPE *srcptr = (TYPE*)allocator.getAddr(vector);
	TYPE *endptr = srcptr + vectorSize;
	while (srcptr < endptr) {
	  *dstptr++ = *srcptr;
	  (*srcptr).~TYPE();
	  srcptr++;
	}
	allocatedSize = s;
	if (vector != 0) {
	  allocator.free(allocator.getAddr(vector));
	}
	vector = allocator.getOffset(newptr);
      }
    }

    void resize(size_t s, SharedMemoryAllocator &allocator, TYPE v = TYPE()) {
      if (s > allocatedSize) {
	size_t asize = allocatedSize == 0 ? 1 : allocatedSize;
	while (asize < s) {
	  asize <<= 1;
	}
	reserve(asize, allocator);
	TYPE *base = (TYPE*)allocator.getAddr(vector);
	TYPE *dstptr = base + vectorSize;
	TYPE *endptr = base + s;
	for (; dstptr < endptr; dstptr++) {
	  *dstptr = v;
	}
      }
      vectorSize = s;
    }

    void serializeAsText(ostream &os, ObjectSpace *objectspace = 0) { 
      unsigned int s = size();
      os << s << " ";
      for (unsigned int i = 0; i < s; i++) {
	Serializer::writeAsText(os, (*this)[i]);
	os << " ";
      }
    }
    void deserializeAsText(istream &is, ObjectSpace *objectspace = 0) { 
      clear();
      size_t s;
      Serializer::readAsText(is, s);
      resize(s);
      for (unsigned int i = 0; i < s; i++) {
	Serializer::readAsText(is, (*this)[i]);
      }
    }
    size_t size() { return vectorSize; }

  public:
    void extend(SharedMemoryAllocator &allocator) {
      extend(vectorSize, allocator);
    }

    void extend(size_t idx, SharedMemoryAllocator &allocator) {
      if (idx >= allocatedSize) {
	if (idx == 0) {
	  reserve(2, allocator);
	} else {
	  uint64_t size = allocatedSize == 0 ? 1 : allocatedSize;
	  do {
	    size <<= 1;
	  } while (size <= idx);
	  if (size > 0xffffffff) {
	    cerr << "Vector is too big. " << size << endl;
	    abort();
	  }
	  reserve(size, allocator);
	}
      }
    }

    off_t vector;
    uint32_t vectorSize;
    uint32_t allocatedSize;
  };

#endif // NGT_SHARED_MEMORY_ALLOCATOR


#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  class ObjectSpace;
  template <class TYPE>
    class PersistentRepository {
  public:
    typedef Vector<off_t> ARRAY;
    typedef TYPE **		iterator;

    PersistentRepository():array(0) {}
    ~PersistentRepository() { close(); }

    void open(const string &mapfile, size_t sharedMemorySize) {
      assert(array == 0);
      SharedMemoryAllocator &allocator = getAllocator();
#ifdef ADVANCED_USE_REMOVED_LIST
      off_t *entryTable = (off_t*)allocator.construct(mapfile, sharedMemorySize);
      if (entryTable == 0) {
	entryTable = (off_t*)construct();
	allocator.setEntry(entryTable);
      }
      assert(entryTable != 0);
      initialize(entryTable);
#else
      void *entry = allocator.construct(mapfile, sharedMemorySize);
      if (entry == 0) {
	array = (ARRAY*)construct();
	allocator.setEntry(array);
      } else {
	array = (ARRAY*)entry;
      }
#endif
    }

    void close() {
      getAllocator().destruct();
    }

#ifdef ADVANCED_USE_REMOVED_LIST
    void *construct() {
      SharedMemoryAllocator &allocator = getAllocator();
      off_t *entryTable = new(allocator) off_t[2];
      entryTable[0] = allocator.getOffset(new(allocator) ARRAY);
      entryTable[1] = allocator.getOffset(new(allocator) Vector<size_t>);
      return entryTable;
    }
    void initialize(void *e) {
      SharedMemoryAllocator &allocator = getAllocator();
      off_t *entryTable = (off_t*)e;
      array = (ARRAY*)allocator.getAddr(entryTable[0]);
      removedList = (Vector<size_t>*)allocator.getAddr(entryTable[1]);
    }
    size_t removedListPop() {
      assert(removedList->size() != 0);
      size_t idx = removedList->back(allocator);
      removedList->pop_back();
      return idx;
    }

    void removedListPush(size_t id) {
      if (removedList->size() == 0) {
	removedList->push_back(id, allocator);
	return;
      }
      Vector<size_t>::iterator rmi
	= std::lower_bound(removedList->begin(allocator), removedList->end(allocator), id, greater<size_t>());
      if ((rmi != removedList->end(allocator)) && ((*rmi) == id)) {
	cerr << "removedListPush: already existed! continue... ID=" << id << endl;
	return;
      }
      removedList->insert(rmi, id, allocator);
    }
#else
    void *construct() {
      SharedMemoryAllocator &allocator = getAllocator();
      return new(allocator) ARRAY;
    }
    void initialize(void *e) {
      assert(array == 0);
      assert(e != 0);
      array = (ARRAY*)e;
    }
#endif
    TYPE *allocate() { return new(allocator) TYPE(allocator); }

    size_t push(TYPE *n) {
      if (size() == 0) {
	push_back(0);
      }
      push_back(n);
      return size() - 1;
    }

    size_t insert(TYPE *n) {
#ifdef ADVANCED_USE_REMOVED_LIST
      if (!removedList->empty()) {
	size_t idx = removedListPop();
	put(idx, n);
	return idx;
      }
#endif
      return push(n);
    }

    bool isEmpty(size_t idx) {
      if (idx < size()) {
	return (*array).at(idx, allocator) == 0;
      } else {
	return true;
      }
    }

    void put(size_t idx, TYPE *n) {
      (*array).extend(idx, allocator);
      if (size() <= idx) {
	resize(idx + 1);
      }
      if ((*this)[idx] != 0) {
	NGTThrowException("put: Not empty");  
      }
      set(idx, n);
    }

    void remove(size_t idx) {
      if (isEmpty(idx)) {
	NGTThrowException("remove: Not in-memory or invalid id");  
      }
      SharedMemoryAllocator &allocator = getAllocator();
      (*this)[idx]->~TYPE();
      allocator.free((*this)[idx]);
      set(idx, (TYPE*)0);
#ifdef ADVANCED_USE_REMOVED_LIST
      removedListPush(idx);
#endif
    }

    inline TYPE *get(size_t idx) {
      if (isEmpty(idx)) {
	stringstream msg;
	msg << "get: Not in-memory or invalid offset of node. " << idx << ":" << array->size();
	NGTThrowException(msg.str());
      }
      return (*this)[idx];
    }

    void serialize(ofstream &os, ObjectSpace *objectspace = 0) {
      NGT::Serializer::write(os, array->size());    
      for (size_t idx = 0; idx < array->size(); idx++) {
	if ((*this)[idx] == 0) {
	  NGT::Serializer::write(os, '-');
	} else {
	  NGT::Serializer::write(os, '+');
	  if (objectspace == 0) {
	    assert(0);
	    //(*this)[idx]->serialize(os, allocator);
	  } else {
	    assert(0);
	    //(*this)[idx]->serialize(os, allocator, objectspace);
	  }
	}
      }
    }

    void deserialize(ifstream &is, ObjectSpace *objectspace = 0) {
      if (!is.is_open()) {
	NGTThrowException("NGT::Common: Not open the specified stream yet.");
      }
      deleteAll();
      (*this).push_back((TYPE*)0);
      size_t s;
      NGT::Serializer::read(is, s);
      for (size_t i = 0; i < s; i++) {
	char type;
	NGT::Serializer::read(is, type);
	switch(type) {
	case '-':
	  {
	    (*this).push_back((TYPE*)0);
#ifdef ADVANCED_USE_REMOVED_LIST
	    if (i != 0) {
	      removedListPush(i);
	    }
#endif
	  }
	  break;
	case '+':
	  {
	    if (objectspace == 0) {
	      TYPE *v = new(allocator) TYPE(allocator);
	      //v->deserialize(is, allocator);
	      assert(0);
	      (*this).push_back(v);
	    } else {
	      TYPE *v = new(allocator) TYPE(allocator, objectspace);
	      //v->deserialize(is, allocator, objectspace);
	      assert(0);
	      (*this).push_back(v);
	    }
	  }
	  break;
	default:
	  {
	    assert(type == '-' || type == '+');
	    break;
	  }
	}
      }
    }

    void serializeAsText(ofstream &os, ObjectSpace *objectspace = 0) {
      os.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
      os << std::setprecision(8);

      os << array->size() << endl;
      for (size_t idx = 0; idx < array->size(); idx++) {
	if ((*this)[idx] == 0) {
	  os << idx << " - " << endl;
	} else {
	  os << idx << " + ";
	  if (objectspace == 0) {
	    (*this)[idx]->serializeAsText(os, allocator);
	  } else {
	    (*this)[idx]->serializeAsText(os, allocator, objectspace);
	  }
	  os << endl;
	}
      }
      os << fixed;
    }


    void deserializeAsText(ifstream &is, ObjectSpace *objectspace = 0) {
      if (!is.is_open()) {
	NGTThrowException("NGT::Common: Not open the specified stream yet.");
      }
      deleteAll();
      size_t s;
      NGT::Serializer::readAsText(is, s);
      (*this).reserve(s);
      for (size_t i = 0; i < s; i++) {
	size_t idx;
	NGT::Serializer::readAsText(is, idx);
	if (i != idx) {
	  cerr << "PersistentRepository: Error. index of a specified import file is invalid. " << idx << ":" << i << endl;
	}
	char type;
	NGT::Serializer::readAsText(is, type);
	switch(type) {
	case '-':
	  {
	    (*this).push_back((TYPE*)0);
#ifdef ADVANCED_USE_REMOVED_LIST
	    if (i != 0) {
	      removedListPush(i);
	    }
#endif
	  }
	  break;
	case '+':
	  {
	    if (objectspace == 0) {
	      TYPE *v = new(allocator) TYPE(allocator);
	      v->deserializeAsText(is, allocator);
	      (*this).push_back(v);
	    } else {
	      TYPE *v = new(allocator) TYPE(allocator, objectspace);
	      v->deserializeAsText(is, allocator, objectspace);
	      (*this).push_back(v);
	    }
	  }
	  break;
	default:
	  {
	    assert(type == '-' || type == '+');
	    break;
	  }
	}
      }
    }
    void deleteAll() {
      for (size_t i = 0; i < array->size(); i++) {
	if ((*array).at(i, allocator) != 0) {
	  allocator.free((*this)[i]);
	}
      }
      array->clear(allocator);
#ifdef ADVANCED_USE_REMOVED_LIST
      while (!removedList->empty()) { removedListPop(); }
#endif
    }

    void set(size_t idx, TYPE *n) {
      array->at(idx, allocator) = allocator.getOffset(n);
    }
    SharedMemoryAllocator &getAllocator() { return allocator; }
    void clear() { array->clear(allocator); }
    iterator begin() { return (iterator)array->begin(allocator); }
    iterator end() { return (iterator)array->end(allocator); }
    bool empty() { return array->empty(); }
    TYPE *operator[](size_t idx) {
      return (TYPE*)allocator.getAddr((*array).at(idx, allocator));
    }
    TYPE *at(size_t idx) {
      return (TYPE*)allocator.getAddr(array->at(idx, allocator));
    }
    void push_back(TYPE *data) { 
      array->push_back(allocator.getOffset(data), allocator); 
    }
    void reserve(size_t s) { array->reserve(s, allocator); }
    void resize(size_t s) { array->resize(s, allocator, (off_t)0); }
    size_t size() { return array->size(); }
    size_t getAllocatedSize() { return array->allocatedSize; }

    ARRAY *array;

    SharedMemoryAllocator allocator;

#ifdef ADVANCED_USE_REMOVED_LIST
    size_t count() { return size() == 0 ? 0 : size() - removedList->size() - 1; }
  protected:
    Vector<size_t>	*removedList;
#endif

  };

#endif	// NGT_SHARED_MEMORY_ALLOCATOR

  class ObjectSpace;

  template <class TYPE>
    class ObjectArray : public vector<TYPE> {
  public:

    size_t push(TYPE &n) {
      if (vector<TYPE>::size() == 0) {
	vector<TYPE>::push_back(0);
      }
      vector<TYPE>::push_back(n);
      return vector<TYPE>::size() - 1;
    }

    size_t insert(TYPE &n) {
#ifdef ADVANCED_USE_REMOVED_LIST
      if (!removedList.empty()) {
	size_t idx = removedList.top();
	removedList.pop();
	put(idx, n);
	return idx;
      }
#endif
      return push(n);
    }

    bool isEmpty(size_t idx) {
      if (idx < vector<TYPE>::size()) {
	return false;
      } else {
	return true;
      }
    }

    void put(size_t idx, TYPE &n) {
      if (vector<TYPE>::size() <= idx) {
	vector<TYPE>::resize(idx + 1);
      }
      (*this)[idx] = n;
    }

    void remove(size_t idx) {
      if (isEmpty(idx)) {
	NGTThrowException("remove: Not in-memory or invalid id");  
      }
    }

    TYPE &get(size_t idx) {
      if (isEmpty(idx)) {
	stringstream msg;
	msg << "get: Not in-memory or invalid offset of node. " << idx << ":" << this->size();
	NGTThrowException(msg.str());
      }
      return (*this)[idx];
    }

    void serialize(ofstream &os, ObjectSpace *objectspace = 0) {
      if (!os.is_open()) {
	NGTThrowException("NGT::Common: Not open the specified stream yet.");
      }
      NGT::Serializer::write(os, vector<TYPE>::size());    
      for (size_t idx = 0; idx < vector<TYPE>::size(); idx++) {
	if ((*this).isEmpty(idx)) {
	  NGT::Serializer::write(os, '-');
	} else {
	  NGT::Serializer::write(os, '+');
	  if (objectspace == 0) {
	    get(idx).serialize(os);
	  } else {
	    get(idx).serialize(os, objectspace);
	  }
	}
      }
    }

    void deserialize(ifstream &is, ObjectSpace *objectspace = 0) {
      if (!is.is_open()) {
	NGTThrowException("NGT::Common: Not open the specified stream yet.");
      }
      vector<TYPE>::clear();
      size_t s;
      NGT::Serializer::read(is, s);
      vector<TYPE>::reserve(s);
      for (size_t i = 0; i < s; i++) {
	char type;
	NGT::Serializer::read(is, type);
	switch(type) {
	case '-':
	  {
	    TYPE dummy;
	    vector<TYPE>::push_back(dummy);
	  }
	  break;
	case '+':
	  {
	    TYPE v;
	    if (objectspace == 0) {
	      v.deserialize(is);
	    } else {
	      v.deserialize(is, objectspace);
	    }
	    vector<TYPE>::push_back(v);
	  }
	  break;
	default:
	  {
	    assert(type == '-' || type == '+');
	    break;
	  }
	}
      }
    }

    void serializeAsText(ofstream &os, ObjectSpace *objectspace = 0) {
      if (!os.is_open()) {
	NGTThrowException("NGT::Common: Not open the specified stream yet.");
      }
      os.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
      os << std::setprecision(8);

      os << vector<TYPE>::size() << endl;
      for (size_t idx = 0; idx < vector<TYPE>::size(); idx++) {
	if ((*this).isEmpty(idx)) {
	  os << idx << " - " << endl;
	} else {
	  os << idx << " + ";
	  if (objectspace == 0) {
	    (*this)[idx].serializeAsText(os);
	  } else {
	    (*this)[idx].serializeAsText(os, objectspace);
	  }
	  os << endl;
	}
      }
      os << fixed;
    }

    void deserializeAsText(ifstream &is, ObjectSpace *objectspace = 0) {
      if (!is.is_open()) {
	NGTThrowException("NGT::Common: Not open the specified stream yet.");
      }
      vector<TYPE>::clear();
      size_t s;
      NGT::Serializer::readAsText(is, s);
      vector<TYPE>::reserve(s);
      for (size_t i = 0; i < s; i++) {
	size_t idx;
	NGT::Serializer::readAsText(is, idx);
	if (i != idx) {
	  cerr << "ObjectArray: Error. index of a specified import file is invalid. " << idx << ":" << i << endl;
	}
	char type;
	NGT::Serializer::readAsText(is, type);
	switch(type) {
	case '-':
	  {
	    TYPE dummy;
	    vector<TYPE>::push_back(dummy);
	  }
	  break;
	case '+':
	  {
	    TYPE v;
	    if (objectspace == 0) {
	      v.deserializeAsText(is);
	    } else {
	      v.deserializeAsText(is, objectspace);
	    }
	    vector<TYPE>::push_back(v);
	  }
	  break;
	default:
	  {
	    assert(type == '-' || type == '+');
	    break;
	  }
	}
      }
    }

    void set(size_t idx, TYPE &n) {
      (*this)[idx] = n;
    }

    void deleteAll() {}

#ifdef ADVANCED_USE_REMOVED_LIST
  protected:
    priority_queue<size_t, vector<size_t>, greater<size_t> >	removedList;
#endif
  };

  template <class TYPE>
    class Repository : public vector<TYPE*> {
  public:

    static TYPE *allocate() { return new TYPE; }

    size_t push(TYPE *n) {
      if (vector<TYPE*>::size() == 0) {
	vector<TYPE*>::push_back(0);
      }
      vector<TYPE*>::push_back(n);
      return vector<TYPE*>::size() - 1;
    }

    size_t insert(TYPE *n) {
#ifdef ADVANCED_USE_REMOVED_LIST
      if (!removedList.empty()) {
	size_t idx = removedList.top();
	removedList.pop();
	put(idx, n);
	return idx;
      }
#endif
      return push(n);
    }

    bool isEmpty(size_t idx) {
      if (idx < vector<TYPE*>::size()) {
	return (*this)[idx] == 0;
      } else {
	return true;
      }
    }

    void put(size_t idx, TYPE *n) {
      if (vector<TYPE*>::size() <= idx) {
	vector<TYPE*>::resize(idx + 1, 0);
      }
      if ((*this)[idx] != 0) {
	NGTThrowException("put: Not empty");  
      }
      (*this)[idx] = n;
    }

    void remove(size_t idx) {
      if (isEmpty(idx)) {
	NGTThrowException("remove: Not in-memory or invalid id");  
      }
      delete (*this)[idx];
      (*this)[idx] = 0;
#ifdef ADVANCED_USE_REMOVED_LIST
      removedList.push(idx);
#endif
    }

    TYPE **getPtr() { return &(*this)[0]; }

    inline TYPE *get(size_t idx) {
      if (isEmpty(idx)) {
	stringstream msg;
	msg << "get: Not in-memory or invalid offset of node. idx=" << idx << " size=" << this->size();
	NGTThrowException(msg.str());
      }
      return (*this)[idx];
    }

    inline TYPE *getWithoutCheck(size_t idx) { return (*this)[idx]; }

    void serialize(ofstream &os, ObjectSpace *objectspace = 0) {
      if (!os.is_open()) {
	NGTThrowException("NGT::Common: Not open the specified stream yet.");
      }
      NGT::Serializer::write(os, vector<TYPE*>::size());    
      for (size_t idx = 0; idx < vector<TYPE*>::size(); idx++) {
	if ((*this)[idx] == 0) {
	  NGT::Serializer::write(os, '-');
	} else {
	  NGT::Serializer::write(os, '+');
	  if (objectspace == 0) {
	    (*this)[idx]->serialize(os);
	  } else {
	    (*this)[idx]->serialize(os, objectspace);
	  }
	}
      }
    }

    void deserialize(ifstream &is, ObjectSpace *objectspace = 0) {
      if (!is.is_open()) {
	NGTThrowException("NGT::Common: Not open the specified stream yet.");
      }
      deleteAll();
      size_t s;
      NGT::Serializer::read(is, s);
      vector<TYPE*>::reserve(s);
      for (size_t i = 0; i < s; i++) {
	char type;
	NGT::Serializer::read(is, type);
	switch(type) {
	case '-':
	  {
	    vector<TYPE*>::push_back(0);
#ifdef ADVANCED_USE_REMOVED_LIST
	    if (i != 0) {
	      removedList.push(i);
	    }
#endif
	  }
	  break;
	case '+':
	  {
	    if (objectspace == 0) {
	      TYPE *v = new TYPE;
	      v->deserialize(is);
	      vector<TYPE*>::push_back(v);
	    } else {
	      TYPE *v = new TYPE(objectspace);
	      v->deserialize(is, objectspace);
	      vector<TYPE*>::push_back(v);
	    }
	  }
	  break;
	default:
	  {
	    assert(type == '-' || type == '+');
	    break;
	  }
	}
      }
    }

    void serializeAsText(ofstream &os, ObjectSpace *objectspace = 0) {
      if (!os.is_open()) {
	NGTThrowException("NGT::Common: Not open the specified stream yet.");
      }
      // The format is almost the same as the default and the best in terms of the string length.
      os.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
      os << std::setprecision(8);

      os << vector<TYPE*>::size() << endl;
      for (size_t idx = 0; idx < vector<TYPE*>::size(); idx++) {
	if ((*this)[idx] == 0) {
	  os << idx << " - " << endl;
	} else {
	  os << idx << " + ";
	  if (objectspace == 0) {
	    (*this)[idx]->serializeAsText(os);
	  } else {
	    (*this)[idx]->serializeAsText(os, objectspace);
	  }
	  os << endl;
	}
      }
      os << fixed;
    }

    void deserializeAsText(ifstream &is, ObjectSpace *objectspace = 0) {
      if (!is.is_open()) {
	NGTThrowException("NGT::Common: Not open the specified stream yet.");
      }
      deleteAll();
      size_t s;
      NGT::Serializer::readAsText(is, s);
      vector<TYPE*>::reserve(s);
      for (size_t i = 0; i < s; i++) {
	size_t idx;
	NGT::Serializer::readAsText(is, idx);
	if (i != idx) {
	  cerr << "Repository: Error. index of a specified import file is invalid. " << idx << ":" << i << endl;
	}
	char type;
	NGT::Serializer::readAsText(is, type);
	switch(type) {
	case '-':
	  {
	    vector<TYPE*>::push_back(0);
#ifdef ADVANCED_USE_REMOVED_LIST
	    if (i != 0) {
	      removedList.push(i);
	    }
#endif
	  }
	  break;
	case '+':
	  {
	    if (objectspace == 0) {
	      TYPE *v = new TYPE;
	      v->deserializeAsText(is);
	      vector<TYPE*>::push_back(v);
	    } else {
	      TYPE *v = new TYPE(objectspace);
	      v->deserializeAsText(is, objectspace);
	      vector<TYPE*>::push_back(v);
	    }
	  }
	  break;
	default:
	  {
	    assert(type == '-' || type == '+');
	    break;
	  }
	}
      }
    }

    void deleteAll() {
      for (size_t i = 0; i < this->size(); i++) {
	if ((*this)[i] != 0) {
	  delete (*this)[i];
	  (*this)[i] = 0;
	}
      }
      this->clear();
#ifdef ADVANCED_USE_REMOVED_LIST
      while(!removedList.empty()){ removedList.pop(); };
#endif
    }

    void set(size_t idx, TYPE *n) {
      (*this)[idx] = n;
    }

#ifdef ADVANCED_USE_REMOVED_LIST
    size_t count() { return vector<TYPE*>::size() == 0 ? 0 : vector<TYPE*>::size() - removedList.size() - 1; }
  protected:
    priority_queue<size_t, vector<size_t>, greater<size_t> >	removedList;
#endif
  };

#pragma pack(2)
  class ObjectDistance {
  public:
    ObjectDistance():id(0), distance(0.0) {}
    ObjectDistance(unsigned int i, float d):id(i), distance(d) {}
    inline bool operator==(const ObjectDistance &o) const {
      return (distance == o.distance) && (id == o.id);
    }
    inline void set(unsigned int i, float d) { id = i; distance = d; }
    inline bool operator<(const ObjectDistance &o) const {
      if (distance == o.distance) {
        return id < o.id;
      } else {
        return distance < o.distance;
      }
    }
    inline bool operator>(const ObjectDistance &o) const {
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
    uint32_t            id;
    float         distance;
  };

#pragma pack()

  class Object;
  class ObjectDistances;

  class Container {
  public:
    Container(Object &o, ObjectID i):object(o), id(i) {}
    Container(Container &c):object(c.object), id(c.id) {}
    Object		&object;
    ObjectID		id;
  };

  typedef priority_queue<ObjectDistance, vector<ObjectDistance>, less<ObjectDistance> > ResultPriorityQueue;

  class SearchContainer : public NGT::Container {
  public:
    SearchContainer(Object &f, ObjectID i): Container(f, i) { initialize(); }
    SearchContainer(Object &f): Container(f, 0) { initialize(); }
    SearchContainer(SearchContainer &sc): Container(sc) { *this = sc; }
    SearchContainer &operator=(SearchContainer &sc) {
      size = sc.size;
      radius = sc.radius;
      explorationCoefficient = sc.explorationCoefficient;
      result = sc.result;
      distanceComputationCount = sc.distanceComputationCount;
      edgeSize = sc.edgeSize;
      workingResult = sc.workingResult;
      useAllNodesInLeaf = sc.useAllNodesInLeaf;  
      return *this;
    }
    virtual ~SearchContainer() {}
    virtual void initialize() {
      size = 10;
      radius = FLT_MAX;
      explorationCoefficient = 1.1;
      result = 0;
      edgeSize = -1;	// dynamically prune the edges during search. -1 means following the index property. 0 means using all edges.
      useAllNodesInLeaf = false;
    }
    void setSize(size_t s) { size = s; };
    void setResults(ObjectDistances *r) { result = r; }
    void setRadius(Distance r) { radius = r; }
    void setEpsilon(float e) { explorationCoefficient = e + 1.0; }
    void setEdgeSize(int e) { edgeSize = e; }

    inline bool resultIsAvailable() { return result != 0; }
    ObjectDistances &getResult() {
      if (result == 0) {
	NGTThrowException("Inner error: results is not set");
      }
      return *result;
    }

    ResultPriorityQueue &getWorkingResult() { return workingResult; }


    size_t		size;
    Distance		radius;
    float		explorationCoefficient;
    int			edgeSize;
    size_t		distanceComputationCount;
    ResultPriorityQueue	workingResult;
    bool		useAllNodesInLeaf;

  private:
    ObjectDistances	*result;
  };


  class InsertContainer : public Container {
  public:
    InsertContainer(Object &f, ObjectID i):Container(f, i) {}
  };

  class Timer {
  public:
  Timer():time(0) {}

    void reset() { time = 0; ntime = 0; }

    void start() {
      struct timespec res;
      clock_getres(CLOCK_REALTIME, &res);
      reset();
      clock_gettime(CLOCK_REALTIME, &startTime);
    }

    void restart() {
      clock_gettime(CLOCK_REALTIME, &startTime);
    }

    void stop() {
      clock_gettime(CLOCK_REALTIME, &stopTime);
      sec = stopTime.tv_sec - startTime.tv_sec;
      nsec = stopTime.tv_nsec - startTime.tv_nsec;
      if (nsec < 0) {
	sec -= 1;
	nsec += 1000000000L;
      }
      time += (double)sec + (double)nsec / 1000000000.0;
      ntime += sec * 1000000000L + nsec;
    }

    friend ostream &operator<<(ostream &os, Timer &t) {
      os << setprecision(6) << t.time << " (sec)";
      return os;
    }

    struct timespec startTime;
    struct timespec stopTime;

    int64_t	sec;
    int64_t	nsec;
    int64_t	ntime;	// nano second
    double      time;	// second
  };

} // namespace NGT

