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
#include	<map>
#include	<set>
#include	<bitset>
#include	<iomanip>
#include	<unordered_set>

#include	<sys/time.h>
#include	<sys/stat.h>
#include	<stdint.h>

#include	"NGT/defines.h"
#include	"NGT/Common.h"
#include	"NGT/Tree.h"
#include	"NGT/Thread.h"
#include	"NGT/Graph.h"


namespace NGT {

  class Property;

  class Index {
  public:

    class Property {
    public:
      typedef ObjectSpace::ObjectType		ObjectType;
      typedef ObjectSpace::DistanceType		DistanceType;
      typedef NeighborhoodGraph::SeedType	SeedType;
      typedef NeighborhoodGraph::GraphType	GraphType;
      enum ObjectAlignment {
	ObjectAlignmentNone	= 0,
	ObjectAlignmentTrue	= 1,
	ObjectAlignmentFalse	= 2
      };
      enum IndexType {
	IndexTypeNone		= 0,
	GraphAndTree		= 1,
	Graph			= 2
      };
      enum DatabaseType {
	DatabaseTypeNone	= 0,
	Memory			= 1,
	MemoryMappedFile	= 2
      };
      Property() { setDefault(); }
      void setDefault() {
	dimension 	= 0;
	threadPoolSize	= 32;
	objectType	= ObjectSpace::ObjectType::Float;
	distanceType	= DistanceType::DistanceTypeL2;
	indexType	= IndexType::GraphAndTree;
	objectAlignment	= ObjectAlignment::ObjectAlignmentFalse;
	pathAdjustmentInterval = 0;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	databaseType	= DatabaseType::MemoryMappedFile;
      	graphSharedMemorySize	= 512; // MB
      	treeSharedMemorySize	= 512; // MB
      	objectSharedMemorySize	= 512; // MB  512 is up to 20M objects.
#else
	databaseType	= DatabaseType::Memory;
#endif
	prefetchOffset	= 0;
      }
      void clear() {
	dimension 	= -1;
	threadPoolSize	= -1;
	objectType	= ObjectSpace::ObjectTypeNone;
	distanceType	= DistanceType::DistanceTypeNone;
	indexType	= IndexTypeNone;
	databaseType	= DatabaseTypeNone;
	objectAlignment	= ObjectAlignment::ObjectAlignmentNone;
	pathAdjustmentInterval	= -1;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      	graphSharedMemorySize	= -1;
      	treeSharedMemorySize	= -1;
      	objectSharedMemorySize	= -1;
#endif
	prefetchOffset	= -1;
      }

      void exportProperty(NGT::PropertySet &p) {
	p.set("Dimension", dimension);
	p.set("ThreadPoolSize", threadPoolSize);
	switch (objectType) {
	case ObjectSpace::ObjectType::Uint8: p.set("ObjectType", "Integer-1"); break;
	case ObjectSpace::ObjectType::Float: p.set("ObjectType", "Float-4"); break;
	default : cerr << "Fatal error. Invalid object type. " << objectType <<endl; abort();
	}
	switch (distanceType) {
	case DistanceType::DistanceTypeNone:			p.set("DistanceType", "None"); break;
	case DistanceType::DistanceTypeL1:			p.set("DistanceType", "L1"); break;
	case DistanceType::DistanceTypeL2:			p.set("DistanceType", "L2"); break;
	case DistanceType::DistanceTypeHamming:			p.set("DistanceType", "Hamming"); break;
	case DistanceType::DistanceTypeAngle:			p.set("DistanceType", "Angle"); break;
	case DistanceType::DistanceTypeCosine:			p.set("DistanceType", "Cosine"); break;
	case DistanceType::DistanceTypeNormalizedAngle:		p.set("DistanceType", "NormalizedAngle"); break;
	case DistanceType::DistanceTypeNormalizedCosine:	p.set("DistanceType", "NormalizedCosine"); break;
	default : cerr << "Fatal error. Invalid distance type. " << distanceType << endl; abort();
	}
	switch (indexType) {      
	case IndexType::GraphAndTree:	p.set("IndexType", "GraphAndTree"); break;
	case IndexType::Graph:		p.set("IndexType", "Graph"); break;
	default : cerr << "Fatal error. Invalid index type. " << indexType << endl; abort();
	}
	switch (databaseType) {
	case DatabaseType::Memory:		p.set("DatabaseType", "Memory"); break;
	case DatabaseType::MemoryMappedFile:	p.set("DatabaseType", "MemoryMappedFile"); break;
	default : cerr << "Fatal error. Invalid database type. " << databaseType << endl; abort();
	}
	switch (objectAlignment) {
	case ObjectAlignment::ObjectAlignmentNone:	p.set("ObjectAlignment", "None"); break;
	case ObjectAlignment::ObjectAlignmentTrue:	p.set("ObjectAlignment", "True"); break;
	case ObjectAlignment::ObjectAlignmentFalse:	p.set("ObjectAlignment", "False"); break;
	default : cerr << "Fatal error. Invalid objectAlignment. " << objectAlignment << endl; abort();
	}
	p.set("PathAdjustmentInterval", pathAdjustmentInterval);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	p.set("GraphSharedMemorySize", graphSharedMemorySize);
	p.set("TreeSharedMemorySize", treeSharedMemorySize);
	p.set("ObjectSharedMemorySize", objectSharedMemorySize);
#endif
	p.set("PrefetchOffset", prefetchOffset);
      }

      void importProperty(NGT::PropertySet &p) {
	setDefault();
	dimension = p.getl("Dimension", dimension);
	threadPoolSize = p.getl("ThreadPoolSize", threadPoolSize);
	PropertySet::iterator it = p.find("ObjectType");
	if (it != p.end()) {
	  if (it->second == "Float-4") {
	    objectType = ObjectSpace::ObjectType::Float;
	  } else if (it->second == "Integer-1") {
	    objectType = ObjectSpace::ObjectType::Uint8;
	  } else {
	    cerr << "Invalid Object Type in the property. " << it->first << ":" << it->second << endl;
	  }
	} else {
	  cerr << "Not found \"ObjectType\"" << endl;
	}
	it = p.find("DistanceType");
	if (it != p.end()) {
	  if (it->second == "None") {
	    distanceType = DistanceType::DistanceTypeNone;
	  } else if (it->second == "L1") {
	    distanceType = DistanceType::DistanceTypeL1;
	  } else if (it->second == "L2") {
	    distanceType = DistanceType::DistanceTypeL2;
	  } else if (it->second == "Hamming") {
	    distanceType = DistanceType::DistanceTypeHamming;
	  } else if (it->second == "Angle") {
	    distanceType = DistanceType::DistanceTypeAngle;
	  } else if (it->second == "Cosine") {
	    distanceType = DistanceType::DistanceTypeCosine;
	  } else if (it->second == "NormalizedAngle") {
	    distanceType = DistanceType::DistanceTypeNormalizedAngle;
	  } else if (it->second == "NormalizedCosine") {
	    distanceType = DistanceType::DistanceTypeNormalizedCosine;
	  } else {
	    cerr << "Invalid Distance Type in the property. " << it->first << ":" << it->second << endl;
	  }
	} else {
	  cerr << "Not found \"DistanceType\"" << endl;
	}
	it = p.find("IndexType");
	if (it != p.end()) {
	  if (it->second == "GraphAndTree") {
	    indexType = IndexType::GraphAndTree;
	  } else if (it->second == "Graph") {
	    indexType = IndexType::Graph;
	  } else {
	    cerr << "Invalid Index Type in the property. " << it->first << ":" << it->second << endl;
	  }
	} else {
	  cerr << "Not found \"IndexType\"" << endl;
	}
	it = p.find("DatabaseType");
	if (it != p.end()) {
	  if (it->second == "Memory") {
	    databaseType = DatabaseType::Memory;
	  } else if (it->second == "MemoryMappedFile") {
	    databaseType = DatabaseType::MemoryMappedFile;
	  } else {
	    cerr << "Invalid Database Type in the property. " << it->first << ":" << it->second << endl;
	  }
	} else {
	  cerr << "Not found \"DatabaseType\"" << endl;
	}
	it = p.find("ObjectAlignment");
	if (it != p.end()) {
	  if (it->second == "None") {
	    objectAlignment = ObjectAlignment::ObjectAlignmentNone;
	  } else if (it->second == "True") {
	    objectAlignment = ObjectAlignment::ObjectAlignmentTrue;
	  } else if (it->second == "False") {
	    objectAlignment = ObjectAlignment::ObjectAlignmentFalse;
	  } else {
	    cerr << "Invalid Object Alignment in the property. " << it->first << ":" << it->second << endl;
	  }
	} else {
	  cerr << "Not found \"ObjectAlignment\"" << endl;
	  objectAlignment = ObjectAlignment::ObjectAlignmentFalse;
	}
	pathAdjustmentInterval  = p.getl("PathAdjustmentInterval", pathAdjustmentInterval);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	graphSharedMemorySize  = p.getl("GraphSharedMemorySize", graphSharedMemorySize);
	treeSharedMemorySize   = p.getl("TreeSharedMemorySize", treeSharedMemorySize);
	objectSharedMemorySize = p.getl("ObjectSharedMemorySize", objectSharedMemorySize);
#endif
	prefetchOffset = p.getl("PrefetchOffset", prefetchOffset);
      }

      void set(NGT::Property &prop);
      void get(NGT::Property &prop);
      int		dimension;
      int		threadPoolSize;
      ObjectSpace::ObjectType	objectType;
      DistanceType	distanceType;
      IndexType		indexType;
      DatabaseType	databaseType;
      ObjectAlignment	objectAlignment;
      int		pathAdjustmentInterval;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      int		graphSharedMemorySize;
      int		treeSharedMemorySize;
      int		objectSharedMemorySize;
#endif
      int		prefetchOffset;
    };

    class InsertionResult {
    public:
      InsertionResult():id(0), identical(false), distance(0.0) {}
      InsertionResult(size_t i, bool tf, Distance d):id(i), identical(tf), distance(d) {}
      size_t	id;
      bool	identical;
      Distance	distance; // the distance between the centroid and the inserted object.
    };

    Index():index(0) {}
    Index(const string &database, bool rdOnly = false):index(0) { open(database, rdOnly); }
    Index(const string &database, NGT::Property &prop):index(0) { open(database, prop);  }
    virtual ~Index() { close(); }

    void open(const string &database, NGT::Property &prop) {
      open(database);
      setProperty(prop);
    }
    void open(const string &database, bool rdOnly = false);

    void close() {
      if (index != 0) { 
	delete index;
	index = 0;
      } 
      path.clear();
    }
    void save() {
      if (path.empty()) {
	NGTThrowException("NGT::Index::saveIndex: path is empty");	
      }
      saveIndex(path);
    }
    static void mkdir(const string &dir) { 
      if (::mkdir(dir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP |  S_IROTH | S_IXOTH) != 0) {
	stringstream msg;
	msg << "NGT::Index::mkdir: Cannot make the specified directory. " << dir;
	NGTThrowException(msg);	
      }
    }
    static void createGraphAndTree(const string &database, NGT::Property &prop, const string &dataFile, size_t dataSize = 0);
    static void createGraphAndTree(const string &database, NGT::Property &prop) { createGraphAndTree(database, prop, ""); }
    static void createGraph(const string &database, NGT::Property &prop, const string &dataFile, size_t dataSize = 0);
    template<typename T> size_t insert(vector<T> &object);
    template<typename T> size_t append(vector<T> &object);
    static void append(const string &database, const string &dataFile, size_t threadSize, size_t dataSize); 
    static void append(const string &database, const float *data, size_t dataSize, size_t threadSize);
    static void remove(const string &database, vector<ObjectID> &objects, bool force = false);
    static void exportIndex(const string &database, const string &file);
    static void importIndex(const string &database, const string &file);
    virtual void load(const string &ifile, size_t dataSize) { getIndex().load(ifile, dataSize); }
    virtual void append(const string &ifile, size_t dataSize) { getIndex().append(ifile, dataSize); }
    virtual void append(const float *data, size_t dataSize) { getIndex().append(data, dataSize); } 
    virtual void append(const double *data, size_t dataSize) { getIndex().append(data, dataSize); } 
    virtual size_t getObjectRepositorySize() { return getIndex().getObjectRepositorySize(); }
    virtual void createIndex(size_t threadNumber) { getIndex().createIndex(threadNumber); }
    virtual void saveIndex(const string &ofile) { getIndex().saveIndex(ofile); }
    virtual void loadIndex(const string &ofile) { getIndex().loadIndex(ofile); }
    virtual Object *allocateObject(const string &textLine, const string &sep) { return getIndex().allocateObject(textLine, sep); }
    virtual Object *allocateObject(const vector<double> &obj) { return getIndex().allocateObject(obj); }
    virtual Object *allocateObject(const vector<float> &obj) { return getIndex().allocateObject(obj); }
    virtual Object *allocateObject(const float *obj, size_t size) { return getIndex().allocateObject(obj, size); }
    virtual size_t getSizeOfElement() { return getIndex().getSizeOfElement(); }
    virtual void setProperty(NGT::Property &prop) { getIndex().setProperty(prop); }
    virtual void getProperty(NGT::Property &prop) { getIndex().getProperty(prop); }
    virtual void deleteObject(Object *po) { getIndex().deleteObject(po); }
    virtual void linearSearch(NGT::SearchContainer &sc) { getIndex().linearSearch(sc); }
    virtual void search(NGT::SearchContainer &sc) { getIndex().search(sc); }
    virtual void search(NGT::SearchContainer &sc, ObjectDistances &seeds) { getIndex().search(sc, seeds); }
    virtual void remove(ObjectID id, bool force = false) { getIndex().remove(id, force); }
    virtual void exportIndex(const string &file) { getIndex().exportIndex(file); }
    virtual void importIndex(const string &file) { getIndex().importIndex(file); }
    virtual bool verify(vector<uint8_t> &status, bool info = false) { return getIndex().verify(status, info); }
    virtual ObjectSpace &getObjectSpace() { return getIndex().getObjectSpace(); }
    virtual size_t getSharedMemorySize(ostream &os, SharedMemoryAllocator::GetMemorySizeType t = SharedMemoryAllocator::GetTotalMemorySize) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      size_t osize = getObjectSpace().getRepository().getAllocator().getMemorySize(t);
#else
      size_t osize = 0;
#endif
      os << "object=" << osize << endl;
      size_t isize = getIndex().getSharedMemorySize(os, t); 
      return osize + isize;
    }
    void searchUsingOnlyGraph(NGT::SearchContainer &sc) { 
      sc.distanceComputationCount = 0;
      ObjectDistances seeds; 
      getIndex().search(sc, seeds); 
    }
    Index &getIndex() {
      if (index == 0) {
	assert(index != 0);
	NGTThrowException("NGT::Index::getIndex: Index is unavailable.");	
      }
      return *index;
    }

    static void destroy(const string &path) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      std::remove(string(path + "/grp").c_str());
      std::remove(string(path + "/grpc").c_str());
      std::remove(string(path + "/trei").c_str());
      std::remove(string(path + "/treic").c_str());
      std::remove(string(path + "/trel").c_str());
      std::remove(string(path + "/trelc").c_str());
      std::remove(string(path + "/objpo").c_str());
      std::remove(string(path + "/objpoc").c_str());
#else
      std::remove(string(path + "/grp").c_str());
      std::remove(string(path + "/tre").c_str());
      std::remove(string(path + "/obj").c_str());
#endif
      std::remove(string(path + "/prf").c_str());
      std::remove(path.c_str());
    }
    
    static void version(ostream &os);
    string getPath(){ return path; }
  protected:
    static void loadAndCreateIndex(Index &index, const string &database, const string &dataFile,
				   size_t threadSize, size_t dataSize);

    Index *index;
    string path;
  };

  class GraphIndex : public Index, 
    public NeighborhoodGraph {
  public:

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    GraphIndex(const string &allocator, bool rdOnly = false);
    GraphIndex(const string &allocator, NGT::Property &prop):readOnly(false) {
      initialize(allocator, prop);
    }
    void initialize(const string &allocator, NGT::Property &prop);
#else // NGT_SHARED_MEMORY_ALLOCATOR
    GraphIndex(const string &database, bool rdOnly = false);
    GraphIndex(NGT::Property &prop):readOnly(false) {
      initialize(prop);
    }

    void initialize(NGT::Property &prop) {
      constructObjectSpace(prop);
      setProperty(prop);
    }

#endif // NGT_SHARED_MEMORY_ALLOCATOR

    virtual ~GraphIndex() {
      destructObjectSpace();
    }
    void constructObjectSpace(NGT::Property &prop);

    void destructObjectSpace() {
      if (objectSpace == 0) {
	return;
      }
      if (property.objectType == NGT::ObjectSpace::ObjectType::Float) {
	ObjectSpaceRepository<float, double> *os = (ObjectSpaceRepository<float, double>*)objectSpace;
#ifndef NGT_SHARED_MEMORY_ALLOCATOR
	os->deleteAll();
#endif
	delete os;
      } else if (property.objectType == NGT::ObjectSpace::ObjectType::Uint8) {
	ObjectSpaceRepository<unsigned char, int> *os = (ObjectSpaceRepository<unsigned char, int>*)objectSpace;
#ifndef NGT_SHARED_MEMORY_ALLOCATOR
	os->deleteAll();
#endif
	delete os;
      } else {
	cerr << "Cannot find Object Type in the property. " << property.objectType << endl;
	return;
      }
      objectSpace = 0;
    }

    virtual void load(const string &ifile, size_t dataSize = 0) {
      if (ifile.empty()) {
	return;
      }
      istream *is;
      ifstream *ifs = 0;
      if (ifile == "-") {
	is = &cin;
      } else {
	ifs = new ifstream;
	ifs->ifstream::open(ifile);
	if (!(*ifs)) {
	  stringstream msg;
	  msg << "Index::load: Cannot open the specified file. " << ifile;
	  NGTThrowException(msg);
	}
	is = ifs;
      }
      try {
	objectSpace->readText(*is, dataSize);
      } catch(Exception &err) {
	if (ifile != "-") {
	  delete ifs;
	}
	throw(err);
      }
      if (ifile != "-") {
	delete ifs;
      }
    }

    virtual void append(const string &ifile, size_t dataSize = 0) {
      ifstream is(ifile.c_str());
      objectSpace->appendText(is, dataSize);
    }

    virtual void append(const float *data, size_t dataSize) { objectSpace->append(data, dataSize); }
    virtual void append(const double *data, size_t dataSize) { objectSpace->append(data, dataSize); }

    virtual void saveIndex(const string &ofile) {
#ifndef NGT_SHARED_MEMORY_ALLOCATOR
      try {
	mkdir(ofile);
      } catch(...) {}
      if (objectSpace != 0) {
	objectSpace->serialize(ofile + "/obj");
      } else {
	cerr << "saveIndex::Warning! ObjectSpace is null. continue saving..." << endl;
      }
      string fname = ofile + "/grp";
      ofstream osg(fname);
      if (!osg.is_open()) {
	stringstream msg;
	msg << "saveIndex:: Cannot open. " << fname;
	NGTThrowException(msg);
      }
      repository.serialize(osg);
#endif
      saveProperty(ofile);
    }

    void saveProperty(const string &file) {
      NGT::PropertySet prop;
      assert(property.dimension != 0);
      GraphIndex::property.exportProperty(prop);
      NeighborhoodGraph::property.exportProperty(prop);
      prop.save(file + "/prf");
    }

    void exportProperty(const string &file) {
      NGT::PropertySet prop;
      assert(property.dimension != 0);
      GraphIndex::property.exportProperty(prop);
      NeighborhoodGraph::property.exportProperty(prop);
      prop.save(file + "/prf");
    }

    virtual void loadIndex(const string &ifile, bool readOnly);

    virtual void exportIndex(const string &ofile) {
      try {
	mkdir(ofile);
      } catch(...) {
	stringstream msg;
	msg << "exportIndex:: Cannot make the directory. " << ofile;
	NGTThrowException(msg);
      }
      objectSpace->serializeAsText(ofile + "/obj");
      ofstream osg(ofile + "/grp");
      repository.serializeAsText(osg);
      exportProperty(ofile);
    }

    virtual void importIndex(const string &ifile) {
      objectSpace->deserializeAsText(ifile + "/obj");
      string fname = ifile + "/grp";
      ifstream isg(fname);
      if (!isg.is_open()) {
	stringstream msg;
	msg << "importIndex:: Cannot open. " << fname;
	NGTThrowException(msg);
      }
      repository.deserializeAsText(isg);
    }

    void linearSearch(NGT::SearchContainer &sc) {
      ObjectSpace::ResultSet results;
      objectSpace->linearSearch(sc.object, sc.radius, sc.size, results);
      ObjectDistances &qresults = sc.getResult();
      qresults.moveFrom(results);
    }

    // GraphIndex
    virtual void search(NGT::SearchContainer &sc) {
      sc.distanceComputationCount = 0;
      ObjectDistances seeds;
      search(sc, seeds);
    }

    // get randomly nodes as seeds.
    template<class REPOSITORY> void getRandomSeeds(REPOSITORY &repo, ObjectDistances &seeds, size_t seedSize) {
      // clear all distances to find the same object as a randomized object.
      for (ObjectDistances::iterator i = seeds.begin(); i != seeds.end(); i++) {
	(*i).distance = 0.0;
      }
      size_t repositorySize = repo.size();
      repositorySize = repositorySize == 0 ? 0 : repositorySize - 1; // Because the head of repository is a dummy.
      seedSize = seedSize > repositorySize ? repositorySize : seedSize;
      vector<ObjectID> deteted;
      while (seedSize > seeds.size()) {
	double random = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
	size_t idx = floor(repositorySize * random) + 1;
	if (repo.isEmpty(idx)) {
	  continue;
	}
	ObjectDistance obj(idx, 0.0);
	if (find(seeds.begin(), seeds.end(), obj) != seeds.end()) {
	  continue;
	}
	seeds.push_back(obj);
      }
    }

    void remove(const ObjectID id, bool force) {
      removeEdgesReliably(id);
      try {
	getObjectRepository().remove(id);
      } catch(Exception &err) {
	cerr << "NGT::GraphIndex::remove:: cannot remove from feature. id=" << id << " " << err.what() << endl;
	throw err;
      }
    }

    virtual void searchForNNGInsertion(Object &po, ObjectDistances &result) {
      NGT::SearchContainer sc(po);
      sc.setResults(&result);
      sc.size = NeighborhoodGraph::property.edgeSizeForCreation;
      sc.radius = FLT_MAX;
      sc.explorationCoefficient = NeighborhoodGraph::property.insertionRadiusCoefficient;
      try {
	GraphIndex::search(sc);
      } catch(Exception &err) {
	throw err;
      }
      if (static_cast<int>(result.size()) < NeighborhoodGraph::property.edgeSizeForCreation && 
	  result.size() < repository.size()) {
	if (sc.edgeSize != 0) {
	  sc.edgeSize = 0;	// not prune edges.
	  try {
	    GraphIndex::search(sc);
	  } catch(Exception &err) {
	    throw err;
	  }
	}
      }
    }

    void searchForKNNGInsertion(Object &po, ObjectID id, ObjectDistances &result) {
      double radius = FLT_MAX;
      size_t size = NeighborhoodGraph::property.edgeSizeForCreation;
      if (id > 0) {
	size = NeighborhoodGraph::property.edgeSizeForCreation + 1;
      }
      ObjectSpace::ResultSet rs;
      objectSpace->linearSearch(po, radius, size, rs);
      result.moveFrom(rs, id);
      if ((size_t)NeighborhoodGraph::property.edgeSizeForCreation != result.size()) {
	cerr << "searchForKNNGInsert::Warning! inconsistency of the sizes. ID=" << id 
	     << " " << NeighborhoodGraph::property.edgeSizeForCreation << ":" << result.size() << endl;
	for (size_t i = 0; i < result.size(); i++) {
	  cerr << result[i].id << ":" << result[i].distance << " ";
	}
	cerr << endl;
      }
    }

    virtual void insert(
			ObjectID id
			) {
      ObjectRepository &fr = objectSpace->getRepository();
      if (fr[id] == 0) {
	cerr << "NGTIndex::insert empty " << id << endl;
	return;
      }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      Object &po = *objectSpace->allocateObject(*fr[id]);
#else
      Object &po = *fr[id];
#endif
      ObjectDistances rs;
      if (NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeANNG) {
	searchForNNGInsertion(po, rs);
      } else {
	searchForKNNGInsertion(po, id, rs);
      }
      insertNode(id, rs);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      objectSpace->deleteObject(&po);
#endif
    }

    virtual void createIndex();
    virtual void createIndex(size_t threadNumber);

    void checkGraph()
    {
      GraphRepository &repo = repository;
      ObjectRepository &fr = objectSpace->getRepository();
      for (size_t id = 0; id < fr.size(); id++){
	if (repo[id] == 0) {
	  cerr << id << " empty" << endl;
	  continue;
	}
	if ((id % 10000) == 0) {
	  cerr << "checkGraph: Processed size=" << id << endl;
	}
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	Object &po = *objectSpace->allocateObject(*fr[id]);
#else
	Object &po = *fr[id];
#endif
	GraphNode *objects = getNode(id);

	ObjectDistances rs;
	NeighborhoodGraph::property.edgeSizeForCreation = objects->size() + 1;
	searchForNNGInsertion(po, rs);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	objectSpace->deleteObject(&po);
#endif

	if (rs.size() != objects->size()) {
	  cerr << "Cannot get the specified number of the results. " << rs.size() << ":" << objects->size() << endl;
	}
	size_t count = 0;
	ObjectDistances::iterator rsi = rs.begin(); 
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	for (GraphNode::iterator ri = objects->begin(repo.allocator);
	     ri != objects->end(repo.allocator) && rsi != rs.end();) {
#else
	for (GraphNode::iterator ri = objects->begin();
	     ri != objects->end() && rsi != rs.end();) {
#endif
	  if ((*ri).distance == (*rsi).distance && (*ri).id == (*rsi).id) {
	    count++;
	    ri++;
	    rsi++;
	  } else if ((*ri).distance < (*rsi).distance) {
	    ri++;
	  } else {
	    rsi++;
	  }
	}
	if (count != objects->size()) {
	  cerr << "id=" << id << " identities=" << count << " " << objects->size() << " " << rs.size() << endl;
	}
      }
    }

    virtual bool verify(vector<uint8_t> &status, bool info)
    {
      bool valid = true;
      cerr << "Started verifying graph and objects" << endl;
      GraphRepository &repo = repository;
      ObjectRepository &fr = objectSpace->getRepository();
      if (repo.size() != fr.size()) {
	if (info) {
	  cerr << "Warning! # of nodes is different from # of objects. " << repo.size() << ":" << fr.size() << endl;
	}
      }
      status.clear();
      status.resize(fr.size(), 0);
      for (size_t id = 1; id < fr.size(); id++){
	status[id] |= repo[id] != 0 ? 0x02 : 0x00;
	status[id] |= fr[id]   != 0 ? 0x01 : 0x00;
      }
      for (size_t id = 1; id < fr.size(); id++) {
	if (fr[id] == 0) {
	  if (id < repo.size() && repo[id] != 0) {
	    cerr << "Error! The node exists in the graph, but the object does not exist. " << id << endl;
	    valid = false;
	  }
	}
	if (fr[id] != 0 && repo[id] == 0) {
	  cerr << "Error. No." << id << " is not registerd in the graph." << endl;
	  valid = false;
	}
	if ((id % 1000000) == 0) {
	  cerr << "  verified " << id << " entries." << endl;
	}
	if (fr[id] != 0) {
	  try {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	    Object *po = objectSpace->allocateObject(*fr[id]);
#else
	    Object *po = fr[id];
#endif
	    if (po == 0) {
	      cerr << "Error! Cannot get the object. " << id << endl;
	      valid = false;
	      continue;
	    }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	    objectSpace->deleteObject(po);
#endif
	  } catch (Exception &err) {
	    cerr << "Error! Cannot get the object. " << id << " " << err.what() << endl;
	    valid = false;
	    continue;
	  }
	}
	if (id >= repo.size()) {
	  cerr << "Error. No." << id << " is not registerd in the object repository. " << repo.size() << endl;
	  valid = false;
	}
	if (id < repo.size() && repo[id] != 0) {
	  try {
	    GraphNode *objects = getNode(id);
	    if (objects == 0) {
	      cerr << "Error! Cannot get the node. " << id << endl;
	      valid = false;
	    }
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    for (GraphNode::iterator ri = objects->begin(repo.allocator);
		 ri != objects->end(repo.allocator); ++ri) {
#else
	    for (GraphNode::iterator ri = objects->begin();
		 ri != objects->end(); ++ri) {
#endif
	      if ((*ri).id == 0 || (*ri).id >= repo.size()) {
		cerr << "Error! Neighbor's ID of the node is wrong." << endl;
		valid = false;
	      }
	      if ((*ri).distance < 0.0) {
		cerr << "Error! Neighbor's distance is munus." << endl;
		valid = false;
	      }
	    }
	  } catch (Exception &err) {
	    cerr << "Error! Cannot get the node. " << id << " " << err.what() << endl;
	    valid = false;
	  }
	}
      }
      return valid;
    }

    static bool showStatisticsOfGraph(NGT::GraphIndex &outGraph, char mode = '-', size_t edgeSize = UINT_MAX)
    {
      long double distance = 0.0;
      size_t numberOfNodes = 0;
      size_t numberOfOutdegree = 0;
      size_t numberOfNodesWithoutEdges = 0;
      size_t maxNumberOfOutdegree = 0;
      size_t minNumberOfOutdegree = SIZE_MAX;
      vector<int64_t> indegreeCount;
      vector<size_t> outdegreeHistogram;
      vector<size_t> indegreeHistogram;
      vector<vector<float> > indegree;
      NGT::GraphRepository &graph = outGraph.repository;
      NGT::ObjectRepository &repo = outGraph.objectSpace->getRepository();
      indegreeCount.resize(graph.size(), 0);
      indegree.resize(graph.size());
      size_t removedObjectCount = 0;
      bool valid = true;
      for (size_t id = 1; id < graph.size(); id++) {
	if (repo[id] == 0) {
	  removedObjectCount++;
	  continue;
	}
	NGT::GraphNode *node = 0;
	try {
	  node = outGraph.getNode(id);
	} catch(NGT::Exception &err) {
	  cerr << "ngt info: Error. Cannot get the node. ID=" << id << ":" << err.what() << endl;
	  valid = false;
	  continue;
	}
	numberOfNodes++;
	if (numberOfNodes % 1000000 == 0) {
	  cerr << "Processed " << numberOfNodes << endl;
	}
	size_t esize = node->size() > edgeSize ? edgeSize : node->size();
	if (esize == 0) {
	  numberOfNodesWithoutEdges++;
	}
	if (esize > maxNumberOfOutdegree) {
	  maxNumberOfOutdegree = esize;
	}
	if (esize < minNumberOfOutdegree) {
	  minNumberOfOutdegree = esize;
	}
	if (outdegreeHistogram.size() <= esize) {
	  outdegreeHistogram.resize(esize + 1);
	}
	outdegreeHistogram[esize]++;
	if (mode == 'e') {
	  cout << id << "," << esize << ": ";
	}
	for (size_t i = 0; i < esize; i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  NGT::ObjectDistance &n = (*node).at(i, graph.allocator);
#else
	  NGT::ObjectDistance &n = (*node)[i];
#endif
	  if (n.id == 0) {
	    cerr << "ngt info: Warning. id is zero." << endl;
	    valid = false;
	  }
	  indegreeCount[n.id]++;
	  indegree[n.id].push_back(n.distance);
	  numberOfOutdegree++;
	  double d = n.distance;
	  if (mode == 'e') {
	    cout << n.id << ":" << d << " ";
	  }
	  distance += d;
	}
	if (mode == 'e') {
	  cout << endl;
	}
      }

      if (mode == 'a') {
	size_t count = 0;
	for (size_t id = 1; id < graph.size(); id++) {
	  if (repo[id] == 0) {
	    continue;
	  }
	  NGT::GraphNode *n = 0;
	  try {
	    n = outGraph.getNode(id);
	  } catch(NGT::Exception &err) {
	    continue;
	  }
	  NGT::GraphNode &node = *n;
	  for (size_t i = 0; i < node.size(); i++) {  
	    NGT::GraphNode *nn = 0;
	    try {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	      nn = outGraph.getNode(node.at(i, graph.allocator).id);
#else
	      nn = outGraph.getNode(node[i].id);
#endif
	    } catch(NGT::Exception &err) {
	      count++;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	      cerr << "Directed edge! " << id << "->" << node.at(i, graph.allocator).id << " no object. " 
		   << node.at(i, graph.allocator).id << endl; 
#else
	      cerr << "Directed edge! " << id << "->" << node[i].id << " no object. " << node[i].id << endl; 
#endif
	      continue;
	    }
	    NGT::GraphNode &nnode = *nn;
	    bool found = false;
	    for (size_t i = 0; i < nnode.size(); i++) {  
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	      if (nnode.at(i, graph.allocator).id == id) {
#else
	      if (nnode[i].id == id) {
#endif
		found = true;
		break;
	      }
	    }
	    if (!found) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	      cerr << "Directed edge! " << id << "->" << node.at(i, graph.allocator).id << " no edge. " 
		   << node.at(i, graph.allocator).id << "->" << id << endl; 
#else
	      cerr << "Directed edge! " << id << "->" << node[i].id << " no edge. " << node[i].id << "->" << id << endl; 
#endif
	      count++;
	    }
	  }
	}
	cerr << "# of not undirected edges=" << count << endl;
      }

      // calculate outdegree distance 10
      size_t d10count = 0;
      long double distance10 = 0.0;
      size_t d10SkipCount = 0;
      const size_t dcsize = 10;
      for (size_t id = 1; id < graph.size(); id++) {
	if (repo[id] == 0) {
	  continue;
	}
	NGT::GraphNode *n = 0;
	try {
	  n = outGraph.getNode(id);
	} catch(NGT::Exception &err) {
	  cerr << "ngt info: Warning. Cannot get the node. ID=" << id << ":" << err.what() << endl;
	  continue;
	}
	NGT::GraphNode &node = *n;
	if (node.size() < dcsize - 1) {
	  d10SkipCount++;
	  continue;
	}
	for (size_t i = 0; i < node.size(); i++) {  
	  if (i >= dcsize) {
	    break;
	  }
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  distance10 += node.at(i, graph.allocator).distance;
#else
	  distance10 += node[i].distance;
#endif
	  d10count++;
	}
      }
      distance10 /= (long double)d10count;

      // calculate indegree distance 10
      size_t ind10count = 0;
      long double indegreeDistance10 = 0.0;
      size_t ind10SkipCount = 0;
      for (size_t id = 1; id < indegree.size(); id++) {
	vector<float> &node = indegree[id];
	if (node.size() < dcsize - 1) {
	  ind10SkipCount++;
	  continue;
	}
	std::sort(node.begin(), node.end());
	for (size_t i = 0; i < node.size(); i++) {  
	  assert(i == 0 || node[i - 1] <= node[i]);
	  if (i >= dcsize) {
	    break;
	  }
	  indegreeDistance10 += node[i];
	  ind10count++;
	}
      }
      indegreeDistance10 /= (long double)ind10count;

      // calculate variance
      double averageNumberOfOutdegree = (double)numberOfOutdegree / (double)numberOfNodes;
      double sumOfSquareOfOutdegree = 0;
      double sumOfSquareOfIndegree = 0;
      for (size_t id = 1; id < graph.size(); id++) {
	if (repo[id] == 0) {
	  continue;
	}
	NGT::GraphNode *node = 0;
	try {
	  node = outGraph.getNode(id);
	} catch(NGT::Exception &err) {
	  cerr << "ngt info: Warning. Cannot get the node. ID=" << id << ":" << err.what() << endl;
	  continue;
	}
	size_t esize = node->size();
	sumOfSquareOfOutdegree += ((double)esize - averageNumberOfOutdegree) * ((double)esize - averageNumberOfOutdegree);
	sumOfSquareOfIndegree += ((double)indegreeCount[id] - averageNumberOfOutdegree) * ((double)indegreeCount[id] - averageNumberOfOutdegree);	
      }

      size_t numberOfNodesWithoutIndegree = 0;
      size_t maxNumberOfIndegree = 0;
      size_t minNumberOfIndegree = INT64_MAX;
      for (size_t id = 1; id < graph.size(); id++) {
	if (graph[id] == 0) {
	  continue;
	}
	if (indegreeCount[id] == 0) {
	  numberOfNodesWithoutIndegree++;
	  cerr << "Error! The node without incoming edges. " << id << endl;
	  valid = false;
	}
	if (indegreeCount[id] > static_cast<int>(maxNumberOfIndegree)) {
	  maxNumberOfIndegree = indegreeCount[id];
	}
	if (indegreeCount[id] < static_cast<int64_t>(minNumberOfIndegree)) {
	  minNumberOfIndegree = indegreeCount[id];
	}
	if (static_cast<int>(indegreeHistogram.size()) <= indegreeCount[id]) {
	  indegreeHistogram.resize(indegreeCount[id] + 1);
	}
	indegreeHistogram[indegreeCount[id]]++;
      }

      size_t count = 0;
      int medianOutdegree = -1;
      size_t modeOutdegree = 0;
      size_t max = 0;
      double c95 = 0.0;
      double c99 = 0.0;
      for (size_t i = 0; i < outdegreeHistogram.size(); i++) {
	count += outdegreeHistogram[i];
	if (medianOutdegree == -1 && count >= numberOfNodes / 2) {
	  medianOutdegree = i;
	}
	if (max < outdegreeHistogram[i]) {
	  max = outdegreeHistogram[i];
	  modeOutdegree = i;
	}
	if (count > numberOfNodes * 0.95) {
	  if (c95 == 0.0) {
	    c95 += i * (count - numberOfNodes * 0.95);
	  } else {
	    c95 += i * outdegreeHistogram[i];
	  }
	}
	if (count > numberOfNodes * 0.99) {
	  if (c99 == 0.0) {
	    c99 += i * (count - numberOfNodes * 0.99);
	  } else {
	    c99 += i * outdegreeHistogram[i];
	  }
	}
      }
      c95 /= (double)numberOfNodes * 0.05;
      c99 /= (double)numberOfNodes * 0.01;

      count = 0;
      int medianIndegree = -1;
      size_t modeIndegree = 0;
      max = 0;
      double c5 = 0.0;
      double c1 = 0.0;
      for (size_t i = 0; i < indegreeHistogram.size(); i++) {
	if (count < numberOfNodes * 0.05) {
	  if (count + indegreeHistogram[i] >= numberOfNodes * 0.05) {
	    c5 += i * (numberOfNodes * 0.05 - count);
	  } else {
	    c5 += i * indegreeHistogram[i];
	  }
	}
	if (count < numberOfNodes * 0.01) {
	  if (count + indegreeHistogram[i] >= numberOfNodes * 0.01) {
	    c1 += i * (numberOfNodes * 0.01 - count);
	  } else {
	    c1 += i * indegreeHistogram[i];
	  }
	}
	count += indegreeHistogram[i];
	if (medianIndegree == -1 && count >= numberOfNodes / 2) {
	  medianIndegree = i;
	}
	if (max < indegreeHistogram[i]) {
	  max = indegreeHistogram[i];
	  modeIndegree = i;
	}
      }
      c5 /= (double)numberOfNodes * 0.05;
      c1 /= (double)numberOfNodes * 0.01;

      cerr << "The size of object array=" << repo.size() << endl;
      cerr << "# of removed objects=" << removedObjectCount << "/" << repo.size() << endl;
      cerr << "# of nodes=" << numberOfNodes << endl;
      cerr << "# of edges=" << numberOfOutdegree << endl;
      cerr << "# of nodes without edges=" << numberOfNodesWithoutEdges << endl;
      cerr << "Max outdegree=" << maxNumberOfOutdegree << endl;
      cerr << "Min outdegree=" << minNumberOfOutdegree << endl;
      cerr << "Average number of edges=" << (double)numberOfOutdegree / (double)numberOfNodes << endl;
      cerr << "Average distance of edges=" << setprecision(10) << distance / (double)numberOfOutdegree << endl;
      cerr << "# of nodes where indegree is 0=" << numberOfNodesWithoutIndegree << endl;
      cerr << "Max indegree=" << maxNumberOfIndegree << endl;
      cerr << "Min indegree=" << minNumberOfIndegree << endl;
      cerr << "#-nodes,#-edges,#-no-indegree,avg-edges,avg-dist,max-out,min-out,v-out,max-in,min-in,v-in,med-out,"
	"med-in,mode-out,mode-in,c95,c5,o-distance(10),o-skip,i-distance(10),i-skip:" 
	   << numberOfNodes << ":" << numberOfOutdegree << ":" << numberOfNodesWithoutIndegree << ":" 
	   << setprecision(10) << (double)numberOfOutdegree / (double)numberOfNodes << ":"
	   << distance / (double)numberOfOutdegree << ":"
	   << maxNumberOfOutdegree << ":" << minNumberOfOutdegree << ":" << sumOfSquareOfOutdegree / (double)numberOfOutdegree<< ":"
	   << maxNumberOfIndegree << ":" << minNumberOfIndegree << ":" << sumOfSquareOfIndegree / (double)numberOfOutdegree << ":"
	   << medianOutdegree << ":" << medianIndegree << ":" << modeOutdegree << ":" << modeIndegree 
	   << ":" << c95 << ":" << c5 << ":" << c99 << ":" << c1 << ":" << distance10 << ":" << d10SkipCount << ":"
	   << indegreeDistance10 << ":" << ind10SkipCount << endl;
      if (mode == 'h') {
	cerr << "#\tout\tin" << endl;
	for (size_t i = 0; i < outdegreeHistogram.size() || i < indegreeHistogram.size(); i++) {
	  size_t out = outdegreeHistogram.size() <= i ? 0 : outdegreeHistogram[i];
	  size_t in = indegreeHistogram.size() <= i ? 0 : indegreeHistogram[i];
	  cerr << i << "\t" << out << "\t" << in << endl;
	}
      }
      return valid;
    }

    size_t getObjectRepositorySize() { return objectSpace->getRepository().size(); }

    size_t getSizeOfElement() { return objectSpace->getSizeOfElement(); }

    Object *allocateObject(const string &textLine, const string &sep) {
      return objectSpace->allocateNormalizedObject(textLine, sep);
    }
    Object *allocateObject(const vector<double> &obj) { 
      return objectSpace->allocateNormalizedObject(obj);
    }
    Object *allocateObject(const vector<float> &obj) { 
      return objectSpace->allocateNormalizedObject(obj);
    }
    Object *allocateObject(const float *obj, size_t size) { 
      return objectSpace->allocateNormalizedObject(obj, size);
    }

    void deleteObject(Object *po) {
      return objectSpace->deleteObject(po);
    }

    ObjectSpace &getObjectSpace() { return *objectSpace; }

    void setProperty(NGT::Property &prop) {
      GraphIndex::property.set(prop);
      NeighborhoodGraph::property.set(prop);
      assert(property.dimension != 0);
    }

    void getProperty(NGT::Property &prop) {
      GraphIndex::property.get(prop);
      NeighborhoodGraph::property.get(prop);
    }

    NeighborhoodGraph::Property &getGraphProperty() { return NeighborhoodGraph::property; }

    virtual size_t getSharedMemorySize(ostream &os, SharedMemoryAllocator::GetMemorySizeType t) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      size_t size = repository.getAllocator().getMemorySize(t);
#else
      size_t size = 0;
#endif
      os << "graph=" << size << endl;
      return size;
    }

    protected:

    template <class REPOSITORY> void getSeedsFromGraph(REPOSITORY &repo, ObjectDistances &seeds) {
      if (repo.size() != 0) {
	size_t seedSize = repo.size() - 1 < (size_t)NeighborhoodGraph::property.seedSize ? 
	  repo.size() - 1 : (size_t)NeighborhoodGraph::property.seedSize;
	if (NeighborhoodGraph::property.seedType == NeighborhoodGraph::SeedTypeRandomNodes ||
	    NeighborhoodGraph::property.seedType == NeighborhoodGraph::SeedTypeNone) {
	  getRandomSeeds(repo, seeds, seedSize);
	} else if (NeighborhoodGraph::property.seedType == NeighborhoodGraph::SeedTypeFixedNodes) {
	  // To check speed using fixed seeds.
	  for (size_t i = 1; i <= seedSize; i++) {
	    ObjectDistance obj(i, 0.0);
	    seeds.push_back(obj);
	  }
	} else if (NeighborhoodGraph::property.seedType == NeighborhoodGraph::SeedTypeFirstNode) {
	  ObjectDistance obj(1, 0.0);
	  seeds.push_back(obj);
	} else {
	  getRandomSeeds(repo, seeds, seedSize);
	}
      }
    }

    // GraphIndex
    virtual void search(NGT::SearchContainer &sc, ObjectDistances &seeds) {
      if (sc.size == 0) {
	while (!sc.workingResult.empty()) sc.workingResult.pop();
	return;
      }
      if (seeds.size() == 0) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR) || !defined(NGT_GRAPH_READ_ONLY_GRAPH)
	getSeedsFromGraph(repository, seeds);
#else
	if (readOnly) {
	  getSeedsFromGraph(searchRepository, seeds);
	} else {
	  getSeedsFromGraph(repository, seeds);
	}
#endif
      }
      NGT::SearchContainer so(sc);
      try {
	if (readOnly) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR) || !defined(NGT_GRAPH_READ_ONLY_GRAPH)
	  NeighborhoodGraph::search(so, seeds);
#else
	  (*searchUnupdatableGraph)(*this, so, seeds);
#endif
	} else {
	  NeighborhoodGraph::search(so, seeds);
	}
	sc.workingResult = std::move(so.workingResult);
	sc.distanceComputationCount = so.distanceComputationCount;
      } catch(Exception &err) {
	cerr << err.what() << endl;
	Exception e(err);
	throw e;
      }
    }

    Index::Property			property;

    bool readOnly;
#ifdef NGT_GRAPH_READ_ONLY_GRAPH
    void (*searchUnupdatableGraph)(NGT::NeighborhoodGraph&, NGT::SearchContainer&, NGT::ObjectDistances&);
#endif
  };

  class GraphAndTreeIndex : public GraphIndex, public DVPTree {
  public:

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    GraphAndTreeIndex(const string &allocator, bool rdOnly = false):GraphIndex(allocator, false) {
      initialize(allocator, 0);
    }
    GraphAndTreeIndex(const string &allocator, NGT::Property &prop);
    void initialize(const string &allocator, size_t sharedMemorySize) {
      DVPTree::objectSpace = GraphIndex::objectSpace;
      DVPTree::open(allocator + "/tre", sharedMemorySize);
    }
#else
    GraphAndTreeIndex(const string &database, bool rdOnly = false) : GraphIndex(database, rdOnly) {
      GraphAndTreeIndex::loadIndex(database, rdOnly);
    }

    GraphAndTreeIndex(NGT::Property &prop) : GraphIndex(prop) {
      DVPTree::objectSpace = GraphIndex::objectSpace;
    }
#endif
    virtual ~GraphAndTreeIndex() {}

    void create() {}

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    void alignObjects()
    {
    }
#else 
    void alignObjects()
    {
      NGT::ObjectSpace &space = getObjectSpace();
      NGT::ObjectRepository &repo = space.getRepository();
      Object **object = repo.getPtr();
      vector<bool> exist(repo.size(), false);
      vector<NGT::Node::ID> leafNodeIDs;
      DVPTree::getAllLeafNodeIDs(leafNodeIDs);
      size_t objectCount = 0;
      for (size_t i = 0; i < leafNodeIDs.size(); i++) {
	ObjectDistances objects;
	DVPTree::getObjectIDsFromLeaf(leafNodeIDs[i], objects);
	for (size_t j = 0; j < objects.size(); j++) {
	  exist[objects[j].id] = true;
	  objectCount++;
	}
      }
      multimap<uint32_t, uint32_t> notexist; 
      if (objectCount != repo.size()) {
        for (size_t id = 1; id < exist.size(); id++) {
	  if (!exist[id]) {
            DVPTree::SearchContainer tso(*object[id]);
            tso.mode = DVPTree::SearchContainer::SearchLeaf;
            tso.radius = 0.0;
            tso.size = 1;
      	    try {
	      DVPTree::search(tso);
      	    } catch (Exception &err) {
	      stringstream msg;
	      msg << "GraphAndTreeIndex::getSeeds: Cannot search for tree.:" << err.what();
	      NGTThrowException(msg);
      	    }
	    notexist.insert(pair<uint32_t, uint32_t>(tso.nodeID.getID(), id));
	    objectCount++;
	  }
	}
      }
      assert(objectCount == repo.size() - 1);

      objectCount = 1;
      vector<pair<uint32_t, uint32_t> > order;  
      for (size_t i = 0; i < leafNodeIDs.size(); i++) {
	ObjectDistances objects;
	DVPTree::getObjectIDsFromLeaf(leafNodeIDs[i], objects);
	for (size_t j = 0; j < objects.size(); j++) {
	  order.push_back(pair<uint32_t, uint32_t>(objects[j].id, objectCount));
	  objectCount++;
	}
	auto nei = notexist.equal_range(leafNodeIDs[i].getID());
	for (auto ii = nei.first; ii != nei.second; ++ii) {
	  order.push_back(pair<uint32_t, uint32_t>((*ii).second, objectCount));
	  objectCount++;
	}
      }
      assert(objectCount == repo.size());
      Object *tmp = space.allocateObject();
      unordered_set<uint32_t> uncopiedObjects;
      for (size_t i = 1; i < repo.size(); i++) {
	uncopiedObjects.insert(i);
      }
      size_t copycount = 0;
      while (!uncopiedObjects.empty()) {
	size_t startID = *uncopiedObjects.begin();
	if (startID == order[startID - 1].first) {
	  uncopiedObjects.erase(startID);
	  copycount++;
	  continue;
	}
	size_t id = startID;
	space.copy(*tmp, *object[id]);
	uncopiedObjects.erase(id);    
	do {
	  space.copy(*object[id], *object[order[id - 1].first]);
	  copycount++;
	  id = order[id - 1].first;
	  uncopiedObjects.erase(id);
	} while (order[id - 1].first != startID);
	space.copy(*object[id], *tmp);
	copycount++;
      }
      space.deleteObject(tmp);

      assert(copycount == repo.size() - 1);

      sort(order.begin(), order.end());
      uncopiedObjects.clear();
      for (size_t i = 1; i < repo.size(); i++) {
	uncopiedObjects.insert(i);
      }
      copycount = 0;
      Object *tmpPtr;
      while (!uncopiedObjects.empty()) {
	size_t startID = *uncopiedObjects.begin();
	if (startID == order[startID - 1].second) {
	  uncopiedObjects.erase(startID);
	  copycount++;
	  continue;
	}
	size_t id = startID;
	tmpPtr = object[id];
	uncopiedObjects.erase(id);    
	do {
	  object[id] = object[order[id - 1].second];
	  copycount++;
	  id = order[id - 1].second;
	  uncopiedObjects.erase(id);
	} while (order[id - 1].second != startID);
	object[id] = tmpPtr;
	copycount++;
      }
      assert(copycount == repo.size() - 1);
    }
#endif // NGT_SHARED_MEMORY_ALLOCATOR

    void load(const string &ifile) {
      GraphIndex::load(ifile);
      DVPTree::objectSpace = GraphIndex::objectSpace;
    }

    void saveIndex(const string &ofile) {
      GraphIndex::saveIndex(ofile);
#ifndef NGT_SHARED_MEMORY_ALLOCATOR
      string fname = ofile + "/tre";
      ofstream ost(fname);
      if (!ost.is_open()) {
	stringstream msg;
	msg << "saveIndex:: Cannot open. " << fname;
	NGTThrowException(msg);
      }
      DVPTree::serialize(ost);
#endif
    }

    void loadIndex(const string &ifile, bool readOnly) {
      DVPTree::objectSpace = GraphIndex::objectSpace;
      ifstream ist(ifile + "/tre");
      DVPTree::deserialize(ist);
#ifdef NGT_GRAPH_READ_ONLY_GRAPH
      if (readOnly) {
	if (property.objectAlignment == NGT::Index::Property::ObjectAlignmentTrue) {
	  alignObjects();
	}
	GraphIndex::NeighborhoodGraph::loadSearchGraph(ifile);
      }
#endif
    }
    
    void exportIndex(const string &ofile) {
      GraphIndex::exportIndex(ofile);
      ofstream ost(ofile + "/tre");
      DVPTree::serializeAsText(ost);
    }

    void importIndex(const string &ifile) {
      string fname = ifile + "/tre";
      ifstream ist(fname);
      if (!ist.is_open()) {
	stringstream msg;
	msg << "importIndex:: Cannot open. " << fname;
	NGTThrowException(msg);
      }
      DVPTree::deserializeAsText(ist);
      GraphIndex::importIndex(ifile);
    }

    void remove(const ObjectID id, bool force = false) {
      Object *obj = 0;
      try {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	obj = GraphIndex::objectSpace->allocateObject(*GraphIndex::objectSpace->getRepository().get(id));
#else
	obj = GraphIndex::objectSpace->getRepository().get(id);
#endif
      } catch (Exception &err) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	GraphIndex::objectSpace->deleteObject(obj);
#endif
	if (force) {
	  try {
	    DVPTree::removeNaively(id);
          } catch(...) {}
	  try {
	    GraphIndex::remove(id, force);
          } catch(...) {}
	  stringstream msg;
	  msg << err.what() << " Even though the object could not be found, the object could be removed from the tree and graph if it existed in them.";
	  NGTThrowException(msg);
        }
	throw err;
      }
      NGT::SearchContainer so(*obj);
      ObjectDistances results;
      so.setResults(&results);
      so.id = 0;
      so.size = 2;
      so.radius = 0.0;
      so.explorationCoefficient = 1.1;
      ObjectDistances	seeds;
      seeds.push_back(ObjectDistance(id, 0.0));
      GraphIndex::search(so, seeds);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      GraphIndex::objectSpace->deleteObject(obj);
#endif
      if (results.size() == 0) {
	NGTThrowException("No found the specified id");
      }
      if (results.size() == 1) {
	try {
	  DVPTree::remove(id);
	} catch(Exception &err) {
	  cerr << "remove:: cannot remove from tree. id=" << id << " " << err.what() << endl;
	  throw err;
	}
      } else {
	ObjectID replaceID = id == results[0].id ? results[1].id : results[0].id;
	try {
	  DVPTree::replace(id, replaceID);
	} catch(Exception &err) {
	}
      }
      GraphIndex::remove(id, force);
    }

    void searchForNNGInsertion(Object &po, ObjectDistances &result) {
      NGT::SearchContainer sc(po);
      sc.setResults(&result);
      sc.size = NeighborhoodGraph::property.edgeSizeForCreation;
      sc.radius = FLT_MAX;
      sc.explorationCoefficient = NeighborhoodGraph::property.insertionRadiusCoefficient;
      sc.useAllNodesInLeaf = true;
      try {
	GraphAndTreeIndex::search(sc);
      } catch(Exception &err) {
	throw err;
      }
      if (static_cast<int>(result.size()) < NeighborhoodGraph::property.edgeSizeForCreation && 
	  result.size() < repository.size()) {
	if (sc.edgeSize != 0) {
	  try {
	    GraphAndTreeIndex::search(sc);
	  } catch(Exception &err) {
	    throw err;
	  }
	}
      }
    }

    void insert(ObjectID id) {
      ObjectRepository &fr = GraphIndex::objectSpace->getRepository();
      if (fr[id] == 0) {
	cerr << "GraphAndTreeIndex::insert empty " << id << endl;
	return;
      }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      Object &po = *GraphIndex::objectSpace->allocateObject(*fr[id]);
#else
      Object &po = *fr[id];
#endif
      ObjectDistances rs;
      if (NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeANNG) {
	searchForNNGInsertion(po, rs);
      } else {
	searchForKNNGInsertion(po, id, rs);
      }

      GraphIndex::insertNode(id, rs);

      if (((rs.size() > 0) && (rs[0].distance != 0.0)) || rs.size() == 0) {
	DVPTree::InsertContainer tiobj(po, id);
	try {
	  DVPTree::insert(tiobj);
	} catch (Exception &err) {
	  cerr << "GraphAndTreeIndex::insert: Fatal error" << endl;
	  cerr << err.what() << endl;
	  return;
	}
      }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      GraphIndex::objectSpace->deleteObject(&po);
#endif
    }

    void createIndex(size_t threadNumber);

    void createIndex(const vector<pair<NGT::Object*, size_t> > &objects, vector<InsertionResult> &ids,
		     double range, size_t threadNumber);

    void createTreeIndex();

    // GraphAndTreeIndex
    void getSeedsFromTree(NGT::SearchContainer &sc, ObjectDistances &seeds) {
      DVPTree::SearchContainer tso(sc.object);
      tso.mode = DVPTree::SearchContainer::SearchLeaf;
      tso.radius = 0.0;
      tso.size = 1;
      tso.distanceComputationCount = 0;
      try {
	DVPTree::search(tso);
      } catch (Exception &err) {
	stringstream msg;
	msg << "GraphAndTreeIndex::getSeeds: Cannot search for tree.:" << err.what();
	NGTThrowException(msg);
      }

      try {
	DVPTree::getObjectIDsFromLeaf(tso.nodeID, seeds);
      } catch (Exception &err) {
	stringstream msg;
	msg << "GraphAndTreeIndex::getSeeds: Cannot get a leaf.:" << err.what();
	NGTThrowException(msg);
      }
      sc.distanceComputationCount += tso.distanceComputationCount;
      if (sc.useAllNodesInLeaf || NeighborhoodGraph::property.seedType == NeighborhoodGraph::SeedTypeAllLeafNodes) {
	return;
      }
      // if seedSize is zero, the result size of the query is used as seedSize.
      size_t seedSize = NeighborhoodGraph::property.seedSize == 0 ? sc.size : NeighborhoodGraph::property.seedSize;
      seedSize = seedSize > sc.size ? sc.size : seedSize;
      if (seeds.size() > seedSize) {
	srand(tso.nodeID.getID());
	// to accelerate thinning data.
	for (size_t i = seeds.size(); i > seedSize; i--) {
	  double random = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
	  size_t idx = floor(i * random);
	  seeds[idx] = seeds[i - 1];
	}
	seeds.resize(seedSize);
      } else if (seeds.size() < seedSize) {
	// A lack of the seeds is compansated by random seeds.
	//getRandomSeeds(seeds, seedSize);
      }
    }

    // GraphAndTreeIndex
    void search(NGT::SearchContainer &sc) {
      sc.distanceComputationCount = 0;
      ObjectDistances	seeds;
      getSeedsFromTree(sc, seeds);
      GraphIndex::search(sc, seeds);
    }

    size_t getSharedMemorySize(ostream &os, SharedMemoryAllocator::GetMemorySizeType t) {
      return GraphIndex::getSharedMemorySize(os, t) + DVPTree::getSharedMemorySize(os, t);
    }

    bool verify(vector<uint8_t> &status, bool info) {
      bool valid = GraphIndex::verify(status, info);
      if (!valid) {
	cerr << "The graph or object is invalid!" << endl;
      }
      valid = valid && DVPTree::verify(GraphIndex::objectSpace->getRepository().size(), status);
      if (!valid) {
	cerr << "The tree is invalid" << endl;
      }
      // status: tree|graph|object
      for (size_t id = 1; id < status.size(); id++) {
	if (status[id] != 0x00 && status[id] != 0x07) {
	  if (status[id] == 0x03) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	    NGT::Object *po = GraphIndex::objectSpace->allocateObject(*GraphIndex::getObjectRepository().get(id));
	    NGT::SearchContainer sc(*po);
#else
	    NGT::SearchContainer sc(*GraphIndex::getObjectRepository().get(id));
#endif
	    NGT::ObjectDistances objects;
	    sc.setResults(&objects);
	    sc.id = 0;
	    sc.radius = 0.0;
	    sc.explorationCoefficient = 1.1;
	    sc.edgeSize = 0;
	    ObjectDistances	seeds;
	    seeds.push_back(ObjectDistance(id, 0.0));
	    objects.clear();
	    GraphIndex::search(sc, seeds);
	    size_t n = 0;
	    bool registeredIdenticalObject = false;
	    for (; n < objects.size(); n++) {
	      if (objects[n].id != id && status[objects[n].id] == 0x07) {
		registeredIdenticalObject = true;
		break;
	      }
	    }
	    if (!registeredIdenticalObject) {
	      cerr << "Warning: not found the registered same objects. id=" << id << " size=" << objects.size() << endl;
	      sc.id = 0;
	      sc.radius = FLT_MAX;
	      sc.explorationCoefficient = 1.2;
	      sc.edgeSize = 0;
	      sc.size = objects.size() < 100 ? 100 : objects.size() * 2;
	      ObjectDistances	seeds;
	      seeds.push_back(ObjectDistance(id, 0.0));
	      objects.clear();
	      GraphIndex::search(sc, seeds);
	      registeredIdenticalObject = false;
	      for (n = 0; n < objects.size(); n++) {
		if (objects[n].distance != 0.0) break;
		if (objects[n].id != id && status[objects[n].id] == 0x07) {
		  registeredIdenticalObject = true;
		  cerr << "info: found by using mode accurate search." << objects[n].id << endl;
		  break;
		}
	      }
	    }
	    if (!registeredIdenticalObject) {
	      cerr << "Warning: not found by using more accurate search." << endl;
	      sc.id = 0;
	      sc.radius = 0.0;
	      sc.explorationCoefficient = 1.1;
	      sc.edgeSize = 0;
	      sc.size = SIZE_MAX;
	      objects.clear();
	      linearSearch(sc);
	      n = 0;
	      registeredIdenticalObject = false;
	      for (; n < objects.size(); n++) {
		if (objects[n].distance != 0.0) break;
		if (objects[n].id != id && status[objects[n].id] == 0x07) {
		  registeredIdenticalObject = true;
		  cerr << "info: found by using linear search. " << objects[n].id << endl;
		  break;
		}
	      }
	    }
	    if (registeredIdenticalObject) {
	      if (info) {
		cerr << "Info ID=" << id << ":" << static_cast<int>(status[id]) << endl;
		cerr << "  found the valid same objects. " << objects[n].id << endl;
	      }

	      GraphNode &node = *GraphIndex::getNode(id);
	      bool found = false;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	      for (auto i = node.begin(GraphIndex::repository.allocator); i != node.end(GraphIndex::repository.allocator); ++i) {
#else
	      for (auto i = node.begin(); i != node.end(); ++i) {
#endif
		if ((*i).id == objects[n].id) {
		  found = true;
		}
	      }
	      if (!found) {
		cerr << "Warning no directed edge from " << id << " to " << objects[n].id << endl;
	      }
	    } else {
	      cerr << "Warning: not found the valid same object even by using linear search." << endl;
	      cerr << "Error! ID=" << id << ":" << static_cast<int>(status[id]) << endl;
	      valid = false;
	    }
	  } else if (status[id] == 0x01) {
	    if (info) {
	      cerr << "Warning! ID=" << id << ":" << static_cast<int>(status[id]) << endl;
	      cerr << "  not inserted into the indexes" << endl;
	    }
	  } else {
	    cerr << "Error! ID=" << id << ":" << static_cast<int>(status[id]) << endl;
	    valid = false;
	  }
	}
      }
      return valid;
    }

  };

  class Property : public Index::Property, public NeighborhoodGraph::Property {
  public:
    void setDefault() {
      Index::Property::setDefault();
      NeighborhoodGraph::Property::setDefault();
    }
    void clear() {
      Index::Property::clear();
      NeighborhoodGraph::Property::clear();
    }
    void set(NGT::Property &p) {
      Index::Property::set(p);
      NeighborhoodGraph::Property::set(p);
    }
    void load(const string &file) {
      NGT::PropertySet prop;
      prop.load(file + "/prf");
      Index::Property::importProperty(prop);
      NeighborhoodGraph::Property::importProperty(prop);
    }
    void importProperty(const string &file) {
      NGT::PropertySet prop;
      prop.load(file + "/prf");
      Index::Property::importProperty(prop);
      NeighborhoodGraph::Property::importProperty(prop);
    }
  };

} // namespace NGT

inline void 
NGT::Index::open(const string &database, bool rdOnly) {
  Index* idx = 0;
  NGT::Property prop;
  prop.load(database);
  if (prop.indexType == NGT::Index::Property::GraphAndTree) {
    idx = new NGT::GraphAndTreeIndex(database, rdOnly);
  } else if (prop.indexType == NGT::Index::Property::Graph) {
    idx = new NGT::GraphIndex(database, rdOnly);
  } else {
    NGTThrowException("Index::Open: Not found IndexType in property file.");
  }
  if (idx == 0) {
    stringstream msg;
    msg << "Index::open: Cannot open. " << database;
    NGTThrowException(msg);
  }
  index = idx;
  path = database;
}

inline void 
  NGT::Index::createGraphAndTree(const string &database, NGT::Property &prop, const string &dataFile,
					   size_t dataSize) {
  if (prop.dimension == 0) {
    NGTThrowException("Index::createGraphAndTree. Dimension is not specified.");
  }
  prop.indexType = NGT::Index::Property::IndexType::GraphAndTree;
  Index *idx = 0;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  mkdir(database);
  idx = new NGT::GraphAndTreeIndex(database, prop);
#else
  idx = new NGT::GraphAndTreeIndex(prop);
#endif
  assert(idx != 0);
  try {
    loadAndCreateIndex(*idx, database, dataFile, prop.threadPoolSize, dataSize);
  } catch(Exception &err) {
    delete idx;
    throw err;
  }
  delete idx;
}
template<typename T>
size_t NGT::Index::append(vector<T> &object) 
{
  if (getObjectSpace().getRepository().size() == 0) {
    getObjectSpace().getRepository().initialize();
  }

  auto *o = getObjectSpace().getRepository().allocateNormalizedPersistentObject(object);
  getObjectSpace().getRepository().push_back(dynamic_cast<PersistentObject*>(o));
  size_t oid = getObjectSpace().getRepository().size() - 1;
  return oid;
}
template<typename T>
size_t NGT::Index::insert(vector<T> &object) 
{
  if (getObjectSpace().getRepository().size() == 0) {
    getObjectSpace().getRepository().initialize();
  }

  auto *o = getObjectSpace().getRepository().allocateNormalizedPersistentObject(object);
  size_t oid = getObjectSpace().getRepository().insert(dynamic_cast<PersistentObject*>(o));
  return oid;
}

inline void 
  NGT::Index::createGraph(const string &database, NGT::Property &prop, const string &dataFile, size_t dataSize) {
  if (prop.dimension == 0) {
    NGTThrowException("Index::createGraphAndTree. Dimension is not specified.");
  }
  prop.indexType = NGT::Index::Property::IndexType::Graph;
  Index *idx = 0;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  mkdir(database);
  idx = new NGT::GraphIndex(database, prop);
#else
  idx = new NGT::GraphIndex(prop);
#endif
  assert(idx != 0);
  try {
    loadAndCreateIndex(*idx, database, dataFile, prop.threadPoolSize, dataSize);
  } catch(Exception &err) {
    delete idx;
    throw err;
  }
  delete idx;
}

inline void 
NGT::Index::loadAndCreateIndex(Index &index, const string &database, const string &dataFile, size_t threadSize, size_t dataSize) {
  NGT::Timer timer;
  timer.start();
  if (dataFile.size() != 0) {
    index.load(dataFile, dataSize);
  } else {
    index.saveIndex(database);
    return;
  }
  timer.stop();
  cerr << "Data loading time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  if (index.getObjectRepositorySize() == 0) {
    NGTThrowException("Index::create: Data file is empty.");
  }
  cerr << "# of objects=" << index.getObjectRepositorySize() - 1 << endl;
  timer.reset();
  timer.start();
  index.createIndex(threadSize);
  timer.stop();
  index.saveIndex(database);
  cerr << "Index creation time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
}

inline void 
NGT::Index::append(const string &database, const string &dataFile, size_t threadSize, size_t dataSize) {
  NGT::Index	index(database);
  NGT::Timer	timer;
  timer.start();
  if (dataFile.size() != 0) {
    index.append(dataFile, dataSize);
  } else {
    NGTThrowException("Index::append: No data file.");
  }
  timer.stop();
  cerr << "Data loading time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << index.getObjectRepositorySize() - 1 << endl;
  timer.reset();
  timer.start();
  index.createIndex(threadSize);
  timer.stop();
  index.saveIndex(database);
  cerr << "Index creation time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  return;
}

inline void 
NGT::Index::append(const string &database, const float *data, size_t dataSize, size_t threadSize) {
  NGT::Index	index(database);
  NGT::Timer	timer;
  timer.start();
  if (data != 0 && dataSize != 0) {
    index.append(data, dataSize);
  } else {
    NGTThrowException("Index::append: No data.");
  }
  timer.stop();
  cerr << "Data loading time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << index.getObjectRepositorySize() - 1 << endl;
  timer.reset();
  timer.start();
  index.createIndex(threadSize);
  timer.stop();
  index.saveIndex(database);
  cerr << "Index creation time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  return;
}

inline void 
NGT::Index::remove(const string &database, vector<ObjectID> &objects, bool force) {
  NGT::Index	index(database);
  NGT::Timer	timer;
  timer.start();
  for (vector<ObjectID>::iterator i = objects.begin(); i != objects.end(); i++) {
    try {
      index.remove(*i, force);
    } catch (Exception &err) {
      cerr << "Warning: Cannot remove the node. ID=" << *i << " : " << err.what() << endl;
      continue;
    }
  }
  timer.stop();
  cerr << "Data removing time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << index.getObjectRepositorySize() - 1 << endl;
  index.saveIndex(database);
  return;
}

inline void 
NGT::Index::importIndex(const string &database, const string &file) {
  Index *idx = 0;
  NGT::Property property;
  property.importProperty(file);
  NGT::Timer	timer;
  timer.start();
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  property.databaseType = NGT::Index::Property::DatabaseType::MemoryMappedFile;
  mkdir(database);
#else
  property.databaseType = NGT::Index::Property::DatabaseType::Memory;
#endif
  if (property.indexType == NGT::Index::Property::IndexType::GraphAndTree) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    idx = new NGT::GraphAndTreeIndex(database, property);
#else
    idx = new NGT::GraphAndTreeIndex(property);
#endif
    assert(idx != 0);
  } else if (property.indexType == NGT::Index::Property::IndexType::Graph) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    idx = new NGT::GraphIndex(database, property);
#else
    idx = new NGT::GraphIndex(property);
#endif
    assert(idx != 0);
  } else {
    NGTThrowException("Index::Open: Not found IndexType in property file.");
  }
  idx->importIndex(file);
  timer.stop();
  cerr << "Data importing time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << idx->getObjectRepositorySize() - 1 << endl;
  idx->saveIndex(database);
  delete idx;
}

inline void 
NGT::Index::exportIndex(const string &database, const string &file) {
  NGT::Index	idx(database);
  NGT::Timer	timer;
  timer.start();
  idx.exportIndex(file);
  timer.stop();
  cerr << "Data exporting time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << idx.getObjectRepositorySize() - 1 << endl;
}

inline void 
NGT::GraphIndex::constructObjectSpace(NGT::Property &prop) {
  assert(prop.dimension != 0);
  switch (prop.objectType) {
  case NGT::ObjectSpace::ObjectType::Float :
    objectSpace = new ObjectSpaceRepository<float, double>(prop.dimension, typeid(float), prop.distanceType);
    break;
  case NGT::ObjectSpace::ObjectType::Uint8 :
    objectSpace = new ObjectSpaceRepository<unsigned char, int>(prop.dimension, typeid(uint8_t), prop.distanceType);
    break;
  default:
    stringstream msg;
    msg << "Invalid Object Type in the property. " << prop.objectType;
    NGTThrowException(msg);	
  }
  prop.prefetchOffset = objectSpace->setPrefetchOffset(prop.prefetchOffset);
#ifdef NGT_GRAPH_READ_ONLY_GRAPH
  searchUnupdatableGraph = NeighborhoodGraph::Search::getMethod(prop.distanceType, prop.objectType);
#endif
}

inline void 
NGT::Index::Property::set(NGT::Property &prop) {
  if (prop.dimension != -1) dimension = prop.dimension;
  if (prop.threadPoolSize != -1) threadPoolSize = prop.threadPoolSize;
  if (prop.objectType != ObjectSpace::ObjectTypeNone) objectType = prop.objectType;
  if (prop.distanceType != DistanceType::DistanceTypeNone) distanceType = prop.distanceType;
  if (prop.indexType != IndexTypeNone) indexType = prop.indexType;
  if (prop.databaseType != DatabaseTypeNone) databaseType = prop.databaseType;
  if (prop.objectAlignment != ObjectAlignmentNone) objectAlignment = prop.objectAlignment;
  if (prop.pathAdjustmentInterval != -1) pathAdjustmentInterval = prop.pathAdjustmentInterval;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  if (prop.graphSharedMemorySize != -1) graphSharedMemorySize = prop.graphSharedMemorySize;
  if (prop.treeSharedMemorySize != -1) treeSharedMemorySize = prop.treeSharedMemorySize;
  if (prop.objectSharedMemorySize != -1) objectSharedMemorySize = prop.objectSharedMemorySize;
#endif
  if (prop.prefetchOffset != -1) prefetchOffset = prop.prefetchOffset;
}

inline void 
NGT::Index::Property::get(NGT::Property &prop) {
  prop.dimension = dimension;
  prop.threadPoolSize = threadPoolSize;
  prop.objectType = objectType;
  prop.distanceType = distanceType;
  prop.indexType = indexType;
  prop.databaseType = databaseType;
  prop.pathAdjustmentInterval = pathAdjustmentInterval;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  prop.graphSharedMemorySize = graphSharedMemorySize;
  prop.treeSharedMemorySize = treeSharedMemorySize;
  prop.objectSharedMemorySize = objectSharedMemorySize;
#endif
  prop.prefetchOffset = prefetchOffset;
}

inline void 
  NGT::GraphIndex::loadIndex(const string &ifile, bool readOnly) {
  objectSpace->deserialize(ifile + "/obj");
#ifdef NGT_GRAPH_READ_ONLY_GRAPH
  if (readOnly && property.indexType == NGT::Index::Property::IndexType::Graph) {
    GraphIndex::NeighborhoodGraph::loadSearchGraph(ifile);
  } else {
    ifstream isg(ifile + "/grp");
    repository.deserialize(isg);
  }
#else
  ifstream isg(ifile + "/grp");
  repository.deserialize(isg);
#endif
}

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
inline
NGT::GraphIndex::GraphIndex(const string &allocator, bool rdonly):readOnly(rdonly) {
  NGT::Property prop;
  prop.load(allocator);
  if (prop.databaseType != NGT::Index::Property::DatabaseType::MemoryMappedFile) {
    NGTThrowException("GraphIndex: Cannot open. Not memory mapped file type.");
  }
  initialize(allocator, prop);
}

inline
NGT::GraphAndTreeIndex::GraphAndTreeIndex(const string &allocator, NGT::Property &prop):GraphIndex(allocator, prop) {
  initialize(allocator, prop.treeSharedMemorySize);
}

inline
void NGT::GraphIndex::initialize(const string &allocator, NGT::Property &prop) {
  constructObjectSpace(prop);
  repository.open(allocator + "/grp", prop.graphSharedMemorySize);
  objectSpace->open(allocator + "/obj", prop.objectSharedMemorySize);
  setProperty(prop);
}
#else // NGT_SHARED_MEMORY_ALLOCATOR
inline
NGT::GraphIndex::GraphIndex(const string &database, bool rdOnly):readOnly(rdOnly) {
  NGT::Property prop;
  prop.load(database);
  if (prop.databaseType != NGT::Index::Property::DatabaseType::Memory) {
    NGTThrowException("GraphIndex: Cannot open. Not memory type.");
  }
  assert(prop.dimension != 0);
  initialize(prop);
  loadIndex(database, readOnly);
}
#endif

