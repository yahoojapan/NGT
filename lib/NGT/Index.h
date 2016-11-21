//
// Copyright (C) 2015-2016 Yahoo Japan Corporation
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

#include	<sys/time.h>
#include	<sys/stat.h>

#include	"NGT/defines.h"
#include	"NGT/Common.h"
#include	"NGT/Tree.h"
#include	"NGT/Thread.h"
#include	"NGT/Graph.h"

#ifdef NGT_EXPERIMENTAL_GRAPH
#include	"NGT/ExperimentalGraph.h"
#endif

namespace NGT {

  class Property;

  class Index {
  public:

    class Property {
    public:
      typedef ObjectSpace::DistanceType		DistanceType;
      typedef NeighborhoodGraph::SeedType	SeedType;
      typedef NeighborhoodGraph::GraphType	GraphType;
      enum IndexType {
	IndexTypeNone	= 0,
	GraphAndTree	= 1,
	Graph		= 2
      };
      enum ObjectType {
	ObjectTypeNone	= 0,
	Uint8		= 1,
	Float		= 2
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
	objectType	= ObjectType::Uint8;
	distanceType	= DistanceType::DistanceTypeL2;
	indexType	= IndexType::GraphAndTree;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	databaseType	= DatabaseType::MemoryMappedFile;
      	graphSharedMemorySize	= 512; // MB
      	treeSharedMemorySize	= 512; // MB
      	objectSharedMemorySize	= 512; // MB  512 is up to 20M objects.
#else
	databaseType	= DatabaseType::Memory;
#endif
      }
      void setNotAvailable() {
	dimension 	= -1;
	threadPoolSize	= -1;
	objectType	= ObjectTypeNone;
	distanceType	= DistanceType::DistanceTypeNone;
	indexType	= IndexTypeNone;
	databaseType	= DatabaseTypeNone;

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      	graphSharedMemorySize	= -1;
      	treeSharedMemorySize	= -1;
      	objectSharedMemorySize	= -1;
#endif
      }

      void exportProperty(NGT::PropertySet &p) {
	p.set("Dimension", dimension);
	p.set("ThreadPoolSize", threadPoolSize);
	switch (objectType) {
	case ObjectType::Uint8: p.set("ObjectType", "Integer-1"); break;
	case ObjectType::Float: p.set("ObjectType", "Float-4"); break;
	default : cerr << "Fatal error. Invalid object type. " << objectType <<endl; abort();
	}
	switch (distanceType) {
	case DistanceType::DistanceTypeNone:	p.set("DistanceType", "None"); break;
	case DistanceType::DistanceTypeL1:	p.set("DistanceType", "L1"); break;
	case DistanceType::DistanceTypeL2:	p.set("DistanceType", "L2"); break;
	case DistanceType::DistanceTypeHamming:	p.set("DistanceType", "Hamming"); break;
	case DistanceType::DistanceTypeAngle:	p.set("DistanceType", "Angle"); break;
	default : cerr << "Fatal error. Invalid distance type. " << distanceType <<endl; abort();
	}
	switch (indexType) {      
	case IndexType::GraphAndTree:	p.set("IndexType", "GraphAndTree"); break;
	case IndexType::Graph:		p.set("IndexType", "Graph"); break;
	default : cerr << "Fatal error. Invalid index type. " << indexType <<endl; abort();
	}
	switch (databaseType) {
	case DatabaseType::Memory:		p.set("DatabaseType", "Memory"); break;
	case DatabaseType::MemoryMappedFile:	p.set("DatabaseType", "MemoryMappedFile"); break;
	default : cerr << "Fatal error. Invalid database type. " << databaseType <<endl; abort();
	}
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	p.set("GraphSharedMemorySize", graphSharedMemorySize);
	p.set("TreeSharedMemorySize", treeSharedMemorySize);
	p.set("ObjectSharedMemorySize", objectSharedMemorySize);
#endif
      }

      void importProperty(NGT::PropertySet &p) {
	setDefault();
	dimension = p.getl("Dimension", dimension);
	threadPoolSize = p.getl("ThreadPoolSize", threadPoolSize);
	PropertySet::iterator it = p.find("ObjectType");
	if (it != p.end()) {
	  if (it->second == "Float-4") {
	    objectType = ObjectType::Float;
	  } else if (it->second == "Integer-1") {
	    objectType = ObjectType::Uint8;
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
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	graphSharedMemorySize  = p.getl("GraphSharedMemorySize", graphSharedMemorySize);
	treeSharedMemorySize   = p.getl("TreeSharedMemorySize", treeSharedMemorySize);
	objectSharedMemorySize = p.getl("ObjectSharedMemorySize", objectSharedMemorySize);
#endif
      }

      void set(NGT::Property &prop);
      void get(NGT::Property &prop);
      int		dimension;
      int		threadPoolSize;
      ObjectType	objectType;
      DistanceType	distanceType;
      IndexType		indexType;
      DatabaseType	databaseType;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      int		graphSharedMemorySize;
      int		treeSharedMemorySize;
      int		objectSharedMemorySize;
#endif
    };


    Index():index(0) {}
    Index(const string &database):index(0) { open(database); }
    Index(const string &database, NGT::Property &prop):index(0) { open(database, prop);  }
    virtual ~Index() { close(); }

    void open(const string &database, NGT::Property &prop) {
      open(database);
      setProperty(prop);	
    }
    void open(const string &database);

    void close() {
      if (index != 0) { 
	delete index;
	index = 0;
      } 
    }
    static void mkdir(const string &dir) { ::mkdir(dir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP |  S_IROTH | S_IXOTH); }
    static void createGraphAndTree(const string &database, NGT::Property &prop, const string &dataFile, size_t dataSize = 0);
    static void createGraph(const string &database, NGT::Property &prop, const string &dataFile, size_t dataSize = 0);
    static void append(const string &database, const string &dataFile, size_t threadSize, size_t dataSize);
    static void remove(const string &database, vector<ObjectID> &objects);
    static void exportIndex(const string &database, const string &file);
    static void importIndex(const string &database, const string &file);
    virtual void load(const string &ifile, size_t dataSize) { getIndex().load(ifile, dataSize); }
    virtual void append(const string &ifile, size_t dataSize) { getIndex().append(ifile, dataSize); }
    virtual size_t getObjectRepositorySize() { return getIndex().getObjectRepositorySize(); }
    virtual void createIndex(size_t threadNumber) { getIndex().createIndex(threadNumber); }
    virtual void saveIndex(const string &ofile) { getIndex().saveIndex(ofile); }
    virtual void loadIndex(const string &ofile) { getIndex().loadIndex(ofile); }
    virtual Object *allocateObject(const string &textLine, const string &sep) { return getIndex().allocateObject(textLine, sep); }
    virtual Object *allocateObject(vector<double> &obj) { return getIndex().allocateObject(obj); }
    virtual size_t getSizeOfElement() { return getIndex().getSizeOfElement(); }
    virtual void setProperty(NGT::Property &prop) { getIndex().setProperty(prop); }
    virtual void getProperty(NGT::Property &prop) { getIndex().getProperty(prop); }
    virtual void deleteObject(Object *po) { getIndex().deleteObject(po); }
    virtual void linearSearch(NGT::SearchContainer &sc) { getIndex().linearSearch(sc); }
    virtual void search(NGT::SearchContainer &sc) { getIndex().search(sc); }
    virtual void search(NGT::SearchContainer &sc, ObjectDistances &seeds) { getIndex().search(sc, seeds); };
    void searchUsingOnlyGraph(NGT::SearchContainer &sc) { ObjectDistances seeds; getIndex().search(sc, seeds); };
    virtual void remove(ObjectID id) { getIndex().remove(id); }
    virtual void exportIndex(const string &file) { getIndex().exportIndex(file); }
    virtual void importIndex(const string &file) { getIndex().importIndex(file); }
    virtual ObjectSpace &getObjectSpace() { return getIndex().getObjectSpace(); };
    Index &getIndex() {
      if (index == 0) {
	assert(index != 0);
	NGTThrowException("NGT::Index::getIndex: Index is unavailable.");	
      }
      return *index;
    }

  protected:
    static void loadAndCreateIndex(Index &index, const string &database, const string &dataFile,
				   size_t threadSize, size_t dataSize);

    Index *index;
  };

  class GraphIndex : public Index, 
#ifdef NGT_EXPERIMENTAL_GRAPH
    public ExperimentalNeighborhoodGraph {
#else
    public NeighborhoodGraph {
#endif
  public:

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    GraphIndex(const string &allocator);
    // When no index exists yet.
    GraphIndex(const string &allocator, NGT::Property &prop) {
      initialize(allocator, prop);
    }
    void initialize(const string &allocator, NGT::Property &prop);
#else // NGT_SHARED_MEMORY_ALLOCATOR
    GraphIndex(const string &database);
    // When no indexs exist yet,
    ////////////
    GraphIndex(NGT::Property &prop) {
      initialize(prop);
    }

    void initialize(NGT::Property &prop) {
      constructObjectSpace(prop);
      setProperty(prop);
    }

#endif // NGT_SHARED_MEMORY_ALLOCATOR

    virtual ~GraphIndex() {
      if (objectSpace != 0) {
	destructObjectSpace();
	objectSpace = 0;
      }
    }
    void constructObjectSpace(NGT::Property &prop);

    void destructObjectSpace() {
      if (objectSpace == 0) {
	return;
      }
      if (property.objectType == NGT::Index::Property::ObjectType::Float) {
	ObjectSpaceT<float, double> *os = (ObjectSpaceT<float, double>*)objectSpace;
#ifndef NGT_SHARED_MEMORY_ALLOCATOR
	os->deleteAll();
#endif
	delete os;
	os = 0;
      } else if (property.objectType == NGT::Index::Property::ObjectType::Uint8) {
	ObjectSpaceT<unsigned char, int> *os = (ObjectSpaceT<unsigned char, int>*)objectSpace;
#ifndef NGT_SHARED_MEMORY_ALLOCATOR
	os->deleteAll();
#endif
	delete os;
	os = 0;
      } else {
	cerr << "Cannot find Object Type in the property. " << property.objectType << endl;
      }
    }

    virtual void load(const string &ifile, size_t dataSize = 0) {
      ifstream is(ifile.c_str());
      objectSpace->readText(is, dataSize);
    }

    virtual void append(const string &ifile, size_t dataSize = 0) {
      ifstream is(ifile.c_str());
      objectSpace->appendText(is, dataSize);
    }

    virtual void saveIndex(const string &ofile) {
#ifndef NGT_SHARED_MEMORY_ALLOCATOR
      try {
	mkdir(ofile);
      } catch(...) {}
      objectSpace->serialize(ofile + "/obj");
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

    virtual void loadIndex(const string &ifile);

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

    virtual void search(NGT::SearchContainer &sc) {
      ObjectDistances seeds;
      search(sc, seeds);
    }

    // GraphIndex
    virtual void search(NGT::SearchContainer &sc, ObjectDistances &seeds) {
      if (seeds.size() == 0 && repository.size() != 0) {
	size_t seedSize = repository.size() - 1 < (size_t)NeighborhoodGraph::property.seedSize ? 
	  repository.size() - 1 : (size_t)NeighborhoodGraph::property.seedSize;
	if (NeighborhoodGraph::property.seedType == NeighborhoodGraph::SeedTypeRandomNodes) {
	  // get randomly nodes as seeds.
	  size_t repositorySize = repository.size();
	  while (seedSize != seeds.size()) {
	    double random = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
	    size_t idx = floor(repositorySize * random) + 1;
	    if (repository.isEmpty(idx)) {
	      continue;
	    }
	    ObjectDistance obj(idx, 0.0);
	    if (find(seeds.begin(), seeds.end(), obj) != seeds.end()) {
	      continue;
	    }
	    seeds.push_back(obj);
	  }
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
	  cerr << "NGT::GraphIndex::search: Fatal Error. Invalid seed type" << endl;
	  abort();
	}
      }

      ObjectDistances::iterator si;
      for (si = seeds.begin(); si != seeds.end(); si++) {
	(*si).distance = -1.0;
      }
      NGT::SearchContainer so(sc.object);
      ObjectDistances &rs = sc.getResult();

      rs.clear();
      so.setResults(&rs);
      so.id = 0;
      so.size = sc.size;
      so.radius = sc.radius;
      so.explorationCoefficient = sc.explorationCoefficient;
      try {
	NeighborhoodGraph::search(so, seeds);
      } catch(Exception &err) {
	cerr << err.what() << endl;
	Exception e(err);
	throw e;
      }
    }

    void remove(ObjectID id) {
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
    }

    void searchForKNNGInsertion(Object &po, ObjectID id, ObjectDistances &result) {
      double radius = FLT_MAX;
      size_t size = NeighborhoodGraph::property.edgeSizeForCreation + 1;
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

    size_t getObjectRepositorySize() { return objectSpace->getRepository().size(); }

    size_t getSizeOfElement() { return objectSpace->getSizeOfElement(); }

    Object *allocateObject(const string &textLine, const string &sep) {
      return objectSpace->allocateObject(textLine, sep);
    }

    Object *allocateObject(vector<double> &obj) {
      return objectSpace->allocateObject(obj);
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

    Index::Property			property;
  };

  class GraphAndTreeIndex : public GraphIndex, public DVPTree {
  public:

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    GraphAndTreeIndex(const string &allocator):GraphIndex(allocator) {
      initialize(allocator, 0);
    }
    GraphAndTreeIndex(const string &allocator, NGT::Property &prop);
    void initialize(const string &allocator, size_t sharedMemorySize) {
      DVPTree::objectSpace = GraphIndex::objectSpace;
      DVPTree::open(allocator + "/tre", sharedMemorySize);
    }
#else
    GraphAndTreeIndex(const string &database) : GraphIndex(database) {
      GraphAndTreeIndex::loadIndex(database);
    }
    // When no indexs exist yet,
    GraphAndTreeIndex(NGT::Property &prop) : GraphIndex(prop) {
      DVPTree::objectSpace = GraphIndex::objectSpace;
    }
#endif

    virtual ~GraphAndTreeIndex() {}

    void create() {}

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

    void loadIndex(const string &ifile) {
      DVPTree::objectSpace = GraphIndex::objectSpace;
      ifstream ist(ifile + "/tre");
      DVPTree::deserialize(ist);
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

    void remove(ObjectID id) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      Object &f = *GraphIndex::objectSpace->allocateObject(*GraphIndex::objectSpace->getRepository().get(id));
#else
      Object &f = *GraphIndex::objectSpace->getRepository().get(id);
#endif
      NGT::SearchContainer so(f);
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
      GraphIndex::objectSpace->deleteObject(&f);
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
	try {
	  DVPTree::replace(id, results[1].id);
	} catch(Exception &err) {
	  cerr << "remove:: warning: cannot replace from tree. id=" << id << " " << err.what() << endl;
	}
      }
      GraphIndex::remove(id);
    }

    void searchForNNGInsertion(Object &po, ObjectDistances &result) {
      NGT::SearchContainer sc(po);
      sc.setResults(&result);
      sc.size = NeighborhoodGraph::property.edgeSizeForCreation;
      sc.radius = FLT_MAX;
      sc.explorationCoefficient = NeighborhoodGraph::property.insertionRadiusCoefficient;
      try {
	GraphAndTreeIndex::search(sc);
      } catch(Exception &err) {
	throw err;
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


    void createTreeIndex();

    void getSeeds(Object &object, ObjectDistances &seeds) {
      DVPTree::SearchContainer tso(object);
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

      try {
	DVPTree::getObjectIDsFromLeaf(tso.nodeID, seeds);
      } catch (Exception &err) {
	stringstream msg;
	msg << "GraphAndTreeIndex::getSeeds: Cannot get a leaf.:" << err.what();
	NGTThrowException(msg);
      }
    }

    // GraphAndTreeIndex
    void search(NGT::SearchContainer &sc) {
      ObjectDistances	seeds;
      getSeeds(sc.object, seeds);
      GraphIndex::search(sc, seeds);
    }

  };

  class Property : public Index::Property, public NeighborhoodGraph::Property {
  public:
    void setDefault() {
      Index::Property::setDefault();
      NeighborhoodGraph::Property::setDefault();
    }
    void setNotAvailable() {
      Index::Property::setNotAvailable();
      NeighborhoodGraph::Property::setNotAvailable();
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
NGT::Index::open(const string &database) {
  Index* idx = 0;
  NGT::Property prop;
  prop.load(database);
  if (prop.indexType == NGT::Index::Property::GraphAndTree) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    idx = new NGT::GraphAndTreeIndex(database);
#else
    idx = new NGT::GraphAndTreeIndex(database);
#endif
  } else if (prop.indexType == NGT::Index::Property::Graph) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    idx = new NGT::GraphIndex(database);
#else
    idx = new NGT::GraphIndex(database);
#endif
  } else {
    NGTThrowException("Index::Open: Not found IndexType in property file.");
  }
  if (idx == 0) {
    stringstream msg;
    msg << "Index::open: Cannot open. " << database;
    NGTThrowException(msg);
  }
  index = idx;
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
    NGTThrowException("Index::create: No data file.");
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
  cerr << "Index creation time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  index.saveIndex(database);
}

inline void 
NGT::Index::append(const string &database, const string &dataFile, size_t threadSize, size_t dataSize) {
  NGT::Index	index(database);
  NGT::Timer	timer;
  timer.start();
  if (dataFile.size() != 0) {
    index.append(dataFile, dataSize);
  } else {
    NGTThrowException("Index::create: No data file.");
  }
  timer.stop();
  cerr << "Data loading time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << index.getObjectRepositorySize() - 1 << endl;
  timer.reset();
  timer.start();
  index.createIndex(threadSize);
  timer.stop();
  cerr << "Index creation time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  index.saveIndex(database);
  return;
}

inline void 
NGT::Index::remove(const string &database, vector<ObjectID> &objects) {
  NGT::Index	index(database);
  NGT::Timer	timer;
  timer.start();
  for (vector<ObjectID>::iterator i = objects.begin(); i != objects.end(); i++) {
    try {
      index.remove(*i);
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
  case NGT::Property::ObjectType::Float :
    objectSpace = new ObjectSpaceT<float, double>(prop.dimension, typeid(float), prop.distanceType);
    break;
  case NGT::Property::ObjectType::Uint8 :
    objectSpace = new ObjectSpaceT<unsigned char, int>(prop.dimension, typeid(uint8_t), prop.distanceType);
    break;
  default:
    cerr << "Invalid Object Type in the property. " << prop.objectType << endl;
  }
}

inline void 
NGT::Index::Property::set(NGT::Property &prop) {
  if (prop.dimension != -1) dimension = prop.dimension;
  if (prop.threadPoolSize != -1) threadPoolSize = prop.threadPoolSize;
  if (prop.objectType != ObjectTypeNone) objectType = prop.objectType;
  if (prop.distanceType != DistanceType::DistanceTypeNone) distanceType = prop.distanceType;
  if (prop.indexType != IndexTypeNone) indexType = prop.indexType;
  if (prop.databaseType != DatabaseTypeNone) databaseType = prop.databaseType;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  if (prop.graphSharedMemorySize != -1) graphSharedMemorySize = prop.graphSharedMemorySize;
  if (prop.treeSharedMemorySize != -1) treeSharedMemorySize = prop.treeSharedMemorySize;
  if (prop.objectSharedMemorySize != -1) objectSharedMemorySize = prop.objectSharedMemorySize;
#endif
}

inline void 
NGT::Index::Property::get(NGT::Property &prop) {
  prop.dimension = dimension;
  prop.threadPoolSize = threadPoolSize;
  prop.objectType = objectType;
  prop.distanceType = distanceType;
  prop.indexType = indexType;
  prop.databaseType = databaseType;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  prop.graphSharedMemorySize = graphSharedMemorySize;
  prop.treeSharedMemorySize = treeSharedMemorySize;
  prop.objectSharedMemorySize = objectSharedMemorySize;
#endif
}

inline void 
NGT::GraphIndex::loadIndex(const string &ifile) {
  NGT::Property prop;
  prop.load(ifile);
  constructObjectSpace(prop);
  objectSpace->deserialize(ifile + "/obj");
  ifstream isg(ifile + "/grp");
  repository.deserialize(isg);
}

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
inline
NGT::GraphIndex::GraphIndex(const string &allocator) {
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
#else
inline
NGT::GraphIndex::GraphIndex(const string &database) {
  NGT::Property prop;
  prop.load(database);
  if (prop.databaseType != NGT::Index::Property::DatabaseType::Memory) {
    NGTThrowException("GraphIndex: Cannot open. Not memory type.");
  }
  assert(prop.dimension != 0);
  initialize(prop);
  loadIndex(database);
}
#endif

