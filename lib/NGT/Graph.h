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

#include	<bitset>

#include	"NGT/defines.h"
#include	"NGT/Common.h"
#include	"NGT/ObjectSpace.h"



#ifndef NGT_GRAPH_CHECK_VECTOR
#include	<unordered_set>
#endif

#ifdef NGT_GRAPH_UNCHECK_STACK
#include	<stack>
#endif

#ifndef NGT_EXPLORATION_COEFFICIENT
#define NGT_EXPLORATION_COEFFICIENT		1.1
#endif

#ifndef NGT_INSERTION_EXPLORATION_COEFFICIENT
#define NGT_INSERTION_EXPLORATION_COEFFICIENT	1.1
#endif

#ifndef NGT_TRUNCATION_THRESHOLD
#define	NGT_TRUNCATION_THRESHOLD		50
#endif

#ifndef NGT_SEED_SIZE
//#define	NGT_SEED_SIZE				10
#define	NGT_SEED_SIZE				50
#endif

#ifndef NGT_CREATION_EDGE_SIZE
#define NGT_CREATION_EDGE_SIZE			10
#endif

namespace NGT {
  class Property;

  typedef GraphNode	GRAPH_NODE;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  class GraphRepository: public PersistentRepository<GRAPH_NODE> {
#else
  class GraphRepository: public Repository<GRAPH_NODE> {
#endif

  public:
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    typedef PersistentRepository<GRAPH_NODE>	VECTOR;
#else
    typedef Repository<GRAPH_NODE>	VECTOR;

    GraphRepository() {
      prevsize = new vector<unsigned short>;
    }
    virtual ~GraphRepository() {
      deleteAll();
      if (prevsize != 0) {
	delete prevsize;
      }
    }
#endif

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    void open(const string &file, size_t sharedMemorySize) {
      SharedMemoryAllocator &allocator = VECTOR::getAllocator();
      off_t *entryTable = (off_t*)allocator.construct(file, sharedMemorySize);
      if (entryTable == 0) {
	entryTable = (off_t*)construct();
	allocator.setEntry(entryTable);
      }
      assert(entryTable != 0);
      this->initialize(entryTable);
    }

    void *construct() {
      SharedMemoryAllocator &allocator = VECTOR::getAllocator();
      off_t *entryTable = new(allocator) off_t[2];
      entryTable[0] = allocator.getOffset(PersistentRepository<GRAPH_NODE>::construct());
      entryTable[1] = allocator.getOffset(new(allocator) Vector<unsigned short>);
      return entryTable;
    }

    void initialize(void *e) {
      SharedMemoryAllocator &allocator = VECTOR::getAllocator();
      off_t *entryTable = (off_t*)e;
      array = (ARRAY*)allocator.getAddr(entryTable[0]);
      PersistentRepository<GRAPH_NODE>::initialize(allocator.getAddr(entryTable[0]));
      prevsize = (Vector<unsigned short>*)allocator.getAddr(entryTable[1]);
    }
#endif

    void insert(ObjectID id, ObjectDistances &objects) {
      GRAPH_NODE *r = allocate();
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      (*r).copy(objects, VECTOR::getAllocator());
#else
      *r = objects;
#endif
      try {
	put(id, r);
      } catch (Exception exp) {
	delete r;
	throw exp;
      }
      if (id >= prevsize->size()) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	prevsize->resize(id + 1, VECTOR::getAllocator(), 0);
#else
	prevsize->resize(id + 1, 0);
#endif
      } else {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	(*prevsize).at(id, VECTOR::getAllocator()) = 0;
#else
	(*prevsize)[id] = 0;
#endif
      }
      return;
    }

    inline GRAPH_NODE *get(ObjectID fid, size_t &minsize) {
      GRAPH_NODE *rs = VECTOR::get(fid);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      minsize = (*prevsize).at(fid, VECTOR::getAllocator());
#else
      minsize = (*prevsize)[fid];
#endif
      return rs;
    }
    void serialize(ofstream &os) {
      VECTOR::serialize(os);
      Serializer::write(os, *prevsize);
    }
    void deserialize(ifstream &is) {
      VECTOR::deserialize(is);      
      Serializer::read(is, *prevsize);
    }
    void show() {
      for (size_t i = 0; i < this->size(); i++) {
	cout << "Show graph " << i << " ";
	if ((*this)[i] == 0) {
	  cout << endl;
	  continue;
	}
	for (size_t j = 0; j < (*this)[i]->size(); j++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  cout << (*this)[i]->at(j, VECTOR::getAllocator()).id << ":" << (*this)[i]->at(j, VECTOR::getAllocator()).distance << " ";
#else
	  cout << (*this)[i]->at(j).id << ":" << (*this)[i]->at(j).distance << " ";
#endif
	}
	cout << endl;
      }
    }

    public:
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    Vector<unsigned short>	*prevsize;
#else
    vector<unsigned short>	*prevsize;
#endif
    };

    class NeighborhoodGraph {
    public:
      enum GraphType {
	GraphTypeNone	= 0,
	GraphTypeANNG	= 1,
	GraphTypeKNNG	= 2,
	GraphTypeBKNNG	= 3,
	GraphTypeDNNG	= 4
      };

      enum SeedType {
	SeedTypeNone		= 0,
	SeedTypeRandomNodes	= 1,
	SeedTypeFixedNodes	= 2,
	SeedTypeFirstNode	= 3,
	SeedTypeAllLeafNodes	= 4
      };

      class Property {
      public:
	Property() { setDefault(); }
	void setDefault() {
	  truncationThreshold		= 0;
	  edgeSizeForCreation 		= NGT_CREATION_EDGE_SIZE;
	  edgeSizeForSearch 		= 0;
	  edgeSizeLimitForCreation 	= 5;
	  insertionRadiusCoefficient 	= NGT_INSERTION_EXPLORATION_COEFFICIENT;
	  seedSize			= NGT_SEED_SIZE;
	  seedType			= SeedTypeNone;
	  truncationThreadPoolSize	= 8;
	  batchSizeForCreation		= 200;
	  graphType			= GraphTypeANNG;
	}
	void clear() {
	  truncationThreshold		= -1;
	  edgeSizeForCreation		= -1;
	  edgeSizeForSearch		= -1;
	  edgeSizeLimitForCreation	= -1;
	  insertionRadiusCoefficient	= -1;
	  seedSize			= -1;
	  seedType			= SeedTypeNone;
	  truncationThreadPoolSize	= -1;
	  batchSizeForCreation		= -1;
	  graphType			= GraphTypeNone;
	}
	void set(NGT::Property &prop);
	void get(NGT::Property &prop);

	void exportProperty(NGT::PropertySet &p) {
	  p.set("IncrimentalEdgeSizeLimitForTruncation", truncationThreshold);
	  p.set("EdgeSizeForCreation", edgeSizeForCreation);
	  p.set("EdgeSizeForSearch", edgeSizeForSearch);
	  p.set("EdgeSizeLimitForCreation", edgeSizeLimitForCreation);
	  assert(insertionRadiusCoefficient >= 1.0);
	  p.set("EpsilonForCreation", insertionRadiusCoefficient - 1.0);
	  p.set("BatchSizeForCreation", batchSizeForCreation);
	  p.set("SeedSize", seedSize);
	  p.set("TruncationThreadPoolSize", truncationThreadPoolSize);
	  switch (graphType) {
	  case NeighborhoodGraph::GraphTypeKNNG: p.set("GraphType", "KNNG"); break;
	  case NeighborhoodGraph::GraphTypeANNG: p.set("GraphType", "ANNG"); break;
	  case NeighborhoodGraph::GraphTypeBKNNG: p.set("GraphType", "BKNNG"); break;
	  default: cerr << "Invalid Graph Type." << endl; abort();
	  }
	  switch (seedType) {
	  case NeighborhoodGraph::SeedTypeRandomNodes: p.set("SeedType", "RandomNodes"); break;
	  case NeighborhoodGraph::SeedTypeFixedNodes: p.set("SeedType", "FixedNodes"); break;
	  case NeighborhoodGraph::SeedTypeFirstNode: p.set("SeedType", "FirstNode"); break;
	  case NeighborhoodGraph::SeedTypeNone: p.set("SeedType", "None"); break;
	  case NeighborhoodGraph::SeedTypeAllLeafNodes: p.set("SeedType", "AllLeafNodes"); break;
	  default: cerr << "Invalid Seed Type." << endl; abort();
	  }
	}
	void importProperty(NGT::PropertySet &p) {
	  setDefault();
	  truncationThreshold = p.getl("IncrimentalEdgeSizeLimitForTruncation", truncationThreshold);
	  edgeSizeForCreation = p.getl("EdgeSizeForCreation", edgeSizeForCreation);
	  edgeSizeForSearch = p.getl("EdgeSizeForSearch", edgeSizeForSearch);
	  edgeSizeLimitForCreation = p.getl("EdgeSizeLimitForCreation", edgeSizeLimitForCreation);
	  insertionRadiusCoefficient = p.getf("EpsilonForCreation", insertionRadiusCoefficient);
	  insertionRadiusCoefficient += 1.0;
	  batchSizeForCreation = p.getl("BatchSizeForCreation", batchSizeForCreation);
	  seedSize = p.getl("SeedSize", seedSize);
	  truncationThreadPoolSize = p.getl("TruncationThreadPoolSize", truncationThreadPoolSize);
	  PropertySet::iterator it = p.find("GraphType");
	  if (it != p.end()) {
	    if (it->second == "KNNG")		graphType = NeighborhoodGraph::GraphTypeKNNG;
	    else if (it->second == "ANNG")	graphType = NeighborhoodGraph::GraphTypeANNG;
	    else if (it->second == "BKNNG")     graphType = NeighborhoodGraph::GraphTypeBKNNG;
	    else { cerr << "Fatal error! Invalid Graph Type. " << it->second << endl; abort(); }
	  }
	  it = p.find("SeedType");
	  if (it != p.end()) {
	    if (it->second == "RandomNodes")		seedType = NeighborhoodGraph::SeedTypeRandomNodes;
	    else if (it->second == "FixedNodes")	seedType = NeighborhoodGraph::SeedTypeFixedNodes;
	    else if (it->second == "FirstNode")		seedType = NeighborhoodGraph::SeedTypeFirstNode;
	    else if (it->second == "None")		seedType = NeighborhoodGraph::SeedTypeNone;
	    else if (it->second == "AllLeafNodes")	seedType = NeighborhoodGraph::SeedTypeAllLeafNodes;
	    else { cerr << "Fatal error! Invalid Seed Type. " << it->second << endl; abort(); }
	  }
	}
	friend ostream & operator<<(ostream& os, const Property& p) {
	  os << "truncationThreshold="		<< p.truncationThreshold << endl;
	  os << "edgeSizeForCreation="		<< p.edgeSizeForCreation << endl;
	  os << "edgeSizeForSearch="		<< p.edgeSizeForSearch << endl;
	  os << "edgeSizeLimitForCreation="	<< p.edgeSizeLimitForCreation << endl;
	  os << "insertionRadiusCoefficient="	<< p.insertionRadiusCoefficient << endl;
	  os << "insertionRadiusCoefficient="	<< p.insertionRadiusCoefficient << endl;
	  os << "seedSize="			<< p.seedSize << endl;
	  os << "seedType="			<< p.seedType << endl;
	  os << "truncationThreadPoolSize="	<< p.truncationThreadPoolSize << endl;
	  os << "batchSizeForCreation="		<< p.batchSizeForCreation << endl;
	  os << "graphType="			<< p.graphType << endl;
	  return os;
	}

	int16_t		truncationThreshold;
	int16_t		edgeSizeForCreation;
	int16_t		edgeSizeForSearch;
	int16_t		edgeSizeLimitForCreation;
	double		insertionRadiusCoefficient;
	int16_t		seedSize;
	SeedType	seedType;
	int16_t		truncationThreadPoolSize;
	int16_t		batchSizeForCreation;
	GraphType	graphType;
      };

      NeighborhoodGraph(): objectSpace(0) {
	property.truncationThreshold = NGT_TRUNCATION_THRESHOLD;
	// initialize random to generate random seeds
#ifdef NGT_DISABLE_SRAND_FOR_RANDOM
	struct timeval randTime;
	gettimeofday(&randTime, 0);
	srand(randTime.tv_usec);
#endif
      }
      inline GraphNode *getNode(ObjectID fid, size_t &minsize) { return repository.get(fid, minsize); }
      inline GraphNode *getNode(ObjectID fid) { return repository.VECTOR::get(fid); }
      void insertNode(ObjectID id,  ObjectDistances &objects) {
	switch (property.graphType) {
	case GraphTypeANNG:
	  insertANNGNode(id, objects);	
	  break;
	case GraphTypeKNNG:
	  insertKNNGNode(id, objects);
	  break;
	case GraphTypeBKNNG:
	  insertBKNNGNode(id, objects);
	  break;
	case GraphTypeNone:
	  NGTThrowException("NGT::insertNode: GraphType is not specified.");
	  break;
	default:
	  NGTThrowException("NGT::insertNode: GraphType is invalid.");
	  break;
	}
      }

      void insertBKNNGNode(ObjectID id, ObjectDistances &results) {
	if (repository.isEmpty(id)) {
	  repository.insert(id, results);
	} else {
	  GraphNode &rs = *getNode(id);
	  for (ObjectDistances::iterator ri = results.begin(); ri != results.end(); ri++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    rs.push_back((*ri), repository.allocator);
#else
	    rs.push_back((*ri));
#endif
	  }
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  std::sort(rs.begin(repository.allocator), rs.end(repository.allocator));
	  ObjectID prev = 0;
	  for (GraphNode::iterator ri = rs.begin(repository.allocator); ri != rs.end(repository.allocator);) {
	    if (prev == (*ri).id) {
	      ri = rs.erase(ri, repository.allocator);
	      continue;
	    }
	    prev = (*ri).id;
	    ri++;
	  }
#else
	  std::sort(rs.begin(), rs.end());
	  ObjectID prev = 0;
	  for (GraphNode::iterator ri = rs.begin(); ri != rs.end();) {
	    if (prev == (*ri).id) {
	      ri = rs.erase(ri);
	      continue;
	    }
	    prev = (*ri).id;
	    ri++;
	  }
#endif
	}
	for (ObjectDistances::iterator ri = results.begin(); ri != results.end(); ri++) {
	  assert(id != (*ri).id);
	  addBKNNGEdge((*ri).id, id, (*ri).distance);
	}
	return;
      }

      void insertKNNGNode(ObjectID id, ObjectDistances &results) {
	repository.insert(id, results);
      }

      void insertANNGNode(ObjectID id, ObjectDistances &results) {
	repository.insert(id, results);
	queue<ObjectID> truncateQueue;
	for (ObjectDistances::iterator ri = results.begin(); ri != results.end(); ri++) {
	  assert(id != (*ri).id);
	  if (addEdge((*ri).id, id, (*ri).distance)) {
	    truncateQueue.push((*ri).id);
	  }
	}
	while (!truncateQueue.empty()) {
	  ObjectID tid = truncateQueue.front();
	  truncateEdges(tid);
	  truncateQueue.pop();
	}
	return;
      }

      void removeEdgesReliably(ObjectID id);

      int truncateEdgesOptimally(ObjectID id, GraphNode &results, size_t truncationSize);

      int truncateEdges(ObjectID id) {
	GraphNode &results = *getNode(id);
	if (results.size() == 0) {
	  return -1;
	}

	size_t truncationSize = NGT_TRUNCATION_THRESHOLD;
	if (truncationSize < (size_t)property.edgeSizeForCreation) {
	  truncationSize = property.edgeSizeForCreation;
	}
	return truncateEdgesOptimally(id, results, truncationSize);
      }

      void search(NGT::SearchContainer &sc, ObjectDistances &seeds);

      void removeEdge(ObjectID fid, ObjectID rmid) {
	// have not been tested yet.
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	GraphNode &rs = *getNode(fid);
	for (GraphNode::iterator ri = rs.begin(repository.allocator); ri != rs.end(repository.allocator); ri++) {
	  if ((*ri).id == rmid) {
	    rs.erase(ri, repository.allocator);
	    break;
	  }
	}
#else
	GraphNode &rs = *getNode(fid);
	for (GraphNode::iterator ri = rs.begin(); ri != rs.end(); ri++) {
	  if ((*ri).id == rmid) {
	    rs.erase(ri);
	    break;
	  }
	}
#endif
      }

      void
	removeNode(ObjectID id) {
	repository.remove(id);
      }
#ifdef NGT_GRAPH_VECTOR_RESULT
      typedef ObjectDistances ResultSet;
#else
      typedef priority_queue<ObjectDistance, vector<ObjectDistance>, less<ObjectDistance> > ResultSet;
#endif

#ifdef NGT_GRAPH_CHECK_VECTOR
#if defined(NGT_GRAPH_CHECK_BITSET)
      typedef bitset<10000001> DistanceCheckedSet;
#elif defined(NGT_GRAPH_CHECK_BOOLEANSET)
      typedef BooleanSet DistanceCheckedSet;
#else
      typedef vector<bool> DistanceCheckedSet;
#endif
#else // NGT_GRAPH_CHECK_VECTOR
      typedef unordered_set<ObjectID> DistanceCheckedSet;
#endif // NGT_GRAPH_CHECK_VECTOR

#ifdef NGT_GRAPH_UNCHECK_STACK
      typedef stack<ObjectDistance> UncheckedSet;
#else
      typedef priority_queue<ObjectDistance, vector<ObjectDistance>, greater<ObjectDistance> > UncheckedSet;
#endif

      void setupSeeds(SearchContainer &sc, ObjectDistances &seeds, ResultSet &results, 
		      UncheckedSet &unchecked, DistanceCheckedSet &distanceChecked);

      int getEdgeSize() {return property.edgeSizeForCreation;}

      ObjectRepository &getObjectRepository() { return objectSpace->getRepository(); }

      ObjectSpace &getObjectSpace() { return *objectSpace; }

      void deleteInMemory() {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	assert(0);
#else
	for (vector<NGT::GraphNode*>::iterator i = repository.begin(); i != repository.end(); i++) {
	  if ((*i) != 0) {
	    delete (*i);
	  }
	}
	repository.clear();
#endif
      }

    protected:
      void
	addBKNNGEdge(ObjectID target, ObjectID addID, Distance addDistance) {
	if (repository.isEmpty(target)) {
	  ObjectDistances objs;
	  objs.push_back(ObjectDistance(addID, addDistance));
	  repository.insert(target, objs);
	  return;
	}
	addEdge(target, addID, addDistance, false);
      }

    public:
      // identityCheck is checking whether the same edge has already added to the node.
      // return whether truncation is needed that means the node has too many edges.
      bool addEdge(ObjectID target, ObjectID addID, Distance addDistance,  bool identityCheck = true) {
	size_t minsize = 0;
	GraphNode &node = property.truncationThreshold == 0 ? *getNode(target) : *getNode(target, minsize);
	ObjectDistance obj(addID, addDistance);
	// this seach ocuppies about 1% of total insertion time.
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	GraphNode::iterator ni = std::lower_bound(node.begin(repository.allocator), node.end(repository.allocator), obj);
	if ((ni != node.end(repository.allocator)) && ((*ni).id == addID)) {
	  if (identityCheck) {
	    stringstream msg;
	    msg << "NGT::addEdge: already existed! " << (*ni).id << ":" << addID;
	    NGTThrowException(msg);
	  }
	  return false;
	}
#else
	GraphNode::iterator ni = std::lower_bound(node.begin(), node.end(), obj);
	if ((ni != node.end()) && ((*ni).id == addID)) {
	  if (identityCheck) {
	    stringstream msg;
	    msg << "NGT::addEdge: already existed! " << (*ni).id << ":" << addID;
	    NGTThrowException(msg);
	  }
	  return false;
	}
#endif
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	node.insert(ni, obj, repository.allocator);
#else
	node.insert(ni, obj);
#endif
	if ((size_t)property.truncationThreshold != 0 && node.size() - minsize > 
	    (size_t)property.truncationThreshold) {
	  return true;
	}
	return false;
      }

    public:

      GraphRepository	repository;
      ObjectSpace	*objectSpace;

      NeighborhoodGraph::Property		property;

    }; // NeighborhoodGraph

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

  } // NGT

