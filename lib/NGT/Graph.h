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

#include	<bitset>

#include	"NGT/defines.h"
#include	"NGT/Common.h"
#include	"NGT/ObjectSpaceRepository.h"



#ifdef NGT_GRAPH_CHECK_HASH_BASED_BOOLEAN_SET
#include	"NGT/HashBasedBooleanSet.h"
#endif

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
#define	NGT_SEED_SIZE				10
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

#ifdef NGT_GRAPH_READ_ONLY_GRAPH
    class ReadOnlyGraphNode : public vector<pair<uint32_t, PersistentObject*>> {
    public:
    };
    class SearchGraphRepository : public vector<ReadOnlyGraphNode> {
    public:
      SearchGraphRepository() {}
      bool isEmpty(size_t idx) { return (*this)[idx].empty(); }

      void deserialize(ifstream &is, ObjectRepository &objectRepository) {
	if (!is.is_open()) {
	  NGTThrowException("NGT::SearchGraph: Not open the specified stream yet.");
	}
	clear();
	size_t s;
	NGT::Serializer::read(is, s);
	resize(s);
	for (size_t id = 0; id < s; id++) {
	  char type;
	  NGT::Serializer::read(is, type);
	  switch(type) {
	  case '-':
	    break;
	  case '+':
	    {
	      ObjectDistances node;
	      node.deserialize(is);
	      ReadOnlyGraphNode &searchNode = at(id);
	      searchNode.reserve(node.size());
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	      for (auto ni = node.begin(); ni != node.end(); ni++) {
		cerr << "not implement" << endl;
		abort();
	      }
#else
	      for (auto ni = node.begin(); ni != node.end(); ni++) {
		searchNode.push_back(pair<uint32_t, Object*>((*ni).id, objectRepository.get((*ni).id)));
	      }
#endif
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

    };
#endif

    class NeighborhoodGraph {
    public:
      enum GraphType {
	GraphTypeNone	= 0,
	GraphTypeANNG	= 1,
	GraphTypeKNNG	= 2,
	GraphTypeBKNNG	= 3,
	GraphTypeONNG	= 4,
	GraphTypeIANNG	= 5,	// Improved ANNG
	GraphTypeDNNG	= 6
      };

      enum SeedType {
	SeedTypeNone		= 0,
	SeedTypeRandomNodes	= 1,
	SeedTypeFixedNodes	= 2,
	SeedTypeFirstNode	= 3,
	SeedTypeAllLeafNodes	= 4
      };

#ifdef NGT_GRAPH_READ_ONLY_GRAPH
      class Search {
      public:
	static void (*getMethod(NGT::ObjectSpace::DistanceType dtype, NGT::ObjectSpace::ObjectType otype))(NGT::NeighborhoodGraph&, NGT::SearchContainer&, NGT::ObjectDistances&)  {
	  switch (otype) {
	    default:
	    case NGT::ObjectSpace::Float:	    
	      switch (dtype) {
	      case NGT::ObjectSpace::DistanceTypeNormalizedCosine : return normalizedCosineSimilarityFloat;
	      case NGT::ObjectSpace::DistanceTypeCosine : 	    return cosineSimilarityFloat;
	      case NGT::ObjectSpace::DistanceTypeNormalizedAngle :  return normalizedAngleFloat;
	      case NGT::ObjectSpace::DistanceTypeAngle : 	    return angleFloat;
	      case NGT::ObjectSpace::DistanceTypeL2 : 		    return l2Float;
	      case NGT::ObjectSpace::DistanceTypeL1 : 		    return l1Float;
	      default:						    return l2Float;
	      }
	      break;
	    case NGT::ObjectSpace::Uint8:
	      switch (dtype) {
	      case NGT::ObjectSpace::DistanceTypeHamming : return hammingUint8;
	      case NGT::ObjectSpace::DistanceTypeL2 : 	   return l2Uint8;
	      case NGT::ObjectSpace::DistanceTypeL1 : 	   return l1Uint8;
	      default : 				   return l2Uint8;
	      }
	    break;
	  }
	  return l1Uint8;
	}
	static void l1Uint8(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds);
	static void l2Uint8(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds);
	static void l1Float(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds);
	static void l2Float(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds);
	static void hammingUint8(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds);
	static void cosineSimilarityFloat(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds);
	static void angleFloat(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds);
	static void normalizedCosineSimilarityFloat(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds);
	static void normalizedAngleFloat(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds);
      };
#endif

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
	  dynamicEdgeSizeBase		= 30;
	  buildTimeLimit		= 0.0;
	  outgoingEdge			= 10;
	  incomingEdge			= 80;
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
	  dynamicEdgeSizeBase		= -1;
	  buildTimeLimit		= -1;
	  outgoingEdge			= -1;
	  incomingEdge			= -1;
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
	  p.set("DynamicEdgeSizeBase", dynamicEdgeSizeBase);
	  p.set("BuildTimeLimit", buildTimeLimit);
	  p.set("OutgoingEdge", outgoingEdge);
	  p.set("IncomingEdge", incomingEdge);
	  switch (graphType) {
	  case NeighborhoodGraph::GraphTypeKNNG: p.set("GraphType", "KNNG"); break;
	  case NeighborhoodGraph::GraphTypeANNG: p.set("GraphType", "ANNG"); break;
	  case NeighborhoodGraph::GraphTypeBKNNG: p.set("GraphType", "BKNNG"); break;
	  case NeighborhoodGraph::GraphTypeONNG: p.set("GraphType", "ONNG"); break;
	  case NeighborhoodGraph::GraphTypeIANNG: p.set("GraphType", "IANNG"); break;
	  default: cerr << "Graph::exportProperty: Fatal error! Invalid Graph Type." << endl; abort();
	  }
	  switch (seedType) {
	  case NeighborhoodGraph::SeedTypeRandomNodes: p.set("SeedType", "RandomNodes"); break;
	  case NeighborhoodGraph::SeedTypeFixedNodes: p.set("SeedType", "FixedNodes"); break;
	  case NeighborhoodGraph::SeedTypeFirstNode: p.set("SeedType", "FirstNode"); break;
	  case NeighborhoodGraph::SeedTypeNone: p.set("SeedType", "None"); break;
	  case NeighborhoodGraph::SeedTypeAllLeafNodes: p.set("SeedType", "AllLeafNodes"); break;
	  default: cerr << "Graph::exportProperty: Fatal error! Invalid Seed Type." << endl; abort();
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
	  dynamicEdgeSizeBase = p.getl("DynamicEdgeSizeBase", dynamicEdgeSizeBase);
	  buildTimeLimit = p.getf("BuildTimeLimit", buildTimeLimit);
	  outgoingEdge = p.getl("OutgoingEdge", outgoingEdge);
	  incomingEdge = p.getl("IncomingEdge", incomingEdge);
	  PropertySet::iterator it = p.find("GraphType");
	  if (it != p.end()) {
	    if (it->second == "KNNG")		graphType = NeighborhoodGraph::GraphTypeKNNG;
	    else if (it->second == "ANNG")	graphType = NeighborhoodGraph::GraphTypeANNG;
	    else if (it->second == "BKNNG")     graphType = NeighborhoodGraph::GraphTypeBKNNG;
	    else if (it->second == "ONNG")      graphType = NeighborhoodGraph::GraphTypeONNG;
	    else if (it->second == "IANNG")	graphType = NeighborhoodGraph::GraphTypeIANNG;
	    else { cerr << "Graph::importProperty: Fatal error! Invalid Graph Type. " << it->second << endl; abort(); }
	  }
	  it = p.find("SeedType");
	  if (it != p.end()) {
	    if (it->second == "RandomNodes")		seedType = NeighborhoodGraph::SeedTypeRandomNodes;
	    else if (it->second == "FixedNodes")	seedType = NeighborhoodGraph::SeedTypeFixedNodes;
	    else if (it->second == "FirstNode")		seedType = NeighborhoodGraph::SeedTypeFirstNode;
	    else if (it->second == "None")		seedType = NeighborhoodGraph::SeedTypeNone;
	    else if (it->second == "AllLeafNodes")	seedType = NeighborhoodGraph::SeedTypeAllLeafNodes;
	    else { cerr << "Graph::importProperty: Fatal error! Invalid Seed Type. " << it->second << endl; abort(); }
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
	  os << "dynamicEdgeSizeBase="		<< p.dynamicEdgeSizeBase << endl;
	  os << "outgoingEdge="			<< p.outgoingEdge << endl;
	  os << "incomingEdge="			<< p.incomingEdge << endl;
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
	int16_t		dynamicEdgeSizeBase;
	float		buildTimeLimit;
	int16_t		outgoingEdge;
	int16_t		incomingEdge;
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
	case GraphTypeIANNG:
	  insertIANNGNode(id, objects);	
	  break;
	case GraphTypeONNG:
	  insertONNGNode(id, objects);	
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

      void insertIANNGNode(ObjectID id, ObjectDistances &results) {
	repository.insert(id, results);
	for (ObjectDistances::iterator ri = results.begin(); ri != results.end(); ri++) {
	  assert(id != (*ri).id);
	  addEdgeDeletingExcessEdges((*ri).id, id, (*ri).distance);
	}
	return;
      }

      void insertONNGNode(ObjectID id, ObjectDistances &results) {
	if (property.truncationThreshold != 0) {
	  stringstream msg;
	  msg << "NGT::insertONNGNode: truncation should be disabled!" << endl;
	  NGTThrowException(msg);
	}
	int count = 0;
	for (ObjectDistances::iterator ri = results.begin(); ri != results.end(); ri++, count++) {
	  assert(id != (*ri).id);
	  if (count >= property.incomingEdge) {
	    break;
	  }
	  addEdge((*ri).id, id, (*ri).distance); 
	}
	if (static_cast<int>(results.size()) > property.outgoingEdge) {
	  results.resize(property.outgoingEdge);
	}
	repository.insert(id, results);
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

      // setup edgeSize
      inline size_t getEdgeSize(NGT::SearchContainer &sc) {
	size_t edgeSize = INT_MAX;
	if (sc.edgeSize < 0) {
	  if (sc.edgeSize == -2) {
	    edgeSize = property.dynamicEdgeSizeBase + pow(10, (sc.explorationCoefficient - 1.0) * 20.0);
	  } else {
	    edgeSize = property.edgeSizeForSearch == 0 ? INT_MAX : property.edgeSizeForSearch;
	  }
	} else {
	  edgeSize = sc.edgeSize == 0 ? INT_MAX : sc.edgeSize;
	}
	return edgeSize;
      }

      void search(NGT::SearchContainer &sc, ObjectDistances &seeds);

#ifdef NGT_GRAPH_READ_ONLY_GRAPH
      template <typename COMPARATOR> void searchReadOnlyGraph(NGT::SearchContainer &sc, ObjectDistances &seeds);
#endif

      void removeEdge(ObjectID fid, ObjectID rmid) {
	GraphNode &rs = *getNode(fid);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	for (GraphNode::iterator ri = rs.begin(repository.allocator); ri != rs.end(repository.allocator); ri++) {
	  if ((*ri).id == rmid) {
	    rs.erase(ri, repository.allocator);
	    break;
	  }
	}
#else
	for (GraphNode::iterator ri = rs.begin(); ri != rs.end(); ri++) {
	  if ((*ri).id == rmid) {
	    rs.erase(ri);
	    break;
	  }
	}
#endif
      }

      void removeEdge(GraphNode &node, ObjectDistance &edge) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	GraphNode::iterator ni = std::lower_bound(node.begin(repository.allocator), node.end(repository.allocator), edge);
	if (ni != node.end(repository.allocator) && (*ni).id == edge.id) {
	  node.erase(ni, repository.allocator);
#else
	GraphNode::iterator ni = std::lower_bound(node.begin(), node.end(), edge);
	if (ni != node.end() && (*ni).id == edge.id) {
	  node.erase(ni);
#endif
	  return;
	}
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	if (ni == node.end(repository.allocator)) {
#else
	if (ni == node.end()) {
#endif
	  stringstream msg;
	  msg << "NGT::removeEdge: Cannot found " << edge.id;
	  NGTThrowException(msg);
	} else {
	  stringstream msg;
	  msg << "NGT::removeEdge: Cannot found " << (*ni).id << ":" << edge.id;
	  NGTThrowException(msg);
	}
      }

      void
	removeNode(ObjectID id) {
	repository.remove(id);
      }

      class BooleanVector : public vector<bool> {
      public:
        inline BooleanVector(size_t s):vector<bool>(s, false) {}
	inline void insert(size_t i) { vector<bool>::operator[](i) = true; }
      };

#ifdef NGT_GRAPH_VECTOR_RESULT
      typedef ObjectDistances ResultSet;
#else
      typedef priority_queue<ObjectDistance, vector<ObjectDistance>, less<ObjectDistance> > ResultSet;
#endif

#if defined(NGT_GRAPH_CHECK_BOOLEANSET)
      typedef BooleanSet DistanceCheckedSet;
#elif defined(NGT_GRAPH_CHECK_VECTOR)
      typedef BooleanVector DistanceCheckedSet;
#elif defined(NGT_GRAPH_CHECK_HASH_BASED_BOOLEAN_SET)
      typedef HashBasedBooleanSet DistanceCheckedSet;
#else
      class DistanceCheckedSet : public unordered_set<ObjectID> {
      public:
	bool operator[](ObjectID id) { return find(id) != end(); }
      };
#endif

      class NodeWithPosition : public ObjectDistance {
       public:
        NodeWithPosition(uint32_t p = 0):position(p){}
        NodeWithPosition(ObjectDistance &o):ObjectDistance(o), position(0){}
	NodeWithPosition &operator=(const NodeWithPosition &n) {
	  ObjectDistance::operator=(static_cast<const ObjectDistance&>(n));
	  position = n.position;
	  assert(id != 0);
	  return *this;
	}
	uint32_t	position;
      };

#ifdef NGT_GRAPH_UNCHECK_STACK
      typedef stack<ObjectDistance> UncheckedSet;
#else
#ifdef NGT_GRAPH_BETTER_FIRST_RESTORE
      typedef priority_queue<NodeWithPosition, vector<NodeWithPosition>, greater<NodeWithPosition> > UncheckedSet;
#else
      typedef priority_queue<ObjectDistance, vector<ObjectDistance>, greater<ObjectDistance> > UncheckedSet;
#endif
#endif
      void setupDistances(NGT::SearchContainer &sc, ObjectDistances &seeds);
      void setupDistances(NGT::SearchContainer &sc, ObjectDistances &seeds, double (&comparator)(const void*, const void*, size_t));

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

      static double (*getComparator())(const void*, const void*, size_t);


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
      void addEdge(GraphNode &node, ObjectID addID, Distance addDistance, bool identityCheck = true) {
	ObjectDistance obj(addID, addDistance);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	GraphNode::iterator ni = std::lower_bound(node.begin(repository.allocator), node.end(repository.allocator), obj);
	if ((ni != node.end(repository.allocator)) && ((*ni).id == addID)) {
	  if (identityCheck) {
	    stringstream msg;
	    msg << "NGT::addEdge: already existed! " << (*ni).id << ":" << addID;
	    NGTThrowException(msg);
	  }
	  return;
	}
#else
	GraphNode::iterator ni = std::lower_bound(node.begin(), node.end(), obj);
	if ((ni != node.end()) && ((*ni).id == addID)) {
	  if (identityCheck) {
	    stringstream msg;
	    msg << "NGT::addEdge: already existed! " << (*ni).id << ":" << addID;
	    NGTThrowException(msg);
	  }
	  return;
	}
#endif
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	node.insert(ni, obj, repository.allocator);
#else
	node.insert(ni, obj);
#endif
      }

      // identityCheck is checking whether the same edge has already added to the node.
      // return whether truncation is needed that means the node has too many edges.
      bool addEdge(ObjectID target, ObjectID addID, Distance addDistance, bool identityCheck = true) {
	size_t minsize = 0;
	GraphNode &node = property.truncationThreshold == 0 ? *getNode(target) : *getNode(target, minsize);
	addEdge(node, addID, addDistance, identityCheck);
	if ((size_t)property.truncationThreshold != 0 && node.size() - minsize > 
	    (size_t)property.truncationThreshold) {
	  return true;
	}
	return false;
      }

      void addEdgeDeletingExcessEdges(ObjectID target, ObjectID addID, Distance addDistance, bool identityCheck = true) {
	GraphNode &node = *getNode(target);
	size_t kEdge = property.edgeSizeForCreation - 1;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	if (node.size() > kEdge && node.at(kEdge, repository.allocator).distance >= addDistance) {
	  GraphNode &linkedNode = *getNode(node.at(kEdge, repository.allocator).id);
	  ObjectDistance linkedNodeEdge(target, node.at(kEdge, repository.allocator).distance);
	  if ((linkedNode.size() > kEdge) && node.at(kEdge, repository.allocator).distance >= 
	    linkedNode.at(kEdge, repository.allocator).distance) {
#else
	if (node.size() > kEdge && node[kEdge].distance >= addDistance) {
	  GraphNode &linkedNode = *getNode(node[kEdge].id);
	  ObjectDistance linkedNodeEdge(target, node[kEdge].distance);
	  if ((linkedNode.size() > kEdge) && node[kEdge].distance >= linkedNode[kEdge].distance) {
#endif
	    try {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	      removeEdge(node, node.at(kEdge, repository.allocator));
#else
	      removeEdge(node, node[kEdge]);
#endif
	    } catch (Exception &exp) {
	      stringstream msg;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	      msg << "addEdge: Cannot remove. (a) " << target << "," << addID << "," << node.at(kEdge, repository.allocator).id << "," << node.at(kEdge, repository.allocator).distance;
#else
	      msg << "addEdge: Cannot remove. (a) " << target << "," << addID << "," << node[kEdge].id << "," << node[kEdge].distance;
#endif
	      msg << ":" << exp.what();
	      NGTThrowException(msg.str());
	    }
	    try {
	      removeEdge(linkedNode, linkedNodeEdge);
	    } catch (Exception &exp) {
	      stringstream msg;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	      msg << "addEdge: Cannot remove. (b) " << target << "," << addID << "," << node.at(kEdge, repository.allocator).id << "," << node.at(kEdge, repository.allocator).distance;
#else
	      msg << "addEdge: Cannot remove. (b) " << target << "," << addID << "," << node[kEdge].id << "," << node[kEdge].distance;
#endif
	      msg << ":" << exp.what();
	      NGTThrowException(msg.str());
	    }
	  }
	}
	addEdge(node, addID, addDistance, identityCheck);
      }


#ifdef NGT_GRAPH_READ_ONLY_GRAPH
      void loadSearchGraph(const string &database) {
	ifstream isg(database + "/grp");
	NeighborhoodGraph::searchRepository.deserialize(isg, NeighborhoodGraph::getObjectRepository());
      }
#endif

    public:

      GraphRepository	repository;
      ObjectSpace	*objectSpace;

#ifdef NGT_GRAPH_READ_ONLY_GRAPH
      SearchGraphRepository searchRepository;
#endif      

      NeighborhoodGraph::Property		property;

    }; // NeighborhoodGraph

  } // NGT

