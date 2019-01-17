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

#include	"defines.h"
#include	"Graph.h"
#include	"Thread.h"
#include	"Index.h"


using namespace NGT;

void 
NeighborhoodGraph::Property::set(NGT::Property &prop) {
  if (prop.truncationThreshold != -1)		truncationThreshold = prop.truncationThreshold;
  if (prop.edgeSizeForCreation != -1)		edgeSizeForCreation = prop.edgeSizeForCreation;
  if (prop.edgeSizeForSearch != -1)		edgeSizeForSearch = prop.edgeSizeForSearch;
  if (prop.edgeSizeLimitForCreation != -1)	edgeSizeLimitForCreation = prop.edgeSizeLimitForCreation;
  if (prop.insertionRadiusCoefficient != -1)	insertionRadiusCoefficient = prop.insertionRadiusCoefficient;
  if (prop.seedSize != -1)			seedSize = prop.seedSize;
  if (prop.seedType != SeedTypeNone)		seedType = prop.seedType;
  if (prop.truncationThreadPoolSize != -1)	truncationThreadPoolSize = prop.truncationThreadPoolSize;
  if (prop.batchSizeForCreation != -1)		batchSizeForCreation = prop.batchSizeForCreation;
  if (prop.dynamicEdgeSizeBase != -1)           dynamicEdgeSizeBase = prop.dynamicEdgeSizeBase;
  if (prop.buildTimeLimit != -1)                buildTimeLimit = prop.buildTimeLimit;
  if (prop.graphType != GraphTypeNone)		graphType = prop.graphType;
}

void 
NeighborhoodGraph::Property::get(NGT::Property &prop) {
  prop.truncationThreshold 		= truncationThreshold;
  prop.edgeSizeForCreation		= edgeSizeForCreation;
  prop.edgeSizeForSearch		= edgeSizeForSearch;
  prop.edgeSizeLimitForCreation		= edgeSizeLimitForCreation;
  prop.insertionRadiusCoefficient	= insertionRadiusCoefficient;
  prop.seedSize				= seedSize;
  prop.seedType				= seedType;
  prop.truncationThreadPoolSize		= truncationThreadPoolSize;
  prop.batchSizeForCreation		= batchSizeForCreation;
  prop.dynamicEdgeSizeBase              = dynamicEdgeSizeBase;
  prop.graphType			= graphType;
  prop.buildTimeLimit                   = buildTimeLimit;
}


#ifdef NGT_GRAPH_READ_ONLY_GRAPH
void 
NeighborhoodGraph::Search::normalizedCosineSimilarityFloat(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds)
{
  graph.searchReadOnlyGraph<PrimitiveComparator::NormalizedCosineSimilarityFloat>(sc, seeds);
}

void 
NeighborhoodGraph::Search::cosineSimilarityFloat(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds)
{
  graph.searchReadOnlyGraph<PrimitiveComparator::CosineSimilarityFloat>(sc, seeds);
}

void 
NeighborhoodGraph::Search::normalizedAngleFloat(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds)
{
  graph.searchReadOnlyGraph<PrimitiveComparator::NormalizedAngleFloat>(sc, seeds);
}

void 
NeighborhoodGraph::Search::angleFloat(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds)
{
  graph.searchReadOnlyGraph<PrimitiveComparator::AngleFloat>(sc, seeds);
}

void 
NeighborhoodGraph::Search::l1Float(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds)
{
  graph.searchReadOnlyGraph<PrimitiveComparator::L1Float>(sc, seeds);
}

void 
NeighborhoodGraph::Search::l2Float(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds)
{
  graph.searchReadOnlyGraph<PrimitiveComparator::L2Float>(sc, seeds);
}

void 
NeighborhoodGraph::Search::l1Uint8(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds)
{
  graph.searchReadOnlyGraph<PrimitiveComparator::L1Uint8>(sc, seeds);
}

void 
NeighborhoodGraph::Search::l2Uint8(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds)
{
  graph.searchReadOnlyGraph<PrimitiveComparator::L2Uint8>(sc, seeds);
}

void 
NeighborhoodGraph::Search::hammingUint8(NeighborhoodGraph &graph, NGT::SearchContainer &sc, ObjectDistances &seeds)
{
  graph.searchReadOnlyGraph<PrimitiveComparator::HammingUint8>(sc, seeds);
}
#endif

void 
NeighborhoodGraph::setupSeeds(NGT::SearchContainer &sc, ObjectDistances &seeds, ResultSet &results, 
			      UncheckedSet &unchecked, DistanceCheckedSet &distanceChecked)
{
  ObjectRepository &objectRepository = getObjectRepository();
  NGT::ObjectSpace::Comparator &comparator = objectSpace->getComparator();
  size_t byteSizeOfObject = objectSpace->getByteSizeOfObject();
  ObjectDistances tmp;
  tmp.reserve(seeds.size());
  size_t seedSize = seeds.size();
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
#ifndef NGT_PREFETCH_DISABLED
  const size_t prefetchPos = MemoryCache::getPrefetchPos(seedSize);
  for (size_t i = 0; i < prefetchPos; i++) {
    MemoryCache::prefetch(reinterpret_cast<unsigned char *>(objectRepository.get(seeds[i].id)), byteSizeOfObject);
  }
#endif
  for (size_t i = 0; i < seedSize; i++) {
#ifndef NGT_PREFETCH_DISABLED
    if (i + prefetchPos < seedSize) {
      MemoryCache::prefetch(reinterpret_cast<unsigned char*>(objectRepository.get(seeds[i + prefetchPos].id)), byteSizeOfObject);
    }
#endif
    ObjectDistance result;
    result.id = seeds[i].id;;
    if (objectRepository.isEmpty(result.id)) {
      cerr << "setupseeds:fatal error " << result.id << endl;
      assert(0);
    }
    Distance d = comparator(sc.object, *objectRepository.get(result.id));
    result.distance = d;
    tmp.push_back(result);
  }
#else
  PersistentObject **objects = objectRepository.getPtr();
#ifndef NGT_PREFETCH_DISABLED
  const size_t prefetchPos = MemoryCache::getPrefetchPos(seedSize);
  for (size_t i = 0; i < prefetchPos; i++) {
    MemoryCache::prefetch(reinterpret_cast<unsigned char *>(objects[seeds[i].id]), byteSizeOfObject);
  }
#endif
  for (size_t i = 0; i < seedSize; i++) {
#ifndef NGT_PREFETCH_DISABLED
    if (i + prefetchPos < seedSize) {
      MemoryCache::prefetch(reinterpret_cast<unsigned char*>(objects[seeds[i + prefetchPos].id]), byteSizeOfObject);
    }
#endif
    ObjectDistance result;
    result.id = seeds[i].id;
    if (objectRepository.isEmpty(result.id)) {
      cerr << "setupseeds:fatal error " << result.id << endl;
      assert(0);
    }
    result.distance = comparator(sc.object, *objects[result.id]);
    tmp.push_back(result);
  }
#endif

#ifdef NGT_DISTANCE_COMPUTATION_COUNT
  sc.distanceComputationCount += seeds.size();
#endif

  std::sort(tmp.begin(), tmp.end());

  for (ObjectDistances::iterator ri = tmp.begin(); ri != tmp.end(); ri++) {
    if ((results.size() < (unsigned int)sc.size) && ((*ri).distance <= sc.radius)) {
      results.push((*ri));
    } else {
      break;
    }
  }

    if (results.size() >= sc.size) {
      sc.radius = results.top().distance;
    }

  for (ObjectDistances::iterator ri = tmp.begin(); ri != tmp.end(); ri++) {
#if !defined(NGT_GRAPH_CHECK_VECTOR) || defined(NGT_GRAPH_CHECK_BOOLEANSET)
    distanceChecked.insert((*ri).id);
#else
    distanceChecked[(*ri).id] = 1;
#endif
    unchecked.push(*ri);
  }
}

#ifdef NGT_GRAPH_READ_ONLY_GRAPH
  template <typename COMPARATOR>
  void
    NeighborhoodGraph::searchReadOnlyGraph(NGT::SearchContainer &sc, ObjectDistances &seeds)
  {
    if (sc.explorationCoefficient == 0.0) {
      sc.explorationCoefficient = NGT_EXPLORATION_COEFFICIENT;
    }

    // setup edgeSize
    size_t edgeSize = getEdgeSize(sc);

    UncheckedSet unchecked;

    DistanceCheckedSet distanceChecked(searchRepository.size());

    ResultSet results;

    setupSeeds(sc, seeds, results, unchecked, distanceChecked);

    Distance explorationRadius = sc.explorationCoefficient * sc.radius;
    const size_t dimension = objectSpace->getPaddedDimension();
    const size_t byteSizeOfObject = objectSpace->getByteSizeOfObject();
    ReadOnlyGraphNode *nodes = &searchRepository.front();
    ReadOnlyGraphNode *neighbors = 0;
    ObjectDistance result;
    ObjectDistance target;
    pair<uint32_t, PersistentObject*> *neighborptr;
    pair<uint32_t, PersistentObject*> *neighborendptr;
    while (!unchecked.empty()) {
      target = unchecked.top();
      unchecked.pop();
      if (target.distance > explorationRadius) {
	break;
      }
      neighbors = &nodes[target.id];
      neighborptr = &(*neighbors)[0];
      neighborendptr = neighborptr;
      size_t neighborSize = neighbors->size();
      neighborendptr += (neighborSize < edgeSize ? neighborSize : edgeSize);

      const size_t prefetchPos = MemoryCache::getPrefetchPos(neighborSize);
      for (size_t i = 0; i < prefetchPos; i++) {
	if (!distanceChecked[(*(neighborptr + i)).first]) {
	  unsigned char *ptr = reinterpret_cast<unsigned char*>((*(neighborptr + i)).second);
	  MemoryCache::prefetch(ptr, byteSizeOfObject);
	}
      }
      for (; neighborptr < neighborendptr; ++neighborptr) {
	if ((neighborptr + prefetchPos < neighborendptr) && !distanceChecked[(*(neighborptr + prefetchPos)).first]) {
	  unsigned char *ptr = reinterpret_cast<unsigned char*>((*(neighborptr + prefetchPos)).second);
	  MemoryCache::prefetch(ptr, byteSizeOfObject);
	}
#ifdef NGT_VISIT_COUNT
	sc.visitCount++;
#endif
	auto &neighbor = *neighborptr;
	if (distanceChecked[neighbor.first]) {
	  continue;
	}
        distanceChecked.insert(neighbor.first);

	Distance distance = COMPARATOR::compare((void*)&sc.object[0], 
						(void*)&(*static_cast<PersistentObject*>(neighbor.second))[0], dimension);
	if (distance <= explorationRadius) {
	  result.set(neighbor.first, distance);
	  unchecked.push(result);
	  if (distance <= sc.radius) {
	    results.push(result);
	    if (results.size() >= sc.size) {
	      if (results.size() > sc.size) {
	        results.pop();
	      }
	      sc.radius = results.top().distance;
	      explorationRadius = sc.explorationCoefficient * sc.radius;
	    } 
	  } 
	} 
      } 
    } 

    if (sc.resultIsAvailable()) { 
      ObjectDistances &qresults = sc.getResult();
      qresults.clear();
      qresults.moveFrom(results);
    } else {
      sc.workingResult = std::move(results);
    }
  }

#endif

  void
    NeighborhoodGraph::search(NGT::SearchContainer &sc, ObjectDistances &seeds)
  {
    if (sc.explorationCoefficient == 0.0) {
      sc.explorationCoefficient = NGT_EXPLORATION_COEFFICIENT;
    }

    // setup edgeSize
    size_t edgeSize = getEdgeSize(sc);

    UncheckedSet unchecked;
#if defined(NGT_GRAPH_CHECK_BITSET)
    DistanceCheckedSet distanceChecked(0);
#elif defined(NGT_GRAPH_CHECK_BOOLEANSET)
    DistanceCheckedSet distanceChecked(repository.size());
#elif defined(NGT_GRAPH_CHECK_HASH_BASED_BOOLEAN_SET)
    DistanceCheckedSet distanceChecked(repository.size());
#elif defined(NGT_GRAPH_CHECK_VECTOR)
    DistanceCheckedSet distanceChecked(repository.size());
#else 
    DistanceCheckedSet distanceChecked;
#endif

    ResultSet results;
    setupSeeds(sc, seeds, results, unchecked, distanceChecked);
    Distance explorationRadius = sc.explorationCoefficient * sc.radius;
    NGT::ObjectSpace::Comparator &comparator = objectSpace->getComparator();
    ObjectRepository &objectRepository = getObjectRepository();
    const size_t byteSizeOfObject = objectSpace->getByteSizeOfObject();
    ObjectDistance result;
#ifdef NGT_GRAPH_BETTER_FIRST_RESTORE
    NodeWithPosition target;
#else
    ObjectDistance target;
#endif
    ObjectDistance *neighborptr;
    ObjectDistance *neighborendptr;
    while (!unchecked.empty()) {
      target = unchecked.top();
      unchecked.pop();
      if (target.distance > explorationRadius) {
	break;
      }
      GraphNode *neighbors = 0;
      try {
	neighbors = getNode(target.id);
      } catch(Exception &err) {
       cerr << "Graph::search: Warning. " << err.what() << "  ID=" << target.id << endl;
       continue;
      }
      if (neighbors->size() == 0) {
	continue;
      }
#ifdef NGT_GRAPH_BETTER_FIRST_RESTORE
      uint32_t position = target.position;
#endif
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
#ifdef NGT_GRAPH_BETTER_FIRST_RESTORE
      neighborptr = &(*neighbors).at(position, repository.allocator);
#else
      neighborptr = &(*neighbors).at(0, repository.allocator);
#endif
#else
#ifdef NGT_GRAPH_BETTER_FIRST_RESTORE
      neighborptr = &(*neighbors)[position];
#else
      neighborptr = &(*neighbors)[0];
#endif
#endif
      neighborendptr = neighborptr;
      if (edgeSize == 0) {
	neighborendptr += neighbors->size();
      } else {
	neighborendptr += (neighbors->size() < edgeSize ? neighbors->size() : edgeSize);
      }
#ifdef NGT_GRAPH_BETTER_FIRST_RESTORE
      neighborendptr -= position;
#endif
#ifdef NGT_GRAPH_BETTER_FIRST_RESTORE
      const size_t prefetchPos = MemoryCache::getPrefetchPos(neighbors->size() - position);
#else
      const size_t prefetchPos = MemoryCache::getPrefetchPos(neighbors->size());
#endif

      for (size_t i = 0; i < prefetchPos; i++) {
	if (!distanceChecked[(*(neighborptr + i)).id]) {
	  unsigned char *ptr = reinterpret_cast<unsigned char*>(objectRepository.get((*(neighborptr + i)).id));
	  MemoryCache::prefetch(ptr, byteSizeOfObject);
	}
      }
#ifdef NGT_GRAPH_BETTER_FIRST_RESTORE
      for (; neighborptr < neighborendptr; ++neighborptr, position++) {
#else
      for (; neighborptr < neighborendptr; ++neighborptr) {
#endif
	if ((neighborptr + prefetchPos < neighborendptr) && !distanceChecked[(*(neighborptr + prefetchPos)).id]) {
	  unsigned char *ptr = reinterpret_cast<unsigned char*>(objectRepository.get((*(neighborptr + prefetchPos)).id));
	  MemoryCache::prefetch(ptr, byteSizeOfObject);
	}
#ifdef NGT_VISIT_COUNT
	sc.visitCount++;
#endif
	ObjectDistance &neighbor = *neighborptr;
	if (distanceChecked[neighbor.id]) {
	  continue;
	}
	distanceChecked.insert(neighbor.id);

#ifdef NGT_EXPLORATION_COEFFICIENT_OPTIMIZATION
	sc.explorationCoefficient = exp(-(double)distanceChecked.size() / 20000.0) / 10.0 + 1.0;
#endif

	Distance distance = comparator(sc.object, *objectRepository.get(neighbor.id));
	sc.distanceComputationCount++;
	if (distance <= explorationRadius) {
	  result.set(neighbor.id, distance);
	  unchecked.push(result);
	  if (distance <= sc.radius) {
	    results.push(result);
	    if (results.size() >= sc.size) {
	      if (results.top().distance >= distance) {
		if (results.size() > sc.size) {
		  results.pop();
		}
		sc.radius = results.top().distance;
		explorationRadius = sc.explorationCoefficient * sc.radius;
	      }
	    }
	  } 
#ifdef NGT_GRAPH_BETTER_FIRST_RESTORE
	  if ((distance < target.distance) && (distance <= explorationRadius) && ((neighborptr + 2) < neighborendptr)) {
	    target.position = position + 1;
	    unchecked.push(target);
	    break;
	  }
#endif
	} 
      } 

    } 
    if (sc.resultIsAvailable()) { 
      ObjectDistances &qresults = sc.getResult();
      qresults.clear();
      qresults.moveFrom(results);
    } else {
      sc.workingResult = std::move(results);
    }
  }


  void
  NeighborhoodGraph::removeEdgesReliably(ObjectID id) {
    GraphNode *nodetmp = 0;
    try {
      nodetmp = getNode(id);
    } catch (Exception &err) {
      stringstream msg;
      msg << "removeEdgesReliably : cannot find a node. ID=" << id;
      msg << ":" << err.what();
      NGTThrowException(msg.str());
    }
    if (nodetmp == 0) {
      stringstream msg;
      msg << "removeEdgesReliably : cannot find a node. ID=" << id;
      NGTThrowException(msg.str());
    }
    GraphNode &node = *nodetmp;
    if (node.size() == 0) {
      cerr << "removeEdgesReliably : Warning! : No edges. ID=" << id << endl;
      try {
	removeNode(id);
      } catch (Exception &err) {
	stringstream msg;
	msg << "removeEdgesReliably : Internal error! : cannot remove a node without edges. ID=" << id;
	msg << ":" << err.what();
	NGTThrowException(msg.str());
      }
      return;
    }

    vector<PersistentObject*> objtbl;
    vector<GraphNode*> nodetbl;
    try {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      for (GraphNode::iterator i = node.begin(repository.allocator); i != node.end(repository.allocator);) {
#else
      for (GraphNode::iterator i = node.begin(); i != node.end();) {
#endif
	if (id == (*i).id) {
	  cerr << "Graph::removeEdgesReliably: Inner error. Destination nodes include a source node. ID="
	       << id << " continue..." << endl;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  i = node.erase(i, repository.allocator);
#else
	  i = node.erase(i);
#endif
	  continue;
	}
	objtbl.push_back(getObjectRepository().get((*i).id));
	GraphNode *n = 0;
	try {
	  n = getNode((*i).id);
	} catch (Exception &err) {
	  cerr << "Graph::removeEdgesReliably: Cannot find edges of a child. ID="
	       << (*i).id << " continue..." << endl;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  i = node.erase(i, repository.allocator);
#else
	  i = node.erase(i);
#endif
	  continue;
	}
	nodetbl.push_back(n);

	ObjectDistance edge;
	edge.id = id;
	edge.distance = (*i).distance;
	{
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  GraphNode::iterator ei = std::lower_bound(n->begin(repository.allocator), n->end(repository.allocator), edge);
	  if (ei != n->end(repository.allocator) && (*ei).id == id) {
	    n->erase(ei, repository.allocator);
#else
	  GraphNode::iterator ei = std::lower_bound(n->begin(), n->end(), edge);
	  if (ei != n->end() && (*ei).id == id) {
	    n->erase(ei);
#endif
	  } else {
	    stringstream msg;
	    msg << "removeEdgesReliably : internal error : cannot find an edge. ID="
		<< id << " d=" << edge.distance << " in " << (*i).id << endl;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    for (GraphNode::iterator ni = n->begin(repository.allocator); ni != n->end(repository.allocator); ni++) {
#else
	    for (GraphNode::iterator ni = n->begin(); ni != n->end(); ni++) {
#endif
	      msg << "check. " << (*ni).id << endl;
	    }
#ifdef NGT_FORCED_REMOVE
	    msg << " anyway continue...";
	    cerr << msg.str() << endl;
#else
	    NGTThrowException(msg.str());
#endif
	  }
	}
	i++;
      }
      for (unsigned int i = 0; i < node.size() - 1; i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	assert(node.at(i, repository.allocator).id != id);
#else
	assert(node[i].id != id);
#endif
	int minj = -1;
	Distance mind = FLT_MAX;
	for (unsigned int j = i + 1; j < node.size(); j++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  assert(node.at(j, repository.allocator).id != id);
#else
	  assert(node[j].id != id);
#endif
	  Distance d = objectSpace->getComparator()(*objtbl[i], *objtbl[j]);
	  if (d < mind) {
	    minj = j;
	    mind = d;
	  }
	}
	assert(minj != -1);
	{
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  ObjectDistance obj = node.at(minj, repository.allocator);
#else
	  ObjectDistance obj = node[minj];
#endif
	  obj.distance = mind;
	  GraphNode &n = *nodetbl[i];
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  GraphNode::iterator ei = std::lower_bound(n.begin(repository.allocator), n.end(repository.allocator), obj);
	  if ((*ei).id != obj.id) {
	    n.insert(ei, obj, repository.allocator);
	  }
#else
	  GraphNode::iterator ei = std::lower_bound(n.begin(), n.end(), obj);
	  if ((*ei).id != obj.id) {
	    n.insert(ei, obj);
	  }
#endif
	}
	{
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  ObjectDistance obj = node.at(i, repository.allocator);
#else
	  ObjectDistance obj = node[i];
#endif
	  obj.distance = mind;
	  GraphNode &n = *nodetbl[minj];
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  GraphNode::iterator ei = std::lower_bound(n.begin(repository.allocator), n.end(repository.allocator), obj);
	  if ((*ei).id != obj.id) {
	    n.insert(ei, obj, repository.allocator);
#else
	  GraphNode::iterator ei = std::lower_bound(n.begin(), n.end(), obj);
	  if ((*ei).id != obj.id) {
	    n.insert(ei, obj);
#endif
	  }
	}
	if ((i + 1 < node.size()) && (i + 1 != (unsigned int)minj)) {
	  ObjectDistance tmpr;
	  PersistentObject *tmpf;
	  GraphNode *tmprs;

#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  tmpr = node.at(i + 1, repository.allocator);
#else
	  tmpr = node[i + 1];
#endif
	  tmpf = objtbl[i + 1];
	  tmprs = nodetbl[i + 1];

#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  node.at(i + 1, repository.allocator) = node.at(minj, repository.allocator);
#else
	  node[i + 1] = node[minj];
#endif
	  objtbl[i + 1] = objtbl[minj];
	  nodetbl[i + 1] = nodetbl[minj];

#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  node.at(minj, repository.allocator) = tmpr;
#else
	  node[minj] = tmpr;
#endif
	  objtbl[minj] = tmpf;
	  nodetbl[minj] = tmprs;
	}
      } 

    } catch(Exception &err) {
      stringstream msg;
      msg << "removeEdgesReliably : Relink error ID=" << id << ":" << err.what();
#ifdef NGT_FORCED_REMOVE
      cerr << msg.str() << " continue..." << endl;
#else
      NGTThrowException(msg.str());
#endif
    }

    try {
      removeNode(id);
    } catch (Exception &err) {
      stringstream msg;
      msg << "removeEdgesReliably : removeEdges error. ID=" << id << ":" << err.what();
      NGTThrowException(msg.str());
    }
  }

class TruncationSearchJob {
public:
  TruncationSearchJob() {}
  TruncationSearchJob &operator=(TruncationSearchJob &d) {
    idx = d.idx;
    object = d.object;
    nearest = d.nearest;
    start = d.start;
    radius = d.radius;
    return *this;
  }
  size_t		idx;
  PersistentObject	*object;
  ObjectDistance	nearest;
  ObjectDistance	start;
  NGT::Distance	radius;
};

class TruncationSearchSharedData {
public:
  TruncationSearchSharedData(NGT::NeighborhoodGraph &g, NGT::ObjectID id, size_t size, NGT::Distance lr) : 
    graphIndex(g), targetID(id), resultSize(size), explorationCoefficient(lr) {}
  NGT::NeighborhoodGraph	&graphIndex;
  NGT::ObjectID		targetID;
  size_t			resultSize;
  NGT::Distance		explorationCoefficient;
};

class TruncationSearchThread : public NGT::Thread {
public:
  TruncationSearchThread() {}
  virtual ~TruncationSearchThread() {}
  virtual int run() {
    NGT::ThreadPool<TruncationSearchJob, TruncationSearchSharedData*, TruncationSearchThread>::Thread &poolThread =
      (NGT::ThreadPool<TruncationSearchJob, TruncationSearchSharedData*, TruncationSearchThread>::Thread&)*this;
    TruncationSearchSharedData &sd = *poolThread.getSharedData();
    for (;;) {
      TruncationSearchJob job;
      try {
	poolThread.getInputJobQueue().popFront(job);
      } catch(NGT::ThreadTerminationException &err) {
	break;
      } catch (NGT::Exception &err) {
	cerr << "TruncationSearchThread::run()::Inner error. continue..." << endl;
      }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      Object *po = sd.graphIndex.objectSpace->allocateObject((Object&)*job.object);
      NGT::SearchContainer ssc(*po);
#else
      NGT::SearchContainer ssc(*job.object);
#endif

      NGT::ObjectDistances srs, results;

      srs.push_back(job.start);
      ssc.setResults(&results);
      ssc.size = sd.resultSize;
      ssc.radius = job.radius;
      ssc.explorationCoefficient = sd.explorationCoefficient;
      ssc.id = 0;
      try {
	sd.graphIndex.search(ssc, srs);
      } catch(...) {
	cerr << "CreateIndexThread::run : Fatal Error!" << endl;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	sd.graphIndex.objectSpace->deleteObject(po);
#endif
      }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      sd.graphIndex.objectSpace->deleteObject(po);
#endif
      job.nearest = results[0];
      poolThread.getOutputJobQueue().pushBack(job);
    }
    return 0;
  }

};

typedef NGT::ThreadPool<TruncationSearchJob, TruncationSearchSharedData*, TruncationSearchThread> TruncationSearchThreadPool;

int
NeighborhoodGraph::truncateEdgesOptimally(
					  ObjectID id,
					  GraphNode &results,
					  size_t truncationSize
					  )
{

  ObjectDistances delNodes;

  size_t osize = results.size();

  size_t resSize = 2;
  TruncationSearchThreadPool threads(property.truncationThreadPoolSize);
  TruncationSearchSharedData sd(*this, id, resSize, 1.1);
  threads.setSharedData(&sd);
  threads.create();

  try {
    for (size_t i = truncationSize; i < results.size(); i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      if (id == results.at(i, repository.allocator).id) {
#else
      if (id == results[i].id) {
#endif
	continue;
      }
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      delNodes.push_back(results.at(i, repository.allocator));
#else
      delNodes.push_back(results[i]);
#endif
    }

#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    GraphNode::iterator ri = results.begin(repository.allocator);
    ri += truncationSize;
    results.erase(ri, results.end(repository.allocator), repository.allocator);
#else
    GraphNode::iterator ri = results.begin();
    ri += truncationSize;
    results.erase(ri, results.end());
#endif

    for (size_t i = 0; i < delNodes.size(); i++) {
      GraphNode::iterator j;
      GraphNode &res = *getNode(delNodes[i].id);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      for (j = res.begin(repository.allocator); j != res.end(repository.allocator); j++) {
#else
      for (j = res.begin(); j != res.end(); j++) {
#endif
	if ((*j).id == id) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  res.erase(j, repository.allocator);
#else
	  res.erase(j);
#endif
	  break;
	}
      }
    }
    bool retry = true;
    size_t maxResSize = osize * 2;
    size_t batchSize = 20;
    TruncationSearchThreadPool::OutputJobQueue &output = threads.getOutputJobQueue();
    TruncationSearchJob job;

    for (; retry == true; resSize = maxResSize) {
      retry = false;
      sd.resultSize = resSize;
      size_t nodeidx = 0;
      for (;;) {
	size_t nodeSize = 0;
	for (; nodeidx < delNodes.size(); nodeidx++) {
	  if (delNodes[nodeidx].id == 0) {
	    continue;
	  }
	  nodeSize++;
	  job.object = getObjectRepository().get(delNodes[nodeidx].id);
	  job.idx = nodeidx;
	  job.start.id = id;
	  job.start.distance = delNodes[nodeidx].distance;
	  job.radius = FLT_MAX;
	  threads.pushInputQueue(job);
	  if (nodeSize >= batchSize) {
	    break;
	  }
	}
	if (nodeSize == 0) {
	  break;
	}
	threads.waitForFinish();

	if (output.size() != nodeSize) {
	  nodeSize = output.size();
	}
	size_t cannotMoveCnt = 0;
	for (size_t i = 0; i < nodeSize; i++) {
	  TruncationSearchJob &ojob = output.front();
	  ObjectID nearestID = ojob.nearest.id;
	  size_t idx = ojob.idx;
	  if (nearestID == delNodes[idx].id) {
	    delNodes[idx].id = 0;
	    output.pop_front();
	    continue;
	  } else if (nearestID == id) {
	    cannotMoveCnt++;
	    if ((resSize < maxResSize) && (cannotMoveCnt > 1)) {
	      retry = true;
	      output.pop_front();
	      continue;
	    }
	  } else {
	  }

	  ObjectID tid = delNodes[idx].id;
	  delNodes[idx].id = 0;

	  GraphNode &delres = *getNode(tid);
	  {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    GraphNode::iterator ei = std::lower_bound(delres.begin(repository.allocator), delres.end(repository.allocator), ojob.nearest);
	    if ((*ei).id != ojob.nearest.id) {
	      delres.insert(ei, ojob.nearest, repository.allocator);
#else
	    GraphNode::iterator ei = std::lower_bound(delres.begin(), delres.end(), ojob.nearest);
	    if ((*ei).id != ojob.nearest.id) {
	      delres.insert(ei, ojob.nearest);
#endif
	    } else {
	      output.pop_front();
	      continue;
	    }
	  }
	  ObjectDistance r;
	  r.distance = ojob.nearest.distance;
	  r.id = tid;
	  if (nearestID != id) {
	    GraphNode &rs = *getNode(nearestID);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    rs.push_back(r, repository.allocator);
	    std::sort(rs.begin(repository.allocator), rs.end(repository.allocator));
#else
	    rs.push_back(r);
	    std::sort(rs.begin(), rs.end());
#endif
	  } else {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    results.push_back(r, repository.allocator);
	    std::sort(results.begin(repository.allocator), results.end(repository.allocator));
#else
	    results.push_back(r);
	    std::sort(results.begin(), results.end());
#endif
	  }
	  output.pop_front();
	}

      } 
    } 

    int cnt = 0;
    for (size_t i = 0; i < delNodes.size(); i++) {
      if (delNodes[i].id != 0) {
	cnt++;
      }
    }
    if (cnt != 0) {
      for (size_t i = 0; i < delNodes.size(); i++) {
	if (delNodes[i].id != 0) {
	}
      }
    }
    threads.terminate();
  } catch (Exception &err) {
    threads.terminate();
    Exception e(err);
    throw e;
  }

  size_t delsize = osize - results.size();

  return delsize;
}

