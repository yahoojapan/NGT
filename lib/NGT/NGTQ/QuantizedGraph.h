//
// Copyright (C) 2016-2020 Yahoo Japan Corporation
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
#include	"NGT/NGTQ/Quantizer.h"

#if !defined(NGT_SHARED_MEMORY_ALLOCATOR) && !defined(NGTQ_SHARED_INVERTED_INDEX)


#define GLOBAL_SIZE	1

namespace NGTQG {
  class SearchContainer : public NGT::SearchContainer {
    public:
      SearchContainer(SearchContainer &sc, NGT::Object &f): NGT::SearchContainer(sc, f), resultExpansion(sc.resultExpansion) {}
      SearchContainer() {}
      
      void setResultExpansion(float re) { resultExpansion = re; }
      float resultExpansion;
  };
  
  class SearchQuery : public NGT::QueryContainer, public NGTQG::SearchContainer {
    public:
      template <typename QTYPE> SearchQuery(const std::vector<QTYPE> &q): NGT::QueryContainer(q) {}
  };
  
  class QuantizedNode {
  public:
    ~QuantizedNode() {
      delete[] static_cast<uint8_t*>(objects);
    }
    std::vector<uint32_t> ids;
    void *objects;
  };
  
  class QuantizedGraphRepository : public std::vector<QuantizedNode> {
    typedef std::vector<QuantizedNode> PARENT;
  public:
    QuantizedGraphRepository(NGTQ::Index &quantizedIndex): numOfSubspaces(quantizedIndex.getQuantizer().property.localDivisionNo) {}
    ~QuantizedGraphRepository() {}

    void *get(size_t id) {
      return PARENT::at(id).objects;
    }

    std::vector<uint32_t>& getIDs(size_t id) {
      return PARENT::at(id).ids;
    }
    

    size_t getNumOfPaddedUint8Objects(size_t noEdges) {
      if (noEdges == 0) {
	return 0;
      }
      return ((noEdges - 1) / (NGTQ_SIMD_BLOCK_SIZE * NGTQ_BATCH_SIZE) + 1) * (NGTQ_SIMD_BLOCK_SIZE * NGTQ_BATCH_SIZE);
    }
    size_t getNumOfPaddedUint4Objects(size_t noUint8Objects) {
      if (noUint8Objects == 0) {
	return 0;
      }
      return (noUint8Objects * numOfSubspaces) / 2 + 1;
    }

    void construct(NGT::Index &ngtindex, NGTQ::Index &quantizedIndex, size_t maxNoOfEdges) {
      NGTQ::InvertedIndexEntry<uint16_t> invertedIndexObjects(numOfSubspaces);
      quantizedIndex.getQuantizer().extractInvertedIndexObject(invertedIndexObjects);

      NGT::GraphAndTreeIndex &index = static_cast<NGT::GraphAndTreeIndex&>(ngtindex.getIndex());
      NGT::NeighborhoodGraph &graph = static_cast<NGT::NeighborhoodGraph&>(index);

      NGT::GraphRepository &graphRepository = graph.repository;
      PARENT::resize(graphRepository.size());

      for (size_t id = 1; id < graphRepository.size(); id++) {
	NGT::GraphNode &node = *graphRepository.VECTOR::get(id);
	size_t noEdges = node.size() < maxNoOfEdges ? node.size() : maxNoOfEdges;
	(*this)[id].ids.reserve(noEdges);
	size_t noObjects = getNumOfPaddedUint8Objects(noEdges);
	uint8_t *objectData = new uint8_t[noObjects * numOfSubspaces]();
	for (auto i = node.begin(); i != node.end(); i++) {
	  if (distance(node.begin(), i) >= static_cast<int64_t>(noEdges)) {
	    break;
	  }
	  if ((*i).id == 0) {
	    std::cerr << "something strange" << std::endl;
	    abort();
	    continue;
	  }
	  (*this)[id].ids.push_back((*i).id);
	  for (size_t idx = 0; idx < numOfSubspaces; idx++) {
            size_t dataNo = distance(node.begin(), i);  
	    size_t blkNo = dataNo / NGTQ_SIMD_BLOCK_SIZE;	
	    size_t oft = dataNo - blkNo * NGTQ_SIMD_BLOCK_SIZE;	
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    abort();
#else
	    objectData[blkNo * (NGTQ_SIMD_BLOCK_SIZE * numOfSubspaces) + NGTQ_SIMD_BLOCK_SIZE * idx + oft] = invertedIndexObjects[(*i).id].localID[idx] - 1;
#endif
	  }
	} 

	{
	  size_t idx = 0;
	  uint8_t *uint4Objects = new uint8_t[getNumOfPaddedUint4Objects(noObjects)]();
	  for (size_t nidx = 0; nidx < noObjects; nidx += NGTQ_SIMD_BLOCK_SIZE * NGTQ_BATCH_SIZE) {
	    for (size_t bcnt = 0; bcnt < NGTQ_BATCH_SIZE; bcnt++) {
	      for (size_t lidx = 0; lidx < numOfSubspaces; lidx++) {
		for (size_t bidx = 0; bidx < NGTQ_SIMD_BLOCK_SIZE; bidx++) {
		  if (idx % 2 == 0) {
		    uint4Objects[idx / 2] = objectData[idx];
		  } else {
		    uint4Objects[idx / 2] |= (objectData[idx] << 4);
		  }
		  idx++;
		}
	      }
	    }
	  }
	  delete[] objectData;
	  (*this)[id].objects = uint4Objects;
        }
      }
    }
    void serialize(std::ofstream &os, NGT::ObjectSpace *objspace = 0) {
      uint64_t n = numOfSubspaces;
      NGT::Serializer::write(os, n);
      n = PARENT::size();
      NGT::Serializer::write(os, n);
      for (auto i = PARENT::begin(); i != PARENT::end(); ++i) {
	NGT::Serializer::write(os, (*i).ids);
	size_t noObjects = getNumOfPaddedUint4Objects(getNumOfPaddedUint8Objects((*i).ids.size()));
	NGT::Serializer::write(os, static_cast<uint8_t*>((*i).objects), noObjects);
      }
    }
    void deserialize(std::ifstream &is, NGT::ObjectSpace *objectspace = 0) {
      try {
	uint64_t n;
	NGT::Serializer::read(is, n);
	numOfSubspaces = n;
	NGT::Serializer::read(is, n);
	PARENT::resize(n);
	for (auto i = PARENT::begin(); i != PARENT::end(); ++i) {
	  NGT::Serializer::read(is, (*i).ids);
	  size_t noObjects = getNumOfPaddedUint4Objects(getNumOfPaddedUint8Objects((*i).ids.size()));
	  uint8_t *objects = new uint8_t[noObjects];
	  NGT::Serializer::read(is, objects, noObjects);
	  (*i).objects = objects;
	}
      } catch(NGT::Exception &err) {
	std::stringstream msg;
	msg << "QuantizedGraph::deserialize: Fatal error. " << err.what();
	NGTThrowException(msg);
      }
    }
    void save(const string &path) {
      const std::string p(path + "/grp");
      std::ofstream os(p);
      serialize(os);
    }
    void load(const string &path) {
      const std::string p(path + "/grp");
      std::ifstream is(p);
      deserialize(is);
    }
    size_t numOfSubspaces;
  };


  class Index : public NGT::Index {
  public:
    Index(const std::string &indexPath, size_t maxNoOfEdges = 128) :
      NGT::Index(indexPath, false),
      path(indexPath),
      quantizedIndex(indexPath + "/qg"),
      quantizedGraph(quantizedIndex)  
      {
	{
	  struct stat st;
	  std::string qgpath(path + "/qg/grp");
	  if (stat(qgpath.c_str(), &st) == 0) {
	    quantizedGraph.load(path + "/qg");
	  } else {
	    quantizedGraph.construct(*this, quantizedIndex, maxNoOfEdges);
	  }
	}
      }

    void save() {
      quantizedGraph.save(path + "/qg");
    }

    
    void searchQuantizedGraph(NGT::NeighborhoodGraph &graph, NGTQG::SearchContainer &sc, NGT::ObjectDistances &seeds) {
      size_t sizeBackup = sc.size;
      if (sc.resultExpansion > 1.0) {
	sc.size *= sc.resultExpansion;
      }

      NGTQ::Quantizer &quantizer = quantizedIndex.getQuantizer();
      NGTQ::QuantizedObjectDistance &quantizedObjectDistance = quantizer.getQuantizedObjectDistance();

      NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8 cache[GLOBAL_SIZE + 1];
      for (int i = 1; i < GLOBAL_SIZE + 1; i++) {
	quantizedObjectDistance.initialize(cache[i]);
	quantizedObjectDistance.createDistanceLookup(sc.object, i, cache[i]);
      }

      if (sc.explorationCoefficient == 0.0) {
	sc.explorationCoefficient = NGT_EXPLORATION_COEFFICIENT;
      }
      NGT::NeighborhoodGraph::UncheckedSet unchecked;
      NGT::NeighborhoodGraph::DistanceCheckedSet distanceChecked(graph.repository.size());
      NGT::NeighborhoodGraph::ResultSet results;

      graph.setupDistances(sc, seeds, NGT::PrimitiveComparator::L2Float::compare);
      graph.setupSeeds(sc, seeds, results, unchecked, distanceChecked);
      NGT::Distance explorationRadius = sc.explorationCoefficient * sc.radius;
      NGT::ObjectDistance result;
      NGT::ObjectDistance target;

      while (!unchecked.empty()) {
	target = unchecked.top();
	unchecked.pop();
	if (target.distance > explorationRadius) {
	  break;
	}
	auto &neighborIDs = quantizedGraph.getIDs(target.id);
	size_t neighborSize = neighborIDs.size();
	float ds[neighborSize + NGTQ_SIMD_BLOCK_SIZE];

#ifdef NGTQG_PREFETCH
	{
	  uint8_t *lid = static_cast<uint8_t*>(quantizedGraph.get(target.id));
	  size_t size = ((neighborSize - 1) / (NGTQ_SIMD_BLOCK_SIZE * NGTQ_BATCH_SIZE) + 1) * (NGTQ_SIMD_BLOCK_SIZE * NGTQ_BATCH_SIZE);
	  size /= 2;
	  size *= quantizedIndex.getQuantizer().divisionNo; 

	  NGT::MemoryCache::prefetch(lid, size);
	}
#endif  //// NGTQG_PREFETCH
	quantizedObjectDistance(quantizedGraph.get(target.id), ds, neighborSize, cache[1]);
	for (size_t idx = 0;idx < neighborSize; idx++) {
	  NGT::Distance distance = ds[idx];
	  auto objid = neighborIDs[idx];
	  if (distance <= explorationRadius) {
	    bool checked = distanceChecked[objid];
	    if (checked) {
	      continue;
	    }
#ifdef NGT_VISIT_COUNT
	    sc.visitCount++;
#endif
	    distanceChecked.insert(objid);
	    result.set(objid, distance);
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
	NGT::ObjectDistances &qresults = sc.getResult();
	qresults.moveFrom(results);
	if (sc.resultExpansion >= 1.0) {
	  {
	    NGT::ObjectRepository &objectRepository = NGT::Index::getObjectSpace().getRepository();
	    NGT::ObjectSpace::Comparator &comparator =  NGT::Index::getObjectSpace().getComparator();
	    for (auto i = qresults.begin(); i != qresults.end(); ++i) {
#ifdef NGTQG_PREFETCH
	      if (static_cast<size_t>(distance(qresults.begin(), i + 10)) < qresults.size()) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)    
		NGT::PersistentObject &o = *objectRepository.get((*(i + 10)).id);
#else
		NGT::Object &o = *objectRepository[(*(i + 10)).id];
#endif
		_mm_prefetch(&o[0], _MM_HINT_T0);
	      }
#endif
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	      NGT::PersistentObject &obj = *objectRepository.get((*i).id);
#else
	      NGT::Object &obj = *objectRepository[(*i).id];
#endif
	      (*i).distance = comparator(sc.object, obj);
	    }
	    std::sort(qresults.begin(), qresults.end());
	  }
	  sc.size = sizeBackup;
	  qresults.resize(sc.size);
	}
      } else {
	if (sc.resultExpansion >= 1.0) {
	  {
	    NGT::ObjectRepository &objectRepository = NGT::Index::getObjectSpace().getRepository();
	    NGT::ObjectSpace::Comparator &comparator =  NGT::Index::getObjectSpace().getComparator();
	    while (!sc.workingResult.empty()) { sc.workingResult.pop(); }
	    while (!results.empty()) {
	      NGT::ObjectDistance obj = results.top();
	      obj.distance = comparator(sc.object, *objectRepository.get(obj.id));
	      results.pop();
	      sc.workingResult.push(obj);
	    }
	    sc.size = sizeBackup;
	    while (sc.workingResult.size() > sc.size) { sc.workingResult.pop(); }
	  }
	} else {
	  sc.workingResult = std::move(results);
	}
      }

    }



    void search(NGT::GraphIndex &index, NGTQG::SearchContainer &sc, NGT::ObjectDistances &seeds) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      NGTThrowException("NGTQG is not available for SHARED.");
#endif
      if (index.getReadOnly()) {
	NGTThrowException("NGTQG is not available for readonly mode.");
      }
      if (sc.size == 0) {
	while (!sc.workingResult.empty()) sc.workingResult.pop();
	return;
      }
      if (seeds.size() == 0) {
	index.getSeedsFromGraph(index.repository, seeds);
      }
      if (sc.expectedAccuracy > 0.0) {
	sc.setEpsilon(getEpsilonFromExpectedAccuracy(sc.expectedAccuracy));
      }
      try {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR) || !defined(NGT_GRAPH_READ_ONLY_GRAPH)
	index.NeighborhoodGraph::search(sc, seeds);
#else
	searchQuantizedGraph(static_cast<NGT::NeighborhoodGraph&>(index), sc, seeds);
#endif
      } catch(NGT::Exception &err) {
	std::cerr << err.what() << std::endl;
	NGT::Exception e(err);
	throw e;
      }
    }

    void search(NGTQG::SearchQuery &sq) {
      NGT::GraphAndTreeIndex &index = static_cast<NGT::GraphAndTreeIndex&>(getIndex());
      NGT::Object *query = Index::allocateObject(sq.getQuery(), sq.getQueryType());
      try {
        NGTQG::SearchContainer sc(sq, *query);
        sc.distanceComputationCount = 0;
        sc.visitCount = 0;
	NGT::ObjectDistances	seeds;
	index.getSeedsFromTree(sc, seeds);
	NGTQG::Index::search(static_cast<NGT::GraphIndex&>(index), sc, seeds);
	sq.workingResult = std::move(sc.workingResult);
	sq.distanceComputationCount = sc.distanceComputationCount;
	sq.visitCount = sc.visitCount;
      } catch(NGT::Exception &err) {
	deleteObject(query);
	throw err;
      }
      deleteObject(query);
    }

    const std::string path;
    NGTQ::Index quantizedIndex;
    NGTQ::Index blobIndex;

    QuantizedGraphRepository quantizedGraph;

  }; 

} 

#endif 
