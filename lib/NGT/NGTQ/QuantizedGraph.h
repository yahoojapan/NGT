//
// Copyright (C) 2020 Yahoo Japan Corporation
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
    
    void construct(NGT::Index &ngtindex, NGTQ::Index &quantizedIndex, size_t maxNoOfEdges) {
      NGTQ::InvertedIndexEntry<uint16_t> invertedIndexObjects(numOfSubspaces);
      quantizedIndex.getQuantizer().extractInvertedIndexObject(invertedIndexObjects);
      
      NGT::GraphAndTreeIndex &index = static_cast<NGT::GraphAndTreeIndex&>(ngtindex.getIndex());
      NGT::NeighborhoodGraph &graph = static_cast<NGT::NeighborhoodGraph&>(index);

      NGT::GraphRepository &graphRepository = graph.repository;
      PARENT::resize(graphRepository.size());

      for (size_t id = 1; id < graphRepository.size(); id++) {
	NGT::GraphNode &node = *graphRepository.VECTOR::get(id);
	size_t numOfEdges = node.size() < maxNoOfEdges ? node.size() : maxNoOfEdges;
	(*this)[id].ids.reserve(numOfEdges);
	NGTQ::QuantizedObjectProcessingStream quantizedStream(quantizedIndex.getQuantizer(), numOfEdges);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	for (auto i = node.begin(graphRepository.allocator); i != node.end(graphRepository.allocator); ++i) {
	  if (distance(node.begin(graphRepository.allocator), i) >= static_cast<int64_t>(numOfEdges)) {
#else
	for (auto i = node.begin(); i != node.end(); i++) {
	  if (distance(node.begin(), i) >= static_cast<int64_t>(numOfEdges)) {
#endif
	    break;
	  }
	  if ((*i).id == 0) {
	    std::cerr << "something strange" << std::endl;
	    abort();
	    continue;
	  }
	  (*this)[id].ids.push_back((*i).id);
	  for (size_t idx = 0; idx < numOfSubspaces; idx++) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
            size_t dataNo = distance(node.begin(graphRepository.allocator), i);  
#else
            size_t dataNo = distance(node.begin(), i);  
#endif
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    abort();
#else
	    if (invertedIndexObjects[(*i).id].localID[idx] < 1 || invertedIndexObjects[(*i).id].localID[idx] > 16) {
	      std::cerr << "Fatal inner error! Invalid local centroid ID. ID=" << (*i).id << ":" << invertedIndexObjects[(*i).id].localID[idx] << std::endl;
	      abort();
	    }
	    quantizedStream.arrangeQuantizedObject(dataNo, idx, invertedIndexObjects[(*i).id].localID[idx] - 1);
#endif
	  }
	} 


	(*this)[id].objects = quantizedStream.compressIntoUint4();
      }
    }

    void serialize(std::ofstream &os, NGT::ObjectSpace *objspace = 0) {
      NGTQ::QuantizedObjectProcessingStream quantizedObjectProcessingStream(numOfSubspaces);
      uint64_t n = numOfSubspaces;
      NGT::Serializer::write(os, n);
      n = PARENT::size();
      NGT::Serializer::write(os, n);
      for (auto i = PARENT::begin(); i != PARENT::end(); ++i) {
	NGT::Serializer::write(os, (*i).ids);
	size_t streamSize = quantizedObjectProcessingStream.getUint4StreamSize((*i).ids.size());
	NGT::Serializer::write(os, static_cast<uint8_t*>((*i).objects), streamSize);
      }
    }

    void deserialize(std::ifstream &is, NGT::ObjectSpace *objectspace = 0) {
      try {
	NGTQ::QuantizedObjectProcessingStream quantizedObjectProcessingStream(numOfSubspaces);
	uint64_t n;
	NGT::Serializer::read(is, n);
	numOfSubspaces = n;
	NGT::Serializer::read(is, n);
	PARENT::resize(n);
	for (auto i = PARENT::begin(); i != PARENT::end(); ++i) {
	  NGT::Serializer::read(is, (*i).ids);
          size_t streamSize = quantizedObjectProcessingStream.getUint4StreamSize((*i).ids.size());
	  uint8_t *objectStream = new uint8_t[streamSize];
	  NGT::Serializer::read(is, objectStream, streamSize);
	  (*i).objects = objectStream;
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

    static size_t getNumberOfSubvectors(size_t dimension, size_t dimensionOfSubvector) {
      if (dimensionOfSubvector == 0) {
	dimensionOfSubvector = dimension > 400 ? 2 : 1;
	dimensionOfSubvector = (dimension % dimensionOfSubvector == 0) ? dimensionOfSubvector : 1;
      }
      if (dimension % dimensionOfSubvector != 0) {
	stringstream msg;
	msg << "Quantizer::getNumOfSubvectors: dimensionOfSubvector is invalid. " << dimension << " : " << dimensionOfSubvector << std::endl;
	NGTThrowException(msg);
      }
      return dimension / dimensionOfSubvector;
    }

    static void buildQuantizedGraph(const std::string indexPath, size_t maxNumOfEdges) {
      NGTQG::Index index(indexPath, maxNumOfEdges);
      index.save();
    }

    static void buildQuantizedObjects(const std::string quantizedIndexPath, NGT::ObjectSpace &objectSpace) {
      NGTQ::Index	quantizedIndex(quantizedIndexPath);
      NGTQ::Quantizer	&quantizer = quantizedIndex.getQuantizer();

      {
	std::vector<float> meanObject(objectSpace.getDimension(), 0);
	quantizedIndex.getQuantizer().globalCodebook.insert(meanObject);
	quantizedIndex.getQuantizer().globalCodebook.createIndex(8);
      }

      vector<pair<NGT::Object*, size_t> > objects;
      for (size_t id = 1; id < objectSpace.getRepository().size(); id++) {
	if (id % 100000 == 0) {
	  std::cerr << "Processed " << id << " objects." << std::endl;
	}
	std::vector<float> object;
	try {
	  objectSpace.getObject(id, object);
	} catch(...) {
	  continue;
	}
	quantizer.insert(object, objects, id);
      }
      if (objects.size() > 0) {
	quantizer.insert(objects);
      }

      quantizedIndex.save();
      quantizedIndex.close();
    }

    static void constructQuantizedGraphFrame(const std::string quantizedIndexPath, size_t dimension, size_t dimensionOfSubvector) {
      NGTQ::Property property;
      NGT::Property globalProperty;
      NGT::Property localProperty;

      property.threadSize = 24;
      property.globalRange = 0;
      property.localRange = 0;
      property.dataType = NGTQ::DataTypeFloat;
      property.distanceType = NGTQ::DistanceTypeL2;
      property.singleLocalCodebook = false;
      property.batchSize = 1000;
      property.centroidCreationMode = NGTQ::CentroidCreationModeStatic; 
      property.localCentroidCreationMode = NGTQ::CentroidCreationModeDynamicKmeans; 

      property.globalCentroidLimit = 1; 
      property.localCentroidLimit = 16; 
      property.localClusteringSampleCoefficient = 100; 
      property.localDivisionNo = NGTQG::Index::getNumberOfSubvectors(dimension, dimensionOfSubvector);
      property.dimension = dimension;

      globalProperty.edgeSizeForCreation = 10;
      globalProperty.edgeSizeForSearch = 40;
      globalProperty.indexType = NGT::Property::GraphAndTree;
      globalProperty.insertionRadiusCoefficient = 1.1;

      localProperty.indexType = NGT::Property::GraphAndTree;
      localProperty.insertionRadiusCoefficient = 1.1;

      NGTQ::Index::create(quantizedIndexPath, property, globalProperty, localProperty);

    }

    static void quantize(const std::string indexPath, float dimensionOfSubvector, size_t maxNumOfEdges) {
      NGT::Index	index(indexPath);
      NGT::ObjectSpace &objectSpace = index.getObjectSpace();


      {
	std::string quantizedIndexPath = indexPath + "/qg";
	struct stat st;
	if (stat(quantizedIndexPath.c_str(), &st) != 0) {
	  NGT::Property ngtProperty;
	  index.getProperty(ngtProperty);
	  //NGTQG::Command::CreateParameters createParameters(args, property.dimension);
	  constructQuantizedGraphFrame(quantizedIndexPath, ngtProperty.dimension, dimensionOfSubvector);
	  buildQuantizedObjects(quantizedIndexPath, objectSpace);
	  if (maxNumOfEdges != 0) {
	    buildQuantizedGraph(indexPath, maxNumOfEdges);
	  }
	}
      }
    }

    const std::string path;
    NGTQ::Index quantizedIndex;
    NGTQ::Index blobIndex;

    QuantizedGraphRepository quantizedGraph;

  }; 

} 

