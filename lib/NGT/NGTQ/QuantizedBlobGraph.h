//
// Copyright (C) 2021 Yahoo Japan Corporation
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

#ifdef NGTQ_QBG
#include	"NGT/NGTQ/QuantizedGraph.h"

#include	<thread>


namespace QBG {
  class SearchContainer : public NGT::SearchContainer {
  public:
    SearchContainer(NGT::Object &q): NGT::SearchContainer(q),
      cutback(0.0), graphExplorationSize(50), exactResultSize(0),
      blobExplorationCoefficient(0.0), numOfProbes(0) {}
    SearchContainer(): NGT::SearchContainer(*reinterpret_cast<NGT::Object*>(0)),
      cutback(0.0), graphExplorationSize(50), exactResultSize(0),
      blobExplorationCoefficient(0.0), numOfProbes(0) {}
    SearchContainer(SearchContainer &sc, NGT::Object &q): NGT::SearchContainer(q) {
      QBG::SearchContainer::operator=(sc);
    }
    SearchContainer &operator=(SearchContainer &sc) {
      NGT::SearchContainer::operator=(sc);
      cutback = sc.cutback;
      graphExplorationSize = sc.graphExplorationSize;
      exactResultSize = sc.exactResultSize;
      blobExplorationCoefficient = sc.blobExplorationCoefficient;
      numOfProbes = sc.numOfProbes;
      objectVector = sc.objectVector;
      return *this;
    }
    void setCutback(float c) { cutback = c; }
    void setGraphExplorationSize(size_t size) { graphExplorationSize = size; }
    void setExactResultSize(size_t esize) { exactResultSize = esize; }
    void setBlobEpsilon(float c) { blobExplorationCoefficient = c + 1.0; }
    void setNumOfProbes(size_t p) { numOfProbes = p; }
    void setObjectVector(std::vector<float> &query) { objectVector = std::move(query); }
    float       cutback;
    size_t      graphExplorationSize;
    size_t      exactResultSize;
    float       blobExplorationCoefficient;
    size_t	numOfProbes;
    std::vector<float>	objectVector;
  };

  class QuantizedBlobGraphRepository : public NGTQG::QuantizedGraphRepository {
  public:
    QuantizedBlobGraphRepository(NGTQ::Index &quantizedIndex): NGTQG::QuantizedGraphRepository(quantizedIndex){
    }
    
    void construct(NGTQ::Index &quantizedIndex) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      std::cerr << "construct: Not implemented" << std::endl;
      abort();
#else
      
      (*this).resize(quantizedIndex.getGlobalCodebookSize());
      NGT::Timer timer;
      timer.start();
      for (size_t gid = 1; gid < quantizedIndex.getGlobalCodebookSize(); gid++) {
	if (gid % 100000 == 0) {
	  timer.stop();
	  std::cerr << "The number of processed blobs=" << gid << " VmSize=" <<  NGT::Common::getProcessVmSizeStr() << " Elapsed time=" << timer << std::endl;
	  timer.restart();
	}
	NGTQ::InvertedIndexEntry<uint16_t> invertedIndexObjects(numOfSubspaces);
	quantizedIndex.getQuantizer().extractInvertedIndexObject(invertedIndexObjects, gid);
	quantizedIndex.getQuantizer().eraseInvertedIndexObject(gid);
	NGTQ::QuantizedObjectProcessingStream quantizedStream(quantizedIndex.getQuantizer(), invertedIndexObjects.size());
	(*this)[gid].ids.reserve(invertedIndexObjects.size());
	for (size_t oidx = 0; oidx < invertedIndexObjects.size(); oidx++) {
	  (*this)[gid].ids.push_back(invertedIndexObjects[oidx].id);
	  for (size_t idx = 0; idx < numOfSubspaces; idx++) {
#ifdef NGTQ_UINT8_LUT
#ifdef NGTQ_SIMD_BLOCK_SIZE
            size_t dataNo = oidx;  
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    abort();
#else
	    quantizedStream.arrangeQuantizedObject(dataNo, idx, invertedIndexObjects[oidx].localID[idx] - 1);
#endif
#else 
	    objectData[idx * noobjs + dataNo] = invertedIndexObjects[oidx].localID[idx] - 1;
#endif 
#else  
	    objectData[idx * noobjs + dataNo] = invertedIndexObjects[oidx].localID[idx];
#endif 
	  }
	}

	(*this)[gid].subspaceID = invertedIndexObjects.subspaceID;
	(*this)[gid].objects = quantizedStream.compressIntoUint4();
      } 
#endif 
    }

  };
  
  class Index : public NGTQ::Index {
  public:
  Index(const std::string &indexPath, bool readOnly = false, bool silence = false) :
    NGTQ::Index(indexPath, readOnly), path(indexPath), quantizedBlobGraph(*this) {
      NGT::StdOstreamRedirector redirector(silence);
      redirector.begin();
      try {
	load();
      } catch (NGT::Exception &err) {
	if (readOnly) {
	} else {
	  quantizedBlobGraph.construct(*this);
	}
      }
      redirector.end();
    }

    ~Index() {}

    bool &getSilence() { return silence; };

    NGT::Object *allocateObject(std::vector<float> &objectVector) {
      auto &globalIndex = getQuantizer().globalCodebookIndex;
      auto dim = getQuantizer().property.dimension;
      objectVector.resize(dim, 0);
      return globalIndex.allocateObject(objectVector);
    }
    
    static void create(const std::string &index, NGTQ::Property &property, 
		       NGT::Property &globalProperty,
#ifdef NGTQ_QBG		       
		       NGT::Property &localProperty,
		       std::vector<float> *rotation,
		       const std::string &objectFile) {
#else
		       NGT::Property &localProperty) {
#endif
      property.quantizerType = NGTQ::QuantizerTypeQBG;
#ifdef NGTQ_QBG
      NGTQ::Index::create(index, property, globalProperty, localProperty, rotation, objectFile);
      try {
	NGT::Index::mkdir(index + getWorkspaceName());
      } catch(...) {}
#else
      NGTQ::Index::create(index, property, globalProperty, localProperty);
#endif
    }

    static void load(const std::string &indexPath, const std::vector<float> &rotation) {
      NGTQ::Index index(indexPath);
      index.getQuantizer().saveRotation(rotation);
    }

    void insert(const size_t id, std::vector<float> &object) {
      getQuantizer().objectList.put(id, object, &getQuantizer().globalCodebookIndex.getObjectSpace());
    }

    NGT::ObjectID append(std::vector<float> &object) {
      NGT::ObjectID id = getQuantizer().objectList.size();
      id = id == 0 ? 1 : id;
      std::cerr << "ID=" << getQuantizer().objectList.size();
      getQuantizer().objectList.put(id, object, &getQuantizer().globalCodebookIndex.getObjectSpace());
      std::cerr << ":" << getQuantizer().objectList.size() << std::endl;
      return id;
    }

    static void append(const std::string &indexName,	// index file
		       const std::string &data,	// data file
		       size_t dataSize = 0,	// data size
		       bool silence = false
		       ) {
      NGT::StdOstreamRedirector redirector(silence);
      redirector.begin();
      QBG::Index index(indexName);
      istream *is;
      if (data == "-") {
	is = &cin;
      } else {
	ifstream *ifs = new ifstream;
	ifs->ifstream::open(data);
	if (!(*ifs)) {
	  cerr << "Cannot open the specified file. " << data << endl;
	  return;
	}
	is = ifs;
      }
      string line;
      vector<pair<NGT::Object*, size_t> > objects;
      size_t count = 0;
      // extract objects from the file and insert them to the object list.
      while(getline(*is, line)) {
	count++;
	std::vector<float>	object;
	NGT::Common::extractVector(line, " ,\t", object);
	if (object.empty()) {
	  cerr << "An empty line or invalid value: " << line << endl;
	  continue;
	}
	index.insert(count, object);

	if (count % 100000 == 0) {
	  cerr << "Processed " << count;
	  cerr << endl;
	}
      }
      if (data != "-") {
	delete is;
      }

      index.save();
      index.close();
      redirector.end();
    }

    static void appendBinary(const std::string &indexName,	// index file
			     const std::string &data,	// data file
			     size_t dataSize = 0,	// data size
			     bool silence = false
		       ) {
      NGT::StdOstreamRedirector redirector(silence);
      redirector.begin();
      QBG::Index index(indexName);
      std::vector<std::string> tokens;
      NGT::Common::tokenize(data, tokens, ".");
      if (tokens.size() < 2) {
	std::stringstream msg;
	msg << "Invalid file name format";
	NGTThrowException(msg);
      }
      StaticObjectFileLoader loader(data, tokens[tokens.size() - 1]);
      size_t idx = 0;
      while (!loader.isEmpty()) {
	idx++;
	if (dataSize > 0 && idx > dataSize) {
	  break;
	}
	auto object = loader.getObject();
	index.insert(idx, object);
	if (idx % 100000 == 0) {
	  std::cerr << "loaded " << static_cast<float>(idx) / 1000000.0 << "M objects." << std::endl;
	}
      }
      index.save();
      index.close();
      redirector.end();
    }

    static void appendFromObjectRepository(const std::string &ngtIndex,	// QG
					   const std::string &qgIndex,	// NGT
					   bool silence = false) {
      NGT::StdOstreamRedirector redirector(silence);
      redirector.begin();

      NGT::Index ngt(ngtIndex);
      QBG::Index qg(qgIndex);
      auto &objectSpace = ngt.getObjectSpace();
      size_t size = objectSpace.getRepository().size();
      for (size_t id = 1; id < size; ++id) {
	std::vector<float> object;
	try {
	  objectSpace.getObject(id, object);
	} catch(...) {
	  std::cerr << "append: Info: removed object. " << id << std::endl;
	}
	qg.insert(id, object);
      }
      cerr << "end of insertion." << endl;
      qg.save();
      qg.close();
      redirector.end();
    }

    void getSeeds(NGT::Index &index, NGT::Object *object, NGT::ObjectDistances &seeds, size_t noOfSeeds) {
      auto &graph = static_cast<NGT::GraphAndTreeIndex&>(index.getIndex());
      NGT::SearchContainer sc(*object);
      sc.setResults(&seeds);
      sc.setSize(noOfSeeds);
      sc.setEpsilon(0.0);
      sc.setEdgeSize(-2);
      graph.search(sc);
    }

    NGT::Distance getDistance(void *objects, std::vector<float> &distances, size_t noOfObjects, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8 &lut
			      ) {
      auto &quantizedObjectDistance = getQuantizer().getQuantizedObjectDistance();
#ifdef NGTQBG_MIN
      auto min = quantizedObjectDistance(objects, distances.data(), noOfObjects, lut);
#else
      quantizedObjectDistance(objects, distances.data(), noOfObjects, lut);
#endif
#ifdef NGTQBG_MIN
      return min;
#endif
    }

    std::tuple<NGT::Distance, NGT::Distance>
      judge(NGTQG::QuantizedNode &ivi, size_t k, NGT::Distance radius,
	    NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8 &lut,
	    NGT::NeighborhoodGraph::ResultSet &result, size_t &foundCount
	    ) {
      auto noOfObjects = ivi.ids.size();
      float distances[NGTQ::QuantizedObjectProcessingStream::getNumOfAlignedObjects(noOfObjects)];
      auto &quantizedObjectDistance = getQuantizer().getQuantizedObjectDistance();
#ifdef NGTQBG_MIN
      float distance = quantizedObjectDistance(ivi.objects, &distances[0], noOfObjects, lut);
#else
      quantizedObjectDistance(ivi.objects, &distances[0], noOfObjects, lut);
#endif

#ifdef NGTQBG_MIN
      if (distance >= radius) {
	return std::make_pair(distance, radius);
      }
#endif
      bool found = false;
      for (size_t i = 0; i < noOfObjects; i++) {
	if (distances[i] <= radius) {
	  result.push(NGT::ObjectDistance(ivi.ids[i], distances[i]));
	  found = true;
	  if (result.size() > k) {
	    result.pop();
	  }
	}
      }
      if (result.size() >= k) {
	radius = result.top().distance;
      }
      if (found) foundCount++;
#ifdef NGTQBG_MIN
      return std::make_pair(distance, radius);
#else
      return std::make_pair(0.0, radius);
#endif
    }

    /////////////////// ///-/
    void searchBlobGraphNaively(QBG::SearchContainer &searchContainer) {
      NGT::Object *query = &searchContainer.object;

      auto &quantizer = getQuantizer();
      auto &globalIndex = quantizer.globalCodebookIndex;
      auto &globalGraph = static_cast<NGT::GraphAndTreeIndex&>(globalIndex.getIndex());
      NGT::ObjectDistances seeds;
      getSeeds(globalIndex, query, seeds, 5);

      if (seeds.empty()) {
	std::cerr << "something wrong." << std::endl;
	return;
      }
      size_t seedBlobID = seeds[0].id;
      auto &quantizedObjectDistance = quantizer.getQuantizedObjectDistance();
#ifdef NGTQG_ZERO_GLOBAL
      NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8 lut;
      quantizedObjectDistance.initialize(lut);
      quantizedObjectDistance.createDistanceLookup(*query, 1, lut); 
#else
#if defined(NGTQG_ROTATION)
      quantizedObjectDistance.rotation->mul(static_cast<float*>(query->getPointer()));
#endif
      uint32_t subspaceID = quantizedBlobGraph[seedBlobID].subspaceID;
      
      std::unordered_map<size_t, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8> luts;
      luts.insert(std::make_pair(subspaceID, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8()));
      auto lutfi = luts.find(subspaceID);
      quantizedObjectDistance.initialize((*lutfi).second);
      quantizedObjectDistance.createDistanceLookup(*query, subspaceID, (*lutfi).second);

#endif
      size_t visitCount = 1;
      size_t foundCount = 0;

      size_t k = searchContainer.size;
      NGT::Distance radius = FLT_MAX;
      NGT::Distance distance;
      NGT::NeighborhoodGraph::ResultSet result;
#ifdef NGTQG_ZERO_GLOBAL
      std::tie(distance, radius) = judge(quantizedBlobGraph[seedBlobID], k, radius, lut, result, foundCount);
#else
      std::tie(distance, radius) = judge(quantizedBlobGraph[seedBlobID], k, radius, (*lutfi).second, result, foundCount);
#endif
      NGT::NeighborhoodGraph::UncheckedSet uncheckedBlobs;
      NGT::NeighborhoodGraph::DistanceCheckedSet distanceCheckedBlobs(globalGraph.repository.size());
      distanceCheckedBlobs.insert(seedBlobID);
      if (globalGraph.searchRepository.size() == 0) {
	std::cerr << "graph is empty! Is it read only?" << std::endl;
	abort();
      }
      auto *nodes = globalGraph.searchRepository.data();
      uncheckedBlobs.push(NGT::ObjectDistance(seedBlobID, distance));
      float explorationRadius = radius * searchContainer.explorationCoefficient;
      while (!uncheckedBlobs.empty()) {
	auto targetBlob = uncheckedBlobs.top();
	if (targetBlob.distance > explorationRadius) {
	  break;
	}
	uncheckedBlobs.pop();
	auto *neighbors = nodes[targetBlob.id].data();
	auto noOfEdges = (searchContainer.edgeSize == 0 || searchContainer.edgeSize > static_cast<int64_t>(nodes[targetBlob.id].size())) ?
	                 nodes[targetBlob.id].size() : searchContainer.edgeSize;
	auto neighborend = neighbors + noOfEdges;
;
	for (auto neighbor = neighbors; neighbor < neighborend; neighbor++) {
	  NGT::ObjectID neighborID = neighbor->first;
	  if (distanceCheckedBlobs[neighborID]) {
	    continue;
	  }
	  visitCount++;
	  distanceCheckedBlobs.insert(neighborID);
#ifdef NGTQG_ZERO_GLOBAL
	  std::tie(distance, radius) = judge(quantizedBlobGraph[neighborID], k, radius, lut, result, foundCount);
#else
          auto luti = luts.find(subspaceID);
          if (luti == luts.end()) {
	    luts.insert(std::make_pair(subspaceID, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8()));
	    luti = luts.find(subspaceID);
	    quantizedObjectDistance.initialize((*luti).second);
	    quantizedObjectDistance.createDistanceLookup(*query, subspaceID, (*luti).second);
	  }
	  std::tie(distance, radius) = judge(quantizedBlobGraph[neighborID], k, radius, (*luti).second, result, foundCount);
#endif 
	  if (static_cast<float>(foundCount) / visitCount < searchContainer.cutback) {
	    uncheckedBlobs = NGT::NeighborhoodGraph::UncheckedSet();
	    break;
	  }
	  if (distance <= explorationRadius) {
	    uncheckedBlobs.push(NGT::ObjectDistance(neighborID, distance));
	  }
	}
      }

      if (searchContainer.resultIsAvailable()) { 
	searchContainer.getResult().clear();
	searchContainer.getResult().moveFrom(result);
      } else {
	searchContainer.workingResult = result;
      }


    }


    void searchBlobNaively(QBG::SearchContainer &searchContainer) {
      NGT::ObjectDistances blobs;
      NGT::SearchContainer sc(searchContainer);
      sc.setResults(&blobs);
      sc.setSize(searchContainer.numOfProbes);

      auto &quantizer = getQuantizer();
      auto &globalIndex = quantizer.globalCodebookIndex;
      auto &objectSpace = globalIndex.getObjectSpace();
      globalIndex.search(sc);


      if (blobs.empty()) {
	std::cerr << "something wrong." << std::endl;
	return;
      }

      auto &quantizedObjectDistance = quantizer.getQuantizedObjectDistance();
      NGT::Object rotatedQuery(&objectSpace);
      objectSpace.copy(rotatedQuery, searchContainer.object);

#if defined(NGTQG_ROTATION)
      quantizedObjectDistance.rotation->mul(static_cast<float*>(rotatedQuery.getPointer()));
#endif
      std::unordered_map<size_t, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8> luts;
      size_t foundCount = 0;

      size_t k = searchContainer.size;
      NGT::Distance radius = FLT_MAX;
      NGT::Distance distance;
      NGT::NeighborhoodGraph::ResultSet result;

      for (size_t idx = 0; idx < blobs.size(); idx++) {
	auto blobID = blobs[idx].id;
	auto subspaceID = quantizedBlobGraph[blobID].subspaceID;
	auto luti = luts.find(subspaceID);
	if (luti == luts.end()) {
	  luts.insert(std::make_pair(subspaceID, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8()));
	  luti = luts.find(subspaceID);
	  quantizedObjectDistance.initialize((*luti).second);
	  quantizedObjectDistance.createDistanceLookup(rotatedQuery, subspaceID, (*luti).second);
	}
	std::tie(distance, radius) = judge(quantizedBlobGraph[blobID], k, radius, (*luti).second, result, foundCount);

      }

      if (searchContainer.resultIsAvailable()) { 
	searchContainer.getResult().clear();
	searchContainer.getResult().moveFrom(result);
      } else {
	searchContainer.workingResult = result;
      }


    }

   void searchBlobGraph(QBG::SearchContainer &searchContainer) {
     auto &globalIndex = getQuantizer().globalCodebookIndex;
     auto &globalGraph = static_cast<NGT::GraphAndTreeIndex&>(globalIndex.getIndex());
     NGT::ObjectDistances	seeds;
     NGT::Object *query = allocateObject(searchContainer.objectVector);
     SearchContainer sc(searchContainer, *query);
     globalGraph.getSeedsFromTree(sc, seeds);
     if (seeds.empty()) {
       globalGraph.getRandomSeeds(globalGraph.repository, seeds, 20);
     }
     searchBlobGraph(sc, seeds);
     globalIndex.deleteObject(query);
     searchContainer.workingResult = std::move(sc.workingResult);
   }

   void searchBlobGraph(QBG::SearchContainer &searchContainer, NGT::ObjectDistances &seeds) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
     std::cerr << "searchBlobGraph: Not implemented. " << std::endl;
     abort();
#else

    auto &quantizer = getQuantizer();
    auto &globalIndex = quantizer.globalCodebookIndex;
    auto &globalGraph = static_cast<NGT::GraphAndTreeIndex&>(globalIndex.getIndex());
    auto &objectSpace = globalIndex.getObjectSpace();

    if (searchContainer.explorationCoefficient == 0.0) {
      searchContainer.explorationCoefficient = NGT_EXPLORATION_COEFFICIENT;
    }

    const auto requestedSize = searchContainer.size;
    searchContainer.size = std::numeric_limits<uint32_t>::max();

    // setup edgeSize
    size_t edgeSize = globalGraph.getEdgeSize(searchContainer);

    NGT::NeighborhoodGraph::UncheckedSet untracedNodes;

    NGT::NeighborhoodGraph::DistanceCheckedSet distanceChecked(globalGraph.searchRepository.size());
    NGT::NeighborhoodGraph::ResultSet results;

    if (objectSpace.getObjectType() == typeid(float)) {
      globalGraph.setupDistances(searchContainer, seeds, NGT::PrimitiveComparator::L2Float::compare);
#ifdef NGT_HALF_FLOAT
    } else if (objectSpace.getObjectType() == typeid(NGT::float16)) {
      globalGraph.setupDistances(searchContainer, seeds, NGT::PrimitiveComparator::L2Float16::compare);
    }
#endif
    std::sort(seeds.begin(), seeds.end());
    NGT::ObjectDistance currentNearestBlob = seeds.front();
    NGT::Distance explorationRadius = searchContainer.blobExplorationCoefficient * currentNearestBlob.distance;
    std::priority_queue<NGT::ObjectDistance, std::vector<NGT::ObjectDistance>, std::greater<NGT::ObjectDistance>> discardedObjects;
    untracedNodes.push(seeds.front());
    distanceChecked.insert(seeds.front().id);
    for (size_t i = 1; i < seeds.size(); i++) {
      untracedNodes.push(seeds[i]);
      distanceChecked.insert(seeds[i].id);
      discardedObjects.push(seeds[i]);
    }
    size_t explorationSize = 1;
    auto &quantizedObjectDistance = quantizer.getQuantizedObjectDistance();
    std::unordered_map<size_t, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8> luts;
    std::vector<float> &rotatedQuery = searchContainer.objectVector;
    if (objectSpace.getObjectType() == typeid(float)) {
      memcpy(rotatedQuery.data(), searchContainer.object.getPointer(), rotatedQuery.size() * sizeof(float));
#ifdef NGT_HALF_FLOAT
    } else if (objectSpace.getObjectType() == typeid(NGT::float16)) {
      auto *ptr = static_cast<NGT::float16*>(searchContainer.object.getPointer());   
      for (size_t i = 0; i < rotatedQuery.size(); i++) {
	rotatedQuery[i] = ptr[i];
      }
#endif
    } else {
      std::cerr << "Fatal inner error! Invalid object type." << std::endl;
    }
    quantizedObjectDistance.rotation->mul(rotatedQuery.data());
    NGT::Distance radius = searchContainer.radius;
    if (requestedSize >= std::numeric_limits<int32_t>::max()) {
      radius *= searchContainer.explorationCoefficient;
    }
    const size_t dimension = objectSpace.getPaddedDimension();
    if (globalGraph.searchRepository.empty()) {
      NGTThrowException("QBG:Index: searchRepository is empty.");
    }
    NGT::ReadOnlyGraphNode *nodes = globalGraph.searchRepository.data();
    NGT::ReadOnlyGraphNode *neighbors = 0;
    NGT::ObjectDistance target;
    const size_t prefetchSize = objectSpace.getPrefetchSize();
    const size_t prefetchOffset = objectSpace.getPrefetchOffset();
    pair<uint64_t, NGT::PersistentObject*> *neighborptr;
    pair<uint64_t, NGT::PersistentObject*> *neighborendptr;
    for (;;) {
      if (untracedNodes.empty() || untracedNodes.top().distance > explorationRadius) {
	explorationSize++;
	auto blobID = currentNearestBlob.id;
	auto subspaceID = quantizedBlobGraph[blobID].subspaceID;
	auto luti = luts.find(subspaceID);
	if (luti == luts.end()) {
	  luts.insert(std::make_pair(subspaceID, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8()));
	  luti = luts.find(subspaceID);
	  quantizedObjectDistance.initialize((*luti).second);
	  quantizedObjectDistance.createDistanceLookup(rotatedQuery.data(), subspaceID, (*luti).second);
	}
	NGT::Distance blobDistance;
	size_t foundCount;
	std::tie(blobDistance, radius) = judge(quantizedBlobGraph[blobID], requestedSize, 
					       radius, (*luti).second, results, foundCount);
#ifdef NGTQBG_MIN
	if (blobDistance > radius * searchContainer.explorationCoefficient) {
	  break;
        }
#endif 
	if (explorationSize > searchContainer.graphExplorationSize) {
	  break;
	}
	if (discardedObjects.empty()) {
	  break;
	}
	currentNearestBlob = discardedObjects.top();
	discardedObjects.pop();
	explorationRadius = searchContainer.blobExplorationCoefficient * currentNearestBlob.distance;
	continue;
      }
      target = untracedNodes.top();
      untracedNodes.pop();

      neighbors = &nodes[target.id];
      neighborptr = &(*neighbors)[0];
      size_t neighborSize = neighbors->size() < edgeSize ? neighbors->size() : edgeSize;
      neighborendptr = neighborptr + neighborSize;

      pair<uint64_t, NGT::PersistentObject*>* nsPtrs[neighborSize];
      size_t nsPtrsSize = 0;
#ifndef PREFETCH_DISABLE
      for (; neighborptr < neighborendptr; ++neighborptr) {
#ifdef NGT_VISIT_COUNT
	searchContainer.visitCount++;
#endif
	if (!distanceChecked[(*(neighborptr)).first]) {
	  distanceChecked.insert((*(neighborptr)).first);
          nsPtrs[nsPtrsSize] = neighborptr;
          if (nsPtrsSize < prefetchOffset) {
            unsigned char *ptr = reinterpret_cast<unsigned char*>((*(neighborptr)).second);
	    NGT::MemoryCache::prefetch(ptr, prefetchSize);
          }
          nsPtrsSize++;
        }
      }
#endif
#ifdef PREFETCH_DISABLE
      for (; neighborptr < neighborendptr; ++neighborptr) {
#else
      for (size_t idx = 0; idx < nsPtrsSize; idx++) {
#endif
#ifdef PREFETCH_DISABLE
	if (distanceChecked[(*(neighborptr)).first]) {
	  continue;
	}
	distanceChecked.insert((*(neighborptr)).first);
#else
	neighborptr = nsPtrs[idx]; 
	if (idx + prefetchOffset < nsPtrsSize) {
	  unsigned char *ptr = reinterpret_cast<unsigned char*>((*(nsPtrs[idx + prefetchOffset])).second);
	  NGT::MemoryCache::prefetch(ptr, prefetchSize);
	}
#endif
#ifdef NGT_DISTANCE_COMPUTATION_COUNT
	searchContainer.distanceComputationCount++;
#endif
	NGT::Distance distance = 0.0;
	if (objectSpace.getObjectType() == typeid(float)) {
	  distance = NGT::PrimitiveComparator::L2Float::compare(searchContainer.object.getPointer(), 
							        neighborptr->second->getPointer(), dimension);
#ifdef NGT_HALF_FLOAT
	} else if (objectSpace.getObjectType() == typeid(NGT::float16)) { 
	  distance = NGT::PrimitiveComparator::L2Float16::compare(searchContainer.object.getPointer(), 
							           neighborptr->second->getPointer(), dimension);
#endif
	} else {
	  assert(false);
	}
	NGT::ObjectDistance r;
	r.set(neighborptr->first, distance);
	untracedNodes.push(r);
	if (distance < currentNearestBlob.distance) {
	  discardedObjects.push(currentNearestBlob);
	  currentNearestBlob = r;
	  explorationRadius = searchContainer.blobExplorationCoefficient * currentNearestBlob.distance;
	} else {
	  discardedObjects.push(r);
	}
      } 
    } 

    if (searchContainer.resultIsAvailable()) { 
      if (searchContainer.exactResultSize > 0) {
	NGT::ObjectDistances &qresults = searchContainer.getResult();
	auto threadid = omp_get_thread_num();
	auto paddedDimension = getQuantizer().globalCodebookIndex.getObjectSpace().getPaddedDimension();
	NGT::ResultPriorityQueue	rs;
	std::vector<float> object;
	qresults.resize(results.size());
	size_t idx = results.size();
	while (!results.empty()) {
	  auto r = results.top();
	  results.pop();
#ifdef MULTIPLE_OBJECT_LISTS
	  quantizer.objectList.get(threadid, r.id, object, &quantizer.globalCodebookIndex.getObjectSpace());
#else
	  quantizer.objectList.get(r.id, object, &quantizer.globalCodebookIndex.getObjectSpace());
#endif
	  if (objectSpace.getObjectType() == typeid(float)) {
	    r.distance = NGT::PrimitiveComparator::compareL2(static_cast<float*>(searchContainer.object.getPointer()),
	    						     static_cast<float*>(object.data()), paddedDimension);
#ifdef NGT_HALF_FLOAT
	  } else if (objectSpace.getObjectType() == typeid(NGT::float16)) {
	    r.distance = NGT::PrimitiveComparator::compareL2(reinterpret_cast<NGT::float16*>(searchContainer.object.getPointer()),
	    						     reinterpret_cast<NGT::float16*>(object.data()), paddedDimension);
#endif
	  }
	  qresults[--idx] = r;
	}
	std::sort(qresults.begin(), qresults.end());
	qresults.resize(searchContainer.exactResultSize);
      } else {
	NGT::ObjectDistances &qresults = searchContainer.getResult();
	qresults.moveFrom(results);
      }
    } else {
      if (searchContainer.exactResultSize > 0) {
	auto threadid = omp_get_thread_num();
	auto paddedDimension = getQuantizer().globalCodebookIndex.getObjectSpace().getPaddedDimension();
	NGT::ResultPriorityQueue	rs;
	std::vector<float> object;
	while (!results.empty()) {
	  auto r = results.top();
	  results.pop();
#ifdef MULTIPLE_OBJECT_LISTS
	  quantizer.objectList.get(threadid, r.id, object, &quantizer.globalCodebookIndex.getObjectSpace());
#else
	  quantizer.objectList.get(r.id, object, &quantizer.globalCodebookIndex.getObjectSpace());
#endif
	  if (objectSpace.getObjectType() == typeid(float)) {
	    r.distance = NGT::PrimitiveComparator::compareL2(static_cast<float*>(searchContainer.object.getPointer()),
	    						     static_cast<float*>(object.data()), paddedDimension);
#ifdef NGT_HALF_FLOAT
	  } else if (objectSpace.getObjectType() == typeid(NGT::float16)) {
	    r.distance = NGT::PrimitiveComparator::compareL2(reinterpret_cast<NGT::float16*>(searchContainer.object.getPointer()),
	    						     reinterpret_cast<NGT::float16*>(object.data()), paddedDimension);
#endif
	  }
	  rs.push(r);
	}
	results = std::move(rs);
      } else {
	searchContainer.workingResult = std::move(results);
      }
    }
#endif 
   }

    void save() {
      quantizedBlobGraph.save(path);
    }

    void load() {
      if (quantizedBlobGraph.stat(path)) {
	quantizedBlobGraph.load(path);
      } else {
	NGTThrowException("Not found the rearranged inverted index. [" + path + "]");
      }
    }

    static void build(const std::string &indexPath, bool silence = false) {
      build(indexPath, "", "", "", 1, 0, silence);
    }

    static void build(const std::string &indexPath,
		      std::string quantizerCodebookFile = "",
		      std::string codebookIndexFile = "",
		      std::string objectIndexFile = "",
		      size_t beginID = 1, size_t endID = 0, bool silence = false) {
      buildNGTQ(indexPath, quantizerCodebookFile, codebookIndexFile, objectIndexFile, beginID, endID, silence);
      buildQBG(indexPath, silence);
      std::cerr << "NGTQ and NGTQBG indices are completed." << std::endl;
      std::cerr << "  vmsize=" << NGT::Common::getProcessVmSizeStr() << std::endl;
      std::cerr << "  peak vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;

    }

    static void build(const std::string &indexPath,
		      std::vector<std::vector<float>> &quantizerCodebook,
		      std::vector<uint32_t> &codebookIndex,
		      std::vector<uint32_t> &objectIndex,
		      size_t beginID = 1, size_t endID = 0) {
      buildNGTQ(indexPath, quantizerCodebook, codebookIndex, objectIndex, beginID, endID);
      buildQBG(indexPath);
      std::cerr << "NGTQ and NGTQBG indices are completed." << std::endl;
      std::cerr << "  vmsize=" << NGT::Common::getProcessVmSizeStr() << std::endl;
      std::cerr << "  peak vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;
    }

    static void buildNGTQ(const std::string &indexPath, bool silence = false) {
      buildNGTQ(indexPath, "", "-", "-", 1, 0, silence);
    }

    static void buildNGTQ(const std::string &indexPath,
			  std::string quantizerCodebookFile = "",
			  std::string codebookIndexFile = "",
			  std::string objectIndexFile = "",
			  size_t beginID = 1, size_t endID = 0, bool silence = false) {
      std::vector<std::vector<float>> quantizerCodebook;
      std::vector<uint32_t> codebookIndex;
      std::vector<uint32_t> objectIndex;
      {
	std::string codebookPath = quantizerCodebookFile;
	if (codebookPath.empty()) {
	  codebookPath = QBG::Index::getQuantizerCodebookFileName(indexPath);
	}
	std::ifstream stream(codebookPath);
	if (!stream) {
	  std::stringstream msg;
	  msg << "Cannot open the codebook. " << codebookPath;
	  NGTThrowException(msg);
	}
	std::string line;
	while (getline(stream, line)) {
	  std::vector<std::string> tokens;
	  NGT::Common::tokenize(line, tokens, " \t");
	  std::vector<float> object;
	  for (auto &token : tokens) {
	    object.push_back(NGT::Common::strtof(token));
	  }
	  if (!quantizerCodebook.empty() && quantizerCodebook[0].size() != object.size()) {
	    std::stringstream msg;
	    msg << "The specified quantizer codebook is invalid. " << quantizerCodebook[0].size()
		<< ":" << object.size() << ":" << quantizerCodebook.size() << ":" << line;
	    NGTThrowException(msg);
	  }
	  if (!object.empty()) {
	    quantizerCodebook.push_back(object);
	  }
	}
      }
      {
        std::string codebookIndexPath = codebookIndexFile;
	if (codebookIndexPath.empty()) {
	  codebookIndexPath = QBG::Index::getCodebookIndexFileName(indexPath);
	}
	if (codebookIndexPath != "-") {
	  cerr << "buildNGTQ: codebook index is " << codebookIndexPath << "." << endl;
	  std::ifstream stream(codebookIndexPath);
	  if (!stream) {
	    std::stringstream msg;
	    msg << "Cannot open the codebook index. " << codebookIndexPath;
	    NGTThrowException(msg);
	  }
	  std::string line;
	  while (getline(stream, line)) {
	    std::vector<std::string> tokens;
	    NGT::Common::tokenize(line, tokens, " \t");
	    std::vector<float> object;
	    if (tokens.size() != 1) {
	      std::stringstream msg;
	      msg << "The specified codebook index is invalid. " << line;
	      NGTThrowException(msg);
	    }
	    codebookIndex.push_back(NGT::Common::strtol(tokens[0]));
	  }
	}
      }
      {
	std::string objectIndexPath = objectIndexFile;
	if (objectIndexPath.empty()) {
	  objectIndexPath = QBG::Index::getObjectIndexFileName(indexPath);
	}
	if (objectIndexPath != "-") {
	  std::ifstream stream(objectIndexPath);
	  if (!stream) {
	    std::stringstream msg;
	    msg << "Cannot open the codebook index. " << objectIndexPath;
	    NGTThrowException(msg);
	  }
	  std::string line;
	  while (getline(stream, line)) {
	    std::vector<std::string> tokens;
	    NGT::Common::tokenize(line, tokens, " \t");
	    std::vector<float> object;
	    if (tokens.size() != 1) {
	      std::stringstream msg;
	      msg << "The specified codebook index is invalid. " << line;
	      NGTThrowException(msg);
	    }
	    objectIndex.push_back(NGT::Common::strtol(tokens[0]));
	  }
        }
      }
      buildNGTQ(indexPath, quantizerCodebook, codebookIndex, objectIndex, beginID, endID, silence);
    }

    static void buildNGTQ(const std::string &indexPath,
		      std::vector<std::vector<float>> &quantizerCodebook,
		      std::vector<uint32_t> &codebookIndex,
		      std::vector<uint32_t> &objectIndex,
			  size_t beginID = 1, size_t endID = 0, bool silence = false) {
      NGT::StdOstreamRedirector redirector(silence);
      //redirector.begin();
      NGT::Timer timer;
      timer.start();
      NGTQ::Index index(indexPath);
      if (quantizerCodebook.size() == 0) {
	index.createIndex(beginID, endID);
      } else {
	if (codebookIndex.size() == 0) {
	  codebookIndex.resize(quantizerCodebook.size());
	}
	if (objectIndex.size() == 0) {
	  size_t size = index.getQuantizer().objectList.size();
	  size = size == 0 ? 0 : size - 1;
	  objectIndex.resize(size);
	}
	if ((quantizerCodebook.size() == 0) || (codebookIndex.size() == 0)) {
	  stringstream msg;
	  msg << "The specified codebooks or indexes are invalild " << quantizerCodebook.size() << ":" << codebookIndex.size();
	  NGTThrowException(msg);
	}
	index.createIndex(quantizerCodebook, codebookIndex, objectIndex, beginID, endID);
      }
      timer.stop();
      std::cerr << "NGTQ index is completed." << std::endl;
      std::cerr << "  time=" << timer << std::endl;
      std::cerr << "  vmsize=" << NGT::Common::getProcessVmSizeStr() << std::endl;
      std::cerr << "  peak vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;
      std::cerr << "saving..." << std::endl;
      index.save();
      redirector.end();
    }

    static void buildQBG(const std::string &indexPath, bool silence = false) {
      std::cerr << "build QBG" << std::endl;
      NGT::Timer timer;
      timer.start();
      QBG::Index index(indexPath, false, silence);
      timer.stop();
      std::cerr << "QBG index is completed." << std::endl;
      std::cerr << "  time=" << timer << std::endl;
      std::cerr << "  vmsize=" << NGT::Common::getProcessVmSizeStr() << std::endl;
      std::cerr << "  peak vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;
      std::cerr << "saving..." << std::endl;
      index.save();
    }

    void extract(std::ostream &os, size_t n, bool random = true) {
      if (n == 0) {
	NGTThrowException("QuantizedBlobGraph::extract # of objects is zero.");
      }
      auto &quantizer = getQuantizer();
      size_t dim = quantizer.property.dimension;
      std::vector<float> object;
      if (random) {
	struct timeval randTime;
	gettimeofday(&randTime, 0);
	srand(randTime.tv_usec);
	for (size_t cnt = 0; cnt < n; cnt++) {
	  double random = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
	  size_t id = 0;
	  while (true) {
	    do {
	      id = floor(quantizer.objectList.size() * random);
	    } while (id >= quantizer.objectList.size());
	    if (quantizer.objectList.get(id, object, &quantizer.globalCodebookIndex.getObjectSpace())) {
	      break;
	    } else {
	      std::cerr << "Cannot get the object. " << id << std::endl;
	    }
	  }
	  if (cnt + 1 % 100000 == 0) {
	    std::cerr << "loaded " << static_cast<float>(cnt + 1) / 1000000.0 << "M objects." << std::endl;
	  }
	  if (dim != 0) {
	    object.resize(dim, 0.0);
	  }
	  for (auto v = object.begin(); v != object.end(); ++v) {
	    if (v + 1 != object.end()) {
	      os << *v << "\t";
	    } else {
	      os << *v << std::endl;;
	    }
	  }
	}
      } else {
	for (size_t cnt = 1; cnt <= n; cnt++) {
	  if (!quantizer.objectList.get(cnt, object, &quantizer.globalCodebookIndex.getObjectSpace())) {
	    std::cerr << "Cannot get the object. " << cnt << std::endl;
	    continue;
	  }
	  if (cnt % 100000 == 0) {
	    std::cerr << "loaded " << static_cast<float>(cnt) / 1000000.0 << "M objects." << std::endl;
	  }
	  if (dim != 0) {
	    object.resize(dim, 0.0);
	  }
	  for (auto v = object.begin(); v != object.end(); ++v) {
	    if (v + 1 != object.end()) {
	      os << *v << "\t";
	    } else {
	      os << *v << std::endl;;
	    }
	  }
	  if (n > 0 && cnt >= n) {
	    break;
	  }
	}
      }
    }


    static void
      load(std::string indexPath, std::string blobs = "", std::string localCodebooks = "", std::string rotationPath = "", int threadSize = 0)
    {
      if (blobs.empty()) {
	blobs = QBG::Index::getBlobFileName(indexPath);
      }
      if (localCodebooks.empty()) {
	localCodebooks = QBG::Index::getPQFileName(indexPath) + "/opt-@.tsv";
      }
      if (rotationPath.empty()) {
	rotationPath = QBG::Index::getRotationFileName(indexPath);
      }

      threadSize = threadSize == 0 ? std::thread::hardware_concurrency() : threadSize;
      assert(threadSize != 0);

      size_t dataSize = 0;
      NGT::Index::append(indexPath + "/global", blobs, threadSize, dataSize);	

      std::cerr << "qbg: loading the local codebooks..." << std::endl;
      NGTQ::Property property;
      property.load(indexPath);

      std::vector<std::string> tokens;
      NGT::Common::tokenize(localCodebooks, tokens, "@");
      if (tokens.size() != 2) {
	NGTThrowException("No @ in the specified local codebook string.");
      }
      for (size_t no = 0; no < property.localDivisionNo; no++) {
	std::stringstream data;
	data << tokens[0] << no << tokens[1];
	std::stringstream localCodebook;
	localCodebook << indexPath << "/local-" << no;
	std::cerr << data.str() << "->" << localCodebook.str() << std::endl;
	NGT::Index::append(localCodebook.str(), data.str(), threadSize, dataSize);
      }

      cerr << "qbg: loading the rotation..." << endl;
#ifdef NGTQ_QBG
      std::vector<float> rotation;

      std::ifstream stream(rotationPath);
      if (!stream) {
	std::stringstream msg;
	msg << "Cannot open the rotation. " << rotationPath;
	NGTThrowException(msg);
      }
      std::string line;
      while (getline(stream, line)) {
	std::vector<std::string> tokens;
	NGT::Common::tokenize(line, tokens, " \t");
	for (auto &token : tokens) {
	  rotation.push_back(NGT::Common::strtof(token));
	}
      }
      std::cerr << "rotation matrix size=" << rotation.size() << std::endl;
      QBG::Index::load(indexPath, rotation);
#endif
    }

    static const std::string getTrainObjectFileName(std::string indexPath) { return indexPath + "/" + getWorkspaceName() + "/objects.tsv"; }
    static const std::string getPQFileName(std::string indexPath) { return indexPath + "/" + getWorkspaceName() + "/kmeans-cluster_opt"; }
    static const std::string getBlobFileName(std::string indexPath) { return indexPath + "/" + getWorkspaceName() + "/kmeans-cluster_centroid.tsv"; }

    static const std::string getQuantizerCodebookFileName(std::string indexPath) { return indexPath + "/" + getWorkspaceName() + "/kmeans-cluster_qcentroid.tsv"; }
    static const std::string getCodebookIndexFileName(std::string indexPath) { return indexPath + "/" + getWorkspaceName() + "/kmeans-cluster_bqindex.tsv"; }
    static const std::string getObjectIndexFileName(std::string indexPath) { return indexPath + "/" + getWorkspaceName() + "/kmeans-cluster_index.tsv"; }
    static const std::string getRotationFileName(std::string indexPath) { return getPQFileName(indexPath) + "/opt_R.tsv"; }

    static const std::string getWorkspaceName() { return "ws"; }

    const std::string path;
    
    QuantizedBlobGraphRepository quantizedBlobGraph;


  };

} 

#endif 
