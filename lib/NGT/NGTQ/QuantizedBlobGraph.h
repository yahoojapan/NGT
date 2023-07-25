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
#include	"NGT/NGTQ/Optimizer.h"
#include	"NGT/NGTQ/HierarchicalKmeans.h"

#include	<thread>


namespace QBG {

  class CreationParameters {
  public:
    CreationParameters() { setDefault(); }
    void setDefault() {
      numOfObjects			= 0;
      threadSize			= 24;
      numOfLocalClusters		= 16;
      dimension				= 0;
#ifdef NGTQ_QBG
      genuineDimension			= 0;
      dimensionOfSubvector		= 1;
      genuineDataType			= ObjectFile::DataTypeFloat;
#endif
      dataType				= NGTQ::DataTypeFloat;
      distanceType			= NGTQ::DistanceType::DistanceTypeL2;
      singleLocalCodebook		= false;
      numOfSubvectors			= 0;
      batchSize				= 1000;
      centroidCreationMode		= NGTQ::CentroidCreationModeStaticLayer;
      localCentroidCreationMode		= NGTQ::CentroidCreationModeStatic;
      localIDByteSize			= 1;
      localClusteringSampleCoefficient	= 10;
      objectListOnMemory		= false;
      
      globalEdgeSizeForCreation		= 10;
      globalEdgeSizeForSearch		= 40;
      globalIndexType			= NGT::Property::GraphAndTree;
      globalInsertionRadiusCoefficient	= 1.1;
      globalGraphType			= NGT::NeighborhoodGraph::GraphTypeANNG;

      localIndexType			= NGT::Property::GraphAndTree;
      localInsertionRadiusCoefficient	= 1.1;
      localGraphType			= NGT::NeighborhoodGraph::GraphTypeANNG;

      verbose				= false;
    }

    static void setProperties(CreationParameters &creation, NGTQ::Property &property, NGT::Property &globalProperty,
			      NGT::Property &localProperty) {
      property.threadSize		= creation.threadSize;
      property.globalCentroidLimit	= 0;
      property.localCentroidLimit	= creation.numOfLocalClusters;
      property.dimension		= creation.dimension;
      property.globalRange		= 0;
      property.localRange		= 0;
      property.localCentroidLimit	= creation.numOfLocalClusters;
#ifdef NGTQ_QBG
      property.genuineDimension		= creation.genuineDimension;
      //-/property.dimensionOfSubvector = creation.dimensionOfSubvector;
      property.genuineDataType		= creation.genuineDataType;
#endif
      property.dataType			= creation.dataType;
      property.distanceType		= creation.distanceType;
      property.singleLocalCodebook	= false;
      property.localDivisionNo		= creation.numOfSubvectors;
      property.batchSize		= creation.batchSize;
      property.centroidCreationMode	= creation.centroidCreationMode;
      property.localCentroidCreationMode = creation.localCentroidCreationMode;
      property.localIDByteSize		= creation.localIDByteSize;
      property.localClusteringSampleCoefficient = creation.localClusteringSampleCoefficient;
      property.objectListOnMemory	= creation.objectListOnMemory;
      globalProperty.edgeSizeForCreation = creation.globalEdgeSizeForCreation;
      globalProperty.edgeSizeForSearch	= creation.globalEdgeSizeForSearch;
      globalProperty.indexType		= creation.globalIndexType;
      globalProperty.insertionRadiusCoefficient = creation.globalInsertionRadiusCoefficient;
      globalProperty.graphType		= creation.globalGraphType;
      localProperty.indexType		= creation.localIndexType;
      localProperty.insertionRadiusCoefficient = creation.localInsertionRadiusCoefficient;
      localProperty.graphType		= creation.localGraphType;
      if (property.localCentroidLimit >= 0xFF) {
	if (property.localIDByteSize < 2) {
	  property.localIDByteSize = 2;
	}
      } else if (property.localCentroidLimit >= 0xFFFF) {
	property.localIDByteSize = 4;
      }
      property.dimension = property.dimension == 0 ? property.genuineDimension : property.dimension;
      property.localDivisionNo = property.localDivisionNo == 0 ? property.dimension : property.localDivisionNo;
    }

    size_t				numOfObjects;
    size_t				threadSize;
    size_t				numOfLocalClusters;
    size_t				dimension;
#ifdef NGTQ_QBG
    size_t				genuineDimension;
    size_t				dimensionOfSubvector;
    ObjectFile::DataType		genuineDataType;
#endif
    NGTQ::DataType			dataType;
    NGTQ::DistanceType			distanceType;
    bool				singleLocalCodebook;
    size_t				numOfSubvectors;
    size_t				batchSize;
    NGTQ::CentroidCreationMode		centroidCreationMode;
    NGTQ::CentroidCreationMode		localCentroidCreationMode;
    size_t				localIDByteSize;
    size_t				localClusteringSampleCoefficient;
    bool				objectListOnMemory;
    
    size_t				globalEdgeSizeForCreation;
    size_t				globalEdgeSizeForSearch;
    NGT::Property::IndexType		globalIndexType;
    float				globalInsertionRadiusCoefficient;
    NGT::Property::GraphType		globalGraphType;

    NGT::Property::IndexType		localIndexType;
    float				localInsertionRadiusCoefficient;
    NGT::Property::GraphType		localGraphType;

    bool				verbose;
  };

  class HierarchicalClusteringParameters {
  public:
    HierarchicalClusteringParameters() { setDefault(); }
    void setDefault() {
      maxSize			= 1000;
      numOfObjects		= 0;
      numOfClusters		= 2;
      numOfTotalClusters	= 0;
      numOfTotalBlobs		= 0;
      clusterID			= -1;
      initMode			= NGT::Clustering::InitializationModeKmeansPlusPlus;
      numOfRandomObjects	= 0;
      numOfFirstObjects		= 0;
      numOfFirstClusters	= 0;
      numOfSecondObjects	= 0;
      numOfSecondClusters	= 0;
      numOfThirdObjects		= 0;
      numOfThirdClusters	= 0;
      extractCentroid		= false;
      clusteringType		= QBG::HierarchicalKmeans::ClusteringTypeThreeLayer;
      epsilonExplorationSize	= 1000;
      expectedRecall		= 0.98;

      verbose			= false;
    }

    size_t		maxSize;
    size_t		numOfObjects;
    size_t		numOfClusters;
    size_t		numOfTotalClusters;
    size_t		numOfTotalBlobs;
    int32_t		clusterID;

    NGT::Clustering::InitializationMode initMode;

    size_t		numOfRandomObjects;

    size_t		numOfFirstObjects;
    size_t		numOfFirstClusters;
    size_t		numOfSecondObjects;
    size_t		numOfSecondClusters;
    size_t		numOfThirdObjects;
    size_t		numOfThirdClusters;
    bool		extractCentroid;

    QBG::HierarchicalKmeans::ClusteringType	clusteringType;
    size_t		epsilonExplorationSize;
    float		expectedRecall;

    bool		verbose;
  };

  class OptimizationParameters {
  public:
    OptimizationParameters() { setDefault(); }
    void setDefault() {
      clusteringType			= NGT::Clustering::ClusteringTypeKmeansWithoutNGT;
      initMode				= NGT::Clustering::InitializationModeHead;
      timelimit				= 24 * 1 * 60.0 * 60.0;
      iteration				= 1000;
      clusterIteration			= 400;
      clusterSizeConstraint		= false;
      clusterSizeConstraintCoefficient	= 10.0;
      convergenceLimitTimes		= 5;
      numOfObjects			= 1000;
      numOfClusters			= 0;
      numOfSubvectors			= 0;
      numOfMatrices			= 1;
      seedNumberOfSteps			= 2;
      seedStep				= 10;
      reject				= 0.9;
      repositioning			= false;
      rotation				= true;
      globalType			= QBG::Optimizer::GlobalTypeNone;
      randomizedObjectExtraction	= true;
      showClusterInfo			= false;

      verbose				= false;
    }
    NGT::Clustering::ClusteringType	clusteringType;
    NGT::Clustering::InitializationMode	initMode;

    float		timelimit;
    size_t		iteration;
    size_t		clusterIteration;		
    bool		clusterSizeConstraint;
    float		clusterSizeConstraintCoefficient;
    size_t		convergenceLimitTimes;		
    size_t		numOfObjects;
    size_t		numOfClusters;
    size_t		numOfSubvectors;
    size_t		numOfMatrices;
    size_t		seedNumberOfSteps;
    size_t		seedStep;
    float		reject;
    bool		repositioning;
    bool		rotation;
    QBG::Optimizer::GlobalType	globalType;
    bool		randomizedObjectExtraction;
    bool		showClusterInfo;

    bool		verbose;
  };

  class BuildParameters {
  public:
    BuildParameters(){ setDefault(); }

    void setDefault() {
      creation.setDefault();
      hierarchicalClustering.setDefault();
      optimization.setDefault();
    }

    void setProperties(NGTQ::Property &property, NGT::Property &globalProperty,
		      NGT::Property &localProperty) {
      CreationParameters::setProperties(creation, property, globalProperty, localProperty);
    }

    void setVerbose(bool s) {
      creation.verbose = s;
      hierarchicalClustering.verbose = s;
      optimization.verbose = s;
      verbose = s;
    }

    CreationParameters			creation;
    HierarchicalClusteringParameters	hierarchicalClustering;
    OptimizationParameters		optimization;

    bool				verbose;
  };


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
      
      (*this).resize(quantizedIndex.getInvertedIndexSize());
      NGT::Timer timer;
      timer.start();
      for (size_t gid = 1; gid < quantizedIndex.getInvertedIndexSize(); gid++) {
	if (gid % 100000 == 0) {
	  timer.stop();
	  std::cerr << "The number of processed blobs=" << gid << " VmSize=" <<  NGT::Common::getProcessVmSizeStr() << " Elapsed time=" << timer << std::endl;
	  timer.restart();
	}
	NGTQ::InvertedIndexEntry<uint16_t> invertedIndexObjects(numOfSubspaces);
	quantizedIndex.getQuantizer().extractInvertedIndexObject(invertedIndexObjects, gid);
	quantizedIndex.getQuantizer().eraseInvertedIndexObject(gid);
	NGTQ::QuantizedObjectProcessingStream quantizedStream(quantizedIndex.getQuantizer().divisionNo, invertedIndexObjects.size());
	rearrange(invertedIndexObjects, (*this)[gid], quantizedStream);
      }
#endif
    }

    static void rearrange(NGTQ::InvertedIndexEntry<uint16_t> &invertedIndexObjects, NGTQG::QuantizedNode &rearrangedObjects, NGTQ::QuantizedObjectProcessingStream &quantizedStream) {
	rearrangedObjects.ids.reserve(invertedIndexObjects.size());
	for (size_t oidx = 0; oidx < invertedIndexObjects.size(); oidx++) {
	  rearrangedObjects.ids.push_back(invertedIndexObjects[oidx].id);
	  for (size_t idx = 0; idx < invertedIndexObjects.numOfSubvectors; idx++) {
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

	rearrangedObjects.subspaceID = invertedIndexObjects.subspaceID;
	rearrangedObjects.objects = quantizedStream.compressIntoUint4();
    }

    static void rearrange(NGTQ::InvertedIndexEntry<uint16_t> &invertedIndexObjects, NGTQG::QuantizedNode &rearrangedObjects) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      std::cerr << "construct: Not implemented" << std::endl;
      abort();
#else
      if (invertedIndexObjects.numOfSubvectors == 0) {
	NGTThrowException("# of subvectors is zero.");
      }
      
      //(*this).resize(quantizedIndex.getInvertedIndexSize());
      NGT::Timer timer;
      timer.start();
      {
	//NGTQ::InvertedIndexEntry<uint16_t> invertedIndexObjects(numOfSubspaces);
	//quantizedIndex.getQuantizer().extractInvertedIndexObject(invertedIndexObjects, gid);
	//quantizedIndex.getQuantizer().eraseInvertedIndexObject(gid);
	NGTQ::QuantizedObjectProcessingStream quantizedStream(invertedIndexObjects.numOfSubvectors, invertedIndexObjects.size());

	rearrange(invertedIndexObjects, rearrangedObjects, quantizedStream);
      }
#endif
    }

    static void rearrange(NGTQ::QuantizedObjectSet &quantizedObjects, NGTQG::QuantizedNode &rearrangedObjects) {
      NGTQ::InvertedIndexEntry<uint16_t> iie;
      iie.set(quantizedObjects);
      rearrange(iie, rearrangedObjects);
    }
  };
  
  class Index : public NGTQ::Index {
  public:
  Index(const std::string &indexPath, bool readOnly = false, bool verbose = false) :
    NGTQ::Index(indexPath, readOnly), path(indexPath), quantizedBlobGraph(*this) {
      searchable = false;
      NGT::StdOstreamRedirector redirector(!verbose);
      redirector.begin();
      try {
	load();
	searchable = true;
      } catch (NGT::Exception &err) {
	if (readOnly) {
	  stringstream msg;
	  msg << "QBG::Index: No quantized blob graph. " << err.what();
	  NGTThrowException(msg);
	} else {
	}
      }
      redirector.end();
    }

    ~Index() {}

    bool &getVerbose() { return verbose; }

#ifdef NGTQ_QBG
    static void create(const std::string &index,
		       BuildParameters &buildParameters,
		       std::vector<float> *rotation = 0,const std::string objectFile = "") {
      create(index, buildParameters.creation, rotation, objectFile);
    }
    static void create(const std::string &index,
		       CreationParameters &creation,
		       std::vector<float> *rotation = 0,const std::string objectFile = "") {
      NGTQ::Property property;
      NGT::Property globalProperty;
      NGT::Property localProperty;
      CreationParameters::setProperties(creation, property, globalProperty, localProperty);
      property.quantizerType = NGTQ::QuantizerTypeQBG;
      NGTQ::Index::create(index, property, globalProperty, localProperty, rotation, objectFile);
    }
#endif
#ifdef NGTQ_QBG
    static void initialize(NGTQ::Property &property, NGT::Property &globalProperty,NGT::Property &localProperty) {
      QBG::CreationParameters params;
      QBG::CreationParameters::setProperties(params, property, globalProperty, localProperty);
    }
#endif
    
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
#else
      NGTQ::Index::create(index, property, globalProperty, localProperty);
#endif
    }

    static void load(const std::string &indexPath, const std::vector<std::vector<float>> &quantizerCodebook, const std::vector<float> &rotation) {
      NGTQ::Index index(indexPath);
      index.getQuantizer().loadQuantizationCodebookAndRotation(quantizerCodebook, rotation);
    }

    void insert(const size_t id, std::vector<float> &object) {
      getQuantizer().objectList.put(id, object, &getQuantizer().globalCodebookIndex.getObjectSpace());
    }
    template<typename T>
    NGT::ObjectID append(std::vector<T> &object) {
      NGT::ObjectID id = getQuantizer().objectList.size();
      id = id == 0 ? 1 : id;
      if (typeid(T) == typeid(float)) {
	auto &obj = *reinterpret_cast<std::vector<float>*>(&object);
	getQuantizer().objectList.put(id, obj, &getQuantizer().globalCodebookIndex.getObjectSpace());
      } else {
	std::vector<float> obj(object.begin(), object.end());
	getQuantizer().objectList.put(id, obj, &getQuantizer().globalCodebookIndex.getObjectSpace());
      }
      return id;
    }

    static void append(const std::string &indexName,	// index file
		       const std::string &data,	// data file
		       size_t dataSize = 0,	// data size
		       bool verbose = false
		       ) {
      NGT::StdOstreamRedirector redirector(!verbose);
      redirector.begin();
      QBG::Index index(indexName);
      auto &quantizer = index.getQuantizer();
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
      size_t idx = quantizer.objectList.size() == 0 ? 0 : quantizer.objectList.size() - 1;
      size_t count = 0;
      // extract objects from the file and insert them to the object list.
      while(getline(*is, line)) {
	idx++;
	count++;
	std::vector<float>	object;
	NGT::Common::extractVector(line, " ,\t", object);
	if (object.empty()) {
	  cerr << "An empty line or invalid value: " << line << endl;
	  continue;
	}
	index.insert(idx, object);

	if (count % 100000 == 0) {
	  std::cerr << "appended " << static_cast<float>(count) / 1000000.0 << "M objects.";
	  if (count != idx) {
	    std::cerr << " # of the total objects=" << static_cast<float>(idx) / 1000000.0 << "M";
	  }
	  cerr << " virtual memory(kbyte)=" << NGT::Common::getProcessVmSize() << std::endl;
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
			     bool verbose = false
		       ) {
      NGT::StdOstreamRedirector redirector(!verbose);
      redirector.begin();
      QBG::Index index(indexName);
      auto &quantizer = index.getQuantizer();
      std::vector<std::string> tokens;
      NGT::Common::tokenize(data, tokens, ".");
      if (tokens.size() < 2) {
	std::stringstream msg;
	msg << "Invalid file name format";
	NGTThrowException(msg);
      }
      StaticObjectFileLoader loader(data, tokens[tokens.size() - 1]);
      size_t idx = quantizer.objectList.size() == 0 ? 0 : quantizer.objectList.size() - 1;
      size_t count = 0;
      while (!loader.isEmpty()) {
	idx++;
	count++;
	if (dataSize > 0 && idx > dataSize) {
	  break;
	}
	auto object = loader.getObject();
	index.insert(idx, object);
	if (count % 1000000 == 0) {
	  std::cerr << "appended " << static_cast<float>(count) / 1000000.0 << "M objects.";
	  if (count != idx) {
	    std::cerr << " # of the total objects=" << static_cast<float>(idx) / 1000000.0 << "M";
	  }
	  cerr << " virtual memory(kbyte)=" << NGT::Common::getProcessVmSize() << std::endl;
	}
      }
      index.save();
      index.close();
      redirector.end();
    }

    float getApproximateDistances(std::vector<float> &query, NGTQG::RearrangedQuantizedObjectSet &quantizedObjects,
				 size_t subspaceID, std::vector<float> &distances) {
      if (query.empty()) {
	NGTThrowException("The specified query is empty.");
      }
      auto &quantizer = this->getQuantizer();
      if (quantizer.getNumOfLocalClusters() != 16) {
	std::stringstream msg;
	msg << "# of the local clusters is not 16. " << quantizer.getNumOfLocalClusters() << std::endl;
	NGTThrowException(msg);
      }
      distances.clear();
      auto noOfObjects = quantizedObjects.ids.size();
      if (noOfObjects == 0) {
	return 0.0;
      }
      auto rotatedQuery = query;
      auto &quantizedObjectDistance = quantizer.getQuantizedObjectDistance();
      quantizedObjectDistance.rotation->mul(rotatedQuery.data());
      NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8 lookupTable;
      quantizedObjectDistance.initialize(lookupTable);
      quantizedObjectDistance.createDistanceLookup(rotatedQuery.data(), subspaceID, lookupTable);
      distances.resize(NGTQ::QuantizedObjectProcessingStream::getNumOfAlignedObjects(noOfObjects));
      auto minDistance = quantizedObjectDistance(quantizedObjects.objects, distances.data(), noOfObjects, lookupTable);
      distances.resize(noOfObjects);
      return minDistance;
    }

    void getApproximateDistances(std::vector<float> &query, NGTQ::QuantizedObjectSet &quantizedObjects,
				 size_t subspaceID, std::vector<float> &distances) {
      if (query.empty()) {
	NGTThrowException("The specified query is empty.");
      }
      auto &quantizer = this->getQuantizer();
      distances.clear();
      auto noOfObjects = quantizedObjects.size();
      if (noOfObjects == 0) {
	return;
      }
      auto rotatedQuery = query;
      auto &quantizedObjectDistance = quantizer.getQuantizedObjectDistance();
      quantizedObjectDistance.rotation->mul(rotatedQuery.data());
      NGTQ::QuantizedObjectDistance::DistanceLookupTable lookupTable;
      quantizedObjectDistance.initialize(lookupTable);
      quantizedObjectDistance.createDistanceLookup(rotatedQuery.data(), subspaceID, lookupTable);
      distances.resize(noOfObjects);
      if (quantizer.localIDByteSize == 1) {
	NGTQ::InvertedIndexEntry<uint8_t> iie;
	iie.set(quantizedObjects);
	for (size_t idx = 0; idx < iie.size(); idx++) {
	  distances[idx] = quantizedObjectDistance(&iie[idx].localID[0], lookupTable);
	}
      } else if (quantizer.localIDByteSize == 2) {
	NGTQ::InvertedIndexEntry<uint16_t> iie;
	iie.set(quantizedObjects);
	for (size_t idx = 0; idx < iie.size(); idx++) {
	  distances[idx] = quantizedObjectDistance(&iie[idx].localID[0], lookupTable);
	}
      } else if (quantizer.localIDByteSize == 4) {
	NGTQ::InvertedIndexEntry<uint32_t> iie;
	iie.set(quantizedObjects);
	for (size_t idx = 0; idx < iie.size(); idx++) {
	  distances[idx] = quantizedObjectDistance(&iie[idx].localID[0], lookupTable);
	}
      }
    }

    static void appendFromObjectRepository(const std::string &ngtIndex,	// QG
					   const std::string &qgIndex,	// NGT
					   bool verbose = false) {
      NGT::StdOstreamRedirector redirector(!verbose);
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

    static void refineDistances(QBG::SearchContainer &searchContainer, NGTQ::Quantizer &quantizer,
				NGT::NeighborhoodGraph::ResultSet &result,
				NGT::ObjectDistances &qresults) {
      auto &objectSpace = quantizer.globalCodebookIndex.getObjectSpace();
      NGT::ResultPriorityQueue qres;
      if (objectSpace.getObjectType() == typeid(float)) {
	refineDistances<float>(searchContainer, quantizer, result, qres);
      } else if (objectSpace.getObjectType() == typeid(uint8_t)) {
	refineDistances<uint8_t>(searchContainer, quantizer, result, qres);
      } else if (objectSpace.getObjectType() == typeid(NGT::float16)) {
	refineDistances<NGT::float16>(searchContainer, quantizer, result, qres);
      } else {
	std::stringstream msg;
	msg << "refineDistances: Fatal error! Invalid datatype. " << objectSpace.getObjectType().name() << std::endl;
	NGTThrowException(msg);
      }
      qresults.resize(qres.size());
      for (int i = qresults.size() - 1; i >= 0; i--) {
	qresults[i] = qres.top();
	qres.pop();
      }
    }

    static void refineDistances(QBG::SearchContainer &searchContainer, NGTQ::Quantizer &quantizer,
				NGT::NeighborhoodGraph::ResultSet &result,
				NGT::ResultPriorityQueue &qresults) {
      auto &objectSpace = quantizer.globalCodebookIndex.getObjectSpace();
      if (objectSpace.getObjectType() == typeid(float)) {
	refineDistances<float>(searchContainer, quantizer, result, qresults);
      } else if (objectSpace.getObjectType() == typeid(uint8_t)) {
	refineDistances<uint8_t>(searchContainer, quantizer, result, qresults);
      } else if (objectSpace.getObjectType() == typeid(NGT::float16)) {
	refineDistances<NGT::float16>(searchContainer, quantizer, result, qresults);
      } else {
	std::stringstream msg;
	msg << "refineDistances: Fatal error! Invalid datatype. " << objectSpace.getObjectType().name() << std::endl;
	NGTThrowException(msg);
      }
    }

    template<typename T>
    static void refineDistances(QBG::SearchContainer &searchContainer, NGTQ::Quantizer &quantizer,
				NGT::NeighborhoodGraph::ResultSet &result,
				NGT::ResultPriorityQueue &qresults) {
      qresults = NGT::ResultPriorityQueue();
      NGT::Object &query = searchContainer.object;
      auto &objectSpace = quantizer.globalCodebookIndex.getObjectSpace();
      auto paddedDimension = objectSpace.getPaddedDimension();
      const size_t prefetchSize = objectSpace.getPrefetchSize();
#ifdef NGTQ_OBJECT_IN_MEMORY
      if (quantizer.objectListOnMemory.size() != 0) {
	while (!result.empty()) {
	  auto r = result.top();
	  result.pop();
	  NGT::Object &object = *quantizer.objectListOnMemory.get(r.id);
	  if (!result.empty()) {
            uint8_t *ptr = static_cast<uint8_t*>(quantizer.objectListOnMemory.get(result.top().id)->getPointer());
	    NGT::MemoryCache::prefetch(ptr, prefetchSize);
	  }
	  r.distance = objectSpace.getComparator()(query, object);
	  qresults.push(r);
	}
      } else {
#endif
	auto threadid = omp_get_thread_num();
	while (!result.empty()) {
	  auto r = result.top();
	  result.pop();
	  std::vector<T> object;
#ifdef MULTIPLE_OBJECT_LISTS
	  quantizer.objectList.get(threadid, r.id, object);
#else
	  quantizer.objectList.get(r.id, object);
#endif
	  r.distance = NGT::PrimitiveComparator::compareL2(static_cast<T*>(query.getPointer()),
							   static_cast<T*>(object.data()), paddedDimension);


	  qresults.push(r);
	}
#ifdef NGTQ_OBJECT_IN_MEMORY
      }
#endif
      while (qresults.size() > searchContainer.exactResultSize) {
	qresults.pop();
      }

    }

    void searchInTwoSteps(QBG::SearchContainer &searchContainer) {
      if (searchContainer.isEmptyObject()) {
	NGT::Object query(searchContainer.objectVector, getQuantizer().globalCodebookIndex.getObjectSpace());
	SearchContainer sc(searchContainer, query);
	searchInTwoSteps(sc);
	searchContainer.workingResult = std::move(sc.workingResult);
	return;
      }
      NGT::ObjectDistances blobs;
      NGT::SearchContainer sc(searchContainer);
      sc.setResults(&blobs);
      sc.setEpsilon(searchContainer.blobExplorationCoefficient - 1.0);
      sc.setSize(searchContainer.numOfProbes);

      auto &quantizer = getQuantizer();
      auto &globalIndex = quantizer.globalCodebookIndex;

      auto &quantizedObjectDistance = quantizer.getQuantizedObjectDistance();
      if (searchContainer.objectVector.size() == 0) {
	NGTThrowException("search: object is null.");
      }
      std::vector<float> rotatedQuery = searchContainer.objectVector;
      {
	NGT::Object *query = allocateObject(searchContainer.objectVector);
	NGT::SearchContainer tsc(sc, *query);
	tsc.setResults(&sc.getResult());
	globalIndex.search(tsc);
	globalIndex.deleteObject(query);
      }
      if (blobs.empty()) {
	std::cerr << "Warning: No blobs can be searched." << std::endl;
	std::cerr << "  global index size=" << globalIndex.getObjectRepositorySize() << std::endl;
	std::cerr << "  size=" << sc.size << std::endl;
	return;
      }

#if defined(NGTQG_ROTATION)
      if (quantizedObjectDistance.rotation != 0) {
	quantizedObjectDistance.rotation->mul(rotatedQuery.data());
      }
#endif
      std::unordered_map<size_t, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8> luts;
      size_t foundCount = 0;
      size_t k = searchContainer.size;
      NGT::Distance radius = FLT_MAX;
      NGT::NeighborhoodGraph::ResultSet result;
#ifdef NGTQBG_COARSE_BLOB
      NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8 lookupTable;
      quantizedObjectDistance.initialize(lookupTable);
#endif
      for (size_t idx = 0; idx < blobs.size(); idx++) {
#ifdef NGTQBG_COARSE_BLOB
	NGT::Distance blobDistance = std::numeric_limits<NGT::Distance>::max();
	auto graphNodeID = blobs[idx].id;
	auto &graphNodeToInvertedIndexEntries = quantizer.getGraphNodeToInvertedIndexEntries();
	auto beginIvtID = graphNodeToInvertedIndexEntries[graphNodeID - 1] + 1;
	auto endIvtID = graphNodeToInvertedIndexEntries[graphNodeID] + 1;
	for (auto blobID = beginIvtID; blobID < endIvtID; blobID++) {
	  auto subspaceID = quantizedBlobGraph[blobID].subspaceID;
	  quantizedObjectDistance.createDistanceLookup(rotatedQuery, subspaceID, lookupTable);
	  NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8 &lut = lookupTable;
#else
	{
	  auto blobID = blobs[idx].id;
	  auto subspaceID = quantizedBlobGraph[blobID].subspaceID;
	  auto luti = luts.find(subspaceID);
	  if (luti == luts.end()) {
	    luts.insert(std::make_pair(subspaceID, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8()));
	    luti = luts.find(subspaceID);
	    quantizedObjectDistance.initialize((*luti).second);
	    quantizedObjectDistance.createDistanceLookup(rotatedQuery.data(), subspaceID, (*luti).second);
          }
	  NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8 &lut = (*luti).second;
#endif

	  NGT::Distance bd;
	  std::tie(bd, radius) = judge(quantizedBlobGraph[blobID], k, radius, lut, result, foundCount);
#ifdef NGTQBG_COARSE_BLOB
	  if (bd < blobDistance) {
	    blobDistance = bd;
	  }
#else
#endif
        }
#ifdef NGTQBG_MIN
#endif
      }
      if (searchContainer.resultIsAvailable()) {
	if (searchContainer.exactResultSize > 0) {
	  NGT::ObjectDistances &qresults = searchContainer.getResult();
	  refineDistances(searchContainer, quantizer, result, qresults);
	} else {
	  searchContainer.getResult().moveFrom(result);
	}
      } else {
	if (searchContainer.exactResultSize > 0) {
	  refineDistances(searchContainer, quantizer, result, searchContainer.workingResult);
	} else {
	  searchContainer.workingResult = std::move(result);
	}
      }
    }

    void searchInOneStep(QBG::SearchContainer &searchContainer) {
      auto &globalIndex = getQuantizer().globalCodebookIndex;
      auto &globalGraph = static_cast<NGT::GraphAndTreeIndex&>(globalIndex.getIndex());
      NGT::ObjectDistances	seeds;
      const size_t dimension = globalIndex.getObjectSpace().getPaddedDimension();
      if (dimension > searchContainer.objectVector.size()) {
	searchContainer.objectVector.resize(dimension);
      }
      NGT::Object query(searchContainer.objectVector, globalIndex.getObjectSpace());
      SearchContainer sc(searchContainer, query);
      globalGraph.getSeedsFromTree(sc, seeds);
      if (seeds.empty()) {
	globalGraph.getRandomSeeds(globalGraph.repository, seeds, 20);
      }
      searchInOneStep(sc, seeds);
      searchContainer.workingResult = std::move(sc.workingResult);
    }

    void searchInOneStep(QBG::SearchContainer &searchContainer, NGT::ObjectDistances &seeds) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      std::cerr << "searchBlobGraph: Not implemented. " << std::endl;
      abort();
#else
      if (!searchable) {
	std::stringstream msg;
	msg << "The specified index is not now searchable. ";
	NGTThrowException(msg);
      }

      auto &quantizer = getQuantizer();
      auto &globalIndex = quantizer.globalCodebookIndex;
      auto &globalGraph = static_cast<NGT::GraphAndTreeIndex&>(globalIndex.getIndex());
      auto &objectSpace = globalIndex.getObjectSpace();

      if (globalGraph.searchRepository.empty()) {
	NGTThrowException("QBG:Index: graph repository is empty.");
      }
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
      } else if (objectSpace.getObjectType() == typeid(uint8_t)) {
	globalGraph.setupDistances(searchContainer, seeds, NGT::PrimitiveComparator::L2Uint8::compare);
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
      std::vector<float> rotatedQuery = searchContainer.objectVector;
      quantizedObjectDistance.rotation->mul(rotatedQuery.data());
      NGT::Distance radius = searchContainer.radius;
      if (requestedSize >= std::numeric_limits<int32_t>::max()) {
	radius *= searchContainer.explorationCoefficient;
      }
      NGT::ReadOnlyGraphNode *nodes = globalGraph.searchRepository.data();
      NGT::ReadOnlyGraphNode *neighbors = 0;
      NGT::ObjectDistance target;
      const size_t prefetchSize = objectSpace.getPrefetchSize();
      const size_t prefetchOffset = objectSpace.getPrefetchOffset();
      pair<uint64_t, NGT::PersistentObject*> *neighborptr;
      pair<uint64_t, NGT::PersistentObject*> *neighborendptr;
#ifdef NGTQBG_COARSE_BLOB
      NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8 lookupTable;
      quantizedObjectDistance.initialize(lookupTable);
#endif
      for (;;) {
	if (untracedNodes.empty() || untracedNodes.top().distance > explorationRadius) {
	  explorationSize++;
	  NGT::Distance blobDistance = std::numeric_limits<NGT::Distance>::max();
#ifdef NGTQBG_COARSE_BLOB
	  auto graphNodeID = currentNearestBlob.id;
	  auto &graphNodeToInvertedIndexEntries = quantizer.getGraphNodeToInvertedIndexEntries();
	  auto beginIvtID = graphNodeToInvertedIndexEntries[graphNodeID - 1] + 1;
	  auto endIvtID = graphNodeToInvertedIndexEntries[graphNodeID] + 1;
	  for (auto blobID = beginIvtID; blobID < endIvtID; blobID++) {
	    auto subspaceID = quantizedBlobGraph[blobID].subspaceID;
	    quantizedObjectDistance.createDistanceLookup(rotatedQuery.data(), subspaceID, lookupTable);
	    NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8 &lut = lookupTable;
#else
	  {
	    auto blobID = currentNearestBlob.id;
	    auto subspaceID = quantizedBlobGraph[blobID].subspaceID;
	    auto luti = luts.find(subspaceID);
	    if (luti == luts.end()) {
	      luts.insert(std::make_pair(subspaceID, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8()));
	      luti = luts.find(subspaceID);
	      quantizedObjectDistance.initialize((*luti).second);
	      quantizedObjectDistance.createDistanceLookup(rotatedQuery.data(), subspaceID, (*luti).second);
	    }
	    NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8 &lut = (*luti).second;
#endif
	    size_t foundCount;
	    NGT::Distance bd;
	    std::tie(bd, radius) = judge(quantizedBlobGraph[blobID], requestedSize,
					 radius, lut, results, foundCount);
#ifdef NGTQBG_COARSE_BLOB
	    if (bd < blobDistance) {
	      blobDistance = bd;
	    }
#else
	    blobDistance = bd;
#endif
	  }

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
	  NGT::Distance distance = objectSpace.getComparator()(searchContainer.object, *neighborptr->second);
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
	  refineDistances(searchContainer, quantizer, results, qresults);
	} else {
	  searchContainer.getResult().moveFrom(results);
	}
      } else {
	if (searchContainer.exactResultSize > 0) {
	  refineDistances(searchContainer, quantizer, results, searchContainer.workingResult);
	} else {
	  searchContainer.workingResult = std::move(results);
	}
      }
#endif
    }

    void search(QBG::SearchContainer &searchContainer) {
      searchInOneStep(searchContainer);
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

    static void buildNGTQ(const std::string &indexPath, bool verbose = false) {
      load(indexPath, QBG::Index::getQuantizerCodebookFile(indexPath), "", "", "", verbose);
      buildNGTQ(indexPath, "", "-", "-", 1, 0, verbose);
      if (verbose) {
	std::cerr << "NGTQ and NGTQBG indices are completed." << std::endl;
	std::cerr << "  vmsize=" << NGT::Common::getProcessVmSizeStr() << std::endl;
	std::cerr << "  peak vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;
      }
    }

    static void build(const std::string &indexPath, bool verbose = false) {
      load(indexPath, "", "", "");
      buildNGTQ(indexPath, "", "", "", 1, 0, verbose);
      buildQBG(indexPath, verbose);
      if (verbose) {
	std::cerr << "NGTQ and NGTQBG indices are completed." << std::endl;
	std::cerr << "  vmsize=" << NGT::Common::getProcessVmSizeStr() << std::endl;
	std::cerr << "  peak vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;
      }
    }

    static void build(const std::string &indexPath,
		      std::string quantizerCodebookFile = "",
		      std::string codebookIndexFile = "",
		      std::string objectIndexFile = "",
		      size_t beginID = 1, size_t endID = 0, bool verbose = false) {
      buildNGTQ(indexPath, quantizerCodebookFile, codebookIndexFile, objectIndexFile, beginID, endID, verbose);
      buildQBG(indexPath, verbose);
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

    static void buildNGTQ(const std::string &indexPath,
			  std::string quantizerCodebookFile = "",
			  std::string codebookIndexFile = "",
			  std::string objectIndexFile = "",
			  size_t beginID = 1, size_t endID = 0, bool verbose = false) {
      std::vector<std::vector<float>> quantizerCodebook;
      std::vector<uint32_t> codebookIndex;
      std::vector<uint32_t> objectIndex;
      {
	std::string codebookPath = quantizerCodebookFile;
	if (codebookPath.empty()) {
	  codebookPath = QBG::Index::getQuantizerCodebookFile(indexPath);
	}
	if (codebookPath != "-") {
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
      }
      {
        std::string codebookIndexPath = codebookIndexFile;
	if (codebookIndexPath.empty()) {
	  codebookIndexPath = QBG::Index::getCodebookIndexFile(indexPath);
	}
	if (codebookIndexPath != "-") {
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
	  objectIndexPath = QBG::Index::getObjectIndexFile(indexPath);
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
	      msg << "The specified object index is invalid. " << line;
	      NGTThrowException(msg);
	    }
	    objectIndex.push_back(NGT::Common::strtol(tokens[0]));
	  }
        }
      }
      buildNGTQ(indexPath, quantizerCodebook, codebookIndex, objectIndex, beginID, endID, verbose);
    }

    static void buildNGTQ(const std::string &indexPath,
			  std::vector<std::vector<float>> &quantizerCodebook,
			  std::vector<uint32_t> &codebookIndex,
			  std::vector<uint32_t> &objectIndex,
			  size_t beginID = 1, size_t endID = 0, bool verbose = false) {
      NGT::StdOstreamRedirector redirector(!verbose);
      redirector.begin();
      NGT::Timer timer;
      timer.start();
      NGTQ::Index index(indexPath);
      if ((quantizerCodebook.size() == 0) && (codebookIndex.size() == 0) && (objectIndex.size() == 0)) {
	index.createIndex(beginID, endID);
      } else {
	if (codebookIndex.size() == 0) {
	  codebookIndex.resize(quantizerCodebook.size());
	}
	if (codebookIndex.size() == 0) {
	  stringstream msg;
	  msg << "The specified codebook indexe invalild " << codebookIndex.size();
	  NGTThrowException(msg);
	}
	if (objectIndex.size() == 0) {
	  size_t size = index.getQuantizer().objectList.size();
	  size = size == 0 ? 0 : size - 1;
	  objectIndex.resize(size);
	}
	index.createIndex(quantizerCodebook, codebookIndex, objectIndex, beginID, endID);
      }

      {
	char *s = getenv("NGT_NOT_REMOVE_WORKSPACE");
	if (s == 0) {
	  const string comrmdir = "rm -rf " + indexPath + "/" + getWorkspaceName();
	  if (system(comrmdir.c_str()) == -1) {
	    std::cerr << "Warning. cannot remove the workspace directory. "
		      << comrmdir << std::endl;
	  }
	}
	const string comrm = "rm -f " + indexPath + "/" + NGTQ::Quantizer::getInvertedIndexFile();
	if (system(comrm.c_str()) == -1) {
	  std::cerr << "Warning. cannot remove the indeverted index. "
		    << comrm << std::endl;
	}
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

    static void buildQBG(const std::string &indexPath, bool verbose = false) {
      NGT::Timer timer;
      timer.start();
      auto readOnly = false;
      QBG::Index index(indexPath, readOnly, verbose);
      try {
	index.load();
	stringstream msg;
	msg << "QBG::Index::buildQBG: The index is already built. ";
	NGTThrowException(msg);
      } catch (...) {}
      index.quantizedBlobGraph.construct(index);

      timer.stop();
      if (verbose) {
	std::cerr << "QBG index is completed." << std::endl;
	std::cerr << "  time=" << timer << std::endl;
	std::cerr << "  vmsize=" << NGT::Common::getProcessVmSizeStr() << std::endl;
	std::cerr << "  peak vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;
	std::cerr << "saving..." << std::endl;
      }
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
	if (n > quantizer.objectList.size() / 2) {
	  if (n > quantizer.objectList.size() - 1) {
	    n = quantizer.objectList.size() - 1;
	  }
	  size_t pickedObjectCount = 0;
	  for (size_t id = 1; id < quantizer.objectList.size(); id++) {
	    double random = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
	    double p = static_cast<double>(n - pickedObjectCount) / static_cast<double>(quantizer.objectList.size() - id);
	    if (p == 0.0) {
	      break;
	    }
	    if (random <= p) {
	      if (!quantizer.objectList.get(id, object, &quantizer.globalCodebookIndex.getObjectSpace())) {
		std::cerr << "Cannot get the object. " << id << std::endl;
		continue;
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
	      pickedObjectCount++;
	      if (pickedObjectCount == n) {
		break;
	      }
	      if (pickedObjectCount % 100000 == 0) {
		std::cerr << "loaded " << static_cast<float>(pickedObjectCount + 1) / 1000000.0 << "M objects." << std::endl;
	      }
	    }
	  }
	} else {
	  std::unordered_set<uint32_t> pickedObjects;
	  for (size_t cnt = 0; cnt < n; cnt++) {
	    size_t id = 0;
	    while (true) {
	      do {
		double random = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
		id = floor(quantizer.objectList.size() * random);
	      } while (pickedObjects.count(id) > 0 || id >= quantizer.objectList.size());
	      if (quantizer.objectList.get(id, object, &quantizer.globalCodebookIndex.getObjectSpace())) {
		pickedObjects.insert(id);
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
    load(std::string indexPath, std::string blobs = "", std::string localCodebooks = "", std::string quantizerCodebook = "", std::string rotationPath = "", bool verbose = false, int threadSize = 0)
    {
      NGT::StdOstreamRedirector redirector(!verbose);
      redirector.begin();
      if (blobs.empty()) {
	blobs = QBG::Index::getBlobFile(indexPath);
      }
      if (localCodebooks.empty()) {
	localCodebooks = QBG::Index::getPQFile(indexPath) + "/" + QBG::Index::getSubvectorPrefix() + "-@";
      }
      if (quantizerCodebook.empty()) {
	quantizerCodebook = QBG::Index::getQuantizerCodebookFile(indexPath);
      }
      if (rotationPath.empty()) {
	rotationPath = QBG::Index::getRotationFile(indexPath);
      }

      threadSize = threadSize == 0 ? std::thread::hardware_concurrency() : threadSize;
      assert(threadSize != 0);

      size_t dataSize = 0;
      {
	const char *ngtDirString = "/tmp/ngt-XXXXXX";
	char ngtDir[strlen(ngtDirString) + 1];
	strcpy(ngtDir, ngtDirString);
	std::string tmpDir = mkdtemp(ngtDir);
	const std::string mvcom = "mv " + indexPath + "/" + NGTQ::Quantizer::getGlobalFile()
	  + " " + tmpDir + "/";
	if (system(mvcom.c_str()) == -1) {
	  std::stringstream msg;
	  msg << "Error! moving is failed. " << mvcom;
	  NGTThrowException(msg);
	}

	NGT::Index::append(tmpDir + "/" + NGTQ::Quantizer::getGlobalFile(), blobs, threadSize, dataSize);

	auto unlog = false;
	NGT::GraphOptimizer graphOptimizer(unlog);
	graphOptimizer.searchParameterOptimization = false;
	graphOptimizer.prefetchParameterOptimization = false;
	graphOptimizer.accuracyTableGeneration = false;
	int numOfOutgoingEdges = 10;
	int numOfIncomingEdges = 120;
	int numOfQueries = 200;
	int numOfResultantObjects = 20;
	graphOptimizer.set(numOfOutgoingEdges, numOfIncomingEdges, numOfQueries, numOfResultantObjects);
	graphOptimizer.execute(tmpDir + "/" + NGTQ::Quantizer::getGlobalFile(), indexPath + "/global");
	const std::string rmcom = "rm -rf " + tmpDir;
	if (system(rmcom.c_str()) == -1) {
	  std::cerr << "Warning. remove is failed. " << rmcom << std::endl;
	}
      }
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
	localCodebook << indexPath << "/" + NGTQ::Quantizer::getLocalPrefix() << no;
	std::cerr << data.str() << "->" << localCodebook.str() << std::endl;
	NGT::Index::append(localCodebook.str(), data.str(), threadSize, dataSize);
      }

#ifdef NGTQ_QBG
      std::vector<std::vector<float>> qCodebook;
      {
	std::ifstream stream(quantizerCodebook);
	if (!stream) {
	  std::stringstream msg;
	  msg << "Cannot open the codebook. " << quantizerCodebook;
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
	  if (!qCodebook.empty() && qCodebook[0].size() != object.size()) {
	    std::stringstream msg;
	    msg << "The specified quantizer codebook is invalid. " << qCodebook[0].size()
		<< ":" << object.size() << ":" << qCodebook.size() << ":" << line;
	    NGTThrowException(msg);
	  }
	  if (!object.empty()) {
	    qCodebook.push_back(object);
	  }
	}
      }
      {
	cerr << "qbg: loading the rotation..." << endl;
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
	QBG::Index::load(indexPath, qCodebook, rotation);
      }
#endif
      redirector.end();
    }

    static const std::string getSubvectorPrefix() { return "sv"; }
    static const std::string getHierarchicalClusteringPrefix() { return "hkc"; }
    static const std::string getSecondCentroidSuffix() { return "_2c"; }
    static const std::string getThirdCentroidSuffix() { return "_3c"; }
    static const std::string get3rdTo2ndSuffix() { return "_3to2"; }
    static const std::string getObjTo3rdSuffix() { return "_oto3"; }
    static const std::string getResidualFile() { return "r"; }
    static const std::string getRotatedResidualFile() { return "Rr"; }
    static const std::string getObjectFile() { return "obj"; }
    static const std::string getRotationFile() { return "R"; }
    static const std::string getWorkSpacePrefix(std::string indexPath) { return indexPath + "/" + getWorkspaceName(); }
    static const std::string getTrainObjectFile(std::string indexPath) { return getWorkSpacePrefix(indexPath) + "/" + getObjectFile(); }
    static const std::string getPrefix(std::string indexPath) { return getWorkSpacePrefix(indexPath) + "/" + getHierarchicalClusteringPrefix(); }
    static const std::string getPQFile(std::string indexPath) { return getPrefix(indexPath) + "_opt"; }
#ifdef NGTQBG_COARSE_BLOB
    static const std::string getBlobFile(std::string indexPath) { return getPrefix(indexPath) + getSecondCentroidSuffix(); }
    static const std::string getQuantizerCodebookFile(std::string indexPath) { return getPrefix(indexPath) + getThirdCentroidSuffix(); }
#else
    static const std::string getBlobFile(std::string indexPath) { return getPrefix(indexPath) + getThirdCentroidSuffix(); }
    static const std::string getQuantizerCodebookFile(std::string indexPath) { return getPrefix(indexPath) + getSecondCentroidSuffix(); }
#endif
    static const std::string getCodebookIndexFile(std::string indexPath) { return getPrefix(indexPath) + get3rdTo2ndSuffix(); }
    static const std::string getObjectIndexFile(std::string indexPath) { return getPrefix(indexPath) + getObjTo3rdSuffix(); }
    static const std::string getRotationFile(std::string indexPath) { return getPQFile(indexPath) + "/" + getRotationFile
(); }

    static const std::string getWorkspaceName() { return "ws"; }

    const std::string path;
    
    QuantizedBlobGraphRepository quantizedBlobGraph;

    bool searchable;


  };

}

#endif
