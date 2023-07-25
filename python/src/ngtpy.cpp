//
// Copyright (C) 2018 Yahoo Japan Corporation
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

#include	"NGT/Index.h"
#include	"NGT/GraphOptimizer.h"
#include	"NGT/version_defs.h"
#include	"NGT/NGTQ/Quantizer.h"
#ifdef NGTQ_QBG
#include	"NGT/NGTQ/QuantizedBlobGraph.h"
#else
#include	"NGT/NGTQ/QuantizedGraph.h"
#endif

#include	<pybind11/pybind11.h>
#include	<pybind11/stl.h>
#include	<pybind11/numpy.h>

namespace py = pybind11;

class Index : public NGT::Index {
public:
  Index(
   const std::string path, 		// ngt index path.
   bool readOnly,			// read only or not.
   bool zeroBasedNumbering,		// object ID numbering.
   bool treeDisabled,
   bool logDisabled			// stderr log is disabled.
  ):NGT::Index(path, readOnly) {
    zeroNumbering = zeroBasedNumbering;
    numOfDistanceComputations = 0;
    if (logDisabled) {
      NGT::Index::disableLog();
    } else {
      NGT::Index::enableLog();
    }
    treeIndex = !treeDisabled;
    defaultNumOfSearchObjects = 20;
    defaultEpsilon = 0.1;
    defaultRadius = FLT_MAX;
    defaultEdgeSize = -1;	// -1: use edge_size_for_search in the profile
    defaultExpectedAccuracy = -1.0;
  }

  static void create(
   const std::string path,
   size_t dimension,
   int edgeSizeForCreation = 10,
   int edgeSizeForSearch = 40,
   const std::string distanceType = "L2",
   const std::string objectType = "Float"
  ) {
    NGT::Property prop;
    prop.dimension = dimension;
    prop.edgeSizeForCreation = edgeSizeForCreation;
    prop.edgeSizeForSearch = edgeSizeForSearch;

    if (objectType == "Float" || objectType == "float") {
      prop.objectType = NGT::Index::Property::ObjectType::Float;
    } else if (objectType == "Byte" || objectType == "byte") {
      prop.objectType = NGT::Index::Property::ObjectType::Uint8;
#ifdef NGT_HALF_FLOAT
    } else if (objectType == "Float16" || objectType == "float16") {
      prop.objectType = NGT::Index::Property::ObjectType::Float16;
#endif
    } else {
      std::stringstream msg;
      msg << "ngtpy::create: invalid object type. " << objectType;
      NGTThrowException(msg);
    }

    if (distanceType == "L1") {
      prop.distanceType = NGT::Index::Property::DistanceType::DistanceTypeL1;
    } else if (distanceType == "L2") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeL2;
    } else if (distanceType == "Normalized L2") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeNormalizedL2;
    } else if (distanceType == "Hamming") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeHamming;
    } else if (distanceType == "Jaccard") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeJaccard;
    } else if (distanceType == "Sparse Jaccard") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeSparseJaccard;
    } else if (distanceType == "Angle") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeAngle;
    } else if (distanceType == "Normalized Angle") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeNormalizedAngle;
    } else if (distanceType == "Cosine") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeCosine;
    } else if (distanceType == "Normalized Cosine") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeNormalizedCosine;
    } else if (distanceType == "Normalized L2") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeNormalizedL2;
    } else {
      std::stringstream msg;
      msg << "ngtpy::create: invalid distance type. " << distanceType;
      NGTThrowException(msg);
    }
    NGT::Index::createGraphAndTree(path, prop);
  }

  void batchInsert(
   py::array_t<double> objects,
   size_t numThreads = 16,
   bool debug = false
  ) {
    py::buffer_info info = objects.request();
    if (debug) {
      std::cerr << info.shape.size() << ":" << info.shape[0] << ":" << info.shape[1] << std::endl;
    }
    if ((objects.flags() & py::detail::npy_api::constants::NPY_ARRAY_C_CONTIGUOUS_) == 0) {
      std::stringstream msg;
      msg << "ngtpy::batchInsert: Error! The array order is not C type. " << static_cast<int>(objects.flags())
	  << ":" << static_cast<int>(py::detail::npy_api::constants::NPY_ARRAY_C_CONTIGUOUS_);
      NGTThrowException(msg);
    }
    auto ptr = static_cast<double *>(info.ptr);
    assert(info.shape.size() == 2);
    NGT::Property prop;
    getProperty(prop);
    if (prop.dimension != info.shape[1]) {
      std::stringstream msg;
      msg << "ngtpy::insert: Error! dimensions are inconsitency. " << prop.dimension << ":" << info.shape[1];
      NGTThrowException(msg);
    }
    NGT::Index::append(ptr, info.shape[0]);
    NGT::Index::createIndex(numThreads);
    numOfDistanceComputations = 0;
  }

  int insert(
   py::array_t<double> object,
   bool debug = false
  ) {
    py::buffer_info info = object.request();
    auto ptr = static_cast<double *>(info.ptr);
    if (debug) {
      std::cerr << info.shape.size() << ":" << info.shape[0] << ":" << info.shape[1] << std::endl;
      for (int i = 0; i < info.shape[0]; i++) {
	std::cerr << *(ptr + i) << " ";
      }
      std::cerr << std::endl;
    }
    std::vector<double> v(ptr, ptr + info.shape[0]);
    int id = NGT::Index::insert(v);
    numOfDistanceComputations = 0;
    return zeroNumbering ? id - 1 : id;
  }

  py::object search(
   py::object query,
   size_t size = 0, 			// the number of resultant objects
   float epsilon = 0.1, 		// search parameter epsilon. the adequate range is from 0.0 to 0.15. negative value is acceptable.
   int edgeSize = -1,			// the number of used edges for each node during the exploration of the graph.
   float expectedAccuracy = -1.0,	// expected accuracy. if this is specified, epsilon that is calculated from this is used for search instead of the specified epsilon.
   bool withDistance = true
  ) {
    py::array_t<float> qobject(query);
    py::buffer_info qinfo = qobject.request();
    NGT::Object *ngtquery = 0;
    try {
      ngtquery = NGT::Index::allocateObject(static_cast<float*>(qinfo.ptr), qinfo.size);
    } catch (NGT::Exception &e) {
      std::cerr << e.what() << std::endl;
      if (!withDistance) {
	return py::array_t<int>();
      } else {
	return py::list();
      }
    }
    NGT::SearchContainer sc(*ngtquery);
    sc.setSize(size == 0 ? defaultNumOfSearchObjects : size);		// the number of resulting objects.
    sc.setRadius(defaultRadius);					// the radius of search.
    if (expectedAccuracy > 0.0) {
      sc.setExpectedAccuracy(expectedAccuracy);
    } else {
      sc.setEpsilon(epsilon <= -1.0 ? defaultEpsilon : epsilon);	// set exploration coefficient.
    }
    sc.setEdgeSize(edgeSize < -2 ? defaultEdgeSize : edgeSize);		// if maxEdge is negative, the specified value in advance is used.

    if (treeIndex) {
      NGT::Index::search(sc);
    } else {
      NGT::Index::searchUsingOnlyGraph(sc);
    }

    numOfDistanceComputations += sc.distanceComputationCount;

    NGT::Index::deleteObject(ngtquery);
    if (!withDistance) {
      NGT::ResultPriorityQueue &r = sc.getWorkingResult();
      py::array_t<int> ids(r.size());
      py::buffer_info idsinfo = ids.request();
      int *endptr = reinterpret_cast<int*>(idsinfo.ptr);
      int *ptr = endptr + (r.size() - 1);
      if (zeroNumbering) {
        while (ptr >= endptr) {
	  *ptr-- = r.top().id - 1;
	  r.pop();
        }
      } else {
        while (ptr >= endptr) {
	  *ptr-- = r.top().id;
	  r.pop();
        }
      }

      return ids;
    }
    py::list results;
    NGT::ObjectDistances r;
    r.moveFrom(sc.getWorkingResult());
    if (zeroNumbering) {
      for (auto ri = r.begin(); ri != r.end(); ++ri) {
	results.append(py::make_tuple((*ri).id - 1, (*ri).distance));
      }
    } else {
      for (auto ri = r.begin(); ri != r.end(); ++ri) {
	results.append(py::make_tuple((*ri).id, (*ri).distance));
      }
    }
    return results;
  }

  py::object linearSearch(
   py::object query,
   size_t size = 0, 			// the number of resultant objects
   bool withDistance = true
  ) {
    py::array_t<float> qobject(query);
    py::buffer_info qinfo = qobject.request();
    NGT::Object *ngtquery = 0;
    try {
      ngtquery = NGT::Index::allocateObject(static_cast<float*>(qinfo.ptr), qinfo.size);
    } catch (NGT::Exception &e) {
      std::cerr << e.what() << std::endl;
      if (!withDistance) {
	return py::array_t<int>();
      } else {
	return py::list();
      }
    }

    NGT::SearchContainer sc(*ngtquery);
    sc.setSize(size == 0 ? defaultNumOfSearchObjects : size);	// the number of resulting objects.
    sc.setRadius(defaultRadius);				// the radius of search.
    NGT::ObjectDistances rs;
    sc.setResults(&rs);

    NGT::Index::linearSearch(sc);

    numOfDistanceComputations += sc.distanceComputationCount;

    NGT::Index::deleteObject(ngtquery);
    if (!withDistance) {
      py::array_t<int> ids(rs.size());
      py::buffer_info idsinfo = ids.request();
      int *ptr = reinterpret_cast<int*>(idsinfo.ptr);
      if (zeroNumbering) {
	for (auto ri = rs.begin(); ri != rs.end(); ++ri) {
	  *ptr++ = (*ri).id - 1;
        }
      } else {
	for (auto ri = rs.begin(); ri != rs.end(); ++ri) {
	  *ptr++ = (*ri).id;
        }
      }
      return ids;
    }
    py::list results;
    if (zeroNumbering) {
      for (auto ri = rs.begin(); ri != rs.end(); ++ri) {
	results.append(py::make_tuple((*ri).id - 1, (*ri).distance));
      }
    } else {
      for (auto ri = rs.begin(); ri != rs.end(); ++ri) {
	results.append(py::make_tuple((*ri).id, (*ri).distance));
      }
    }
    return results;
  }

  void remove(size_t id) {
    id = zeroNumbering ? id + 1 : id;
    NGT::Index::remove(id);
  }

  void refineANNG(
    float epsilon,		// epsilon for search
    float accuracy,		// expected accuracy for search
    int numOfEdges,		// k to build AKNNG
    int numOfExploredEdges,	// # of explored edges for search
    size_t batchSize)		// batch size to search at the same time
  {
    bool unlog = NGT::Index::redirect;
    NGT::GraphReconstructor::refineANNG(*this, unlog, epsilon, accuracy, numOfEdges, numOfExploredEdges, batchSize);
  }

  std::vector<float> getObject(size_t id) {
    id = zeroNumbering ? id + 1 : id;
    NGT::Property prop;
    NGT::Index::getProperty(prop);
    std::vector<float> object;
    object.reserve(prop.dimension);
    switch (prop.objectType) {
    case NGT::ObjectSpace::ObjectType::Uint8:
      {
	auto *obj = static_cast<uint8_t*>(NGT::Index::getObjectSpace().getObject(id));
	for (int i = 0; i < prop.dimension; i++) {
	  object.push_back(*obj++);
	}
	break;
      }
    default:
    case NGT::ObjectSpace::ObjectType::Float:
      {
	auto *obj = static_cast<float*>(NGT::Index::getObjectSpace().getObject(id));
	for (int i = 0; i < prop.dimension; i++) {
	  object.push_back(*obj++);
	}
	break;
      }
    }
    return object;
  }

  void set(
   size_t numOfSearchObjects, 		// the number of resultant objects
   NGT::Distance radius,		// search radius.
   float epsilon,	 		// search parameter epsilon. the adequate range is from 0.0 to 0.05.
   int edgeSize,			// the number of edges for each node
   float expectedAccuracy		// Expected accuracy
  ) {
    defaultNumOfSearchObjects = numOfSearchObjects > 0 ? numOfSearchObjects : defaultNumOfSearchObjects;
    defaultEpsilon	      = epsilon > -1.0 ? epsilon : defaultEpsilon;
    defaultRadius	      = radius >= 0.0 ? radius : defaultRadius;
    defaultEdgeSize	      = edgeSize >= -2 ? edgeSize : defaultEdgeSize;
    defaultExpectedAccuracy   = expectedAccuracy > 0.0 ? expectedAccuracy : defaultExpectedAccuracy;
  }

  size_t getNumOfDistanceComputations() { return numOfDistanceComputations; }

  bool		zeroNumbering;	    // for object ID numbering. zero-based or one-based numbering.
  size_t	numOfDistanceComputations;
  size_t	numOfSearchObjects; // k
  bool		treeIndex;
  size_t	defaultNumOfSearchObjects; // k
  float		defaultEpsilon;
  float		defaultRadius;
  int64_t	defaultEdgeSize;
  float		defaultExpectedAccuracy;
};

class Optimizer : public NGT::GraphOptimizer {
public:
  using NGT::GraphOptimizer::GraphOptimizer;

  int optimizeNumberOfEdgesForANNG(
    const std::string path,		// anng index path
    int numOfQueries,
    int numOfResults,
    int numOfThreads,
    float targetAccuracy,
    int targetNoOfObjects,
    int numOfSampleObjects,
    int maxNoOfEdges) {

    NGT::GraphOptimizer::ANNGEdgeOptimizationParameter p;

    p.noOfQueries	= numOfQueries > 0 ? numOfQueries : p.noOfQueries;
    p.noOfResults	= numOfResults > 0 ? numOfResults : p.noOfResults;
    p.noOfThreads	= numOfThreads >= 0 ? numOfThreads : p.noOfThreads;
    p.targetAccuracy	= targetAccuracy > 0.0 ? targetAccuracy : p.targetAccuracy;
    p.targetNoOfObjects	= targetNoOfObjects >= 0 ? targetNoOfObjects : p.targetNoOfObjects;
    p.noOfSampleObjects	= numOfSampleObjects > 0 ? numOfSampleObjects : p.noOfSampleObjects;
    p.maxNoOfEdges	= maxNoOfEdges > 0 ? maxNoOfEdges : p.maxNoOfEdges;

    auto edge = NGT::GraphOptimizer::optimizeNumberOfEdgesForANNG(path, p);
    if (!logDisabled) {
      std::cerr << "the optimized number of edges is" << edge.first << "(" << edge.second << ")" << std::endl;
    }
    return edge.first;
  }

};

class QuantizedIndex : public NGTQG::Index {
public:
  QuantizedIndex(
   const  std::string path, 		// ngt index path.
   size_t maxNoOfEdges,			// the maximum number of quantized graph edges.
   bool   zeroBasedNumbering,		// object ID numbering.
   bool   treeDisabled,			// not use the tree index.
   bool   logDisabled			// stderr log is disabled.
  ):NGTQG::Index(path, maxNoOfEdges) {
    zeroNumbering = zeroBasedNumbering;
    numOfDistanceComputations = 0;
    treeIndex = !treeDisabled;
    withDistance = true;;
    defaultNumOfSearchObjects = 20;
    defaultEpsilon = 0.02;
    defaultRadius = FLT_MAX;
    defaultResultExpansion = 3.0;
    defaultEdgeSize = 0; // not used
    if (logDisabled) {
      NGT::Index::disableLog();
    } else {
      NGT::Index::enableLog();
    }
  }

  py::object search(
   py::object query,
   size_t size, 		// the number of resultant objects
   float epsilon, 		// search parameter epsilon. the adequate range is from 0.0 to 0.05.
   float resultExpansion,	// the number of inner resultant objects
   int edgeSize		// the number of used edges for each node during the exploration of the graph.
  ) {
    py::array_t<float> qobject(query);
    py::buffer_info qinfo = qobject.request();
    std::vector<float> qvector(static_cast<float*>(qinfo.ptr), static_cast<float*>(qinfo.ptr) + qinfo.size);
    try {
      NGTQG::SearchQuery sc(qvector);
      size		= size > 0 ? size : defaultNumOfSearchObjects;
      epsilon		= epsilon > -1.0 ? epsilon : defaultEpsilon;
      resultExpansion	= resultExpansion >= 0.0 ? resultExpansion : defaultResultExpansion;
      edgeSize		= edgeSize >= -2 ? edgeSize : defaultEdgeSize;
      sc.setSize(size);				// the number of resulting objects.
      sc.setRadius(defaultRadius);			// the radius of search.
      sc.setEpsilon(epsilon);			// set exploration coefficient.
      sc.setResultExpansion(resultExpansion);	// set result expansion.
      sc.setEdgeSize(edgeSize);			// if maxEdge is minus, the specified value in advance is used.

      NGTQG::Index::search(sc);

      numOfDistanceComputations += sc.distanceComputationCount;

      if (!withDistance) {
	NGT::ResultPriorityQueue &r = sc.getWorkingResult();
	py::array_t<int> ids(r.size());
	py::buffer_info idsinfo = ids.request();
	int *endptr = reinterpret_cast<int*>(idsinfo.ptr);
	int *ptr = endptr + (r.size() - 1);
	if (zeroNumbering) {
	  while (ptr >= endptr) {
	    *ptr-- = r.top().id - 1;
	    r.pop();
	  }
	} else {
	  while (ptr >= endptr) {
	    *ptr-- = r.top().id;
	    r.pop();
	  }
	}

	return ids;
      }
      py::list results;
      NGT::ObjectDistances r;
      r.moveFrom(sc.getWorkingResult());
      if (zeroNumbering) {
	for (auto ri = r.begin(); ri != r.end(); ++ri) {
	  results.append(py::make_tuple((*ri).id - 1, (*ri).distance));
	}
      } else {
	for (auto ri = r.begin(); ri != r.end(); ++ri) {
	  results.append(py::make_tuple((*ri).id, (*ri).distance));
	}
      }
      return results;
    } catch (NGT::Exception &e) {
      std::cerr << e.what() << std::endl;
      if (!withDistance) {
	return py::array_t<int>();
      } else {
	return py::list();
      }
    }
  }

  void setWithDistance(bool v) { withDistance = v; }

  void set(
   size_t numOfSearchObjects, 		// the number of resultant objects
   NGT::Distance radius,		// search radius.
   float epsilon,	 		// search parameter epsilon. the adequate range is from 0.0 to 0.05.
   float resultExpansion,		// the number of inner resultant objects
   int edgeSize				// not used
  ) {
    defaultNumOfSearchObjects = numOfSearchObjects > 0 ? numOfSearchObjects : defaultNumOfSearchObjects;
    defaultEpsilon	      = epsilon > -1.0 ? epsilon : defaultEpsilon;
    defaultRadius	      = radius >= 0.0 ? radius : defaultRadius;
    defaultResultExpansion    = resultExpansion >= 0.0 ? resultExpansion : defaultResultExpansion;
    defaultEdgeSize	      = edgeSize >= -2 ? edgeSize : defaultEdgeSize;
  }

  bool		zeroNumbering;	    // for object ID numbering. zero-based or one-based numbering.
  size_t	numOfDistanceComputations;
  bool		treeIndex;
  bool		withDistance;
  size_t	defaultNumOfSearchObjects; // k
  float		defaultEpsilon;
  float		defaultRadius;
  float		defaultResultExpansion;
  int64_t	defaultEdgeSize;
};


class BatchResults {
public:
  BatchResults() {}
  size_t getSize() { return size; }
  void convert() {
    if (results.size() == 0) {
      return;
    }
    resultList.clear();
    resultList.resize(results.size());
    for (size_t idx = 0; idx < size; idx++) {
      if (resultList[idx].size() != results[idx].size() && results[idx].size() != 0) {
	resultList[idx].clear();
	resultList[idx].resize(results[idx].size());
	size_t rank = results[idx].size();
	resultList[idx].resize(rank);
	rank--;
	while (!results[idx].empty()) {
	  resultList[idx][rank] = results[idx].top();
	  results[idx].pop();
	  rank--;
	}
      }
    }
    results.clear();
  }
  py::object get(size_t idx) {
    convert();
    if (idx >= size) {
      py::list result;
      return result;
    }
    py::list result;
    for (auto ri = resultList[idx].begin(); ri != resultList[idx].end(); ++ri) {
      result.append(py::make_tuple((*ri).id - 1, (*ri).distance));
    }
    return result;
  }
  py::array_t<int> getIDs() {
    convert();
    if (size == 0 || resultList[0].size() == 0) {
      NGTThrowException("ngtpy::BatchResults::get: empty.");
    }
    size_t nobjects = resultList[0].size();
    py::array_t<uint32_t> r({size, nobjects});
    auto wr = r.mutable_unchecked<2>();
    for (size_t idx = 0; idx < size; idx++) {
      if (resultList[idx].size() != nobjects) {
	NGTThrowException("ngtpy::BatchResults::get: not knn results.");
      }
      for (auto ri = resultList[idx].begin(); ri != resultList[idx].end(); ++ri) {
	wr(idx, std::distance(resultList[idx].begin(), ri)) = (*ri).id - 1;
      }
    }
    return r;
  }

  py::array_t<int> getIndexedIDs() {
    convert();
    size_t count = 0;
    for (size_t idx = 0; idx < size; idx++) {
      count += resultList[idx].size();
    }
    py::array_t<int> results(count);
    auto wresults = results.mutable_unchecked<1>();
    size_t pos = 0;
    for (size_t idx = 0; idx < size; idx++) {
      for (auto ri = resultList[idx].begin(); ri != resultList[idx].end(); ++ri) {
	wresults(pos++) = (*ri).id - 1;
      }
    }
    return results;
  }

  py::array_t<float> getIndexedDistances() {
    convert();
    size_t count = 0;
    for (size_t idx = 0; idx < size; idx++) {
      count += resultList[idx].size();
    }
    py::array_t<float> results(count);
    auto wresults = results.mutable_unchecked<1>();
    size_t pos = 0;
    for (size_t idx = 0; idx < size; idx++) {
      for (auto ri = resultList[idx].begin(); ri != resultList[idx].end(); ++ri) {
	wresults(pos++) = (*ri).distance;
      }
    }
    return results;
  }

  py::array_t<uint32_t> getIndex() {
    convert();
    py::array_t<int> results(size + 1);
    auto wresults = results.mutable_unchecked<1>();
    size_t count = 0;
    wresults(0) = 0;
    for (size_t idx = 0; idx < size; idx++) {
      count += resultList[idx].size();
      wresults(idx + 1) = count;
    }
    return results;
  }

  std::vector<NGT::ResultPriorityQueue> results;
  std::vector<NGT::ObjectDistances> resultList;
  size_t size;
};

class QuantizedBlobIndex : public QBG::Index {
public:
  QuantizedBlobIndex(
   const  std::string path, 		// ngt index path.
   size_t maxNoOfEdges,			// the maximum number of quantized graph edges.
   bool   zeroBasedNumbering,		// object ID numbering.
   bool   treeDisabled,			// not use the tree index.
   bool   logDisabled,			// stderr log is disabled.
   bool   readOnly			// open mode
  ):QBG::Index(path, readOnly) {
    zeroNumbering = zeroBasedNumbering;
    numOfDistanceComputations = 0;
    treeIndex = !treeDisabled;
    withDistance = true;;
    defaultNumOfSearchObjects = 20;
    defaultEpsilon = 0.02;
    defaultBlobEpsilon = 0.0;
    defaultResultExpansion = 3.0;
    defaultEdgeSize = -2;
    defaultExplorationSize = 200;
    defaultExactResultExpansion = 0.0;
    defaultNumOfProbes = 0;
  }

  void batchInsert(
   py::array_t<double> objects,
   bool debug = false
  ) {
    py::buffer_info info = objects.request();
    if (debug) {
      std::cerr << info.shape.size() << ":" << info.shape[0] << ":" << info.shape[1] << std::endl;
    }
    if ((objects.flags() & py::detail::npy_api::constants::NPY_ARRAY_C_CONTIGUOUS_) == 0) {
      std::stringstream msg;
      msg << "ngtpy::batchInsert: Error! The array order is not C type. " << static_cast<int>(objects.flags())
	  << ":" << static_cast<int>(py::detail::npy_api::constants::NPY_ARRAY_C_CONTIGUOUS_);
      NGTThrowException(msg);
    }
    auto ptr = static_cast<double *>(info.ptr);
    assert(info.shape.size() == 2);
    for (int idx = 0; idx < info.shape[0]; idx++) {
      if (debug) {
        for (int i = 0; i < info.shape[1]; i++) {
          std::cerr << *(ptr + i) << " ";
        }
        std::cerr << std::endl;
      }
      std::vector<float> v(ptr, ptr + info.shape[1]);
      ptr += info.shape[1];
      QBG::Index::append(v);
    }
  }

  py::array_t<uint32_t> batchSearchTmp(
   py::array_t<float> queries,
   size_t size	 		
  ) {
    const py::buffer_info &qinfo = queries.request();
    const std::vector<long int> &qshape = qinfo.shape;
    auto nOfQueries = qshape[0];
    auto dimension = qshape[1];
    auto *queryPtr = static_cast<float*>(qinfo.ptr);

    size	= size > 0 ? size : defaultNumOfSearchObjects;

    py::array_t<uint32_t> results({nOfQueries, static_cast<long>(size)});
    auto wresults = results.mutable_unchecked<2>();

#pragma omp parallel for schedule(dynamic)
    for (int idx = 0; idx < nOfQueries; idx++) {
      NGT::Object queryObject(dimension * sizeof(float));
      float *qptr = queryPtr + idx * dimension;
      memcpy(queryObject.getPointer(), qptr, dimension * sizeof(float));
      QBG::SearchContainer sc(queryObject);
      sc.setSize(size);
      sc.setEpsilon(defaultEpsilon);
      sc.setBlobEpsilon(defaultBlobEpsilon);
      sc.setEdgeSize(defaultEdgeSize);
      sc.setGraphExplorationSize(defaultExplorationSize);
      QBG::Index::search(sc);
      NGT::ResultPriorityQueue &r = sc.getWorkingResult();
      if (r.size() != size) {
	std::cerr << "result size is invalid? " << r.size() << ":" << size << std::endl;
      }
      size_t rank = size;
      while (!r.empty()) {
	rank--;
	wresults(idx, rank) = r.top().id - 1;
	r.pop();
      }
    }
    return results;
  }

  void batchSearchInTwoSteps(
    py::array_t<float> queries,
    BatchResults &results,
    size_t size
  ) {
    const py::buffer_info &qinfo = queries.request();
    const std::vector<long int> &qshape = qinfo.shape;
    auto nOfQueries = qshape[0];
    size_t dimension = qshape[1];
    size_t psedoDimension = QBG::Index::getQuantizer().globalCodebookIndex.getObjectSpace().getPaddedDimension();
    auto *queryPtr = static_cast<float*>(qinfo.ptr);

    size	= size > 0 ? size : defaultNumOfSearchObjects;

    results.results.clear();
    results.resultList.clear();
    results.results.resize(nOfQueries);

#pragma omp parallel for schedule(dynamic)
    for (int idx = 0; idx < nOfQueries; idx++) {
      float *qptr = queryPtr + idx * dimension;
      vector<float> query(psedoDimension, 0);
      memcpy(query.data(), qptr, dimension * sizeof(float));
      QBG::SearchContainer sc;
      sc.setObjectVector(query);
      sc.setSize(size);
      sc.setEpsilon(defaultEpsilon);
      sc.setBlobEpsilon(defaultBlobEpsilon);
      sc.setEdgeSize(defaultEdgeSize);
      sc.setNumOfProbes(defaultNumOfProbes);
#ifdef NGTQBG_FUNCTION_SELECTOR
      sc.functionSelector = defaultFunctionSelector;	/////////////////
#endif
      QBG::Index::searchInTwoSteps(sc);
      results.results[idx] = std::move(sc.getWorkingResult());
    }
    results.size = results.results.size();
    return;
  }

  void batchSearchInOneStep(
    py::array_t<float> queries,
    BatchResults &results,
    size_t size
  ) {
    const py::buffer_info &qinfo = queries.request();
    const std::vector<long int> &qshape = qinfo.shape;
    auto nOfQueries = qshape[0];
    size_t dimension = qshape[1];
    size_t psedoDimension = QBG::Index::getQuantizer().globalCodebookIndex.getObjectSpace().getPaddedDimension();
    auto *queryPtr = static_cast<float*>(qinfo.ptr);

    size	= size > 0 ? size : defaultNumOfSearchObjects;

    results.results.clear();
    results.resultList.clear();
    results.results.resize(nOfQueries);

    size_t searchSize = size;
    size_t searchExactResultSize = 0;
    if (defaultExactResultExpansion >= 1.0) {
      searchSize = static_cast<float>(size) * defaultExactResultExpansion;
      searchExactResultSize = size;
    }
#pragma omp parallel for schedule(dynamic)
    for (int idx = 0; idx < nOfQueries; idx++) {
      float *qptr = queryPtr + idx * dimension;
      vector<float> query(psedoDimension, 0);
      memcpy(query.data(), qptr, dimension * sizeof(float));
      QBG::SearchContainer sc;
      sc.setObjectVector(query);
      sc.setSize(searchSize);
      sc.setExactResultSize(searchExactResultSize);
      sc.setEpsilon(defaultEpsilon);
      sc.setBlobEpsilon(defaultBlobEpsilon);
      sc.setEdgeSize(defaultEdgeSize);
      sc.setGraphExplorationSize(defaultExplorationSize);
      QBG::Index::search(sc);
      results.results[idx] = std::move(sc.getWorkingResult());
      //QBG::Index::getQuantizer().globalCodebookIndex.deleteObject(queryObject);
    }
    results.size = results.results.size();
    return;
  }

  void batchSearch(
    py::array_t<float> queries,
    BatchResults &results,
    size_t size
  ) {
    if ((queries.flags() & py::detail::npy_api::constants::NPY_ARRAY_C_CONTIGUOUS_) == 0) {
      std::stringstream msg;
      msg << "ngtpy::batchSearch: Error! The array order is not C type. " << static_cast<int>(queries.flags())
	  << ":" << static_cast<int>(py::detail::npy_api::constants::NPY_ARRAY_C_CONTIGUOUS_);
      NGTThrowException(msg);
    }
    if (defaultNumOfProbes == 0) {
      batchSearchInOneStep(queries, results, size);
    } else {
      batchSearchInTwoSteps(queries, results, size);
    }
    return;
  }

  void batchRangeSearch(
    py::array_t<float> queries,
    BatchResults &results,
    float radius
  ) {
    if ((queries.flags() & py::detail::npy_api::constants::NPY_ARRAY_C_CONTIGUOUS_) == 0) {
      std::stringstream msg;
      msg << "ngtpy::batchRangeSearch: Error! The array order is not C type. " << static_cast<int>(queries.flags())
	  << ":" << static_cast<int>(py::detail::npy_api::constants::NPY_ARRAY_C_CONTIGUOUS_);
      NGTThrowException(msg);
    }
    const py::buffer_info &qinfo = queries.request();
    const std::vector<long int> &qshape = qinfo.shape;
    auto nOfQueries = qshape[0];
    size_t dimension = qshape[1];
    size_t psedoDimension = QBG::Index::getQuantizer().globalCodebookIndex.getObjectSpace().getPaddedDimension();
    auto *queryPtr = static_cast<float*>(qinfo.ptr);
    radius = radius >= 0 ? radius : defaultRadius;
    radius = sqrt(radius);

    results.results.clear();
    results.resultList.clear();
    results.results.resize(nOfQueries);

#pragma omp parallel for schedule(dynamic)
    for (int idx = 0; idx < nOfQueries; idx++) {
      float *qptr = queryPtr + idx * dimension;
      vector<float> query(psedoDimension, 0);
      memcpy(query.data(), qptr, dimension * sizeof(float));
      QBG::SearchContainer sc;
      sc.setObjectVector(query);
      sc.setSize(std::numeric_limits<uint32_t>::max());
      sc.setRadius(radius);
      sc.setEpsilon(defaultEpsilon);
      sc.setBlobEpsilon(defaultBlobEpsilon);
      sc.setEdgeSize(defaultEdgeSize);
      sc.setGraphExplorationSize(defaultExplorationSize);
      QBG::Index::search(sc);
      results.results[idx] = std::move(sc.getWorkingResult());
    }
    results.size = results.results.size();
    return;
  }

  py::object search(
    py::object query,
    size_t size, 		// the number of resultant objects
    float epsilon 		// search parameter epsilon. the adequate range is from 0.0 to 0.05.
  ) {
    py::array_t<float> qobject(query);
    py::buffer_info qinfo = qobject.request();
    std::vector<float> qvector(static_cast<float*>(qinfo.ptr), static_cast<float*>(qinfo.ptr) + qinfo.size);
    try {
      QBG::SearchContainer sc;
      sc.setObjectVector(qvector);
      size		= size > 0 ? size : defaultNumOfSearchObjects;
      epsilon		= epsilon > -1.0 ? epsilon : defaultEpsilon;
      if (defaultExactResultExpansion >= 1.0) {
	sc.setSize(static_cast<float>(size) * defaultExactResultExpansion);
	sc.setExactResultSize(size);
      } else {
	sc.setSize(size);				// the number of resulting objects.
      }
      sc.setEpsilon(epsilon);			// set exploration coefficient.
      sc.setBlobEpsilon(defaultBlobEpsilon);
      sc.setEdgeSize(defaultEdgeSize);
      sc.setGraphExplorationSize(defaultExplorationSize);
      sc.setNumOfProbes(defaultNumOfProbes);
#ifdef NGTQBG_FUNCTION_SELECTOR
      sc.functionSelector = defaultFunctionSelector;	/////////////////
#endif
      QBG::Index::searchInTwoSteps(sc);

      numOfDistanceComputations += sc.distanceComputationCount;

      if (!withDistance) {
	//-/std::cerr << "without distances" << std::endl;
	NGT::ResultPriorityQueue &r = sc.getWorkingResult();
	py::array_t<int> ids(r.size());
	py::buffer_info idsinfo = ids.request();
	int *endptr = reinterpret_cast<int*>(idsinfo.ptr);
	int *ptr = endptr + (r.size() - 1);
	if (zeroNumbering) {
	  while (ptr >= endptr) {
	    *ptr-- = r.top().id - 1;
	    r.pop();
	  }
	} else {
	  while (ptr >= endptr) {
	    *ptr-- = r.top().id;
	    r.pop();
	  }
	}

	return ids;
      }
      py::list results;
      NGT::ObjectDistances r;
      r.moveFrom(sc.getWorkingResult());
      if (zeroNumbering) {
	for (auto ri = r.begin(); ri != r.end(); ++ri) {
	  results.append(py::make_tuple((*ri).id - 1, (*ri).distance));
	}
      } else {
	for (auto ri = r.begin(); ri != r.end(); ++ri) {
	  results.append(py::make_tuple((*ri).id, (*ri).distance));
	}
      }
      return results;
    } catch (NGT::Exception &e) {
      std::cerr << e.what() << std::endl;
      if (!withDistance) {
	return py::array_t<int>();
      } else {
	return py::list();
      }
    }
  }

  void setWithDistance(bool v) { withDistance = v; }

  void set(
   size_t numOfSearchObjects, 		// the number of resultant objects
   float epsilon,	 		// search parameter epsilon. the adequate range is from 0.0 to 0.05.
   float blobEpsilon,	 		// search parameter blob epsilon. the recommended value is 0.0.
   float resultExpansion,		// the number of inner resultant objects
   float radius,
   int edgeSize,
   int explorationSize,			// the maximum number of nodes that are searched
   float exactResultExpansion,
   int numOfProbes
#ifdef NGTQBG_FUNCTION_SELECTOR
   , size_t functionSelector ///////////////////
#endif
  ) {
    defaultNumOfSearchObjects = numOfSearchObjects > 0 ? numOfSearchObjects : defaultNumOfSearchObjects;
    defaultEpsilon	      = epsilon > -1.0 ? epsilon : defaultEpsilon;
    defaultBlobEpsilon	      = blobEpsilon > -1.0 ? blobEpsilon : defaultBlobEpsilon;
    defaultResultExpansion    = resultExpansion >= 0.0 ? resultExpansion : defaultResultExpansion;
    defaultEdgeSize	      = edgeSize >= -2 ? edgeSize : defaultEdgeSize;
    defaultExplorationSize    = explorationSize > 0 ? explorationSize : defaultExplorationSize;
    defaultRadius	      = radius >= 0.0 ? radius : defaultRadius;
    defaultExactResultExpansion = exactResultExpansion > 0.0 ? exactResultExpansion : defaultExactResultExpansion;
    defaultNumOfProbes	      = numOfProbes > 0 ? numOfProbes : defaultNumOfProbes;
#ifdef NGTQBG_FUNCTION_SELECTOR
    defaultFunctionSelector   = functionSelector; ////////////////
#endif
  }

  bool		zeroNumbering;	    // for object ID numbering. zero-based or one-based numbering.
  size_t	numOfDistanceComputations;
  bool		treeIndex;
  bool		withDistance;
  size_t	defaultNumOfSearchObjects; // k
  float		defaultEpsilon;

  float		defaultBlobEpsilon;
  float		defaultResultExpansion;
  int64_t	defaultEdgeSize;
  size_t	defaultExplorationSize;
  float		defaultRadius;
  float		defaultExactResultExpansion;
  size_t	defaultNumOfProbes;
#ifdef NGTQBG_FUNCTION_SELECTOR
  size_t	defaultFunctionSelector; //////////////////////////
#endif
};

PYBIND11_MODULE(ngtpy, m) {
    m.doc() = "ngt python";

    m.attr("__version__") = NGT_VERSION;

    m.def("create", &::Index::create,
          py::arg("path"),
          py::arg("dimension"),
          py::arg("edge_size_for_creation") = 10,
          py::arg("edge_size_for_search") = 40,
          py::arg("distance_type") = "L2",
          py::arg("object_type") = "Float");

    py::class_<Index>(m, "Index")
      .def(py::init<const std::string &, bool, bool, bool, bool>(),
           py::arg("path"),
           py::arg("read_only") = false,
           py::arg("zero_based_numbering") = true,
	   py::arg("tree_disabled") = false,
           py::arg("log_disabled") = false)
      .def("search", &::Index::search,
           py::arg("query"),
           py::arg("size") = 0,
           py::arg("epsilon") = -FLT_MAX,
           py::arg("edge_size") = INT_MIN,
           py::arg("expected_accuracy") = -FLT_MAX,
           py::arg("with_distance") = true)
      .def("linear_search", &::Index::linearSearch,
           py::arg("query"),
           py::arg("size") = 0,
           py::arg("with_distance") = true)
      .def("get_num_of_distance_computations", &::Index::getNumOfDistanceComputations)
      .def("save", (void (NGT::Index::*)()) &NGT::Index::save)
      .def("close", &NGT::Index::close)
      .def("remove", &::Index::remove,
           py::arg("object_id"))
      .def("build_index", &NGT::Index::createIndex,
           py::arg("num_threads") = 8,
           py::arg("target_size_of_graph") = 0)
      .def("get_num_of_objects", &::Index::getNumberOfObjects)
      .def("get_size_of_object_repository", &::Index::getObjectRepositorySize)
      .def("get_size_of_graph_repository", &::Index::getGraphRepositorySize)
      .def("get_object", &::Index::getObject,
           py::arg("object_id"))
      .def("batch_insert", &::Index::batchInsert,
           py::arg("objects"),
           py::arg("num_threads") = 8,
           py::arg("debug") = false)
      .def("insert", &::Index::insert,
           py::arg("object"),
           py::arg("debug") = false)
      .def("refine_anng", &::Index::refineANNG,
	   py::arg("epsilon") = 0.1,
	   py::arg("expected_accuracy") = 0.0,
	   py::arg("num_of_edges") = 0,
	   py::arg("num_of_explored_edges") = INT_MIN,
	   py::arg("batch_size") = 10000)
      .def("set", &::Index::set,
           py::arg("num_of_search_objects") = 0,
	   py::arg("search_radius") = -FLT_MAX,
	   py::arg("epsilon") = -FLT_MAX,
	   py::arg("edge_size") = INT_MIN,
	   py::arg("expected_accuracy") = -FLT_MAX)
      .def("export_index", (void (NGT::Index::*)(const std::string&)) &NGT::Index::exportIndex,
           py::arg("path"))
      .def("import_index", (void (NGT::Index::*)(const std::string&)) &NGT::Index::importIndex,
           py::arg("path"));

    py::class_<Optimizer>(m, "Optimizer")
      .def(py::init<int, int, int, int, float, float, float, float, double, double, bool>(),
	   py::arg("num_of_outgoings") = -1,
	   py::arg("num_of_incomings") = -1,
	   py::arg("num_of_queries") = -1,
	   py::arg("num_of_objects") = -1,
	   py::arg("low_accuracy_from") = -1.0,
	   py::arg("low_accuracy_to") = -1.0,
	   py::arg("high_accuracy_from") = -1.0,
	   py::arg("high_accuracy_to") = -1.0,
	   py::arg("gt_epsilon") = -DBL_MAX,
	   py::arg("margin") = -1.0,
	   py::arg("log_disabled") = false)
      .def("execute", &NGT::GraphOptimizer::execute,
	   py::arg("in_path"),
	   py::arg("out_path"))
      .def("adjust_search_coefficients", &NGT::GraphOptimizer::adjustSearchCoefficients,
	   py::arg("path"))
      .def("set", (void (NGT::GraphOptimizer::*)(int, int, int, int, float, float, float, float,
						 double, double)) &NGT::GraphOptimizer::set,
	   py::arg("num_of_outgoings") = -1,
	   py::arg("num_of_incomings") = -1,
	   py::arg("num_of_queries") = -1,
	   py::arg("num_of_objects") = -1,
	   py::arg("low_accuracy_from") = -1.0,
	   py::arg("low_accuracy_to") = -1.0,
	   py::arg("high_accuracy_from") = -1.0,
	   py::arg("high_accuracy_to") = -1.0,
	   py::arg("gt_epsilon") = -DBL_MAX,
	   py::arg("margin") = -1.0)
      .def("set_processing_modes", &NGT::GraphOptimizer::setProcessingModes,
	   py::arg("shortcut_reduction") = true,
	   py::arg("search_parameter_optimization") = true,
	   py::arg("prefetch_parameter_optimization") = true,
	   py::arg("accuracy_table_generation") = true)
      .def("optimize_search_parameters", &NGT::GraphOptimizer::optimizeSearchParameters,
	   py::arg("path"))
      .def("optimize_number_of_edges_for_anng", &::Optimizer::optimizeNumberOfEdgesForANNG,
	   py::arg("path"),
	   py::arg("num_of_queries") = -1,
	   py::arg("num_of_results") = -1,
	   py::arg("num_of_threads") = -1,
	   py::arg("target_accuracy") = -1,
	   py::arg("target_num_of_objects") = -1,
	   py::arg("num_of_sample_objects") = -1,
	   py::arg("max_num_of_edges") = -1);

    py::class_<QuantizedIndex>(m, "QuantizedIndex")
      .def(py::init<const std::string &, size_t, bool, bool, bool>(),
           py::arg("path"),
	   py::arg("max_no_of_edges") = 128,
           py::arg("zero_based_numbering") = true,
	   py::arg("tree_disabled") = false,
           py::arg("log_disabled") = false)
      .def("search", &::QuantizedIndex::search,
           py::arg("query"),
           py::arg("size") = 0,
           py::arg("epsilon") = -FLT_MAX,
           py::arg("result_expansion") = -FLT_MAX,
           py::arg("edge_size") = INT_MIN)
      .def("set_with_distance", &::QuantizedIndex::setWithDistance,
           py::arg("boolean") = true)
      .def("set", &::QuantizedIndex::set,
           py::arg("num_of_search_objects") = 0,
	   py::arg("search_radius") = -FLT_MAX,
           py::arg("epsilon") = -FLT_MAX,
           py::arg("result_expansion") = -FLT_MAX,
           py::arg("edge_size") = INT_MIN)
      // set_defaults is deprecated
      .def("set_defaults", &::QuantizedIndex::set,
           py::arg("size") = 0,
	   py::arg("search_radius") = -FLT_MAX,
           py::arg("epsilon") = -FLT_MAX,
           py::arg("result_expansion") = -FLT_MAX,
           py::arg("edge_size") = INT_MIN);


    py::class_<QuantizedBlobIndex>(m, "QuantizedBlobIndex")
      .def(py::init<const std::string &, size_t, bool, bool, bool, bool>(),
           py::arg("path"),
	   py::arg("max_no_of_edges") = 128,
           py::arg("zero_based_numbering") = true,
	   py::arg("tree_disabled") = false,
           py::arg("log_disabled") = false,
           py::arg("read_only") = true)
      .def("save", (void (QBG::Index::*)()) &QBG::Index::save)
      .def("batch_insert", &::QuantizedBlobIndex::batchInsert,
           py::arg("objects"),
           py::arg("debug") = false)
      .def("batch_search_tmp", &::QuantizedBlobIndex::batchSearchTmp,
    	   py::arg("query"),
           py::arg("size") = 0)
      .def("batch_search", &::QuantizedBlobIndex::batchSearch,
    	   py::arg("query"),
	   py::arg("results"),
           py::arg("size") = 0)
      .def("batch_range_search", &::QuantizedBlobIndex::batchRangeSearch,
    	   py::arg("query"),
	   py::arg("results"),
           py::arg("radius") = -FLT_MAX)
      .def("search", &::QuantizedBlobIndex::search,
           py::arg("query"),
           py::arg("size") = 0,
           py::arg("epsilon") = -FLT_MAX)
      .def("set_with_distance", &::QuantizedBlobIndex::setWithDistance,
           py::arg("boolean") = true)
      .def("set", &::QuantizedBlobIndex::set,
           py::arg("num_of_search_objects") = 0,
           py::arg("epsilon") = -FLT_MAX,
           py::arg("blob_epsilon") = -FLT_MAX,
           py::arg("result_expansion") = -FLT_MAX,
           py::arg("radius") = -FLT_MAX,
	   py::arg("edge_size") = INT_MIN,
           py::arg("exploration_size") = 0,
           py::arg("exact_result_expansion") = 0.0,
           py::arg("num_of_probes") = INT_MIN
#ifdef NGTQBG_FUNCTION_SELECTOR
	   , py::arg("function_selector") = 0 ///////////
#endif
	   );

    py::class_<BatchResults>(m, "BatchResults")
      .def(py::init<>())
      .def("get", &::BatchResults::get,
	   py::arg("position"))
      .def("get_ids", &::BatchResults::getIDs)
      .def("get_indexed_ids", &::BatchResults::getIndexedIDs)
      .def("get_indexed_distances", &::BatchResults::getIndexedDistances)
      .def("get_index", &::BatchResults::getIndex)
      .def("get_size", &::BatchResults::getSize);

}
