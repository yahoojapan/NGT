//
// Copyright (C) 2018-2019 Yahoo Japan Corporation
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

#include	<pybind11/pybind11.h>
#include	<pybind11/stl.h>
#include	<pybind11/numpy.h>

namespace py = pybind11;

class Index : public NGT::Index {
public:
  Index(
   const string path, 			// ngt index path.
   bool readOnly,			// read only or not.
   bool zeroBasedNumbering,		// object ID numbering.
   bool logDisabled			// stderr log is disabled.
  ):NGT::Index(path, readOnly) {
    zeroNumbering = zeroBasedNumbering;
    numOfDistanceComputations = 0;
    numOfSearchObjects = 10;
    searchRadius = FLT_MAX;
    if (logDisabled) {
      NGT::Index::disableLog();
    } else {
      NGT::Index::enableLog();
    }
  }

  static void create(
   const string path,
   size_t dimension,
   int edgeSizeForCreation = 10,
   int edgeSizeForSearch = 40,
   const string distanceType = "L2",
   const string objectType = "Float"
  ) {
    NGT::Property prop;
    prop.dimension = dimension;
    prop.edgeSizeForCreation = edgeSizeForCreation;
    prop.edgeSizeForSearch = edgeSizeForSearch;

    if (objectType == "Float" || objectType == "float") {
      prop.objectType = NGT::Index::Property::ObjectType::Float;
    } else if (objectType == "Byte" || objectType == "byte") {
      prop.objectType = NGT::Index::Property::ObjectType::Uint8;
    } else {
      std::stringstream msg;
      msg << "ngtpy::create: invalid object type. " << objectType;
      NGTThrowException(msg);
    }

    if (distanceType == "L1") {
      prop.distanceType = NGT::Index::Property::DistanceType::DistanceTypeL1;
    } else if (distanceType == "L2") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeL2;
    } else if (distanceType == "Hamming") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeHamming;
    } else if (distanceType == "Jaccard") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeJaccard;
    } else if (distanceType == "Angle") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeAngle;
    } else if (distanceType == "Normalized Angle") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeNormalizedAngle;
    } else if (distanceType == "Cosine") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeCosine;
    } else if (distanceType == "Normalized Cosine") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeNormalizedCosine;
    } else {
      std::stringstream msg;
      msg << "ngtpy::create: invalid distance type. " << distanceType;
      NGTThrowException(msg);
    }
    NGT::Index::createGraphAndTree(path, prop);
  }

  void batchInsert(
   py::array_t<double> objects, 
   size_t numThreads = 8,
   bool debug = false
  ) {
    py::buffer_info info = objects.request();
    if (debug) {
      std::cerr << info.shape.size() << ":" << info.shape[0] << ":" << info.shape[1] << std::endl;
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
    return NGT::Index::insert(v);
    numOfDistanceComputations = 0;
  }

  py::object search(
   py::object query,
   size_t size = 0, 			// the number of resultant objects
   float epsilon = 0.1, 		// search parameter epsilon. the adequate range is from 0.0 to 0.15. minus value is acceptable.
   int edgeSize = -1,			// the number of used edges for each node during the exploration of the graph.
   bool withDistance = true
  ) {
    py::array_t<float> qobject(query);
    py::buffer_info qinfo = qobject.request();
    NGT::Object *ngtquery = 0;
    try {
      ngtquery = NGT::Index::allocateObject(static_cast<float*>(qinfo.ptr), qinfo.size);
    } catch (NGT::Exception &e) {
      std::cerr << e.what() << endl;
      if (!withDistance) {
	return py::array_t<int>();
      } else {
	return py::list();
      }
    }
    NGT::SearchContainer sc(*ngtquery);
    if (size == 0) {
      sc.setSize(numOfSearchObjects);		// the number of resulting objects.
    } else {
      sc.setSize(size);				// the number of resulting objects.
    }
    sc.setRadius(searchRadius);			// the radius of search.
    sc.setEpsilon(epsilon);			// set exploration coefficient.
    sc.setEdgeSize(edgeSize);			// if maxEdge is minus, the specified value in advance is used.

    NGT::Index::search(sc);

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
      std::cerr << e.what() << endl;
      if (!withDistance) {
	return py::array_t<int>();
      } else {
	return py::list();
      }
    }

    NGT::SearchContainer sc(*ngtquery);
    if (size == 0) {
      sc.setSize(numOfSearchObjects);		// the number of resulting objects.
    } else {
      sc.setSize(size);				// the number of resulting objects.
    }
    sc.setRadius(searchRadius);			// the radius of search.
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

  vector<float> getObject(size_t id) {
    id = zeroNumbering ? id + 1 : id;
    NGT::Property prop;
    NGT::Index::getProperty(prop);
    vector<float> object;
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

  void set(size_t k, NGT::Distance r) {
    if (k > 0) {
      numOfSearchObjects = k;
    }
    if (r >= 0.0) {
      searchRadius = r; 
    }
  }

  size_t getNumOfDistanceComputations() { return numOfDistanceComputations; }

  bool		zeroNumbering;	    // for object ID numbering. zero-based or one-based numbering.
  size_t	numOfDistanceComputations;
  size_t	numOfSearchObjects; // k
  NGT::Distance	searchRadius;
};

PYBIND11_MODULE(ngtpy, m) {
    m.doc() = "ngt python";

    m.def("create", &::Index::create, 
          py::arg("path"), 
          py::arg("dimension"), 
          py::arg("edge_size_for_creation") = 10, 
          py::arg("edge_size_for_search") = 40, 
          py::arg("distance_type") = "L2", 
          py::arg("object_type") = "Float");

    py::class_<Index>(m, "Index")
      .def(py::init<const std::string &, bool, bool, bool>(), 
           py::arg("path"),
           py::arg("read_only") = false,
           py::arg("zero_based_numbering") = true,
           py::arg("log_disabled") = false)
      .def("search", &::Index::search, 
           py::arg("query"), 
           py::arg("size") = 0, 
           py::arg("epsilon") = 0.1, 
           py::arg("edge_size") = -1,
           py::arg("with_distance") = true)
      .def("linear_search", &::Index::linearSearch, 
           py::arg("query"), 
           py::arg("size") = 0, 
           py::arg("with_distance") = true)
      .def("get_num_of_distance_computations", &::Index::getNumOfDistanceComputations)
      .def("save", &NGT::Index::save)
      .def("close", &NGT::Index::close)
      .def("remove", &::Index::remove, 
           py::arg("object_id"))
      .def("build_index", &NGT::Index::createIndex, 
           py::arg("num_threads") = 8)
      .def("get_object", &::Index::getObject, 
           py::arg("object_id"))
      .def("batch_insert", &::Index::batchInsert, 
           py::arg("objects"),
           py::arg("num_threads") = 8, 
           py::arg("debug") = false)
      .def("insert", &::Index::insert, 
           py::arg("object"),
           py::arg("debug") = false)
      .def("set", &::Index::set, 
           py::arg("num_of_search_objects") = 0,
	   py::arg("search_radius") = -1.0);

    py::class_<NGT::GraphOptimizer>(m, "Optimizer")
      .def(py::init<int, int, int, float, float, float, float, double, double, bool>(),
	   py::arg("num_of_outgoings") = -1,
	   py::arg("num_of_incomings") = -1,
	   py::arg("num_of_queries") = -1,
	   py::arg("low_accuracy_from") = -1.0,
	   py::arg("low_accuracy_to") = -1.0,
	   py::arg("high_accuracy_from") = -1.0,
	   py::arg("high_accuracy_to") = -1.0,
	   py::arg("gt_epsilon") = DBL_MIN,
	   py::arg("merge") = -1.0,
	   py::arg("log_disabled") = false)
      .def("execute", &NGT::GraphOptimizer::execute, 
	   py::arg("in_index_path"),
	   py::arg("out_index_path"))
      .def("adjust_search_coefficients", &NGT::GraphOptimizer::adjustSearchCoefficients, 
	   py::arg("index_path"))
      .def("set", &NGT::GraphOptimizer::set, 
	   py::arg("num_of_outgoings") = -1,
	   py::arg("num_of_incomings") = -1,
	   py::arg("num_of_queries") = -1,
	   py::arg("low_accuracy_from") = -1.0,
	   py::arg("low_accuracy_to") = -1.0,
	   py::arg("high_accuracy_from") = -1.0,
	   py::arg("high_accuracy_to") = -1.0,
	   py::arg("gt_epsilon") = DBL_MIN,
	   py::arg("merge") = -1.0);
}

