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

#include	<pybind11/pybind11.h>
#include	<pybind11/stl.h>
#include	<pybind11/numpy.h>

namespace py = pybind11;

class Index : public NGT::Index {
public:
  Index(
   const string path, 			// ngt index path.
   bool zeroBasedNumbering = true	// object ID numbering.
  ):NGT::Index(path) {
    indexDecrement = zeroBasedNumbering ? 1 : 0;
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
      std::cerr << "ngtpy::create: invalid object type. " << objectType << std::endl;
      return;
    }

    if (distanceType == "L1") {
      prop.distanceType = NGT::Index::Property::DistanceType::DistanceTypeL1;
    } else if (distanceType == "L2") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeL2;
    } else if (distanceType == "Hamming") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeHamming;
    } else if (distanceType == "Angle") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeAngle;
    } else if (distanceType == "Cosine") {
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeCosine;
    } else {
      std::cerr << "ngtpy::create: invalid distance type. " << distanceType << std::endl;
      return;
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
      std::cerr << "ngtpy::insert: Error! dimensions are inconsitency. " << prop.dimension << ":" << info.shape[1] << std::endl;
      return;
    }
///-</ delete for open release
#if 0
    for (int idx = 0; idx < info.shape[0]; idx++) {
      if (debug) {
        for (int i = 0; i < info.shape[1]; i++) {
          std::cerr << *(ptr + i) << " ";
        }
        std::cerr << std::endl;
      }
      std::vector<double> v(ptr, ptr + info.shape[1]);
      ptr += info.shape[1];
      NGT::Index::insert(v);
    }
#else
///->/
    NGT::Index::append(ptr, info.shape[0]);
///-</ delete for open release
#endif
///->/
    NGT::Index::createIndex(numThreads);
  }

  std::vector<int> search(
   std::vector<double> query,		// query
   size_t size = 10, 			// the number of resultant objects
   float epsilon = 0.1, 		// search parameter epsilon. the adequate range is from 0.0 to 0.15. minus value is acceptable.
   int edgeSize = -1			// the number of used edges for each node during the exploration of the graph.
  ) {
    NGT::Object *ngtquery = NGT::Index::allocateObject(query);
    NGT::SearchContainer sc(*ngtquery);
    NGT::ObjectDistances objects;
    sc.setResults(&objects);			// set the result set.
    sc.setSize(size);				// the number of resultant objects.
    sc.setEpsilon(epsilon);			// set exploration coefficient.
    sc.setEdgeSize(edgeSize);			// if maxEdge is minus, the specified value in advance is used.

    NGT::Index::search(sc);

    std::vector<int> ids;
    for (size_t i = 0; i < objects.size(); i++) {
      ids.push_back(objects[i].id - indexDecrement);
    }
///-</ delete for open release
#if 0
    for (size_t i = 0; i < objects.size(); i++) {
      std::cerr << objects[i].id - indexDecrement << ":" << objects[i].distance << " ";
    }
    std::cerr << std::endl;
#endif
///->/
    NGT::Index::deleteObject(ngtquery);

    return ids;
  }

  size_t indexDecrement;	// for object ID numbering. zero-based or one-based numbering.
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
      .def(py::init<const std::string &, bool>(), 
           py::arg("path"), 
           py::arg("zero_based_numbering") = true)
      .def("search", &::Index::search, 
           py::arg("query"), 
           py::arg("size") = 10, 
           py::arg("epsilon") = 0.1, 
           py::arg("edge_size") = -1)
      .def("save", &NGT::Index::saveIndex,
           py::arg("path"))
      .def("close", &NGT::Index::close)
      .def("batch_insert", &::Index::batchInsert, 
           py::arg("objects"),
           py::arg("num_threads") = 8, 
           py::arg("debug") = false);
}

