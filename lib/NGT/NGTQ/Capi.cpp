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

#include <string>
#include <iostream>
#include <sstream>

#include "NGT/Capi.h"
#include "NGT/NGTQ/Quantizer.h"
#include "NGT/NGTQ/Capi.h"
#include "NGT/NGTQ/QuantizedGraph.h"
#include "NGT/NGTQ/QuantizedBlobGraph.h"
#include "NGT/NGTQ/Optimizer.h"
#include "NGT/NGTQ/HierarchicalKmeans.h"

#ifdef NGTQ_QBG

static bool operate_error_string_(const std::stringstream &ss, NGTError error) {
  if (error != NULL) {
    try {
      std::string *error_str = static_cast<std::string *>(error);
      *error_str             = ss.str();
    } catch (std::exception &err) {
      std::cerr << ss.str() << " > " << err.what() << std::endl;
      return false;
    }
  } else {
    std::cerr << ss.str() << std::endl;
  }
  return true;
}

void ngtqg_initialize_query_parameters(NGTQGQueryParameters *parameters) {
  parameters->size             = 20;
  parameters->epsilon          = 0.03;
  parameters->result_expansion = 3.0;
  parameters->radius           = FLT_MAX;
}

void ngtqg_initialize_query(NGTQGQuery *query) {
  query->query            = 0;
  query->size             = 20;
  query->epsilon          = 0.03;
  query->result_expansion = 3.0;
  query->radius           = FLT_MAX;
}

NGTQGIndex ngtqg_open_index(const char *index_path, NGTError error) {
  try {
    std::string index_path_str(index_path);
    auto *index = new NGTQG::Index(index_path_str);
    index->disableLog();
    return static_cast<NGTQGIndex>(index);
  } catch (std::exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return NULL;
  }
}

void ngtqg_close_index(NGTQGIndex index) {
  if (index == NULL) return;
  (static_cast<NGTQG::Index *>(index))->close();
  delete static_cast<NGTQG::Index *>(index);
}

static bool ngtqg_search_index_(NGTQG::Index *pindex, NGTQG::SearchQuery *sq, NGTQGQueryParameters &params,
                                NGTObjectDistances results) {
  sq->setResults(static_cast<NGT::ObjectDistances *>(results)); // set the result set.
  sq->setSize(params.size);                                     // the number of resultant objects.
  sq->setRadius(params.radius);                                 // search radius.
  sq->setEpsilon(params.epsilon);                               // exploration coefficient.
  sq->setResultExpansion(params.result_expansion);              // result expansion.

  pindex->search(*sq);

  return true;
}

template <typename T>
static bool ngtqg_search_index_(NGTQGIndex index, T *query, NGTQGQueryParameters &params,
                                NGTObjectDistances results, NGTError error) {
  if (index == NULL || query == NULL || results == NULL) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : parametor error: index = " << index << " query = " << query
       << " results = " << results;
    operate_error_string_(ss, error);
    return false;
  }

  NGTQG::Index *pindex = static_cast<NGTQG::Index *>(index);
  int32_t dim          = pindex->getObjectSpace().getDimension();

  NGT::Object *ngtquery = NULL;

  if (params.radius < 0.0) {
    params.radius = FLT_MAX;
  }

  try {
    std::vector<T> vquery(&query[0], &query[dim]);
    NGTQG::SearchQuery sq(vquery);
    ngtqg_search_index_(pindex, &sq, params, results);
  } catch (std::exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    if (ngtquery != NULL) {
      pindex->deleteObject(ngtquery);
    }
    return false;
  }
  return true;
}

bool ngtqg_search_index(NGTQGIndex index, NGTQGQuery query, NGTObjectDistances results, NGTError error) {
  NGTQGQueryParameters qparams = {
      .size             = query.size,
      .epsilon          = query.epsilon,
      .result_expansion = query.result_expansion,
      .radius           = query.radius,
  };
  return ngtqg_search_index_(index, query.query, qparams, results, error);
}

bool ngtqg_search_index_float(NGTQGIndex index, NGTQGQueryFloat query, NGTObjectDistances results,
                              NGTError error) {
  return ngtqg_search_index_(index, query.query, query.params, results, error);
}

bool ngtqg_search_index_uint8(NGTQGIndex index, NGTQGQueryUint8 query, NGTObjectDistances results,
                              NGTError error) {
  return ngtqg_search_index_(index, query.query, query.params, results, error);
}

bool ngtqg_search_index_float16(NGTQGIndex index, NGTQGQueryFloat16 query, NGTObjectDistances results,
                                NGTError error) {
  return ngtqg_search_index_(index, static_cast<NGT::float16 *>(query.query), query.params, results, error);
}

void ngtqg_initialize_quantization_parameters(NGTQGQuantizationParameters *parameters) {
  parameters->dimension_of_subvector = 0;
  parameters->max_number_of_edges    = 128;
}

bool ngtqg_quantize(const char *indexPath, NGTQGQuantizationParameters parameters, NGTError error) {
  try {
    NGTQG::Index::quantize(indexPath, parameters.dimension_of_subvector, parameters.max_number_of_edges,
                           false);
    return true;
  } catch (std::exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return false;
  }
}

uint32_t qbg_get_result_size(QBGObjectDistances results, NGTError error) {
  return ngt_get_result_size(results, error);
}

NGTObjectDistance qbg_get_result(const QBGObjectDistances results, const uint32_t idx, NGTError error) {
  return ngt_get_result(results, idx, error);
}

void qbg_destroy_results(QBGObjectDistances results) { ngt_destroy_results(results); }

void qbg_initialize_construction_parameters(QBGConstructionParameters *parameters) {
  parameters->extended_dimension   = 0;
  parameters->dimension            = 0;
  parameters->number_of_subvectors = 1;
  parameters->number_of_blobs      = 0;
  parameters->internal_data_type   = NGTQ::DataTypeFloat;
  parameters->data_type            = NGTQ::DataTypeFloat;
  parameters->distance_type        = NGTQ::DistanceType::DistanceTypeL2;
}

bool qbg_create(const char *indexPath, QBGConstructionParameters *parameters, NGTQGError error) {

  try {
    std::vector<float> r;
    NGTQ::Property property;
    NGT::Property globalProperty;
    NGT::Property localProperty;
    property.dimension = parameters->extended_dimension;
    if (property.dimension == 0) {
      property.dimension = parameters->dimension;
    }
    property.genuineDimension          = parameters->dimension;
    property.globalRange               = 0;
    property.localRange                = 0;
    property.globalCentroidLimit       = parameters->number_of_blobs;
    property.localCentroidLimit        = 16;
    property.localDivisionNo           = parameters->number_of_subvectors;
    property.singleLocalCodebook       = false;
    property.centroidCreationMode      = NGTQ::CentroidCreationModeStaticLayer;
    property.localCentroidCreationMode = NGTQ::CentroidCreationModeStatic;
    property.localIDByteSize           = 1;
    property.dataType                  = static_cast<NGTQ::DataType>(parameters->internal_data_type);
    property.genuineDataType           = static_cast<ObjectFile::DataType>(parameters->data_type);
    property.distanceType              = static_cast<NGTQ::DistanceType>(parameters->distance_type);

    globalProperty.edgeSizeForCreation        = 10;
    globalProperty.edgeSizeForSearch          = 40;
    globalProperty.indexType                  = NGT::Property::GraphAndTree;
    globalProperty.insertionRadiusCoefficient = 1.1;

    localProperty.indexType                  = globalProperty.indexType;
    localProperty.insertionRadiusCoefficient = globalProperty.insertionRadiusCoefficient;

    std::vector<float> *rotation = 0;
    const std::string objectPath;
    QBG::Index::create(indexPath, property, globalProperty, localProperty, rotation, objectPath);
  } catch (NGT::Exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return false;
  }

  return true;
}

QBGIndex qbg_open_index(const char *index_path, bool prebuilt, QBGError error) {
  try {
    std::string index_path_str(index_path);
    auto *index = new QBG::Index(index_path_str, prebuilt);
    return static_cast<QBGIndex>(index);
  } catch (std::exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return NULL;
  }
}

void qbg_close_index(QBGIndex index) {
  if (index == NULL) return;
  (static_cast<QBG::Index *>(index))->close();
  delete static_cast<QBG::Index *>(index);
  index = 0;
}

bool qbg_save_index(QBGIndex index, QBGError error) {
  if (index == NULL) return false;
  try {
    (static_cast<QBG::Index *>(index))->save();
  } catch (std::exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return false;
  }
  return true;
}

template <typename T>
static ObjectID qbg_append_object_(QBGIndex index, T *obj, uint32_t obj_dim, QBGError error) {
  if (index == NULL || obj == NULL || obj_dim == 0) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : parametor error: index = " << index << " obj = " << obj
       << " obj_dim = " << obj_dim;
    operate_error_string_(ss, error);
    return 0;
  }

  try {
    auto *pindex = static_cast<QBG::Index *>(index);
    std::vector<T> vobj(&obj[0], &obj[obj_dim]);
    return pindex->append(vobj);
  } catch (std::exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return 0;
  }
}

ObjectID qbg_append_object(QBGIndex index, float *obj, uint32_t obj_dim, QBGError error) {
  return qbg_append_object_(index, obj, obj_dim, error);
}

ObjectID qbg_append_object_as_uint8(QBGIndex index, uint8_t *obj, uint32_t obj_dim, QBGError error) {
  return qbg_append_object_(index, obj, obj_dim, error);
}

ObjectID qbg_append_object_as_float16(QBGIndex index, NGTFloat16 *obj, uint32_t obj_dim, QBGError error) {
  return qbg_append_object_(index, static_cast<NGT::float16 *>(obj), obj_dim, error);
}

template <typename T>
static ObjectID qbg_insert_object_(QBGIndex index, T *obj, uint32_t obj_dim, QBGError error) {
  if (index == NULL || obj == NULL || obj_dim == 0) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : parametor error: index = " << index << " obj = " << obj
       << " obj_dim = " << obj_dim;
    operate_error_string_(ss, error);
    return 0;
  }

  try {
    auto *pindex = static_cast<QBG::Index *>(index);
    std::vector<float> vobj(&obj[0], &obj[obj_dim]);
    return pindex->insert(vobj);
  } catch (std::exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return 0;
  }
}

ObjectID qbg_insert_object(QBGIndex index, float *obj, uint32_t obj_dim, QBGError error) {
  return qbg_insert_object_(index, obj, obj_dim, error);
}

ObjectID qbg_insert_object_as_uint8(QBGIndex index, uint8_t *obj, uint32_t obj_dim, QBGError error) {
  return qbg_insert_object_(index, obj, obj_dim, error);
}

ObjectID qbg_insert_object_as_float16(QBGIndex index, NGTFloat16 *obj, uint32_t obj_dim, QBGError error) {
  return qbg_insert_object_(index, static_cast<NGT::float16 *>(obj), obj_dim, error);
}

template <typename T>
uint32_t qbg_insert_objects_(QBGIndex index, T *obj, uint32_t n_of_objects, ObjectID *ids, NGTError error) {
  if (index == NULL || obj == NULL || n_of_objects == 0 || ids == NULL) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : parametor error: index = " << index << " obj = " << obj
       << " n_of_objects = " << n_of_objects << " ids=" << ids;
    operate_error_string_(ss, error);
    return 0;
  }

  try {
    auto *pindex = static_cast<QBG::Index *>(index);

    int32_t dim = pindex->getQuantizer().property.genuineDimension;
    std::vector<std::vector<float>> objs;
    objs.reserve(n_of_objects);
    for (size_t idx = 0; idx < n_of_objects; idx++) {
      std::vector<float> vobj(&obj[idx * dim], &obj[idx * dim + dim]);
      objs.emplace_back(vobj);
    }
    std::vector<uint32_t> insertedIDs;
    pindex->insert(objs, insertedIDs);
    for (size_t idx = 0; idx < insertedIDs.size(); idx++) {
      ids[idx] = insertedIDs[idx];
    }
    return insertedIDs.size();
  } catch (std::exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return 0;
  }
}

uint32_t qbg_insert_objects(QBGIndex index, float *objects, uint32_t n_of_objects, ObjectID *ids,
                            NGTError error) {
  return qbg_insert_objects_(index, objects, n_of_objects, ids, error);
}

uint32_t qbg_insert_objects_as_uint8(QBGIndex index, uint8_t *objects, uint32_t n_of_objects, ObjectID *ids,
                                     NGTError error) {
  return qbg_insert_objects_(index, objects, n_of_objects, ids, error);
}

uint32_t qbg_insert_objects_as_float16(QBGIndex index, NGTFloat16 *objects, uint32_t n_of_objects,
                                       ObjectID *ids, NGTError error) {
  return qbg_insert_objects_(index, static_cast<NGT::float16 *>(objects), n_of_objects, ids, error);
}

bool qbg_remove_objects(QBGIndex index, ObjectID *ids, uint32_t n_of_ids, QBGError error) {
  if (index == NULL || n_of_ids == 0 || ids == NULL) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : parametor error: index=" << index << " n_of_ids=" << n_of_ids
       << " ids=" << ids;
    operate_error_string_(ss, error);
    return false;
  }

  try {
    auto *pindex = static_cast<QBG::Index *>(index);
    std::vector<uint32_t> removedIDs(&ids[0], &ids[n_of_ids]);
    pindex->remove(removedIDs);
  } catch (std::exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return false;
  }
  return true;
}

bool qbg_remove_object(QBGIndex index, ObjectID id, QBGError error) {
  if (index == NULL || id == 0) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : parametor error: index=" << index << " id=" << id;
    operate_error_string_(ss, error);
    return false;
  }

  try {
    auto *pindex = static_cast<QBG::Index *>(index);
    pindex->remove(id);
  } catch (std::exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return false;
  }
  return true;
}

void qbg_initialize_build_parameters(QBGBuildParameters *parameters) {
  parameters->hierarchical_clustering_init_mode =
      static_cast<int>(NGT::Clustering::InitializationModeKmeansPlusPlus);
  parameters->number_of_first_objects   = 0;
  parameters->number_of_first_clusters  = 0;
  parameters->number_of_second_objects  = 0;
  parameters->number_of_second_clusters = 0;
  parameters->number_of_third_clusters  = 0;

  parameters->number_of_objects    = 1000;
  parameters->number_of_subvectors = 1;
  parameters->optimization_clustering_init_mode =
      static_cast<int>(NGT::Clustering::InitializationModeKmeansPlusPlus);
  parameters->rotation_iteration  = 2000;
  parameters->subvector_iteration = 400;
  parameters->number_of_matrices  = 3;
  parameters->rotation            = true;
  parameters->repositioning       = false;
}

bool qbg_build_index(const char *index_path, QBGBuildParameters *parameters, QBGError error) {

  QBG::HierarchicalKmeans hierarchicalKmeans;

  hierarchicalKmeans.maxSize            = 1000;
  hierarchicalKmeans.numOfClusters      = 2;
  hierarchicalKmeans.numOfTotalClusters = 0;
  hierarchicalKmeans.numOfTotalBlobs    = 0;
  hierarchicalKmeans.clusterID          = -1;
  hierarchicalKmeans.initMode =
      static_cast<NGT::Clustering::InitializationMode>(parameters->hierarchical_clustering_init_mode);
  hierarchicalKmeans.numOfRandomObjects  = 0;
  hierarchicalKmeans.extractCentroid     = false;
  hierarchicalKmeans.numOfFirstObjects   = parameters->number_of_first_objects;
  hierarchicalKmeans.numOfFirstClusters  = parameters->number_of_first_clusters;
  hierarchicalKmeans.numOfSecondObjects  = parameters->number_of_second_objects;
  hierarchicalKmeans.numOfSecondClusters = parameters->number_of_second_clusters;
  hierarchicalKmeans.numOfThirdClusters  = parameters->number_of_third_clusters;
  hierarchicalKmeans.numOfObjects        = 0;
  //-/hierarchicalKmeans.threeLayerClustering = true;
  hierarchicalKmeans.clusteringType = QBG::HierarchicalKmeans::ClusteringTypeThreeLayer;
  hierarchicalKmeans.verbose        = false;

  try {
    hierarchicalKmeans.clustering(index_path);
  } catch (NGT::Exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return false;
  }

  QBG::Optimizer optimizer;

  optimizer.numberOfObjects    = parameters->number_of_objects;
  optimizer.numberOfClusters   = 16;
  optimizer.numberOfSubvectors = 0;
  optimizer.clusteringType     = NGT::Clustering::ClusteringTypeKmeansWithNGT;
  optimizer.initMode =
      static_cast<NGT::Clustering::InitializationMode>(parameters->optimization_clustering_init_mode);
  optimizer.convergenceLimitTimes = 5;
  optimizer.iteration             = parameters->rotation_iteration;
  optimizer.clusterIteration      = parameters->subvector_iteration;
  optimizer.clusterSizeConstraint = false;
  optimizer.numberOfMatrices      = parameters->number_of_matrices;
  optimizer.seedNumberOfSteps     = 2;
  optimizer.seedStep              = 10;
  optimizer.reject                = 0.9;
  optimizer.timelimit             = 24 * 2;
  optimizer.timelimit *= 60.0 * 60.0;
  optimizer.rotation      = parameters->rotation;
  optimizer.repositioning = parameters->repositioning;
  optimizer.globalType    = QBG::Optimizer::GlobalTypeNone;
  optimizer.verbose       = false;

  try {
    auto nthreads = omp_get_max_threads();
    optimizer.optimize(index_path, nthreads);
  } catch (NGT::Exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return false;
  }

  try {
    auto verbose = false;
    QBG::Index::build(index_path, verbose);
  } catch (NGT::Exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return false;
  }
  return true;
}

void qbg_initialize_query_parameters(QBGQueryParameters *parameters) {
  parameters->number_of_results        = 20;
  parameters->epsilon                  = 0.1;
  parameters->blob_epsilon             = 0.0;
  parameters->result_expansion         = 3.0;
  parameters->number_of_explored_blobs = 256;
  parameters->number_of_edges          = 0;
  parameters->radius                   = 0;
}

void qbg_initialize_query(QBGQuery *parameters) {
  parameters->query                    = 0;
  parameters->number_of_results        = 20;
  parameters->epsilon                  = 0.1;
  parameters->blob_epsilon             = 0.0;
  parameters->result_expansion         = 3.0;
  parameters->number_of_explored_blobs = 256;
  parameters->number_of_edges          = 0;
  parameters->radius                   = 0;
}

static bool qbg_search_index_(QBG::Index *pindex, std::vector<float> &query, QBGQueryParameters &param,
                              NGTObjectDistances results) {
  // set search parameters.
  if (param.radius < 0.0) {
    param.radius = FLT_MAX;
  }

  QBG::SearchContainer sc;
  sc.setObjectVector(query);
  sc.setResults(static_cast<NGT::ObjectDistances *>(results));
  if (param.result_expansion >= 1.0) {
    sc.setSize(static_cast<float>(param.number_of_results) * param.result_expansion);
    sc.setExactResultSize(param.number_of_results);
  } else {
    sc.setSize(param.number_of_results);
    sc.setExactResultSize(0);
  }
  sc.setEpsilon(param.epsilon);
  sc.setBlobEpsilon(param.blob_epsilon);
  sc.setEdgeSize(param.number_of_edges);
  sc.setGraphExplorationSize(param.number_of_explored_blobs);

  pindex->search(sc);

  return true;
}

template <typename T>
static bool qbg_search_index_(QBGIndex index, T *query, QBGQueryParameters &params,
                              NGTObjectDistances results, QBGError error) {
  if (index == NULL || query == NULL || results == NULL) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : parametor error: index = " << index << " query = " << query
       << " results = " << results;
    operate_error_string_(ss, error);
    return false;
  }

  auto *pindex = static_cast<QBG::Index *>(index);
  int32_t dim  = pindex->getQuantizer().property.genuineDimension;

  try {
    std::vector<T> oquery(&query[0], &query[dim]);
    std::vector<float> vquery(oquery.begin(), oquery.end());
    qbg_search_index_(pindex, vquery, params, results);
  } catch (std::exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return false;
  }

  return true;
}

bool qbg_search_index(QBGIndex index, QBGQuery query, NGTObjectDistances results, QBGError error) {
  QBGQueryParameters qparams = {
      .number_of_results        = query.number_of_results,
      .epsilon                  = query.epsilon,
      .blob_epsilon             = query.blob_epsilon,
      .result_expansion         = query.result_expansion,
      .number_of_explored_blobs = query.number_of_explored_blobs,
      .number_of_edges          = query.number_of_edges,
      .radius                   = query.radius,
  };
  return qbg_search_index_(index, query.query, qparams, results, error);
}

bool qbg_search_index_float(QBGIndex index, QBGQueryFloat query, NGTObjectDistances results, QBGError error) {
  return qbg_search_index_(index, query.query, query.params, results, error);
}

bool qbg_search_index_uint8(QBGIndex index, QBGQueryUint8 query, NGTObjectDistances results, QBGError error) {
  return qbg_search_index_(index, query.query, query.params, results, error);
}

bool qbg_search_index_float16(QBGIndex index, QBGQueryFloat16 query, NGTObjectDistances results,
                              QBGError error) {
  return qbg_search_index_(index, static_cast<NGT::float16 *>(query.query), query.params, results, error);
}

template <typename T> static T *qbg_get_object_(QBGIndex index, ObjectID id, QBGError error) {
  if (index == NULL) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : parametor error: index = " << index;
    operate_error_string_(ss, error);
    return 0;
  }

  auto *pindex = static_cast<QBG::Index *>(index);

  try {
    if (typeid(T) == typeid(float)) {
      auto object = pindex->getObject(id);
      size_t size = sizeof(float) * object.size();
      auto obj    = malloc(size);
      if (obj == 0) {
        std::stringstream ss;
        ss << "Capi : " << __FUNCTION__ << "() : Error: Cannot allocate memory.";
        operate_error_string_(ss, error);
        return 0;
      }
      memcpy(obj, object.data(), size);
      return static_cast<T *>(obj);
    } else {
      auto o = pindex->getObject(id);
      std::vector<T> object(o.begin(), o.end());
      size_t size = sizeof(T) * object.size();
      auto obj    = malloc(size);
      if (obj == 0) {
        std::stringstream ss;
        ss << "Capi : " << __FUNCTION__ << "() : Error: Cannot allocate memory.";
        operate_error_string_(ss, error);
        return 0;
      }
      memcpy(obj, object.data(), size);
      return static_cast<T *>(obj);
    }
  } catch (std::exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return 0;
  }
}

float *qbg_get_object(QBGIndex index, ObjectID id, QBGError error) {
  return qbg_get_object_<float>(index, id, error);
}

uint8_t *qbg_get_object_as_uint8(QBGIndex index, ObjectID id, QBGError error) {
  return qbg_get_object_<uint8_t>(index, id, error);
}

NGTFloat16 *qbg_get_object_as_float16(QBGIndex index, ObjectID id, QBGError error) {
  return qbg_get_object_<NGT::float16>(index, id, error);
}

size_t qbg_get_dimension(QBGIndex index, QBGError error) {
  if (index == NULL) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : parametor error: index = " << index;
    operate_error_string_(ss, error);
    return 0;
  }
  auto *pindex = static_cast<QBG::Index *>(index);
  return pindex->getQuantizer().property.genuineDimension;
}

#endif
