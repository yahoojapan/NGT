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

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#include "NGT/Capi.h"

  typedef void* NGTQGIndex;
  typedef NGTObjectDistance NGTObjectDistance;
  typedef NGTError NGTQGError;

  typedef struct {
    float	*query;
    size_t	size;		// # of returned objects
    float	epsilon;
    float	result_expansion;
    float	radius;
  } NGTQGQuery;

  typedef struct {
    size_t	size;		// # of returned objects
    float	epsilon;
    float	result_expansion;
    float	radius;
  } NGTQGQueryParameters;

  typedef struct {
    float					*query;
    NGTQGQueryParameters	params;
  } NGTQGQueryFloat;

  typedef struct {
    uint8_t					*query;
    NGTQGQueryParameters	params;
  } NGTQGQueryUint8;

  typedef struct {
    NGTFloat16				*query;
    NGTQGQueryParameters	params;
  } NGTQGQueryFloat16;

  typedef struct {
    float	dimension_of_subvector;
    size_t	max_number_of_edges;
  } NGTQGQuantizationParameters;

  NGTQGIndex ngtqg_open_index(const char *, NGTQGError);

  void ngtqg_close_index(NGTQGIndex);

  void ngtqg_initialize_quantization_parameters(NGTQGQuantizationParameters *);

  bool ngtqg_quantize(const char *, NGTQGQuantizationParameters, NGTQGError);

  void ngtqg_initialize_query_parameters(NGTQGQueryParameters *);

  void ngtqg_initialize_query(NGTQGQuery *);

  bool ngtqg_search_index(NGTQGIndex, NGTQGQuery, NGTObjectDistances, NGTQGError);

  bool ngtqg_search_index_float(NGTQGIndex, NGTQGQueryFloat, NGTObjectDistances, NGTQGError);

  bool ngtqg_search_index_uint8(NGTQGIndex, NGTQGQueryUint8, NGTObjectDistances, NGTQGError);

  bool ngtqg_search_index_float16(NGTQGIndex, NGTQGQueryFloat16, NGTObjectDistances, NGTQGError);

  // QBG CAPI

  typedef void* QBGIndex;
  typedef NGTError QBGError;
  typedef NGTObjectDistances QBGObjectDistances;

  uint32_t qbg_get_result_size(QBGObjectDistances results, NGTError error);

  NGTObjectDistance qbg_get_result(const QBGObjectDistances results, const uint32_t idx, NGTError error);

  void qbg_destroy_results(QBGObjectDistances results);

  typedef struct {
    size_t	extended_dimension;
    size_t	dimension;
    size_t	number_of_subvectors;
    size_t	number_of_blobs;
    int		internal_data_type;
    int		data_type;
    int		distance_type;
  } QBGConstructionParameters;

  typedef struct {
    // hierarchical kmeans
    int		hierarchical_clustering_init_mode;
    size_t	number_of_first_objects;
    size_t	number_of_first_clusters;
    size_t	number_of_second_objects;
    size_t	number_of_second_clusters;
    size_t	number_of_third_clusters;
    // optimization
    size_t	number_of_objects;
    size_t	number_of_subvectors;
    int		optimization_clustering_init_mode;
    size_t	rotation_iteration;
    size_t	subvector_iteration;
    size_t	number_of_matrices;
    bool	rotation;
    bool	repositioning;
  } QBGBuildParameters;

  typedef struct {
    size_t	number_of_results;
    float	epsilon;
    float	blob_epsilon;
    float	result_expansion;
    size_t	number_of_explored_blobs;
    size_t	number_of_edges;
    float	radius;
  } QBGQueryParameters;

  typedef struct {
    float	*query;
    size_t	number_of_results;
    float	epsilon;
    float	blob_epsilon;
    float	result_expansion;
    size_t	number_of_explored_blobs;
    size_t	number_of_edges;
    float	radius;
  } QBGQuery;

  typedef struct {
    float		   		*query;
    QBGQueryParameters	params;
  } QBGQueryFloat;

  typedef struct {
    uint8_t		   		*query;
    QBGQueryParameters	params;
  } QBGQueryUint8;

  typedef struct {
    NGTFloat16	   		*query;
    QBGQueryParameters	params;
  } QBGQueryFloat16;

  void qbg_initialize_construction_parameters(QBGConstructionParameters *parameters);

  bool qbg_create(const char *indexPath, QBGConstructionParameters *parameters, QBGError error);

  QBGIndex qbg_open_index(const char *index_path, bool read_only, QBGError error);

  void qbg_close_index(QBGIndex index);

  bool qbg_save_index(QBGIndex index, QBGError error);

  ObjectID qbg_append_object(QBGIndex index, float *obj, uint32_t obj_dim, QBGError error);

  ObjectID qbg_append_object_as_uint8(QBGIndex index, uint8_t *obj, uint32_t obj_dim, QBGError error);

  ObjectID qbg_append_object_as_float16(QBGIndex index, NGTFloat16 *obj, uint32_t obj_dim, QBGError error);

  void qbg_initialize_build_parameters(QBGBuildParameters *parameters);

  bool qbg_build_index(const char *index_path, QBGBuildParameters *parameters, QBGError error);

  void qbg_initialize_query_parameters(QBGQueryParameters *parameters);

  void qbg_initialize_query(QBGQuery *parameters);

  bool qbg_search_index(QBGIndex index, QBGQuery query, NGTObjectDistances results, QBGError error);

  bool qbg_search_index_float(QBGIndex index, QBGQueryFloat query, NGTObjectDistances results, QBGError error);

  bool qbg_search_index_uint8(QBGIndex index, QBGQueryUint8 query, NGTObjectDistances results, QBGError error);

  bool qbg_search_index_float16(QBGIndex index, QBGQueryFloat16 query, NGTObjectDistances results, QBGError error);

  float* qbg_get_object(QBGIndex index, ObjectID id,  QBGError error);

  uint8_t* qbg_get_object_as_uint8(QBGIndex index, ObjectID id,  QBGError error);

  NGTFloat16* qbg_get_object_as_float16(QBGIndex index, ObjectID id, QBGError error);

  size_t qbg_get_dimension(QBGIndex index, QBGError error);

#ifdef __cplusplus
}
#endif
