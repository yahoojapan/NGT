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

/***
  {
    // simple quantization and search example

    std::string	indexPath = "onng_index";	// ONNG
    std::string	queryPath = "query.tsv";	// Query file.
    NGTError err = ngt_create_error_object();
    
    // quantize the specified existing index
    //   build quantized objects and a quantized graph
    NGTQGQuantizationPrameters quantizationParameters;
    ngtqg_initialize_quantization_parameters(&quantizationParameters);
    ngtqg_quantize(indexPath.c_str(), quantizationParameters, err);

    // open the index (ANNG or ONNG).
    index = ngtqg_open_index(indexPath.c_str(), err);
    if (index == NULL) {
      std::cerr << ngt_get_error_string(err) << std::endl;
      return false;
    }

    std::ifstream is(queryPath);		// open a query file.
    if (!is) {
      std::cerr << "Cannot open the specified file. " << queryPath << std::endl;
      return false;
    }

    // get the dimension of the index to check the dimension of the query
    NGTProperty property = ngt_create_property(err);
    ngt_get_property(index, property, err);
    size_t dimension = ngt_get_property_dimension(property, err);
    ngt_destroy_property(property);

    std::string line;
    float queryVector[dimension];
    if (!getline(is, line)) {	// read a query object from the query file.
      std::cerr << "no data" << std::endl;
    }
    std::vector<std::string> tokens;
    NGT::Common::tokenize(line, tokens, " \t");      // split a string into words by the separators.
    // create a query vector from the tokens.
    if (tokens.size() != dimension) {
      std::cerr << "dimension of the query is invalid. dimesion=" << tokens.size() << ":" << dimension << std::endl;
      return false;
    }
    for (std::vector<std::string>::iterator ti = tokens.begin(); ti != tokens.end(); ++ti) {
      queryVector[distance(tokens.begin(), ti)] = NGT::Common::strtod(*ti);
    }
    // set search prameters.
    NGTObjectDistances result = ngt_create_empty_results(err);
    NGTQGQuery query;
    ngtqg_initialize_query(&query);
    query.query = queryVector;
    query.size = 20;
    query.epsilon = 0.03;
    query.result_expansion = 2;

    // search with the quantized graph
    bool status = ngtqg_search_index(index, query, result, err);
    auto rsize = ngt_get_result_size(result, err);

    // show resultant objects.
    std::cout << "Rank\tID\tDistance" << std::endl;
    for (size_t i = 0; i < rsize; i++) {
      NGTObjectDistance object = ngt_get_result(result, i, err);
      std::cout << i + 1 << "\t" << object.id << "\t" << object.distance << std::endl;
    }
    ngt_destroy_results(result);
    ngtqg_close_index(index);
  }
***/

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
  float		*query;
  size_t	size;		// # of returned objects
  float		epsilon;
  float		result_expansion;
  float		radius;
} NGTQGQuery;

typedef struct {
  float dimension_of_subvector;
  size_t max_number_of_edges;
} NGTQGQuantizationPrameters;

NGTQGIndex ngtqg_open_index(const char *, NGTError);

void ngtqg_close_index(NGTQGIndex);

void ngtqg_initialize_quantization_parameters(NGTQGQuantizationPrameters *);

void ngtqg_quantize(const char *, NGTQGQuantizationPrameters, NGTError);

void ngtqg_initialize_query(NGTQGQuery *);

bool ngtqg_search_index(NGTQGIndex, NGTQGQuery, NGTObjectDistances, NGTError);

#ifdef __cplusplus
}
#endif
