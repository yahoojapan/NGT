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
#include "NGT/NGTQ/Capi.h"
#include "NGT/NGTQ/QuantizedGraph.h"

static bool operate_error_string_(const std::stringstream &ss, NGTError error){
  if(error != NULL){  
    try{
      std::string *error_str = static_cast<std::string*>(error);
      *error_str = ss.str();
    }catch(std::exception &err){
      std::cerr << ss.str() << " > " << err.what() << std::endl;
      return false;
    }
  }else{
    std::cerr << ss.str() << std::endl;
  }
  return true;
}

void ngtqg_initialize_query(NGTQGQuery *query) {
  query->query = 0;
  query->size = 20;
  query->epsilon = 0.03;
  query->result_expansion = 3.0;
  query->radius = FLT_MAX;
}

NGTQGIndex ngtqg_open_index(const char *index_path, NGTError error) {
  try{
    std::string index_path_str(index_path);
    auto *index = new NGTQG::Index(index_path_str);
    index->disableLog();
    return static_cast<NGTQGIndex>(index);
  }catch(std::exception &err){
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    return NULL;
  }
}

void ngtqg_close_index(NGTQGIndex index) {
    if(index == NULL) return;
    (static_cast<NGTQG::Index*>(index))->close();
    delete static_cast<NGTQG::Index*>(index);
}

static bool ngtqg_search_index_(NGTQG::Index* pindex, std::vector<float> &query, NGTQGQuery &param, NGTObjectDistances results) {
  // set search prameters.
  NGTQG::SearchQuery sq(query);  // Query.

  sq.setResults(static_cast<NGT::ObjectDistances*>(results));          // set the result set.
  sq.setSize(param.size);                        // the number of resultant objects.
  sq.setRadius(param.radius);                    // search radius.
  sq.setEpsilon(param.epsilon);                  // exploration coefficient.
  sq.setResultExpansion(param.result_expansion); // result expansion.

  auto tmp = static_cast<NGT::ObjectDistances*>(results);

  pindex->search(sq);

  return true;
}

bool ngtqg_search_index(NGTQGIndex index, NGTQGQuery query, NGTObjectDistances results, NGTError error) {
  if(index == NULL || query.query == NULL || results == NULL){
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : parametor error: index = " << index << " query = " << query.query << " results = " << results;
    operate_error_string_(ss, error);
    return false;
  }
  
  NGTQG::Index* pindex = static_cast<NGTQG::Index*>(index);   
  int32_t dim = pindex->getObjectSpace().getDimension();

  NGT::Object *ngtquery = NULL;

  if(query.radius < 0.0){
    query.radius = FLT_MAX;
  }

  try{
    std::vector<float> vquery(&query.query[0], &query.query[dim]);
    ngtqg_search_index_(pindex, vquery, query, results);
  }catch(std::exception &err) {
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);
    if(ngtquery != NULL){
      pindex->deleteObject(ngtquery);
    }
    return false;
  }
  return true;
}

void ngtqg_initialize_quantization_parameters(NGTQGQuantizationPrameters *parameters) {
  parameters->dimension_of_subvector = 0;
  parameters->max_number_of_edges = 128;
}

void ngtqg_quantize(const char *indexPath, NGTQGQuantizationPrameters parameters, NGTError error) {
  try{
    NGTQG::Index::quantize(indexPath, parameters.dimension_of_subvector, parameters.max_number_of_edges);
  }catch(std::exception &err){
    std::stringstream ss;
    ss << "Capi : " << __FUNCTION__ << "() : Error: " << err.what();
    operate_error_string_(ss, error);    
    return;
  }
}

