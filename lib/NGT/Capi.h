//
// Copyright (C) 2015-2019 Yahoo Japan Corporation
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

typedef unsigned int ObjectID;
typedef void* NGTIndex;
typedef void* NGTProperty;
typedef void* NGTObjectSpace;
typedef void* NGTObjectDistances;
typedef void* NGTError;

typedef struct {
  ObjectID id;
  float distance;
} NGTObjectDistance;

NGTIndex ngt_open_index(const char *, NGTError);

NGTIndex ngt_create_graph_and_tree(const char *, NGTProperty, NGTError);

NGTProperty ngt_create_property(NGTError);

bool ngt_save_index(const NGTIndex, const char *, NGTError);

bool ngt_get_property(const NGTIndex, NGTProperty, NGTError);

int32_t ngt_get_property_dimension(NGTProperty, NGTError);

bool ngt_set_property_dimension(NGTProperty, int32_t, NGTError);

bool ngt_set_property_edge_size_for_creation(NGTProperty, int16_t, NGTError);

bool ngt_set_property_edge_size_for_search(NGTProperty, int16_t, NGTError);

int32_t ngt_get_property_object_type(NGTProperty, NGTError);

bool ngt_is_property_object_type_float(int32_t);

bool ngt_is_property_object_type_integer(int32_t);

bool ngt_set_property_object_type_float(NGTProperty, NGTError);

bool ngt_set_property_object_type_integer(NGTProperty, NGTError);

bool ngt_set_property_distance_type_l1(NGTProperty, NGTError);

bool ngt_set_property_distance_type_l2(NGTProperty, NGTError);

bool ngt_set_property_distance_type_angle(NGTProperty, NGTError);

bool ngt_set_property_distance_type_hamming(NGTProperty, NGTError);

bool ngt_set_property_distance_type_cosine(NGTProperty, NGTError);

NGTObjectDistances ngt_create_empty_results(NGTError);

bool ngt_search_index(NGTIndex, double*, int32_t, size_t, float, float, NGTObjectDistances, NGTError);
  
int32_t ngt_get_size(NGTObjectDistances, NGTError); // deprecated
  
uint32_t ngt_get_result_size(NGTObjectDistances, NGTError); 

NGTObjectDistance ngt_get_result(const NGTObjectDistances, const uint32_t, NGTError);

ObjectID ngt_insert_index(NGTIndex, double*, uint32_t, NGTError);

ObjectID ngt_append_index(NGTIndex, double*, uint32_t, NGTError);

bool ngt_batch_append_index(NGTIndex, float*, uint32_t, NGTError);

bool ngt_create_index(NGTIndex, uint32_t, NGTError);

bool ngt_remove_index(NGTIndex, ObjectID, NGTError);

NGTObjectSpace ngt_get_object_space(NGTIndex, NGTError);

float* ngt_get_object_as_float(NGTObjectSpace, ObjectID, NGTError);

uint8_t* ngt_get_object_as_integer(NGTObjectSpace, ObjectID, NGTError);

void ngt_destroy_results(NGTObjectDistances);

void ngt_destroy_property(NGTProperty);

void ngt_close_index(NGTIndex);

int16_t ngt_get_property_edge_size_for_creation(NGTProperty, NGTError);

int16_t ngt_get_property_edge_size_for_search(NGTProperty, NGTError);

int32_t ngt_get_property_distance_type(NGTProperty, NGTError);

NGTError ngt_create_error_object();

const char *ngt_get_error_string(const NGTError);

void ngt_clear_error_string(NGTError);
  
void ngt_destroy_error_object(NGTError);
  

#ifdef __cplusplus
}
#endif
