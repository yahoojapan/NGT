# -*- coding: utf-8 -*-

#
# Copyright (C) 2015 Yahoo Japan Corporation
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from ctypes import (Structure,
                    POINTER,
                    c_char,
                    c_uint,
                    c_uint8,
                    c_void_p,
                    c_char_p,
                    c_int,
                    c_double,
                    c_float,
                    cdll)
import ctypes.util
import csv


class NativeError(Exception):
    pass


class APIError(Exception):
    pass


class Index(object):
    '''
    NGT: Neighborhood Graph and Tree for Indexing High-dimensional Data
      NGT provides the functionality of searching for approximate nearest neighbors in high-dimensional data.

    Example:
    ```python
      from ngt import base as ngt
      import random

      dim = 10
      objects = []
      for i in range(0, 100) :
          vector = random.sample(range(100), dim)
          objects.append(vector)

      query = objects[0]
      index = ngt.Index.create(b"tmp", dim)
      index.insert(objects)
      # You can also insert objects from a file like this.
      # index.insert_from_tsv('list.dat')

      index.save()
      # You can load the saved index like this.
      # index = ngt.Index(b"tmp")

      result = index.search(query, 3)

      for i, o in enumerate(result) :
          print(str(i) + ": " + str(o.id) + ", " + str(o.distance))
          object = index.get_object(o.id)
          print(object)
    ```
    '''
    class ObjectDistance(Structure):
        _fields_ = [
            ('id', c_uint),
            ('distance', c_float)
        ]

        def __repr__(self):
            return "ObjectDistances({},{})".format(
                self.id,
                self.distance)

    __ngt = cdll.LoadLibrary(ctypes.util.find_library("ngt"))

    __ngt.ngt_open_index.argtypes = [c_char_p, c_void_p]
    __ngt.ngt_open_index.restype = c_void_p

    __ngt.ngt_create_graph_and_tree.argtypes = [c_char_p, c_void_p, c_void_p]
    __ngt.ngt_create_graph_and_tree.restype = c_void_p

    __ngt.ngt_create_property.argtypes = [c_void_p]
    __ngt.ngt_create_property.restype = c_void_p

    __ngt.ngt_save_index.argtypes = [c_void_p, c_char_p, c_void_p]

    __ngt.ngt_get_property.argtypes = [c_void_p, c_void_p, c_void_p]

    __ngt.ngt_get_property_dimension.argtypes = [c_void_p, c_void_p]

    __ngt.ngt_set_property_dimension.argtypes = [c_void_p, c_int, c_void_p]

    __ngt.ngt_set_property_edge_size_for_creation.argtypes = [
                                                c_void_p, c_int, c_void_p]

    __ngt.ngt_set_property_edge_size_for_search.argtypes = [
                                                c_void_p, c_int, c_void_p]

    __ngt.ngt_get_property_object_type.argtypes = [c_void_p]

    __ngt.ngt_is_property_object_type_float.argtypes = [c_int]

    __ngt.ngt_is_property_object_type_integer.argtypes = [c_int]

    __ngt.ngt_set_property_object_type_float.argtypes = [c_void_p, c_void_p]

    __ngt.ngt_set_property_object_type_integer.argtypes = [c_void_p, c_void_p]

    __ngt.ngt_set_property_distance_type_l1.argtypes = [c_void_p, c_void_p]

    __ngt.ngt_set_property_distance_type_l2.argtypes = [c_void_p, c_void_p]

    __ngt.ngt_set_property_distance_type_angle.argtypes = [c_void_p, c_void_p]

    __ngt.ngt_set_property_distance_type_hamming.argtypes = [c_void_p, c_void_p]

    __ngt.ngt_set_property_distance_type_jaccard.argtypes = [c_void_p, c_void_p]

    __ngt.ngt_set_property_distance_type_cosine.argtypes = [c_void_p, c_void_p]

    __ngt.ngt_create_empty_results.argtype = [c_void_p]
    __ngt.ngt_create_empty_results.restype = c_void_p

    __ngt.ngt_get_size.argtypes = [c_void_p, c_void_p]

    __ngt.ngt_get_result.argtypes = [c_void_p, c_int, c_void_p]
    __ngt.ngt_get_result.restype = ObjectDistance

    __ngt.ngt_create_index.argtypes = [c_void_p, c_int, c_void_p]

    __ngt.ngt_remove_index.argtypes = [c_void_p, c_uint, c_void_p]

    __ngt.ngt_get_object_space.argtypes = [c_void_p]
    __ngt.ngt_get_object_space.restype = c_void_p

    __ngt.ngt_get_object_as_float.argtypes = [c_void_p, c_int, c_void_p]
    __ngt.ngt_get_object_as_float.restype = POINTER(c_float)

    __ngt.ngt_get_object_as_integer.argtypes = [c_void_p, c_int, c_void_p]
    __ngt.ngt_get_object_as_integer.restype = POINTER(c_uint8)

    __ngt.ngt_destroy_results.argtypes = [c_void_p]

    __ngt.ngt_destroy_property.argtypes = [c_void_p]

    __ngt.ngt_close_index.argtypes = [c_void_p]

    __ngt.ngt_get_property_edge_size_for_creation.argtypes = [
                                                        c_void_p, c_void_p]

    __ngt.ngt_get_property_edge_size_for_search.argtypes = [c_void_p, c_void_p]

    __ngt.ngt_get_property_distance_type.argtypes = [c_void_p, c_void_p]

    __ngt.ngt_create_error_object.restype = c_void_p

    __ngt.ngt_get_error_string.argtypes = [c_void_p]
    __ngt.ngt_get_error_string.restype = c_char_p

    __ngt.ngt_clear_error_string.argtypes = [c_void_p]

    __ngt.ngt_destroy_error_object.argtypes = [c_void_p]

    @staticmethod
    def _check_error_num(status, error):
        if status == 0:
            message = Index.__ngt.ngt_get_error_string(error)
            Index.__ngt.ngt_clear_error_string(error)
            raise NativeError(message)

    @staticmethod
    def _check_error_pnum(status, error):
        if status < 0:
            message = Index.__ngt.ngt_get_error_string(error)
            Index.__ngt.ngt_clear_error_string(error)
            raise NativeError(message)

    @staticmethod
    def _check_error_obj(status, error):
        if status is None:
            message = Index.__ngt.ngt_get_error_string(error)
            Index.__ngt.ngt_clear_error_string(error)
            raise NativeError(message)

    @staticmethod
    def create(path, dimension,
               edge_size_for_creation=10, edge_size_for_search=40,
               object_type="Float", distance_type="L2"):
        '''
        create an empty index with the specified parameters.
          edge_size_for_creation : Number of edges for each node in the graph.
          edge_size_for_search   : Number of edges to search.
          object_type            : Type of the data object. (Float, Integer [Integer is 1 byte])
          distance_type          : Type of the distance function. (L1,L2,Angle,Hamming,Cosine)
        '''
        err = None
        prop = None
        index = None
        try:
            err = Index.__ngt.ngt_create_error_object()

            prop = Index.__ngt.ngt_create_property(err)
            Index._check_error_obj(prop, err)

            stat = Index.__ngt.ngt_set_property_dimension(prop, dimension, err)
            Index._check_error_num(stat, err)

            stat = Index.__ngt.ngt_set_property_edge_size_for_creation(
                                     prop, edge_size_for_creation, err)
            Index._check_error_num(stat, err)
            stat = Index.__ngt.ngt_set_property_edge_size_for_search(
                                     prop, edge_size_for_search, err)
            Index._check_error_num(stat, err)

            if object_type == "Float":
                stat = Index.__ngt.ngt_set_property_object_type_float(
                                                                prop, err)
                Index._check_error_num(stat, err)
            elif object_type == "Integer":
                stat = Index.__ngt.ngt_set_property_object_type_integer(
                                                                prop, err)
                Index._check_error_num(stat, err)

            if distance_type == "L1":
                stat = Index.__ngt.ngt_set_property_distance_type_l1(
                                                                prop, err)
                Index._check_error_num(stat, err)
            elif distance_type == "L2":
                stat = Index.__ngt.ngt_set_property_distance_type_l2(
                                                                prop, err)
                Index._check_error_num(stat, err)
            elif distance_type == "Angle":
                stat = Index.__ngt.ngt_set_property_distance_type_angle(
                                                                prop, err)
                Index._check_error_num(stat, err)
            elif distance_type == "Hamming":
                stat = Index.__ngt.ngt_set_property_distance_type_hamming(
                                                                prop, err)
                Index._check_error_num(stat, err)
            elif distance_type == "Jaccard":
                stat = Index.__ngt.ngt_set_property_distance_type_jaccard(
                                                                prop, err)
                Index._check_error_num(stat, err)
            elif distance_type == "Cosine":
                stat = Index.__ngt.ngt_set_property_distance_type_cosine(
                                                                prop, err)
                Index._check_error_num(stat, err)

            index = Index.__ngt.ngt_create_graph_and_tree(path, prop, err)
            Index._check_error_obj(index, err)

            # close and delete the index object
            Index.__ngt.ngt_close_index(index)
            index = None

            obj = Index(path)

        finally:
            if err is not None:
                Index.__ngt.ngt_destroy_error_object(err)
            if prop is not None:
                Index.__ngt.ngt_destroy_property(prop)
            if index is not None:
                Index.__ngt.ngt_close_index(index)

        return obj

    def __init__(self, path):
        '''
        open the specified index by the path.

            path : Path of the opened index.
        '''
        self.path = path

        # initial object and get dimension
        self.err = self.__ngt.ngt_create_error_object()
        self.index = self.__ngt.ngt_open_index(path, self.err)
        self._check_error_obj(self.index, self.err)

        self.prop = self.__ngt.ngt_create_property(self.err)
        self._check_error_obj(self.prop, self.err)

        stat = self.__ngt.ngt_get_property(self.index, self.prop, self.err)
        self._check_error_num(stat, self.err)

        self.dim = self.__ngt.ngt_get_property_dimension(self.prop, self.err)
        self._check_error_pnum(self.dim, self.err)

        self.otype = self.__ngt.ngt_get_property_object_type(
                                                    self.prop, self.err)
        self._check_error_obj(self.otype, self.err)
        self.is_float = self.__ngt.ngt_is_property_object_type_float(
                                                    self.otype) > 0
        self.ospace = self.__ngt.ngt_get_object_space(self.index, self.err)
        self._check_error_obj(self.ospace, self.err)

        # c method define (depends dimension)
        self.__ngt.ngt_search_index.argtypes = [
                                                c_void_p,
                                                (c_double * self.dim),
                                                c_int,
                                                c_int,
                                                c_float,
                                                c_float,
                                                c_void_p,
                                                c_void_p
                                               ]
        self.__ngt.ngt_insert_index.argtypes = [
            c_void_p, (c_double * self.dim), c_int, c_void_p]
        self.__ngt.ngt_batch_append_index.argtypes = [
                                                c_void_p,
                                                c_void_p,
                                                c_int,
                                                c_void_p
                                               ]
        self.__ngt.ngt_insert_index.restype = c_uint

    def search(self, query, k=20, epsilon=0.1):
        '''
        search for the k nearest neighbors of the specifiecd query object.

            query   : Query object.
            k       : Number of searched objects.
            epsilon : Epsilon defines a search range.
        '''
        try:
            err = self.__ngt.ngt_create_error_object()
            results = self.__ngt.ngt_create_empty_results(err)
            self._check_error_obj(results, err)

            cvec = (c_double * len(query))(*query)
            stat = self.__ngt.ngt_search_index(
                self.index, cvec, self.dim, k, epsilon, -1.0, results, err)
            self._check_error_num(stat, err)
            rsize = self.__ngt.ngt_get_size(results, err)
            self._check_error_pnum(rsize, err)

            ret = []
            for i in range(rsize):
                result = self.__ngt.ngt_get_result(results, i, err)
                self._check_error_obj(result, err)
                ret.append(result)

        finally:
            if err is not None:
                self.__ngt.ngt_destroy_error_object(err)
            if results is not None:
                self.__ngt.ngt_destroy_results(results)
        return ret

    def insert_object(self, object):
        '''
        insert the specified object into the index.
        must call build_index after call this method.

            object : Inserted object.
            return : The ID of the inserted object.
        '''
        cvec = (c_double * len(object))(*object)
        id = self.__ngt.ngt_insert_index(self.index, cvec, self.dim, self.err)
        self._check_error_pnum(id, self.err)
        return id

    def insert(self, objects, num_threads=8):
        '''
        insert the specified objects into the index and build the index.

            objects : Inserted objects.
            return  : List of the IDs of the inserted objects.
        '''
        idList = []
        for object in objects:
            idList.append(self.insert_object(object))
        self.build_index(num_threads)
        return idList

    def insert_blob(self, objects, num_threads=8):
        '''
        insert the specified objects into the index and build the index.
        Although this is the same as the fucntion insert(), both implementations are different.

            objects : Inserted objects.
            num_threads : Number of threads in building index.
        '''
        data_count = len(objects)
        if len(objects[0]) != self.dim:
            message = "insert_blob: Inconsistent dimensionality. " \
                + "The expected is {}. The specified is {}".format(
                                         self.dim, len(objects[0]))
            raise APIError(message)
        merged_vectors = []
        for object in objects:
            merged_vectors.extend(object)
        cvec = (c_float * len(merged_vectors))(*merged_vectors)
        stat = self.__ngt.ngt_batch_append_index(
            self.index, cvec, data_count, self.err)
        self._check_error_num(stat, self.err)
        self.build_index(num_threads)
        return

    def insert_from_tsv(self, path, num_threads=8, dlmt='\t'):
        '''
        insert objects in the specified file and build the index.

            path : Path of the object file.
            num_threads : Number of threads in building index.
            dlmt : Delimiter to sepalate each element in the object file.
        '''
        idList = []
        with open(path) as f:
            reader = csv.reader(f, delimiter=dlmt)
            for row in reader:
                object = [float(x) for x in row]
                idList.append(self.insert_object(object))
        self.build_index(num_threads)
        return idList

    def build_index(self, num_threads=8):
        '''
        build the inserted objects into the index.

            num_threads : Number of threads in building index.
        '''
        stat = self.__ngt.ngt_create_index(self.index, num_threads, self.err)
        self._check_error_num(stat, self.err)

    def remove(self, id):
        '''
        remove the specified object by id

            id : Object id.
        '''
        stat = self.__ngt.ngt_remove_index(self.index, id, self.err)
        self._check_error_num(stat, self.err)

    def get_object(self, id):
        '''
        get the specfied object by id.

            id : Object id.
        '''
        vec = None
        try:
            err = self.__ngt.ngt_create_error_object()
            if self.is_float:
                cvec = self.__ngt.ngt_get_object_as_float(self.ospace, id, err)
            else:
                cvec = self.__ngt.ngt_get_object_as_integer(
                                                        self.ospace, id, err)
            vec = []
            try:
                for i in range(self.dim):
                    vec.append(cvec[i])
            except:
                # NULL that NGT cpai returns cannot be handled with python NONE.
                # That is why some lines below try to handle with the error.
                message = Index.__ngt.ngt_get_error_string(err)
                Index.__ngt.ngt_clear_error_string(err)
                raise NativeError(message)

        finally:
            if err is not None:
                self.__ngt.ngt_destroy_error_object(err)

        return vec

    def save(self, path=None):
        '''
        save the index.

            path : Path to save the index. default overwrite the files.
        '''
        if path is None:
            path = self.path
        stat = self.__ngt.ngt_save_index(self.index, path, self.err)
        self._check_error_num(stat, self.err)

    def __del__(self):
        '''
            close this index.
        '''
        if self.err:
            self.__ngt.ngt_destroy_error_object(self.err)
        if self.prop:
            self.__ngt.ngt_destroy_property(self.prop)
        if self.index:
            self.__ngt.ngt_close_index(self.index)
