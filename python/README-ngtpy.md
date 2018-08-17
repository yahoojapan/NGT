ngtpy
=====

NGT python

Class Index
===========

## Member Functions

### \_\_init\_\_
Opens the specified index and creates the index object for the index.

      __init__(self: ngtpy.Index, path: str, read_only: bool=False, zero_based_numbering: bool=True)

**Returns**  
None.

**path**   
Specifies the path of the index to open.

**read\_only**   
Opens the index as read-only. False means opening as read-write.

**zero\_based\_numbering**   
Specifies zero-based numbering for object IDs. False means one-based numbering.

### close
Closes the index.

      close(self: ngtpy.Index)

**Returns**   
None.

### insert
Inserts the specified object into the index. Not build the index for the object. The function build_index below should be called to build the index after inserting objects by using this function.

      int insert(self: ngtpy.Index, object: numpy.ndarray[float64])

**Returns**   
ID for the inserted object.

**object**   
Specifies the inserted object.

### build_index
Build the index for the objects that have been inserted by using the function insert.

      build_index(self: ngtpy.Index, num_threads: int=8)

**Returns**   
None.

**num_thread**  
Specifies the number of threads to build the index.

### batch_insert
Inserts the specified objects and builds the index for the objects. This function is almost the same as executing the ngt command "ngt append". You may execute the ngt command instead of calling this function.

      batch_insert(self: ngtpy.Index, objects: numpy.ndarray[float64], num_threads: int=8)

**Returns**  
None.

**objects**   
Specifies the inserted objects.

**num_thread**   
Specifies the number of threads to build the index.


### remove
Removes the specified object.

      remove(self: ngtpy.Index, object_id: int)

**Returns**  
None.

**object\_id**   
Specifies the removed object ID.

### save
Saves the index.

      save(self: ngtpy.Index)

**Returns**  
None.

### get_object
Gets the specified object.

      List[float] get_object(self: ngtpy.Index, object_id: int)

**Returns**   
The specified object.


### search
Searches the nearest objects to the specified query object.

      object search(self: ngtpy.Index, query: object, size: int=10, epsilon: float=0.1, edge_size: int=-1, with_distance: bool=True)

**Returns**   
The list of tuples(object ID, distance) as the search result. 

**query**   
Specifies the query object.

**size**   
Specifies the number of the objects as the search result.

**epsilon**   
Specifies epsilon which defines the explored range for the graph.

**edge\_size**   
Specifies the number of edges for each node to explore the graph.

**with\_distance**   
Specifies object IDs with distances as the result. False means that the result is a list of only object IDs.

FUNCTIONS
=========

### create
Creates an empty index. This function is almost the same as executing the ngt command "ngt create" with an empty inserted object file. You may execute the ngt command instead of calling this function.

      create(path: str, dimension: int, edge_size_for_creation: int=10, edge_size_for_search: int=40, distance_type: str='L2', object_type: str='Float')


**Returns**   
None.

**path**   
Specifies the path of the index.

**dimension**  
Specifies the dimensionality of the inserted object.

**edge\_size\_for\_creation**   
Specifies the initial number of edges for each node.

**edge\_size\_for\_search**   
Specifies the number of edges for each node to explore the graph for the search processing.

**distance\_type**   
Specifies the distance function for the objects.
- __L1__: L1 distance
- __L2__: L2 distance (default)
- __Angle__: Angle distance
- __Normalized Angle__: Normalized angle distance. The specified data are automatically normalized to be appended to the index.
- __Cosine__: Cosine similarity
- __Normalized Cosine__: Normalized cosine similarity. The specified data are automatically normalized to be appended to the index.
- __Hamming__: Hamming distance

**object\_type**  
Specifies the data type of the objects.
- __Float__: 4 byte floating point number
- __Byte__: 1 byte unsigned integer

