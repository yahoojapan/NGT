QBG
===

Command-line interface for NGT with Quantization for indexing high-dimensional data

Command
=======



**qbg** - proximity search for high dimensional data with quantization



      $ qbg command [option] index [data]
        
**Note:**

When the environment variable POSIXLY_CORERECT is set on some platforms such as Cygwin, you should specifiy options 
before the command as follows.

      $ qbg [option] command index [data]



**qbg** handles two types of graphs with quantization: Quantized Graph (QG) and Quantized Blob Graph (QBG).

**command** for the quantized graph is one of:

-   *[create-qg](#create-qg)*
-   *[build-qg](#build-qg)*
-   *[search-qg](#search-qg)*

**command** for the quantized blob graph is one of:

-   *[create](#create)*
-   *[append](#append)*
-   *[build](#build)*
-   *[search](#search)*

### CREATE-QG
Make and initialize a QG directory for the quantized graph in the specified NGT index directory, and insert the data in the NGT index into the QG index.

      $ qbg create-qg [-P number_of_extended_dimensions] [-Q number_of_subvector_dimensions] index

*index*  
Specify the name of the directory for the existing index such as ANNG or ONNG to be quantized. The index only with L2 distance and normalized cosine similarity distance can be quantized. You should build the ANNG or ONNG with normalized cosine similarity in order to use cosine similarity for the quantized graph.

**-P** *number_of_extended_dimensions*  
Specify the number of the extended dimensions. The number should be greater than or equal to the number of the genuine dimensions, and also should be a multiple of 4. When this option is not specified, the smallest multiple of 4 that is greater than the dimension is set to the number of the extended dimensions.

**-Q** *number_of_subvector_dimension*  
Specify the number of the subvector dimensions. The number should be less than or equal to the the number of the extended dimensions, and also should be a divisor of the number of the extended dimensions. When this option is not specified, one is set to the number of the subvector dimensions.


### BUILD-QG

Quantize the objects of the specified index and build a quantized graph into the index.

      $ qbg build-qg [-o number_of_objects_for_quantization] [-E max_number_of_edges] [-M number_of_trials] index

*index*  
Specify the name of the directory for the existing index such as ANNG or ONNG to be quantized. The index only with L2 distance and normalized cosine similarity distance can be quantized. You should build the ANNG or ONNG with normalized cosine similarity in order to use cosine similarity for the quantized graph.

**-o** *number_of_objects_for_quantization*  
Specify the number of object for quantization and optimization. The number should be less than or equal to the number of the registered objects.

**-E** *max_number_of_edges*  
Specify the maximum number of edges to build a qunatized graph. Since every 16 objects that are associated with edges of each node are processed, the number should be a multiple of 16.

**-M** *number_of_trials*  
Specify the number of trials to optimize the subvector quantization.

### SEARCH-QG

Search the index using the specified query data.

      $ qbg search-qg [-n number_of_search_objects] [-e search_range_coefficient] [-p result_expansion]
          [-r search_radius] index query_data
        

*index*  
Specify the path of the existing quantized index.

*query_data*  
Specify the path of the file containing query data. This file shall consist of one item of query data per line and each dimensional element of that data item shall be delimited by a space or tab. Each search shall be sequentially performed when providing multiple queries.

**-n** *number_of_search_objects* (default: 20)  
Specify the number of search objects.

**-e** *search_range_coefficient* (default = 0.02)  
Specify the magnification coefficient (epsilon) of the search range. A larger value means higher accuracy but slower searching, while a smaller value means a drop in accuracy but faster searching. While it is desirable to adjust this value within the range of 0 - 0.1, a negative value (> -1.0) may also be specified.

**-p** *result_expansion* (default = 3.0)   
Specify the expansion ratio of the number of approximate inner search objects to the number of search objects. For example, when the ratio is 10 and the number of search objects is 20, the number of the approximate search objects is set to 200 inside the search processing. A larger value brings higher accuracy but slower searching.

### CREATE
Make and initialize a QBG directory for the quantized blob graph.

      $ qbg create [-d number_of_dimension] [-P number_of_extended_dimensions] [-O object_type] [-D distance_function] [-C number_of_blobs] index

*index*  
Specify the name of the directory for QBG.

**-d** *number_of_dimensions*  
Specify the number of dimensions of registration data.

**-P** *number_of_extended_dimensions*  
Specify the number of the extended dimensions. The number should be greater than or equal to the number of the genuine dimensions, and also should be a multiple of 4. When this option is not specified, the smallest multiple of 4 that is greater than the dimension is set to the number of the extended dimensions.

**-O** *object_type*  
Specify the data object type.
- __c__: 1 byte unsigned integer
- __f__: 4 byte floating point number (default)

**-D** *distance_function*  
Specify the distance function as follows.
- __2__: L2 distance (default)
- __c__: Cosine similarity

**-C** *number_of_blobs*  
Specify the number of blobs that should be less than or equal to the number of quantization clusters.

**-N** *number_of_subvectors*  
Specify the number of subvectors that should be a divisor of the number of the extended dimensions.

### APPEND

Append the specified data to the specified index.

      $ qbg append index registration_data

### BUILD

Quantize the objects of the specified index and build a quantized graph into the index.

      $ qbg build [-o number_of_objects_for_quantization] [-E max_number_of_edges] [-M number_of_trials] [-P rotation] index

*index*  
Specify the name of the directory for the existing index such as ANNG or ONNG to be quantized. The index only with L2 distance and normalized cosine similarity distance can be quantized. You should build the ANNG or ONNG with normalized cosine similarity in order to use cosine similarity for the quantized graph.

**-o** *number_of_objects_for_quantization*  
Specify the number of object for quantization and optimization. The number should be less than or equal to the number of the registered objects.

**-P** *rotation*  
Specify the transform matrix type for the inserted and query object to optimize the subvector quantization.
- __r__: Rotation matrix.
- __R__: Rotation and repositioning matrix.
- __p__: Repositioning matrix.
- __n__: No matrix.

**-M** *number_of_trials*  
Specify the number of trials to optimize the subvector quantization.

### SEARCH

Search the index using the specified query data.

      $ qbg search [-n number_of_search_objects] [-e search_range_coefficient] [-p result_expansion]
          index query_data
        

*index*  
Specify the path of the existing quantized index.

*query_data*  
Specify the path of the file containing query data. This file shall consist of one item of query data per line and each dimensional element of that data item shall be delimited by a space or tab. Each search shall be sequentially performed when providing multiple queries.

**-n** *number_of_search_objects* (default: 20)  
Specify the number of search objects.

**-e** *search_range_coefficient* (default = 0.02)  
Specify the magnification coefficient (epsilon) of the search range. A larger value means higher accuracy but slower searching, while a smaller value means a drop in accuracy but faster searching. While it is desirable to adjust this value within the range of 0 - 0.1, a negative value (> -1.0) may also be specified.

**-B** *blob_search_range_coefficient* (default = 0.0)  
Specify the magnification coefficient (epsilon) of the search range for the quantized blob graph.

**-N** *number_of_explored_nodes* (default = 256)  
Specify the number of the explored nodes in the graph. When the number of the explored nodes reached the specified number, the search is terminated. 

**-p** *result_expansion* (default = 0.0)   
Specify the expansion ratio of the number of approximate inner search objects to the number of search objects. For example, when the ratio is 10 and the number of search objects is 20, the number of the approximate search objects is set to 200 inside the search processing. A larger value brings higher accuracy but slower searching.


Examples of using the quantized graph
-------------------------------------

### Setup data

      $ curl -L -O https://github.com/yahoojapan/NGT/raw/main/tests/datasets/ann-benchmarks/sift-128-euclidean.tsv
      $ curl -L -O https://github.com/yahoojapan/NGT/raw/main/tests/datasets/ann-benchmarks/sift-128-euclidean_query.tsv
      $ head -1 sift-128-euclidean_query.tsv > query.tsv

### Build the quantized graph

Build an ANNG for 128-dimensional, floating point data:

      $ ngt create -d 128 -o f -D 2 anng sift-128-euclidean.tsv
      Data loading time=15.4804 (sec) 15480.4 (msec)
      # of objects=1000000
      Processed 100000 objects. time= 4.26452 (sec)
      ...
      Processed 1000000 objects. time= 7.06745 (sec)
      Index creation time=63.3504 (sec) 63350.4 (msec)

Create and initialize the quantized graph:

      $ qbg create-qg anng
      creating...
      appending...

Build the quantized graph:

      $ qbg build-qg anng
      optimizing...
      building the inverted index...
      building the quantized graph... 

### Search with the quantized graph

Search k nearest neighbors with the quantized graph:

      $ qbg search-qg -n 20 -e 0.02 anng query.tsv
      Query No.1
      Rank	ID	Distance
      1	932086	232.871
      2	934877	234.715
      3	561814	243.99
      ...
      20	2177	276.781
      Query Time= 0.0005034 (sec), 0.5034 (msec)
      Average Query Time= 0.0005034 (sec), 0.5034 (msec), (0.0005034/1)

Examples of building the quantized graph for higher performance
------------------------------------------------------------

Build an ANNG having more edges for higher performance:

      $ ngt create -d 128 -o f -D 2 -E 40 anng-40 sift-128-euclidean.tsv

Build an ONNG:

      $ ngt reconstruct-graph -m S -E 64 -o 64 -i 120 anng-40 onng-40

Create and initialize the quantized graph:

      $ qbg create-qg onng-40

Build the quantized graph:

      $ qbg build-qg onng-40

Search k nearest neighbors with the quantized graph:

      $ qbg search -n 20 -e 0.02 onng-40 query.tsv


Examples of using the quantized blob graph
-------------------------------------

### Setup data

      $ curl -L -O https://github.com/yahoojapan/NGT/raw/main/tests/datasets/ann-benchmarks/sift-128-euclidean.tsv
      $ curl -L -O https://github.com/yahoojapan/NGT/raw/main/tests/datasets/ann-benchmarks/sift-128-euclidean_query.tsv
      $ head -1 sift-128-euclidean_query.tsv > query.tsv

### Build the quantized blob graph

Create and initialize the quantized blob graph:

      $ qbg create -d 128 -D 2 -N 128 qbg-index

Append objects:

      $ qbg append qbg-index sift-128-euclidean.tsv

Build the quantized graph:

      $ qbg build qbg-index

### Search with the quantized blob graph

Search k nearest neighbors with the quantized blob graph:

      $ qbg search -n 20 -e 0.02 qbg-index query.tsv
      

