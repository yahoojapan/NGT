NGTQG
===

Neighborhood Graph and Tree for Indexing High-dimensional Data with Quantized Graph

Command
=======

## Name

**ngtqg** - proximity search for high dimensional data with quantized graph

## Synopsis

      $ ngtqg command [option] index [data]
        
**Note:**

When the environment variable POSIXLY_CORERECT is set on some platforms such as Cygwin, you should specifiy options 
before the command as follows.

      $ ngtqg [option] command index [data]

## Description

**ngtqg** provides high-speed nearest neighbor searches against high dimensional data.

**command** is one of:

-   *[quantize](#quantize)*
-   *[search](#search)*

### QUANTIZE

Quantize the objects of the specified index and build a quantized graph into the index.

      $ ngtqg quantize [-E max_no_of_edges] index

*index*  
Specify the name of the directory for the existing index such as ANNG or ONNG to be quantized. The index only with L2 distance and normalized cosine similarity distance can be quantized. You should build the ANNG or ONNG with normalized cosine similarity in order to use cosine similarity for the quantized graph.

**-E** *max_no_of_edges*  
Specify the maximum number of edges to build a qunatized graph. Since every 16 objects that are associated with edges of each node are processed, the number should be a multiple of 16.

### SEARCH

Search the index using the specified query data.

      $ ngtqg search [-n no_of_search_objects] [-e search_range_coefficient] [-p result_expansion]
          [-r search_radius] index query_data
        

*index*  
Specify the path of the existing quantized index.

*query_data*  
Specify the path of the file containing query data. This file shall consist of one item of query data per line and each dimensional element of that data item shall be delimited by a space or tab. Each search shall be sequentially performed when providing multiple queries.

**-n** *no_of_search_objects* (default: 20)  
Specify the number of search objects.

**-e** *search_range_coefficient* (default = 0.02)  
Specify the magnification coefficient (epsilon) of the search range. A larger value means higher accuracy but slower searching, while a smaller value means a drop in accuracy but faster searching. While it is desirable to adjust this value within the range of 0 - 0.1, a negative value (> -1.0) may also be specified.

**-p** *result_expansion*  (default = 3.0)
Specify the expansion ratio of the number of approximate inner search objects to the number of search objects. For example, when the ratio is 10 and the number of search objects is 20, the number of the approximate search objects is set to 200 inside the search processing. A larger value brings higher accuracy but slower searching.

**-r** *search_radius* (default = infinite circle)  
Specify the search range in terms of the radius of a circle.

Examples of using the quantized graph
-------------------------------------

### Setup data

      $ curl -L -O https://github.com/yahoojapan/NGT/raw/master/tests/datasets/ann-benchmarks/sift-128-euclidean.tsv
      $ curl -L -O https://github.com/yahoojapan/NGT/raw/master/tests/datasets/ann-benchmarks/sift-128-euclidean_query.tsv
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

Quantize the objects in the ANNG and build the quantized graph:

      $ ngtqg quantize anng
      Clustering of the subvector is complete. anng/qg/local-0:17
      ...
      Clustering of the subvector is complete. anng/qg/local-127:17
      Processed 100000 objects.
      ...
      Processed 1000000 objects.

### Search with the quantize graph

Search k nearest neighbors with the quantized graph:

      $ ngtqg search -n 20 -e 0.02 anng query.tsv
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

Quantize the objects and build a quantized graph from the ONNG:

      $ ngtqg quantize onng-40

Search k nearest neighbors with the quantized graph:

      $ ngtqg search -n 20 -e 0.02 onng-40 query.tsv

