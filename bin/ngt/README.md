NGT
===

Neighborhood Graph and Tree for Indexing High-dimensional Data

[Home](/README.md) / [Build](/README.md#build) / [Command](/bin/ngt/README.md#command) / [License](/README.md#license) / [Publications](/README.md#publications) / [About Us](http://research-lab.yahoo.co.jp/)

Command
=======



**ngt** - proximity search for high dimensional data



      $ ngt command [option] index [data]
        
**Note:**

When the environment variable POSIXLY_CORERECT is set on some platforms such as Cygwin, you should specifiy options 
before the command as follows.

      $ ngt [option] command index [data]



**ngt** provides high-speed nearest neighbor searches against a large volume of data (several million to several 10 million items of data) in high dimensional vector data space (several ten to several thousand dimensions).

**command** is one of:

-   *[create](#create)*
-   *[append](#append)*
-   *[search](#search)*
-   *[remove](#remove)*
-   *[prune](#prune)*
-   *[reconstruct graph](#reconstruct-graph)*

### CREATE

Constructs the specified index with the specified data.

      $ ngt create -d no_of_dimensions [-p no_of_threads] [-b no_of_batch_processes] 
          [-i index_type] [-g graph_type] [-t edge_reduction_threshold] 
          [-e search_range_coefficient] [-E no_of_edges] [-S no_of_edges_at_search_time] 
          [-o object_type] [-D distance_function] [-n no_of_registration_data] 
          index registration_data
        

*index*  
Specifies the name of the directory for the index to be generated. After data registration, the directory consists of multiple files for the index.

*registration\_data*  
Specifies the vector data to be registered. These data shall consist of one object (data item) per line and each dimensional element shall be delimited by a space or tab.

**-d** *no\_of\_dimensions*  
Specifies the number of dimensions of registration data. Specification is unnecessary if each row of the registration data file consists only of dimensional elements. However, if attribute information or other types of data follow the dimensional elements, such subsequent data will be ignored based on the number of dimensions specified here.

**-p** *no\_of\_threads* (default = 24, recomended value = number of cores)   
Specifies the number of threads to be used for parallel processing at generation time.

**-b** *no\_of\_batch\_processes* (default = recomended value = 200)  
Specifies the number of objects for batch processing, which is normally performed for a fixed number of objects. It is generally not necessary to specify this option.

**-i** *index\_type* (__g__|__t__)  
Selects between creation of graph only or graph and tree.
- __g__: Generates only a graph index. At search time, the search start node in the graph is randomly determined.
- __t__: Generates a tree index in addition to a graph index. At search time, the search start node in the graph is determined using the tree. (default/recommended)

**-g** *graph\_type*  
Specifies the type of graph index.
- __a__: Generates ANNG. Enables high-speed registration. (default/recommended)
- __k__: Generates KNNG. Registration is very slow. (for experimental use&mdash;not recommended)
- __b__: Generates BKNNG. Registration is very slow, but search has high level of performance. (for experimental use&mdash;not recommended)

**-t** *edge\_reduction\_threshold* (default = recomended value = 0)  
Specifies the increase in number of edges that acts as a criterion for executing excess-edge reduction processing. Although excess-edge reduction processing can only be used when selecting ANNG, it involves heavy processing and is essential unnecessary unless there is a need to reduce the amount of consumed memory as much as possible. Specifying 0 here prevents the execution of excess-edge reduction processing.

**-e** *search\_range\_coefficient* (default = recomended value = 0.1)  
When specifying ANNG or BKNNG, neighboring nodes connected to an item of registration data (node) by edges are obtained by searching and combined by edges. This option specifies the magnification coefficient of the search range at search time.

**-E** *no\_of\_edges* (default = recomended value = 10)  
Specifies the number of initial edges of each node at graph generation time. Once an index has been generated, the number of edges will be equal to or greater than this specified number in the case of ANNG or BKNNG and equal to this specified number in the case of KNNG.

**-S** *no\_of\_edges\_at\_search\_time* (default = recomended value = 40)  
Specifies the number of edges at search time accompanying or following index generation. This value is used when not specifying the number of edges by the search command. It is specified to conduct searches by a number of edges less than the actual number of edges of each node in the graph. Since a large number of edges may be generated in the case of ANNG or BKNNG, limiting the number of edges can help improve search performance. Specifying 0 here indicates that the number of edges is not to be limited (that all actual edges are to be used). If 0 is specified, index generation is relatively slow, but top performance can be obtained at search time.

**-o** *object\_type*  
Specifies the data object type.
- __c__: 1 byte unsigned integer
- __f__: 4 byte floating point number (default)

**-D** *distance\_function*  
Specifies the distance function as follows.
- __1__: L1 distance
- __2__: L2 distance (default)
- __a__: Angle distance
- __A__: Normalized angle distance. The specified data are automatically normalized to be appended to the index.
- __c__: Cosine similarity
- __C__: Normalized cosine similarity. The specified data are automatically normalized to be appended to the index.
- __h__: Hamming distance. 1 byte unsigned integer should be specified for the data object type.

**-n** *no\_of\_registration\_data*  
Specifies the number of data items to be registered. If not specified, all data in the specified file will be registered.

### APPEND

Adds the specified data to the specified index.

      $ ngt append [-p no_of_threads] [-d no_of_dimensions] [-n no_of_registration_data]
          index registration_data
        

*index*  
Specifies the name of the existing index.

*registration\_data*  
Specifies the vector data to be registered. These data shall consist of one object (data item) per line and each dimensional element shall be delimited by a space or tab.

**-p** *no\_of\_threads* (default = recomended value = 24)   
Specifies the number of threads to be used for parallel processing at generation time.

**-d** *no\_of\_dimensions*  
Specifies the number of dimensions of registration data. Specification is unnecessary if each row of the registration data file consists only of dimensional elements. However, if attribute information or other types of data follow the dimensional elements, such subsequent data will be ignored based on the number of dimensions specified here.

**-n** *no\_of\_registration\_data*  
Specifies the number of data items to be registered. If not specified, all data in the specified file will be registered.

### SEARCH

Searches the index using the specified query data.

      $ ngt search [-i index_type] [-e search_range_coefficient] [-n no_of_search_results] 
          [-E max_no_of_edges] [-r search_radius] index query_data
        

*index*  
Specifies the name of the existing index.

*query\_data*  
Specifies the name of the file containing query data. This file shall consist of one item of query data per line and each dimensional element of that data item shall be delimited by a space or tab the same as registration data. Each search shall be sequentially performed when providing multiple queries.

**-i** *index\_type* (__g__|__t__|__s__)  
Selects between using of graph only or graph and tree.
- __g__: Uses only a graph index. May be specified even if a tree exists. The search start point in the graph is randomly determined.
- __t__: Uses a tree index in addition to a graph index. An error occurs if no tree exists. The search start point in the graph is determined using the tree. Enables high-speed searches compared with graph only. (default/recommended)
- __s__: Performs a linear search using no index.

**-e** *search\_range\_coefficient* (default = recomended value = 0.1)  
Specifies the magnification coefficient of the search range. A larger value means greater accuracy but slower searching, while a smaller value means a drop in accuracy but faster searching. While it is desirable to adjust this value within the range of 0 - 0.3, a negative value may also be specified.

**-n** *no\_of\_search\_results* (default: 20)  
Specifies the number of search results.

**-E** *max\_no\_of\_edges* (default = value specified with the create command or 40; recomended value = 40)   
Specifies the maximum number of edges to be used in the search. This option is specified when conducting a search with fewer edges than the number of edges of each node on the graph. Since a large number of edges can be generated in the case of ANNG or BKNNG, limiting the number of edges in this way tends to improve search performance. Specifying zero here indicates no limiting of number of edges (use all actual edges).

**-r** *search\_radius* (default = infinite circle)  
Specifies the search range in terms of the radius of a circle.

### REMOVE

Removes the specified object from the index.

      $ ngt remove [-d object_id_specification_method] index object_id
        

*index*  
Specifies the name of the existing index.

*object\_id*  
Specifies the ID of the object to be removed.

**-d** *object\_id\_specification\_method* (__f__|__d__) (default = f)  
Specifies the method for specifying the ID of the object to be removed. Specifying __f__ indicates that the following object-ID specification is to be treated as a file name. That file shall consist of one entry per line, each indicting the ID of an object to be removed. Specifying __d__ indicates that the following object-ID specification is to be treated simply as an object-ID referring to the object to be removed.

### PRUNE

Prunes long edges in the graph of the index. Although this command shortens the query time, to further shorten the query time, the path adjustment of the following command reconstruct graph is recommended.

      $ ngt prune -e no_of_forcedly_pruned_edges -s no_of_selectively_pruned_edge index

*index*  
Specifies the name of the existing index.

**-e** *no\_of\_forcedly\_pruned\_edges*   
Specifies the maximum number of edges for each node. Excess edges over the specified number are removed for each node in descending of edge length.

**-s** *no\_of\_selectively\_pruned\_edge*   
Specifies the number (rank) of edges that should not be removed for each node. If the rank of an edge in ascending order of edge length for each node exceeds the specified number and the edge has alternative paths, the edge is removed.

no\_of\_forcedly\_pruned\_edges should be greater than no\_of\_selectively\_pruned\_edges.

### RECONSTRUCT GRAPH

constructs the index with the reconstructed graph from the specified index.

      $ ngt reconstruct-graph [-m mode] -o no_of_original_edges -i no_of_reverse_edge input_index reconstructed_index

*input_index*  
Specifies the name of the existing index.

*reconstructed_index*  
Specifies the name of the reconstructed index.

**-o** *no_of_original_edges*   
Specifies the number of edges for each node to add to the reconstructed graph from the original graph. The specified number also means the lower bound of the outdegrees of the reconstructed graph.

**-i** *no_of_reverse_edges*  
Specifies the number of edges for each node to add to the reconstructed graph from the original graph. Unlike *no_of_original_edges*, after the direction of the edges are reveresed, the edges are added to the reconstructed graph. The specified number also means the lower bound of the indegrees of the reconstructed graph.

**-m** *mode*   
Specified the mode of the path adjustment.
- __s__: no path adjustment.
- __S__: path adjustment.




### Create

Construct an index for 128-dimensional, 1-byte-integer data:

      $ cd (NGT_TOP_DIR)
      $ ngt create -d 128 -o c index ./data/sift-dataset-5k.tsv
      Data loading time=0.160748 (sec) 160.748 (msec)
      # of objects=5000
      Index creation time=0.379659 (sec) 379.659 (msec)

### Search

Perform a neighborhood search by three queries specified in a file:

      $ cd (NGT_TOP_DIR)
      $ ngt search -n 20 index ./data/sift-query-3.tsv
      Query No.1
      Rank	ID	Distance
      1	3031	239.332
      2	4079	240.002
      3	3164	244.504
      4	3718	246.763
      5	157	251.094
      6	2422	251.185
      7	1313	251.34
      8	379	252.446
      9	3521	260.158
      10	2594	261.132
      11	4627	262.381
      12	2159	263.471
      13	3519	264.909
      14	1764	265.136
      15	4400	266.156
      16	2717	266.914
      17	3168	269.637
      18	4236	270.673
      19	4700	272.725
      20	679	272.973
      Query Time= 0.000472 (sec), 0.472 (msec)
      Query No.2
      Rank	ID	Distance
      1	2726	291.983
      2	924	296.987
      3	3638	298.585
      4	858	300.376
      5	1453	306.805
      6	174	307.789
      7	2992	308.485
      8	2980	308.93
      9	1525	309.49
      10	244	309.816
      11	910	310.446
      12	3310	310.585
      13	2433	311.482
      14	1633	311.735
      15	3761	312.44
      16	407	313.252
      17	4546	313.876
      18	697	315.108
      19	34	315.563
      20	2189	316.193
      Query Time= 0.000478 (sec), 0.478 (msec)
      Query No.3
      Rank	ID	Distance
      1	762	194.286
      2	1046	212.695
      3	4906	215.244
      4	2905	216.539
      5	4142	219.479
      6	1879	219.616
      7	4398	223.352
      8	3842	223.468
      9	233	224.127
      10	2794	224.366
      11	2476	224.804
      12	1848	225.803
      13	3364	226.561
      14	4098	226.74
      15	3023	228.884
      16	4113	229.325
      17	1036	232.852
      18	1740	233.144
      19	2302	233.818
      20	2440	233.91
      Query Time= 0.00018 (sec), 0.18 (msec)
      Average Query Time= 0.000376667 (sec), 0.376667 (msec), (0.00113/3)



#### [ONNG](/README.md#onng)
```
$ ngt create -i t -g a -S 0 -e 0.0 -E no_of_edges -d dimensionality_of_data -o data_type -D distatnce_type anng-index vector-data.dat
$ ngt reconstruct-graph -m S -o outdegree -i indegree anng-index onng-index
```
e.g.  
```
$ ngt create -i t -g a -S 0 -e 0.0 -E 100 -d 128 -o c -D 2 anng-index vector-data.dat
$ ngt reconstruct-graph -m S -o 10 -i 120 anng-index onng-index
```
#### [PANNG](/README.md#panng)
```
$ ngt create -i g|t -g a -S 0 -e 0.1 -E no_of_edges -d dimensionality_of_data -o data_type -D distatnce_type panng-index vector-data.dat
$ ngt prune -e no_of_forcedly_pruned_edges -s no_of_selectively_pruned_edges panng-index
```
e.g.  
```
$ ngt create -i t -g a -S 0 -e 0.1 -E 10 -d 128 -o c -D 2 panng-index vector-data.dat
$ ngt prune -e 60 -s 30 panng-index
```
#### [ANNGT](/README.md#anngt)
```
$ ngt create -i t -g a -S 0 -e 0.1 -E no_of_edges(k) -d dimensionality_of_data -o data_type -D distance_type anngt-index vector-data.dat
```
e.g.
```  
$ ngt create -i t -g a -S 0 -e 0.1 -E 16 -d 128 -o c -D 2 anngt-index vector-data.dat
```
#### [ANNG](/README.md#anng)
```
$ ngt create -i g -g a -S 0 -e 0.1 -E no_of_edges(k) -d dimensionality_of_data -o data_type -D distance_type anng-index vector-data.dat
```
e.g.
```
$ ngt create -i g -g a -S 0 -e 0.1 -E 16 -d dimensionality_of_data -o data_type -D distance_type anng-index vector-data.dat
```
#### KNNG  
```
$ ngt create -i g -g k -S 0 -E no_of_edges(k) -d dimensionality_of_data -o data_type -D distance_type knng-index vector-data.dat
```
e.g.
```  
$ ngt create -i g -g k -S 0 -E 20 -d 128 -o c -D 2 knng-index vector-data.dat
```
