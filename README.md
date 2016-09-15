NGT
===

Neighborhood Graph and Tree for Indexing High-dimensional Data

[Build](#build) / [Command](#command) / [License](#license) / [Publications](#publications) / [About Us](http://research-lab.yahoo.co.jp/)

**NGT** provides commands and a library for performing high-speed approximate nearest neighbor searches against a large volume of data (several million to several 10 million items of data) in high dimensional vector data space (several ten to several thousand dimensions).

The experimental code used in our [ACL 2016](http://acl2016.org/) paper will be available by the end of October 2016.

#### Download

- [NGT GitHub](https://github.com/yahoojapan/NGT/)

Build
-----------------

      $ unzip NGT-x.x.x.zip
      $ cd NGT-x.x.x
      $ mkdir build
      $ cd build 
      $ cmake ..
      $ make 
      $ make install
      $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib64

### Shared memory use

The index can be placed in shared memory. Using shared memory can reduce the amount of memory needed when multiple processes are using the same index. It can also improve the boot-up speed of an index for a large volume of registration data. Since changes become necessary at build time, please add the following parameter when executing "cmake" in order to use shared memory.

      $ cmake -DNGT_SHARED_MEMORY_ALLOCATOR=ON ..

Note: Since there is no lock function, the index should be used only for reference when multiple processes are using the same index,


Command
-------

### Name

**ngt** - proximity search for high dimensional data

### Synopsis

      $ ngt command [option] index [data]
        

### Description

**ngt** provides high-speed nearest neighbor searches against a large volume of data (several million to several 10 million items of data) in high dimensional vector data space (several ten to several thousand dimensions).

**command** is one of:

-   *[create](#create)*
-   *[append](#append)*
-   *[search](#search)*
-   *[remove](#remove)*

#### CREATE

Constructs the specified index with the specified data.

      $ ngt create -d no_of_dimensions [-p no_of_threads] [-b no_of_batch_processes] 
          [-i index_type] [-g graph_type] [-t edge_reduction_threshold] 
          [-e search_range_coefficient] [-E no_of_edges] [-S no_of_edges_at_search_time] 
          [-o object_type] [-D distance_function] [-n no_of_registration_data] 
          index registration_data
        

*index*  
Specifies the name of the index to be generated. After data registration, an index is generated consisting of multiple files each with an extension attached to this index name.

*registration\_data*  
Specifies the vector data to be registered. These data shall consist of one object (data item) per line and each dimensional element shall be delimited by a space or tab.

**-d** *no\_of\_dimensions*  
Specifies the number of dimensions of registration data. Specification is unnecessary if each row of the registration data file consists only of dimensional elements. However, if attribute information or other types of data follow the dimensional elements, such subsequent data will be ignored based on the number of dimensions specified here.

**-p** *no\_of\_threads* (default = 24, recommend value = number of cores)   
Specifies the number of threads to be used for parallel processing at generation time.

**-b** *no\_of\_batch\_processes* (default = recommend value = 200)  
Specifies the number of objects for batch processing, which is normally performed for a fixed number of objects. It is generally not necessary to specify this option.

**-i** *index\_type* (__g__|__t__)  
Selects between creation of graph only or graph and tree.
- __g__: Generates only a graph index. At search time, the search start point in the graph is randomly determined.
- __t__: Generates a tree index in addition to a graph index. At search time, the search start point in the graph is determined using the tree. (default/recommended)

**-g** *graph\_type*  
Specifies the type of graph index.
- __a__: Generates ANNG. Enables high-speed registration. (default/recommended)
- __k__: Generates KNNG. Registration is very slow. (for experimental use&mdash;not recommended)
- __b__: Generates BKNNG. Registration is very slow, but search has highest level of performance. (for experimental use&mdash;not recommended)

**-t** *edge\_reduction\_threshold* (default = recommend value = 0)  
Specifies the increase in number of edges that acts as a criterion for executing excess-edge reduction processing. Although excess-edge reduction processing can only be used when selecting ANNG, it involves heavy processing and is essential unnecessary unless there is a need to reduce the amount of consumed memory as much as possible. Specifying 0 here prevents the execution of excess-edge reduction processing.

**-e** *search\_range\_coefficient* (default = recommend value = 0.1)  
When specifying ANNG or BKNNG, neighboring nodes connected to an item of registration data (node) by edges are obtained by searching and combined by edges. This option specifies the magnification coefficient of the search range at search time.

**-E** *no\_of\_edges* (default = recommend value = 10)  
Specifies the number of initial edges of each node at graph generation time. Once an index has been generated, the number of edges will be equal to or greater than this specified number in the case of ANNG or BKNNG and equal to this specified number in the case of KNNG.

**-S** *no\_of\_edges\_at\_search\_time* (default = recommend value = 40)  
Specifies the number of edges at search time accompanying or following index generation. This value is used when not specifying the number of edges by the search command. It is specified to conduct searches by a number of edges less than the actual number of edges of each node in the graph. Since a large number of edges may be generated in the case of ANNG or BKNNG, limiting the number of edges can help improve search performance. Specifying 0 here indicates that the number of edges is not to be limited (that all actual edges are to be used). If 0 is specified, index generation is relatively slow, but top performance can be obtained at search time.

**-o** *object\_type*  
Specifies the data object type.
- __c__: 1 byte integer
- __f__: 4 byte floating point (default)

**-D** *distance\_function*  
Specifies the distance function as follows.
- __1__: L1 distance
- __2__: L2 distance (default)
- __a__: angle
- __h__: Hamming distance

**-n** *no\_of\_registration\_data*  
Specifies the number of data items to be registered. If not specified, all data in the specified file will be registered.

#### APPEND

Adds the specified data to the specified index.

      $ ngt append [-p no_of_threads] [-d no_of_dimensions] [-n no_of_registration_data]
          index registration_data
        

*index*  
Specifies name of existing index (the portion excluding index file extension).

*registration\_data*  
Specifies the vector data to be registered. These data shall consist of one object (data item) per line and each dimensional element shall be delimited by a space or tab.

**-p** *no\_of\_threads* (default = recommend value = 24)   
Specifies the number of threads to be used for parallel processing at generation time.

**-d** *no\_of\_dimensions*  
Specifies the number of dimensions of registration data. Specification is unnecessary if each row of the registration data file consists only of dimensional elements. However, if attribute information or other types of data follow the dimensional elements, such subsequent data will be ignored based on the number of dimensions specified here.

**-n** *no\_of\_registration\_data*  
Specifies the number of data items to be registered. If not specified, all data in the specified file will be registered.

#### SEARCH

Searches the index using the specified query data.

      $ ngt search [-i index_type] [-e search_range_coefficient] [-n no_of_searches] 
          [-E max_no_of_edges] [-r search_radius] index query_data
        

*index*  
Specifies name of existing index (the portion excluding index file extension).

*query\_data*  
Specifies the name of the file containing query data. This file shall consist of one item of query data per line and each dimensional element of that data item shall be delimited by a space or tab the same as registration data. A sequential search shall be performed when providing multiple queries.

**-i** *index\_type* (__g__|__t__|__s__)  
Selects between creation of graph only or graph and tree.
- __g__: Uses only a graph index. May be specified even if a tree exists. The search start point in the graph is randomly determined.
- __t__: Uses a tree index in addition to a graph index. An error occurs if no tree exists. The search start point in the graph is determined using the tree. Enables high-speed searches compared with graph only. (default/recommended)
- __s__: Performs a linear search using no index.

**-e** *search\_range\_coefficient* (default = recommend value = 0.1)  
Specifies the magnification coefficient of the search range. A larger value means greater accuracy but slower searching, while a smaller value means a drop in accuracy but faster searching. While it is desirable to adjust this value within the range of 0 - 0.3, a negative value may also be specified.

**-n** *no\_of\_searches* (default: 20)  
Specifies the number of search results.

**-E** *max\_no\_of\_edges* (default = value specified with the create command or 40; recommend value = 40)   
Specifies the maximum number of edges to be used in the search. This option is specified when conducting a search with fewer edges than the number of edges of each node on the graph. Since a large number of edges can be generated in the case of ANNG or BKNNG, limiting the number of edges in this way tends to improve search performance. Specifying zero here indicates no limiting of number of edges (use all actual edges).

**-r** *search\_radius* (default = infinite circle)  
Specifies the search range in terms of the radius of a circle.

#### REMOVE

Removes the specified object from the index.

      $ ngt remove [-d object_id_specification_method] index object_id
        

*index*  
Specifies name of existing index (the portion excluding index file extension).

*object\_id*  
Specifies the ID of the object to be removed.

**-d** *object\_id\_specification\_method* (__f__|__d__) (default = f)  
Specifies the method for specifying the ID of the object to be removed. Specifying __f__ indicates that the following object-ID specification is to be treated as a file name. That file shall consist of one entry per line, each indicting the ID of an object to be removed. Specifying __d__ indicates that the following object-ID specification is to be treated simply as an object-ID referring to the object to be removed.

### Examples of using ngt command

#### Create

Construct an index for 128-dimensional, 1-byte-integer data:

      $ cd (NGT_TOP_DIR)
      $ ngt create -d 128 -o c index ./data/sift-dataset-5k.tsv
      Data loading time=0.160748 (sec) 160.748 (msec)
      # of objects=5000
      Index creation time=0.379659 (sec) 379.659 (msec)

#### Search

Perform a neighborhood search by three queries specified in a file:

      $ cd (NGT_TOP_DIR)
      $ ngt search -n 20 index ./data/sift-query-3.tsv
      Query No.1
      Rank	ID		Distance
      1		3031	239.332
      2		4079	240.002
      3		3164	244.504
      4		3718	246.763
      5		157		251.094
      6		2422	251.185
      7		1313	251.34
      8		379		252.446
      9		3521	260.158
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
      20	679		272.973
      Query Time= 0.000472 (sec), 0.472 (msec)
      Query No.2
      Rank	ID		Distance
      1		2726	291.983
      2		924		296.987
      3		3638	298.585
      4		858		300.376
      5		1453	306.805
      6		174		307.789
      7		2992	308.485
      8		2980	308.93
      9		1525	309.49
      10	244		309.816
      11	910		310.446
      12	3310	310.585
      13	2433	311.482
      14	1633	311.735
      15	3761	312.44
      16	407		313.252
      17	4546	313.876
      18	697		315.108
      19	34		315.563
      20	2189	316.193
      Query Time= 0.000478 (sec), 0.478 (msec)
      Query No.3
      Rank	ID		Distance
      1		762		194.286
      2		1046	212.695
      3		4906	215.244
      4		2905	216.539
      5		4142	219.479
      6		1879	219.616
      7		4398	223.352
      8		3842	223.468
      9		233		224.127
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

License
-------

This software in part uses patents held by Yahoo Japan Corporation and therefore 
adopts the license shown below. For this reason, please refrain from commercial 
use.

This work is licensed under the Creative Commons  Attribution-NonCommercial-ShareAlike 
4.0 International License. To view a copy of this license, visit 
http://creativecommons.org/licenses/by-nc-sa/4.0/.

Publications
------------

-   Sugawara, K., Kobayashi, H. and Iwasaki, M.: On Approximately Searching for Similar Word Embeddings. Proc. of ACL2016 (2016) 2265-2275. ([pdf](https://aclweb.org/anthology/P/P16/P16-1214.pdf))
-   Iwasaki, M.: Applying a Graph-Structured Index to Product Image Search (in Japanese). IIEEJ Journal 42(5) (2013) 633-641. ([pdf](http://i.yimg.jp/i/docs/research_lab/articles/miwasaki-iieej-jnl-2013.pdf))
-   Iwasaki, M.: Proximity search using approximate k nearest neighbor graph with a tree structured index (in Japanese). IPSJ Journal 52(2) (2011) 817-828. ([pdf](http://i.yimg.jp/i/docs/research_lab/articles/miwasaki-ipsj-jnl-2011.pdf))
-   Iwasaki, M.: Proximity search in metric spaces using approximate k nearest neigh-bor graph (in Japanese). IPSJ Trans. on Database 3(1) (2010) 18-28. ([pdf](http://i.yimg.jp/i/docs/research_lab/articles/miwasaki-ipsj-tod-2010.pdf))

Copyright &copy; 2015-2016 Yahoo Japan Corporation All Rights Reserved.

