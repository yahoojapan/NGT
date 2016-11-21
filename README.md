NGT
===

Neighborhood Graph and Tree for Indexing High-dimensional Data

[Home](/README.md) / [Build](/README.md#build) / [Command](/bin/ngt/README.md#command) / [License](/README.md#license) / [Publications](/README.md#publications) / [About Us](http://research-lab.yahoo.co.jp/)

**NGT** provides commands and a library for performing high-speed approximate nearest neighbor searches against a large volume of data (several million to several 10 million items of data) in high dimensional vector data space (several ten to several thousand dimensions).

#### Downloads

- [Source code](https://github.com/yahoojapan/NGT/)
- [Experimental software](https://s.yimg.jp/dl/docs/research_lab/ngt-experimental-1.0.zip) : This software was used in our [ACL 2016](http://acl2016.org/) [paper](https://aclweb.org/anthology/P/P16/P16-1214.pdf).

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

License
-------

Copyright (C) 2015-2016 Yahoo Japan Corporation

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this software except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Publications
------------

##### PANNG
- Iwasaki, M.: Pruned Bi-directed K-nearest Neighbor Graph for Proximity Search. Proc. of SISAP2016 (2016) 20-33.
- Sugawara, K., Kobayashi, H. and Iwasaki, M.: On Approximately Searching for Similar Word Embeddings. Proc. of ACL2016 (2016) 2265-2275. ([pdf](https://aclweb.org/anthology/P/P16/P16-1214.pdf))

##### ANNGT
- Iwasaki, M.: Applying a Graph-Structured Index to Product Image Search (in Japanese). IIEEJ Journal 42(5) (2013) 633-641. ([pdf](http://i.yimg.jp/i/docs/research_lab/articles/miwasaki-iieej-jnl-2013.pdf))
- Iwasaki, M.: Proximity search using approximate k nearest neighbor graph with a tree structured index (in Japanese). IPSJ Journal 52(2) (2011) 817-828. ([pdf](http://i.yimg.jp/i/docs/research_lab/articles/miwasaki-ipsj-jnl-2011.pdf))

##### ANNG
- Iwasaki, M.: Proximity search in metric spaces using approximate k nearest neigh-bor graph (in Japanese). IPSJ Trans. on Database 3(1) (2010) 18-28. ([pdf](http://i.yimg.jp/i/docs/research_lab/articles/miwasaki-ipsj-tod-2010.pdf))

Copyright &copy; 2015-2016 Yahoo Japan Corporation All Rights Reserved.

