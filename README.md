NGT
===

Neighborhood Graph and Tree for Indexing High-dimensional Data

[Home](/README.md) / [Installation](/README.md#Installation) / [Command](/bin/ngt/README.md#command) / [License](/README.md#license) / [Publications](/README.md#publications) / [About Us](http://research-lab.yahoo.co.jp/en/) / [日本語](/README-jp.md)

**NGT** provides commands and a library for performing high-speed approximate nearest neighbor searches against a large volume of data (several million to several 10 million items of data) in high dimensional vector data space (several ten to several thousand dimensions).

News
----

- 01/17/2019 Python NGT can be installed via pip from PyPI
- 12/14/2018 [NGTQ](bin/ngtq/README.md) (NGT with Quantization) is now available (v1.5.0)
- 08/08/2018 [ONNG](README.md#onng) is now available. (v1.4.0)

Installation
------------

### Downloads

- [Releases](https://github.com/yahoojapan/NGT/releases)

### Build

#### On Linux

      $ unzip NGT-x.x.x.zip
      $ cd NGT-x.x.x
      $ mkdir build
      $ cd build 
      $ cmake ..
      $ make 
      $ make install
      $ ldconfig /usr/local/lib

#### On Mac using homebrew

      $ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
      $ brew install cmake
      $ brew install gcc
      $ ln -s ./gcc-8 /usr/local/bin/gcc
      $ ln -s ./g++-8 /usr/local/bin/g++
      $ export CXX=g++
      $ export CC=gcc
      $ unzip NGT-x.x.x.zip
      $ cd NGT-x.x.x
      $ mkdir build
      $ cd build 
      $ cmake ..
      $ make 
      $ make install

#### Shared memory use

The index can be placed in shared memory with memory mapped files. Using shared memory can reduce the amount of memory needed when multiple processes are using the same index. In addition, it can not only handle an index with a large number of objects that cannot be loaded into memory, but also reduce time to open it. Since changes become necessary at build time, please add the following parameter when executing "cmake" in order to use shared memory.

      $ cmake -DNGT_SHARED_MEMORY_ALLOCATOR=ON ..

Note: Since there is no lock function, the index should be used only for reference when multiple processes are using the same index.

#### Large-scale data use

When you insert more than about 5 million objects, please add the following parameter to improve the search time.

      $ cmake -DNGT_LARGE_DATASET=ON ..

Utilities
---------

- Command : [ngt](/bin/ngt/README.md#command), [ngtq](bin/ngtq/README.md)
- Server : [ngtd](https://github.com/yahoojapan/ngtd)

Supported Programming Languages
-------------------------------

- [Python](/python/README.md)
- [Go](https://github.com/yahoojapan/gongt)
- C
- C++([sample code](bin/search/search.cpp))

License
-------

Copyright (C) 2015-2019 Yahoo Japan Corporation

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this software except in compliance with the License. You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and limitations under the License.

Yahoo Japan Corporation has aquired several patents on the technologies utilized in this software. However, the patent rights shall not be exercised under Apache License Version 2.0, only when the patented techniques are used with this software.

Contributor License Agreement
-----------------------------

This project requires contributors to accept the terms in the [Contributor License Agreement (CLA)](https://gist.github.com/ydnjp/3095832f100d5c3d2592).

Please note that contributors to the NGT repository on GitHub (https://github.com/yahoojapan/NGT) shall be deemed to have accepted the CLA without individual written agreements.

Publications
------------
##### [ONNG](bin/ngt/README.md#onng)
- Iwasaki, M., Miyazaki, D.: Optimization of Indexing Based on k-Nearest Neighbor Graph for Proximity. arXiv:1810.07355 [cs] (2018). ([pdf](https://arxiv.org/abs/1810.07355))

##### [PANNG](bin/ngt/README.md#panng)
- Iwasaki, M.: Pruned Bi-directed K-nearest Neighbor Graph for Proximity Search. Proc. of SISAP2016 (2016) 20-33. ([pdf](https://link.springer.com/chapter/10.1007/978-3-319-46759-7_2))
- Sugawara, K., Kobayashi, H. and Iwasaki, M.: On Approximately Searching for Similar Word Embeddings. Proc. of ACL2016 (2016) 2265-2275. ([pdf](https://aclweb.org/anthology/P/P16/P16-1214.pdf))

##### [ANNGT](bin/ngt/README.md#anngt)
- Iwasaki, M.: Applying a Graph-Structured Index to Product Image Search (in Japanese). IIEEJ Journal 42(5) (2013) 633-641. ([pdf](https://s.yimg.jp/i/docs/research_lab/articles/miwasaki-iieej-jnl-2013.pdf))
- Iwasaki, M.: Proximity search using approximate k nearest neighbor graph with a tree structured index (in Japanese). IPSJ Journal 52(2) (2011) 817-828. ([pdf](https://s.yimg.jp/i/docs/research_lab/articles/miwasaki-ipsj-jnl-2011.pdf))

##### [ANNG](bin/ngt/README.md#anng)
- Iwasaki, M.: Proximity search in metric spaces using approximate k nearest neighbor graph (in Japanese). IPSJ Trans. on Database 3(1) (2010) 18-28. ([pdf](https://s.yimg.jp/i/docs/research_lab/articles/miwasaki-ipsj-tod-2010.pdf))

Copyright &copy; 2015-2019 Yahoo Japan Corporation All Rights Reserved.

