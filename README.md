NGT
===

Neighborhood Graph and Tree for Indexing High-dimensional Data

[Home](/README.md) / [Build](/README.md#build) / [Command](/bin/ngt/README.md#command) / [License](/README.md#license) / [Publications](/README.md#publications) / [About Us](http://research-lab.yahoo.co.jp/)

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

##### PANNG
- Sugawara, K., Kobayashi, H. and Iwasaki, M.: On Approximately Searching for Similar Word Embeddings. Proc. of ACL2016 (2016) 2265-2275. ([pdf](https://aclweb.org/anthology/P/P16/P16-1214.pdf))

##### ANNGT
- Iwasaki, M.: Applying a Graph-Structured Index to Product Image Search (in Japanese). IIEEJ Journal 42(5) (2013) 633-641. ([pdf](http://i.yimg.jp/i/docs/research_lab/articles/miwasaki-iieej-jnl-2013.pdf))
- Iwasaki, M.: Proximity search using approximate k nearest neighbor graph with a tree structured index (in Japanese). IPSJ Journal 52(2) (2011) 817-828. ([pdf](http://i.yimg.jp/i/docs/research_lab/articles/miwasaki-ipsj-jnl-2011.pdf))

##### ANNG
- Iwasaki, M.: Proximity search in metric spaces using approximate k nearest neigh-bor graph (in Japanese). IPSJ Trans. on Database 3(1) (2010) 18-28. ([pdf](http://i.yimg.jp/i/docs/research_lab/articles/miwasaki-ipsj-tod-2010.pdf))

Copyright &copy; 2015-2016 Yahoo Japan Corporation All Rights Reserved.

