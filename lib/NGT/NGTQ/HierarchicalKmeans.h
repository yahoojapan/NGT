//
// Copyright (C) 2021 Yahoo Japan Corporation
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

#include "Quantizer.h"

namespace QBG {
  class HierarchicalKmeans {
  public:
    typedef NGTQ::Quantizer::ObjectList QBGObjectList;

    class HKNode {
    public:
      bool leaf;
    };

    class HKLeafNode : public HKNode {
    public:
    HKLeafNode():id(0){ leaf = true; }
      std::vector<uint32_t> members;
      uint32_t id;
    };

    class HKInternalNode : public HKNode {
    public:
      HKInternalNode() { leaf = false; }
      std::vector<std::pair<uint32_t, std::vector<float>>> children;
    };

    HierarchicalKmeans() {
      silence = true;
    }

    HierarchicalKmeans(QBG::BuildParameters &param) {
#ifdef NGTQ_QBG
      maxSize			= param.hierarchicalClustering.maxSize;
      numOfObjects		= param.hierarchicalClustering.numOfObjects;
      numOfClusters		= param.hierarchicalClustering.numOfClusters;
      numOfTotalClusters	= param.hierarchicalClustering.numOfTotalClusters;
      numOfTotalBlobs		= param.hierarchicalClustering.numOfTotalBlobs;
      clusterID			= param.hierarchicalClustering.clusterID;

      initMode			= param.hierarchicalClustering.initMode;

      numOfRandomObjects	= param.hierarchicalClustering.numOfRandomObjects;

      numOfFirstObjects		= param.hierarchicalClustering.numOfFirstObjects;
      numOfFirstClusters	= param.hierarchicalClustering.numOfFirstClusters;
      numOfSecondObjects	= param.hierarchicalClustering.numOfSecondObjects;
      numOfSecondClusters	= param.hierarchicalClustering.numOfSecondClusters;
      numOfThirdClusters	= param.hierarchicalClustering.numOfThirdClusters;
      extractCentroid		= param.hierarchicalClustering.extractCentroid;

      threeLayerClustering	= param.hierarchicalClustering.threeLayerClustering;
      silence			= param.silence;
#endif
    }

    static int32_t searchLeaf(std::vector<HKNode*> &nodes, int32_t rootID, float *object) {
      auto nodeID = rootID;
      while (true) {
	auto *node = nodes[nodeID];
	if (node->leaf) {
	  return nodeID;
	} else {
	  HKInternalNode &internalNode = static_cast<HKInternalNode&>(*node);      
	  float min = std::numeric_limits<float>::max();
	  int32_t minid = 0;
	  for (auto &c : internalNode.children) {
	    auto d = NGT::PrimitiveComparator::compareL2(reinterpret_cast<float*>(&object[0]),
							 c.second.data(), c.second.size());
	    if (d < min) {
	      min = d;
	      minid = c.first;
	    }
	  }
	  nodeID = minid;
	}
      }
      return -1;
    }

    static void aggregateObjects(HKLeafNode &leafNode, std::vector<std::vector<float>> &vectors,
			  NGT::ObjectSpace &objectSpace, QBGObjectList &objectList)
    {
      vectors.reserve(leafNode.members.size() + 1);
      std::vector<float> obj;
      for (auto &m : leafNode.members) {
	objectList.get(m, obj, &objectSpace);
	vectors.push_back(obj);
      }
    }

    static void aggregateObjects(HKLeafNode &leafNode, std::vector<std::vector<float>> &vectors,
			  NGT::ObjectSpace &objectSpace, QBGObjectList &objectList,
			  std::vector<float> &object)
    {
      aggregateObjects(leafNode, vectors, objectSpace, objectList);
      vectors.push_back(std::move(object));
    }

    static void split(uint32_t id, std::vector<std::vector<float>> &vectors,
				 std::vector<HKNode*> &nodes, int32_t leafNodeID, NGT::Clustering &clustering)
    {
      HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[leafNodeID]);
      std::vector<NGT::Clustering::Cluster> clusters;
      clustering.kmeans(vectors, clusters);
      auto *newNode = new HKInternalNode;
      for (auto &cluster : clusters) {
	auto centroid = std::move(cluster.centroid);
	centroid.resize(vectors[0].size());
	newNode->children.push_back(std::make_pair(nodes.size(), std::move(centroid)));
	auto *cnode = new HKLeafNode;
	nodes.push_back(cnode);
	for (auto &member : cluster.members) {
	  if (member.vectorID > leafNode.members.size()) {
	    std::cerr << "Fatal error. member:" << member.vectorID << ":" << leafNode.members.size() << std::endl;
	    abort();
	  }
	  if (member.vectorID == leafNode.members.size()) {
	    cnode->members.push_back(id);
	  } else {
	    cnode->members.push_back(leafNode.members[member.vectorID]);
	  }
	}
      }
      delete nodes[leafNodeID]; 
      nodes[leafNodeID] = newNode; 
    }

    static double computeError(std::vector<HKNode*> &nodes, NGT::ObjectSpace &objectSpace, QBGObjectList &objectList) {
      std::cerr << "node size=" << nodes.size() << std::endl;
      double distance = 0.0;
      size_t dcount = 0;
      for (auto *node : nodes) {
	if (node->leaf) {
	} else {
	  HKInternalNode &internalNode = static_cast<HKInternalNode&>(*node);
	  std::vector<float> obj;
	  for (auto &child : internalNode.children) {
	    if (nodes[child.first]->leaf) {
	      if (dcount % 100000 == 0) {
		std::cerr << "Processed leaves=" << dcount << std::endl;
	      }
	      auto centroid = child.second;
	      HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[child.first]);
	      for (auto &m : leafNode.members) {
		objectList.get(m, obj, &objectSpace);
		distance += NGT::Clustering::distanceL2(centroid, obj);
		dcount++;
	      }
	    }
	  }
	} 
      }
      distance /= dcount;
      std::cout << "# of vectors=" << dcount << std::endl;
      std::cout << "Quantization error=" << distance << std::endl;
      return distance;
    }

    static size_t extractCentroids(std::ostream &oStream, std::vector<HKNode*> &nodes) {
      std::cerr << "node size=" << nodes.size() << std::endl;
      size_t clusterCount = 0;
      size_t objectCount = 0;
      size_t leafID = 0;
      for (auto *node : nodes) {
	if (node->leaf) {
	  HKLeafNode &leafNode = static_cast<HKLeafNode&>(*node);
	  objectCount += leafNode.members.size();
	} else {
	  HKInternalNode &internalNode = static_cast<HKInternalNode&>(*node);
	  for (auto &child : internalNode.children) {
	    if (nodes[child.first]->leaf) {
	      if (static_cast<HKLeafNode*>(nodes[child.first])->id == 0) {
		static_cast<HKLeafNode*>(nodes[child.first])->id = leafID;
	      } else if (static_cast<HKLeafNode*>(nodes[child.first])->id != leafID) {
		std::cerr << "leaf ID is invalid?" << std::endl;
	      }
	      leafID++;
	      size_t count = 0;
	      clusterCount++;
	      for (auto &v : child.second) {
		oStream << v;
		if (++count == child.second.size()) {
		  oStream << std::endl;
		} else {
		  oStream << "\t";
		}
	      }
	    }
	  }
	} 
      }
      std::cerr << "# of clusters=" << clusterCount << std::endl;
      return objectCount;
    }

    static size_t extractIndex(std::ostream &oStream, std::vector<HKNode*> &nodes, size_t numOfObjects) {
      std::vector<int32_t> clusterID(numOfObjects, -1);
      std::cerr << "numOfObjects=" << numOfObjects << std::endl;
      std::cerr << "node size=" << nodes.size() << std::endl;
      for (auto *node : nodes) {
	if (node->leaf) {
	  HKLeafNode &leafNode = static_cast<HKLeafNode&>(*node);
	  for (auto &member : leafNode.members) {
	    if (member > numOfObjects) {
	      std::cerr << "output index: Internal fatal error. " << member << ":" << numOfObjects - 1 << std::endl;
	      abort();
	    }
	    if (member == 0) {
	      std::cerr << "output index: Internal fatal error. Invalid ID" << std::endl;
	      abort();
	    }
	    clusterID[member - 1] = leafNode.id;
	  }
	}
      }
      std::cerr << "clusterID.size=" << clusterID.size() << std::endl;
      size_t count = 0;
      for (auto cid : clusterID) {
	count++;
	oStream << cid << std::endl;
      }
      std::cerr << "# of id=" << count << std::endl;
      return count;
    }

    static void extractBtoQAndQCentroid(std::ostream &btoqStream, std::ostream &qStream,
				 std::vector<HKNode*> &nodes, size_t numOfThirdClusters) {
      std::cerr << "extractBtoQ" << std::endl;
      std::vector<int32_t> btoq(numOfThirdClusters);
      std::cerr << "numOfThirdClusters=" << numOfThirdClusters << std::endl;
      std::cerr << "node size=" << nodes.size() << std::endl;
      size_t rootID = 0;
      HKInternalNode &root = static_cast<HKInternalNode&>(*nodes[rootID]);
      std::cerr << "first=" << root.children.size() << std::endl;
      size_t secondCount = 0;	
      size_t thirdCount = 0;	
      size_t objectCount = 0;
      size_t leafID = 0;
      size_t qID = 0;
      for (auto &c1 : root.children) {
	HKInternalNode &node1 = static_cast<HKInternalNode&>(*nodes[c1.first]);
	std::cerr << "second=" << node1.children.size() << std::endl;
	secondCount += node1.children.size();
	for (auto &c2 : node1.children) {
	  HKInternalNode &node2 = static_cast<HKInternalNode&>(*nodes[c2.first]);
	  std::cerr << "third=" << node2.children.size() << std::endl;
	  thirdCount += node2.children.size();
	  size_t count = 0;
	  for (auto &v : c2.second) {
	    qStream << v;
	    if (++count == c2.second.size()) {
	      qStream << std::endl;
	    } else {
	      qStream << "\t";
	    }
	  }
	  for (auto &c3 : node2.children) {
	    btoqStream << qID << std::endl;
	    HKLeafNode &leaf = static_cast<HKLeafNode&>(*nodes[c3.first]);
	    objectCount += leaf.members.size();
	    if (leaf.id != leafID++) {
	      std::cerr << "leaf is invalid" << leaf.id << ":" << leafID << std::endl;
	      abort();
	    }
	  }
	  qID++;
	}
      }
      std::cerr << "second=" << secondCount << std::endl;
      std::cerr << "third=" << thirdCount << std::endl;
      std::cerr << "object=" << objectCount << std::endl;
    }

    static void extractRandomObjectsFromEachBlob(std::ostream &oStream, std::vector<HKNode*> &nodes, size_t numOfObjects,
						 size_t numOfRandomObjects, NGTQ::QuantizerInstance<uint8_t>& quantizer, bool extractCentroid) {
      std::cerr << "node size=" << nodes.size() << std::endl;
      std::vector<std::vector<std::vector<float>>> randomObjects(numOfObjects);
      std::vector<std::vector<float>> centroids(numOfObjects);
      for (auto *node : nodes) {
	if (node->leaf) {
	  HKLeafNode &leafNode = static_cast<HKLeafNode&>(*node);
	  std::vector<uint32_t> randomObjectIDXs;
	  if (numOfRandomObjects >= leafNode.members.size()) {
	    randomObjectIDXs = leafNode.members;
	    while (randomObjectIDXs.size() < numOfRandomObjects) {
	      double random = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
	      uint32_t idx = floor(leafNode.members.size() * random);
	      if (idx >= leafNode.members.size()) {
		std::cerr << "Internal error. " << idx << ":" << leafNode.members.size() << std::endl;
		abort();
	      }
	      randomObjectIDXs.push_back(leafNode.members[idx]);	  
	    }
	  } else {
	    srand(leafNode.id);
	    while (randomObjectIDXs.size() < numOfRandomObjects) {
	      uint32_t idx = 0;
	      do {
		double random = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
		idx = floor(leafNode.members.size() * random);
		if (idx >= leafNode.members.size()) {
		  std::cerr << "Internal error. " << idx << ":" << leafNode.members.size() << std::endl;
		  abort();
		}
	      } while (std::find(randomObjectIDXs.begin(), randomObjectIDXs.end(), leafNode.members[idx]) != randomObjectIDXs.end());
	      std::cerr << "IDX=" << idx << "/" << leafNode.members.size() << std::endl;
	      randomObjectIDXs.push_back(leafNode.members[idx]);
	    }
	  }
	  std::cerr << "randomObjectIDXs=" << randomObjectIDXs.size() << std::endl;
	  for (auto member : randomObjectIDXs) {
	    if (member == 0) {
	      std::cerr << "output index: Internal fatal error. Invalid ID. " <<  member << std::endl;
	      abort();
	    }
	    std::vector<float> object;
	    quantizer.objectList.get(member, object, &quantizer.globalCodebookIndex.getObjectSpace());
	    if (leafNode.id >= numOfObjects) {
	      std::cerr << "Internal error! Wrong leaf ID. " << leafNode.id << ":" << numOfObjects << std::endl;
	      abort();
	    }
	    randomObjects[leafNode.id].push_back(object);
	  }
	} else { 
	  if (extractCentroid) {
	    HKInternalNode &internalNode = static_cast<HKInternalNode&>(*node);
	    for (auto &child : internalNode.children) {
	      if (nodes[child.first]->leaf) {
		HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[child.first]);
		centroids[leafNode.id] = child.second;
	      }
	    }
	  }
	} 
      }
      for (size_t idx = 0; idx < centroids.size(); idx++) {
	auto &c = centroids[idx];
	if (extractCentroid && c.empty()) {
	  std::cerr << "qbg: Fatal error! The centroid is empty." << std::endl;
	  abort();
	}
	for (size_t i = 0; i < c.size(); i++) {
	  oStream << c[i];
	  if (i + 1 != c.size()) {
	    oStream << "\t";
	  } else {
	    oStream << std::endl;;
	  }
	}
	auto &ros = randomObjects[idx];
	for (auto &ro : ros) {
	  if (ro.empty()) {
	    std::cerr << "qbg: Fatal error! The random object vector is empty." << std::endl;
	    abort();
	  }
	  for (size_t i = 0; i < ro.size(); i++) {
	    oStream << ro[i];
	    if (i + 1 != ro.size()) {
	      oStream << "\t";
	    } else {
	      oStream << std::endl;;
	    }
	  }
	}
      }
    }

    static void extractBtoQIndex(std::ofstream &of, std::vector<HKNode*> &nodes, std::vector<uint32_t> &qNodeIDs) {
      size_t leafID = 0;
      for (size_t qnidx = 0; qnidx < qNodeIDs.size(); qnidx++) {
	if (nodes[qNodeIDs[qnidx]]->leaf) {
	  std::cerr << "Fatal error. this should be an internal node." << std::endl;
	  abort();
	}
	HKInternalNode &inode = static_cast<HKInternalNode&>(*nodes[qNodeIDs[qnidx]]);
	for (auto &c : inode.children) {
	  if (!nodes[c.first]->leaf) {
	    std::cerr << "Fatal error. this should be a leaf." << std::endl;
	    abort();
	  }
	  HKLeafNode &leaf = static_cast<HKLeafNode&>(*nodes[c.first]);
	  if (leaf.id == 0) {
	    leaf.id = leafID;
	  }
	  of << qnidx << std::endl;
	  leafID++;
	}
      }
}


    static void hierarchicalKmeans(uint32_t id, int32_t rootID, std::vector<float> &object,
				   QBGObjectList &objectList, NGT::ObjectSpace &objectSpace,
				   std::vector<HKNode*> &nodes,   NGT::Clustering &clustering, size_t maxSize) {
      NGT::Timer timer;
      objectList.get(id, object, &objectSpace);
      int32_t nodeID = searchLeaf(nodes, rootID, reinterpret_cast<float*>(&object[0]));
      if (nodeID < 0) {
	std::cerr << "Fatal inner error! node ID=" << nodeID << std::endl;
	exit(1);
      }
      auto *node = nodes[nodeID];
      HKLeafNode &leafNode = static_cast<HKLeafNode&>(*node);
      if (leafNode.members.size() >= maxSize) {
	NGT::Timer subtimer;
	subtimer.start();
	std::vector<std::vector<float>> vectors;
	aggregateObjects(leafNode, vectors, objectSpace, objectList, object);
	subtimer.stop();
	std::cerr << "aggregate time=" << subtimer << std::endl;
	subtimer.start();
	split(id, vectors, nodes, nodeID, clustering);
	subtimer.stop();
	std::cerr << "split time=" << subtimer << std::endl;
      } else {
	leafNode.members.push_back(id);
      }
    }

    static void hierarchicalKmeansBatch(std::vector<uint32_t> &batch, std::vector<pair<uint32_t, uint32_t>> &exceededLeaves,
					int32_t rootID, std::vector<float> &object,
					QBGObjectList &objectList, NGT::ObjectSpace &objectSpace,
					std::vector<HKNode*> &nodes, NGT::Clustering &clustering, size_t maxSize, size_t &nleaves,
					size_t maxExceededLeaves) {

      if (batch.size() == 0) {
	return;
      }

      int32_t nodeIDs[batch.size()];

#pragma omp parallel for
      for (size_t idx = 0; idx < batch.size(); idx++) {
	auto id = batch[idx];
#pragma omp critical
	objectList.get(id, object, &objectSpace);
	int32_t nodeID = searchLeaf(nodes, rootID, reinterpret_cast<float*>(&object[0]));
	if (nodeID < 0) {
	  std::cerr << "Fatal inner error! node ID=" << nodeID << std::endl;
	  exit(1);
        }
	nodeIDs[idx] = nodeID;
      }

  
      for (size_t idx = 0; idx < batch.size(); idx++) {
	auto id = batch[idx];
	HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[nodeIDs[idx]]);
	leafNode.members.push_back(id);
	if (leafNode.members.size() > maxSize) {
	  auto i = exceededLeaves.begin();
	  for (; i != exceededLeaves.end(); i++) {
	    if (static_cast<int32_t>((*i).second) == nodeIDs[idx]) break;
	  }
	  if (i == exceededLeaves.end()) {
	    exceededLeaves.push_back(std::make_pair(batch[idx], nodeIDs[idx]));
	  }
	}
      }

      batch.clear();

      if (exceededLeaves.size() < maxExceededLeaves) {
	return;
      }
	
      std::vector<std::vector<NGT::Clustering::Cluster>> clusters(exceededLeaves.size());
#pragma omp parallel for
      for (size_t idx = 0; idx < exceededLeaves.size(); idx++) {
	HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[exceededLeaves[idx].second]);
	std::vector<std::vector<float>> vectors;
#pragma omp critical
	aggregateObjects(leafNode, vectors, objectSpace, objectList);
	clustering.kmeans(vectors, clusters[idx]);
      }

      std::cerr << "exceeded leaves=" << exceededLeaves.size() << std::endl;
      for (size_t idx = 0; idx < exceededLeaves.size(); idx++) {
	auto leafNodeID = exceededLeaves[idx].second;
	HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[leafNodeID]);
	auto *newNode = new HKInternalNode;
	for (auto &cluster : clusters[idx]) {
	  newNode->children.push_back(std::make_pair(nodes.size(), std::move(cluster.centroid)));
	  auto *cnode = new HKLeafNode;
	  nodes.push_back(cnode);
	  for (auto &member : cluster.members) {
	    cnode->members.push_back(leafNode.members[member.vectorID]);
	  }
	}
	nleaves += clusters[idx].size() - 1;
	delete nodes[leafNodeID]; 
	nodes[leafNodeID] = newNode; 
      }
      exceededLeaves.clear();

    }

    static void hierarchicalKmeansWithNumberOfClusters(size_t numOfTotalClusters, size_t numOfObjects, size_t numOfLeaves, 
						       QBGObjectList &objectList, NGT::ObjectSpace &objectSpace, 
						       std::vector<HKNode*> &nodes, NGT::Clustering::InitializationMode initMode){
      std::cerr << "numOfTotalClusters=" << numOfTotalClusters << std::endl;
      std::cerr << "numOfLeaves=" << numOfLeaves << std::endl;
      if (numOfLeaves > numOfTotalClusters) {
	std::cerr << "# of clusters is invalid. " << numOfLeaves << ":" << numOfTotalClusters << std::endl;
	abort();
      }
      auto numOfRemainingClusters = numOfTotalClusters;
      auto numOfRemainingVectors = numOfObjects;
      size_t leafCount = 0;
      size_t nodeSize = nodes.size(); 
      for (size_t nidx = 0; nidx < nodeSize; nidx++) {
	if (nodes[nidx]->leaf) {
	  leafCount++;
	  if (numOfLeaves >= 100 && leafCount % (numOfLeaves / 100) == 0) {
	    std::cerr << "Processed leaves: " << leafCount << " " << leafCount * 100 / numOfLeaves << "%" << std::endl;
	  }
	  HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[nidx]);
	  std::vector<std::vector<float>> vectors;
	  aggregateObjects(leafNode, vectors, objectSpace, objectList);
	  size_t nClusters = round(static_cast<float>(leafNode.members.size()) / numOfRemainingVectors * numOfRemainingClusters);
	  nClusters = nClusters == 0 ? 1 : nClusters;
	  numOfRemainingVectors -= leafNode.members.size();
	  numOfRemainingClusters -= nClusters;
	  NGT::Clustering clustering(initMode, NGT::Clustering::ClusteringTypeKmeansWithoutNGT, 1000, nClusters);
	  NGT::Timer timer;
	  timer.start();
	  split(0, vectors, nodes, nidx, clustering);
	  timer.stop();
	  if (nodes[nidx]->leaf) {
	    std::cerr << "At this moment, the second node should be an internal" << std::endl;
	    abort();
	  }
	}
      } 
    }

    static void hierarchicalKmeansWithNumberOfClustersInParallel(size_t numOfTotalClusters, size_t numOfObjects, size_t numOfLeaves, 
								 QBGObjectList &objectList, NGT::ObjectSpace &objectSpace, 
								 std::vector<HKNode*> &nodes, NGT::Clustering::InitializationMode initMode){
      NGT::Timer timer;
      timer.start();
      auto numOfRemainingClusters = numOfTotalClusters;
      auto numOfRemainingVectors = numOfObjects;
      size_t leafCount = 0;

      std::vector<pair<uint32_t, size_t>> leafNodes;
      leafNodes.reserve(numOfLeaves);
      for (size_t nidx = 0; nidx < nodes.size(); nidx++) {
	if (nodes[nidx]->leaf) {
	  leafCount++;
	  {
	    size_t step = 10;
	    if (numOfLeaves >= step && leafCount % (numOfLeaves / step) == 0) {
	      std::cerr << "Processed leaves: " << leafCount << " " << leafCount * step / numOfLeaves << "%" << std::endl;
	    }
	  }
	  HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[nidx]);
	  size_t nClusters = round(static_cast<float>(leafNode.members.size()) / numOfRemainingVectors * numOfRemainingClusters);
	  nClusters = nClusters == 0 ? 1 : nClusters;
	  numOfRemainingVectors -= leafNode.members.size();
	  numOfRemainingClusters -= nClusters;
	  leafNodes.push_back(std::make_pair(nidx, nClusters));
	}
      }
      timer.stop();
      std::cerr << "hierarchicalKmeansWithNumberOfClustersInParallel: extract leaves. Time=" << timer << std::endl;
      timer.start();

      std::cerr << "start kmeans..." << std::endl;
      std::vector<std::vector<NGT::Clustering::Cluster>> clusters(leafNodes.size());
#pragma omp parallel for
      for (size_t nidx = 0; nidx < leafNodes.size(); nidx++) {
	HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[leafNodes[nidx].first]);
	std::vector<std::vector<float>> vectors;
#pragma omp critical
	aggregateObjects(leafNode, vectors, objectSpace, objectList);
	NGT::Clustering clustering(initMode, NGT::Clustering::ClusteringTypeKmeansWithoutNGT, 1000, leafNodes[nidx].second);
	clustering.kmeans(vectors, clusters[nidx]);
      }

      timer.stop();
      std::cerr << "hierarchicalKmeansWithNumberOfClustersInParallel: kmeans. Time=" << timer << std::endl;
      timer.start();

      std::cerr << "add nodes..." << std::endl;
      for (size_t idx = 0; idx < leafNodes.size(); idx++) {
	auto leafNodeID = leafNodes[idx].first;
	HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[leafNodeID]);
	auto *newNode = new HKInternalNode;
	for (auto &cluster : clusters[idx]) {
	  newNode->children.push_back(std::make_pair(nodes.size(), std::move(cluster.centroid)));
	  auto *cnode = new HKLeafNode;
	  nodes.push_back(cnode);
	  for (auto &member : cluster.members) {
	    cnode->members.push_back(leafNode.members[member.vectorID]);
	  }
	}
	delete nodes[leafNodeID]; 
	nodes[leafNodeID] = newNode; 
      }
      timer.stop();
      std::cerr << "hierarchicalKmeansWithNumberOfClustersInParallel: add nodes. Time=" << timer << std::endl;

    }

    static void flattenClusters(std::vector<NGT::Clustering::Cluster> &upperClusters, 
				std::vector<std::vector<NGT::Clustering::Cluster>> &lowerClusters,
				size_t numOfLowerClusters,
				std::vector<NGT::Clustering::Cluster> &flatClusters) {


      flatClusters.clear();
      flatClusters.reserve(numOfLowerClusters);

      for (size_t idx1 = 0; idx1 < lowerClusters.size(); idx1++) {
	for (size_t idx2 = 0; idx2 < lowerClusters[idx1].size(); idx2++) {
	  for (auto &m : lowerClusters[idx1][idx2].members) {
	    m.vectorID = upperClusters[idx1].members[m.vectorID].vectorID;
	  }
	  flatClusters.push_back(lowerClusters[idx1][idx2]);
	}
      }

    }

#ifndef MULTIPLE_OBJECT_LISTS
    void subclustering(std::vector<NGT::Clustering::Cluster> &upperClusters, size_t numOfLowerClusters, size_t numOfObjects,
		       NGT::ObjectSpace &objectSpace, QBGObjectList &objectList, 
		       NGT::Clustering::InitializationMode initMode, std::vector<std::vector<NGT::Clustering::Cluster>> &lowerClusters) {
      std::vector<uint32_t> nPartialClusters(upperClusters.size());
      auto numOfRemainingClusters = numOfLowerClusters;
      auto numOfRemainingVectors = numOfObjects;
      size_t ts = 0;
      for (size_t idx = 0; idx < upperClusters.size(); idx++) {
	size_t ncs = round(static_cast<float>(upperClusters[idx].members.size()) / numOfRemainingVectors * 
			   numOfRemainingClusters);
	ncs = ncs == 0 ? 1 : ncs;
	numOfRemainingVectors -= upperClusters[idx].members.size();
	if (numOfRemainingClusters >= ncs) {
	  numOfRemainingClusters -= ncs;
	}
	nPartialClusters[idx] = ncs;
	ts += ncs;
      }
      std::cerr << "numOfRemainingClusters=" << numOfRemainingClusters << std::endl;
      std::cerr << "numOfRemainingVectors=" << numOfRemainingVectors << std::endl;
      std::cerr << "upperClusters=" << upperClusters.size() << std::endl;
      std::cerr << "total=" << ts << ":" << numOfLowerClusters << std::endl;
      if (ts < numOfLowerClusters || numOfRemainingClusters != 0) {
	std::cerr << "subclustering: Internal error! " << std::endl;
	exit(1);
      }

      lowerClusters.resize(upperClusters.size());
#pragma omp parallel for schedule(dynamic)  
      for (size_t idx = 0; idx < upperClusters.size(); idx++) {
	std::vector<std::vector<float>> partialVectors;
	partialVectors.reserve(upperClusters[idx].members.size());
	std::vector<float> obj;
#pragma omp critical
	{
	  for (auto &m : upperClusters[idx].members) {
	    objectList.get(m.vectorID + 1, obj, &objectSpace);
	    partialVectors.push_back(obj);
	  }
	}
	if (upperClusters[idx].members.size() != partialVectors.size()) {
	  std::cerr << "the sizes of members are not consistent" << std::endl;
	  abort();
	}
	NGT::Clustering lowerClustering(initMode, NGT::Clustering::ClusteringTypeKmeansWithoutNGT, 1000);
	lowerClustering.kmeans(partialVectors, nPartialClusters[idx], lowerClusters[idx]);
	if (nPartialClusters[idx] != lowerClusters[idx].size()) {
	  std::cerr << "the sizes of cluster members are not consistent" << std::endl;	    
	  abort();
	}
      }
      size_t nc = 0;
      size_t mc = 0;
      for (auto &cs : lowerClusters) {
	nc += cs.size();
	for (auto &c : cs) {
	  mc += c.members.size();
	}
      }
      std::cerr << "# of clusters=" << nc << " # of members=" << mc << std::endl;
    }
    void subclustering(std::vector<NGT::Clustering::Cluster> &upperClusters, size_t numOfLowerClusters, size_t numOfObjects,
		       NGT::ObjectSpace &objectSpace, QBGObjectList &objectList, 
		       NGT::Clustering::InitializationMode initMode, std::vector<NGT::Clustering::Cluster> &flatLowerClusters) {

      std::vector<std::vector<NGT::Clustering::Cluster>> lowerClusters;
      subclustering(upperClusters, numOfLowerClusters, numOfObjects, objectSpace, objectList, initMode, lowerClusters);

      flattenClusters(upperClusters, lowerClusters, numOfLowerClusters, flatLowerClusters);

    }

#else 
    static void subclustering(std::vector<NGT::Clustering::Cluster> &upperClusters, size_t numOfLowerClusters, size_t numOfObjects,
			      NGT::ObjectSpace &objectSpace, QBGObjectList &objectList, 
			      NGT::Clustering::InitializationMode initMode, std::vector<std::vector<NGT::Clustering::Cluster>> &lowerClusters) {
      std::vector<uint32_t> nPartialClusters(upperClusters.size());
      auto numOfRemainingClusters = numOfLowerClusters;
      auto numOfRemainingVectors = numOfObjects;
      size_t ts = 0;
      for (size_t idx = 0; idx < upperClusters.size(); idx++) {
	size_t ncs = round(static_cast<float>(upperClusters[idx].members.size()) / numOfRemainingVectors * 
			   numOfRemainingClusters);
	ncs = ncs == 0 ? 1 : ncs;
	numOfRemainingVectors -= upperClusters[idx].members.size();
	if (numOfRemainingClusters >= ncs) {
	  numOfRemainingClusters -= ncs;
	}
	nPartialClusters[idx] = ncs;
	ts += ncs;
      }

      std::cerr << "numOfRemainingClusters=" << numOfRemainingClusters << std::endl;
      std::cerr << "numOfRemainingVectors=" << numOfRemainingVectors << std::endl;
      std::cerr << "upperClusters=" << upperClusters.size() << std::endl;
      std::cerr << "total=" << ts << ":" << numOfLowerClusters << std::endl;
      if (ts < numOfLowerClusters || numOfRemainingClusters != 0) {
	std::cerr << "subclustering: Internal error! " << std::endl;
	exit(1);
      }

      auto nthreads = omp_get_max_threads();
      if (!objectList.openMultipleStreams(nthreads)) {
	std::cerr << "Cannot open multiple streams." << std::endl;
	abort();
      }

      lowerClusters.resize(upperClusters.size());
#pragma omp parallel for schedule(dynamic)
      for (size_t idx = 0; idx < upperClusters.size(); idx++) {
	std::vector<std::vector<float>> partialVectors;
	partialVectors.reserve(upperClusters[idx].members.size());
	std::vector<float> obj;
	auto threadid = omp_get_thread_num();
	//#pragma omp critical
	{
	  for (auto &m : upperClusters[idx].members) {
	    if (threadid >= nthreads) {
	      std::cerr << "inner fatal error. # of threads=" << nthreads << ":" << threadid << std::endl;
	      exit(1);
	    }
	    if (!objectList.get(threadid, m.vectorID + 1, obj, &objectSpace)) {
	      std::cerr << "subclustering: Fatal error! cannot get!!!! " << m.vectorID + 1 << std::endl;
	      abort();
	    }
	    partialVectors.push_back(obj);
	  }
	}
	if (upperClusters[idx].members.size() != partialVectors.size()) {
	  std::cerr << "the sizes of members are not consistent" << std::endl;
	  abort();
	}
	NGT::Clustering lowerClustering(initMode, NGT::Clustering::ClusteringTypeKmeansWithoutNGT, 1000);
	lowerClustering.kmeans(partialVectors, nPartialClusters[idx], lowerClusters[idx]);
	if (nPartialClusters[idx] != lowerClusters[idx].size()) {
	  std::cerr << "the sizes of cluster members are not consistent" << std::endl;	    
	  abort();
	}
      }
      size_t nc = 0;
      size_t mc = 0;
      for (auto &cs : lowerClusters) {
	nc += cs.size();
	for (auto &c : cs) {
	  mc += c.members.size();
	}
      }
      std::cerr << "# of clusters=" << nc << " # of members=" << mc << std::endl;
    }

    static void subclustering(std::vector<NGT::Clustering::Cluster> &upperClusters, size_t numOfLowerClusters, size_t numOfObjects,
			      NGT::ObjectSpace &objectSpace, QBGObjectList &objectList,
			      NGT::Clustering::InitializationMode initMode, std::vector<NGT::Clustering::Cluster> &flatLowerClusters) {

      std::vector<std::vector<NGT::Clustering::Cluster>> lowerClusters;
      subclustering(upperClusters, numOfLowerClusters, numOfObjects, objectSpace, objectList, initMode, lowerClusters);

      flattenClusters(upperClusters, lowerClusters, numOfLowerClusters, flatLowerClusters);

    }

#endif 


    static void assign(std::vector<NGT::Clustering::Cluster> &clusters, size_t beginID, size_t endID,
		       NGT::ObjectSpace &objectSpace, QBGObjectList &objectList) {

#ifdef MULTIPLE_OBJECT_LISTS
      if (!objectList.openMultipleStreams(omp_get_max_threads())) {
	std::cerr << "Cannot open multiple streams." << std::endl;
	abort();
      }
#endif

      size_t count = 0;
#pragma omp parallel for
      for (size_t id = beginID; id <= endID; id++) {
	std::vector<float> obj;
	//#pragma omp critical
#ifdef MULTIPLE_OBJECT_LISTS
	objectList.get(omp_get_thread_num(), id, obj, &objectSpace);
#else
	objectList.get(id, obj, &objectSpace);
#endif
	float min = std::numeric_limits<float>::max();
	int minidx = -1;
	for (size_t cidx = 0; cidx != clusters.size(); cidx++) {
	  auto d = NGT::PrimitiveComparator::compareL2(reinterpret_cast<float*>(obj.data()),
						       clusters[cidx].centroid.data(), obj.size());
	  if (d < min) {
	    min = d;
	    minidx = cidx;
	  }
	}
	if (minidx < 0) {
	  std::cerr << "assign: Fatal error!" << std::endl;
	  abort();
	}
#pragma omp critical
	{
	  clusters[minidx].members.push_back(NGT::Clustering::Entry(id - 1, minidx, min));
	  count++;
	  if (count % 1000000 == 0) {
	    std::cerr << "# of assigned objects=" << count << std::endl;
	  }
	}
      }

    }

    static void assignWithNGT(std::vector<NGT::Clustering::Cluster> &clusters, size_t beginID, size_t endID,
			      NGT::ObjectSpace &objectSpace, QBGObjectList &objectList) {
      if (beginID > endID) {
	std::cerr << "assignWithNGT::Warning. beginID:" << beginID << " > endID:" << endID << std::endl;
	return;
      }

      NGT::Property prop;
      prop.dimension = objectSpace.getDimension();
      prop.objectType = NGT::Index::Property::ObjectType::Float;
      prop.distanceType = NGT::Property::DistanceType::DistanceTypeL2;
      prop.edgeSizeForCreation = 10;
      prop.edgeSizeForSearch = 40;

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      NGT::Index index(prop, "dummy");
#else
      NGT::Index index(prop);
#endif
      for (size_t cidx = 0; cidx < clusters.size(); cidx++) {
	if (cidx % 100000 == 0) {
	  std::cerr << "# of appended cluster objects=" << cidx << std::endl;
	}
	index.append(clusters[cidx].centroid);
      }
      std::cerr << "createIndex..." << std::endl;
      index.createIndex(500);

      std::cerr << "assign with NGT..." << std::endl;
      endID++;
#ifdef MULTIPLE_OBJECT_LISTS
      if (!objectList.openMultipleStreams(omp_get_max_threads())) {
	std::cerr << "Cannot open multiple streams." << std::endl;
	abort();
      }
#endif
      std::vector<pair<uint32_t, float>> clusterIDs(endID - beginID);
#pragma omp parallel for
      for (size_t id = beginID; id < endID; id++) {
	std::vector<float> obj;
#ifdef MULTIPLE_OBJECT_LISTS
	objectList.get(omp_get_thread_num(), id, obj, &objectSpace);
#else
	objectList.get(id, obj, &objectSpace);
#endif
	NGT::SearchQuery sc(obj);
	NGT::ObjectDistances	objects;
	sc.setResults(&objects);
	sc.setSize(10);
	sc.setEpsilon(0.12);
	index.search(sc);
	//index.linearSearch(sc);
	clusterIDs[id - beginID] = make_pair(objects[0].id - 1, objects[0].distance);
      }
      std::cerr << "pushing..." << std::endl;
      for (size_t id = beginID; id < endID; id++) {
	auto cid = clusterIDs[id - beginID].first;
	auto cdistance = clusterIDs[id - beginID].second;
	clusters[cid].members.push_back(NGT::Clustering::Entry(id - 1, cid, cdistance));
      }
    }

#ifdef NGTQ_QBG
    void treeBasedTopdownClustering(std::string prefix, QBG::Index &index, uint32_t rootID, std::vector<float> &object, std::vector<HKNode*> &nodes, NGT::Clustering &clustering) {
      auto &quantizer = static_cast<NGTQ::QuantizerInstance<uint8_t>&>(index.getQuantizer());
      auto &objectSpace = quantizer.globalCodebookIndex.getObjectSpace();
      QBGObjectList &objectList = quantizer.objectList;
      NGT::Timer timer;
      timer.start();
      std::vector<uint32_t> batch;
      std::vector<pair<uint32_t, uint32_t>> exceededLeaves;
      size_t nleaves = 1;
      size_t nOfThreads = 32;
      for (size_t id = 1; id <= numOfObjects; id++) {
	if (id % (numOfObjects / 100) == 0) {
	  timer.stop();
	  std::cerr << "# of processed objects=" << id << " " << id * 100 / numOfObjects << "% " << timer << " # of leaves=" << nleaves << std::endl;
	  timer.start();
	}
	batch.push_back(id);
	if (batch.size() > 100000) {
	  size_t kmeansBatchSize = nleaves < nOfThreads ? nleaves : nOfThreads;
	  hierarchicalKmeansBatch(batch, exceededLeaves, rootID, object, objectList, objectSpace, nodes, 
				  clustering, maxSize, nleaves, kmeansBatchSize);

	}
      }
      hierarchicalKmeansBatch(batch, exceededLeaves, rootID, object, objectList, objectSpace, nodes,
			      clustering, maxSize, nleaves, 0);

      if (numOfTotalClusters != 0) {
	NGT::Timer timer;
	timer.start();
	size_t numOfLeaves = 0;
	for (auto node : nodes) {
	  if (node->leaf) {
	    numOfLeaves++;
	  }
	}
	std::cerr << "# of nodes=" << nodes.size() << std::endl;
	std::cerr << "# of leaves=" << numOfLeaves << std::endl;
	std::cerr << "clustering for quantization." << std::endl;
	hierarchicalKmeansWithNumberOfClustersInParallel(numOfTotalClusters, numOfObjects, numOfLeaves, 
							 objectList, objectSpace, nodes, initMode);
	if (numOfTotalBlobs != 0) {
	  NGT::Timer timer;
	  timer.start();
	  size_t numOfLeaves = 0;
	  for (auto node : nodes) {
	    if (node->leaf) {
	      numOfLeaves++;
	    }
	  }
	  std::cerr << "# of leaves=" << numOfLeaves << ":" << numOfTotalClusters << std::endl;
	  if (numOfLeaves != numOfTotalClusters) {
	    std::cerr << "# of leaves is invalid " << numOfLeaves << ":" << numOfTotalClusters << std::endl;
	    abort();
	  }
	  {
	    std::ofstream of(prefix + QBG::Index::getSecondCentroidSuffix());
	    extractCentroids(of, nodes);
	  }
	  std::vector<uint32_t> qNodeIDs;
	  for (uint32_t nid = 0; nid < nodes.size(); nid++) {
	    if (nodes[nid]->leaf) {
	      qNodeIDs.push_back(nid);
	    }
	  }
	  std::cerr << "clustering to make blobs." << std::endl;
	  hierarchicalKmeansWithNumberOfClustersInParallel(numOfTotalBlobs, numOfObjects, numOfTotalClusters, 
							   objectList, objectSpace, nodes, initMode);
	  {
	    std::ofstream of(prefix + QBG::Index::get3rdTo2ndSuffix());
	    extractBtoQIndex(of, nodes, qNodeIDs);
	  }
	}
      } 

    }

    void multilayerClustering(std::string prefix, QBG::Index &index) {
      auto &quantizer = static_cast<NGTQ::QuantizerInstance<uint8_t>&>(index.getQuantizer());
      auto &objectSpace = quantizer.globalCodebookIndex.getObjectSpace();
      {
	std::cerr << "Three layer clustering..." << std::endl;
	std::cerr << "HiearchicalKmeans::clustering: # of clusters=" << numOfThirdClusters << ":" << index.getQuantizer().property.globalCentroidLimit << std::endl;
	if (index.getQuantizer().objectList.size() <= 1) {
	  NGTThrowException("optimize: No objects");
	}
	if (numOfThirdClusters == 0) {
	  if (index.getQuantizer().property.globalCentroidLimit == 0) {
	    numOfThirdClusters = index.getQuantizer().objectList.size() / 1000;
	    numOfThirdClusters = numOfThirdClusters == 0 ? 1 : numOfThirdClusters;
	    numOfThirdClusters = numOfThirdClusters > 1000000 ? 1000000 : numOfThirdClusters;
	  } else {
	    numOfThirdClusters = index.getQuantizer().property.globalCentroidLimit;
	  }
	}
	if (numOfThirdClusters != 0 && index.getQuantizer().property.globalCentroidLimit != 0 &&
	    numOfThirdClusters != index.getQuantizer().property.globalCentroidLimit) {
	}
	auto &quantizer = static_cast<NGTQ::QuantizerInstance<uint8_t>&>(index.getQuantizer());
	QBGObjectList &objectList = quantizer.objectList;
	if (numOfObjects == 0) {
	  numOfObjects = objectList.size() - 1;
	}

	std::cerr << "The first layer. " << numOfFirstClusters << ":" << numOfFirstObjects << std::endl;
	if (numOfThirdClusters == 0 || numOfObjects == 0) {
	  NGTThrowException("numOfThirdClusters or numOfObjects are zero");
	}
	numOfSecondClusters = numOfSecondClusters == 0 ? numOfThirdClusters : numOfSecondClusters;
	numOfFirstClusters = numOfFirstClusters == 0 ? static_cast<size_t>(sqrt(numOfSecondClusters)) : numOfFirstClusters;
	numOfSecondObjects = numOfSecondClusters * 100;
	numOfSecondObjects = numOfSecondObjects > numOfObjects ? numOfObjects : numOfSecondObjects;
	numOfFirstObjects = numOfFirstClusters * 2000;
	numOfFirstObjects = numOfFirstObjects > numOfSecondObjects ? numOfSecondObjects : numOfFirstObjects;
	if (numOfFirstObjects < numOfFirstClusters) {
	  std::stringstream msg;
	  msg << "# of objects for the first should be larger than # of the first clusters. " << numOfFirstObjects << ":" << numOfFirstClusters;
	  NGTThrowException(msg);
	}
	if (numOfFirstClusters > numOfSecondClusters) {
	  std::stringstream msg;
	  msg << "# of the first clusters should be larger than or equal to # of the second clusters. " << numOfFirstClusters << ":" << numOfSecondClusters;
	  NGTThrowException(msg);
	}
	if (numOfSecondClusters > numOfThirdClusters) {
	  std::stringstream msg;
	  msg << "# of the third clusters should be larger than or equal to # of the second clusters. " << numOfSecondClusters << ":" << numOfThirdClusters;
	  NGTThrowException(msg);
	}
	if (numOfFirstClusters > numOfSecondClusters) {
	  std::stringstream msg;
	  msg << "# of the second clusters should be larger than # of the first clusters. " << numOfFirstClusters << ":" << numOfSecondClusters;
	  NGTThrowException(msg);
	}

	std::cerr << "Three layer clustering:" << numOfFirstClusters << ":" << numOfFirstObjects << "," << numOfSecondClusters << ":" << numOfSecondObjects << "," << numOfThirdClusters << ":" << numOfObjects << std::endl;

	NGT::Clustering firstClustering(initMode, NGT::Clustering::ClusteringTypeKmeansWithoutNGT, 300);
	float clusterSizeConstraint = 5.0;
	firstClustering.setClusterSizeConstraintCoefficient(clusterSizeConstraint);
	std::cerr << "size constraint=" << clusterSizeConstraint << std::endl;
	std::vector<std::vector<float>> vectors;
	vectors.reserve(numOfFirstObjects);
	std::vector<float> obj;
	for (size_t id = 1; id <= numOfFirstObjects; id++) {
	  if (id % 1000000 == 0) {
	    std::cerr << "# of prcessed objects is " << id << std::endl;
	  }
	  if (!objectList.get(id, obj, &objectSpace)) {
	    std::stringstream msg;
	    msg << "qbg: Cannot get object. ID=" << id;
	    NGTThrowException(msg);
	  }
	  vectors.push_back(obj);
	}
	std::cerr << "Kmeans... " << vectors.size() << " to " << numOfFirstClusters << std::endl;
	std::vector<NGT::Clustering::Cluster> firstClusters;
	NGT::Timer timer;

	timer.start();
	firstClustering.kmeans(vectors, numOfFirstClusters, firstClusters);
	timer.stop();
	std::cerr << "# of clusters=" << firstClusters.size() << " time=" << timer << std::endl;

	std::vector<std::vector<float>> otherVectors;
	timer.start();
	std::cerr << "Assign for the second. (" << numOfFirstObjects << "-" << numOfSecondObjects << ")..." << std::endl;
	assign(firstClusters, numOfFirstObjects + 1, numOfSecondObjects, objectSpace, objectList);
	timer.stop();
	std::cerr << "Assign(1) time=" << timer << std::endl;

	std::cerr << "subclustering for the second." << std::endl;
	std::vector<NGT::Clustering::Cluster> secondClusters;
	timer.start();
	subclustering(firstClusters, numOfSecondClusters, numOfSecondObjects, objectSpace, objectList, initMode, secondClusters);
	timer.stop();
	std::cerr << "subclustering(1) time=" << timer << std::endl;
	std::cerr << "save quantization centroid" << std::endl;
	NGT::Clustering::saveClusters(prefix + QBG::Index::getSecondCentroidSuffix(), secondClusters);
	timer.start();
	std::cerr << "Assign for the third. (" << numOfSecondObjects << "-" << numOfObjects << ")..." << std::endl;
	assignWithNGT(secondClusters, numOfSecondObjects + 1, numOfObjects, objectSpace, objectList);
	timer.stop();
	std::cerr << "Assign(2) time=" << timer << std::endl;
	std::cerr << "subclustering for the third." << std::endl;
	std::vector<std::vector<NGT::Clustering::Cluster>> thirdClusters;
	timer.start();
	subclustering(secondClusters, numOfThirdClusters, numOfObjects, objectSpace, objectList, initMode, thirdClusters);
	timer.stop();
	std::cerr << "subclustering(2) time=" << timer << std::endl;
	{
	  std::vector<size_t> bqindex;
	  for (size_t idx1 = 0; idx1 < thirdClusters.size(); idx1++) {
	    for (size_t idx2 = 0; idx2 < thirdClusters[idx1].size(); idx2++) {
	      bqindex.push_back(idx1);
	    }
	  }
	  std::cerr << "save bqindex..." << std::endl;
	  NGT::Clustering::saveVector(prefix + QBG::Index::get3rdTo2ndSuffix(), bqindex);
	}
	
	std::vector<NGT::Clustering::Cluster> thirdFlatClusters;
	flattenClusters(secondClusters, thirdClusters, numOfThirdClusters, thirdFlatClusters);

	std::cerr << "save centroid..." << std::endl;
	NGT::Clustering::saveClusters(prefix + QBG::Index::getThirdCentroidSuffix(), thirdFlatClusters);

	{
	  std::vector<size_t> cindex(numOfObjects);
	  for (size_t cidx = 0; cidx < thirdFlatClusters.size(); cidx++) {
	    for (auto mit = thirdFlatClusters[cidx].members.begin(); mit != thirdFlatClusters[cidx].members.end(); ++mit) {
	      size_t vid = (*mit).vectorID;
	      cindex[vid] = cidx;
	    }
	  }
	  std::cerr << "save index... " << cindex.size() << std::endl;
	  NGT::Clustering::saveVector(prefix + QBG::Index::getObjTo3rdSuffix(), cindex);
	}
	std::cerr << "end of clustering" << std::endl;
	return;
      }

    }

    void clustering(std::string indexPath, std::string prefix = "", std::string objectIDsFile = "") {
      NGT::StdOstreamRedirector redirector(silence);
      redirector.begin();

      bool readOnly = false;
      QBG::Index index(indexPath, readOnly);
      index.getQuantizer().objectList.size();
      std::cerr << "clustering... " << std::endl;
      if (threeLayerClustering) {
	try {
	  if (numOfObjects == 0) {
	    numOfObjects = index.getQuantizer().objectList.size() - 1;
	  }
	  if (numOfObjects != index.getQuantizer().objectList.size() - 1) {
	    std::cerr << "HierarchicalKmeans::clustering: Warning! # of objects is invalid." << std::endl;
	    std::cerr << "     " << index.getQuantizer().objectList.size() - 1 << " is set to # of object instead of " << numOfObjects << std::endl;
	    numOfObjects = index.getQuantizer().objectList.size() - 1;
	  }
	  if (prefix.empty()) {
	    std::cerr << "Prefix is not specified." << std::endl;
	    prefix = indexPath + "/" + QBG::Index::getWorkspaceName();
	    try {
	      NGT::Index::mkdir(prefix);
	    } catch(...) {}
	    prefix +="/" + QBG::Index::getHierarchicalClusteringPrefix();
	    std::cerr << prefix << " is used" << std::endl;
	  }
	  auto &quantizer = static_cast<NGTQ::QuantizerInstance<uint8_t>&>(index.getQuantizer());
	  auto &objectSpace = quantizer.globalCodebookIndex.getObjectSpace();
	  size_t paddedDimension = objectSpace.getPaddedDimension();
	  size_t dimension = objectSpace.getDimension();
	  if (paddedDimension != dimension) {
	    std::cerr << "HierarachicalKmeans: Warning! Dimensions are inconsistent. Dimension=" << paddedDimension << ":" << dimension << std::endl;
	  }
	  multilayerClustering(prefix, index);
	} catch(NGT::Exception &err) {
	  redirector.end();
	  throw err;
	}
	redirector.end();
	return;
      }
      try {
	NGT::Clustering::ClusteringType clusteringType = NGT::Clustering::ClusteringTypeKmeansWithoutNGT;

	uint32_t rootID = 0;
	std::vector<HKNode*> nodes;
	nodes.push_back(new HKLeafNode);

	std::vector<float> object;
	size_t iteration = 1000;
	NGT::Clustering clustering(initMode, clusteringType, iteration, numOfClusters);
	auto &quantizer = static_cast<NGTQ::QuantizerInstance<uint8_t>&>(index.getQuantizer());
	QBGObjectList &objectList = quantizer.objectList;
	if (objectIDsFile.empty()) {
	  treeBasedTopdownClustering(prefix, index, rootID, object, nodes, clustering);
	} else {
	  std::cerr << "Cluster ID=" << clusterID << std::endl;
	  if (clusterID < 0) {
	    std::stringstream msg;
	    msg << "Any target cluster ID is not specified.";
	    NGTThrowException(msg);
	  }
	  std::ifstream objectIDs(objectIDsFile);
	  if (!objectIDs) {
	    std::stringstream msg;
	    msg << "Cannot open the object id file. " << objectIDsFile;
	    NGTThrowException(msg);
	  }
	  auto &objectSpace = quantizer.globalCodebookIndex.getObjectSpace();
	  uint32_t id = 1;
	  int32_t cid;
	  size_t ccount = 0;
	  while (objectIDs >> cid) {
	    std::cerr << cid << std::endl;
	    if (id % 100000 == 0) {
	      std::cerr << "# of processed objects=" << id << std::endl;
	    }
	    if (cid == -1) {
	      continue;
	    }
	    if (cid == clusterID) {
	      ccount++;
	      hierarchicalKmeans(id, rootID, object, objectList, objectSpace, nodes, clustering, maxSize);
	    }
	    id++;
	  }
	}
	size_t objectCount = 0;
	if (prefix.empty()) {
	  objectCount = extractCentroids(std::cout, nodes);
	} else {
	  {
	    std::ofstream of(prefix + QBG::Index::getThirdCentroidSuffix());
	    objectCount = extractCentroids(of, nodes);
	  }
	  {
	    std::ofstream of(prefix + QBG::Index::getObjTo3rdSuffix());
	    extractIndex(of, nodes, numOfObjects);
	  }
	  if (numOfFirstObjects > 0) {
	    std::ofstream btoqof(prefix + QBG::Index::get3rdTo2ndSuffix());
	    std::ofstream qcof(prefix + QBG::Index::getSecondCentroidSuffix());
	    extractBtoQAndQCentroid(btoqof, qcof, nodes, numOfThirdClusters);
	  }
	  if (numOfRandomObjects > 0) {
	    std::ofstream of(prefix + "_random_object.tsv");
	    if (extractCentroid) {
	      extractRandomObjectsFromEachBlob(of, nodes, numOfObjects, numOfRandomObjects - 1, quantizer, extractCentroid);
	    } else {
	      extractRandomObjectsFromEachBlob(of, nodes, numOfObjects, numOfRandomObjects, quantizer, extractCentroid);
	    }
	  }
	}
	if (objectCount != numOfObjects) {
	  std::cerr << "# of objects is invalid. " << objectCount << ":" << numOfObjects << std::endl;
	}
      } catch(NGT::Exception &err) {
	redirector.end();
	throw err;
      }
      redirector.end();
    }
#endif 

    size_t	maxSize;
    size_t	numOfObjects;
    size_t	numOfClusters;
    size_t	numOfTotalClusters;
    size_t	numOfTotalBlobs;
    int32_t	clusterID;

    NGT::Clustering::InitializationMode initMode;

    size_t	numOfRandomObjects;

    size_t	numOfFirstObjects;
    size_t	numOfFirstClusters;
    size_t	numOfSecondObjects;
    size_t	numOfSecondClusters;
    size_t	numOfThirdClusters;
    bool	extractCentroid;

    bool	threeLayerClustering;
    bool	silence;
  };
} 
