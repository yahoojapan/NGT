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


#include "NGT/NGTQ/QuantizedBlobGraph.h"
#include "NGT/NGTQ/HierarchicalKmeans.h"

QBG::HierarchicalKmeans::HierarchicalKmeans(QBG::BuildParameters &param) {
#ifdef NGTQ_QBG
  setParameters(param.hierarchicalClustering);
#endif
}
QBG::HierarchicalKmeans::HierarchicalKmeans(QBG::HierarchicalClusteringParameters &hierarchicalClustering) {
  setParameters(hierarchicalClustering);
}

void QBG::HierarchicalKmeans::setParameters(QBG::HierarchicalClusteringParameters &hierarchicalClustering) {
#ifdef NGTQ_QBG
  maxSize		= hierarchicalClustering.maxSize;
  numOfObjects		= hierarchicalClustering.numOfObjects;
  numOfClusters		= hierarchicalClustering.numOfClusters;
  numOfTotalClusters	= hierarchicalClustering.numOfTotalClusters;
  numOfTotalBlobs	= hierarchicalClustering.numOfTotalBlobs;
  clusterID		= hierarchicalClustering.clusterID;

  initMode		= hierarchicalClustering.initMode;

  numOfRandomObjects	= hierarchicalClustering.numOfRandomObjects;

  numOfFirstObjects	= hierarchicalClustering.numOfFirstObjects;
  numOfFirstClusters	= hierarchicalClustering.numOfFirstClusters;
  numOfSecondObjects	= hierarchicalClustering.numOfSecondObjects;
  numOfSecondClusters	= hierarchicalClustering.numOfSecondClusters;
  numOfThirdClusters	= hierarchicalClustering.numOfThirdClusters;
  numOfThirdObjects	= hierarchicalClustering.numOfThirdObjects;
  extractCentroid	= hierarchicalClustering.extractCentroid;

  clusteringType	= hierarchicalClustering.clusteringType;
  epsilonExplorationSize = hierarchicalClustering.epsilonExplorationSize;
  expectedRecall	= hierarchicalClustering.expectedRecall;
  verbose		= hierarchicalClustering.verbose;
#endif
}

void QBG::HierarchicalKmeans::initialize() {
#ifdef NGTQ_QBG
  HierarchicalClusteringParameters params;
  setParameters(params);
#endif
}

#ifdef NGTQ_QBG
void QBG::HierarchicalKmeans::treeBasedTopdownClustering(std::string prefix, QBG::Index &index, uint32_t rootID, std::vector<float> &object, std::vector<HKNode*> &nodes, NGT::Clustering &clustering) {
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

void QBG::HierarchicalKmeans::threeLayerClustering(std::string prefix, QBG::Index &index) {
  auto &quantizer = static_cast<NGTQ::QuantizerInstance<uint8_t>&>(index.getQuantizer());
  auto &objectSpace = quantizer.globalCodebookIndex.getObjectSpace();
  {
    std::cerr << "Three layer clustering..." << std::endl;
    std::cerr << "HiearchicalKmeans::clustering: # of clusters=" << numOfThirdClusters << ":" << index.getQuantizer().property.globalCentroidLimit << std::endl;
    if (index.getQuantizer().objectList.size() <= 1) {
      NGTThrowException("HierarchicelKmeans: No objects");
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
    if (numOfThirdObjects > numOfObjects) {
      numOfThirdObjects = numOfObjects;
    }

    if (numOfThirdClusters == 0 || numOfObjects == 0) {
      NGTThrowException("numOfThirdClusters or numOfObjects are zero");
    }
    numOfThirdObjects = numOfThirdObjects == 0 ? numOfObjects : numOfThirdObjects;
    numOfSecondClusters = numOfSecondClusters == 0 ? numOfThirdClusters : numOfSecondClusters;
    numOfFirstClusters = numOfFirstClusters == 0 ? static_cast<size_t>(sqrt(numOfSecondClusters)) : numOfFirstClusters;
    numOfSecondObjects = numOfSecondObjects == 0 ? numOfSecondClusters * 100 : numOfSecondObjects;
    numOfSecondObjects = numOfSecondObjects > numOfObjects ? numOfObjects : numOfSecondObjects;
    numOfFirstObjects = numOfFirstObjects == 0 ? numOfFirstClusters * 2000 : numOfFirstObjects;
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

    NGT::Clustering firstClustering(initMode, NGT::Clustering::ClusteringTypeKmeansWithoutNGT, 300);
    float clusterSizeConstraint = 5.0;
    firstClustering.setClusterSizeConstraintCoefficient(clusterSizeConstraint);
    std::cerr << "size constraint=" << clusterSizeConstraint << std::endl;
    std::vector<std::vector<float>> vectors;
    vectors.reserve(numOfFirstObjects);
    std::cerr << "loading objects..." << std::endl;
    NGT::Timer timer;
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
    std::cerr << "loading objects time=" << timer << std::endl;
    std::cerr << "clustering for the first... (" << vectors.size() << "->" << numOfFirstClusters << ") " << std::endl;
    std::vector<NGT::Clustering::Cluster> firstClusters;

    timer.start();
    firstClustering.kmeans(vectors, numOfFirstClusters, firstClusters);
    timer.stop();
    std::cerr << "end of clustering for the first. # of clusters=" << firstClusters.size() << " time=" << timer << std::endl;

    std::vector<std::vector<float>> otherVectors;
    timer.start();
    std::cerr << "assign for the second. (" << numOfFirstObjects << "->" << numOfSecondObjects << ")..." << std::endl;
    assign(firstClusters, numOfFirstObjects + 1, numOfSecondObjects, objectSpace, objectList);
    timer.stop();
    std::cerr << "end of assign for the second. time=" << timer << ", vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;

    std::cerr << "subclustering for the second (" << numOfSecondClusters << ")..." << std::endl;
    std::vector<NGT::Clustering::Cluster> secondClusters;
    timer.start();
    subclustering(firstClusters, numOfSecondClusters, numOfSecondObjects, objectSpace, objectList, initMode, secondClusters);
    timer.stop();
    std::cerr << "end of subclustering for the second. time=" << timer << ", vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;
    timer.start();
    NGT::Clustering::clearMembers(secondClusters);
    std::cerr << "assign for the third. (" << 1 << "->" << numOfThirdObjects << ")..." << std::endl;
    assignWithNGT(secondClusters, 1, numOfThirdObjects, objectSpace, objectList, epsilonExplorationSize, expectedRecall);
    {
      size_t noOfRemovedClusters = NGT::Clustering::removeEmptyClusters(secondClusters);
      if (noOfRemovedClusters != 0) {
	std::cerr << "Clustering: Warning. # of removed clusters=" << noOfRemovedClusters << std::endl;
      }
    }
    timer.stop();
    std::cerr << "end of assign for the third. time=" << timer << ", vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;
    std::cerr << "subclustering for the third (" << numOfThirdClusters << ")..." << std::endl;
    std::vector<std::vector<NGT::Clustering::Cluster>> thirdClusters;
    timer.start();
    subclustering(secondClusters, numOfThirdClusters, numOfThirdObjects, objectSpace, objectList, initMode, thirdClusters);
    timer.stop();
    std::cerr << "end of subclustering for the third. time=" << timer << ", vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;

    std::vector<NGT::Clustering::Cluster> thirdFlatClusters;
    flattenClusters(secondClusters, thirdClusters, numOfThirdClusters, thirdFlatClusters);

    timer.start();
    NGT::Clustering::clearMembers(thirdFlatClusters);
    std::cerr << "assign all for the third (" << 1 << "-" << numOfObjects << ")..." << std::endl;
    assignWithNGT(thirdFlatClusters, 1, numOfObjects, objectSpace, objectList, epsilonExplorationSize, expectedRecall);
    timer.stop();
    std::cerr << "end of assign all for the third. time=" << timer << ", vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;

    {
      std::vector<size_t> bqindex;
      size_t idx = 0;
      for (size_t idx1 = 0; idx1 < thirdClusters.size(); idx1++) {
	for (size_t idx2 = 0; idx2 < thirdClusters[idx1].size(); idx2++, idx++) {
	  if (thirdClusters[idx1][idx2].members.size() == 0) {
	    std::stringstream msg;
	    msg << "Fatal error! found an empty cluster in thirdClusters.";
	    NGTThrowException(msg);
	  }
	  if (thirdFlatClusters[idx].members.size() == 0) {
	    std::cerr << "warning. found an empty cluster in thirdFlatClusters. " << idx << std::endl;
	  } else {
	    bqindex.push_back(idx1);
	  }
	}
      }
      std::cerr << "save the 3rd to the 2nd index..." << std::endl;
      NGT::Clustering::saveVector(prefix + QBG::Index::get3rdTo2ndSuffix(), bqindex);
    }

    std::cerr << "save quantization centroid" << std::endl;
    NGT::Clustering::saveClusters(prefix + QBG::Index::getSecondCentroidSuffix(), secondClusters);

    std::cerr << "save the third centroid..." << std::endl;
    auto skipEmptyClusters = true;
    NGT::Clustering::saveClusters(prefix + QBG::Index::getThirdCentroidSuffix(), thirdFlatClusters, skipEmptyClusters);
    {
      std::vector<size_t> cindex(numOfObjects);
      size_t idx = 0;
      for (size_t cidx = 0; cidx < thirdFlatClusters.size(); cidx++) {
	if (thirdFlatClusters[cidx].members.size() == 0) {
	  continue;
	}
	for (auto mit = thirdFlatClusters[cidx].members.begin(); mit != thirdFlatClusters[cidx].members.end(); ++mit) {
	  size_t vid = (*mit).vectorID;
	  cindex[vid] = idx;
	}
	idx++;
      }
      std::cerr << "save index... " << cindex.size() << std::endl;
      NGT::Clustering::saveVector(prefix + QBG::Index::getObjTo3rdSuffix(), cindex);
    }
    std::cerr << "end of clustering" << std::endl;
    return;
  }

}

void QBG::HierarchicalKmeans::twoPlusLayerClustering(std::string prefix, QBG::Index &index) {
  auto &quantizer = static_cast<NGTQ::QuantizerInstance<uint8_t>&>(index.getQuantizer());
  auto &objectSpace = quantizer.globalCodebookIndex.getObjectSpace();
  {
    std::cerr << "Two layer clustering..." << std::endl;
    std::cerr << "HiearchicalKmeans::clustering: # of clusters=" << numOfThirdClusters << ":" << index.getQuantizer().property.globalCentroidLimit << std::endl;
    if (index.getQuantizer().objectList.size() <= 1) {
      NGTThrowException("HierarchicelKmeans: No objects");
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
    if (numOfThirdObjects > numOfObjects) {
      numOfThirdObjects = numOfObjects;
    }

    std::cerr << "The first layer. " << numOfFirstClusters << ":" << numOfFirstObjects << std::endl;
    if (numOfThirdClusters == 0 || numOfObjects == 0) {
      NGTThrowException("numOfThirdClusters or numOfObjects are zero");
    }
    numOfThirdObjects = numOfThirdObjects == 0 ? numOfObjects : numOfThirdObjects;
    numOfSecondClusters = numOfSecondClusters == 0 ? numOfThirdClusters : numOfSecondClusters;
    numOfFirstClusters = numOfFirstClusters == 0 ? static_cast<size_t>(sqrt(numOfSecondClusters)) : numOfFirstClusters;
    numOfSecondObjects = numOfSecondObjects == 0 ? numOfSecondClusters * 100 : numOfSecondObjects;
    numOfSecondObjects = numOfSecondObjects > numOfObjects ? numOfObjects : numOfSecondObjects;
    numOfFirstObjects = numOfFirstObjects == 0 ? numOfFirstClusters * 2000 : numOfFirstObjects;
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

    std::cerr << "Two layer clustering. " << numOfFirstClusters << ":" << numOfFirstObjects << "," << numOfSecondClusters << ":" << numOfSecondObjects << "," << numOfThirdClusters << ":" << numOfThirdObjects << " " << numOfObjects << std::endl;

    NGT::Clustering firstClustering(initMode, NGT::Clustering::ClusteringTypeKmeansWithoutNGT, 300);
    float clusterSizeConstraint = 5.0;
    firstClustering.setClusterSizeConstraintCoefficient(clusterSizeConstraint);
    std::cerr << "size constraint=" << clusterSizeConstraint << std::endl;
    std::vector<std::vector<float>> vectors;
    size_t reserveSize = numOfThirdClusters > numOfFirstObjects ? numOfThirdClusters : numOfFirstObjects;
    vectors.reserve(reserveSize);
    std::cerr << "loading objects..." << std::endl;
    NGT::Timer timer;
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
    std::cerr << "loading objects time=" << timer << " vmsize=" << NGT::Common::getProcessVmSize() << std::endl;
    std::cerr << "kmeans clustering... " << vectors.size() << " to " << numOfFirstClusters << std::endl;
    std::vector<NGT::Clustering::Cluster> firstClusters;

    timer.start();
    firstClustering.kmeans(vectors, numOfFirstClusters, firstClusters);
    timer.stop();
    std::cerr << "kmeans clustering. # of clusters=" << firstClusters.size() << " time=" << timer << " vmsize=" << NGT::Common::getProcessVmSize() << std::endl;
    timer.start();
    std::cerr << "assign for the third (" << numOfFirstObjects << "-" << numOfThirdObjects << ")..." << std::endl;
    assign(firstClusters, numOfFirstObjects + 1, numOfThirdObjects, objectSpace, objectList);
    timer.stop();
    std::cerr << "assign for the third. time=" << timer << " vmsize=" << NGT::Common::getProcessVmSize() << std::endl;

    std::cerr << "subclustering for the third (" << numOfThirdClusters << ", " << numOfThirdObjects << ")..." << std::endl;
    std::vector<NGT::Clustering::Cluster> thirdFlatClusters;
    timer.start();
    subclustering(firstClusters, numOfThirdClusters, numOfThirdObjects, objectSpace, objectList, initMode, thirdFlatClusters);
    timer.stop();
    std::cerr << "subclustering for the third. time=" << timer << ", vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;

    std::cerr << "assign all for the third (" << numOfThirdObjects << "-" << numOfObjects << ")..." << std::endl;
    timer.start();
    assignWithNGT(thirdFlatClusters, numOfThirdObjects + 1, numOfObjects, objectSpace, objectList);
    timer.stop();
    std::cerr << "assign all for the third. time=" << timer << ", vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;

    std::vector<NGT::Clustering::Cluster> secondClusters;

    std::cerr << "clustering for the second. # of objects=" << thirdFlatClusters.size() << " # of clusters=" << numOfSecondClusters << std::endl;
    timer.start();
    vectors.clear();
    vectors.reserve(thirdFlatClusters.size());
    for (auto cit = thirdFlatClusters.begin(); cit != thirdFlatClusters.end(); ++cit) {
      std::vector<float> &v = (*cit).centroid;
      vectors.push_back(v);
    }
    if (clusteringType == QBG::HierarchicalKmeans::ClusteringTypeTwoPlusOneLayer) {
      std::cerr << "two + one layer clustering..." << std::endl;
      NGT::Clustering clustering(initMode, NGT::Clustering::ClusteringTypeKmeansWithoutNGT, 1000);
      clustering.kmeans(vectors, numOfSecondClusters, secondClusters);
    } else if (clusteringType == QBG::HierarchicalKmeans::ClusteringTypeTwoPlusOneLayer) {
      auto n = numOfSecondObjects / numOfSecondClusters;
      std::cerr << "two + one layer clustering with NGT... " << n << std::endl;
      NGT::Clustering clustering(initMode, NGT::Clustering::ClusteringTypeKmeansWithNGT, n);
      clustering.kmeans(vectors, numOfSecondClusters, secondClusters);
    } else if (clusteringType == QBG::HierarchicalKmeans::ClusteringTypeTwoPlusTwoLayer) {
      std::cerr << "two + two layer clustering..." << std::endl;
      twoLayerClustering(vectors, numOfSecondClusters, secondClusters, 100);
    } else {
      std::stringstream msg;
      msg << "Invalid clustering type:" << clusteringType << std::endl;
      NGTThrowException(msg);
    }
    timer.stop();
    std::cerr << "clustering for the second. time=" << timer << " vmsize=" << NGT::Common::getProcessVmSize() << std::endl;
    NGT::Clustering::saveClusters(prefix + QBG::Index::getSecondCentroidSuffix(), secondClusters);

    timer.start();
    std::vector<NGT::Clustering::Cluster> thirdPermutedClusters;
    thirdPermutedClusters.reserve(thirdFlatClusters.size());
    for (size_t sidx = 0; sidx < secondClusters.size(); sidx++) {
      for (auto mit = secondClusters[sidx].members.begin(); mit != secondClusters[sidx].members.end(); ++mit) {
	thirdPermutedClusters.emplace_back(thirdFlatClusters[(*mit).vectorID]);
	(*mit).vectorID = thirdPermutedClusters.size() - 1;
      }
    }
    timer.stop();
    std::cerr << "permuting cluster time=" << timer << " vmsize=" << NGT::Common::getProcessVmSize() << std::endl;

    std::vector<size_t> cindex(numOfObjects);
    for (size_t cidx = 0; cidx < thirdPermutedClusters.size(); cidx++) {
      for (auto mit = thirdPermutedClusters[cidx].members.begin(); mit != thirdPermutedClusters[cidx].members.end(); ++mit) {
	size_t vid = (*mit).vectorID;
	cindex[vid] = cidx;
      }
    }
    std::cerr << "save index... " << cindex.size() << std::endl;
    NGT::Clustering::saveVector(prefix + QBG::Index::getObjTo3rdSuffix(), cindex);

    std::vector<size_t> bqindex(thirdPermutedClusters.size());
    for (size_t idx1 = 0; idx1 < secondClusters.size(); idx1++) {
      for (size_t idx2 = 0; idx2 < secondClusters[idx1].members.size(); idx2++) {
	bqindex[secondClusters[idx1].members[idx2].vectorID] = idx1;
      }
    }
    std::cerr << "save bqindex..." << std::endl;
    NGT::Clustering::saveVector(prefix + QBG::Index::get3rdTo2ndSuffix(), bqindex);
    std::cerr << "save the third  centroid" << std::endl;
    NGT::Clustering::saveClusters(prefix + QBG::Index::getThirdCentroidSuffix(), thirdPermutedClusters);

  }

  std::cerr << "end of clustering" << std::endl;
  return;
}

void QBG::HierarchicalKmeans::multiLayerClustering(QBG::Index &index, std::string prefix, std::string objectIDsFile) {
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
}

void QBG::HierarchicalKmeans::clustering(std::string indexPath, std::string prefix, std::string objectIDsFile) {
  NGT::StdOstreamRedirector redirector(!verbose);
  redirector.begin();

  std::cerr << "The specified params=FC:" << numOfFirstClusters << ":FO:" << numOfFirstObjects
	    << ",SC:" << numOfSecondClusters << ":SO:" << numOfSecondObjects
	    << ",TC:" << numOfThirdClusters << ":TO:" << numOfThirdObjects << ",O:" << numOfObjects << std::endl;

  bool readOnly = false;
  QBG::Index index(indexPath, readOnly);
  if (index.getQuantizer().objectList.size() <= 1) {
    NGTThrowException("No objects in the index.");
  }
  if (clusteringType == QBG::HierarchicalKmeans::ClusteringTypeMultiLayer) {
    try {
      multiLayerClustering(index, prefix, objectIDsFile);
    } catch(NGT::Exception &err) {
      redirector.end();
      throw err;
    }
  } else {
    std::cerr << "three layer clustering... " << std::endl;
    try {
      if (numOfObjects == 0) {
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
      if (clusteringType == QBG::HierarchicalKmeans::ClusteringTypeThreeLayer) {
	threeLayerClustering(prefix, index);
      } else {
	twoPlusLayerClustering(prefix, index);
      }
    } catch(NGT::Exception &err) {
      redirector.end();
      throw err;
    }
  }
  redirector.end();
}

#endif

