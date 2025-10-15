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

void QBG::Index::append(const std::string &indexName, // index file
                        const std::string &data,      // data file
                        size_t dataSize,              // data size
                        bool verbose) {
  NGT::StdOstreamRedirector redirector(!verbose);
  redirector.begin();
  QBG::Index index(indexName);
  auto &quantizer = index.getQuantizer();
  istream *is;
  if (data == "-") {
    is = &cin;
  } else {
    ifstream *ifs = new ifstream;
    ifs->ifstream::open(data);
    if (!(*ifs)) {
      std::stringstream msg;
      msg << "Cannot open the specified file. " << data;
      NGTThrowException(msg);
    }
    is = ifs;
  }
  string line;
  size_t idx   = quantizer.objectList.size() == 0 ? 0 : quantizer.objectList.size() - 1;
  size_t count = 0;
  // extract objects from the file and insert them to the object list.
  while (getline(*is, line)) {
    idx++;
    count++;
    std::vector<float> object;
    NGT::Common::extractVector(line, " ,\t", object);
    if (object.empty()) {
      cerr << "Empty line or invalid value: " << line << endl;
      continue;
    }
    if ((quantizer.property.distanceType == NGTQ::DistanceType::DistanceTypeInnerProduct) &&
        (object.size() + 1 == quantizer.objectList.genuineDimension)) {
      object.emplace_back(0);
    }
    if (object.size() != quantizer.objectList.genuineDimension) {
      std::stringstream msg;
      msg << "The dimension of the specified object is inconsistent with the dimension of the index. "
          << object.size() << ":" << quantizer.objectList.genuineDimension;
      NGTThrowException(msg);
    }
    index.insert(idx, object);

    if (count % 100000 == 0) {
      std::cerr << "appended " << static_cast<float>(count) / 1000000.0 << "M objects.";
      if (count != idx) {
        std::cerr << " # of the total objects=" << static_cast<float>(idx) / 1000000.0 << "M";
      }
      cerr << " virtual memory(kbyte)=" << NGT::Common::getProcessVmSize() << std::endl;
    }
  }
  if (data != "-") {
    delete is;
  }

  index.save();
  index.close();
  redirector.end();
}

void QBG::Index::append(const std::string &indexName,  // index file
                        NGT::ObjectSpace &objectSpace, // object space including objects
                        bool verbose) {
  NGT::StdOstreamRedirector redirector(!verbose);
  redirector.begin();
  QBG::Index index(indexName);
  auto &quantizer = index.getQuantizer();
  auto &repo      = objectSpace.getRepository();

  size_t count = 0;
  // extract objects from the file and insert them to the object list.
  for (size_t idx = 1; idx < repo.size(); idx++) {
    std::vector<float> object;
    object.clear();
    if (repo.isEmpty(idx)) {
      continue;
    }
    objectSpace.getObject(idx, object);
    count++;
    if (object.size() != quantizer.objectList.genuineDimension) {
      std::stringstream msg;
      msg << "The dimension of the specified object is inconsistent with the dimension of the index. "
          << object.size() << ":" << quantizer.objectList.genuineDimension;
      NGTThrowException(msg);
    }
    index.insert(idx, object);

    if (count % 100000 == 0) {
      std::cerr << "appended " << static_cast<float>(count) / 1000000.0 << "M objects.";
      if (count != idx) {
        std::cerr << " # of the total objects=" << static_cast<float>(idx) / 1000000.0 << "M";
      }
      cerr << " virtual memory(kbyte)=" << NGT::Common::getProcessVmSize() << std::endl;
    }
  }

  redirector.end();
}

void QBG::Index::appendBinary(const std::string &indexName, // index file
                              const std::string &data,      // data file
                              size_t dataSize,              // data size
                              bool verbose) {
  NGT::StdOstreamRedirector redirector(!verbose);
  redirector.begin();
  QBG::Index index(indexName);
  std::vector<std::string> tokens;
  NGT::Common::tokenize(data, tokens, ".");
  if (tokens.size() < 2) {
    std::stringstream msg;
    msg << "Invalid file name format. " << data;
    NGTThrowException(msg);
  }
  auto &quantizer = index.getQuantizer();
  StaticObjectFileLoader loader(data, tokens[tokens.size() - 1]);
  size_t idx   = quantizer.objectList.size() == 0 ? 0 : quantizer.objectList.size() - 1;
  size_t count = 0;
  while (!loader.isEmpty()) {
    idx++;
    count++;
    if (dataSize > 0 && idx > dataSize) {
      break;
    }
    auto object = loader.getObject();
    if ((quantizer.property.distanceType == NGTQ::DistanceType::DistanceTypeInnerProduct) &&
        (object.size() + 1 == quantizer.objectList.genuineDimension)) {
      object.emplace_back(0);
    }
    if (object.size() != quantizer.objectList.genuineDimension) {
      std::stringstream msg;
      msg << "The dimension of the specified object is inconsistent with the dimension of the index. "
          << object.size() << ":" << quantizer.objectList.genuineDimension;
      NGTThrowException(msg);
    }
    index.insert(idx, object);
    if (count % 1000000 == 0) {
      std::cerr << "appended " << static_cast<float>(count) / 1000000.0 << "M objects.";
      if (count != idx) {
        std::cerr << " # of the total objects=" << static_cast<float>(idx) / 1000000.0 << "M";
      }
      cerr << " virtual memory(kbyte)=" << NGT::Common::getProcessVmSize() << std::endl;
    }
  }
  index.save();
  index.close();
  redirector.end();
}

void QBG::Index::preprocessingForNGT(std::string &indexPath, std::string &objectPath, bool verbose) {
  NGT::Property prop;
  {
    if (verbose) {
      std::cerr << "Opening the NGT index to get the property..." << std::endl;
    }
    NGT::Index index(indexPath);
    index.getProperty(prop);
  }
  if (verbose) {
    std::cerr << "object type:" << prop.objectType << std::endl;
    std::cerr << "dimension:" << prop.dimension << std::endl;
    std::cerr << "distance:" << prop.distanceType << std::endl;
  }
  if (prop.objectType != NGT::ObjectSpace::Qint4) {
    std::stringstream msg;
    msg << "Not qint4 object type. " << prop.objectType << std::endl;
    NGTThrowException(msg);
  }
  if (prop.distanceType != NGT::ObjectSpace::DistanceTypeL2) {
    std::stringstream msg;
    msg << "Not L2 distance type. " << prop.distanceType << std::endl;
    NGTThrowException(msg);
  }

  QBG::BuildParameters buildParameters;
  buildParameters.creation.localClusterDataType = NGTQ::ClusterDataTypePQ4;
  buildParameters.creation.genuineDimension     = prop.dimension;
  buildParameters.creation.dimension        = ((buildParameters.creation.genuineDimension + 15) / 16) * 16;
  buildParameters.creation.numOfSubvectors  = prop.dimension;
  buildParameters.creation.distanceType     = NGTQ::DistanceType::DistanceTypeL2;
  buildParameters.creation.genuineDataType  = ObjectFile::DataTypeFloat;
  buildParameters.creation.globalObjectType = NGT::ObjectSpace::ObjectType::Float;

  std::string qbgIndexPath = indexPath + "/" + NGT::Quantizer::getQbgIndex();
  if (verbose) {
    std::cerr << "qbg: creating.." << std::endl;
  }
  QBG::Index::create(qbgIndexPath, buildParameters, 0, objectPath);

  if (verbose) {
    std::cerr << "qbg: appending..." << std::endl;
  }
  size_t dataSize  = 0;
  std::string mode = "";
  QBG::Index::append(qbgIndexPath, objectPath, dataSize, verbose);

  QBG::Optimizer optimizer;

  optimizer.unifiedPQ     = true;
  optimizer.rotation      = false;
  optimizer.repositioning = false;
  optimizer.globalType    = QBG::Optimizer::GlobalTypeZero;

  if (verbose) {
    std::cerr << "qbg: optimizing..." << std::endl;
  }
  optimizer.optimize(qbgIndexPath);

  if (verbose) {
    std::cerr << "qbg: building..." << std::endl;
  }
  QBG::Index::buildNGTQ(qbgIndexPath, verbose);
}

#ifdef NGTQ_OBGRAPH

void QBG::Index::constructONNGForSetBlobID(size_t outEdge, size_t inEdge, size_t inRank,
                                           std::vector<NGT::ObjectDistances> &results,
                                           std::vector<std::vector<NGT::ObjectID>> &edges) {

  std::vector<std::vector<std::tuple<NGT::ObjectID, uint32_t, float>>> inEdges(results.size());
  if (inEdge > 0) {
    for (size_t id = 0; id < results.size(); id++) {
      for (size_t rank = 0; rank < results[id].size(); rank++) {
        if (inRank != 0 && inRank < rank) continue;
        inEdges[results[id][rank].id].emplace_back(std::make_tuple(id, rank, results[id][rank].distance));
      }
    }
  }

  edges.reserve(results.size());
  for (size_t id = 0; id < results.size(); id++) {
    std::vector<NGT::ObjectID> egs;
    std::vector<std::tuple<NGT::ObjectID, uint32_t, float>> iedges;
    if (inEdge > 0) {
      iedges = std::move(inEdges[id]);
      std::sort(iedges.begin(), iedges.end(),
                [](const std::tuple<NGT::ObjectID, uint32_t, float> &a,
                   const std::tuple<NGT::ObjectID, uint32_t, float> &b) -> bool {
                  if (std::get<1>(a) < std::get<1>(b)) {
                    return true;
                  } else if (std::get<1>(a) == std::get<1>(b)) {
                    return std::get<2>(a) < std::get<2>(b);
                  } else {
                    return false;
                  }
                });
    }
    auto ie = iedges.begin();
    for (size_t e = 0; e < results[id].size(); e++) {
      if (outEdge != 0 && e >= outEdge && ie == iedges.end()) break;
      if (outEdge == 0 || e < outEdge) {
        bool found = false;
        for (auto &ee : egs) {
          if (ee == results[id][e].id) {
            found = true;
            break;
          }
        }
        if (!found) {
          egs.emplace_back(results[id][e].id);
        }
      }
      while (ie != iedges.end() && std::get<1>((*ie)) == e) {
        if (inEdge > 0 && static_cast<size_t>(std::distance(iedges.begin(), ie)) >= inEdge) break;
        bool found = false;
        for (auto &ee : egs) {
          if (ee == std::get<0>((*ie))) {
            found = true;
            break;
          }
        }
        if (!found) {
          egs.emplace_back(std::get<0>(*ie));
        }
        ie++;
      }
    }
    for (; ie != iedges.end(); ++ie) {
      if (inEdge > 0 && static_cast<size_t>(std::distance(iedges.begin(), ie)) >= inEdge) break;
      bool found = false;
      for (auto &ee : egs) {
        if (ee == std::get<0>((*ie))) {
          found = true;
          break;
        }
      }
      if (found) continue;
      egs.emplace_back(std::get<0>((*ie)));
    }
    edges.emplace_back(std::move(egs));
    NGT::ObjectDistances().swap(results[id]);
    std::vector<std::tuple<NGT::ObjectID, uint32_t, float>>().swap(inEdges[id]);
  }
}

void QBG::Index::constructExperimentalGraphForSetBlobID(vector<NGT::ObjectDistances> &results, size_t outEdge,
                                                        float threshold,
                                                        std::vector<std::vector<NGT::ObjectID>> &edges) {
  size_t maxRank = 0;
  for (size_t srcID = 0; srcID < results.size(); srcID++) {
    std::cerr << "ID=" << srcID << " size=" << results[srcID].size() << std::endl;
    if (results[srcID].size() == 0) continue;
    if (maxRank == 0) {
      maxRank = results[srcID].size() - 1;
    } else if (maxRank != (results[srcID].size() - 1)) {
      std::stringstream msg;
      msg << "Different number of results. " << srcID << " " << maxRank << ":" << results[srcID].size()
          << std::endl;
      ;
      for (size_t r = 0; r < results[srcID].size(); r++) {
        msg << results[srcID][r].id << "," << results[srcID][r].distance << " ";
      }
      NGTThrowException(msg);
    }
    for (size_t srcRank = 0; srcRank < results[srcID].size(); srcRank++) {
      auto dstID = results[srcID][srcRank].id;
      if (dstID == 0) continue;
      auto distance = results[srcID][srcRank].distance;
      if (distance < 0.0) continue;
      size_t dstRank;
      for (dstRank = 0; dstRank < results[dstID].size(); dstRank++) {
        if (results[dstID][dstRank].id == srcID) {
          break;
        }
      }
      if (dstRank >= results[dstID].size()) {
        dstRank = maxRank + 1;
      } else {
      }
      float srcScore                   = static_cast<float>(srcRank) / (maxRank + 1);
      float dstScore                   = static_cast<float>(dstRank) / (maxRank + 1);
      results[srcID][srcRank].distance = -1;
      if (dstRank <= maxRank) {
        results[dstID][dstRank].distance = -1;
      }
      if (srcScore * dstScore > threshold) {
        results[srcID][srcRank].id = 0;
        if (dstRank <= maxRank) {
          results[dstID][dstRank].id = 0;
        }
      }
    }
  }

  for (size_t srcID = 0; srcID < results.size(); srcID++) {
    std::vector<NGT::ObjectID> egs;
    for (size_t srcRank = 0; srcRank < results[srcID].size(); srcRank++) {
      auto dstID = results[srcID][srcRank].id;
      if (dstID == 0) continue;
      auto distance = results[srcID][srcRank].distance;
      if (distance > 0.0) {
        std::cerr << "something wrong" << std::endl;
      }
      egs.emplace_back(dstID);
    }
    edges.emplace_back(std::move(egs));
    NGT::ObjectDistances().swap(results[srcID]);
  }
}

void QBG::Index::setBlobIDs(QBG::Index &index, std::vector<std::vector<NGT::ObjectID>> &edges) {
  auto &quantizedBlobGraph = index.quantizedBlobGraph;
  auto &quantizer          = index.getQuantizer();

  std::cerr << "object list size=" << quantizer.objectList.size();
  std::vector<uint32_t> oid2bid(quantizer.objectList.size());
  for (auto blobi = quantizedBlobGraph.begin(); blobi != quantizedBlobGraph.end(); ++blobi) {
    uint32_t bidx = distance(quantizedBlobGraph.begin(), blobi);
    if (bidx % 10000 == 0) {
      std::cerr << "# of extructed blobs is " << bidx << std::endl;
    }
    for (const auto id : (*blobi).ids) {
      oid2bid[id] = bidx;
    }
  }

  for (auto blobi = quantizedBlobGraph.begin(); blobi != quantizedBlobGraph.end(); ++blobi) {
    auto &blob    = *blobi;
    uint32_t bidx = distance(quantizedBlobGraph.begin(), blobi);
    if (bidx % 10000 == 0) {
      std::cerr << "# of processed blobs is " << bidx << std::endl;
    }
    if (blob.blobIDs.size() < blob.ids.size()) {
      blob.blobIDs.resize(blob.ids.size());
    }
    for (size_t idx = 0; idx < blob.ids.size(); idx++) {
      const auto id = blob.ids[idx];
      std::unordered_set<uint32_t> uniqueBlobIDs;

      for (auto &oid : edges[id]) {
        auto bi = oid2bid[oid];
        if (bi != bidx) {
          uniqueBlobIDs.insert(bi);
        }
      }

      std::vector<uint32_t> tempVector(uniqueBlobIDs.begin(), uniqueBlobIDs.end());
      blob.blobIDs[idx] = std::move(tempVector);
    }
  }
  quantizedBlobGraph.setGraphType(NGTQ::GraphTypeObjectBlobGraph);
  quantizer.property.graphType = NGTQ::GraphTypeObjectBlobGraph;
  index.saveProperty();
  index.save();
}

void QBG::Index::extractEdgesFromGraph(std::string &ngtPath, std::vector<std::vector<NGT::ObjectID>> &edges) {
  NGT::Index ngtIndex(ngtPath, false);
  NGT::GraphIndex &graph = static_cast<NGT::GraphIndex &>(ngtIndex.getIndex());
  edges.resize(graph.repository.size());
#pragma omp parallel for
  for (size_t id = 0; id < graph.repository.size(); id++) {
    if (id % 100000 == 0) {
      std::cerr << "# of processed objects:" << id << std::endl;
    }
    std::vector<NGT::ObjectID> egs;
    try {
      auto &node = *graph.getNode(id);
      egs.reserve(node.size());
      for (size_t rank = 0; rank < node.size(); rank++) {
        egs.emplace_back(node[rank].id);
      }
    } catch (NGT::Exception &err) {
      std::cerr << "Error: Cannot get the node for ID " << id << ". " << err.what() << std::endl;
      continue;
    }
    edges[id] = std::move(egs);
  }
}

void QBG::Index::extractEdgesBySearch(QBG::Index &index, std::vector<NGT::ObjectDistances> &results,
                                      size_t outEdge, size_t inEdge, QBG::SearchContainer &sc) {
  size_t size     = std::max(inEdge, outEdge);
  auto &quantizer = index.getQuantizer();
  results.resize(quantizer.objectList.size());
  const size_t batchSize = 1000;
  for (size_t batchID = 0; batchID < quantizer.objectList.size(); batchID += batchSize) {
    auto id = batchID;
    std::vector<std::vector<float>> objects;
    objects.reserve(batchSize);
    for (size_t idx = 0; id < quantizer.objectList.size() && idx < batchSize; id++, idx++) {
      if (id % 100000 == 0) {
        std::cerr << "# of processed objects:" << id << std::endl;
      }
      std::vector<float> object;
      if (id != 0) {
        quantizer.objectList.get(id, object);
      }
      objects.emplace_back(std::move(object));
    }
#pragma omp parallel for
    for (size_t idx = 0; idx < objects.size(); idx++) {
      if (objects[idx].size() == 0) continue;
      size_t oid = batchID + idx;
      QBG::SearchContainer searchContainer;
      NGT::ObjectDistances rObjects;
      searchContainer = sc;
      searchContainer.setObjectVector(objects[idx]);
      searchContainer.setResults(&rObjects);
      searchContainer.setSize(size + 1);
      searchContainer.setRefinementExpansion(sc.refinementExpansion);
      searchContainer.setNumOfProbes(sc.numOfProbes);
      searchContainer.setEpsilon(sc.explorationCoefficient - 1.0);
      searchContainer.setBlobEpsilon(sc.blobExplorationCoefficient - 1.0);
      searchContainer.setGraphExplorationSize(sc.graphExplorationSize);
      index.searchInTwoSteps(searchContainer);
      results[oid].reserve(rObjects.size() - 1);
      bool found = false;
      for (size_t i = 0; i < rObjects.size(); i++) {
        if (rObjects[i].id == oid) {
          found = true;
          continue;
        }
        results[oid].emplace_back(rObjects[i]);
      }
      if (results[oid].size() != size) {
        if (found) {
          std::cerr << "Something wrong? " << oid << std::endl;
        }
        while (results[oid].size() > size) {
          results[oid].pop_back();
        }
        if (results[oid].size() < size) {
          std::cerr << "Too few results! " << oid << std::endl;
        }
      }
    }
  }
}

#endif
