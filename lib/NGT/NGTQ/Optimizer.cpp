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

#include "QuantizedBlobGraph.h"
#include "Optimizer.h"

QBG::Optimizer::Optimizer(QBG::BuildParameters &param) {
#ifdef NGTQ_QBG
  clusteringType		= param.optimization.clusteringType;
  initMode			= param.optimization.initMode;

  timelimit			= param.optimization.timelimit;
  iteration			= param.optimization.iteration;
  clusterIteration		= param.optimization.clusterIteration;
  clusterSizeConstraint		= param.optimization.clusterSizeConstraint;
  clusterSizeConstraintCoefficient	= param.optimization.clusterSizeConstraintCoefficient;
  convergenceLimitTimes		= param.optimization.convergenceLimitTimes;
  numberOfObjects		= param.optimization.numberOfObjects;
  numberOfClusters		= param.optimization.numberOfClusters;
  numberOfSubvectors		= param.optimization.numberOfSubvectors;
  nOfMatrices			= param.optimization.nOfMatrices;
  seedStartObjectSizeRate	= param.optimization.seedStartObjectSizeRate;
  seedStep			= param.optimization.seedStep;
  reject			= param.optimization.reject;
  repositioning			= param.optimization.repositioning;
  rotation			= param.optimization.rotation;
  globalType			= param.optimization.globalType;
  randomizedObjectExtraction	= param.optimization.randomizedObjectExtraction;
  showClusterInfo		= param.optimization.showClusterInfo;
  silence			= param.silence;
#endif
}


#ifdef NGTQ_QBG
void QBG::Optimizer::evaluate(string global, vector<vector<float>> &vectors, char clusteringType, string &ofile, size_t &numberOfSubvectors, size_t &subvectorSize)
{
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
  std::cerr << "evaluate: Not implemented." << std::endl;
  abort();
#else
  vector<vector<float>> residualVectors;
  {
    // compute residual vectors by global centroids.
#ifdef NGT_CLUSTERING
    vector<NGT::Clustering::Cluster> globalCentroid;
#else
    vector<Cluster> globalCentroid;
#endif
    if (!global.empty()) {
      cerr << "generate residual vectors." << endl;
      try {
#ifdef NGT_CLUSTERING
	NGT::Clustering::loadClusters(global, globalCentroid);
#else
	loadClusters(global, globalCentroid);
#endif
      } catch (...) {
	cerr << "Cannot load vectors. " << global << endl;
	return;
      }
      if (clusteringType == 'k') {
#ifdef NGT_CLUSTERING
	NGT::Clustering::assign(vectors, globalCentroid);
#else
	assign(vectors, globalCentroid);
#endif
      } else {
	cerr << "Using NGT" << endl;
#ifdef NGT_CLUSTERING
	std::cerr << "Not implemented" << std::endl;
	abort();
#else
	assignWithNGT(vectors, globalCentroid);
#endif
      }
      residualVectors.resize(vectors.size());
      cerr << "global centroid size=" << globalCentroid.size() << endl;
      for (size_t cidx = 0; cidx < globalCentroid.size(); ++cidx) {
	for (auto mit = globalCentroid[cidx].members.begin(); mit != globalCentroid[cidx].members.end(); ++mit) {
	  size_t vid = (*mit).vectorID;
	  residualVectors[vid] = vectors[vid];
#ifdef NGT_CLUSTERING
	  NGT::Clustering::subtract(residualVectors[vid], globalCentroid[cidx].centroid);
#else
	  subtract(residualVectors[vid], globalCentroid[cidx].centroid);
#endif
	}
      }
    }
  }

  Matrix<float> R;
  Matrix<float>::load(ofile + QBG::Index::getRotationFile(), R);    
  vector<vector<float>> qv(vectors.size());	// quantized vector
  vector<vector<float>> xp;	// residual vector
  if (residualVectors.empty()) {
    xp = vectors;
  } else {
    xp = residualVectors;
  }
  Matrix<float>::mulSquare(xp, R);
  for (size_t m = 0; m < numberOfSubvectors; m++) {
#ifdef NGT_CLUSTERING
    vector<NGT::Clustering::Cluster> subClusters;
#else
    vector<Cluster> subClusters;
#endif
    stringstream str;
    str << ofile << "-" << m;
#ifdef NGT_CLUSTERING
    NGT::Clustering::loadClusters(str.str(), subClusters);
#else
    loadClusters(str.str(), subClusters);
#endif
    vector<vector<float>> subVectors;
    extractSubvector(xp, subVectors, m * subvectorSize, subvectorSize);
    if (clusteringType == 'k') {
#ifdef NGT_CLUSTERING
      NGT::Clustering::assign(subVectors, subClusters);
#else
      assign(subVectors, subClusters);
#endif
    } else {
      cerr << "Using NGT for subvector" << endl;
#ifdef NGT_CLUSTERING
      std::cerr << "not implemented" << std::endl;
      abort();
#else
      assignWithNGT(subVectors, subClusters);
#endif
    }
#ifdef NGT_CLUSTERING
    double distortion = NGT::Clustering::calculateML2(subVectors, subClusters);
#else
    double distortion = calculateML2(subVectors, subClusters);
#endif
    cout << "distortion[" << m << "]=" << distortion << endl;
    vector<vector<float>> subCentroids(vectors.size());
    for (size_t cidx = 0; cidx < subClusters.size(); ++cidx) {
#ifdef NGT_CLUSTERING
      vector<NGT::Clustering::Entry> &members = subClusters[cidx].members;
#else
      vector<Entry> &members = subClusters[cidx].members;
#endif
      for (size_t eidx = 0; eidx < members.size(); ++eidx) {
#ifdef NGT_CLUSTERING
	NGT::Clustering::Entry &entry = members[eidx];
#else
	Entry &entry = members[eidx];
#endif
	assert(cidx == entry.centroidID);
	subCentroids[entry.vectorID] = subClusters[cidx].centroid;
      }
    }
    catSubvector(qv, subCentroids);
  }
#ifdef NGT_CLUSTERING
  double distortion = NGT::Clustering::distanceL2(qv, xp);
#else
  double distortion = distanceL2(qv, xp);
#endif
  cout << "distortion=" << distortion << endl;
  return;
#endif 
}
#endif

#ifdef NGTQ_QBG
void QBG::Optimizer::evaluate(vector<vector<float>> &vectors, string &ofile, size_t &numberOfSubvectors, size_t &subvectorSize) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
  std::cerr << "evaluate: Not implemented." << std::endl;
  abort();
#else
  cerr << "Evaluate" << endl;
  Matrix<float> R;
  Matrix<float>::load(ofile + QBG::Index::getRotationFile(), R);    
  vector<vector<float>> xp = vectors;
  Matrix<float>::mulSquare(xp, R);
  for (size_t m = 0; m < numberOfSubvectors; m++) {
#ifdef NGT_CLUSTERING
    vector<NGT::Clustering::Cluster> subClusters;
#else
    vector<Cluster> subClusters;
#endif
    stringstream str;
    str << ofile << "-" << m;
#ifdef NGT_CLUSTERING
    NGT::Clustering::loadClusters(str.str(), subClusters);
#else
    loadClusters(str.str(), subClusters);
#endif
    vector<vector<float>> subVectors;
    extractSubvector(xp, subVectors, m * subvectorSize, subvectorSize);
#ifdef NGT_CLUSTERING
    NGT::Clustering::assign(subVectors, subClusters);
    double distortion = NGT::Clustering::calculateML2(subVectors, subClusters);
#else
    assign(subVectors, subClusters);
    double distortion = calculateML2(subVectors, subClusters);
#endif
    cout << "distortion[" << m << "]=" << distortion << endl;
    for (size_t cidx = 0; cidx < subClusters.size(); ++cidx) {
      cout << "  members[" << cidx << "]=" << subClusters[cidx].members.size() << endl;
    }
  }
  return;
#endif 
}
#endif

#ifdef NGTQ_QBG
void QBG::Optimizer::optimize(const std::string indexPath, size_t threadSize) {
  NGT::StdOstreamRedirector redirector(silence);
  redirector.begin();
  try {
    QBG::Index index(indexPath);
    if (index.getQuantizer().objectList.size() <= 1) {
      NGTThrowException("optimize: No objects");
    }
    std::cerr << "optimize: # of objects=" << numberOfObjects << std::endl;
    if (numberOfObjects == 0) {
      numberOfObjects = 1000;
      numberOfObjects = index.getQuantizer().objectList.size() - 1 < numberOfObjects ? index.getQuantizer().objectList.size() - 1 : numberOfObjects;
    }
    std::cerr << "optimize: updated # of objects=" << numberOfObjects << std::endl;
    std::cerr << "optimize: # of clusters=" << index.getQuantizer().property.localCentroidLimit << ":" << numberOfClusters << std::endl;
    if (index.getQuantizer().property.localCentroidLimit == 0 && numberOfClusters == 0) {
      std::stringstream msg;
      msg << "optimize: # of clusters is illegal. " << index.getQuantizer().property.localCentroidLimit << ":" << numberOfClusters;
      NGTThrowException(msg);
    }
    if (index.getQuantizer().property.localCentroidLimit != 0 && numberOfClusters != 0 && 
	index.getQuantizer().property.localCentroidLimit != numberOfClusters) {
      std::cerr << "optimize: warning! # of clusters is already specified. " << index.getQuantizer().property.localCentroidLimit << ":" << numberOfClusters << std::endl;
    }
    if (numberOfClusters == 0) {
      numberOfClusters = index.getQuantizer().property.localCentroidLimit;
    }

    if (numberOfSubvectors == 0 && index.getQuantizer().property.localDivisionNo == 0) {
      std::stringstream msg;
      msg << "optimize: # of subvectors is illegal. " << numberOfSubvectors << ":" << index.getQuantizer().property.localDivisionNo;
      NGTThrowException(msg);
    }
    if (numberOfSubvectors != 0 && index.getQuantizer().property.localDivisionNo != 0 &&
	numberOfSubvectors != index.getQuantizer().property.localDivisionNo) {
      std::cerr << "optimize: warning! # of subvectros is already specified. " << numberOfSubvectors << ":" << index.getQuantizer().property.localDivisionNo << std::endl;
    }
    if (numberOfSubvectors == 0) {
      numberOfSubvectors = index.getQuantizer().property.localDivisionNo;
    }

    const std::string ws = indexPath + "/" + QBG::Index::getWorkspaceName();
    try {
      NGT::Index::mkdir(ws);
    } catch(...) {}
    const std::string object = QBG::Index::getTrainObjectFile(indexPath);
    std::ofstream ofs;
    ofs.open(object);
    index.extract(ofs, numberOfObjects, randomizedObjectExtraction);
    if (globalType == GlobalTypeZero) {
      assert(index.getQuantizer().objectList.pseudoDimension != 0);
      std::vector<std::vector<float>> global(1);
      global[0].resize(index.getQuantizer().property.dimension, 0.0);
      NGT::Clustering::saveVectors(QBG::Index::getQuantizerCodebookFile(indexPath), global);
    } else if (globalType == GlobalTypeMean) {
      std::vector<std::vector<float>> vectors;
      std::string objects = QBG::Index::getTrainObjectFile(indexPath);
#ifdef NGT_CLUSTERING
      NGT::Clustering::loadVectors(objects, vectors);
#else
      loadVectors(objects, vectors);
#endif
      if (vectors.size() == 0 || vectors[0].size() == 0) {
	NGTThrowException("Optimizer::optimize: invalid input vectors");
      }
      std::vector<std::vector<float>> global(1);
      global[0].resize(index.getQuantizer().property.dimension, 0);
      for (auto v = vectors.begin(); v != vectors.end(); ++v) {
	for (size_t i = 0; i < (*v).size(); i++) {
	  global[0][i] += (*v)[i];
	}
      }
      for (size_t i = 0; i < global[0].size(); i++) {
	global[0][i] /= vectors.size();
      }
      NGT::Clustering::saveVectors(QBG::Index::getQuantizerCodebookFile(indexPath), global);
    }

    optimizeWithinIndex(indexPath);

  } catch(NGT::Exception &err) {
    redirector.end();
    throw err;
  }

  redirector.end();
}
#endif 

#ifdef NGTQ_QBG
void QBG::Optimizer::optimizeWithinIndex(std::string indexPath) {
  std::string object;
  std::string pq;
  std::string global;
  {
    object = QBG::Index::getTrainObjectFile(indexPath);
    pq = QBG::Index::getPQFile(indexPath);
    global = QBG::Index::getQuantizerCodebookFile(indexPath);
  }

  try {
    NGT::Index::mkdir(pq);
  } catch(...) {}
  pq += "/";
  optimize(object, pq, global);
}
#endif 



#ifdef NGTQ_QBG
void QBG::Optimizer::optimize(std::string invector, std::string ofile, std::string global) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
  std::cerr << "optimize: Not implemented." << std::endl;
  abort();
#else
  vector<vector<float>> vectors;


#ifdef NGT_CLUSTERING
  NGT::Clustering::loadVectors(invector, vectors);
#else
  loadVectors(invector, vectors);
#endif

  if (vectors.size() == 0) {
    std::stringstream msg;
    msg << "Optimizer: error! the specified vetor file is empty. " << invector << ". the size=" << vectors.size();
    NGTThrowException(msg);
  }

  dim = vectors[0].size();
  subvectorSize = dim / numberOfSubvectors;
  if (dim % numberOfSubvectors != 0) {
    std::stringstream msg;
    msg << "# of subspaces (m) is illegal. " << dim << ":" << numberOfSubvectors;
    NGTThrowException(msg);
  }


  timelimitTimer.start();

  Matrix<float> reposition;
  if (repositioning) {
    reposition.zero(dim, dim);
    size_t dstidx = 0;
    for (size_t didx = 0; didx < numberOfSubvectors; didx++) {
      for (size_t sdidx = 0; sdidx < subvectorSize; sdidx++) {
	size_t srcidx = didx + sdidx * numberOfSubvectors;
	auto col = dstidx;
	auto row = srcidx;
	reposition.set(row, col, 1.0);
	dstidx++;
      }
    }
    std::cerr << "Optimizer: Each axis was repositioned." << std::endl;
  }


  vector<vector<vector<NGT::Clustering::Cluster>>> localClusters;
  vector<double> errors;

  bool useEye = false;
  nOfMatrices = nOfMatrices == 0 ? 1 : nOfMatrices;
  if (!rotation) {
    iteration = 1;
    seedStartObjectSizeRate = 1.0;
    seedStep = 2;
  }
  useEye = !rotation;
  vector<Matrix<float>> rs(nOfMatrices);
  for (auto &r: rs) {
    if (useEye) {
      r.eye(dim);
    } else {
      r.randomRotation(dim);
    }
  }

  for (size_t vsize = static_cast<float>(vectors.size()) * seedStartObjectSizeRate; ; vsize *= seedStep) {
    auto partialVectors = vectors;
    if (vsize < vectors.size()) {
      partialVectors.resize(vsize);
    }

    optimize(partialVectors,
	     global,
	     ofile,
	     reposition,
	     rs,
	     localClusters,
	     errors);
    if (rs.size() > 1) {
      nOfMatrices = static_cast<float>(nOfMatrices) * (1.0 - reject);
      nOfMatrices = nOfMatrices == 0 ? 1 : nOfMatrices;
      vector<pair<double, pair<Matrix<float>*, vector<vector<NGT::Clustering::Cluster>>*>>> sortedErrors;
      for (size_t idx = 0; idx < errors.size(); idx++) {
	sortedErrors.emplace_back(make_pair(errors[idx], make_pair(&rs[idx], &localClusters[idx])));
      }
      sort(sortedErrors.begin(), sortedErrors.end());
      vector<Matrix<float>> tmpMatrix;
      vector<vector<vector<NGT::Clustering::Cluster>>> tmpLocalClusters;
      for (size_t idx = 0; idx < nOfMatrices; idx++) {
	tmpMatrix.emplace_back(*sortedErrors[idx].second.first);
	tmpLocalClusters.emplace_back(*sortedErrors[idx].second.second);
      }
      if (tmpMatrix.size() != nOfMatrices) {
	std::cerr << "something strange. " << tmpMatrix.size() << ":" << nOfMatrices << std::endl;
      }
      rs = std::move(tmpMatrix);
      localClusters = std::move(tmpLocalClusters);
    }
    if (vsize >= vectors.size()) {
      break;
    }
  }

  if (rs.size() != 1) {
    std::cerr << "Optimizer: Warning. rs.size=" << rs.size() << std::endl;
  }
  auto minR = std::move(rs[0]);
  auto minLocalClusters = std::move(localClusters[0]);
  //-/size_t pos = std::distance(std::find(ofile.rbegin(), ofile.rend(), '.'), ofile.rend()) - 1;
  if (repositioning) {
    Matrix<float> repositionedR(reposition);
    repositionedR.mul(minR);
    Matrix<float>::save(ofile + QBG::Index::getRotationFile(), repositionedR);
  } else {
    Matrix<float>::save(ofile + QBG::Index::getRotationFile(), minR);
  }
  if (showClusterInfo) {
    if (minLocalClusters.size() != numberOfSubvectors) {
      std::stringstream msg;
      msg << "Fatal error. minLocalClusters.size() != numberOfSubvectors " <<
	minLocalClusters.size() << ":" << numberOfSubvectors;
      NGTThrowException(msg);
    }
    float totalRate = 0.0;
    for (size_t m = 0; m < minLocalClusters.size(); m++) {
      size_t min = std::numeric_limits<size_t>::max();
      size_t max = 0;
      size_t nOfVectors = 0;
      for (size_t i = 0; i < minLocalClusters[m].size(); i++) {
	nOfVectors += minLocalClusters[m][i].members.size();
	if (minLocalClusters[m][i].members.size() < min) {
	  min = minLocalClusters[m][i].members.size();
	}
	if (minLocalClusters[m][i].members.size() > max) {
	  max = minLocalClusters[m][i].members.size();
	}
      }
      float rate = static_cast<float>(max - min) / static_cast<float>(nOfVectors);
      totalRate += rate;
      std::cout << "cluster " << m << " " << rate << "," << max - min << "," << min << "," << max << " : ";
      for (size_t i = 0; i < minLocalClusters[m].size(); i++) {
	std::cout << minLocalClusters[m][i].members.size() << " ";
      }
      std::cout << std::endl;
    }
    totalRate /= minLocalClusters.size();
    std::cout << "Range rate=" << totalRate << std::endl;
    std::cout << "Error=" << errors[0] << std::endl;
  }
  for (size_t m = 0; m < numberOfSubvectors; m++) {
    stringstream str;
    str << ofile << QBG::Index::getSubvectorPrefix() << "-" << m;
#ifdef NGT_CLUSTERING
    NGT::Clustering::saveClusters(str.str(), minLocalClusters[m]);
#else
    saveClusters(str.str(), minLocalClusters[m]);
#endif
  }
#endif 
}
#endif
