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
  setParameters(param.optimization);
#endif
}
QBG::Optimizer::Optimizer(QBG::OptimizationParameters &param) {
  setParameters(param);
}

void QBG::Optimizer::initialize() {
#ifdef NGTQ_QBG
  OptimizationParameters params;
  setParameters(params);
#endif
}

void QBG::Optimizer::setParameters(QBG::OptimizationParameters &param) {
#ifdef NGTQ_QBG
  clusteringType		= param.clusteringType;
  initMode			= param.initMode;

  timelimit			= param.timelimit;
  iteration			= param.iteration;
  clusterIteration		= param.clusterIteration;
  clusterSizeConstraint		= param.clusterSizeConstraint;
  clusterSizeConstraintCoefficient	= param.clusterSizeConstraintCoefficient;
  convergenceLimitTimes		= param.convergenceLimitTimes;
  numberOfObjects		= param.numOfObjects;
  numberOfClusters		= param.numOfClusters;
  numberOfSubvectors		= param.numOfSubvectors;
  numberOfMatrices			= param.numOfMatrices;
  seedNumberOfSteps		= param.seedNumberOfSteps;
  seedStep			= param.seedStep;
  reject			= param.reject;
  repositioning			= param.repositioning;
  rotation			= param.rotation;
  globalType			= param.globalType;
  randomizedObjectExtraction	= param.randomizedObjectExtraction;
  showClusterInfo		= param.showClusterInfo;
  verbose			= param.verbose;
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
  NGT::StdOstreamRedirector redirector(!verbose);
  redirector.begin();
  if (threadSize == 0) {
    threadSize = omp_get_max_threads();
  }
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
      std::stringstream msg;
      msg << "optimize: # of clusters is already specified. " << index.getQuantizer().property.localCentroidLimit << ":" << numberOfClusters;
      NGTThrowException(msg);

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
      size_t count = 0;
      {
	ifstream ifs(QBG::Index::getCodebookIndexFile(indexPath));
	if (!ifs) {
	  count = 1;
	  std::cerr << "the codebook index file is missing. this index must be QG." << std::endl;
	} else {
	  size_t id;
	  while (ifs >> id) {
	    count++;
	  }
	}
      }
      ofstream ofs(QBG::Index::getCodebookIndexFile(indexPath));
      if (!ofs) {
	std::stringstream msg;
	msg << "Cannot open the file. " << QBG::Index::getCodebookIndexFile(indexPath);
	NGTThrowException(msg);
      }
      for (size_t i = 0; i < count; i++) {
	ofs << "0" << std::endl;
      }
    } else if (globalType == GlobalTypeMean) {
      std::vector<std::vector<float>> vectors;
      std::string objects = QBG::Index::getTrainObjectFile(indexPath);
#ifdef NGT_CLUSTERING
      NGT::Clustering::loadVectors(objects, vectors);
#else
      loadVectors(objects, vectors);
#endif
      if (vectors.size() == 0 || vectors[0].size() == 0) {
	NGTThrowException("optimize::optimize: invalid input vectors");
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
void QBG::Optimizer::optimize(vector<vector<float>> &vectors, vector<vector<float>> &globalCentroid, Matrix<float> &r, vector<vector<NGT::Clustering::Cluster>> &localClusters, vector<double> &errors) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
  std::cerr << "optimize: Not implemented." << std::endl;
  abort();
#else

  if (vectors.size() == 0) {
    std::stringstream msg;
    msg << "optimize: error! the specified vetor is empty. the size=" << vectors.size();
    NGTThrowException(msg);
  }

  auto dim = vectors[0].size();
  if (numberOfSubvectors == 0) {
    std::stringstream msg;
    msg << "# of subspaces (m) is zero.";
    NGTThrowException(msg);
  }
  subvectorSize = dim / numberOfSubvectors;
  if (dim % numberOfSubvectors != 0) {
    std::stringstream msg;
    msg << "# of subspaces (m) is illegal. " << dim << ":" << numberOfSubvectors;
    NGTThrowException(msg);
  }

  generateResidualObjects(globalCentroid, vectors);


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
    std::cerr << "optimize: Each axis was repositioned." << std::endl;
  }


  vector<vector<vector<NGT::Clustering::Cluster>>> localClustersSet;

  numberOfMatrices = numberOfMatrices == 0 ? 1 : numberOfMatrices;
  if (!rotation) {
    iteration = 1;
    seedNumberOfSteps = 1;
    seedStep = 2;
    numberOfMatrices = 1;
  }

  seedNumberOfSteps = numberOfMatrices == 1 ? 1 : seedNumberOfSteps;
  seedNumberOfSteps = seedNumberOfSteps < 1 ? 1 : seedNumberOfSteps;
  if (numberOfMatrices > 100) {
    std::stringstream msg;
    msg << "# of matrices is too large. Should be less than 10. " << numberOfMatrices << " was specified.";
    NGTThrowException(msg);
  }
  if (seedStep > 100 || seedStep == 0) {
    std::stringstream msg;
    msg << "the seed step is illegal. Should be less than 100. " << seedStep << " was specified.";
    NGTThrowException(msg);
  }
  vector<Matrix<float>> rs(numberOfMatrices);
  for (auto &r: rs) {
    if (!rotation) {
      r.eye(dim);
    } else {
      r.randomRotation(dim);
    }
  }

  NGT::Timer timer;
  timer.start();
  for (int step = seedNumberOfSteps - 1; ; step--) {
    size_t vsize = vectors.size() / pow(seedStep, step);
    std::cerr << "optimize: # of vectors=" << vsize << "/" << vectors.size()
	      << ", # of matrices=" << numberOfMatrices << std::endl;
    if (vsize <= 1) {
      std::stringstream msg;
      msg << "# of partial vectors is too small, because # of vectors is too small or seedStep is too large."
	  << " # of partial vectors=" << vsize << " # of vectors=" <<  vectors.size() << " seedStep=" << seedStep << std::endl;
      NGTThrowException(msg);
    }
    auto partialVectors = vectors;
    if (vsize < vectors.size()) {
      partialVectors.resize(vsize);
    }

    optimize(partialVectors,
	     reposition,
	     rs,
	     localClustersSet,
	     errors);
    if (rs.size() > 1) {
      numberOfMatrices = static_cast<float>(numberOfMatrices) * (1.0 - reject);
      numberOfMatrices = numberOfMatrices == 0 ? 1 : numberOfMatrices;
      vector<pair<double, pair<Matrix<float>*, vector<vector<NGT::Clustering::Cluster>>*>>> sortedErrors;
      for (size_t idx = 0; idx < errors.size(); idx++) {
	sortedErrors.emplace_back(make_pair(errors[idx], make_pair(&rs[idx], &localClustersSet[idx])));
      }
      sort(sortedErrors.begin(), sortedErrors.end());
      vector<Matrix<float>> tmpMatrix;
      vector<vector<vector<NGT::Clustering::Cluster>>> tmpLocalClusters;
      for (size_t idx = 0; idx < numberOfMatrices; idx++) {
	tmpMatrix.emplace_back(*sortedErrors[idx].second.first);
	tmpLocalClusters.emplace_back(*sortedErrors[idx].second.second);
      }
      if (tmpMatrix.size() != numberOfMatrices) {
	std::cerr << "something strange. " << tmpMatrix.size() << ":" << numberOfMatrices << std::endl;
      }
      rs = std::move(tmpMatrix);
      localClustersSet = std::move(tmpLocalClusters);
    }
    timer.stop();
    std::cerr << "optimize: time=" << timer << std::endl;
    timer.restart();
    if (vsize >= vectors.size()) {
      break;
    }
  }

  if (rs.size() != 1) {
    std::cerr << "optimize: Warning. rs.size=" << rs.size() << std::endl;
  }

  localClusters = std::move(localClustersSet[0]);
  if (repositioning) {
    Matrix<float> repositionedR(reposition);
    repositionedR.mul(rs[0]);
    r = std::move(repositionedR);
  } else {
    r = std::move(rs[0]);
  }

  //-/size_t pos = std::distance(std::find(ofile.rbegin(), ofile.rend(), '.'), ofile.rend()) - 1;

#endif
}

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

  vector<vector<float>> globalCentroid;
  NGT::Clustering::loadVectors(global, globalCentroid);

  Matrix<float> r;
  vector<vector<NGT::Clustering::Cluster>> localClusters;
  vector<double> errors;

  optimize(vectors, globalCentroid, r, localClusters, errors);

  Matrix<float>::save(ofile + QBG::Index::getRotationFile(), r);

  if (showClusterInfo) {
    if (localClusters.size() != numberOfSubvectors) {
      std::stringstream msg;
      msg << "Fatal error. localClusters.size() != numberOfSubvectors "
	  << localClusters.size() << ":" << numberOfSubvectors;
      NGTThrowException(msg);
    }
    float totalRate = 0.0;
    for (size_t m = 0; m < localClusters.size(); m++) {
      size_t min = std::numeric_limits<size_t>::max();
      size_t max = 0;
      size_t nOfVectors = 0;
      for (size_t i = 0; i < localClusters[m].size(); i++) {
	nOfVectors += localClusters[m][i].members.size();
	if (localClusters[m][i].members.size() < min) {
	  min = localClusters[m][i].members.size();
	}
	if (localClusters[m][i].members.size() > max) {
	  max = localClusters[m][i].members.size();
	}
      }
      float rate = static_cast<float>(max - min) / static_cast<float>(nOfVectors);
      totalRate += rate;
      std::cout << "cluster " << m << " " << rate << "," << max - min << "," << min << "," << max << " : ";
      for (size_t i = 0; i < localClusters[m].size(); i++) {
	std::cout << localClusters[m][i].members.size() << " ";
      }
      std::cout << std::endl;
    }
    totalRate /= localClusters.size();
    std::cout << "Range rate=" << totalRate << std::endl;
    std::cout << "Error=" << errors[0] << std::endl;
  }
  for (size_t m = 0; m < numberOfSubvectors; m++) {
    stringstream str;
    str << ofile << QBG::Index::getSubvectorPrefix() << "-" << m;
#ifdef NGT_CLUSTERING
    NGT::Clustering::saveClusters(str.str(), localClusters[m]);
#else
    saveClusters(str.str(), localClusters[m]);
#endif
  }
#endif
}
#endif
