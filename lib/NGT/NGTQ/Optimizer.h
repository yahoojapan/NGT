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

#define NGT_CLUSTERING

#include "QuantizedBlobGraph.h"
#ifdef NGT_CLUSTERING
#include "NGT/Clustering.h"
#else
#include "Cluster.h"
#endif

#include "Matrix.h"

namespace NGTQ {
  class Optimizer {
  public:
    enum GlobalType {
      GlobalTypeNone = 0,
      GlobalTypeZero = 1,
      GlobalTypeMean = 2
    };

    Optimizer() {
      numberOfClusters = 0;
      numberOfSubvectors = 0;
      clusteringType = NGT::Clustering::ClusteringTypeKmeansWithNGT;
      initMode = NGT::Clustering::InitializationModeRandom;
      convergenceLimitTimes = 5;
      iteration = 100;
      clusterIteration = 100;
      clusterSizeConstraint = false;
      iteration = 100;
      repositioning = false;
      rotation = true;
      globalType = GlobalTypeNone;
      silence = false;
    }

    static void
      extractSubvector(vector<vector<float>> &vectors, vector<vector<float>> &subvectors, size_t start , size_t size)
    {
      size_t vsize = vectors.size();
      subvectors.clear();
      subvectors.resize(vsize);
      for (size_t vidx = 0; vidx < vsize; vidx++) {
	subvectors[vidx].reserve(size);
	for (size_t i = 0; i < size; i++) {
	  subvectors[vidx].push_back(vectors[vidx][start + i]);
	}
      }
    }

    static void
      catSubvector(vector<vector<float>> &vectors, vector<vector<float>> &subvectors)
    {
      if (vectors.size() == 0) {
	vectors.resize(subvectors.size());
      }
      assert(vectors.size() == subvectors.size());
      size_t vsize = vectors.size();
      size_t subvsize = subvectors[0].size();
      size_t size = vectors[0].size() + subvsize;
      for (size_t vidx = 0; vidx < vsize; vidx++) {
	subvectors[vidx].reserve(size);
	for (size_t i = 0; i < subvsize; i++) {
	  vectors[vidx].push_back(subvectors[vidx][i]);
	}
      }
    }

    static void
#ifdef NGT_CLUSTERING
      extractQuantizedVector(vector<vector<float>> &qvectors, vector<NGT::Clustering::Cluster> &clusters) 
#else
      extractQuantizedVector(vector<vector<float>> &qvectors, vector<Cluster> &clusters) 
#endif
    {
      for (size_t cidx = 0; cidx < clusters.size(); ++cidx) {
	for (auto mit = clusters[cidx].members.begin(); mit != clusters[cidx].members.end(); ++mit) {
	  size_t vid = (*mit).vectorID;
	  if (vid >= qvectors.size()) {
	    qvectors.resize(vid + 1);
	  }
	  qvectors[vid] = clusters[cidx].centroid;
	}
      }
    }

    static double
      squareDistance(vector<vector<float>> &va, vector<vector<float>> &vb)
    {
      assert(va.size() == vb.size());
      size_t vsize = va.size();
      assert(va[0].size() == vb[0].size());
      size_t dim = va[0].size();
      double distance = 0;
      for (size_t vidx = 0; vidx < vsize; vidx++) {
	for (size_t i = 0; i < dim; i++) {
	  double d = va[vidx][i] - vb[vidx][i];
	  distance += d * d;
	}
      }
      distance /= dim * vsize;
      return distance;
    }

    void evaluate(string global, vector<vector<float>> &vectors, char clusteringType, string &ofile, size_t &numberOfSubvectors, size_t &subvectorSize)
    {
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
      Matrix<float>::load(ofile + "_R.tsv", R);    
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
	str << ofile << "-" << m << ".tsv";
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
    }

    void evaluate(vector<vector<float>> &vectors, string &ofile, size_t &numberOfSubvectors, size_t &subvectorSize) {
      cerr << "Evaluate" << endl;
      Matrix<float> R;
      Matrix<float>::load(ofile + "_R.tsv", R);    
      vector<vector<float>> xp = vectors;
      Matrix<float>::mulSquare(xp, R);
      for (size_t m = 0; m < numberOfSubvectors; m++) {
#ifdef NGT_CLUSTERING
	vector<NGT::Clustering::Cluster> subClusters;
#else
	vector<Cluster> subClusters;
#endif
	stringstream str;
	str << ofile << "-" << m << ".tsv";
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
    }

    void
      generateResidualObjects(string global, vector<vector<float>> &vectors)
    {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      std::cerr << "generateResidualObjects: Not implemented." << std::endl;
      abort();
#else
      if (global.empty()) {
	NGTThrowException("A global codebook is not specified!");
      }
      vector<vector<float>> residualVectors;
      vector<vector<float>> globalCentroid;
      try {
	NGT::Clustering::loadVectors(global, globalCentroid);
      } catch (...) {
	std::stringstream msg;
	msg << "Optimizer::generateResidualObjects: Cannot load global vectors. " << global;
	NGTThrowException(msg);
      }

      NGT::Property property;
      property.objectType = NGT::Index::Property::ObjectType::Float;
      property.distanceType = NGT::Index::Property::DistanceType::DistanceTypeL2;
      if (globalCentroid[0].size() != vectors[0].size()) {
	std::cerr << "optimizer: Warning. The dimension is inconsistency. " << globalCentroid[0].size() << ":" << vectors[0].size() << std::endl;
      }
      property.dimension = vectors[0].size();
      NGT::Index index(property);
      for (auto &c : globalCentroid) {
	if (static_cast<int>(c.size()) != property.dimension) {
	  c.resize(property.dimension);
	}
	index.append(c);
      }
      index.createIndex(100);

      for (size_t idx = 0; idx < vectors.size(); idx++) {
	auto &v = vectors[idx];
	if ((idx + 1) % 10000 == 0) {
	}
	NGT::ObjectDistances gc;
	NGT::SearchQuery query(v);
	query.setResults(&gc);
	query.setSize(10);
	query.setEpsilon(0.1);
	index.search(query);
	if (gc.empty()) {
	  std::cerr << "inner fatal error. no results! something wrong." << std::endl;
	  abort();
	}
	if (gc[0].id == 0 || gc[0].id > globalCentroid.size()) {
	  std::cerr << "wrong id " << gc[0].id << ":" <<  globalCentroid.size() << std::endl;
	  abort();
	}
	auto gcidx = gc[0].id - 1;
	try {
	  NGT::Clustering::subtract(v, globalCentroid[gcidx]);
	} catch (NGT::Exception &err) {
	  std::cerr << err.what() << ":" << v.size() << "x" << globalCentroid[gcidx].size() << std::endl;
	  abort();
	}
      }
#endif 
    }

    static void optimizeRotation(
				 size_t iteration,
				 vector<vector<float>> &vectors,
				 Matrix<float> &xt,
				 Matrix<float> &R,
				 Matrix<float> &minR,
				 vector<vector<NGT::Clustering::Cluster>> &minLocalClusters,
				 NGT::Clustering::ClusteringType clusteringType,
				 NGT::Clustering::InitializationMode initMode,
				 size_t numberOfClusters,
				 size_t numberOfSubvectors,
				 size_t subvectorSize,
				 size_t clusterIteration,
				 bool clusterSizeConstraint,
				 size_t convergenceLimitTimes,
				 double &minDistortion,
				 NGT::Timer &timelimitTimer, float timelimit,
				 bool rotation
				 ) {

      minDistortion = DBL_MAX;

      int minIt = 0;
      for (size_t it = 0; it < iteration; it++) {
	vector<vector<float>> xp = vectors;
	Matrix<float>::mulSquare(xp, R);
	float distance = 0.0;
#ifdef NGT_CLUSTERING
	vector<vector<NGT::Clustering::Cluster>> localClusters(numberOfSubvectors);
#else
	vector<vector<Cluster>> localClusters(numberOfSubvectors);
#endif
	vector<vector<float>> subQuantizedVectors[numberOfSubvectors];
#define ERROR_CALCULATION  
#ifdef ERROR_CALCULATION
	vector<float> subvectorDistances(numberOfSubvectors);
#endif
#pragma omp parallel for
	for (size_t m = 0; m < numberOfSubvectors; m++) {
	  vector<vector<float>> subVectors;
	  extractSubvector(xp, subVectors, m * subvectorSize, subvectorSize);
#ifdef NGT_CLUSTERING
	  vector<NGT::Clustering::Cluster> &clusters = localClusters[m];
#else
	  vector<Cluster> &clusters = localClusters[m];
#endif
#ifdef NGT_CLUSTERING
	  NGT::Clustering clustering(initMode, clusteringType, clusterIteration);
	  //clustering.setClusterSizeConstraintCoefficient(1.2);
	  clustering.clusterSizeConstraint = clusterSizeConstraint;
	  clustering.kmeans(subVectors, numberOfClusters, clusters);
#else
	  size_t reassign;
	  if (clusteringType == 'k') {
	    reassign = kmeansClustering(initMode, subVectors, numberOfClusters, clusters, clusterSizeConstraint, 0);
	  } else {
	    reassign = kmeansClusteringWithNGT(initMode, subVectors, numberOfClusters, clusters, clusterSizeConstraint, 0);    
	  }
#endif
	  extractQuantizedVector(subQuantizedVectors[m], clusters);
	  assert(subQuantizedVectors[m].size() == vectors.size());
	  assert(subQuantizedVectors[m][0].size() == subvectorSize);
	  // 入力部分ベクトルと量子化部分ベクトルを比較して量子化誤差の計算
#ifdef ERROR_CALCULATION
	  double d = squareDistance(subQuantizedVectors[m], subVectors);
	  subvectorDistances[m] = d;
#endif
	}
	distance = 0.0;
#ifdef ERROR_CALCULATION
	for (size_t m = 0; m < numberOfSubvectors; m++) {
	  distance += subvectorDistances[m];
	}
#endif
	vector<vector<float>> quantizedVectors;
	for (size_t m = 0; m < numberOfSubvectors; m++) {
	  catSubvector(quantizedVectors, subQuantizedVectors[m]);
	}
	distance = sqrt(distance / numberOfSubvectors);
	if (minDistortion > distance) {
	  minDistortion = distance;
	  minR = R;
	  minIt = it;
	  minLocalClusters = localClusters;
	}
	if (it + 1 > iteration || it - minIt > convergenceLimitTimes) {
	  break;
	}
	timelimitTimer.stop();
	if (timelimitTimer.time > timelimit) {
	  std::cerr << "Optimizer: Warning. The elapsed time exceeded the limit-time. " << timelimit << ":" << timelimitTimer.time << std::endl;
	  timelimitTimer.restart();
	  break;
	}
	timelimitTimer.restart();
	if (rotation) {
	  Matrix<float> a(xt);
	  a.mul(quantizedVectors);
	  Matrix<float> u, s, v;
	  Matrix<float>::svd(a, u, s, v);
	  v.transpose();
	  u.mul(v);
	  R = u;
	}
      }
    }

    void optimize(vector<vector<float>> &vectors,
		  string global,
		  string ofile,
		  Matrix<float> &reposition,
		  vector<Matrix<float>> &rs,
		  vector<vector<vector<NGT::Clustering::Cluster>>> &localClusters,
		  vector<double> &errors
		  ) {
      if (vectors.size() == 0) {
	NGTThrowException("the vector is empty");
      }
      generateResidualObjects(global, vectors);
      if (!reposition.isEmpty()) {
	Matrix<float>::mulSquare(vectors, reposition);
      }
      Matrix<float> xt(vectors);
#ifndef NGTQ_BLAS_MATRIX
      xt.transpose();
#endif
      localClusters.resize(rs.size());
      errors.resize(rs.size());
      NGT::Timer timer;
      for (size_t ri = 0; ri < rs.size(); ri++) {
	auto imode = initMode;
	if (imode == NGT::Clustering::InitializationModeBest) {
	  imode = ri % 2 == 0 ? NGT::Clustering::InitializationModeRandom : NGT::Clustering::InitializationModeKmeansPlusPlus;
	}
	timer.start();
	Matrix<float> optr;
	optimizeRotation(
			 iteration,		
			 vectors,		
			 xt,			
			 rs[ri],		
			 optr,			
			 localClusters[ri],
			 clusteringType,
			 imode,
			 numberOfClusters,
			 numberOfSubvectors,	
			 subvectorSize,		
			 clusterIteration,	
			 clusterSizeConstraint,
			 convergenceLimitTimes,
			 errors[ri],
			 timelimitTimer, timelimit,
			 rotation
			 );
	timer.stop();
	rs[ri] = optr;
      }
    }

#ifdef NGTQ_QBG
    void optimize(const std::string indexPath, bool random, size_t threadSize = 0) {
      NGT::StdOstreamRedirector redirector(silence);
      redirector.begin();
      {
	QBG::Index index(indexPath);
	if (index.getQuantizer().objectList.size() <= 1) {
	  NGTThrowException("optimize: No objects");
	}
	if (numberOfObjects == 0) {
	  numberOfObjects = index.getQuantizer().objectList.size() - 1;
	}
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
	const std::string object = QBG::Index::getTrainObjectFileName(indexPath);
	std::ofstream ofs;
	ofs.open(object);
	index.extract(ofs, numberOfObjects, random);
	if (globalType == GlobalTypeZero) {
	  assert(index.getQuantizer().objectList.pseudoDimension != 0);
	  std::vector<std::vector<float>> global(1);
	  global[0].resize(index.getQuantizer().property.dimension, 0.0);
	  NGT::Clustering::saveVectors(QBG::Index::getQuantizerCodebookFileName(indexPath), global);
	} else if (globalType == GlobalTypeMean) {
	  std::vector<std::vector<float>> vectors;
	  std::string objects = QBG::Index::getTrainObjectFileName(indexPath);
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
	  NGT::Clustering::saveVectors(QBG::Index::getQuantizerCodebookFileName(indexPath), global);
	}
      }

      optimize(indexPath);

      if (globalType == GlobalTypeNone) {
	QBG::Index::load(indexPath, "", "", "", threadSize);
      } else {
	QBG::Index::load(indexPath, QBG::Index::getQuantizerCodebookFileName(indexPath), "", "", threadSize);
      }
      redirector.end();
    }
#endif 

#ifdef NGTQ_QBG
    void optimize(std::string indexPath) {
      std::string object;
      std::string pq;
      std::string global;
      {
	object = QBG::Index::getTrainObjectFileName(indexPath);
	pq = QBG::Index::getPQFileName(indexPath);
	global = QBG::Index::getQuantizerCodebookFileName(indexPath);
      }

      try {
	NGT::Index::mkdir(pq);
      } catch(...) {}
      pq += "/opt.tsv";
      optimize(object, pq, global);
    }
#endif 

    void optimize(std::string invector, std::string ofile, std::string global) {
      vector<vector<float>> vectors;


#ifdef NGT_CLUSTERING
      NGT::Clustering::loadVectors(invector, vectors);
#else
      loadVectors(invector, vectors);
#endif

      if (vectors.size() == 0) {
	std::cerr << "Optimizer: error! the specified vetor file is empty. " << invector << ". the size=" << vectors.size() << std::endl;
	exit(1);
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
      size_t pos = std::distance(std::find(ofile.rbegin(), ofile.rend(), '.'), ofile.rend()) - 1;
      string of = ofile.substr(0, pos);
      if (repositioning) {
	Matrix<float> repositionedR(reposition);
	repositionedR.mul(minR);
	Matrix<float>::save(of + "_R.tsv", repositionedR);
      } else {
	Matrix<float>::save(of + "_R.tsv", minR);
      }
      for (size_t m = 0; m < numberOfSubvectors; m++) {
	stringstream str;
	str << of << "-" << m << ".tsv";
#ifdef NGT_CLUSTERING
	NGT::Clustering::saveClusters(str.str(), minLocalClusters[m]);
#else
	saveClusters(str.str(), minLocalClusters[m]);
#endif
      }

    }

    NGT::Clustering::ClusteringType	clusteringType;
    NGT::Clustering::InitializationMode	initMode;

    NGT::Timer		timelimitTimer;
    float		timelimit;
    size_t		iteration;
    size_t		clusterIteration;		
    bool		clusterSizeConstraint;
    size_t		convergenceLimitTimes;		
    size_t		dim;
    size_t		numberOfObjects;
    size_t		numberOfClusters;
    size_t		numberOfSubvectors;
    size_t		subvectorSize;
    size_t		nOfMatrices;
    float		seedStartObjectSizeRate;
    size_t		seedStep;
    float		reject;
    bool		repositioning;
    bool		rotation;
    GlobalType		globalType;
    bool		silence;
  };
} // namespace NGTQ
