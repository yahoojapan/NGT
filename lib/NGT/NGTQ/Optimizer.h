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

#ifdef NGT_CLUSTERING
#include "NGT/Clustering.h"
#else
#include "Cluster.h"
#endif

#ifndef NGTQ_BLAS_FOR_ROTATION
#define NGT_DISABLE_BLAS
#endif

#include "Matrix.h"


namespace QBG {
  class BuildParameters;
  class OptimizationParameters;

  class Optimizer {
  public:
    enum GlobalType {
      GlobalTypeNone = 0,
      GlobalTypeZero = 1,
      GlobalTypeMean = 2
    };

    Optimizer() { initialize(); }

    Optimizer(QBG::BuildParameters &param);
    Optimizer(QBG::OptimizationParameters &param);

    void setParameters(QBG::OptimizationParameters &param);

    void initialize();

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
      for (size_t vidx = 0; vidx < vsize; vidx++) {
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

#ifdef NGTQ_QBG
    void evaluate(string global, vector<vector<float>> &vectors, char clusteringType, string &ofile, size_t &numberOfSubvectors, size_t &subvectorSize);

    void evaluate(vector<vector<float>> &vectors, string &ofile, size_t &numberOfSubvectors, size_t &subvectorSize);
#endif

    void
      generateResidualObjects(vector<vector<float>> &globalCentroid, vector<vector<float>> &vectors)
    {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      std::cerr << "generateResidualObjects: Not implemented." << std::endl;
      abort();
#else
      vector<vector<float>> residualVectors;

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

    void
      generateResidualObjects(string global, vector<vector<float>> &vectors)
    {
      if (global.empty()) {
	NGTThrowException("A global codebook is not specified!");
      }
      vector<vector<float>> globalCentroid;
      try {
	NGT::Clustering::loadVectors(global, globalCentroid);
      } catch (...) {
	std::stringstream msg;
	msg << "Optimizer::generateResidualObjects: Cannot load global vectors. " << global;
	NGTThrowException(msg);
      }
      generateResidualObjects(globalCentroid, vectors);
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
				 float clusterSizeConstraintCoefficient,
				 size_t convergenceLimitTimes,
				 double &minDistortion,
				 NGT::Timer &timelimitTimer, float timelimit,
				 bool rotation
				 ) {

      if (numberOfClusters <= 1) {
	std::stringstream msg;
	msg << "Optimizer::optimize: # of clusters is zero or one. " << numberOfClusters;
	NGTThrowException(msg);
      }
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
	  clustering.setClusterSizeConstraintCoefficient(clusterSizeConstraintCoefficient);
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
	quantizedVectors.resize(subQuantizedVectors[0].size());
	for (size_t i = 0; i < quantizedVectors.size(); i++) {
	  quantizedVectors[i].reserve(xp[0].size());
	}
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
		  Matrix<float> &reposition,
		  vector<Matrix<float>> &rs,
		  vector<vector<vector<NGT::Clustering::Cluster>>> &localClusters,
		  vector<double> &errors
		  ) {
      if (vectors.size() == 0) {
	NGTThrowException("the vector is empty");
      }
      if (!reposition.isEmpty()) {
	Matrix<float>::mulSquare(vectors, reposition);
      }
      Matrix<float> xt(vectors);
#ifndef NGTQ_BLAS_FOR_ROTATION
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
			 clusterSizeConstraintCoefficient,
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
    void optimize(const std::string indexPath, size_t threadSize = 0);
    void optimizeWithinIndex(std::string indexPath);
    void optimize(std::string invector, std::string ofile, std::string global);
    void optimize(vector<vector<float>> &vectors, vector<vector<float>> &globalCentroid, Matrix<float> &r, vector<vector<NGT::Clustering::Cluster>> &localClusters, vector<double> &errors);
#endif
    NGT::Timer		timelimitTimer;
    size_t		subvectorSize;

    NGT::Clustering::ClusteringType	clusteringType;
    NGT::Clustering::InitializationMode	initMode;
    size_t		iteration;
    size_t		clusterIteration;		
    bool		clusterSizeConstraint;
    float		clusterSizeConstraintCoefficient;
    size_t		convergenceLimitTimes;		
    size_t		numberOfObjects;
    size_t		numberOfClusters;
    size_t		numberOfSubvectors;
    size_t		numberOfMatrices;
    size_t		seedNumberOfSteps;
    size_t		seedStep;
    float		reject;
    bool		repositioning;
    bool		rotation;
    GlobalType		globalType;
    bool		randomizedObjectExtraction;
    bool		verbose;
    bool		showClusterInfo;
    float		timelimit;
  };
} // namespace QBG
