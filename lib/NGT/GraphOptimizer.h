//
// Copyright (C) 2015-2020 Yahoo Japan Corporation
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

#include	"GraphReconstructor.h"
#include	"Optimizer.h"

namespace NGT {
  class GraphOptimizer {
  public:
    GraphOptimizer(bool unlog = false) {
      init();
      logDisabled = unlog;
    }

    GraphOptimizer(int outgoing, int incoming, int nofqs, 
		   float baseAccuracyFrom, float baseAccuracyTo,
		   float rateAccuracyFrom, float rateAccuracyTo,
		   double qte, double m,
		   bool unlog		// stderr log is disabled.
		   ) {
      init();
      set(outgoing, incoming, nofqs, baseAccuracyFrom, baseAccuracyTo,
	  rateAccuracyFrom, rateAccuracyTo, qte, m);
      logDisabled = unlog;
    }

    void init() {
      numOfOutgoingEdges = 10;
      numOfIncomingEdges= 120;
      numOfQueries = 100;
      baseAccuracyRange = std::pair<float, float>(0.30, 0.50);
      rateAccuracyRange = std::pair<float, float>(0.80, 0.90);
      gtEpsilon = 0.1;
      margin = 0.2;
      logDisabled = false;
    }

    void adjustSearchCoefficients(const std::string indexPath){
      NGT::Index		index(indexPath);
      NGT::GraphIndex	&graph = static_cast<NGT::GraphIndex&>(index.getIndex());
      NGT::Optimizer	optimizer(index);
      if (logDisabled) {
	optimizer.disableLog();
      } else {
	optimizer.enableLog();
      }
      try {
	auto coefficients = optimizer.adjustSearchEdgeSize(baseAccuracyRange, rateAccuracyRange, numOfQueries, gtEpsilon, margin);
	NGT::NeighborhoodGraph::Property &prop = graph.getGraphProperty();
	prop.dynamicEdgeSizeBase = coefficients.first;
	prop.dynamicEdgeSizeRate = coefficients.second;
      } catch(NGT::Exception &err) {
	std::stringstream msg;
	msg << "Optimizer::adjustSearchCoefficients: Cannot adjust the search coefficients. " << err.what();
	NGTThrowException(msg);      
      }
      graph.saveIndex(indexPath);
    }

    static double measureQueryTime(NGT::Index &index, size_t start) {
      NGT::ObjectSpace &objectSpace = index.getObjectSpace();
      NGT::ObjectRepository &objectRepository = objectSpace.getRepository();
      size_t nQueries = 200;
      nQueries = objectRepository.size() - 1 < nQueries ? objectRepository.size() - 1 : nQueries;

      size_t step = objectRepository.size() / nQueries;
      assert(step != 0);
      std::vector<size_t> ids;
      for (size_t startID = start; startID < step; startID++) {
	for (size_t id = startID; id < objectRepository.size(); id += step) {
	  if (!objectRepository.isEmpty(id)) {
	    ids.push_back(id);
	  }
	}
	if (ids.size() >= nQueries) {
	  ids.resize(nQueries);
	  break;
	}
      }
      if (nQueries > ids.size()) {
	std::cerr << "# of Queries is not enough." << std::endl;
	return DBL_MAX;
      }

      NGT::Timer timer;
      timer.reset();
      for (auto id = ids.begin(); id != ids.end(); id++) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	NGT::Object *obj = objectSpace.allocateObject(*objectRepository.get(*id));
	NGT::SearchContainer searchContainer(*obj);
#else
	NGT::SearchContainer searchContainer(*objectRepository.get(*id));
#endif
	NGT::ObjectDistances objects;
	searchContainer.setResults(&objects);
	searchContainer.setSize(10);
	searchContainer.setEpsilon(0.1);
	timer.restart();
	index.search(searchContainer);
	timer.stop();
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	objectSpace.deleteObject(obj);
#endif
      }
      return timer.time * 1000.0;
    }

    static std::pair<size_t, double> searchMinimumQueryTime(NGT::Index &index, size_t prefetchOffset, 
							    int maxPrefetchSize, size_t seedID) {
      NGT::ObjectSpace &objectSpace = index.getObjectSpace();
      int step = 256;
      int prevPrefetchSize = 64;
      size_t minPrefetchSize = 0;
      double minTime = DBL_MAX;
      for (step = 256; step != 32; step /= 2) {
        double prevTime = DBL_MAX;
        for (int prefetchSize = prevPrefetchSize - step < 64 ? 64 : prevPrefetchSize - step; prefetchSize <= maxPrefetchSize; prefetchSize += step) {
          objectSpace.setPrefetchOffset(prefetchOffset);
          objectSpace.setPrefetchSize(prefetchSize);
          double time = measureQueryTime(index, seedID);
          if (prevTime < time) {
            break;
          }
          prevTime = time;
          prevPrefetchSize = prefetchSize;
        }
        if (minTime > prevTime) {
          minTime = prevTime;
          minPrefetchSize = prevPrefetchSize;
        }
      }
      return std::make_pair(minPrefetchSize, minTime);
    }

    static std::pair<size_t, size_t> adjustPrefetchParameters(NGT::Index &index) {

      bool gridSearch = false;
      {
	double time = measureQueryTime(index, 1);
	if (time < 500.0) {
	  gridSearch = true;
	}
      }

      size_t prefetchOffset = 0;
      size_t prefetchSize = 0;
      std::vector<std::pair<size_t, size_t>> mins;
      NGT::ObjectSpace &objectSpace = index.getObjectSpace();
      int maxSize = objectSpace.getByteSizeOfObject() * 4;
      maxSize = maxSize < 64 * 28 ? maxSize : 64 * 28; 
      for (int trial = 0; trial < 10; trial++) {
	size_t minps = 0;
	size_t minpo = 0;
	if (gridSearch) {
	  double minTime = DBL_MAX;
	  for (size_t po = 1; po <= 10; po++) {
	    auto min = searchMinimumQueryTime(index, po, maxSize, trial + 1);
	    if (minTime > min.second) {
	      minTime = min.second;
	      minps = min.first;
	      minpo = po;
	    }
	  }
	} else {
	  double prevTime = DBL_MAX;
	  for (size_t po = 1; po <= 10; po++) {
	    auto min = searchMinimumQueryTime(index, po, maxSize, trial + 1);
	    if (prevTime < min.second) {
	      break;
	    }
	    prevTime = min.second;
	    minps = min.first;
	    minpo = po;
	  }
	}
	if (std::find(mins.begin(), mins.end(), std::make_pair(minpo, minps)) != mins.end()) {
	  prefetchOffset = minpo;
	  prefetchSize = minps;
	  mins.push_back(std::make_pair(minpo, minps));
	  break;
	}
	mins.push_back(std::make_pair(minpo, minps));
      }
      return std::make_pair(prefetchOffset, prefetchSize);
    }

    void execute(
		 const std::string inIndexPath,
		 const std::string outIndexPath
		 ){
      if ((numOfOutgoingEdges < 0 && numOfIncomingEdges >= 0) ||
	  (numOfOutgoingEdges >= 0 && numOfIncomingEdges < 0)) {
	NGTThrowException("Optimizer::execute: Specified any of the number of edges is invalid.");
      }
      {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	if (access(outIndexPath.c_str(), 0) == 0) {
	  std::stringstream msg;
	  msg << "Optimizer::execute: The specified index exists. " << outIndexPath;
	  NGTThrowException(msg);
	}
	const std::string com = "cp -r " + inIndexPath + " " + outIndexPath;
	system(com.c_str());
	NGT::Index	outIndex(outIndexPath);
#else
	NGT::Index	outIndex(inIndexPath);
#endif
	NGT::GraphIndex	&outGraph = static_cast<NGT::GraphIndex&>(outIndex.getIndex());
	NGT::Timer timer;
	timer.start();
	std::vector<NGT::ObjectDistances> graph;
	NGT::StdOstreamRedirector redirector(logDisabled);
	redirector.begin();
	try {
	  std::cerr << "Optimizer::execute: Extract the graph data." << std::endl;
	  // extract only edges from the index to reduce the memory usage.
	  NGT::GraphReconstructor::extractGraph(graph, outIndex);
	  if (numOfOutgoingEdges >= 0) {
	    NGT::GraphReconstructor::convertToANNG(graph);
	    NGT::GraphReconstructor::reconstructGraph(graph, outIndex, numOfOutgoingEdges, numOfIncomingEdges);
	  }
	  timer.stop();
	  std::cerr << "Optimizer::execute: Graph reconstruction time=" << timer.time << " (sec) " << std::endl;
	  timer.reset();
	  timer.start();
	  NGT::GraphReconstructor::adjustPathsEffectively(outIndex);
	  timer.stop();
	  std::cerr << "Optimizer::execute: Path adjustment time=" << timer.time << " (sec) " << std::endl;
	} catch (NGT::Exception &err) {
	  redirector.end();
	  throw(err);
	}
	redirector.end();
	NGT::Optimizer optimizer(outIndex);
	if (logDisabled) {
	  optimizer.disableLog();
	} else {
	  optimizer.enableLog();
	}
	try {
	  auto coefficients = optimizer.adjustSearchEdgeSize(baseAccuracyRange, rateAccuracyRange, numOfQueries, gtEpsilon, margin);
	  NGT::NeighborhoodGraph::Property &prop = outGraph.getGraphProperty();
	  prop.dynamicEdgeSizeBase = coefficients.first;
	  prop.dynamicEdgeSizeRate = coefficients.second;
	} catch(NGT::Exception &err) {
	  std::stringstream msg;
	  msg << "Optimizer::execute: Cannot adjust the search coefficients. " << err.what();
	  NGTThrowException(msg);      
	}

	outGraph.saveIndex(outIndexPath);
      }

      try {
	NGT::Index	outIndex(outIndexPath, true);
	auto prefetch = adjustPrefetchParameters(outIndex);
	NGT::Property prop;
	outIndex.getProperty(prop);
	prop.prefetchOffset = prefetch.first;
	prop.prefetchSize = prefetch.second;
	outIndex.setProperty(prop);
	static_cast<NGT::GraphIndex&>(outIndex.getIndex()).saveProperty(outIndexPath);
      } catch(NGT::Exception &err) {
	std::stringstream msg;
	msg << "Optimizer::execute: Cannot adjust prefetch parameters. " << err.what();
	NGTThrowException(msg);
      }

    }

    void set(int outgoing, int incoming, int nofqs, 
	     float baseAccuracyFrom, float baseAccuracyTo,
	     float rateAccuracyFrom, float rateAccuracyTo,
	     double qte, double m
	     ) {
      if (outgoing >= 0) {
	numOfOutgoingEdges = outgoing;
      }
      if (incoming >= 0) {
	numOfIncomingEdges = incoming;
      }
      if (nofqs > 0) {
	numOfQueries = nofqs;
      }
      if (baseAccuracyFrom > 0.0) {
	baseAccuracyRange.first = baseAccuracyFrom;
      }
      if (baseAccuracyTo > 0.0) {
	baseAccuracyRange.second = baseAccuracyTo;
      }
      if (rateAccuracyFrom > 0.0) {
	rateAccuracyRange.first = rateAccuracyFrom;
      }
      if (rateAccuracyTo > 0.0) {
	rateAccuracyRange.second = rateAccuracyTo;
      }
      if (qte >= -1.0) {
	gtEpsilon = qte;
      }
      if (m > 0.0) {
	margin = m;
      }
    }

    size_t numOfOutgoingEdges;
    size_t numOfIncomingEdges;
    std::pair<float, float> baseAccuracyRange;
    std::pair<float, float> rateAccuracyRange;
    size_t numOfQueries;
    double gtEpsilon;
    double margin;
    bool logDisabled;
  };

}; // NGT

