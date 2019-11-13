//
// Copyright (C) 2015-2019 Yahoo Japan Corporation
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
      mergin = 0.2;
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
	auto coefficients = optimizer.adjustSearchEdgeSize(baseAccuracyRange, rateAccuracyRange, numOfQueries, gtEpsilon, mergin);
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

    void execute(
		 const std::string inIndexPath,
		 const std::string outIndexPath
		 ){
      if ((numOfOutgoingEdges < 0 && numOfIncomingEdges >= 0) ||
	  (numOfOutgoingEdges >= 0 && numOfIncomingEdges < 0)) {
	NGTThrowException("Optimizer::execute: Specified any of the number of edges is invalid.");
      }
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
	auto coefficients = optimizer.adjustSearchEdgeSize(baseAccuracyRange, rateAccuracyRange, numOfQueries, gtEpsilon, mergin);
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
	mergin = m;
      }
    }

    size_t numOfOutgoingEdges;
    size_t numOfIncomingEdges;
    std::pair<float, float> baseAccuracyRange;
    std::pair<float, float> rateAccuracyRange;
    size_t numOfQueries;
    double gtEpsilon;
    double mergin;
    bool logDisabled;
  };

}; // NGT

