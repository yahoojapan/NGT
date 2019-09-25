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

#pragma once

#include "Command.h"

#define NGT_LOG_BASED_OPTIMIZATION

namespace NGT {
  class Optimizer {
  public:


    Optimizer(NGT::Index &i, size_t n = 10):index(i), nOfResults(n) { 
    }
    ~Optimizer() {}

    class MeasuredValue {
    public:
    MeasuredValue():keyValue(0.0), totalCount(0), meanAccuracy(0.0), meanTime(0.0), meanDistanceCount(0.0), meanVisitCount(0.0) {}
      double	keyValue;
      size_t	totalCount;
      float	meanAccuracy;
      double	meanTime;
      double	meanDistanceCount;
      double	meanVisitCount;
    };

    void enableLog() { redirector.disable(); }
    void disableLog() { redirector.enable(); }

    static void search(NGT::Index &index, istream &gtStream, Command::SearchParameter &sp, vector<MeasuredValue> &acc) {
      ifstream		is(sp.query);
      if (!is) {
	stringstream msg;
	msg << "Cannot open the specified file. " << sp.query;
	NGTThrowException(msg);
      }

      search(index, gtStream, sp, acc);
    }

    static void search(NGT::Index &index, istream &queries, istream &gtStream, Command::SearchParameter &sp, vector<MeasuredValue> &acc) {
      sp.stepOfEpsilon = 1.0;
      stringstream resultStream;
      NGT::Command::search(index, sp, queries, resultStream);
      resultStream.clear();
      resultStream.seekg(0, ios_base::beg);
      string type;
      size_t actualResultSize = 0;
      gtStream.seekg(0, ios_base::end);      
      auto pos = gtStream.tellg();
      if (pos == 0) {
	evaluate(resultStream, acc, type, actualResultSize);
      } else {
	gtStream.clear();
	gtStream.seekg(0, ios_base::beg);
	evaluate(gtStream, resultStream, acc, type, actualResultSize);
      }

      assert(acc.size() == 1);
    }

    static void
      evaluate(istream &resultStream, vector<MeasuredValue> &accuracies, string &type, 
	       size_t &resultDataSize, size_t specifiedResultSize = 0, size_t groundTruthSize = 0, bool recall = false)
    {

      resultDataSize = 0;

      if (recall) {
	if (specifiedResultSize == 0) {
	  stringstream msg;
	  msg << "For calculating recalls, the result size should be specified.";
	  NGTThrowException(msg);
	}
	resultDataSize = specifiedResultSize;
      } else {
	checkAndGetSize(resultStream, resultDataSize);
      }

      string line;
      size_t queryNo = 1;
      map<double, double> totalAccuracy;
      map<double, double> totalTime;
      map<double, size_t> totalDistanceCount;
      map<double, size_t> totalVisitCount;
      map<double, size_t> totalCount;

      resultStream.clear();
      resultStream.seekg(0, ios_base::beg);

      do {
	unordered_set<size_t> gt;
	double furthestDistance = 0.0;
	sumup(resultStream, queryNo, totalAccuracy, totalTime, totalDistanceCount, totalVisitCount, totalCount, 
	      gt, resultDataSize, type, recall, furthestDistance);
	queryNo++;
      } while (!resultStream.eof());

      accuracies.clear();
      for (auto it = totalAccuracy.begin(); it != totalAccuracy.end(); ++it) {
	MeasuredValue a;
	a.keyValue = (*it).first;
	a.totalCount = totalCount[a.keyValue];
	a.meanAccuracy = totalAccuracy[a.keyValue] / (double)totalCount[a.keyValue];
	a.meanTime = totalTime[a.keyValue] / (double)totalCount[a.keyValue];
	a.meanDistanceCount = (double)totalDistanceCount[a.keyValue] / (double)totalCount[a.keyValue];
	a.meanVisitCount = (double)totalVisitCount[a.keyValue] / (double)totalCount[a.keyValue];
	accuracies.push_back(a);
      }
    }

    static void
      evaluate(istream &gtStream, istream &resultStream, vector<MeasuredValue> &accuracies, string &type, 
	       size_t &resultDataSize, size_t specifiedResultSize = 0, size_t groundTruthSize = 0, bool recall = false)
    {

      resultDataSize = 0;

      if (recall) {
	if (specifiedResultSize == 0) {
	  stringstream msg;
	  msg << "For calculating recalls, the result size should be specified.";
	  NGTThrowException(msg);
	}
	resultDataSize = specifiedResultSize;
      } else {
	checkAndGetSize(resultStream, resultDataSize);
      }

      string line;
      size_t queryNo = 1;
      map<double, double> totalAccuracy;
      map<double, double> totalTime;
      map<double, size_t> totalDistanceCount;
      map<double, size_t> totalVisitCount;
      map<double, size_t> totalCount;

      resultStream.clear();
      resultStream.seekg(0, ios_base::beg);

      while (getline(gtStream, line)) {
	vector<string> tokens;
	NGT::Common::tokenize(line, tokens, "=");
	if (tokens.size() == 0) {
	  continue;
	}
	if (tokens[0] == "# Query No.") {
	  if (tokens.size() > 1 && (size_t)NGT::Common::strtol(tokens[1]) == queryNo) {
	    unordered_set<size_t> gt;
	    double furthestDistance;
	    if (groundTruthSize == 0) {
	      loadGroundTruth(gtStream, gt, resultDataSize, furthestDistance);
	    } else {
	      loadGroundTruth(gtStream, gt, groundTruthSize, furthestDistance);
	    }
	    sumup(resultStream, queryNo, totalAccuracy, totalTime, totalDistanceCount, totalVisitCount, totalCount, 
		  gt, resultDataSize, type, recall, furthestDistance);
	    queryNo++;
	  }
	}
      }

      accuracies.clear();
      for (auto it = totalAccuracy.begin(); it != totalAccuracy.end(); ++it) {
	MeasuredValue a;
	a.keyValue = (*it).first;
	a.totalCount = totalCount[a.keyValue];
	a.meanAccuracy = totalAccuracy[a.keyValue] / (double)totalCount[a.keyValue];
	a.meanTime = totalTime[a.keyValue] / (double)totalCount[a.keyValue];
	a.meanDistanceCount = (double)totalDistanceCount[a.keyValue] / (double)totalCount[a.keyValue];
	a.meanVisitCount = (double)totalVisitCount[a.keyValue] / (double)totalCount[a.keyValue];
	accuracies.push_back(a);
      }
    }


    static void
      loadGroundTruth(istream & gtf, unordered_set<size_t> & gt, size_t resultDataSize, double &distance) {
      string line;
      size_t dataCount = 0;
      size_t searchCount = 0;
      while (getline(gtf, line)) {
	if (line.size() != 0 && line.at(0) == '#') {
	  vector<string> gtf;
	  NGT::Common::tokenize(line, gtf, "=");
	  if (gtf.size() >= 1) {
	    if (gtf[0] == "# End of Search") {
	      searchCount++;
	    }
	    if (gtf[0] == "# End of Query") {
	      if (searchCount != 1) {
		stringstream msg;
		msg << "Error: gt has not just one search result.";
		NGTThrowException(msg);
	      }
	      if (dataCount < resultDataSize) {
		stringstream msg;
		msg << "Error: gt data is less than result size! " << dataCount << ":" << resultDataSize;
		NGTThrowException(msg);
	      }
	      return;
	    }
	    continue;
	  }
	}
	dataCount++;
	if (dataCount > resultDataSize) {
	  continue;
	}
	vector<string> result;      
	NGT::Common::tokenize(line, result, " \t");
	if (result.size() < 3) {
	  stringstream msg;
	  msg << "result format is wrong. ";
	  NGTThrowException(msg);
	}
	size_t id = NGT::Common::strtol(result[1]);
	distance = NGT::Common::strtod(result[2]);
	try {
	  gt.insert(id);
	} catch(...) {
	  stringstream msg;
	  msg << "Cannot insert id into the gt. " << id;
	  NGTThrowException(msg);
	}
      } 
    }

    static void checkAndGetSize(istream &resultStream, size_t &resultDataSize)
    {
      size_t lineCount = 0;
      size_t prevDataCount = 0;
      string line;
      bool warn = false;

      while (getline(resultStream, line)) {
	lineCount++;
	if (line.size() != 0 && line.at(0) == '#') {
	  vector<string> tf;
	  NGT::Common::tokenize(line, tf, "=");
	  if (tf.size() >= 1 && tf[0] == "# Query No.") {
	    size_t dataCount = 0;
	    string lastDataLine;
	    while (getline(resultStream, line)) {
	      lineCount++;
	      if (line.size() != 0 && line.at(0) == '#') {
		vector<string> gtf;
		NGT::Common::tokenize(line, gtf, "=");
		if (gtf.size() >= 1 && gtf[0] == "# End of Search") {
		  if (prevDataCount == 0) {
		    prevDataCount = dataCount;
		  } else {
		    if (prevDataCount != dataCount) {
		      warn = true;
		      cerr << "Warning!: Result sizes are inconsistent! $prevDataCount:$dataCount" << endl;;
		      cerr << "  Line No." << lineCount << ":"  << lastDataLine << endl;
		      if (prevDataCount < dataCount) {
			prevDataCount = dataCount;
		      }
		    }
		  }
		  dataCount = 0;
		  break;
		}
		continue;
	      }
	      lastDataLine = line;
	      vector<string> result;      
	      NGT::Common::tokenize(line, result, " \t");
	      if (result.size() < 3) {
		stringstream msg;
		msg << "result format is wrong. ";
		NGTThrowException(msg);
	      }
	      size_t rank = NGT::Common::strtol(result[0]);
	      dataCount++;
	      if (rank != dataCount) {
		stringstream msg;
		msg << "check: inner error! " << rank << ":" << dataCount;
		NGTThrowException(msg);
	      }
	    }
	  }
	}
      }
      resultDataSize = prevDataCount;
      if (warn) {
	cerr << "Warning! ****************************************************************************" << endl;
	cerr << " Check if the result number $$resultDataSize is correct." << endl;
	cerr << "Warning! ****************************************************************************" << endl;
      }
    }

    static void sumup(istream &resultStream, 
		      size_t queryNo, 
		      map<double, double> &totalAccuracy, 
		      map<double, double> &totalTime,
		      map<double, size_t> &totalDistanceCount,
		      map<double, size_t> &totalVisitCount,
		      map<double, size_t> &totalCount,
		      unordered_set<size_t> &gt,
		      const size_t resultDataSize,
		      string &keyValue,
		      bool recall,
		      double furthestDistance)
    {
      string line;
      size_t lineNo = 0;
      while (getline(resultStream, line)) {
	lineNo++;
	size_t resultNo = 0;
	if (line.size() != 0 && line.at(0) == '#') {
	  vector<string> tf;
	  NGT::Common::tokenize(line, tf, "=");
	  if (tf.size() >= 1 && tf[0] == "# Query No.") {
	    if (tf.size() >= 2 && (size_t)NGT::Common::strtol(tf[1]) == queryNo) {
	      size_t relevantCount = 0;
	      size_t dataCount = 0;
	      string epsilon;
	      string expansion;  
	      double queryTime = 0.0;
	      size_t distanceCount = 0;
	      size_t visitCount = 0;
	      while (getline(resultStream, line)) {
		lineNo++;
		if (line.size() != 0 && line.at(0) == '#') {
		  vector<string> gtf;
		  NGT::Common::tokenize(line, gtf, "=");
		  if (gtf.size() >= 2 && (gtf[0] == "# Epsilon" || gtf[0] == "# Factor")) {
		    epsilon = gtf[1];
		  } else if (gtf.size() >= 2 && gtf[0] == "# Result expansion") {
		    expansion = gtf[1];
		  } else if (gtf.size() >= 2 && gtf[0] == "# Query Time (msec)") {
		    queryTime = NGT::Common::strtod(gtf[1]);
		  } else if (gtf.size() >= 2 && gtf[0] == "# Distance Computation") {
		    distanceCount = NGT::Common::strtol(gtf[1]);
		  } else if (gtf.size() >= 2 && gtf[0] == "# Visit Count") {
		    visitCount = NGT::Common::strtol(gtf[1]);
		  } else if (gtf.size() >= 1 && gtf[0] == "# End of Query") {
		    return;
		  } else if (gtf.size() >= 1 && gtf[0] == "# End of Search") {
		    resultNo++;
		    if (recall == false && resultDataSize != dataCount) {
		      cerr << "Warning! ****************************************************************************" << endl;
		      cerr << "  Use $resultDataSize instead of $dataCount as the result size to compute accuracy. " <<  endl;
		      cerr << "    # of the actual search resultant objects=" << dataCount << endl;
		      cerr << "    the specified # of search objects or # of the ground truth data=" << resultDataSize << endl;
		      cerr << "    Line No.=" << lineNo << " Query No.=" << queryNo << " Result No.=" << resultNo << endl;
		      cerr << "Warning! ****************************************************************************" << endl;
		    }
		    double accuracy = (double)relevantCount / (double)resultDataSize;
		    double key;
		    if (epsilon != "") {
		      key = NGT::Common::strtod(epsilon);
		      keyValue = "Factor (Epsilon)";
		    } else if (expansion != "") {
		      key = NGT::Common::strtod(expansion);
		      keyValue = "Expansion";
		    } else {
		      stringstream msg;
		      msg << "check: inner error! " << epsilon;
		      cerr << "Cannot find epsilon.";
		      NGTThrowException(msg);
		    }
		    {
		      auto di = totalAccuracy.find(key);
		      if (di != totalAccuracy.end()) {
			(*di).second += accuracy;
		      } else {
			totalAccuracy.insert(std::make_pair(key, accuracy));
		      }
		    }
		    {
		      auto di = totalTime.find(key);
		      if (di != totalTime.end()) {
			(*di).second += queryTime;
		      } else {
			totalTime.insert(std::make_pair(key, queryTime));
		      }
		    }
		    {
		      auto di = totalDistanceCount.find(key);
		      if (di != totalDistanceCount.end()) {
			(*di).second += distanceCount;
		      } else {
			totalDistanceCount.insert(std::make_pair(key, distanceCount));
		      }
		    }
		    {
		      auto di = totalVisitCount.find(key);
		      if (di != totalVisitCount.end()) {
			(*di).second += visitCount;
		      } else {
			totalVisitCount.insert(std::make_pair(key, visitCount));
		      }
		    }
		    {
		      auto di = totalCount.find(key);
		      if (di != totalCount.end()) {
			(*di).second ++;
		      } else {
			totalCount.insert(std::make_pair(key, 1));
		      }
		    }
		    relevantCount = 0;
		    dataCount = 0;
		  } 
		  continue;
		} 
		vector<string> result;      
		NGT::Common::tokenize(line, result, " \t");
		if (result.size() < 3) {
		  cerr << "result format is wrong. " << endl;
		  abort();
		}
		size_t rank = NGT::Common::strtol(result[0]);
		size_t id = NGT::Common::strtol(result[1]);
		double distance = NGT::Common::strtod(result[2]);
		if (gt.count(id) != 0) {
		  relevantCount++;
		} else {
		  if (furthestDistance > 0.0 && distance <= furthestDistance) {
		    relevantCount++;
		    if (distance < furthestDistance) {
		      //cerr << "Optimizer:Warning!. The ground truth has a missing object. " << id << ":" << distance << ":" << furthestDistance << endl;
		    }
		  }
		}
		dataCount++;
		if (rank != dataCount) {
		  cerr << "inner error! $rank $dataCount !!" << endl;;
		  abort();
		}
	      } 
	    } else { 
	      cerr << "Fatal error! : Cannot find query No. " << queryNo << endl;
	      abort();
	    } 
	  } 
	} 
      } 
    }

    static void exploreEpsilonForAccuracy(NGT::Index &index, istream &queries, istream &gtStream, 
					  Command::SearchParameter &sp, pair<float, float> accuracyRange, double mergin) 
    {
      double fromUnder = 0.0;
      double fromOver = 1.0;
      double toUnder = 0.0;
      double toOver = 1.0;
      float fromUnderEpsilon = -0.9;
      float fromOverEpsilon = -0.9;
      float toUnderEpsilon = -0.9;
      float toOverEpsilon = -0.9;

      float accuracyRangeFrom = accuracyRange.first;
      float accuracyRangeTo = accuracyRange.second;

      double range = accuracyRangeTo - accuracyRangeFrom;

      vector<MeasuredValue> acc;

      {
	float startEpsilon = -0.6;
	float epsilonStep = 0.02;
	size_t count;
	for (count = 0;; count++) {
	  float epsilon = round((startEpsilon + epsilonStep * count) * 100.0F) / 100.0F; 
	  if (epsilon > 0.25F) {
	    stringstream msg;
	    msg << "exploreEpsilonForAccuracy:" << endl;
	    msg << "Error!! Epsilon (lower bound) is too large. " << epsilon << "," << startEpsilon << "," << epsilonStep << "," << count;
	    NGTThrowException(msg);
	  }
	  acc.clear();
	  sp.beginOfEpsilon = sp.endOfEpsilon = fromOverEpsilon = epsilon;
	  queries.clear();
	  queries.seekg(0, ios_base::beg);
	  search(index, queries, gtStream, sp, acc);
	  if (acc[0].meanAccuracy >= accuracyRangeFrom) {
	    break;
	  }
	}
	if (fromOverEpsilon == startEpsilon) {
	  stringstream msg;
	  msg << "exploreEpsilonForAccuracy:" << endl;
	  msg << "Error! startEpsilon should be reduced for the specified range.";
	  NGTThrowException(msg);
	}
	fromOver = acc[0].meanAccuracy;

	if (fromOver < accuracyRangeTo) {
	  startEpsilon = fromOverEpsilon;
	  for (count = 0;; count++) {
	    float epsilon = round((startEpsilon + epsilonStep * count) * 100.0F) / 100.0F; 
	    sp.beginOfEpsilon = sp.endOfEpsilon = toOverEpsilon = epsilon;
	    if (epsilon > 0.25F) {
	      stringstream msg;
	      msg << "exploreEpsilonForAccuracy:" << endl;
	      msg << "Error!! Epsilon (upper bound) is too large. " << epsilon << "," << startEpsilon << "," << epsilonStep << "," << count;
	      NGTThrowException(msg);
	    }
	    acc.clear();
	    queries.clear();
	    queries.seekg(0, ios_base::beg);
	    search(index, queries, gtStream, sp, acc);
	    epsilon += epsilonStep;
	    if (acc[0].meanAccuracy >= accuracyRangeTo) {
	      break;
	    }
	  }
	  toOver = acc[0].meanAccuracy;
	} else {
	  toOver = fromOver;
	  toOverEpsilon = fromOverEpsilon;
	}
	fromUnderEpsilon = fromOverEpsilon - epsilonStep;
      }
      sp.beginOfEpsilon = sp.endOfEpsilon = fromUnderEpsilon;
      while (true) {
	acc.clear();
	queries.clear();
	queries.seekg(0, ios_base::beg);
	search(index, queries, gtStream, sp, acc);
	if (acc[0].meanAccuracy >= fromUnder && acc[0].meanAccuracy <= accuracyRangeFrom) {
	  fromUnder = acc[0].meanAccuracy;
	  fromUnderEpsilon = acc[0].keyValue;
	}
	if (acc[0].meanAccuracy <= fromOver && acc[0].meanAccuracy > accuracyRangeFrom) {
	  fromOver = acc[0].meanAccuracy;
	  fromOverEpsilon = acc[0].keyValue;
	}
	if (acc[0].meanAccuracy <= toOver && acc[0].meanAccuracy > accuracyRangeTo) {
	  toOver = acc[0].meanAccuracy;
	  toOverEpsilon = acc[0].keyValue;
	}
	if (acc[0].meanAccuracy >= toUnder && acc[0].meanAccuracy <= accuracyRangeTo) {
	  toUnder = acc[0].meanAccuracy;
	  toUnderEpsilon = acc[0].keyValue;
	}

	if (fromUnder < accuracyRangeFrom - range * mergin) {
	  if ((fromUnderEpsilon + fromOverEpsilon) / 2.0 == sp.beginOfEpsilon) {
	    stringstream msg;
	    msg << "exploreEpsilonForAccuracy:" << endl;
	    msg << "Error!! Not found proper under epsilon for mergin=" << mergin << " and the number of queries." << endl;
	    msg << "        Should increase mergin or the number of queries to get the proper epsilon. ";
	    NGTThrowException(msg);
	  } else {
	    sp.beginOfEpsilon = sp.endOfEpsilon = (fromUnderEpsilon + fromOverEpsilon) / 2.0;
	  }
	} else if (toOver > accuracyRangeTo + range * mergin) {
	  if ((toUnderEpsilon + toOverEpsilon) / 2.0 == sp.beginOfEpsilon) {
	    stringstream msg;
	    msg << "exploreEpsilonForAccuracy:" << endl;
	    msg << "Error!! Not found proper over epsilon for mergin=" << mergin << " and the number of queries." << endl;
	    msg << "        Should increase mergin or the number of queries to get the proper epsilon. ";
	    NGTThrowException(msg);
	  } else {
	    sp.beginOfEpsilon = sp.endOfEpsilon = (toUnderEpsilon + toOverEpsilon) / 2.0;
	  }
	} else {
	  if (fromUnderEpsilon == toOverEpsilon) {
	    stringstream msg;
	    msg << "exploreEpsilonForAccuracy:" << endl;
	    msg << "Error!! From and to epsilons are the same. Cannot continue.";
	    NGTThrowException(msg);
	  }
	  sp.beginOfEpsilon = fromUnderEpsilon;
	  sp.endOfEpsilon = toOverEpsilon;
	  return;
	}
      }
      stringstream msg;
      msg << "Something wrong!";
      NGTThrowException(msg);
    }

    MeasuredValue measure(istream &queries, istream &gtStream, Command::SearchParameter &searchParameter, pair<float, float> accuracyRange, double mergin) {

      exploreEpsilonForAccuracy(index, queries, gtStream, searchParameter, accuracyRange, mergin);
    
      stringstream resultStream;
      queries.clear();
      queries.seekg(0, ios_base::beg);
      NGT::Command::search(index, searchParameter, queries, resultStream);
      gtStream.clear();
      gtStream.seekg(0, ios_base::beg);
      resultStream.clear();
      resultStream.seekg(0, ios_base::beg);
      string type;
      vector<MeasuredValue> accuracies;
      size_t actualResultSize = 0;
      evaluate(gtStream, resultStream, accuracies, type, actualResultSize);
      size_t size;
      double distanceCount, visitCount, time;
      calculateMeanValues(accuracies, accuracyRange.first, accuracyRange.second, size, distanceCount, visitCount, time);
      if (distanceCount == 0) {
	stringstream msg;
	msg << "measureDistance: Error! Distance count is zero.";
	NGTThrowException(msg);
      }
      MeasuredValue v;
      v.meanVisitCount = visitCount;
      v.meanDistanceCount = distanceCount;
      v.meanTime = time;
      return v;
    }

    pair<size_t, double> adjustBaseSearchEdgeSize(stringstream &queries, Command::SearchParameter &searchParameter, stringstream &gtStream, pair<float, float> accuracyRange, float merginInit = 0.2, size_t prevBase = 0) {
      searchParameter.edgeSize = -2;
      size_t minimumBase = 4;
      size_t minimumStep = 2;
      size_t baseStartInit = 1;
      while (prevBase != 0) {
	prevBase >>= 1;
	baseStartInit <<= 1;
      }
      baseStartInit >>= 2;
      baseStartInit = baseStartInit < minimumBase ? minimumBase : baseStartInit;
      while(true) {
	try {
	  float mergin = merginInit;
	  size_t baseStart = baseStartInit;
	  double minTime = DBL_MAX;
	  size_t minBase = 0;
	  map<size_t, double> times;
	  cerr << "adjustBaseSearchEdgeSize::explore for the mergin " << mergin << ", " << baseStart << "..." << endl;
	  for (size_t baseStep = 16; baseStep != 1; baseStep /= 2) {
	    double prevTime = DBL_MAX;
	    for (size_t base = baseStart; ; base += baseStep) {
	      if (base > 1000) {
		stringstream msg;
		msg << "base is too large! " << base;
		NGTThrowException(msg);
	      }
	      searchParameter.step = 10;
	      NGT::GraphIndex &graphIndex = static_cast<GraphIndex&>(index.getIndex());
	      NeighborhoodGraph::Property &prop = graphIndex.getGraphProperty();
	      prop.dynamicEdgeSizeBase = base;
	      double time;
	      if (times.count(base) == 0) {
		for (;;) {
		  try {
		    auto values = measure(queries, gtStream, searchParameter, accuracyRange, mergin);
		    time = values.meanTime;
		    break;
		  } catch(NGT::Exception &err) {
		    if (err.getMessage().find("Error!! Epsilon") != std::string::npos &&
			err.getMessage().find("is too large") != std::string::npos) {
		      cerr << "Warning: Cannot adjust the base edge size." << err.what() << endl;
		      cerr << "Try again with the next base" << endl;
		      NGTThrowException("**Retry**"); 
		    }
		    if (mergin > 0.4) {
		      cerr << "Warning: Cannot adjust the base even for the widest mergin " << mergin << ". " << err.what();
		      NGTThrowException("**Retry**"); 
		    } else {
		      cerr << "Warning: Cannot adjust the base edge size for mergin " << mergin << ". " << err.what() << endl;
		      cerr << "Try again for the next mergin." << endl;
		      mergin += 0.05;
		    }
		  }
		}
		times.insert(std::make_pair(base, time));
		cerr << "adjustBaseSearchEdgeSize::base=" << base << ", query time=" << time << endl;
	      } else {
		time = times.at(base);
	      }
	      if (prevTime <= time) {
		if (baseStep == minimumStep) {
		  return std::make_pair(minBase, minTime);
		} else {
		  baseStart = static_cast<int>(minBase) - static_cast<int>(baseStep) < static_cast<int>(baseStart) ? baseStart : minBase - baseStep;
		  break;
		}
	      }
	      prevTime = time;
	      if (time < minTime) {
		minTime = time;
		minBase = base;
	      }
	    }
	  }
	} catch(NGT::Exception &err) {
	  if (err.getMessage().find("**Retry**") != std::string::npos) {
	    baseStartInit += minimumStep;
	  } else {
	    throw err;
	  }
	}
      }
    }

    size_t adjustBaseSearchEdgeSize(pair<float, float> accuracyRange, size_t querySize, double epsilon, float mergin = 0.2) {
      cerr << "adjustBaseSearchEdgeSize::Extract queries for GT..." << endl;
      stringstream queries;
      extractQueries(querySize, queries);

      cerr << "adjustBaseSearchEdgeSize::create GT..." << endl;
      Command::SearchParameter searchParameter;
      stringstream gtStream;
      createGroundTruth(index, epsilon, searchParameter, queries, gtStream);

      auto base = adjustBaseSearchEdgeSize(queries, searchParameter, gtStream, accuracyRange, mergin);
      return base.first;
    }


    pair<size_t, double> adjustRateSearchEdgeSize(stringstream &queries, Command::SearchParameter &searchParameter, stringstream &gtStream, pair<float, float> accuracyRange, float merginInit = 0.2, size_t prevRate = 0) {
      searchParameter.edgeSize = -2;
      size_t minimumRate = 2;
      size_t minimumStep = 4;
      size_t rateStartInit = 1;
      while (prevRate != 0) {
	prevRate >>= 1;
	rateStartInit <<= 1;
      }
      rateStartInit >>= 2;
      rateStartInit = rateStartInit < minimumRate ? minimumRate : rateStartInit;
      while (true) {
	try {
	  float mergin = merginInit;
	  size_t rateStart = rateStartInit;
	  double minTime = DBL_MAX;
	  size_t minRate = 0;
	  map<size_t, double> times;
	  cerr << "adjustRateSearchEdgeSize::explore for the mergin " << mergin << ", " << rateStart << "..." << endl;
	  for (size_t rateStep = 16; rateStep != 1; rateStep /= 2) {
	    double prevTime = DBL_MAX;
	    for (size_t rate = rateStart; rate < 2000; rate += rateStep) {
	      if (rate > 1000) {
		stringstream msg;
		msg << "rate is too large! " << rate;
		NGTThrowException(msg);
	      }
	      searchParameter.step = 10;
	      NGT::GraphIndex &graphIndex = static_cast<GraphIndex&>(index.getIndex());
	      NeighborhoodGraph::Property &prop = graphIndex.getGraphProperty();
	      prop.dynamicEdgeSizeRate = rate;
	      double time;
	      if (times.count(rate) == 0) {
		for (;;) {
		  try {
		    auto values = measure(queries, gtStream, searchParameter, accuracyRange, mergin);
		    time = values.meanTime;
		    break;
		  } catch(NGT::Exception &err) {
		    if (err.getMessage().find("Error!! Epsilon") != std::string::npos &&
			err.getMessage().find("is too large") != std::string::npos) {
		      cerr << "Warning: Cannot adjust the rate of edge size." << err.what() << endl;
		      cerr << "Try again with the next rate" << endl;
		      NGTThrowException("**Retry**");
		    }
		    if (mergin > 0.4) {
		      cerr << "Error: Cannot adjust the rate even for the widest mergin " << mergin << ". " << err.what();
		      NGTThrowException("**Retry**"); 
		    } else {
		      cerr << "Warning: Cannot adjust the rate of edge size for mergin " << mergin << ". " << err.what() << endl;
		      cerr << "Try again for the next mergin." << endl;
		      mergin += 0.05;
		    }
		  }
		}
		times.insert(std::make_pair(rate, time));
		cerr << "adjustRateSearchEdgeSize::rate=" << rate << ", query time=" << time << endl;
	      } else {
		time = times.at(rate);
	      }
	      if (prevTime <= time) {
		if (rateStep == minimumStep) {
		  return std::make_pair(minRate, minTime);
		} else {
		  rateStart = static_cast<int>(minRate) - static_cast<int>(rateStep) < static_cast<int>(rateStart) ? rateStart : minRate - rateStep;
		  break;
		}
	      }
	      prevTime = time;
	      if (time < minTime) {
		minTime = time;
		minRate = rate;
	      }
	    }
	  }
	} catch(NGT::Exception &err) {
	  if (err.getMessage().find("**Retry**") != std::string::npos) {
	    rateStartInit += minimumStep;
	  } else {
	    throw err;
	  }
	}
      }
    }



    pair<size_t, size_t> adjustSearchEdgeSize(pair<float, float> baseAccuracyRange, pair<float, float> rateAccuracyRange, size_t querySize, double epsilon, float mergin = 0.2) {


      stringstream queries;
      stringstream gtStream;

      Command::SearchParameter searchParameter;
      NGT::GraphIndex &graphIndex = static_cast<GraphIndex&>(index.getIndex());
      NeighborhoodGraph::Property &prop = graphIndex.getGraphProperty();
      searchParameter.size = nOfResults;
      redirector.begin();
      try {
	cerr << "adjustSearchEdgeSize::Extract queries for GT..." << endl;
	extractQueries(querySize, queries);
	cerr << "adjustSearchEdgeSize::create GT..." << endl;
	createGroundTruth(index, epsilon, searchParameter, queries, gtStream);
      } catch (NGT::Exception &err) {
	cerr << "adjustSearchEdgeSize::Error!! Cannot adjust. " << err.what() << endl;
	redirector.end();
	return pair<size_t, size_t>(0, 0);
      }
      redirector.end();

      auto prevBase = pair<size_t, double>(0, 0);
      auto prevRate = pair<size_t, double>(0, 0);
      auto base = pair<size_t, double>(0, 0);
      auto rate = pair<size_t, double>(20, 0);

      map<pair<size_t, size_t>, double> history;
      redirector.begin();
      for(;;) {
	try {
	  prop.dynamicEdgeSizeRate = rate.first;
	  cerr << "adjustRateSearchEdgeSize::Base: rate=" << prop.dynamicEdgeSizeRate << endl;
	  prevBase = base;
	  base = adjustBaseSearchEdgeSize(queries, searchParameter, gtStream, baseAccuracyRange, mergin, prevBase.first);
	  cerr << "adjustRateSearchEdgeSize::Base: base=" << prevBase.first << "->" << base.first << ",rate=" << prevRate.first << "->" << rate.first << endl;
	  if (prevBase.first == base.first) {
	    break;
	  }
	  prop.dynamicEdgeSizeBase = base.first;
	  cerr << "adjustRateSearchEdgeSize::Rate: base=" << prop.dynamicEdgeSizeBase << endl;
	  prevRate = rate;
	  rate = adjustRateSearchEdgeSize(queries, searchParameter, gtStream, rateAccuracyRange, mergin, prevRate.first);
	  cerr << "adjustRateSearchEdgeSize::Rate base=" << prevBase.first << "->" << base.first << ",rate=" << prevRate.first << "->" << rate.first << endl;
	  if (prevRate.first == rate.first) {
	    break;
	  }
	  if (history.count(std::make_pair(base.first, rate.first)) != 0) {
	    cerr << "adjustRateSearchEdgeSize::Warning! Found an infinite loop." << endl;
	    double minTime = rate.second;
	    pair<size_t, size_t> min = std::make_pair(base.first, rate.first);
	    for (auto i = history.begin(); i != history.end(); ++i) {
	      double dc = (*i).second;
	      if (dc < minTime) {
		minTime = dc;
		min = (*i).first;
	      }
	    }
	    return min;
	  }
	  // store parameters here to prioritize high accuracy
	  history.insert(std::make_pair(std::make_pair(base.first, rate.first), rate.second));
	} catch (NGT::Exception &err) {
	  cerr << "adjustRateSearchEdgeSize::Error!! Cannot adjust. " << err.what() << endl;
	  redirector.end();
	  return pair<size_t, size_t>(0, 0);
	}
      }
      redirector.end();
      return std::make_pair(base.first, rate.first);
    }

    static void adjustSearchEdgeSize(Args &args)
    {
      const string usage = "Usage: ngt adjust-edge-size [-m mergin] [-e epsilon-for-ground-truth] [-q #-of-queries] [-n #-of-results] index";

      string indexName;
      try {
	indexName = args.get("#1");
      } catch (...) {
	cerr << "ngt: Error: DB is not specified" << endl;
	cerr << usage << endl;
	return;
      }

      pair<float, float> baseAccuracyRange = pair<float, float>(0.30, 0.50);
      pair<float, float> rateAccuracyRange = pair<float, float>(0.80, 0.90);

      string opt = args.getString("A", "");
      if (opt.size() != 0) {
	vector<string> tokens;
	NGT::Common::tokenize(opt, tokens, ":");
	if (tokens.size() >= 1) { baseAccuracyRange.first = NGT::Common::strtod(tokens[0]); }
	if (tokens.size() >= 2) { baseAccuracyRange.second = NGT::Common::strtod(tokens[1]); }
	if (tokens.size() >= 3) { rateAccuracyRange.first = NGT::Common::strtod(tokens[2]); }
	if (tokens.size() >= 4) { rateAccuracyRange.second = NGT::Common::strtod(tokens[3]); }
      }

      double mergin = args.getf("m", 0.2);
      double epsilon = args.getf("e", 0.1);
      size_t querySize = args.getl("q", 100);
      size_t nOfResults = args.getl("n", 10);

      cerr << "adjustRateSearchEdgeSize::range= " << baseAccuracyRange.first << "-" << baseAccuracyRange.second 
	   << "," << rateAccuracyRange.first << "-" << rateAccuracyRange.second << endl;
      cerr << "adjustRateSearchEdgeSize::# of queries=" << querySize << endl;

      NGT::Index	index(indexName);

      Optimizer		optimizer(index, nOfResults);
      try {
	auto v = optimizer.adjustSearchEdgeSize(baseAccuracyRange, rateAccuracyRange, querySize, epsilon, mergin);
	NGT::GraphIndex &graphIndex = static_cast<GraphIndex&>(index.getIndex());
	NeighborhoodGraph::Property &prop = graphIndex.getGraphProperty();
	if (v.first > 0) {
	  prop.dynamicEdgeSizeBase = v.first;
	}
	if (v.second > 0) {
	  prop.dynamicEdgeSizeRate = v.second;
	}
	if (prop.dynamicEdgeSizeRate > 0 && prop.dynamicEdgeSizeBase > 0) {
	  graphIndex.saveProperty(indexName);
	}
      } catch (NGT::Exception &err) {
	cerr << "adjustRateSearchEdgeSize::Error!! Cannot adjust. " << err.what() << endl;
	return;
      }
    }


    void outputObject(ostream &os, size_t id1, size_t id2, NGT::Property &prop) {
      switch (prop.objectType) {
      case NGT::ObjectSpace::ObjectType::Uint8:
	{
	  auto *obj1 = static_cast<uint8_t*>(index.getObjectSpace().getObject(id1));
	  auto *obj2 = static_cast<uint8_t*>(index.getObjectSpace().getObject(id2));
	  for (int i = 0; i < prop.dimension; i++) {
	    int d = (*obj1++ + *obj2++) / 2;
	    os << d;
	    if (i + 1 != prop.dimension) {
	      os << "\t";
	    }
	  }
	  os << endl;
	}
	break;
      default:
      case NGT::ObjectSpace::ObjectType::Float:
	{
	  auto *obj1 = static_cast<float*>(index.getObjectSpace().getObject(id1));
	  auto *obj2 = static_cast<float*>(index.getObjectSpace().getObject(id2));
	  for (int i = 0; i < prop.dimension; i++) {
	    os << (*obj1++ + *obj2++) / 2.0F;
	    if (i + 1 != prop.dimension) {
	      os << "\t";
	    }
	  }
	  os << endl;
	}
	break;
      }
    }

    void
      extractQueries(size_t nqueries, ostream &os, bool similarObject = false) {

      NGT::Property prop;
      index.getProperty(prop);

      size_t osize = index.getObjectRepositorySize();
      size_t interval = osize / nqueries;
      size_t count = 0;
      for (size_t id1 = 1; id1 < osize && count < nqueries; id1 += interval, count++) {
	size_t oft = 0;
	while (index.getObjectSpace().getRepository().isEmpty(id1 + oft)) {
	  oft++;
	  if (id1 + oft >= osize) {
	    stringstream msg;
	    msg << "Too many empty entries to extract.";
	    NGTThrowException(msg);
	  }
	}
	if (similarObject) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	  NGT::Object *query = index.getObjectSpace().allocateObject(*index.getObjectSpace().getRepository().get(id1 + oft));
#else
	  NGT::Object *query = index.getObjectSpace().getRepository().get(id1 + oft);
#endif
	  NGT::SearchContainer sc(*query);
	  NGT::ObjectDistances results;
	  sc.setResults(&results);
	  sc.setSize(nOfResults);
	  index.search(sc);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	  index.getObjectSpace().deleteObject(query);
#endif
	  if (results.size() < 2) {
	    stringstream msg;
	    msg << "Cannot get even two results for queries.";
	    NGTThrowException(msg);
	  }
	  size_t id2 = 1;
	  for (size_t i = 1; i < results.size(); i++) {
	    if (results[i].distance > 0.0) {
	      id2 = results[i].id;
	      break;
	    }
	  }
	  outputObject(os, id1 + oft, id2, prop);
	} else {
	  size_t id2 = id1 + oft + 1;
	  while (index.getObjectSpace().getRepository().isEmpty(id2)) {
	    id2++;
	    if (id2 >= osize) {
	      stringstream msg;
	      msg << "Too many empty entries to extract.";
	      NGTThrowException(msg);
	    }
	  }
	  outputObject(os, id1 + oft, id2, prop);
	}
      }
      assert(count == nqueries);
    
    }


    static void
      extractQueries(Args &args)
    {
      const string usage = "Usage: ngt eval-query -n #-of-queries index";

      string indexName;
      try {
	indexName = args.get("#1");
      } catch (...) {
	cerr << "ngt: Error: DB is not specified" << endl;
	cerr << usage << endl;
	return;
      }
      size_t nqueries = args.getl("n", 1000);

      NGT::Index	index(indexName);
      NGT::Optimizer	optimizer(index);
      optimizer.extractQueries(nqueries, cout);
    }

    static void createGroundTruth(NGT::Index &index, double epsilon, Command::SearchParameter &searchParameter, stringstream &queries, stringstream &gtStream){
      queries.clear();
      queries.seekg(0, ios_base::beg);

      Args args;
      args.insert(pair<string, string>("#1", "dummy"));
      args.insert(pair<string, string>("#2", "dummy"));
      searchParameter.parse(args);
      searchParameter.outputMode = 'e';
      searchParameter.edgeSize = -1;
      searchParameter.beginOfEpsilon = searchParameter.endOfEpsilon = epsilon;

      NGT::Command::search(index, searchParameter, queries, gtStream);
    }

    static int 
      calculateMeanValues(vector<MeasuredValue> &accuracies, double accuracyRangeFrom, double accuracyRangeTo, 
			  size_t &size, double &meanDistanceCount, double &meanVisitCount, double &meanTime) {
      int stat = 0;
      size = 0;
      if (accuracies.front().meanAccuracy > accuracyRangeFrom) {
	stat = 0x1;
      }
      if (accuracies.back().meanAccuracy < accuracyRangeTo) {
	stat |= 0x2;
      }
      if (stat != 0) {
	return stat;
      }
      vector<MeasuredValue> acc;
      acc = accuracies;
      for (auto start = acc.rbegin(); start != acc.rend(); ++start) {
	if ((*start).meanAccuracy <= accuracyRangeFrom) {
	  ++start;
	  acc.erase(acc.begin(), start.base());
	  break;
	}
      }
      for (auto end = acc.begin(); end != acc.end(); ++end) {
	if ((*end).meanAccuracy >= accuracyRangeTo) {
	  end++;
	  acc.erase(end, acc.end());
	  break;
	}
      }
      vector<pair<double, double>> distance;
      vector<pair<double, double>> visit;
      vector<pair<double, double>> time;
      for (auto i = acc.begin(); i != acc.end(); ++i) {
#ifdef NGT_LOG_BASED_OPTIMIZATION
	if ((*i).meanDistanceCount > 0.0) {
	  (*i).meanDistanceCount = log10((*i).meanDistanceCount);
	}
	if ((*i).meanVisitCount > 0.0) {
	  (*i).meanVisitCount = log10((*i).meanVisitCount);
	}
#endif
	distance.push_back(make_pair((*i).meanDistanceCount, (*i).meanAccuracy));
	visit.push_back(make_pair((*i).meanVisitCount, (*i).meanAccuracy));
	time.push_back(make_pair((*i).meanTime, (*i).meanAccuracy));
      }
      {
	size_t last = distance.size() - 1;
	double xfrom = (distance[1].second * distance[0].first - distance[0].second * distance[1].first + 
			accuracyRangeFrom * (distance[1].first - distance[0].first)) / 
	  (distance[1].second - distance[0].second);
	double xto = (distance[last].second * distance[last - 1].first - distance[last - 1].second * distance[last].first + 
		      accuracyRangeTo * (distance[last].first - distance[last - 1].first)) / 
	  (distance[last].second - distance[last - 1].second);
	distance[0].first = xfrom;
	distance[0].second = accuracyRangeFrom;
	distance[last].first = xto;
	distance[last].second = accuracyRangeTo;
	double area = 0.0;
	for (size_t i = 0; i < distance.size() - 1; ++i) {
	  area += ((distance[i].first + distance[i + 1].first) * (distance[i + 1].second - distance[i].second)) / 2.0;
	}
	meanDistanceCount = area / (distance[last].second - distance[0].second);
      }
      {
	size_t last = visit.size() - 1;
	double xfrom = (visit[1].second * visit[0].first - visit[0].second * visit[1].first + 
			accuracyRangeFrom * (visit[1].first - visit[0].first)) / 
	  (visit[1].second - visit[0].second);
	double xto = (visit[last].second * visit[last - 1].first - visit[last - 1].second * visit[last].first + 
		      accuracyRangeTo * (visit[last].first - visit[last - 1].first)) / 
	  (visit[last].second - visit[last - 1].second);
	visit[0].first = xfrom;
	visit[0].second = accuracyRangeFrom;
	visit[last].first = xto;
	visit[last].second = accuracyRangeTo;
	double area = 0.0;
	for (size_t i = 0; i < visit.size() - 1; ++i) {
	  area += ((visit[i].first + visit[i + 1].first) * (visit[i + 1].second - visit[i].second)) / 2.0;
	}
	meanVisitCount = area / (visit[last].second - visit[0].second);
      }
      {
	size_t last = time.size() - 1;
	double xfrom = (time[1].second * time[0].first - time[0].second * time[1].first + 
			accuracyRangeFrom * (time[1].first - time[0].first)) / 
	  (time[1].second - time[0].second);
	double xto = (time[last].second * time[last - 1].first - time[last - 1].second * time[last].first + 
		      accuracyRangeTo * (time[last].first - time[last - 1].first)) / 
	  (time[last].second - time[last - 1].second);
	time[0].first = xfrom;
	time[0].second = accuracyRangeFrom;
	time[last].first = xto;
	time[last].second = accuracyRangeTo;
	double area = 0.0;
	for (size_t i = 0; i < time.size() - 1; ++i) {
	  area += ((time[i].first + time[i + 1].first) * (time[i + 1].second - time[i].second)) / 2.0;
	}
	meanTime = area / (time[last].second - time[0].second);
      }
      assert(distance.size() == time.size());
      size = distance.size();
      return 0;
    }

    static void evaluate(Args &args)
    {
      const string usage = "Usage: ngt eval [-n number-of-results] [-m mode(r=recall)] [-g ground-truth-size] [-o output-mode] ground-truth search-result\n"
	"   Make a ground truth list (linear search): \n"
	"       ngt search -i s -n 20 -o e index query.list > ground-truth.list";

      string gtFile, resultFile;
      try {
	gtFile = args.get("#1");
      } catch (...) {
	cerr << "ground truth is not specified." << endl;
	cerr << usage << endl;
	return;
      }
      try {
	resultFile = args.get("#2");
      } catch (...) {
	cerr << "result file is not specified." << endl;
	cerr << usage << endl;
	return;
      }

      size_t resultSize = args.getl("n", 0);
      if (resultSize != 0) {
	cerr << "The specified number of results=" << resultSize << endl;
      }

      size_t groundTruthSize = args.getl("g", 0);

      bool recall = false;
      if (args.getChar("m", '-') == 'r') {
	cerr << "Recall" << endl;
	recall = true;
      }
      char omode = args.getChar("o", '-');
    
      ifstream	resultStream(resultFile);
      if (!resultStream) {
	cerr << "Cannot open the specified target file. " << resultFile << endl;
	cerr << usage << endl;
	return;
      }

      ifstream	gtStream(gtFile);
      if (!gtStream) {
	cerr << "Cannot open the specified GT file. " << gtFile << endl;
	cerr << usage << endl;
	return;
      }

      vector<MeasuredValue> accuracies;
      string type;
      size_t actualResultSize = 0;
      evaluate(gtStream, resultStream, accuracies, type, actualResultSize, resultSize, groundTruthSize, recall);

      cout << "# # of evaluated resultant objects per query=" << actualResultSize << endl;
      if (recall) {
	cout << "# " << type << "\t# of Queries\tRecall\t";
      } else {
	cout << "# " << type << "\t# of Queries\tPrecision\t";
      }
      if (omode == 'd') {
	cout << "# of computations\t# of visted nodes" << endl;
	for (auto it = accuracies.begin(); it != accuracies.end(); ++it) {
	  cout << (*it).keyValue << "\t" << (*it).totalCount << "\t" << (*it).meanAccuracy << "\t" 
	       << (*it).meanDistanceCount << "\t" << (*it).meanVisitCount << endl;
	}
      } else {
	cout << "Time(msec)\t# of computations\t# of visted nodes" << endl;
	for (auto it = accuracies.begin(); it != accuracies.end(); ++it) {
	  cout << (*it).keyValue << "\t" << (*it).totalCount << "\t" << (*it).meanAccuracy << "\t" << (*it).meanTime << "\t" 
	       << (*it).meanDistanceCount << "\t" << (*it).meanVisitCount << endl;
	}
      }

    }

    NGT::Index &index;
    size_t nOfResults;
    StdOstreamRedirector redirector;
  };

}; // NGT


