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

#define NGT_LOG_BASED_OPTIMIZATION

namespace NGT {
  class Optimizer {
  public:
    Optimizer() {}
    ~Optimizer() {}

    class Accuracy {
    public:
      double	keyValue;
      size_t	totalCount;
      float	averageAccuracy;
      double	averageTime;
      double	averageDistanceCount;
      double	averageVisitCount;
    };

    static void search(NGT::Index &index, istream &gtStream, Command::SearchParameter &sp, vector<Accuracy> &acc) {
      ifstream		is(sp.query);
      if (!is) {
	stringstream msg;
	msg << "Cannot open the specified file. " << sp.query;
	NGTThrowException(msg);
      }

      search(index, gtStream, sp, acc);
    }

    static void search(NGT::Index &index, istream &queries, istream &gtStream, Command::SearchParameter &sp, vector<Accuracy> &acc) {
      sp.stepOfEpsilon = 1.0;
      stringstream resultStream;
      NGT::Command::search(index, sp, queries, resultStream);
      gtStream.clear();
      gtStream.seekg(0, ios_base::beg);
      resultStream.clear();
      resultStream.seekg(0, ios_base::beg);
      string type;
      size_t actualResultSize = 0;
      evaluate(gtStream, resultStream, acc, type, actualResultSize);
      assert(acc.size() == 1);
    }

    static void
      evaluate(istream &gtStream, istream &resultStream, vector<Accuracy> &accuracies, string &type, 
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
	    if (groundTruthSize == 0) {
	      loadGroundTruth(gtStream, gt, resultDataSize);
	    } else {
	      loadGroundTruth(gtStream, gt, groundTruthSize);
	    }
	    sumup(resultStream, queryNo, totalAccuracy, totalTime, totalDistanceCount, totalVisitCount, totalCount, 
		  gt, resultDataSize, type, recall);
	    queryNo++;
	  }
	}
      }

      accuracies.clear();
      for (auto it = totalAccuracy.begin(); it != totalAccuracy.end(); ++it) {
	Accuracy a;
	a.keyValue = (*it).first;
	a.totalCount = totalCount[a.keyValue];
	a.averageAccuracy = totalAccuracy[a.keyValue] / (double)totalCount[a.keyValue];
	a.averageTime = totalTime[a.keyValue] / (double)totalCount[a.keyValue];
	a.averageDistanceCount = (double)totalDistanceCount[a.keyValue] / (double)totalCount[a.keyValue];
	a.averageVisitCount = (double)totalVisitCount[a.keyValue] / (double)totalCount[a.keyValue];
	accuracies.push_back(a);
      }
    }

    static void
      loadGroundTruth(istream & gtf, unordered_set<size_t> & gt, size_t resultDataSize) {
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
		      bool recall)
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
		if (gt.count(id) != 0) {
		  relevantCount++;
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
					  Command::SearchParameter &sp, pair<float, float> interval, double mergin) 
    {
      double fromUnder = 0.0;
      double fromOver = 1.0;
      double toUnder = 0.0;
      double toOver = 1.0;
      float fromUnderEpsilon = -0.9;
      float fromOverEpsilon = -0.9;
      float toUnderEpsilon = -0.9;
      float toOverEpsilon = -0.9;

      float intervalFrom = interval.first;
      float intervalTo = interval.second;

      double range = intervalTo - intervalFrom;

      vector<Accuracy> acc;

      {
	float startEpsilon = -0.4;
	float epsilonStep = 0.02;
	size_t count;
	for (count = 0;; count++) {
	  float epsilon = startEpsilon + epsilonStep * count;
	  if (epsilon > 0.2) {
	    stringstream msg;
	    msg << "exploreEpsilonForAccuracy:" << endl;
	    msg << "Error!! Epsilon is too large. " << epsilon;
	    NGTThrowException(msg);
	  }
	  acc.clear();
	  sp.beginOfEpsilon = sp.endOfEpsilon = fromOverEpsilon = epsilon;
	  queries.clear();
	  queries.seekg(0, ios_base::beg);
	  search(index, queries, gtStream, sp, acc);
	  if (acc[0].averageAccuracy >= intervalFrom) {
	    break;
	  }
	}
	if (fromOverEpsilon == startEpsilon) {
	  stringstream msg;
	  msg << "exploreEpsilonForAccuracy:" << endl;
	  msg << "Error! startEpsilon should be reduced for the specified range.";
	  NGTThrowException(msg);
	}
	fromOver = acc[0].averageAccuracy;

	if (fromOver < intervalTo) {
	  startEpsilon = fromOverEpsilon;
	  for (count = 0;; count++) {
	    float epsilon = startEpsilon + epsilonStep * count;
	    sp.beginOfEpsilon = sp.endOfEpsilon = toOverEpsilon = epsilon;
	    if (epsilon > 0.2) {
	      stringstream msg;
	      msg << "exploreEpsilonForAccuracy:" << endl;
	      msg << "Error!! Epsilon is too large. " << epsilon;
	      NGTThrowException(msg);
	    }
	    acc.clear();
	    queries.clear();
	    queries.seekg(0, ios_base::beg);
	    search(index, queries, gtStream, sp, acc);
	    epsilon += epsilonStep;
	    if (acc[0].averageAccuracy >= intervalTo) {
	      break;
	    }
	  }
	  toOver = acc[0].averageAccuracy;
	} else {
	  toOver = fromOver;
	  toOverEpsilon = fromOverEpsilon;
	}
	if (fromOverEpsilon == toOverEpsilon) {
	  cerr << "Warning!! fromOverEpsilon equals toOverEpsilon  " << fromOverEpsilon << ". This might cause some problems." << endl;
	}
	fromUnderEpsilon = fromOverEpsilon - epsilonStep;
      }
      sp.beginOfEpsilon = sp.endOfEpsilon = fromUnderEpsilon;
      while (true) {
	acc.clear();
	queries.clear();
	queries.seekg(0, ios_base::beg);
	search(index, queries, gtStream, sp, acc);
	if (acc[0].averageAccuracy >= fromUnder && acc[0].averageAccuracy <= intervalFrom) {
	  fromUnder = acc[0].averageAccuracy;
	  fromUnderEpsilon = acc[0].keyValue;
	}
	if (acc[0].averageAccuracy <= fromOver && acc[0].averageAccuracy > intervalFrom) {
	  fromOver = acc[0].averageAccuracy;
	  fromOverEpsilon = acc[0].keyValue;
	}
	if (acc[0].averageAccuracy <= toOver && acc[0].averageAccuracy > intervalTo) {
	  toOver = acc[0].averageAccuracy;
	  toOverEpsilon = acc[0].keyValue;
	}
	if (acc[0].averageAccuracy >= toUnder && acc[0].averageAccuracy <= intervalTo) {
	  toUnder = acc[0].averageAccuracy;
	  toUnderEpsilon = acc[0].keyValue;
	}

	if (fromUnder < intervalFrom - range * mergin) {
	  if ((fromUnderEpsilon + fromOverEpsilon) / 2.0 == sp.beginOfEpsilon) {
	    stringstream msg;
	    msg << "exploreEpsilonForAccuracy:" << endl;
	    msg << "Error!! Not found proper under epsilon for mergin=" << mergin << " and the number of queries." << endl;
	    msg << "        Should increase mergin or the number of queries to get the proper epsilon.";
	    NGTThrowException(msg);
	  } else {
	    sp.beginOfEpsilon = sp.endOfEpsilon = (fromUnderEpsilon + fromOverEpsilon) / 2.0;
	  }
	} else if (toOver > intervalTo + range * mergin) {
	  if ((toUnderEpsilon + toOverEpsilon) / 2.0 == sp.beginOfEpsilon) {
	    stringstream msg;
	    msg << "exploreEpsilonForAccuracy:" << endl;
	    msg << "Error!! Not found proper over epsilon for mergin=" << mergin << " and the number of queries." << endl;
	    msg << "        Should increase mergin or the number of queries to get the proper epsilon.";
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

    static double measureDistance(NGT::Index &index, istream &queries, istream &gtStream, Command::SearchParameter &searchParameter, pair<float, float> interval, double mergin) {

      exploreEpsilonForAccuracy(index, queries, gtStream, searchParameter, interval, mergin);
    
      stringstream resultStream;
      queries.clear();
      queries.seekg(0, ios_base::beg);
      NGT::Command::search(index, searchParameter, queries, resultStream);
      gtStream.clear();
      gtStream.seekg(0, ios_base::beg);
      resultStream.clear();
      resultStream.seekg(0, ios_base::beg);
      string type;
      vector<Accuracy> accuracies;
      size_t actualResultSize = 0;
      evaluate(gtStream, resultStream, accuracies, type, actualResultSize);
      size_t size;
      double distanceCount;
      calculateAverageDistanceCount(accuracies, interval.first, interval.second, size, distanceCount);
      if (distanceCount == 0) {
	stringstream msg;
	msg << "measureDistance: Error! Distance count is zero.";
	NGTThrowException(msg);
      }
      return distanceCount;
    }

    static size_t adjustBaseSearchEdgeSize(NGT::Index &index, pair<float, float> interval, size_t querySize, double epsilon, float mergin = 0.2) {

      cerr << "adjustBaseSearchEdgeSize::Extract queries for GT..." << endl;
      stringstream queries;
      extractQueries(index, querySize, queries);

      queries.clear();
      queries.seekg(0, ios_base::beg);

      Args args;
      args.insert(pair<string, string>("#1", "dummy"));
      args.insert(pair<string, string>("#2", "dummy"));
      Command::SearchParameter searchParameter(args);
      searchParameter.outputMode = 'e';
      searchParameter.edgeSize = -2;
      searchParameter.beginOfEpsilon = searchParameter.endOfEpsilon = epsilon;

      cerr << "adjustBaseSearchEdgeSize::create GT..." << endl;
      stringstream gtStream;
      NGT::Command::search(index, searchParameter, queries, gtStream);

      while(true) {
	double prevDistanceComputation = INT_MAX;
	size_t prevEdgeBase = 0;
	size_t edgeSizeBaseStart = 10;
	size_t edgeSizeBaseStep = 10;
	cerr << "adjustBaseSearchEdgeSize::explore for the mergin " << mergin << "..." << endl;
	for (size_t edgeSizeBase = edgeSizeBaseStart; edgeSizeBase < 120; edgeSizeBase += edgeSizeBaseStep) {
	  searchParameter.step = 10;
	  NGT::GraphIndex &graphIndex = static_cast<GraphIndex&>(index.getIndex());
	  NeighborhoodGraph::Property &prop = graphIndex.getGraphProperty();
	  prop.dynamicEdgeSizeBase = edgeSizeBase;
	  try {
	    double distanceComputation = measureDistance(index, queries, gtStream, searchParameter, interval, mergin);
	    cerr << "adjustBaseSearchEdgeSize::Base edge size=" << edgeSizeBase << ", distance computation=" << distanceComputation << endl;
	    if (prevDistanceComputation < distanceComputation) {
	      return prevEdgeBase;
	    }
	    prevDistanceComputation = distanceComputation;
	    prevEdgeBase = edgeSizeBase;
	  } catch(NGT::Exception &err) {
	    if (mergin > 0.4) {
	      stringstream msg;
	      msg << "Error: Cannot adjust the base edge size even for the widest mergin " << mergin << ". " << err.what();
	      NGTThrowException(msg);
	    }
	    cerr << "Warning: Cannot adjust the base edge size for mergin " << mergin << ". " << err.what() << endl;
	    cerr << "Try again for the next mergin." << endl;
	    mergin += 0.05;
	    break;
	  }

	}
      }
    }

    static void adjustBaseSearchEdgeSize(Args &args)
    {
      const string usage = "Usage: ngt performance [-m mergin] [-e epsilon-for-ground-truth] [-n #-of queries] index";

      string indexName;
      try {
	indexName = args.get("#1");
      } catch (...) {
	cerr << "ngt: Error: DB is not specified" << endl;
	cerr << usage << endl;
	return;
      }

      float intervalFrom = 0.6;
      float intervalTo = 0.8;
      string opt = args.getString("A", "");
      if (opt.size() != 0) {
	vector<string> tokens;
	NGT::Common::tokenize(opt, tokens, ":");
	if (tokens.size() >= 1) { intervalFrom = NGT::Common::strtod(tokens[0]); }
	if (tokens.size() >= 2) { intervalTo = NGT::Common::strtod(tokens[1]); }
      }

      double mergin = args.getf("m", 0.2);
      double epsilon = args.getf("e", 0.1);
      size_t querySize = args.getl("n", 100);

      NGT::Index	index(indexName);

      size_t baseEdgeSize = 0;
      try {
	baseEdgeSize = adjustBaseSearchEdgeSize(index, pair<float, float>(intervalFrom, intervalTo), querySize, epsilon, mergin);
      } catch (NGT::Exception &err) {
	cerr << "adjustBaseSearchEdgeSize::Error!! Cannot adjust. " << err.what() << endl;
	return;
      }

      cerr << "adjustBaseSearchEdgeSize::The best base edge size=" << baseEdgeSize << endl;
      NGT::GraphIndex &graphIndex = static_cast<GraphIndex&>(index.getIndex());
      NeighborhoodGraph::Property &prop = graphIndex.getGraphProperty();
      prop.dynamicEdgeSizeBase = baseEdgeSize;
      graphIndex.saveProperty(indexName);
      cerr << "adjustBaseSearchEdgeSize::Set the base edge size to the index." << endl;

    }

    static void
      extractQueries(NGT::Index &index, size_t nqueries, ostream &os) {

      NGT::Property prop;
      index.getProperty(prop);

      size_t osize = index.getObjectRepositorySize();
      size_t interval = osize / nqueries;
      size_t count = 0;
      for (size_t id = 1; id < osize && count < nqueries; id += interval, count++) {
	size_t oft = 0;
	while (index.getObjectSpace().getRepository().isEmpty(id + oft)) {
	  oft++;
	  if (id + oft >= osize) {
	    cerr << "Too many empty entries to extract." << endl;
	    return;
	  }
	}
	switch (prop.objectType) {
	case NGT::ObjectSpace::ObjectType::Uint8:
	  {
	    auto *obj = static_cast<uint8_t*>(index.getObjectSpace().getObject(id + oft));
	    for (int i = 0; i < prop.dimension; i++) {
	      os << static_cast<int>(*obj++);
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
	    auto *obj = static_cast<float*>(index.getObjectSpace().getObject(id + oft));
	    for (int i = 0; i < prop.dimension; i++) {
	      os << *obj++;
	      if (i + 1 != prop.dimension) {
		os << "\t";
	      }
	    }
	    os << endl;
	  }
	  break;
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

      extractQueries(index, nqueries, cout);
    }

    static int 
      calculateAverageDistanceCount(vector<Accuracy> &accuracies, double intervalFrom, double intervalTo, size_t &size, 
				    double &averageDistanceCount) {
      int stat = 0;
      size = 0;
      averageDistanceCount = DBL_MAX;
      if (accuracies.front().averageAccuracy > intervalFrom) {
	stat = 0x1;
      }
      if (accuracies.back().averageAccuracy < intervalTo) {
	stat |= 0x2;
      }
      if (stat != 0) {
	return stat;
      }
      vector<Accuracy> acc;
      acc = accuracies;
      for (auto start = acc.rbegin(); start != acc.rend(); ++start) {
	if ((*start).averageAccuracy <= intervalFrom) {
	  ++start;
	  acc.erase(acc.begin(), start.base());
	  break;
	}
      }
      for (auto end = acc.begin(); end != acc.end(); ++end) {
	if ((*end).averageAccuracy >= intervalTo) {
	  end++;
	  acc.erase(end, acc.end());
	  break;
	}
      }
      vector<pair<double, double>> data;
      for (auto i = acc.begin(); i != acc.end(); ++i) {
#ifdef NGT_LOG_BASED_OPTIMIZATION
	if ((*i).averageDistanceCount > 0.0) {
	  (*i).averageDistanceCount = log10((*i).averageDistanceCount);
	}
	if ((*i).averageVisitCount > 0.0) {
	  (*i).averageVisitCount = log10((*i).averageVisitCount);
	}
#endif
	data.push_back(make_pair((*i).averageDistanceCount, (*i).averageAccuracy));
      }
      size_t last = data.size() - 1;
      double xfrom = (data[1].second * data[0].first - data[0].second * data[1].first + 
		      intervalFrom * (data[1].first - data[0].first)) / 
	(data[1].second - data[0].second);
      double xto = (data[last].second * data[last - 1].first - data[last - 1].second * data[last].first + 
		    intervalTo * (data[last].first - data[last - 1].first)) / 
	(data[last].second - data[last - 1].second);
      data[0].first = xfrom;
      data[0].second = intervalFrom;
      data[last].first = xto;
      data[last].second = intervalTo;
      double area = 0.0;
      for (size_t i = 0; i < data.size() - 1; ++i) {
	area += ((data[i].first + data[i + 1].first) * (data[i + 1].second - data[i].second)) / 2.0;
      }
      averageDistanceCount = area / (data[last].second - data[0].second);
      size = data.size();
      return 0;
    }

  };

}; // NGT


