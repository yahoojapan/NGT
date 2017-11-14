//
// Copyright (C) 2015-2017 Yahoo Japan Corporation
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

#include	"NGT/Index.h"


namespace NGT {

class Args : public map<string, string>
{
public:
  Args(int argc, char **argv):
    option("a:b:c:d:e:f:g:hi:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:"
	   "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z:")
  {
    int opt;
    while ((opt = getopt(argc, argv, option)) != -1) {
      if ((char)opt == 'h') {
	string str;
	str.append(1, (char)opt);
	insert(pair<string, string>(str, ""));
	continue;
      }
      string str;
      str.append(1, (char)opt);
      insert(pair<string, string>(str, string(optarg)));
    }
    for (int i = 0; optind < argc; optind++, i++) {
      stringstream ss;
      ss << "#" << i;
      insert(pair<string, string>(ss.str(), string(argv[optind])));
    }
  }
  string &find(const char *s) { return get(s); }
  char getChar(const char *s, char v) {
    try {
      return get(s)[0];
    } catch (...) {
      return v;
    }
  }
  string getString(const char *s, const char *v) {
    try {
      return get(s);
    } catch (...) {
      return v;
    }
  }
  string &get(const char *s) {
    Args::iterator ai;
    ai = map<string, string>::find(string(s));
    if (ai == this->end()) {
      stringstream msg;
      msg << s << ": Not specified" << endl;
      NGTThrowException(msg.str());
    }
    return ai->second;
  }
  long getl(const char *s, long v) {
    char *e;
    long val;
    try {
      val = strtol(get(s).c_str(), &e, 10);
    } catch (...) {
      return v;
    }
    if (*e != 0) {
      stringstream msg;
      msg << "ARGS::getl: Illegal string. Option=-" << s << " Specified value=" << get(s) 
	  << " Illegal string=" << e << endl;
      NGTThrowException(msg.str());
    }
    return val;
  }
  float getf(const char *s, float v) {
    char *e;
    float val;
    try {
      val = strtof(get(s).c_str(), &e);
    } catch (...) {
      return v;
    }
    if (*e != 0) {
      stringstream msg;
      msg << "ARGS::getf: Illegal string. Option=-" << s << " Specified value=" << get(s) 
	  << " Illegal string=" << e << endl;
      NGTThrowException(msg.str());
    }
    return val;
  }
  const char *option;
};


class Command {
public:
  Command():debugLevel(0) {}

  void 
  create(Args &args)
  {
    const string usage = "Usage: ngt create "
      "-d dimension [-p #-of-thread] [-i index-type(t|g)] [-g graph-type(a|k|b)] "
      "[-t truncation-edge-limit] [-E edge-size] [-S edge-size-for-search] [-L edge-size-limit] "
      "[-e epsilon] [-o object-type(f|c)] [-D distance-function] [-n data-size] "
      "index(output) data.tsv(input)";
    string database;
    try {
      database = args.get("#1");
    } catch (...) {
      cerr << "ngt: Error: DB is not specified." << endl;
      cerr << usage << endl;
      return;
    }
    string data;
    try {
      data = args.get("#2");
    } catch (...) {
      cerr << "ngt: Error: Data is not specified." << endl;
    }

    NGT::Property property;

    property.edgeSizeForCreation = args.getl("E", 10);
    property.edgeSizeForSearch = args.getl("S", 40);
    property.batchSizeForCreation = args.getl("b", 200);
    property.insertionRadiusCoefficient = args.getf("e", 0.1) + 1.0;
    property.truncationThreshold = args.getl("t", 0);
    property.dimension = args.getl("d", 0);
    property.threadPoolSize = args.getl("p", 24);

    if (property.dimension <= 0) {
      cerr << "ngt: Error: Specify greater than 0 for # of your data dimension by a parameter -d." << endl;
      cerr << usage << endl;
      return;
    }

    char graphType = args.getChar("g", 'a');
    switch(graphType) {
    case 'a': property.graphType = NGT::Property::GraphType::GraphTypeANNG; break;
    case 'k': property.graphType = NGT::Property::GraphType::GraphTypeKNNG; break;
    case 'b': property.graphType = NGT::Property::GraphType::GraphTypeBKNNG; break;
    case 'd': property.graphType = NGT::Property::GraphType::GraphTypeDNNG; break;
    default:
      cerr << "ngt: Error: Invalid graph type. " << graphType << endl;
      cerr << usage << endl;
      return;
    }    

    char seedType = args.getChar("s", 'r');
    switch(seedType) {
    case 'f': property.seedType = NGT::Property::SeedType::SeedTypeFixedNodes; break;
    case '1': property.seedType = NGT::Property::SeedType::SeedTypeFirstNode; break;
    default:
    case 'r': property.seedType = NGT::Property::SeedType::SeedTypeRandomNodes; break;
    }

    char objectType = args.getChar("o", 'f');
    char distanceType = args.getChar("D", '2');

    size_t dataSize = args.getl("n", 0);
    char indexType = args.getChar("i", 't');

    if (debugLevel >= 1) {
      cerr << "edgeSizeForCreation=" << property.edgeSizeForCreation << endl;
      cerr << "edgeSizeForSearch=" << property.edgeSizeForSearch << endl;
      cerr << "edgeSizeLimit=" << property.edgeSizeLimitForCreation << endl;
      cerr << "batch size=" << property.batchSizeForCreation << endl;
      cerr << "graphType=" << property.graphType << endl;
      cerr << "epsilon=" << property.insertionRadiusCoefficient - 1.0 << endl;
      cerr << "thread size=" << property.threadPoolSize << endl;
      cerr << "dimension=" << property.dimension << endl;
      cerr << "indexType=" << indexType << endl;
    }

    switch (objectType) {
    case 'f': 
      property.objectType = NGT::Index::Property::ObjectType::Float;
      break;
    case 'c':
      property.objectType = NGT::Index::Property::ObjectType::Uint8;
      break;
    default:
      cerr << "ngt: Error: Invalid object type. " << objectType << endl;
      cerr << usage << endl;
      return;
    }

    switch (distanceType) {
    case '1': 
      property.distanceType = NGT::Index::Property::DistanceType::DistanceTypeL1;
      break;
    case '2':
      property.distanceType = NGT::Index::Property::DistanceType::DistanceTypeL2;
      break;
    case 'a':
      property.distanceType = NGT::Index::Property::DistanceType::DistanceTypeAngle;
      break;
    case 'h':
      property.distanceType = NGT::Index::Property::DistanceType::DistanceTypeHamming;
      break;
    case 'c':
      property.distanceType = NGT::Index::Property::DistanceType::DistanceTypeCosine;
      break;
    default:
      cerr << "ngt: Error: Invalid distance type. " << distanceType << endl;
      cerr << usage << endl;
      return;
    }

    switch (indexType) {
    case 't':
      NGT::Index::createGraphAndTree(database, property, data, dataSize);
      break;
    case 'g':
      NGT::Index::createGraph(database, property, data, dataSize);	
      break;
    }
  }

  void 
  append(Args &args)
  {
    const string usage = "Usage: ngt append [-p #-of-thread] [-d dimension] [-n data-size] "
      "index(output) data.tsv(input)";
    string database;
    try {
      database = args.get("#1");
    } catch (...) {
      cerr << "ngt: Error: DB is not specified." << endl;
      cerr << usage << endl;
      return;
    }
    string data;
    try {
      data = args.get("#2");
    } catch (...) {
      cerr << "ngt: Error: Data is not specified." << endl;
      cerr << usage << endl;
      return;
    }

    int threadSize = args.getl("p", 50);
    size_t dimension = args.getl("d", 0);
    size_t dataSize = args.getl("n", 0);

    if (debugLevel >= 1) {
      cerr << "thread size=" << threadSize << endl;
      cerr << "dimension=" << dimension << endl;
    }

    try {
      NGT::Index::append(database, data, threadSize, dataSize);	
    } catch (NGT::Exception &err) {
      cerr << "ngt: Error " << err.what() << endl;
      cerr << usage << endl;
    } catch (...) {
      cerr << "ngt: Error" << endl;
      cerr << usage << endl;
    }
  }

  class SearchParameter {
  public:
    SearchParameter() {}
    SearchParameter(Args &args) {
      try {
	query = args.get("#2");
      } catch (...) {
	NGTThrowException("ngt: Error: Query is not specified");
      }
      querySize = args.getl("Q", 0);
      indexType	= args.getChar("i", 't');
      size		= args.getl("n", 20);
      // edgeSize
      // -1 : using the size which was specified at the index creation.
      //  0 : no limitation for the edge size.
      edgeSize	= args.getl("E", -1);
      outputMode	= args.getString("o", "-");
      radius		= args.getf("r", FLT_MAX);
      trial		= args.getl("t", 1);
      {
	beginOfEpsilon = endOfEpsilon = stepOfEpsilon = 0.1;
	string epsilon = args.getString("e", "0.1");
	vector<string> tokens;
	NGT::Common::tokenize(epsilon, tokens, ":");
	if (tokens.size() >= 1) { beginOfEpsilon = endOfEpsilon = NGT::Common::strtod(tokens[0]); }
	if (tokens.size() >= 2) { endOfEpsilon = NGT::Common::strtod(tokens[1]); }
	if (tokens.size() >= 3) { stepOfEpsilon = NGT::Common::strtod(tokens[2]); }
	step = 0;
	if (tokens.size() >= 4) { step = NGT::Common::strtol(tokens[3]); }
      }
    }
    string	query;
    size_t	querySize;
    char	indexType;
    int		size;
    long	edgeSize;
    string	outputMode;
    float	radius;
    float	beginOfEpsilon;
    float	endOfEpsilon;
    float	stepOfEpsilon;
    size_t	step;
    size_t	trial;
  };

  static void
  search(NGT::Index &index, SearchParameter &searchParameter, ostream &stream)
  {

    if (searchParameter.outputMode[0] == 'e') { 
      stream << "# Beginning of Evaluation" << endl; 
    }

    ifstream		is(searchParameter.query);
    if (!is) {
      cerr << "Cannot open the specified file. " << searchParameter.query << endl;
      return;
    }
    string line;
    double totalTime	= 0;
    size_t queryCount	= 0;
    while(getline(is, line)) {
      if (searchParameter.querySize > 0 && queryCount >= searchParameter.querySize) {
	break;
      }
      NGT::Object *object = index.allocateObject(line, " \t");
      queryCount++;
      size_t step = searchParameter.step == 0 ? UINT_MAX : searchParameter.step;
      for (size_t n = 0; n <= step; n++) {
	NGT::SearchContainer sc(*object);
	double epsilon;
	if (searchParameter.step != 0) {
	  epsilon = searchParameter.beginOfEpsilon + (searchParameter.endOfEpsilon - searchParameter.beginOfEpsilon) * n / step; 
	} else {
	  epsilon = searchParameter.beginOfEpsilon + searchParameter.stepOfEpsilon * n;
	  if (epsilon > searchParameter.endOfEpsilon) {
	    break;
	  }
	}
	NGT::ObjectDistances objects;
	sc.setResults(&objects);
	sc.setSize(searchParameter.size);
	sc.setRadius(searchParameter.radius);
	sc.setEpsilon(epsilon);
	sc.setEdgeSize(searchParameter.edgeSize);
	NGT::Timer timer;
	try {
	  if (searchParameter.outputMode[0] == 'e') {
	    double time = 0.0;
	    uint64_t ntime = 0;
	    double minTime = DBL_MAX;
	    size_t trial = searchParameter.trial <= 1 ? 2 : searchParameter.trial;
	    for (size_t t = 0; t < trial; t++) {
	      switch (searchParameter.indexType) {
	      case 't': timer.start(); index.search(sc); timer.stop(); break;
	      case 'g': timer.start(); index.searchUsingOnlyGraph(sc); timer.stop(); break;
	      case 's': timer.start(); index.linearSearch(sc); timer.stop(); break;
	      }
	      if (minTime > timer.time) {
		minTime = timer.time;
	      }
	      time += timer.time;
	      ntime += timer.ntime;
	    }
	    time /= (double)searchParameter.trial;
	    ntime /= searchParameter.trial;
	    timer.time = minTime;
	    timer.ntime = ntime;
	  } else {
	    switch (searchParameter.indexType) {
	    case 't': timer.start(); index.search(sc); timer.stop(); break;
	    case 'g': timer.start(); index.searchUsingOnlyGraph(sc); timer.stop(); break;
	    case 's': timer.start(); index.linearSearch(sc); timer.stop(); break;
	    }
	  }
	} catch (NGT::Exception &err) {
	  throw err;
	}
	totalTime += timer.time;
	if (searchParameter.outputMode[0] == 'e') {
	  stream << "# Query No.=" << queryCount << endl;
	  stream << "# Query=" << line.substr(0, 20) + " ..." << endl;
	  stream << "# Index Type=" << searchParameter.indexType << endl;
	  stream << "# Size=" << searchParameter.size << endl;
	  stream << "# Radius=" << searchParameter.radius << endl;
	  stream << "# Epsilon=" << epsilon << endl;
	  stream << "# Query Time (msec)=" << timer.time * 1000.0 << endl;
	  stream << "# Distance Computation=" << sc.distanceComputationCount << endl;
	} else {
	  stream << "Query No." << queryCount << endl;
	  stream << "Rank\tID\tDistance" << endl;
	}
	for (size_t i = 0; i < objects.size(); i++) {
	  stream << i + 1 << "\t" << objects[i].id << "\t";
	  stream << objects[i].distance << endl;
	}
	if (searchParameter.outputMode[0] == 'e') {
	  stream << "# End of Search" << endl;
	} else {
	  stream << "Query Time= " << timer.time << " (sec), " << timer.time * 1000.0 << " (msec)" << endl;
	}
      } // for
      index.deleteObject(object);
      if (searchParameter.outputMode[0] == 'e') {
	stream << "# End of Query" << endl;
      }
    } // while
    if (searchParameter.outputMode[0] == 'e') {
      stream << "# Average Query Time (msec)=" << totalTime * 1000.0 / (double)queryCount << endl;
      stream << "# Number of queries=" << queryCount << endl;
      stream << "# End of Evaluation" << endl;

      if (searchParameter.outputMode == "e+") {
	// show graph information
	size_t esize = searchParameter.edgeSize;
	long double distance = 0.0;
	size_t numberOfNodes = 0;
	size_t numberOfEdges = 0;

	NGT::GraphIndex	&graph = (NGT::GraphIndex&)index.getIndex();
	for (size_t id = 1; id < graph.repository.size(); id++) {
	  NGT::GraphNode *node = 0;
	  try {
	    node = graph.getNode(id);
	  } catch(NGT::Exception &err) {
	    cerr << "Graph::search: Warning. Cannot get the node. ID=" << id << ":" << err.what() << " If the node was removed, no problem." << endl;
	    continue;
	  }
	  numberOfNodes++;
	  if (numberOfNodes % 1000000 == 0) {
	    cerr << "Processed " << numberOfNodes << endl;
	  }
	  for (size_t i = 0; i < node->size(); i++) {
	    if (esize != 0 && i >= esize) {
	      break;
	    }
	    numberOfEdges++;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    distance += (*node).at(i, graph.repository.allocator).distance;
#else
	    distance += (*node)[i].distance;
#endif
	  }
	}

	stream << "# # of nodes=" << numberOfNodes << endl;
	stream << "# # of edges=" << numberOfEdges << endl;
	stream << "# Average number of edges=" << (double)numberOfEdges / (double)numberOfNodes << endl;
	stream << "# Average distance of edges=" << setprecision(10) << distance / (double)numberOfEdges << endl;
      }
    } else {
      stream << "Average Query Time= " << totalTime / (double)queryCount  << " (sec), " 
	   << totalTime * 1000.0 / (double)queryCount << " (msec), (" 
	   << totalTime << "/" << queryCount << ")" << endl;
    }
  }


  void
  search(Args &args) {
    const string usage = "Usage: ngt search [-i g|t|s] [-n result-size] [-e epsilon] [-E edge-size] "
      "[-o output-mode] index(input) query.tsv(input)";

    string database;
    try {
      database = args.get("#1");
    } catch (...) {
      cerr << "ngt: Error: DB is not specified" << endl;
      cerr << usage << endl;
      return;
    }

    SearchParameter searchParameter(args);

    if (debugLevel >= 1) {
      cerr << "indexType=" << searchParameter.indexType << endl;
      cerr << "size=" << searchParameter.size << endl;
      cerr << "edgeSize=" << searchParameter.edgeSize << endl;
      cerr << "epsilon=" << searchParameter.beginOfEpsilon << "<->" << searchParameter.endOfEpsilon << "," 
	   << searchParameter.stepOfEpsilon << endl;
    }

    try {
      NGT::Property property;
      property.clear();
      NGT::Index	index(database, property);
      search(index, searchParameter, cout);
    } catch (NGT::Exception &err) {
      cerr << "ngt: Error " << err.what() << endl;
      cerr << usage << endl;
    } catch (...) {
      cerr << "ngt: Error" << endl;
      cerr << usage << endl;
    }

  }


  void
  remove(Args &args)
  {
    const string usage = "Usage: ngt remove [-d object-ID-type(f|d)] index(input) object-ID(input)";
    string database;
    try {
      database = args.get("#1");
    } catch (...) {
      cerr << "ngt: Error: DB is not specified" << endl;
      cerr << usage << endl;
      return;
    }
    try {
      args.get("#2");
    } catch (...) {
      cerr << "ngt: Error: ID is not specified" << endl;
      cerr << usage << endl;
      return;
    }
    char dataType = args.getChar("d", 'f');
    if (debugLevel >= 1) {
      cerr << "dataType=" << dataType << endl;
    }

    try {
      vector<NGT::ObjectID> objects;
      if (dataType == 'f') {
	string ids;
	try {
	  ids = args.get("#2");
	} catch (...) {
	  cerr << "ngt: Error: Data file is not specified" << endl;
	  cerr << usage << endl;
	  return;
	}
	ifstream is(ids);
	if (!is) {
	  cerr << "ngt: Error: Cannot open the specified file. " << ids << endl;
	  cerr << usage << endl;
	  return;
	}
	string line;
	int count = 0;
	while(getline(is, line)) {
	  count++;
	  vector<string> tokens;
	  NGT::Common::tokenize(line, tokens, "\t ");
	  if (tokens.size() == 0 || tokens[0].size() == 0) {
	    continue;
	  }
	  char *e;
	  size_t id;
	  try {
	    id = strtol(tokens[0].c_str(), &e, 10);
	    objects.push_back(id);
	  } catch (...) {
	    cerr << "Illegal data. " << tokens[0] << endl;
	  }
	  if (*e != 0) {
	    cerr << "Illegal data. " << e << endl;
	  }
	  cerr << "removed ID=" << id << endl;	
	}
      } else {
	size_t id = args.getl("#2", 0);
	cerr << "removed ID=" << id << endl;
	objects.push_back(id);
      }
      NGT::Index::remove(database, objects);
    } catch (NGT::Exception &err) {
      cerr << "ngt: Error " << err.what() << endl;
      cerr << usage << endl;
    } catch (...) {
      cerr << "ngt: Error" << endl;
      cerr << usage << endl;
    }
  }

  void
  exportIndex(Args &args)
  {
    const string usage = "Usage: ngt export index(input) export-file(output)";
    string database;
    try {
      database = args.get("#1");
    } catch (...) {
      cerr << "ngt: Error: DB is not specified" << endl;
      cerr << usage << endl;
      return;
    }
    string exportFile;
    try {
      exportFile = args.get("#2");
    } catch (...) {
      cerr << "ngt: Error: ID is not specified" << endl;
      cerr << usage << endl;
      return;
    }
    try {
      NGT::Index::exportIndex(database, exportFile);
    } catch (NGT::Exception &err) {
      cerr << "ngt: Error " << err.what() << endl;
      cerr << usage << endl;
    } catch (...) {
      cerr << "ngt: Error" << endl;
      cerr << usage << endl;
    }
  }

  void
  importIndex(Args &args)
  {
    const string usage = "Usage: ngt import index(output) import-file(input)";
    string database;
    try {
      database = args.get("#1");
    } catch (...) {
      cerr << "ngt: Error: DB is not specified" << endl;
      cerr << usage << endl;
      return;
    }
    string importFile;
    try {
      importFile = args.get("#2");
    } catch (...) {
      cerr << "ngt: Error: ID is not specified" << endl;
      cerr << usage << endl;
      return;
    }

    try {
      NGT::Index::importIndex(database, importFile);
    } catch (NGT::Exception &err) {
      cerr << "ngt: Error " << err.what() << endl;
      cerr << usage << endl;
    } catch (...) {
      cerr << "ngt: Error" << endl;
      cerr << usage << endl;
    }

  }

  void
  prune(Args &args)
  {
    const string usage = "Usage: ngt prune -e #-of-forcedly-pruned-edges -s #-of-selecively-pruned-edge";
    string indexName;
    try {
      indexName = args.get("#1");
    } catch (...) {
      cerr << "Index is not specified" << endl;
      cerr << usage << endl;
      return;
    }

    // the number of forcedly pruned edges
    size_t forcedlyPrunedEdgeSize	= args.getl("e", 0);
    // the number of selectively pruned edges
    size_t selectivelyPrunedEdgeSize	= args.getl("s", 0);

    cerr << "forcedly pruned edge size=" << forcedlyPrunedEdgeSize << endl;
    cerr << "selectively pruned edge size=" << selectivelyPrunedEdgeSize << endl;

    if (selectivelyPrunedEdgeSize == 0 && forcedlyPrunedEdgeSize == 0) {
      cerr << "prune: Error! Either of selective edge size or remaining edge size should be specified." << endl;
      cerr << usage << endl;
      return;
    }

    if (forcedlyPrunedEdgeSize != 0 && selectivelyPrunedEdgeSize != 0 && selectivelyPrunedEdgeSize >= forcedlyPrunedEdgeSize) {
      cerr << "prune: Error! selective edge size is less than remaining edge size." << endl;
      cerr << usage << endl;
      return;
    }

    NGT::Index	index(indexName);
    cerr << "loaded the input index." << endl;

    NGT::GraphIndex	&graph = (NGT::GraphIndex&)index.getIndex();

    for (size_t id = 1; id < graph.repository.size(); id++) {
      try {
	NGT::GraphNode &node = *graph.getNode(id);
	if (id % 1000000 == 0) {
	  cerr << "Processed " << id << endl;
	}
	if (forcedlyPrunedEdgeSize > 0 && node.size() >= forcedlyPrunedEdgeSize) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	  node.resize(forcedlyPrunedEdgeSize, graph.repository.allocator);
#else
	  node.resize(forcedlyPrunedEdgeSize);
#endif
	}
	if (selectivelyPrunedEdgeSize > 0 && node.size() >= selectivelyPrunedEdgeSize) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	  cerr << "not implemented" << endl;
	  abort();
#else
	  size_t rank = 0;
	  for (NGT::GraphNode::iterator i = node.begin(); i != node.end(); ++rank) {
	    if (rank >= selectivelyPrunedEdgeSize) {
	      bool found = false;
	      for (size_t t1 = 0; t1 < node.size() && found == false; ++t1) {
		if (t1 >= selectivelyPrunedEdgeSize) {
		  break;
		}
		if (rank == t1) {
		  continue;
		}
		NGT::GraphNode &node2 = *graph.getNode(node[t1].id);
		for (size_t t2 = 0; t2 < node2.size(); ++t2) {		
		  if (t2 >= selectivelyPrunedEdgeSize) {
		    break;
		  }
		  if (node2[t2].id == (*i).id) {
		    found = true;
		    break;
		  }
		} // for
	      } // for
	      if (found) {
		//remove
		i = node.erase(i);
		continue;
	      }
	    }
	    i++;
	  } // for
#endif
	}
	  
      } catch(NGT::Exception &err) {
	cerr << "Graph::search: Warning. Cannot get the node. ID=" << id << ":" << err.what() << endl;
	continue;
      }
    }

    graph.saveIndex(indexName);

  }



  void
  info(Args &args)
  {
    const string usage = "Usage: ngt info [-E #-of-edges] [-m h|e] index";

    string database;
    try {
      database = args.get("#1");
    } catch (...) {
      cerr << "ngt: Error: DB is not specified" << endl;
      cerr << usage << endl;
      return;
    }

    size_t edgeSize = args.getl("E", UINT_MAX);
    char mode = args.getChar("m", '-');

    try {
      NGT::GraphIndex	index(database);
      NGT::GraphIndex::showStatisticsOfGraph(index, mode, edgeSize);
    } catch (NGT::Exception &err) {
      cerr << "ngt: Error " << err.what() << endl;
      cerr << usage << endl;
    } catch (...) {
      cerr << "ngt: Error" << endl;
      cerr << usage << endl;
    }
  }


  //////////////////////////////////////////////////////////////

  void setDebugLevel(int level) { debugLevel = level; }
  int getDebugLevel() { return debugLevel; }

protected:
  int debugLevel;

};

}; // NGT
