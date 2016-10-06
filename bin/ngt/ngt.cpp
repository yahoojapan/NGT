
//
// Copyright (C) 2015-2016 Yahoo! JAPAN Research
//
// This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
// To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
//

#include	<sys/mman.h>

#include	"NGT/Index.h"


static float roundFloat(float f, int digit)
{
  return roundf(f * pow(10.0, digit)) / pow(10.0, digit);
}

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

class NGTCommand {
public:
  NGTCommand():debugLevel(0) {}

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

  void
  search(Args &args)
  {
    const string usage = "Usage: ngt search [-i g|t|s] [-n result-size] [-e epsilon] [-E edge-size] [-o output-mode] index(input) query.tsv(input)";
    string database;
    try {
      database = args.get("#1");
    } catch (...) {
      cerr << "ngt: Error: DB is not specified" << endl;
      cerr << usage << endl;
      return;
    }
    string query;
    try {
      query = args.get("#2");
    } catch (...) {
      cerr << "ngt: Error: Query is not specified" << endl;
      cerr << usage << endl;
      return;
    }

    char   indexType	= args.getChar("i", 't');
    int    size		= args.getl("n", 20);
    size_t edgeSize	= args.getl("E", 0);
    char   outputMode	= args.getChar("o", '-');
    float  radius	= args.getf("r", FLT_MAX);

    float beginOfEpsilon, endOfEpsilon, stepOfEpsilon;
    {
      beginOfEpsilon = endOfEpsilon = stepOfEpsilon = 0.1;
      string epsilon = args.getString("e", "0.1");
      vector<string> tokens;
      NGT::Common::tokenize(epsilon, tokens, ":");
      if (tokens.size() >= 1) { beginOfEpsilon = endOfEpsilon = NGT::Common::strtod(tokens[0]); }
      if (tokens.size() >= 2) { endOfEpsilon = NGT::Common::strtod(tokens[1]); }
      if (tokens.size() >= 3) { stepOfEpsilon = NGT::Common::strtod(tokens[2]); }
    }

    if (debugLevel >= 1) {
      cerr << "indexType=" << indexType << endl;
      cerr << "size=" << size << endl;
      cerr << "edgeSize=" << edgeSize << endl;
      cerr << "epsilon=" << beginOfEpsilon << "<->" << endOfEpsilon << "," << stepOfEpsilon << endl;
    }
    
    try {
      NGT::Property property;
      property.load(database);
      if (edgeSize != 0) {
	property.edgeSizeForSearch = edgeSize;
      }
      NGT::Index	index(database, property);
      ifstream		is(query);
      if (!is) {
	cerr << "Cannot open the specified file. " << query << endl;
	return;
      }
      if (outputMode == 's') { cout << "# Beginning of Evaluation" << endl; }
      string line;
      double totalTime	= 0;
      int    queryCount	= 0;
      while(getline(is, line)) {
	NGT::Object *object = index.allocateObject(line, " \t");
	queryCount++;
	for (float epsilon = beginOfEpsilon; epsilon <= endOfEpsilon; epsilon += stepOfEpsilon) {
	  epsilon = roundFloat(epsilon, 6);
	  NGT::SearchContainer sc(*object);
	  NGT::ObjectDistances objects;
	  sc.setResults(&objects);
	  sc.setSize(size);
	  sc.setRadius(radius);
	  sc.setEpsilon(epsilon);
	  if (debugLevel >= 1) {
	    cerr << "size=" << sc.size << endl;
	    cerr << "explorationCoefficient=" << sc.explorationCoefficient << endl;
	  }
	  NGT::Timer timer;
	  try {
	    if (outputMode == 'e') {
	      switch (indexType) {
	      case 't': index.search(sc); break;
	      case 'g': index.searchUsingOnlyGraph(sc); break;
	      case 's': index.linearSearch(sc); break;
	      }
	    }
	    switch (indexType) {
	    case 't': timer.start(); index.search(sc); timer.stop(); break;
	    case 'g': timer.start(); index.searchUsingOnlyGraph(sc); timer.stop(); break;
	    case 's': timer.start(); index.linearSearch(sc); timer.stop(); break;
	    }
	  } catch (NGT::Exception &err) {
	    throw err;
	  }
	  totalTime += timer.time;
	  if (outputMode == 'e') {
	    cout << "# Query No.=" << queryCount << endl;
	    cout << "# Query=" << line.substr(0, 20) + " ..." << endl;
	    cout << "# Index Type=" << indexType << endl;
	    cout << "# Size=" << size << endl;
	    cout << "# Radius=" << radius << endl;
	    cout << "# Epsilon=" << epsilon << endl;
	    cout << "# Query Time (msec)=" << timer.time * 1000.0 << endl;
	  } else {
	    cout << "Query No." << queryCount << endl;
	    cout << "Rank\tID\tDistance" << endl;
	  }
	  for (size_t i = 0; i < objects.size(); i++) {
	    cout << i + 1 << "\t" << objects[i].id << "\t";
	    cout << objects[i].distance << endl;
	  }
	  if (outputMode == 'e') {
	    cout << "# End of Search" << endl;
	  } else {
	    cout << "Query Time= " << timer.time << " (sec), " << timer.time * 1000.0 << " (msec)" << endl;
	  }
	} // for
	if (outputMode == 'e') {
	  cout << "# End of Query" << endl;
	}
	index.deleteObject(object);
      } // while
      if (outputMode == 'e') {
	cout << "# Average Query Time (msec)=" << totalTime * 1000.0 / (double)queryCount << endl;
	cout << "# Number of queries=" << queryCount << endl;
	cout << "# End of Evaluation" << endl;

	// show graph information
	size_t esize = edgeSize;
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

	cerr << "# of nodes=" << numberOfNodes << endl;
	cerr << "# of edges=" << numberOfEdges << endl;
	cerr << "Average number of edges=" << (double)numberOfEdges / (double)numberOfNodes << endl;
	cerr << "Average distance of edges=" << setprecision(10) << distance / (double)numberOfEdges << endl;

      } else {
	cout << "Average Query Time= " << totalTime / (double)queryCount  << " (sec), " 
	     << totalTime * 1000.0 / (double)queryCount << " (msec), (" 
	     << totalTime << "/" << queryCount << ")" << endl;
      }


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
	  //for (size_t i = 0; i < node.size(); ++i) {
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
    const string usage = "Usage: ngt info [-E #-of edges] [-g virtually-created-graph-type] [-m (h)] index(output)";
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
      long double distance = 0.0;
      size_t numberOfNodes = 0;
      size_t numberOfOutdegree = 0;
      size_t numberOfNodesWithoutEdges = 0;
      size_t maxNumberOfOutdegree = 0;
      size_t minNumberOfOutdegree = SIZE_MAX;
      vector<uint64_t> indegreeCount;
      vector<size_t> outdegreeHistogram;
      vector<size_t> indegreeHistogram;
      indegreeCount.resize(index.repository.size(), 0);
      for (size_t id = 1; id < index.repository.size(); id++) {
	NGT::GraphNode *node = 0;
	try {
	  node = index.getNode(id);
	} catch(NGT::Exception &err) {
	  cerr << "ngt info: Warning. Cannot get the node. ID=" << id << ":" << err.what() << endl;
	  continue;
	}
	numberOfNodes++;
	if (numberOfNodes % 1000000 == 0) {
	  cerr << "Processed " << numberOfNodes << endl;
	}
	size_t esize = node->size() > edgeSize ? edgeSize : node->size();
	if (esize == 0) {
	  numberOfNodesWithoutEdges++;
	}
	if (esize > maxNumberOfOutdegree) {
	  maxNumberOfOutdegree = esize;
	}
	if (esize < minNumberOfOutdegree) {
	  minNumberOfOutdegree = esize;
	}
	if (outdegreeHistogram.size() <= esize) {
	  outdegreeHistogram.resize(esize + 1);
	}
	outdegreeHistogram[esize]++;
	for (size_t i = 0; i < esize; i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  NGT::ObjectDistance &n = (*node).at(i, index.repository.allocator);
#else
	  NGT::ObjectDistance &n = (*node)[i];
#endif
	  if (n.id == 0) {
	    cerr << "ngt info: Warning. id is zero." << endl;
	  }
	  indegreeCount[n.id]++;
	  numberOfOutdegree++;
	  distance += n.distance;
	}
      }
      size_t numberOfNodesWithoutIndegree = 0;
      size_t maxNumberOfIndegree = 0;
      size_t minNumberOfIndegree = SIZE_MAX;
      for (size_t id = 1; id < index.repository.size(); id++) {
	if (indegreeCount[id] == 0) {
	  numberOfNodesWithoutIndegree++;
	}
	if (indegreeCount[id] > maxNumberOfIndegree) {
	  maxNumberOfIndegree = indegreeCount[id];
	}
	if (indegreeCount[id] < minNumberOfIndegree) {
	  minNumberOfIndegree = indegreeCount[id];
	}
	if (indegreeHistogram.size() <= indegreeCount[id]) {
	  indegreeHistogram.resize(indegreeCount[id] + 1);
	}
	indegreeHistogram[indegreeCount[id]]++;
      }
      cout << "# of nodes=" << numberOfNodes << endl;
      cout << "# of edges=" << numberOfOutdegree << endl;
      cout << "# of nodes without edges=" << numberOfNodesWithoutEdges << endl;
      cout << "Max outdegree=" << maxNumberOfOutdegree << endl;
      cout << "Min outdegree=" << minNumberOfOutdegree << endl;
      cout << "Average number of edges=" << (double)numberOfOutdegree / (double)numberOfNodes << endl;
      cout << "Average distance of edges=" << setprecision(10) << distance / (double)numberOfOutdegree << endl;
      cout << "# of nodes where indegree is 0=" << numberOfNodesWithoutIndegree << endl;
      cout << "Max indegree=" << maxNumberOfIndegree << endl;
      cout << "Min indegree=" << minNumberOfIndegree << endl;

      if (mode == 'h') {
	cout << "#\tout\tin" << endl;
	for (size_t i = 0; i < outdegreeHistogram.size() || i < indegreeHistogram.size(); i++) {
	  size_t out = outdegreeHistogram.size() <= i ? 0 : outdegreeHistogram[i];
	  size_t in = indegreeHistogram.size() <= i ? 0 : indegreeHistogram[i];
	  cout << i << "\t" << out << "\t" << in << endl;
	}
      }

    } catch (NGT::Exception &err) {
      cerr << "ngt: Error " << err.what() << endl;
      cerr << usage << endl;
    } catch (...) {
      cerr << "ngt: Error" << endl;
      cerr << usage << endl;
    }
  }


  //////////////////////////////////////////////////////////////

  void help() {
    cerr << "Usage : ngt command database [data]" << endl;
    cerr << "           commands : create search remove append export import prune" << endl;
  }

  void execute(Args args) {
    string command;
    try {
      command = args.get("#0");
    } catch(...) {
      help();
      return;
    }

    debugLevel = args.getl("X", 0);

    try {
      if (debugLevel >= 1) {
	cerr << "ngt: command=" << command << endl;
      }
      if (command == "search") {
	search(args);
      } else if (command == "create") {
	create(args);
      } else if (command == "append") {
	append(args);
      } else if (command == "remove") {
	remove(args);
      } else if (command == "export") {
	exportIndex(args);
      } else if (command == "import") {
	importIndex(args);
      } else if (command == "prune") {
	prune(args);
      } else if (command == "info") {
	info(args);
      } else {
	cerr << "ngt: Error: Illegal command. " << command << endl;
      }
    } catch(NGT::Exception &err) {
      cerr << "ngt: Error: " << err.what() << endl;
    }
  }

  int debugLevel;

};

int
main(int argc, char **argv)
{
  Args args(argc, argv);

  NGTCommand ngt;

  ngt.execute(args);

}


