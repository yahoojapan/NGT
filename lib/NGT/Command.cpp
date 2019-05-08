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

#include	"NGT/Command.h"
#include	"NGT/GraphReconstructor.h"

#include	"NGT/Optimizer.h"


  void 
  NGT::Command::create(Args &args)
  {
    const string usage = "Usage: ngt create "
      "-d dimension [-p #-of-thread] [-i index-type(t|g)] [-g graph-type(a|k|b|o|i)] "
      "[-t truncation-edge-limit] [-E edge-size] [-S edge-size-for-search] [-L edge-size-limit] "
      "[-e epsilon] [-o object-type(f|c)] [-D distance-function(1|2|a|A|h|c|C)] [-n #-of-inserted-objects] "
      "[-P path-adjustment-interval] [-B dynamic-edge-size-base] [-A object-alignment(t|f)] "
      "[-T build-time-limit] [-O outgoing x incoming] "
      "index(output) [data.tsv(input)]";
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
    } catch (...) {}

    NGT::Property property;

    property.edgeSizeForCreation = args.getl("E", 10);
    property.edgeSizeForSearch = args.getl("S", 40);
    property.batchSizeForCreation = args.getl("b", 200);
    property.insertionRadiusCoefficient = args.getf("e", 0.1) + 1.0;
    property.truncationThreshold = args.getl("t", 0);
    property.dimension = args.getl("d", 0);
    property.threadPoolSize = args.getl("p", 24);
    property.pathAdjustmentInterval = args.getl("P", 0);
    property.dynamicEdgeSizeBase = args.getl("B", 30);
    property.buildTimeLimit = args.getf("T", 0.0);

    if (property.dimension <= 0) {
      cerr << "ngt: Error: Specify greater than 0 for # of your data dimension by a parameter -d." << endl;
      cerr << usage << endl;
      return;
    }

    property.objectAlignment = args.getChar("A", 'f') == 't' ? NGT::Property::ObjectAlignmentTrue : NGT::Property::ObjectAlignmentFalse;

    char graphType = args.getChar("g", 'a');
    switch(graphType) {
    case 'a': property.graphType = NGT::Property::GraphType::GraphTypeANNG; break;
    case 'k': property.graphType = NGT::Property::GraphType::GraphTypeKNNG; break;
    case 'b': property.graphType = NGT::Property::GraphType::GraphTypeBKNNG; break;
    case 'd': property.graphType = NGT::Property::GraphType::GraphTypeDNNG; break;
    case 'o': property.graphType = NGT::Property::GraphType::GraphTypeONNG; break;
    case 'i': property.graphType = NGT::Property::GraphType::GraphTypeIANNG; break;
    default:
      cerr << "ngt: Error: Invalid graph type. " << graphType << endl;
      cerr << usage << endl;
      return;
    }    

    if (property.graphType == NGT::Property::GraphType::GraphTypeONNG) {
      property.outgoingEdge = 10;
      property.incomingEdge = 80;
      string str = args.getString("O", "-");
      if (str != "-") {
	vector<string> tokens;
	NGT::Common::tokenize(str, tokens, "x");
	if (str != "-" && tokens.size() != 2) {
	  cerr << "ngt: Error: outgoing/incoming edge size specification is invalid. (out)x(in) " << str << endl;
	  cerr << usage << endl;
	  return;
	}
	property.outgoingEdge = NGT::Common::strtod(tokens[0]);
	property.incomingEdge = NGT::Common::strtod(tokens[1]);
	cerr << "ngt: ONNG out x in=" << property.outgoingEdge << "x" << property.incomingEdge << endl;
      }
    }

    char seedType = args.getChar("s", '-');
    switch(seedType) {
    case 'f': property.seedType = NGT::Property::SeedType::SeedTypeFixedNodes; break;
    case '1': property.seedType = NGT::Property::SeedType::SeedTypeFirstNode; break;
    case 'r': property.seedType = NGT::Property::SeedType::SeedTypeRandomNodes; break;
    case 'l': property.seedType = NGT::Property::SeedType::SeedTypeAllLeafNodes; break;
    default:
    case '-': property.seedType = NGT::Property::SeedType::SeedTypeNone; break;
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
    case 'A':
      property.distanceType = NGT::Index::Property::DistanceType::DistanceTypeNormalizedAngle;
      break;
    case 'h':
      property.distanceType = NGT::Index::Property::DistanceType::DistanceTypeHamming;
      break;
    case 'c':
      property.distanceType = NGT::Index::Property::DistanceType::DistanceTypeCosine;
      break;
    case 'C':
      property.distanceType = NGT::Index::Property::DistanceType::DistanceTypeNormalizedCosine;
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
  NGT::Command::append(Args &args)
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
  NGT::Command::search(NGT::Index &index, NGT::Command::SearchParameter &searchParameter, istream &is, ostream &stream)
  {

    if (searchParameter.outputMode[0] == 'e') { 
      stream << "# Beginning of Evaluation" << endl; 
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
	  if (searchParameter.outputMode != "ei") {
	    // not ignore exceptions
	    throw err;
	  }
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
  NGT::Command::search(Args &args) {
    const string usage = "Usage: ngt search [-i index-type(g|t|s)] [-n result-size] [-e epsilon] [-E edge-size] "
      "[-m open-mode(r|w)] [-o output-mode] index(input) query.tsv(input)";

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
      NGT::Index	index(database, searchParameter.openMode == 'r');
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
  NGT::Command::remove(Args &args)
  {
    const string usage = "Usage: ngt remove [-d object-ID-type(f|d)] [-m f] index(input) object-ID(input)";
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
    char mode = args.getChar("m", '-');
    bool force = false;
    if (mode == 'f') {
      force = true;
    }
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
      NGT::Index::remove(database, objects, force);
    } catch (NGT::Exception &err) {
      cerr << "ngt: Error " << err.what() << endl;
      cerr << usage << endl;
    } catch (...) {
      cerr << "ngt: Error" << endl;
      cerr << usage << endl;
    }
  }

  void
  NGT::Command::exportIndex(Args &args)
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
  NGT::Command::importIndex(Args &args)
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
  NGT::Command::prune(Args &args)
  {
    const string usage = "Usage: ngt prune -e #-of-forcedly-pruned-edges -s #-of-selecively-pruned-edge index(in/out)";
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
  NGT::Command::reconstructGraph(Args &args)
  {
    const string usage = "Usage: ngt reconstruct-graph [-i index-type] [-m mode] [-P path-adjustment-mode] -o #-of-original-edges -i #-of-reverse-edges index(input) index(output)\n"
      "\t-i input indes type\n"
      "\t\ta: anng\n"
      "\t\tn: not anng\n"
      "\t-m mode\n"
      "\t\ts: Edge adjustment. (default)\n"
      "\t\tS: Edge adjustment and path adjustment.\n"
      "\t\tc: Edge adjustment with the constraint.\n"
      "\t\tC: Edge adjustment with the constraint and path adjustment.\n"
      "\t\tP: Path adjustment.\n"
      "\t-P path-adjustment-mode\n"
      "\t\ta: Advanced method. High-speed. Not guarantee the paper method. (default)\n"
      "\t\tothers: Slow and less memory usage, but guarantee the paper method.\n";

    string inIndexName;
    try {
      inIndexName = args.get("#1");
    } catch (...) {
      cerr << "ngt::reconstructGraph: Input index is not specified." << endl;
      cerr << usage << endl;
      return;
    }
    string outIndexName;
    try {
      outIndexName = args.get("#2");
    } catch (...) {
      cerr << "ngt::reconstructGraph: Output index is not specified." << endl;
      cerr << usage << endl;
      return;
    }

    // the number (rank) of original edges
    int originalEdgeSize	= args.getl("o", -1);
    // the number (rank) of reverse edges
    int reverseEdgeSize		= args.getl("i", -1);
    if ((originalEdgeSize < 0 && reverseEdgeSize >= 0) ||
	(originalEdgeSize >= 0 && reverseEdgeSize < 0)) {
      cerr << "ngt::reconstructGraph: specified both of the edges(-i -o) or neither of them." << endl;
      cerr << usage << endl;
      return;
    }
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    if (access(outIndexName.c_str(), 0) == 0) {
      cerr << "ngt:reconstructGraph: the specified index exists. " << outIndexName << endl;
      cerr << usage << endl;
      return;
    }
    const string com = "cp -r " + inIndexName + " " + outIndexName;
    system(com.c_str());
    NGT::Index	outIndex(outIndexName);
#else
    NGT::Index	outIndex(inIndexName);
#endif
    cerr << "ngt::reconstructGraph: Extract the graph data." << endl;
    // extract only edges from the index to reduce the memory usage.
    NGT::GraphIndex	&outGraph = (NGT::GraphIndex&)outIndex.getIndex();
    Timer timer;
    timer.start();
    vector<NGT::ObjectDistances> graph;
    graph.reserve(outGraph.repository.size());
    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      if (id % 1000000 == 0) {
	cerr << "Processed " << id << " objects." << endl;
      }
      try {
	NGT::GraphNode &node = *outGraph.getNode(id);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	NGT::ObjectDistances nd;
	nd.reserve(node.size());
	for (auto n = node.begin(outGraph.repository.allocator); n != node.end(outGraph.repository.allocator); ++n) {
	  nd.push_back(ObjectDistance((*n).id, (*n).distance));
        }
	graph.push_back(nd);
#else
	graph.push_back(node);
#endif
	if (graph.back().size() != graph.back().capacity()) {
	  cerr << "ngt::reconstructGraph: Warning! The graph size must be the same as the capacity. " << id << endl;
	}
      } catch(NGT::Exception &err) {
	cerr << "ngt::reconstructGraph: Warning! Cannot get the node. ID=" << id << ":" << err.what() << endl;
	continue;
      }
    }

    char mode = args.getChar("m", 's');
    char pamode = args.getChar("P", 'a');
    char indexType = args.getChar("I", 'a');

    if (originalEdgeSize >= 0) {
      switch (mode) {
      case 's': // SA
      case 'S': // SA and path adjustment
	if (indexType != 'a') {
	  NGT::GraphReconstructor::convertToANNG(graph);
	}
	NGT::GraphReconstructor::reconstructGraph(graph, outIndex, originalEdgeSize, reverseEdgeSize);
	break;
      case 'c': // SAC
      case 'C': // SAC and path adjustment
	if (indexType != 'a') {
	  NGT::GraphReconstructor::convertToANNG(graph);
	}
	NGT::GraphReconstructor::reconstructGraphWithConstraint(graph, outIndex, originalEdgeSize, reverseEdgeSize);
	break;
      case 'P':
	break;
      default:
	cerr << "ngt::reconstructGraph: Invalid mode. " << mode << endl;
	return;
      }
    }
    timer.stop();
    cerr << "ngt::Graph reconstruction time=" << timer.time << " (sec) " << endl;

    if (mode == 'S' || mode == 'C'|| mode == 'P') {
      timer.reset();
      timer.start();
      if (pamode == 'a') {
        GraphReconstructor::adjustPathsEffectively(outIndex);
      } else {
        GraphReconstructor::adjustPaths(outIndex);
      }
      timer.stop();
      cerr << "ngt::Path adjustment time=" << timer.time << " (sec) " << endl;
    }

    float intervalFrom = 0.6;
    float intervalTo = 0.8;
    size_t querySize = 100;
    double gtEpsilon = 0.1;

    try {
      size_t baseEdgeSize = NGT::Optimizer::adjustBaseSearchEdgeSize(outIndex, pair<float, float>(intervalFrom, intervalTo), querySize, gtEpsilon);
      NeighborhoodGraph::Property &prop = outGraph.getGraphProperty();
      prop.dynamicEdgeSizeBase = baseEdgeSize;
      cerr << "Reconstruct Graph: adjust the base search edge size. " << baseEdgeSize << endl;
    } catch(NGT::Exception &err) {
      cerr << "Warning: Cannot adjust the base edge size. " << err.what() << endl;
    }

    outGraph.saveIndex(outIndexName);

  }



  void
  NGT::Command::info(Args &args)
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
      NGT::Index	index(database);
      stringstream msg;
      size_t smsize = index.getSharedMemorySize(msg);
      cerr << "SharedMemorySize=" << smsize << endl;
      NGT::GraphIndex::showStatisticsOfGraph(static_cast<NGT::GraphIndex&>(index.getIndex()), mode, edgeSize);
      if (mode == 'v') {
	vector<uint8_t> status;
	index.verify(status);
      }
    } catch (NGT::Exception &err) {
      cerr << "ngt: Error " << err.what() << endl;
      cerr << usage << endl;
    } catch (...) {
      cerr << "ngt: Error" << endl;
      cerr << usage << endl;
    }
  }


