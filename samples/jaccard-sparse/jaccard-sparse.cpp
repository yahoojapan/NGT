
// sort -R sparse_binary.tsv |head -10 > sparse_binary_query_10.tsv
// ./jaccard-sparse create -d 100 -D J sparse
// ./jaccard-sparse append sparse sparse_binary.tsv
// ./jaccard-sparse search sparse sparse_binary_query_10.tsv
//

#include	"NGT/Command.h"

using namespace std;

void help() {
  cerr << "Usage : jaccard-sparse command index [data]" << endl;
  cerr << "           command : info create search append" << endl;
}

void 
append(NGT::Args &args)
{
  const string usage = "Usage: jaccard-sparse append [-p #-of-thread] [-n data-size] "
    "index(output) [data.tsv(input)]";
  string database;
  try {
    database = args.get("#1");
  } catch (...) {
    cerr << "jaccard-sparse: Error: DB is not specified." << endl;
    cerr << usage << endl;
    return;
  }
  string data;
  try {
    data = args.get("#2");
  } catch (...) {
    cerr << "jaccard-sparse: Warning: No specified object file. Just build an index for the existing objects." << endl;
  }

  int threadSize = args.getl("p", 50);
  size_t dataSize = args.getl("n", 0);

  std::istream *is;
  std::ifstream *ifs = 0;

  try {
    NGT::Index	index(database);
    if (data == "-") {
      is = &std::cin;
    } else {
      ifs = new std::ifstream;
      ifs->std::ifstream::open(data);
      if (!(*ifs)) {
	cerr << "Cannot open the specified data file. " << data << endl;
	return;
      }
      is = ifs;
    }
    string line;
    size_t count = 0;
    while(getline(*is, line)) {
      if (dataSize > 0 && count >= dataSize) {
	break;
      }
      count++;
      vector<uint32_t> object;
      stringstream	linestream(line);
      while (!linestream.eof()) {
	uint32_t value;
	linestream >> value;
	if (linestream.fail()) {
	  object.clear();
	  break;
	}
	object.push_back(value);
      }
      if (object.empty()) {
	std::cerr << "jaccard-sparse: Empty line or invalid value. " << count << ":" << line << std::endl;
	continue;
      }
    }
    if (data != "-") {
      delete ifs;
    }
    index.createIndex(threadSize);
    index.saveIndex(database);
  } catch (NGT::Exception &err) {
    if (data != "-") {
      delete ifs;
    }
    cerr << "jaccard-sparse: Error " << err.what() << endl;
    cerr << usage << endl;
  }
  return;
}


void
search(NGT::Index &index, NGT::Command::SearchParameters &searchParameters, ostream &stream)
{

  std::ifstream		is(searchParameters.query);
  if (!is) {
    std::cerr << "Cannot open the specified file. " << searchParameters.query << std::endl;
    return;
  }

  if (searchParameters.outputMode[0] == 'e') { 
    stream << "# Beginning of Evaluation" << endl; 
  }

  string line;
  double totalTime	= 0;
  size_t queryCount	= 0;
  double epsilon	= searchParameters.beginOfEpsilon;

  while(getline(is, line)) {
    if (searchParameters.querySize > 0 && queryCount >= searchParameters.querySize) {
      break;
    }
    vector<uint32_t>	query;
    stringstream	linestream(line);
    while (!linestream.eof()) {
      uint32_t value;
      linestream >> value;
      query.push_back(value);
    }
    auto sparseQuery = index.makeSparseObject(query);
    queryCount++;
    NGT::SearchQuery	sc(sparseQuery);
    NGT::ObjectDistances objects;
    sc.setResults(&objects);
    sc.setSize(searchParameters.size);
    sc.setRadius(searchParameters.radius);
    if (searchParameters.accuracy > 0.0) {
      sc.setExpectedAccuracy(searchParameters.accuracy);
    } else {
      sc.setEpsilon(epsilon);
    }
    sc.setEdgeSize(searchParameters.edgeSize);
    NGT::Timer timer;
    switch (searchParameters.indexType) {
    case 't': timer.start(); index.search(sc); timer.stop(); break;
    case 'g': timer.start(); index.searchUsingOnlyGraph(sc); timer.stop(); break;
    case 's': timer.start(); index.linearSearch(sc); timer.stop(); break;
    }
    totalTime += timer.time;
    if (searchParameters.outputMode[0] == 'e') {
      stream << "# Query No.=" << queryCount << endl;
      stream << "# Query=" << line.substr(0, 20) + " ..." << endl;
      stream << "# Index Type=" << searchParameters.indexType << endl;
      stream << "# Size=" << searchParameters.size << endl;
      stream << "# Radius=" << searchParameters.radius << endl;
      stream << "# Epsilon=" << epsilon << endl;
      stream << "# Query Time (msec)=" << timer.time * 1000.0 << endl;
      stream << "# Distance Computation=" << sc.distanceComputationCount << endl;
      stream << "# Visit Count=" << sc.visitCount << endl;
    } else {
      stream << "Query No." << queryCount << endl;
      stream << "Rank\tID\tDistance" << endl;
    }
    for (size_t i = 0; i < objects.size(); i++) {
      stream << i + 1 << "\t" << objects[i].id << "\t";
      stream << objects[i].distance << endl;
    }
    if (searchParameters.outputMode[0] == 'e') {
      stream << "# End of Search" << endl;
    } else {
      stream << "Query Time= " << timer.time << " (sec), " << timer.time * 1000.0 << " (msec)" << endl;
    }
    if (searchParameters.outputMode[0] == 'e') {
      stream << "# End of Query" << endl;
    }
  }
  if (searchParameters.outputMode[0] == 'e') {
    stream << "# Average Query Time (msec)=" << totalTime * 1000.0 / (double)queryCount << endl;
    stream << "# Number of queries=" << queryCount << endl;
    stream << "# End of Evaluation" << endl;
  } else {
    stream << "Average Query Time= " << totalTime / (double)queryCount  << " (sec), " 
	   << totalTime * 1000.0 / (double)queryCount << " (msec), (" 
	   << totalTime << "/" << queryCount << ")" << endl;
  }
}

void
search(NGT::Args &args) {
  const string usage = "Usage: ngt search [-i index-type(g|t|s)] [-n result-size] [-e epsilon] [-E edge-size] "
    "[-m open-mode(r|w)] [-o output-mode] index(input) query.tsv(input)";

  string database;
  try {
    database = args.get("#1");
  } catch (...) {
    cerr << "jaccard-sparse: Error: DB is not specified" << endl;
    cerr << usage << endl;
    return;
  }

  NGT::Command::SearchParameters searchParameters(args);

  try {
    NGT::Index	index(database, searchParameters.openMode == 'r');
    search(index, searchParameters, cout);
  } catch (NGT::Exception &err) {
    cerr << "jaccard-sparse: Error " << err.what() << endl;
    cerr << usage << endl;
  } catch (...) {
    cerr << "jaccard-sparse: Error" << endl;
    cerr << usage << endl;
  }

}

int
main(int argc, char **argv)
{

  NGT::Args args(argc, argv);

  NGT::Command ngt;

  string command;
  try {
    command = args.get("#0");
  } catch(...) {
    help();
    return 0;
  }

  try {
    if (command == "create") {
      ngt.create(args);
    } else if (command == "append") {
      append(args);
    } else if (command == "search") {
      search(args);
    } else {
      cerr << "jaccard-sparse: Error: Illegal command. " << command << endl;
      help();
    }
  } catch(NGT::Exception &err) {
    cerr << "jaccard-sparse: Error: " << err.what() << endl;
    help();
    return 0;
  }
  return 0;
}


