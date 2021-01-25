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


#include	"NGT/NGTQ/NGTQGCommand.h"

#if !defined(NGT_SHARED_MEMORY_ALLOCATOR) && !defined(NGTQ_SHARED_INVERTED_INDEX)

using namespace std;


void 
NGTQG::Command::create(NGT::Args &args)
{

  const string usage = "Usage: ngtqg create "
    "[-D distance-function] "
    "[-p #-of-thread] [-d dimension] [-R global-codebook-range] [-r local-codebook-range] "
    "[-C global-codebook-size-limit] [-c local-codebook-size-limit] "
    "[-Q quantization-ratio] [-i index-type (t:Tree|g:Graph)] "
    "[-M global-centroid-creation-mode (d|s)] [-l global-centroid-creation-mode (d|k|s)] "
    "[-s local-sample-coefficient] "
    "index(output)";

  try {
    NGT::Command::CreateParameters createParameters(args);

    switch (createParameters.indexType) {
    case 't':
      NGT::Index::createGraphAndTree(createParameters.index, createParameters.property, createParameters.objectPath, createParameters.numOfObjects);
      break;
    case 'g':
      NGT::Index::createGraph(createParameters.index, createParameters.property, createParameters.objectPath, createParameters.numOfObjects);
      break;
    }
  } catch(NGT::Exception &err) {
    std::cerr << err.what() << std::endl;
    cerr << usage << endl;
  }

  NGTQG::Command::CreateParameters createParameters(args);

  try {
    char localCentroidCreationMode = args.getChar("l", 'd');
    switch(localCentroidCreationMode) {
    case 'd': createParameters.property.localCentroidCreationMode = NGTQ::CentroidCreationModeDynamic; break;
    case 's': createParameters.property.localCentroidCreationMode = NGTQ::CentroidCreationModeStatic; break;
    case 'k': createParameters.property.localCentroidCreationMode = NGTQ::CentroidCreationModeDynamicKmeans; break;
    default:
      cerr << "ngt: Invalid centroid creation mode. " << localCentroidCreationMode << endl;
      cerr << usage << endl;
      return;
    }
  
    createParameters.index += "/qg";
  
    cerr << "ngtqg: Create" << endl;
    NGTQ::Index::create(createParameters.index, createParameters.property, createParameters.globalProperty, createParameters.localProperty);
  } catch(NGT::Exception &err) {
    std::cerr << err.what() << std::endl;
    cerr << usage << endl;
  }


}

void 
NGTQG::Command::build(NGT::Args &args)
{

  NGT::Command::append(args);

}


void 
NGTQG::Command::quantize(NGT::Args &args)
{
  const std::string usage = "Usage: ngtqg quantize "
    "[-m quantization-mode(q|g|a)] [-E max-number-of-edges] [creation parameters] index\n"
    "\t-m mode\n"
    "\t\ta: quantize the objects and build and save a quantized graph. (default)\n"
    "\t\tq: just quantize the objects but not build and save a quantized graph.\n"
    "\t\tg: not quantize the objects but build and save a quantized graph.\n";
  string indexPath;
  try {
    indexPath = args.get("#1");
  } catch (...) {
    cerr << "An index is not specified." << endl;
    cerr << usage << endl;
    return;
  }
  char mode = args.getChar("m", 'a');
  size_t maxNumOfEdges = args.getl("E", 128);
  if (mode == 'q') {
    maxNumOfEdges = 0;
  }
  if (mode != 'g') {
    NGT::Index	index(indexPath);
    NGT::ObjectSpace &objectSpace = index.getObjectSpace();


    {
      std::string quantizedIndexPath = indexPath + "/qg";
      struct stat st;
      if (stat(quantizedIndexPath.c_str(), &st) != 0) {
	NGT::Property property;
	index.getProperty(property);
	NGTQG::Command::CreateParameters createParameters(args, property.dimension);

	try {
	  createParameters.property.centroidCreationMode = NGTQ::CentroidCreationModeStatic;
	  NGTQ::Index::create(quantizedIndexPath, createParameters.property, createParameters.globalProperty, createParameters.localProperty);
	} catch(NGT::Exception &err) {
	  std::cerr << err.what() << std::endl;
	  //cerr << usage << endl;
	}
      }

      NGTQ::Index	quantizedIndex(quantizedIndexPath);
      NGTQ::Quantizer	&quantizer = quantizedIndex.getQuantizer();

      {
	std::vector<float> meanObject(objectSpace.getDimension(), 0);
	quantizedIndex.getQuantizer().globalCodebook.insert(meanObject);
	quantizedIndex.getQuantizer().globalCodebook.createIndex(8);
      }
    
      vector<pair<NGT::Object*, size_t> > objects;
      for (size_t id = 1; id < objectSpace.getRepository().size(); id++) {
	if (id % 100000 == 0) {
	  std::cerr << "Processed " << id << " objects." << std::endl;
	}
	std::vector<float> object;
	try {
	  objectSpace.getObject(id, object);
	} catch(...) {
	  continue;
	}
	quantizer.insert(object, objects, id);
      }
      if (objects.size() > 0) {
	quantizer.insert(objects);
      }

      quantizedIndex.save();
      quantizedIndex.close();
    }
  }
  if (maxNumOfEdges != 0) {
    NGTQG::Index index(indexPath, maxNumOfEdges);
    index.save();
  }
}

void
search(NGTQG::Index &index, NGTQG::Command::SearchParameters &searchParameters, ostream &stream)
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

  while(getline(is, line)) {
    if (searchParameters.querySize > 0 && queryCount >= searchParameters.querySize) {
      break;
    }
    vector<float>	query;
    stringstream	linestream(line);
    while (!linestream.eof()) {
      float value;
      linestream >> value;
      if (linestream.fail()) {
	NGTThrowException("NGTQG: invalid stream.");
      }
      query.push_back(value);
    }
    queryCount++;
    float beginOfParam = 0;
    float endOfParam = 0;
    float stepOfParam = FLT_MAX;
    bool resultExpansionEnabled = true;
    if (searchParameters.beginOfResultExpansion != searchParameters.endOfResultExpansion) {
      beginOfParam = searchParameters.beginOfResultExpansion;
      endOfParam   = searchParameters.endOfResultExpansion;
      stepOfParam  = searchParameters.stepOfResultExpansion;
      resultExpansionEnabled = true;
    } else if (searchParameters.beginOfEpsilon != searchParameters.endOfEpsilon) {
      beginOfParam = searchParameters.beginOfEpsilon;
      endOfParam   = searchParameters.endOfEpsilon;
      stepOfParam  = searchParameters.stepOfEpsilon;
      resultExpansionEnabled = false;
    }
    for (float param = beginOfParam; param <= endOfParam; param += stepOfParam) {
      NGTQG::SearchQuery	searchQuery(query);
      NGT::ObjectDistances	objects;
      searchQuery.setResults(&objects);
      searchQuery.setSize(searchParameters.size);
      searchQuery.setRadius(searchParameters.radius);
      float epsilon, resultExpansion;
      if (stepOfParam == FLT_MAX) {
	resultExpansion = searchParameters.beginOfResultExpansion;
	epsilon = searchParameters.beginOfEpsilon;
      } else if (resultExpansionEnabled) {
	resultExpansion = param;
	epsilon = searchParameters.beginOfEpsilon;
      } else {
	resultExpansion = searchParameters.beginOfResultExpansion;
	epsilon = param;
      }
      searchQuery.setResultExpansion(resultExpansion);
      searchQuery.setEpsilon(epsilon);
      searchQuery.setEdgeSize(searchParameters.edgeSize);
      NGT::Timer timer;
      switch (searchParameters.indexType) {
      case 't': timer.start(); index.NGTQG::Index::search(searchQuery); timer.stop(); break;
      case 's': timer.start(); index.linearSearch(searchQuery); timer.stop(); break;
      }
      totalTime += timer.time;
      if (searchParameters.outputMode[0] == 'e') {
	stream << "# Query No.=" << queryCount << endl;
	stream << "# Query=" << line.substr(0, 20) + " ..." << endl;
	stream << "# Index Type=" << searchParameters.indexType << endl;
	stream << "# Size=" << searchParameters.size << endl;
	stream << "# Radius=" << searchParameters.radius << endl;
	stream << "# Epsilon*=" << epsilon << endl;
	stream << "# Result Expansion*=" << resultExpansion << endl;
	stream << "# Factor=" << param << endl;
	stream << "# Query Time (msec)=" << timer.time * 1000.0 << endl;
	stream << "# Distance Computation=" << searchQuery.distanceComputationCount << endl;
	stream << "# Visit Count=" << searchQuery.visitCount << endl;
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
NGTQG::Command::search(NGT::Args &args) {
  const string usage = "Usage: ngtqg search [-i index-type(g|t|s)] [-n result-size] [-e epsilon] [-E edge-size] "
    "[-o output-mode] [-p result-expansion] index(input) query.tsv(input)";

  string indexPath;
  try {
    indexPath = args.get("#1");
  } catch (...) {
    cerr << "An index is not specified." << endl;
    cerr << usage << endl;
    return;
  }
  NGTQG::Command::SearchParameters searchParameters(args);

  NGTQG::Index index(indexPath, 128);

  if (debugLevel >= 1) {
    cerr << "indexType=" << searchParameters.indexType << endl;
    cerr << "size=" << searchParameters.size << endl;
    cerr << "edgeSize=" << searchParameters.edgeSize << endl;
    cerr << "epsilon=" << searchParameters.beginOfEpsilon << "<->" << searchParameters.endOfEpsilon << "," 
	 << searchParameters.stepOfEpsilon << endl;
  }

  try {
    ::search(index, searchParameters, std::cout);
  } catch (NGT::Exception &err) {
    cerr << "ngtqg: Error " << err.what() << endl;
    cerr << usage << endl;
  } catch (std::exception &err) {
    cerr << "ngtqg: Error " << err.what() << endl;
    cerr << usage << endl;
  } catch (...) {
    cerr << "ngtqg: Error" << endl;
    cerr << usage << endl;
  }

}


void
NGTQG::Command::info(NGT::Args &args)
{
  const string usage = "Usage: ngt info [-E #-of-edges] [-m h|e] index";

  cerr << "NGT version: " << NGT::Index::getVersion() << endl;

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

#endif 
