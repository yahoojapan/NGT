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

#include	"NGT/defines.h"
#include	"NGT/Common.h"
#include	"NGT/ObjectSpaceRepository.h"
#include	"NGT/Index.h"
#include	"NGT/Thread.h"
#include	"NGT/GraphReconstructor.h"
#include	"NGT/Version.h"

using namespace std;
using namespace NGT;


void 
Index::version(ostream &os) 
{
  os << "libngt:" << endl;
  Version::get(os);
}

string
Index::getVersion()
{
  return Version::getVersion();
}

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
NGT::Index::Index(NGT::Property &prop, const string &database) {
  if (prop.dimension == 0) {
    NGTThrowException("Index::Index. Dimension is not specified.");
  }
  Index* idx = 0;
  mkdir(database);
  if (prop.indexType == NGT::Index::Property::GraphAndTree) {
    idx = new NGT::GraphAndTreeIndex(database, prop);
  } else if (prop.indexType == NGT::Index::Property::Graph) {
    idx = new NGT::GraphIndex(database, prop);
  } else {
    NGTThrowException("Index::Index: Not found IndexType in property file.");
  }
  if (idx == 0) {
    stringstream msg;
    msg << "Index::Index: Cannot construct. ";
    NGTThrowException(msg);
  }
  index = idx;
  path = "";
}
#else
NGT::Index::Index(NGT::Property &prop) {
  if (prop.dimension == 0) {
    NGTThrowException("Index::Index. Dimension is not specified.");
  }
  Index* idx = 0;
  if (prop.indexType == NGT::Index::Property::GraphAndTree) {
    idx = new NGT::GraphAndTreeIndex(prop);
  } else if (prop.indexType == NGT::Index::Property::Graph) {
    idx = new NGT::GraphIndex(prop);
  } else {
    NGTThrowException("Index::Index: Not found IndexType in property file.");
  }
  if (idx == 0) {
    stringstream msg;
    msg << "Index::Index: Cannot construct. ";
    NGTThrowException(msg);
  }
  index = idx;
  path = "";
}
#endif

void 
NGT::Index::open(const string &database, bool rdOnly) {
  NGT::Property prop;
  prop.load(database);
  Index* idx = 0;
  if (prop.indexType == NGT::Index::Property::GraphAndTree) {
    idx = new NGT::GraphAndTreeIndex(database, rdOnly);
  } else if (prop.indexType == NGT::Index::Property::Graph) {
    idx = new NGT::GraphIndex(database, rdOnly);
  } else {
    NGTThrowException("Index::Open: Not found IndexType in property file.");
  }
  if (idx == 0) {
    stringstream msg;
    msg << "Index::open: Cannot open. " << database;
    NGTThrowException(msg);
  }
  index = idx;
  path = database;
}

void 
NGT::Index::createGraphAndTree(const string &database, NGT::Property &prop, const string &dataFile,
			       size_t dataSize, bool redirect) {
  if (prop.dimension == 0) {
    NGTThrowException("Index::createGraphAndTree. Dimension is not specified.");
  }
  prop.indexType = NGT::Index::Property::IndexType::GraphAndTree;
  Index *idx = 0;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  mkdir(database);
  idx = new NGT::GraphAndTreeIndex(database, prop);
#else
  idx = new NGT::GraphAndTreeIndex(prop);
#endif
  assert(idx != 0);
  StdOstreamRedirector redirector(redirect);
  redirector.begin();
  try {
    loadAndCreateIndex(*idx, database, dataFile, prop.threadPoolSize, dataSize);
  } catch(Exception &err) {
    delete idx;
    redirector.end();
    throw err;
  }
  delete idx;
  redirector.end();
}

void 
NGT::Index::createGraph(const string &database, NGT::Property &prop, const string &dataFile, size_t dataSize, bool redirect) {
  if (prop.dimension == 0) {
    NGTThrowException("Index::createGraphAndTree. Dimension is not specified.");
  }
  prop.indexType = NGT::Index::Property::IndexType::Graph;
  Index *idx = 0;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  mkdir(database);
  idx = new NGT::GraphIndex(database, prop);
#else
  idx = new NGT::GraphIndex(prop);
#endif
  assert(idx != 0);
  StdOstreamRedirector redirector(redirect);
  redirector.begin();
  try {
    loadAndCreateIndex(*idx, database, dataFile, prop.threadPoolSize, dataSize);
  } catch(Exception &err) {
    delete idx;
    redirector.end();
    throw err;
  }
  delete idx;
  redirector.end();
}

void 
NGT::Index::loadAndCreateIndex(Index &index, const string &database, const string &dataFile, size_t threadSize, size_t dataSize) {
  NGT::Timer timer;
  timer.start();
  if (dataFile.size() != 0) {
    index.load(dataFile, dataSize);
  } else {
    index.saveIndex(database);
    return;
  }
  timer.stop();
  cerr << "Data loading time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  if (index.getObjectRepositorySize() == 0) {
    NGTThrowException("Index::create: Data file is empty.");
  }
  cerr << "# of objects=" << index.getObjectRepositorySize() - 1 << endl;
  timer.reset();
  timer.start();
  index.createIndex(threadSize);
  timer.stop();
  index.saveIndex(database);
  cerr << "Index creation time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
}

void 
NGT::Index::append(const string &database, const string &dataFile, size_t threadSize, size_t dataSize) {
  NGT::Index	index(database);
  NGT::Timer	timer;
  timer.start();
  if (dataFile.size() != 0) {
    index.append(dataFile, dataSize);
  } else {
    NGTThrowException("Index::append: No data file.");
  }
  timer.stop();
  cerr << "Data loading time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << index.getObjectRepositorySize() - 1 << endl;
  timer.reset();
  timer.start();
  index.createIndex(threadSize);
  timer.stop();
  index.saveIndex(database);
  cerr << "Index creation time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  return;
}

void 
NGT::Index::append(const string &database, const float *data, size_t dataSize, size_t threadSize) {
  NGT::Index	index(database);
  NGT::Timer	timer;
  timer.start();
  if (data != 0 && dataSize != 0) {
    index.append(data, dataSize);
  } else {
    NGTThrowException("Index::append: No data.");
  }
  timer.stop();
  cerr << "Data loading time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << index.getObjectRepositorySize() - 1 << endl;
  timer.reset();
  timer.start();
  index.createIndex(threadSize);
  timer.stop();
  index.saveIndex(database);
  cerr << "Index creation time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  return;
}

void 
NGT::Index::remove(const string &database, vector<ObjectID> &objects, bool force) {
  NGT::Index	index(database);
  NGT::Timer	timer;
  timer.start();
  for (vector<ObjectID>::iterator i = objects.begin(); i != objects.end(); i++) {
    try {
      index.remove(*i, force);
    } catch (Exception &err) {
      cerr << "Warning: Cannot remove the node. ID=" << *i << " : " << err.what() << endl;
      continue;
    }
  }
  timer.stop();
  cerr << "Data removing time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << index.getObjectRepositorySize() - 1 << endl;
  index.saveIndex(database);
  return;
}

void 
NGT::Index::importIndex(const string &database, const string &file) {
  Index *idx = 0;
  NGT::Property property;
  property.importProperty(file);
  NGT::Timer	timer;
  timer.start();
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  property.databaseType = NGT::Index::Property::DatabaseType::MemoryMappedFile;
  mkdir(database);
#else
  property.databaseType = NGT::Index::Property::DatabaseType::Memory;
#endif
  if (property.indexType == NGT::Index::Property::IndexType::GraphAndTree) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    idx = new NGT::GraphAndTreeIndex(database, property);
#else
    idx = new NGT::GraphAndTreeIndex(property);
#endif
    assert(idx != 0);
  } else if (property.indexType == NGT::Index::Property::IndexType::Graph) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    idx = new NGT::GraphIndex(database, property);
#else
    idx = new NGT::GraphIndex(property);
#endif
    assert(idx != 0);
  } else {
    NGTThrowException("Index::Open: Not found IndexType in property file.");
  }
  idx->importIndex(file);
  timer.stop();
  cerr << "Data importing time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << idx->getObjectRepositorySize() - 1 << endl;
  idx->saveIndex(database);
  delete idx;
}

void 
NGT::Index::exportIndex(const string &database, const string &file) {
  NGT::Index	idx(database);
  NGT::Timer	timer;
  timer.start();
  idx.exportIndex(file);
  timer.stop();
  cerr << "Data exporting time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << idx.getObjectRepositorySize() - 1 << endl;
}

void 
NGT::GraphIndex::constructObjectSpace(NGT::Property &prop) {
  assert(prop.dimension != 0);
  switch (prop.objectType) {
  case NGT::ObjectSpace::ObjectType::Float :
    objectSpace = new ObjectSpaceRepository<float, double>(prop.dimension, typeid(float), prop.distanceType);
    break;
  case NGT::ObjectSpace::ObjectType::Uint8 :
    objectSpace = new ObjectSpaceRepository<unsigned char, int>(prop.dimension, typeid(uint8_t), prop.distanceType);
    break;
  default:
    stringstream msg;
    msg << "Invalid Object Type in the property. " << prop.objectType;
    NGTThrowException(msg);	
  }
  prop.prefetchOffset = objectSpace->setPrefetchOffset(prop.prefetchOffset);
  prop.prefetchSize = objectSpace->setPrefetchSize(prop.prefetchSize);
}

void 
NGT::Index::Property::set(NGT::Property &prop) {
  if (prop.dimension != -1) dimension = prop.dimension;
  if (prop.threadPoolSize != -1) threadPoolSize = prop.threadPoolSize;
  if (prop.objectType != ObjectSpace::ObjectTypeNone) objectType = prop.objectType;
  if (prop.distanceType != DistanceType::DistanceTypeNone) distanceType = prop.distanceType;
  if (prop.indexType != IndexTypeNone) indexType = prop.indexType;
  if (prop.databaseType != DatabaseTypeNone) databaseType = prop.databaseType;
  if (prop.objectAlignment != ObjectAlignmentNone) objectAlignment = prop.objectAlignment;
  if (prop.pathAdjustmentInterval != -1) pathAdjustmentInterval = prop.pathAdjustmentInterval;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  if (prop.graphSharedMemorySize != -1) graphSharedMemorySize = prop.graphSharedMemorySize;
  if (prop.treeSharedMemorySize != -1) treeSharedMemorySize = prop.treeSharedMemorySize;
  if (prop.objectSharedMemorySize != -1) objectSharedMemorySize = prop.objectSharedMemorySize;
#endif
  if (prop.prefetchOffset != -1) prefetchOffset = prop.prefetchOffset;
  if (prop.prefetchSize != -1) prefetchSize = prop.prefetchSize;
}

void 
NGT::Index::Property::get(NGT::Property &prop) {
  prop.dimension = dimension;
  prop.threadPoolSize = threadPoolSize;
  prop.objectType = objectType;
  prop.distanceType = distanceType;
  prop.indexType = indexType;
  prop.databaseType = databaseType;
  prop.pathAdjustmentInterval = pathAdjustmentInterval;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  prop.graphSharedMemorySize = graphSharedMemorySize;
  prop.treeSharedMemorySize = treeSharedMemorySize;
  prop.objectSharedMemorySize = objectSharedMemorySize;
#endif
  prop.prefetchOffset = prefetchOffset;
  prop.prefetchSize = prefetchSize;
}

class CreateIndexJob {
public:
  CreateIndexJob() {}
  CreateIndexJob &operator=(const CreateIndexJob &d) {
    id = d.id;
    results = d.results;
    object = d.object;
    batchIdx = d.batchIdx;
    return *this;
  }
  friend bool operator<(const CreateIndexJob &ja, const CreateIndexJob &jb) { return ja.batchIdx < jb.batchIdx; }
  NGT::ObjectID		id;
  NGT::Object		*object;	// this will be a node of the graph later.
  NGT::ObjectDistances	*results;
  size_t		batchIdx;
};

class CreateIndexSharedData {
public:
  CreateIndexSharedData(NGT::GraphIndex &nngt) : graphIndex(nngt) {}
  NGT::GraphIndex &graphIndex;
};

class CreateIndexThread : public NGT::Thread {
public:
  CreateIndexThread() {}
  virtual ~CreateIndexThread() {}
  virtual int run();

};

 typedef NGT::ThreadPool<CreateIndexJob, CreateIndexSharedData*, CreateIndexThread> CreateIndexThreadPool;

int
CreateIndexThread::run() {

  NGT::ThreadPool<CreateIndexJob, CreateIndexSharedData*, CreateIndexThread>::Thread &poolThread =
    (NGT::ThreadPool<CreateIndexJob, CreateIndexSharedData*, CreateIndexThread>::Thread&)*this;

  CreateIndexSharedData &sd = *poolThread.getSharedData();
  NGT::GraphIndex &graphIndex = sd.graphIndex;

  for(;;) {
    CreateIndexJob job;
    try {
      poolThread.getInputJobQueue().popFront(job);
    } catch(NGT::ThreadTerminationException &err) {
      break;
    } catch(NGT::Exception &err) {
      cerr << "CreateIndex::search:Error! popFront " << err.what() << endl;
      break;
    }
    ObjectDistances *rs = new ObjectDistances;
    Object &obj = *job.object;
    try {
      if (graphIndex.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeKNNG) {
	graphIndex.searchForKNNGInsertion(obj, job.id, *rs);	// linear search
      } else {
	graphIndex.searchForNNGInsertion(obj, *rs);
      }
    } catch(NGT::Exception &err) {
      cerr << "CreateIndex::search:Fatal error! ID=" << job.id << " " << err.what() << endl;
      abort();
    } 
    job.results = rs;
    poolThread.getOutputJobQueue().pushBack(job);
  }

  return 0;

}

class BuildTimeController {
public:
  BuildTimeController(GraphIndex &graph, NeighborhoodGraph::Property &prop):property(prop) {
    noOfInsertedObjects = graph.objectSpace->getRepository().size() - graph.repository.size();
    interval = 10000;
    count = interval;
    edgeSizeSave = property.edgeSizeForCreation;
    insertionRadiusCoefficientSave = property.insertionRadiusCoefficient;
    buildTimeLimit = property.buildTimeLimit;
    time = 0.0;
    timer.start();
  }
  ~BuildTimeController() {
    property.edgeSizeForCreation = edgeSizeSave;
    property.insertionRadiusCoefficient = insertionRadiusCoefficientSave;
  }
  void adjustEdgeSize(size_t c) {
    if (buildTimeLimit > 0.0 && count <= c) {
      timer.stop();
      double estimatedTime = time + timer.time / interval * (noOfInsertedObjects - count);
      estimatedTime /= 60 * 60;	// hour
      const size_t edgeInterval = 5;
      const int minimumEdge = 5;
      const float radiusInterval = 0.02;
      if (estimatedTime > buildTimeLimit) {
	if (property.insertionRadiusCoefficient - radiusInterval >= 1.0) {
	  property.insertionRadiusCoefficient -= radiusInterval;
	} else {
	  property.edgeSizeForCreation -= edgeInterval;
	  if (property.edgeSizeForCreation < minimumEdge) {
	    property.edgeSizeForCreation = minimumEdge;
	  }
	}
      }
      time += timer.time;
      count += interval;
      timer.start();
    }
  }

  size_t	noOfInsertedObjects;
  size_t	interval;
  size_t	count ;
  size_t	edgeSizeSave;
  double	insertionRadiusCoefficientSave;
  Timer		timer;
  double	time;
  double	buildTimeLimit;
  NeighborhoodGraph::Property &property;
};

void 
GraphAndTreeIndex::createTreeIndex() 
{
  ObjectRepository &fr = GraphIndex::objectSpace->getRepository();
  for (size_t id = 0; id < fr.size(); id++){
    if (id % 100000 == 0) {
      cerr << " Processed id=" << id << endl;
    }
    if (fr.isEmpty(id)) {
      continue;
    }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    Object *f = GraphIndex::objectSpace->allocateObject(*fr[id]);
    DVPTree::InsertContainer tiobj(*f, id);
#else
    DVPTree::InsertContainer tiobj(*fr[id], id);
#endif
    try {
      DVPTree::insert(tiobj);
    } catch (Exception &err) {
      cerr << "GraphAndTreeIndex::createTreeIndex: Warning. ID=" << id << ":";
      cerr << err.what() << " continue.." << endl;
    }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    GraphIndex::objectSpace->deleteObject(f);
#endif
  }
}

void 
NGT::GraphIndex::loadIndex(const string &ifile, bool readOnly) {
  objectSpace->deserialize(ifile + "/obj");
#ifdef NGT_GRAPH_READ_ONLY_GRAPH
  if (readOnly && property.indexType == NGT::Index::Property::IndexType::Graph) {
    GraphIndex::NeighborhoodGraph::loadSearchGraph(ifile);
  } else {
    ifstream isg(ifile + "/grp");
    repository.deserialize(isg);
  }
#else
  ifstream isg(ifile + "/grp");
  repository.deserialize(isg);
#endif
}

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
NGT::GraphIndex::GraphIndex(const string &allocator, bool rdonly):readOnly(rdonly) {
  NGT::Property prop;
  prop.load(allocator);
  if (prop.databaseType != NGT::Index::Property::DatabaseType::MemoryMappedFile) {
    NGTThrowException("GraphIndex: Cannot open. Not memory mapped file type.");
  }
  initialize(allocator, prop);
#ifdef NGT_GRAPH_READ_ONLY_GRAPH
  searchUnupdatableGraph = NeighborhoodGraph::Search::getMethod(prop.distanceType, prop.objectType,
								objectSpace->getRepository().size());
#endif
}

NGT::GraphAndTreeIndex::GraphAndTreeIndex(const string &allocator, NGT::Property &prop):GraphIndex(allocator, prop) {
  initialize(allocator, prop.treeSharedMemorySize);
}

void 
NGT::GraphIndex::initialize(const string &allocator, NGT::Property &prop) {
  constructObjectSpace(prop);
  repository.open(allocator + "/grp", prop.graphSharedMemorySize);
  objectSpace->open(allocator + "/obj", prop.objectSharedMemorySize);
  setProperty(prop);
}
#else // NGT_SHARED_MEMORY_ALLOCATOR
NGT::GraphIndex::GraphIndex(const string &database, bool rdOnly):readOnly(rdOnly) {
  NGT::Property prop;
  prop.load(database);
  if (prop.databaseType != NGT::Index::Property::DatabaseType::Memory) {
    NGTThrowException("GraphIndex: Cannot open. Not memory type.");
  }
  assert(prop.dimension != 0);
  initialize(prop);
  loadIndex(database, readOnly);
#ifdef NGT_GRAPH_READ_ONLY_GRAPH
  if (prop.searchType == "Large") {
    searchUnupdatableGraph = NeighborhoodGraph::Search::getMethod(prop.distanceType, prop.objectType, 10000000);
  } else if (prop.searchType == "Small") {
    searchUnupdatableGraph = NeighborhoodGraph::Search::getMethod(prop.distanceType, prop.objectType, 0);
  } else {
    searchUnupdatableGraph = NeighborhoodGraph::Search::getMethod(prop.distanceType, prop.objectType,
                                                                  objectSpace->getRepository().size());
  }
#endif
}
#endif

void
GraphIndex::createIndex()
{
  GraphRepository &anngRepo = repository;
  ObjectRepository &fr = objectSpace->getRepository();
  size_t	pathAdjustCount = property.pathAdjustmentInterval;
  NGT::ObjectID id = 1;
  size_t count = 0;
  BuildTimeController buildTimeController(*this, NeighborhoodGraph::property);
  for (; id < fr.size(); id++) {
    if (id < anngRepo.size() && anngRepo[id] != 0) {
      continue;
    }
    insert(id);
    buildTimeController.adjustEdgeSize(++count);
    if (pathAdjustCount > 0 && pathAdjustCount <= id) {
      GraphReconstructor::adjustPathsEffectively(static_cast<GraphIndex&>(*this));
      pathAdjustCount += property.pathAdjustmentInterval;
    }
  }
}

size_t
searchMultipleQueryForCreation(GraphIndex &neighborhoodGraph, 
				 NGT::ObjectID &id, 
				 CreateIndexJob &job, 
				 CreateIndexThreadPool &threads)
{
  ObjectRepository &repo = neighborhoodGraph.objectSpace->getRepository();
  GraphRepository &anngRepo = neighborhoodGraph.repository;
  size_t cnt = 0;
  for (; id < repo.size(); id++) {
    if (repo[id] == 0) {
      continue;
    }
    if (neighborhoodGraph.NeighborhoodGraph::property.graphType != NeighborhoodGraph::GraphTypeBKNNG) {
      if (id < anngRepo.size() && anngRepo[id] != 0) {
	continue;
      }
    }
    job.id = id;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    job.object = neighborhoodGraph.objectSpace->allocateObject(*repo[id]);
#else
    job.object = repo[id];
#endif
    job.batchIdx = cnt;
    threads.pushInputQueue(job);
    cnt++;
    if (cnt >= (size_t)neighborhoodGraph.NeighborhoodGraph::property.batchSizeForCreation) {
      id++;
      break;
    }
  } // for
  return cnt;
}

void
insertMultipleSearchResults(GraphIndex &neighborhoodGraph, 
			    CreateIndexThreadPool::OutputJobQueue &output, 
			    size_t dataSize)
{
  // compute distances among all of the resultant objects
  if (neighborhoodGraph.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeANNG ||
      neighborhoodGraph.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeIANNG ||
      neighborhoodGraph.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeONNG ||
      neighborhoodGraph.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeDNNG) {
    // This processing occupies about 30% of total indexing time when batch size is 200.
    // Only initial batch objects should be connected for each other.
    // The number of nodes in the graph is checked to know whether the batch is initial.
    //size_t size = NeighborhoodGraph::property.edgeSizeForCreation;
    size_t size = neighborhoodGraph.NeighborhoodGraph::property.edgeSizeForCreation;
    // add distances from a current object to subsequence objects to imitate of sequential insertion.

    sort(output.begin(), output.end());	// sort by batchIdx

    for (size_t idxi = 0; idxi < dataSize; idxi++) {
      // add distances
      ObjectDistances &objs = *output[idxi].results;
      for (size_t idxj = 0; idxj < idxi; idxj++) {
	ObjectDistance	r;
	r.distance = neighborhoodGraph.objectSpace->getComparator()(*output[idxi].object, *output[idxj].object);
	r.id = output[idxj].id;
	objs.push_back(r);
      }
      // sort and cut excess edges	    
      std::sort(objs.begin(), objs.end());
      if (objs.size() > size) {
	objs.resize(size);
      }
    } // for (size_t idxi ....
  } // if (neighborhoodGraph.graphType == NeighborhoodGraph::GraphTypeUDNNG)
  // insert resultant objects into the graph as edges
  for (size_t i = 0; i < dataSize; i++) {
    CreateIndexJob &gr = output[i];
    if ((*gr.results).size() == 0) {
    }
    if (static_cast<int>(gr.id) > neighborhoodGraph.NeighborhoodGraph::property.edgeSizeForCreation &&
	static_cast<int>(gr.results->size()) < neighborhoodGraph.NeighborhoodGraph::property.edgeSizeForCreation) {
      cerr << "CreateIndex: The specified number of edges could not be acquired, because the pruned parameter [-S] might be set." << endl;
      cerr << "  The node id=" << gr.id << endl;
      cerr << "  The number of edges for the node=" << gr.results->size() << endl;
      cerr << "  The pruned parameter (edgeSizeForSearch [-S])=" << neighborhoodGraph.NeighborhoodGraph::property.edgeSizeForSearch << endl;
    }
    neighborhoodGraph.insertNode(gr.id, *gr.results);
  }
}

void 
GraphIndex::createIndex(size_t threadPoolSize) 
{
  if (threadPoolSize <= 1) {
    createIndex();
  } else {
    Timer		timer;
    size_t	timerInterval = 100000;
    size_t	timerCount = timerInterval;
    size_t	count = 0;
    timer.start();

    size_t	pathAdjustCount = property.pathAdjustmentInterval;
    CreateIndexThreadPool threads(threadPoolSize);
    CreateIndexSharedData sd(*this);

    threads.setSharedData(&sd);
    threads.create();
    CreateIndexThreadPool::OutputJobQueue &output = threads.getOutputJobQueue();

    BuildTimeController buildTimeController(*this, NeighborhoodGraph::property);

    try {
      CreateIndexJob job;
      NGT::ObjectID id = 1;
      for (;;) {
	// search for the nearest neighbors
	size_t cnt = searchMultipleQueryForCreation(*this, id, job, threads);
	if (cnt == 0) {
	  break;
	}
	// wait for the completion of the search
	threads.waitForFinish();
	if (output.size() != cnt) {
	  cerr << "NNTGIndex::insertGraphIndexByThread: Warning!! Thread response size is wrong." << endl;
	  cnt = output.size();
	}
	// insertion
	insertMultipleSearchResults(*this, output, cnt);

	while (!output.empty()) {
	  delete output.front().results;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	  GraphIndex::objectSpace->deleteObject(output.front().object);
#endif
	  output.pop_front();
	}

	count += cnt;
	if (timerCount <= count) {
	  timer.stop();
	  cerr << "Processed " << timerCount << " time= " << timer << endl;
	  timerCount += timerInterval;
	  timer.start();
	}
	buildTimeController.adjustEdgeSize(count);
	if (pathAdjustCount > 0 && pathAdjustCount <= count) {
	  GraphReconstructor::adjustPathsEffectively(static_cast<GraphIndex&>(*this));
	  pathAdjustCount += property.pathAdjustmentInterval;
	}
      }
    } catch(Exception &err) {
      threads.terminate();
      throw err;
    }
    threads.terminate();
  }

}

void 
GraphAndTreeIndex::createIndex(size_t threadPoolSize) 
{
  assert(threadPoolSize > 0);

  Timer	timer;
  size_t	timerInterval = 100000;
  size_t	timerCount = timerInterval;
  size_t	count = 0;
  timer.start();

  size_t	pathAdjustCount = property.pathAdjustmentInterval;
  CreateIndexThreadPool threads(threadPoolSize);

  CreateIndexSharedData sd(*this);

  threads.setSharedData(&sd);
  threads.create();
  CreateIndexThreadPool::OutputJobQueue &output = threads.getOutputJobQueue();

  BuildTimeController buildTimeController(*this, NeighborhoodGraph::property);

  try {
    CreateIndexJob job;
    NGT::ObjectID id = 1;
    for (;;) {
      size_t cnt = searchMultipleQueryForCreation(*this, id, job, threads);

      if (cnt == 0) {
	break;
      }
      threads.waitForFinish();

      if (output.size() != cnt) {
	cerr << "NNTGIndex::insertGraphIndexByThread: Warning!! Thread response size is wrong." << endl;
	cnt = output.size();
      }

      insertMultipleSearchResults(*this, output, cnt);

      for (size_t i = 0; i < cnt; i++) {
	CreateIndexJob &job = output[i];
	if (((job.results->size() > 0) && ((*job.results)[0].distance != 0.0)) ||
	    (job.results->size() == 0)) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	  Object *f = GraphIndex::objectSpace->allocateObject(*job.object);
	  DVPTree::InsertContainer tiobj(*f, job.id);
#else
	  DVPTree::InsertContainer tiobj(*job.object, job.id);
#endif
	  try {
	    DVPTree::insert(tiobj);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	    GraphIndex::objectSpace->deleteObject(f);
#endif
	  } catch (Exception &err) {
	    cerr << "NGT::createIndex: Fatal error. ID=" << job.id << ":";
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	    GraphIndex::objectSpace->deleteObject(f);
#endif
	    if (NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeKNNG) {
	      cerr << err.what() << " continue.." << endl;
	    } else {
	      throw err;
	    }
	  }
	}
      } // for

      while (!output.empty()) {
	delete output.front().results;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	GraphIndex::objectSpace->deleteObject(output.front().object);
#endif
	output.pop_front();
      }

      count += cnt;
      if (timerCount <= count) {
	timer.stop();
	cerr << "Processed " << timerCount << " objects. time= " << timer << endl;
	timerCount += timerInterval;
	timer.start();
      }
      buildTimeController.adjustEdgeSize(count);
      if (pathAdjustCount > 0 && pathAdjustCount <= count) {
	GraphReconstructor::adjustPathsEffectively(static_cast<GraphIndex&>(*this));
	pathAdjustCount += property.pathAdjustmentInterval;
      }
    }
  } catch(Exception &err) {
    threads.terminate();
    throw err;
  }
  threads.terminate();
}


void 
GraphAndTreeIndex::createIndex(const vector<pair<NGT::Object*, size_t> > &objects, 
			       vector<InsertionResult> &ids, 
			       double range, size_t threadPoolSize)
{
  Timer		timer;
  size_t	timerInterval = 100000;
  size_t	timerCount = timerInterval;
  size_t	count = 0;
  timer.start();
  if (threadPoolSize <= 0) {
    cerr << "Not implemented!!" << endl;
    abort();
  } else {
    CreateIndexThreadPool threads(threadPoolSize);
    CreateIndexSharedData sd(*this);
    threads.setSharedData(&sd);
    threads.create();
    CreateIndexThreadPool::OutputJobQueue &output = threads.getOutputJobQueue();
    try {
      CreateIndexJob job;
      size_t	idx = 0;
      for (;;) {
	size_t	cnt = 0;
	{
	  for (; idx < objects.size(); idx++) {
	    if (objects[idx].first == 0) {
	      ids.push_back(InsertionResult());
	      continue;
	    }
	    job.id = 0;
	    job.results = 0;
	    job.object = objects[idx].first;
	    job.batchIdx = ids.size();
	    // insert an empty entry to prepare.
	    ids.push_back(InsertionResult(job.id, false, 0.0));
	    threads.pushInputQueue(job);
	    cnt++;
	    if (cnt >= (size_t)NeighborhoodGraph::property.batchSizeForCreation) {
	      idx++;
	      break;
	    }
	  } 
	}
	if (cnt == 0) {
	  break;
	}
	threads.waitForFinish();
	if (output.size() != cnt) {
	  cerr << "NNTGIndex::insertGraphIndexByThread: Warning!! Thread response size is wrong." << endl;
	  cnt = output.size();
	}
	{
	  // This processing occupies about 30% of total indexing time when batch size is 200.
	  // Only initial batch objects should be connected for each other.
	  // The number of nodes in the graph is checked to know whether the batch is initial.
	  size_t size = NeighborhoodGraph::property.edgeSizeForCreation;
	  // add distances from a current object to subsequence objects to imitate of sequential insertion.

	  sort(output.begin(), output.end());	
	  for (size_t idxi = 0; idxi < cnt; idxi++) {
	    // add distances
	    ObjectDistances &objs = *output[idxi].results;
	    for (size_t idxj = 0; idxj < idxi; idxj++) {
	      if (output[idxj].id == 0) {
		// unregistered object
		continue;
	      }
	      ObjectDistance	r;
	      r.distance = GraphIndex::objectSpace->getComparator()(*output[idxi].object, *output[idxj].object);
	      r.id = output[idxj].id;
	      objs.push_back(r);
	    }
	    // sort and cut excess edges	    
	    std::sort(objs.begin(), objs.end());
	    if (objs.size() > size) {
	      objs.resize(size);
	    }
	    if ((objs.size() > 0) && (range < 0.0 || ((double)objs[0].distance <= range + FLT_EPSILON))) {
	      // The line below was replaced by the line above to consider EPSILON for float comparison. 170702
	      // if ((objs.size() > 0) && (range < 0.0 || (objs[0].distance <= range))) {
	      // An identical or similar object already exits
	      ids[output[idxi].batchIdx].identical = true;
	      ids[output[idxi].batchIdx].id = objs[0].id;
	      ids[output[idxi].batchIdx].distance = objs[0].distance;
	      output[idxi].id = 0;
	    } else {
	      assert(output[idxi].id == 0);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	      PersistentObject *obj = GraphIndex::objectSpace->allocatePersistentObject(*output[idxi].object);
	      output[idxi].id = GraphIndex::objectSpace->insert(obj);
#else
	      output[idxi].id = GraphIndex::objectSpace->insert(output[idxi].object);
#endif
	      ids[output[idxi].batchIdx].id = output[idxi].id;
	    }
	  } 
	}
	// insert resultant objects into the graph as edges
	for (size_t i = 0; i < cnt; i++) {
	  CreateIndexJob &job = output.front();
	  if (job.id != 0) {
	    if (property.indexType == NGT::Property::GraphAndTree) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	      Object *f = GraphIndex::objectSpace->allocateObject(*job.object);
	      DVPTree::InsertContainer tiobj(*f, job.id);
#else
	      DVPTree::InsertContainer tiobj(*job.object, job.id);
#endif
	      try {
		DVPTree::insert(tiobj);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
		GraphIndex::objectSpace->deleteObject(f);
#endif
	      } catch (Exception &err) {
		cerr << "NGT::createIndex: Fatal error. ID=" << job.id << ":" << err.what();
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
		GraphIndex::objectSpace->deleteObject(f);
#endif
		if (NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeKNNG) {
		  cerr << err.what() << " continue.." << endl;
		} else {
		  throw err;
		}
	      }
	    }
	    if (((*job.results).size() == 0) && (job.id != 1)) {
	      cerr  << "insert warning!! No searched nodes!. If the first time, no problem. " << job.id << endl;
	    }
	    GraphIndex::insertNode(job.id, *job.results);
	  } 
	  if (job.results != 0) {
	    delete job.results;
	  }
	  output.pop_front();
	}
	
	count += cnt;
	if (timerCount <= count) {
	  timer.stop();
	  cerr << "Processed " << timerCount << " time= " << timer << endl;
	  timerCount += timerInterval;
	  timer.start();
	}
      }
    } catch(Exception &err) {
      cerr << "thread terminate!" << endl;
      threads.terminate();
      throw err;
    }
    threads.terminate();
  }
}

static bool 
findPathAmongIdenticalObjects(GraphAndTreeIndex &graph, size_t srcid, size_t dstid) {
  stack<size_t> nodes;
  unordered_set<size_t> done;
  nodes.push(srcid);
  while (!nodes.empty()) {
    auto tid = nodes.top();
    nodes.pop();
    done.insert(tid);
    GraphNode &node = *graph.GraphIndex::getNode(tid);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    for (auto i = node.begin(graph.repository.allocator); i != node.end(graph.GraphIndex::repository.allocator); ++i) {
#else
    for (auto i = node.begin(); i != node.end(); ++i) {
#endif
      if ((*i).distance != 0.0) {
	break;
      }
      if ((*i).id == dstid) {
	return true;
      }
      if (done.count((*i).id) == 0) {
        nodes.push((*i).id);
      }
    }
  }
  return false;
}

bool 
GraphAndTreeIndex::verify(vector<uint8_t> &status, bool info, char mode) {
  bool valid = GraphIndex::verify(status, info);
  if (!valid) {
    cerr << "The graph or object is invalid!" << endl;
  }
  bool treeValid = DVPTree::verify(GraphIndex::objectSpace->getRepository().size(), status);
  if (!treeValid) {
    cerr << "The tree is invalid" << endl;
  }
  valid = valid && treeValid;
  // status: tree|graph|object
  cerr << "Started checking consistency..." << endl;
  for (size_t id = 1; id < status.size(); id++) {
    if (id % 100000 == 0) {
      cerr << "The number of processed objects=" << id << endl;
    }
    if (status[id] != 0x00 && status[id] != 0x07) {
      if (status[id] == 0x03) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	NGT::Object *po = GraphIndex::objectSpace->allocateObject(*GraphIndex::getObjectRepository().get(id));
	NGT::SearchContainer sc(*po);
#else
	NGT::SearchContainer sc(*GraphIndex::getObjectRepository().get(id));
#endif
	NGT::ObjectDistances objects;
	sc.setResults(&objects);
	sc.id = 0;
	sc.radius = 0.0;
	sc.explorationCoefficient = 1.1;
	sc.edgeSize = 0;
	ObjectDistances	seeds;
	seeds.push_back(ObjectDistance(id, 0.0));
	objects.clear();
	try {
	  GraphIndex::search(sc, seeds);
	} catch(Exception &err) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	  GraphIndex::objectSpace->deleteObject(po);
#endif
	  cerr << "Fatal Error!: Cannot search! " << err.what() << endl;
	  objects.clear();
	}
	size_t n = 0;
	bool registeredIdenticalObject = false;
	for (; n < objects.size(); n++) {
	  if (objects[n].id != id && status[objects[n].id] == 0x07) {
	    registeredIdenticalObject = true;
	    break;
	  }
	}
	if (!registeredIdenticalObject) {
	  if (info) {
	    cerr << "info: not found the registered same objects. id=" << id << " size=" << objects.size() << endl;
	  }
	  sc.id = 0;
	  sc.radius = FLT_MAX;
	  sc.explorationCoefficient = 1.2;
	  sc.edgeSize = 0;
	  sc.size = objects.size() < 100 ? 100 : objects.size() * 2;
	  ObjectDistances	seeds;
	  seeds.push_back(ObjectDistance(id, 0.0));
	  objects.clear();
	  try {
	    GraphIndex::search(sc, seeds);
	  } catch(Exception &err) {
	    cerr << "Fatal Error!: Cannot search! " << err.what() << endl;
	    objects.clear();
	  }
	  registeredIdenticalObject = false;
	  for (n = 0; n < objects.size(); n++) {
	    if (objects[n].distance != 0.0) break;
	    if (objects[n].id != id && status[objects[n].id] == 0x07) {
	      registeredIdenticalObject = true;
	      if (info) {
		cerr << "info: found by using mode accurate search. " << objects[n].id << endl;
	      }
	      break;
	    }
	  }
	}
	if (!registeredIdenticalObject && mode != 's') {
	  if (info) {
	    cerr << "info: not found by using more accurate search." << endl;
	  }
	  sc.id = 0;
	  sc.radius = 0.0;
	  sc.explorationCoefficient = 1.1;
	  sc.edgeSize = 0;
	  sc.size = SIZE_MAX;
	  objects.clear();
	  linearSearch(sc);
	  n = 0;
	  registeredIdenticalObject = false;
	  for (; n < objects.size(); n++) {
	    if (objects[n].distance != 0.0) break;
	    if (objects[n].id != id && status[objects[n].id] == 0x07) {
	      registeredIdenticalObject = true;
	      if (info) {
		cerr << "info: found by using linear search. " << objects[n].id << endl;
	      }
	      break;
	    }
	  }
	}
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	GraphIndex::objectSpace->deleteObject(po);
#endif
	if (registeredIdenticalObject) {
	  if (info) {
	    cerr << "Info ID=" << id << ":" << static_cast<int>(status[id]) << endl;
	    cerr << "  found the valid same objects. " << objects[n].id << endl;
	  }
	  GraphNode &fromNode = *GraphIndex::getNode(id);
	  bool fromFound = false;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	  for (auto i = fromNode.begin(GraphIndex::repository.allocator); i != fromNode.end(GraphIndex::repository.allocator); ++i) {
#else
	  for (auto i = fromNode.begin(); i != fromNode.end(); ++i) {
#endif
	    if ((*i).id == objects[n].id) {
	      fromFound = true;
	    }
	  }
	  GraphNode &toNode = *GraphIndex::getNode(objects[n].id);
	  bool toFound = false;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	  for (auto i = toNode.begin(GraphIndex::repository.allocator); i != toNode.end(GraphIndex::repository.allocator); ++i) {
#else
	  for (auto i = toNode.begin(); i != toNode.end(); ++i) {
#endif
	    if ((*i).id == id) {
	      toFound = true;
	    }
	  }
	  if (!fromFound || !toFound) {
	    if (info) {
	      if (!fromFound && !toFound) {
		cerr << "Warning no undirected edge between " << id << "(" << fromNode.size() << ") and " 
		     << objects[n].id << "(" << toNode.size() << ")." << endl;
	      } else if (!fromFound) {
		cerr << "Warning no directed edge from " << id << "(" << fromNode.size() << ") to " 
		     << objects[n].id << "(" << toNode.size() << ")." << endl;
	      } else if (!toFound) {
		cerr << "Warning no reverse directed edge from " << id << "(" << fromNode.size() << ") to " 
		     << objects[n].id << "(" << toNode.size() << ")." << endl;
	      }
	    }
	    if (!findPathAmongIdenticalObjects(*this, id, objects[n].id)) {
	      cerr << "Warning no path from " << id << " to " << objects[n].id << endl;
	    }
	    if (!findPathAmongIdenticalObjects(*this, objects[n].id, id)) {
	      cerr << "Warning no reverse path from " << id << " to " << objects[n].id << endl;
	    }
	  }
	} else {
	  if (mode == 's') {
	    cerr << "Warning: not found the valid same object, but not try to use linear search." << endl;
	    cerr << "Error! ID=" << id << ":" << static_cast<int>(status[id]) << endl;
	  } else {
	    cerr << "Warning: not found the valid same object even by using linear search." << endl;
	    cerr << "Error! ID=" << id << ":" << static_cast<int>(status[id]) << endl;
	    valid = false;
	  }
	}
      } else if (status[id] == 0x01) {
	  if (info) {
	    cerr << "Warning! ID=" << id << ":" << static_cast<int>(status[id]) << endl;
	    cerr << "  not inserted into the indexes" << endl;
	  }
      } else {
	  cerr << "Error! ID=" << id << ":" << static_cast<int>(status[id]) << endl;
	  valid = false;
      }
    }
  }
  return valid;
}



