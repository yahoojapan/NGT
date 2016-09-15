
//
// Copyright (C) 2015-2016 Yahoo! JAPAN Research
//
// This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
// To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
//

#include	"NGT/defines.h"
#include	"NGT/Common.h"
#include	"NGT/ObjectSpace.h"
#include	"NGT/Index.h"
#include	"NGT/Thread.h"

using namespace NGT;


class CreateIndexJob {
public:
  CreateIndexJob() {}
  CreateIndexJob &operator=(CreateIndexJob &d) {
    id = d.id;
    results = d.results;
    object = d.object;
    return *this;
  }
  NGT::ObjectID		id;
  NGT::Object		*object;	// this will be a node of the graph later.
  NGT::ObjectDistances	*results;
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
      cerr << "CreateIndex::search:popFront " << err.what() << endl;
      break;
    }
    ObjectDistances *rs = new ObjectDistances;
    Object &obj = *job.object;
    if (graphIndex.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeKNNG) {
      graphIndex.searchForKNNGInsertion(obj, job.id, *rs);
    } else {
      graphIndex.searchForNNGInsertion(obj, *rs);
    }
    job.results = rs;
    poolThread.getOutputJobQueue().pushBack(job);
  }

  return 0;

}

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
GraphIndex::createIndex()
{
  GraphRepository &anngRepo = repository;
  ObjectRepository &fr = objectSpace->getRepository();
  NGT::ObjectID id = 1;
  for (; id < fr.size(); id++) {
    if (id < anngRepo.size() && anngRepo[id] != 0) {
      continue;
    }
    insert(id);
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
      neighborhoodGraph.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeDNNG) {
    // This processing occupies about 30% of total indexing time when batch size is 200.
    // Only initial batch objects should be connected for each other.
    // The number of nodes in the graph is checked to know whether the batch is initial.
    size_t size = neighborhoodGraph.NeighborhoodGraph::property.edgeSizeForCreation;
    vector<Distance> distances;
    distances.reserve(dataSize * dataSize);
    for (size_t i = 0; i < dataSize; i++) {
      for (size_t j = 0; j < dataSize; j++) {
	Distance d;
	if (i == j) {
	  d = 0;
	  distances.push_back(d);
	  continue;
	} else if (i > j) {
	  d = distances[j * dataSize + i];
	  distances.push_back(d);
	  continue;
	} else {
	  d = neighborhoodGraph.objectSpace->getComparator()(*output[i].object, *output[j].object);
	  distances.push_back(d);
	}
      }
    }
    // add distances from a current object to subsequence objects to imitate of sequential insertion.
    for (size_t i = 0; i < dataSize; i++) {
      for (size_t j = i + 1; j < dataSize; j++) {
	ObjectDistance	r;
	r.distance = distances[i * dataSize + j];
	r.id = output[i].id;
	output[j].results->push_back(r);
      }
    }
    // sort and cut excess edges
    for (size_t i = 0; i < dataSize; i++) {
      ObjectDistances &objs = *output[i].results;
      std::sort(objs.begin(), objs.end());
      if (objs.size() > size) {
	objs.resize(size);
      }
    }
  } // if (neighborhoodGraph.graphType == NeighborhoodGraph::GraphTypeUDNNG)
  // insert resultant objects into the graph as edges
  for (size_t i = 0; i < dataSize; i++) {
    CreateIndexJob &gr = output.front();
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    neighborhoodGraph.objectSpace->deleteObject(gr.object);
#endif
    if ((*gr.results).size() == 0) {
    }
    neighborhoodGraph.insertNode(gr.id, *gr.results);
    delete gr.results;
    output.pop_front();
  }
}

void 
GraphIndex::createIndex(size_t threadPoolSize) 
{
  Timer		timer;
  size_t	timerInterval = 100000;
  size_t	timerCount = timerInterval;
  size_t	count = 0;
  timer.start();
  if (threadPoolSize <= 1) {
    createIndex();
  } else {
    CreateIndexThreadPool threads(threadPoolSize);
    CreateIndexSharedData sd(*this);

    threads.setSharedData(&sd);
    threads.create();
    CreateIndexThreadPool::OutputJobQueue &output = threads.getOutputJobQueue();
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
	count += cnt;
	if (timerCount <= count) {
	  timer.stop();
	  cerr << "Processed " << timerCount << " time= " << timer << endl;
	  timerCount += timerInterval;
	  timer.start();
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
  Timer		timer;
  size_t	timerInterval = 100000;
  size_t	timerCount = timerInterval;
  size_t	count = 0;
  timer.start();
  if (threadPoolSize <= 1) {
    GraphIndex::createIndex();
  } else {
    CreateIndexThreadPool threads(threadPoolSize);

    CreateIndexSharedData sd(*this);

    threads.setSharedData(&sd);
    threads.create();
    CreateIndexThreadPool::OutputJobQueue &output = threads.getOutputJobQueue();
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
	insertMultipleSearchResults(*this, output, cnt);

	count += cnt;
	if (timerCount <= count) {
	  timer.stop();
	  cerr << "Processed " << timerCount << " time= " << timer << endl;
	  timerCount += timerInterval;
	  timer.start();
	}
      }
    } catch(Exception &err) {
      threads.terminate();
      throw err;
    }
    threads.terminate();
  }

}



