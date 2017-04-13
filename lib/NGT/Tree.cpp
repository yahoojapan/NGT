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

#include	"NGT/defines.h"

#include	"NGT/Tree.h"
#include	"NGT/Node.h"

#include	<vector>

using namespace NGT;

void
DVPTree::insert(InsertContainer &iobj) {
  SearchContainer q(iobj.object);
  q.mode = SearchContainer::SearchLeaf;
  q.vptree = this;
  q.radius = 0.0;

  search(q);

  iobj.vptree = this;

  assert(q.nodeID.getType() == Node::ID::Leaf);
  LeafNode *ln = (LeafNode*)getNode(q.nodeID);
  insert(iobj, ln);

  return;
}

void
DVPTree::insert(InsertContainer &iobj,  LeafNode *leafNode) 
{
  LeafNode &leaf = *leafNode;
  size_t fsize = leaf.getObjectSize();
  if (fsize != 0) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    Distance d = objectSpace->getComparator()(iobj.object, leaf.getPivot(*objectSpace));
#else
    Distance d = objectSpace->getComparator()(iobj.object, leaf.getPivot());
#endif

    for (size_t i = 0; i < fsize; i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      if (leaf.getObjectIDs(leafNodes.allocator)[i].distance == d) {
#else
      if (leaf.getObjectIDs()[i].distance == d) {
#endif
	Distance idd = 0.0;
	ObjectID loid;
        try {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
	  loid = leaf.getObjectIDs(leafNodes.allocator)[i].id;
	  idd = objectSpace->getComparator()(iobj.object, *getObjectRepository().get(loid));
#else
	  loid = leaf.objectIDs[i].id;
	  idd = objectSpace->getComparator()(iobj.object, *getObjectRepository().get(loid));
#endif
        } catch (Exception &e) {
          stringstream msg;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
          msg << "LeafNode::insert: Cannot find object which belongs to a leaf node. id="
              << leaf.getObjectIDs(leafNodes.allocator)[i].id << ":" << e.what() << endl;
#else
          msg << "LeafNode::insert: Cannot find object which belongs to a leaf node. id="
              << leaf.getObjectIDs()[i].id << ":" << e.what() << endl;
#endif
          NGTThrowException(msg.str());
        }
        if (idd == 0.0) {
	  if (loid == iobj.id) {
	    stringstream msg;
	    msg << "DVPTree::insert:already existed. " << iobj.id;
	    NGTThrowException(msg);
	  }
	  return;
        }
      }
    }
  }

  if (leaf.getObjectSize() >= leafObjectsSize) {
    split(iobj, leaf);
  } else {
    insertObject(iobj, leaf);
  }

  return;
}
Node::ID 
DVPTree::split(InsertContainer &iobj, LeafNode &leaf)
{
  Node::Objects *fs = getObjects(leaf, iobj);
  int pv = DVPTree::MaxVariance;
  switch (splitMode) {
  case DVPTree::MaxVariance:
    pv = LeafNode::selectPivotByMaxVariance(iobj, *fs);
    break;
  case DVPTree::MaxDistance:
    pv = LeafNode::selectPivotByMaxDistance(iobj, *fs);
    break;
  }
  LeafNode::splitObjects(iobj, *fs, pv);
  Node::ID nid = recombineNodes(iobj, *fs, leaf);
  delete fs;

  return nid;
}

Node::ID
DVPTree::recombineNodes(InsertContainer &ic, Node::Objects &fs, LeafNode &leaf)
{
  LeafNode *ln[internalChildrenSize];
  Node::ID targetParent = leaf.parent;
  Node::ID targetId = leaf.id;
  ln[0] = &leaf;
  ln[0]->objectSize = 0;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  for (size_t i = 1; i < internalChildrenSize; i++) {
    ln[i] = new(leafNodes.allocator) LeafNode(leafNodes.allocator);
  }
#else
  for (size_t i = 1; i < internalChildrenSize; i++) {
    ln[i] = new LeafNode;
  }
#endif
  InternalNode *in = createInternalNode();
  Node::ID inid = in->id;
  try {
    if (targetParent.getID() != 0) {
      InternalNode &pnode = *(InternalNode*)getNode(targetParent);
      for (size_t i = 0; i < internalChildrenSize; i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	if (pnode.getChildren(internalNodes.allocator)[i] == targetId) {
	  pnode.getChildren(internalNodes.allocator)[i] = inid;
#else
	if (pnode.getChildren()[i] == targetId) {
	  pnode.getChildren()[i] = inid;
#endif
	  break;
	}
      }
    }
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    in->setPivot(*getObjectRepository().get(fs[0].id), *objectSpace, internalNodes.allocator);
#else
    in->setPivot(*getObjectRepository().get(fs[0].id), *objectSpace);
#endif

    in->parent = targetParent;

    int fsize = fs.size();
    int cid = fs[0].clusterID;
#ifdef NGT_NODE_USE_VECTOR
    LeafNode::ObjectIDs fid;
    fid.id = fs[0].id;
    fid.distance = 0.0;
    ln[cid]->objectIDs.push_back(fid);
#else
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    ln[cid]->getObjectIDs(leafNodes.allocator)[ln[cid]->objectSize].id = fs[0].id;
    ln[cid]->getObjectIDs(leafNodes.allocator)[ln[cid]->objectSize++].distance = 0.0;
#else
    ln[cid]->getObjectIDs()[ln[cid]->objectSize].id = fs[0].id;
    ln[cid]->getObjectIDs()[ln[cid]->objectSize++].distance = 0.0;
#endif
#endif
    if (fs[0].leafDistance == Node::Object::Pivot) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      ln[cid]->setPivot(*getObjectRepository().get(fs[0].id), *objectSpace, leafNodes.allocator);
#else
      ln[cid]->setPivot(*getObjectRepository().get(fs[0].id), *objectSpace);
#endif
    } else {
      NGTThrowException("recombineNodes: internal error : illegal pivot.");
    }
    ln[cid]->parent = inid;
    for (int i = 1; i < fsize; i++) {
      int clusterID = fs[i].clusterID;
      Distance ld;
      if (fs[i].leafDistance == Node::Object::Pivot) {
        // pivot
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	ln[clusterID]->setPivot(*getObjectRepository().get(fs[i].id), *objectSpace, leafNodes.allocator);
#else
	ln[clusterID]->setPivot(*getObjectRepository().get(fs[i].id), *objectSpace);
#endif
        ld = 0.0;
      } else {
        ld = fs[i].leafDistance;
      }

#ifdef NGT_NODE_USE_VECTOR
      fid.id = fs[i].id;
      fid.distance = ld;
      ln[clusterID]->objectIDs.push_back(fid);
#else
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      ln[clusterID]->getObjectIDs(leafNodes.allocator)[ln[clusterID]->objectSize].id = fs[i].id;
      ln[clusterID]->getObjectIDs(leafNodes.allocator)[ln[clusterID]->objectSize++].distance = ld;
#else
      ln[clusterID]->getObjectIDs()[ln[clusterID]->objectSize].id = fs[i].id;
      ln[clusterID]->getObjectIDs()[ln[clusterID]->objectSize++].distance = ld;
#endif
#endif
      ln[clusterID]->parent = inid;
      if (clusterID != cid) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
        in->getBorders(internalNodes.allocator)[cid] = fs[i].distance;
#else
        in->getBorders()[cid] = fs[i].distance;
#endif
        cid = fs[i].clusterID;
      }
    }
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    in->getChildren(internalNodes.allocator)[0] = targetId;
#else
    in->getChildren()[0] = targetId;
#endif
    for (size_t i = 1; i < internalChildrenSize; i++) {
      insertNode(ln[i]);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      in->getChildren(internalNodes.allocator)[i] = ln[i]->id;
#else
      in->getChildren()[i] = ln[i]->id;
#endif
    }
  } catch(Exception &e) {
    throw e;
  }
  return inid;
}

void
DVPTree::insertObject(InsertContainer &ic, LeafNode &leaf) {
  if (leaf.getObjectSize() == 0) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    leaf.setPivot(*getObjectRepository().get(ic.id), *objectSpace, leafNodes.allocator);
#else
    leaf.setPivot(*getObjectRepository().get(ic.id), *objectSpace);
#endif
#ifdef NGT_NODE_USE_VECTOR
    LeafNode::ObjectIDs fid;
    fid.id = ic.id;
    fid.distance = 0;
    leaf.objectIDs.push_back(fid);
#else
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    leaf.getObjectIDs(leafNodes.allocator)[leaf.objectSize].id = ic.id;
    leaf.getObjectIDs(leafNodes.allocator)[leaf.objectSize++].distance = 0;
#else
    leaf.getObjectIDs()[leaf.objectSize].id = ic.id;
    leaf.getObjectIDs()[leaf.objectSize++].distance = 0;
#endif
#endif
  } else {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    Distance d = objectSpace->getComparator()(ic.object, leaf.getPivot(*objectSpace));
#else
    Distance d = objectSpace->getComparator()(ic.object, leaf.getPivot());
#endif

#ifdef NGT_NODE_USE_VECTOR
    LeafNode::ObjectIDs fid;
    fid.id = ic.id;
    fid.distance = d;
    leaf.objectIDs.push_back(fid);
    std::sort(leaf.objectIDs.begin(), leaf.objectIDs.end(), LeafNode::ObjectIDs());
#else
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    leaf.getObjectIDs(leafNodes.allocator)[leaf.objectSize].id = ic.id;
    leaf.getObjectIDs(leafNodes.allocator)[leaf.objectSize++].distance = d;
#else
    leaf.getObjectIDs()[leaf.objectSize].id = ic.id;
    leaf.getObjectIDs()[leaf.objectSize++].distance = d;
#endif
#endif
  }
}

Node::Objects *
DVPTree::getObjects(LeafNode &n, Container &iobj)
{
  int size = n.getObjectSize() + 1;

  Node::Objects *fs = new Node::Objects(size);
  for (size_t i = 0; i < n.getObjectSize(); i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    (*fs)[i].object = getObjectRepository().get(n.getObjectIDs(leafNodes.allocator)[i].id);
    (*fs)[i].id = n.getObjectIDs(leafNodes.allocator)[i].id;
#else
    (*fs)[i].object = getObjectRepository().get(n.getObjectIDs()[i].id);
    (*fs)[i].id = n.getObjectIDs()[i].id;
#endif
  }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  (*fs)[n.getObjectSize()].object = getObjectRepository().get(iobj.id);
#else
  (*fs)[n.getObjectSize()].object = &iobj.object;
#endif
  (*fs)[n.getObjectSize()].id = iobj.id;
  return fs;
}

void
DVPTree::removeEmptyNodes(InternalNode &inode) {

  int csize = internalChildrenSize;

  InternalNode *target = &inode;
  for(;;) {
    for (int i = 0; i < csize; i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      if (target->getChildren(internalNodes.allocator)[i].getType() == Node::ID::Internal) {
#else
      if (target->getChildren()[i].getType() == Node::ID::Internal) {
#endif
	return;
      }
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      LeafNode &ln = *static_cast<LeafNode*>(getNode(target->getChildren(internalNodes.allocator)[i]));
#else
      LeafNode &ln = *static_cast<LeafNode*>(getNode(target->getChildren()[i]));
#endif
      if (ln.getObjectSize() != 0) {
	return;
      }
    }

    for (int i = 0; i < csize; i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      removeNode(target->getChildren(internalNodes.allocator)[i]);
#else
      removeNode(target->getChildren()[i]);
#endif
    }
    if (target->parent.getID() == 0) {
      removeNode(target->id);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      LeafNode *root = new(leafNodes.allocator) LeafNode(leafNodes.allocator);
#else
      LeafNode *root = new LeafNode;
#endif
      insertNode(root);
      if (root->id.getID() != 1) {
	NGTThrowException("Root id Error");
      }
      return;
    }

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    LeafNode *ln = new(leafNodes.allocator) LeafNode(leafNodes.allocator);
#else
    LeafNode *ln = new LeafNode;
#endif
    ln->parent = target->parent;
    insertNode(ln);

    InternalNode &in = *(InternalNode*)getNode(ln->parent);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    in.updateChild(*this, target->id, ln->id, internalNodes.allocator);
#else
    in.updateChild(*this, target->id, ln->id);
#endif
    removeNode(target->id);
    target = &in;
  }

  return;
}


void
DVPTree::search(SearchContainer &sc, InternalNode &node, UncheckedNode &uncheckedNode)
{
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
  Distance d = objectSpace->getComparator()(sc.object, node.getPivot(*objectSpace));
#else
  Distance d = objectSpace->getComparator()(sc.object, node.getPivot());
#endif

  int bsize = internalChildrenSize - 1;

  vector<ObjectDistance> regions;

  ObjectDistance child;
  int mid;
  for (mid = 0; mid < bsize; mid++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    if (d < node.getBorders(internalNodes.allocator)[mid]) {
      if (d + sc.radius < node.getBorders(internalNodes.allocator)[mid]) {
#else
    if (d < node.getBorders()[mid]) {
      if (d + sc.radius < node.getBorders()[mid]) {
#endif
        child.id = mid;
        child.distance = 0.0;
        regions.push_back(child);
        break;
      } else {
        child.id = mid;
        child.distance = 0.0;
        regions.push_back(child);
        continue;
      }
    } else {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      if (d < node.getBorders(internalNodes.allocator)[mid] + sc.radius) {
#else
      if (d < node.getBorders()[mid] + sc.radius) {
#endif
        child.id = mid;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
        child.distance = d - node.getBorders(internalNodes.allocator)[mid];
#else
        child.distance = d - node.getBorders()[mid];
#endif
        regions.push_back(child);
        continue;
      } else {
        continue;
      }
    }
  }

  if (mid == bsize) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    if (d >= node.getBorders(internalNodes.allocator)[mid - 1]) {
#else
    if (d >= node.getBorders()[mid - 1]) {
#endif
      child.id = mid;
      child.distance = 0.0;
      regions.push_back(child);
    } else {
      child.id = mid;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      child.distance = node.getBorders(internalNodes.allocator)[mid - 1] - d;
#else
      child.distance = node.getBorders()[mid - 1] - d;
#endif
      regions.push_back(child);
    }
  }

  sort(regions.begin(), regions.end());

  vector<ObjectDistance>::iterator i;
  if (sc.mode == DVPTree::SearchContainer::SearchLeaf) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    if (node.getChildren(internalNodes.allocator)[regions.front().id].getType() == Node::ID::Leaf) {
      sc.nodeID.setRaw(node.getChildren(internalNodes.allocator)[regions.front().id].get());
#else
    if (node.getChildren()[regions.front().id].getType() == Node::ID::Leaf) {
      sc.nodeID.setRaw(node.getChildren()[regions.front().id].get());
#endif
      assert(uncheckedNode.empty());
    } else {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      uncheckedNode.push(node.getChildren(internalNodes.allocator)[regions.front().id]);
#else
      uncheckedNode.push(node.getChildren()[regions.front().id]);
#endif
    }
  } else {
    for (i = regions.begin(); i != regions.end(); i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      uncheckedNode.push(node.getChildren(internalNodes.allocator)[i->id]);
#else
      uncheckedNode.push(node.getChildren()[i->id]);
#endif
    }
  }
  
}

void
DVPTree::search(SearchContainer &so, LeafNode &node, UncheckedNode &uncheckedNode)
{
  DVPTree::SearchContainer &q = (DVPTree::SearchContainer&)so;

  if (node.getObjectSize() == 0) {
    return;
  }
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
  Distance pq = objectSpace->getComparator()(q.object, node.getPivot(*objectSpace));
#else
  Distance pq = objectSpace->getComparator()(q.object, node.getPivot());
#endif

  ObjectDistance r;

  for (size_t i = 0; i < node.getObjectSize(); i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    if ((node.getObjectIDs(leafNodes.allocator)[i].distance <= pq + q.radius) &&
        (node.getObjectIDs(leafNodes.allocator)[i].distance >= pq - q.radius)) {
#else
    if ((node.getObjectIDs()[i].distance <= pq + q.radius) &&
        (node.getObjectIDs()[i].distance >= pq - q.radius)) {
#endif
      Distance d = 0;
      try {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	d = objectSpace->getComparator()(q.object, *q.vptree->getObjectRepository().get(node.getObjectIDs(leafNodes.allocator)[i].id));
#else
	d = objectSpace->getComparator()(q.object, *q.vptree->getObjectRepository().get(node.getObjectIDs()[i].id));
#endif
      } catch(...) {
        NGTThrowException("VpTree::LeafNode::search: Internal fatal error : Cannot get object");
      }
      if (d <= q.radius) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
        r.id = node.getObjectIDs(leafNodes.allocator)[i].id;
#else
        r.id = node.getObjectIDs()[i].id;
#endif
        r.distance = d;
	so.getResult().push_back(r);
	std::sort(so.getResult().begin(), so.getResult().end());
	if (so.getResult().size() > q.size) {
	  so.getResult().resize(q.size);
	}
      }
    }
  }
}

void 
DVPTree::search(SearchContainer &sc) {
  ((SearchContainer&)sc).vptree = this;
  Node *root = getRootNode();
  assert(root != 0);
  if (sc.mode == DVPTree::SearchContainer::SearchLeaf) {
    if (root->id.getType() == Node::ID::Leaf) {
      sc.nodeID.setRaw(root->id.get());
      return;
    }
  }

  UncheckedNode uncheckedNode;
  uncheckedNode.push(root->id);

  while (!uncheckedNode.empty()) {
    Node::ID nodeid = uncheckedNode.top();
    uncheckedNode.pop();
    Node *cnode = getNode(nodeid);
    if (cnode == 0) {
      cerr << "Error! child node is null. but continue." << endl;
      continue;
    }
    if (cnode->id.getType() == Node::ID::Internal) {
      search(sc, (InternalNode&)*cnode, uncheckedNode);
    } else if (cnode->id.getType() == Node::ID::Leaf) {
      search(sc, (LeafNode&)*cnode, uncheckedNode);
    } else {
      abort();
    }
  }
}

