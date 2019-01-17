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
#include	"NGT/Optimizer.h"

#ifndef BUILD_DATE
#define BUILD_DATE	"-"
#endif
#ifndef GIT_HASH
#define GIT_HASH	"-"
#endif
#ifndef GIT_DATE
#define GIT_DATE	"-"
#endif
#ifndef GIT_TAG
#define GIT_TAG	"-"
#endif

static void
version(ostream &os)
{
  os << "ngt:" << endl;
  os << "  Built date:" << BUILD_DATE << endl;
  os << "  The last git tag:" << GIT_TAG << endl;
  os << "  The last git commit hash:" << GIT_HASH << endl;
  os << "  The last git commit date:" << GIT_DATE << endl;
}

void help() {
  cerr << "Usage : ngt command database data" << endl;
  cerr << "           command : create search remove append export import prune reconstruct-graph" << endl;
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

  ngt.setDebugLevel(args.getl("X", 0));

  try {
    if (ngt.getDebugLevel() >= 1) {
      cerr << "ngt: command=" << command << endl;
    }
    if (command == "search") {
      ngt.search(args);
    } else if (command == "create") {
      ngt.create(args);
    } else if (command == "append") {
      ngt.append(args);
    } else if (command == "remove") {
      ngt.remove(args);
    } else if (command == "export") {
      ngt.exportIndex(args);
    } else if (command == "import") {
      ngt.importIndex(args);
    } else if (command == "prune") {
      ngt.prune(args);
    } else if (command == "reconstruct-graph") {
      ngt.reconstructGraph(args);
#ifndef NGT_SHARED_MEMORY_ALLOCATOR
    } else if (command == "extract-query") {
      NGT::Optimizer::extractQueries(args);
    } else if (command == "adjust-edge-size") {
      NGT::Optimizer::adjustBaseSearchEdgeSize(args);
#endif

    } else if (command == "info") {
      version(cerr);
      NGT::Index::version(cerr);
      ngt.info(args);
    } else {
      cerr << "ngt: Error: Illegal command. " << command << endl;
    }
  } catch(NGT::Exception &err) {
    cerr << "ngt: Error: " << err.what() << endl;
  }

}
