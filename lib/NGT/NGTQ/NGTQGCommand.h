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
#pragma once

#include	"NGT/NGTQ/QuantizedGraph.h"
#include	"NGT/Command.h"
#include	"NGT/NGTQ/NGTQCommand.h"

namespace NGTQG {
  
  class Command : public NGT::Command {
  public:
    class CreateParameters : public NGTQ::Command::CreateParameters {
    public:
      CreateParameters(NGT::Args &args): NGTQ::Command::CreateParameters(args, 's', 1, 16) {
        setup(args, property.dimension);
      }
      CreateParameters(NGT::Args &args, int dimension): NGTQ::Command::CreateParameters(args, 's', 1, 16) {
        setup(args, dimension);
      }
      void setup(NGT::Args &args, size_t dimension) {
        size_t quantizationRatio = args.getl("Q", 0);
        if (quantizationRatio == 0) {
	  if ((dimension > 400) && (dimension % 2 == 0)) {
	    property.localDivisionNo = dimension / 2;
	  } else {
	    property.localDivisionNo = dimension;
	  }
	} else {
	  property.localDivisionNo = dimension;
	}
	property.dimension = dimension;
	char localCentroidCreationMode = args.getChar("l", 'k');
	switch(localCentroidCreationMode) {
	case 'd': property.localCentroidCreationMode = NGTQ::CentroidCreationModeDynamic; break;
	case 's': property.localCentroidCreationMode = NGTQ::CentroidCreationModeStatic; break;
	case 'k': property.localCentroidCreationMode = NGTQ::CentroidCreationModeDynamicKmeans; break;
	default:
	  std::stringstream msg;
	  msg << "NGTQG::Command::CreateParameters: Error: Invalid centroid creation mode. " << localCentroidCreationMode;
	  NGTThrowException(msg);
	}
      }
    };

    class SearchParameters : public NGT::Command::SearchParameters {
    public:
      SearchParameters(NGT::Args &args): NGT::Command::SearchParameters(args, "0.02") {
        stepOfResultExpansion = 2;
	std::string resultExpansion = args.getString("p", "3.0");
	std::vector<std::string> tokens;
	NGT::Common::tokenize(resultExpansion, tokens, ":");
	if (tokens.size() >= 1) { beginOfResultExpansion = endOfResultExpansion = NGT::Common::strtod(tokens[0]); }
	if (tokens.size() >= 2) { endOfResultExpansion = NGT::Common::strtod(tokens[1]); }
	if (tokens.size() >= 3) { stepOfResultExpansion = NGT::Common::strtod(tokens[2]); }
      }
      float	beginOfResultExpansion;
      float	endOfResultExpansion;
      float	stepOfResultExpansion;
  };

#if defined(NGT_SHARED_MEMORY_ALLOCATOR) || defined(NGTQ_SHARED_INVERTED_INDEX)
    void create(NGT::Args &args) {};
    void build(NGT::Args &args) {};
    void quantize(NGT::Args &args) {};
    void search(NGT::Args &args) {};
    void info(NGT::Args &args) {};
#else
    void create(NGT::Args &args);
    void build(NGT::Args &args);
    void quantize(NGT::Args &args);
    void search(NGT::Args &args);
    void info(NGT::Args &args);
#endif
    
    void setDebugLevel(int level) { debugLevel = level; }
    int getDebugLevel() { return debugLevel; }

    void help() {
      cerr << "Usage : ngtqg command database [data]" << endl;
      cerr << "           command : create build quantize search" << endl;
    }

    void execute(NGT::Args args) {
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
	  cerr << "ngt::command=" << command << endl;
	}
	if (command == "search") {
	  search(args);
	} else if (command == "create") {
	  create(args);
	} else if (command == "build") {
	  build(args);
	} else if (command == "quantize") {
	  quantize(args);
	} else if (command == "info") {
	  info(args);
	} else {
	  cerr << "Illegal command. " << command << endl;
	}
      } catch(NGT::Exception &err) {
	cerr << "ngtqg: Fatal error: " << err.what() << endl;
      }
    }

  };

}; // NGTQG
