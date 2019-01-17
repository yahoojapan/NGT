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
#pragma once

#include	"NGT/Index.h"

namespace NGT {


class Command {
public:
  class SearchParameter {
  public:
    SearchParameter() {}
    SearchParameter(Args &args) {
      openMode	= args.getChar("m", 'r');
      try {
	query = args.get("#2");
      } catch (...) {
	NGTThrowException("ngt: Error: Query is not specified");
      }
      querySize = args.getl("Q", 0);
      indexType	= args.getChar("i", 't');
      size	= args.getl("n", 20);
      // edgeSize
      // -1 : using the size which was specified at the index creation.
      //  0 : no limitation for the edge size.
      // e(-2) : automatically set it according to epsilon.
      if (args.getChar("E", '-') == 'e') {
	edgeSize	= -2;
      } else {
	edgeSize	= args.getl("E", -1);
      }
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
    char	openMode;
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

  Command():debugLevel(0) {}

  void create(Args &args);
  void append(Args &args);
  static void search(NGT::Index &index, SearchParameter &searchParameter, ostream &stream)
  {
    ifstream		is(searchParameter.query);
    if (!is) {
      cerr << "Cannot open the specified file. " << searchParameter.query << endl;
      return;
    }
    search(index, searchParameter, is, stream);
  }
  static void search(NGT::Index &index, SearchParameter &searchParameter, istream &is, ostream &stream);
  void search(Args &args);
  void remove(Args &args);
  void exportIndex(Args &args);
  void importIndex(Args &args);
  void prune(Args &args);
  void reconstructGraph(Args &args);


  void info(Args &args);
  void setDebugLevel(int level) { debugLevel = level; }
  int getDebugLevel() { return debugLevel; }

protected:
  int debugLevel;

};

}; // NGT
