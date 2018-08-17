//
// Copyright (C) 2015-2018 Yahoo Japan Corporation
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

class Args : public map<string, string>
{
public:
  Args() {}
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
