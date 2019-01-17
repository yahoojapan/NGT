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

#include	"NGT/Index.h"

int
main(int argc, char **argv)
{
  string	indexFile	= "index";	// Index.
  string	query		= "./data/sift-query-3.tsv";	// Query file.
  int		size		= 20;		// The number of resultant objects.
  float		radius		= FLT_MAX;	// Radius of search range.
  float		epsilon		= 0.1;		// Epsilon to expand explored range.

  try {
    NGT::Index	index(indexFile);		// open the specified index.
    ifstream	is(query);			// open a query file.
    if (!is) {
      cerr << "Cannot open the specified file. " << query << endl;
      return 1;
    }
    string line;
    if (getline(is, line)) {    		// read  a query object from a query file.
      NGT::Object *query = 0;
      {
	vector<string> tokens;
	NGT::Common::tokenize(line, tokens, " \t");      // split a string into words by the separators.
	// create a vector from the words.
	vector<double> obj;
	for (vector<string>::iterator ti = tokens.begin(); ti != tokens.end(); ++ti) {
	  obj.push_back(NGT::Common::strtol(*ti));
	}
	// allocate query object.
	query = index.allocateObject(obj);
      }
      // set search prameters.
      NGT::SearchContainer sc(*query);		// search parametera container.
      NGT::ObjectDistances objects;		// a result set.
      sc.setResults(&objects);			// set the result set.
      sc.setSize(size);				// the number of resultant objects.
      sc.setRadius(radius);			// search radius.
      sc.setEpsilon(epsilon);			// set exploration coefficient.

      index.search(sc);

      // output resultant objects.
      cout << "Rank\tID\tDistance" << endl;
      for (size_t i = 0; i < objects.size(); i++) {
	cout << i + 1 << "\t" << objects[i].id << "\t" << objects[i].distance << endl;
      }
      // delete the query object.
      index.deleteObject(query);
    }
  } catch (NGT::Exception &err) {
    cerr << "Error " << err.what() << endl;
    return 1;
  } catch (...) {
    cerr << "Error" << endl;
    return 1;
  }

  return 0;
}


