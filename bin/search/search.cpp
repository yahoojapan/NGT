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
  std::string	indexFile	= "index";	// Index.
  std::string	query		= "./data/sift-query-3.tsv";	// Query file.
  int		size		= 20;		// The number of resultant objects.
  float		radius		= FLT_MAX;	// Radius of search range.
  float		epsilon		= 0.1;		// Exploration coefficient(epsilon) to expand explored range

  try {
    NGT::Index	index(indexFile);		// open the specified index.
    std::ifstream	is(query);			// open a query file.
    if (!is) {
      std::cerr << "Cannot open the specified file. " << query << std::endl;
      return 1;
    }
    std::string line;
    if (getline(is, line)) {    		// read  a query object from a query file.
      std::vector<double> queryObject;
      {
	std::vector<std::string> tokens;
	NGT::Common::tokenize(line, tokens, " \t");      // split a string into words by the separators.
	// create a vector from the words.
	for (std::vector<std::string>::iterator ti = tokens.begin(); ti != tokens.end(); ++ti) {
	  queryObject.push_back(NGT::Common::strtol(*ti));
	}
      }
      // set search prameters.
      NGT::SearchQuery query(queryObject);	// search query
      NGT::ObjectDistances objects;		// result objects
      query.setResults(&objects);		// set the result objects.
      query.setSize(size);			// the number of result objects
      query.setRadius(radius);			// search radius
      query.setEpsilon(epsilon);		// set exploration coefficient.

      index.search(query);

      // output resultant objects.
      std::cout << "Rank\tID\tDistance" << std::endl;
      for (size_t i = 0; i < objects.size(); i++) {
	std::cout << i + 1 << "\t" << objects[i].id << "\t" << objects[i].distance << std::endl;
      }
    }
  } catch (NGT::Exception &err) {
    std::cerr << "Error " << err.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Error" << std::endl;
    return 1;
  }

  return 0;
}


