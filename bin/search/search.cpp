
//
// Copyright (C) 2015-2016 Yahoo! JAPAN Research
//
// This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
// To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
//

#include	"NGT/Index.h"

int
main(int argc, char **argv)
{
  string	database	= "index";	// Index.
  string	query		= "./data/sift-query-3.tsv";	// Query file.
  int		size		= 20;		// The number of resultant objects.
  float		radius		= FLT_MAX;	// Radius of search range.
  float		epsilon		= 0.1;		// Epsilon to expand explored range.

  try {
    NGT::Index	index(database);		// open the specified index.
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
      cout << "Rank\tIN-ID\tID\tDistance" << endl;
      for (size_t i = 0; i < objects.size(); i++) {
	cout << i + 1 << "\t" << objects[i].id << "\t";
	cout << "(nul)\t";
	cout << objects[i].distance << endl;
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


