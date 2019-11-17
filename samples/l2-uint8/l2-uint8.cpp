
#include	"NGT/Index.h"

using namespace std;

int
main(int argc, char **argv)
{
  string	indexFile	= "index";
  string	objectFile	= "./data/sift-dataset-5k.tsv";
  string	queryFile	= "./data/sift-query-3.tsv";
  // index construction
  try {
    NGT::Property	property;
    property.dimension		= 128;
    property.objectType		= NGT::ObjectSpace::ObjectType::Float;
    property.distanceType	= NGT::Index::Property::DistanceType::DistanceTypeL2;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    NGT::Index	index(property, indexFile);
#else
    NGT::Index	index(property);
#endif
    ifstream	is(objectFile);
    string	line;
    while (getline(is, line)) {
      vector<string>	tokens;
      NGT::Common::tokenize(line, tokens, " \t");      // split an object string into values by separators.
      vector<float>	obj;
      for (vector<string>::iterator ti = tokens.begin(); ti != tokens.end(); ++ti) {
	obj.push_back(NGT::Common::strtol(*ti));
      }
      obj.resize(property.dimension);  // cut off additional data in the file.
      index.append(obj);
    }
    index.createIndex(16);
    index.saveIndex(indexFile);
  } catch (NGT::Exception &err) {
    cerr << "Error " << err.what() << endl;
    return 1;
  } catch (...) {
    cerr << "Error" << endl;
    return 1;
  }

  // nearest neighbor search
  try {
    NGT::Index		index(indexFile);
    NGT::Property	property;
    index.getProperty(property);
    ifstream		is(queryFile);
    string		line;
    while (getline(is, line)) {
      vector<float>	query;
      {
	vector<string>	tokens;
	NGT::Common::tokenize(line, tokens, " \t");
	for (vector<string>::iterator ti = tokens.begin(); ti != tokens.end(); ++ti) {
	  query.push_back(NGT::Common::strtol(*ti));
	}
	query.resize(property.dimension);
	cout << "Query : ";
	for (size_t i = 0; i < 5; i++) {
	  cout << query[i] << " ";
	}
	cout << "...";
      }

      NGT::SearchQuery		sc(query);
      NGT::ObjectDistances	objects;
      sc.setResults(&objects);
      sc.setSize(10);
      sc.setEpsilon(0.1);

      index.search(sc);
      cout << endl << "Rank\tID\tDistance" << std::showbase << endl;
      for (size_t i = 0; i < objects.size(); i++) {
	cout << i + 1 << "\t" << objects[i].id << "\t" << objects[i].distance << "\t: ";
	NGT::ObjectSpace &objectSpace = index.getObjectSpace();
	float *object = static_cast<float*>(objectSpace.getObject(objects[i].id));
	for (size_t idx = 0; idx < 5; idx++) {
	  cout << object[idx] << " ";
	}
	cout << "..." << endl;
      }
      cout << endl;
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


