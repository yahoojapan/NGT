
#include	"NGT/Index.h"

using namespace std;

int
main(int argc, char **argv)
{
  string	indexFile	= "index";
  string	objectFile	= "./data/sift-dataset-5k.tsv";
  string	queryFile	= "./data/sift-query-3.tsv";
  size_t	bitSize		= 64;
  // index construction
  try {
    NGT::Property	property;
    property.dimension		= bitSize / 8;
    property.objectType		= NGT::ObjectSpace::ObjectType::Uint8;
    property.distanceType	= NGT::Index::Property::DistanceType::DistanceTypeHamming;
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
      vector<uint8_t>	obj;
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
    string	line;
    while (getline(is, line)) {
      vector<uint8_t>	query;
      {
	vector<string>	tokens;
	NGT::Common::tokenize(line, tokens, " \t");
	for (vector<string>::iterator ti = tokens.begin(); ti != tokens.end(); ++ti) {
	  query.push_back(NGT::Common::strtol(*ti));
	}
	query.resize(property.dimension);
	cout << "Query : ";
	for (size_t i = 0; i < query.size(); i++) {
	  cout << std::bitset<8>(query[i]) << " ";
	}
      }

      NGT::SearchQuery		sc(query);
      NGT::ObjectDistances	objects;
      sc.setResults(&objects);
      sc.setSize(10);
      sc.setEpsilon(0.2);

      index.search(sc);
      cout << endl << "Rank\tID\tDistance" << std::showbase << endl;
      for (size_t i = 0; i < objects.size(); i++) {
	cout << i + 1 << "\t" << objects[i].id << "\t" << objects[i].distance << "\t: ";
	NGT::ObjectSpace &objectSpace = index.getObjectSpace();
	uint8_t *object = static_cast<uint8_t*>(objectSpace.getObject(objects[i].id));
	for (size_t idx = 0; idx < objectSpace.getDimension(); idx++) {
	  cout << std::bitset<8>(object[idx]) << " ";
	}
	cout << endl;
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


