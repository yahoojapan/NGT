
#include	"NGT/NGTQ/QuantizedGraph.h"

int
main(int argc, char **argv)
{
#ifdef NGTQ_QBG
  string	indexPath	= "index";
  string	objectFile	= "./data/sift-dataset-5k.tsv";
  string	queryFile	= "./data/sift-query-3.tsv";

  // index construction
  try {
    NGT::Property	property;
    property.dimension		= 128;
    property.objectType		= NGT::ObjectSpace::ObjectType::Uint8;
    property.distanceType	= NGT::Index::Property::DistanceType::DistanceTypeL2;
    std::cout << "creating the index framework..." << std::endl;
    NGT::Index::create(indexPath, property);
    NGT::Index	index(indexPath);
    ifstream	is(objectFile);
    string	line;
    std::cout << "appending the objects..." << std::endl;
    while (getline(is, line)) {
      vector<float>	obj;
      stringstream	linestream(line);
      while (!linestream.eof()) {
	int value;
	linestream >> value;
	if (linestream.fail()) {
	  obj.clear();
	  break;
	}
	obj.push_back(value);
      }
      if (obj.empty()) {
	cerr << "An empty line or invalid value: " << line << endl;
	continue;
      }
      obj.resize(property.dimension);  // cut off additional data in the file.
      index.insert(obj);
    }
    std::cout << "building the index..." << std::endl;
    index.createIndex(16);
    index.save();
  } catch (NGT::Exception &err) {
    cerr << "Error " << err.what() << endl;
    return 1;
  } catch (...) {
    cerr << "Error" << endl;
    return 1;
  }

  // quantization
  size_t dimensionOfSubvector = 1;
  size_t maxNumberOfEdges = 50;
  try {
    std::cout << "quantizing the index..." << std::endl;
    NGTQG::Index::quantize(indexPath, dimensionOfSubvector, maxNumberOfEdges, true);
  } catch (NGT::Exception &err) {
    cerr << "Error " << err.what() << endl;
    return 1;
  } catch (...) {
    cerr << "Error" << endl;
    return 1;
  }

  // nearest neighbor search
  try {
    NGT::Index		index(indexPath);
    NGT::Property	property;
    index.getProperty(property);
    ifstream		is(queryFile);
    string		line;
    std::cout << "searching the index..." << std::endl;
    while (getline(is, line)) {
      vector<uint8_t>	query;
      {
	stringstream	linestream(line);
	while (!linestream.eof()) {
	  int value;
	  linestream >> value;
	  query.push_back(value);
	}
	query.resize(property.dimension);
	cout << "Query : ";
	for (size_t i = 0; i < 5; i++) {
	  cout << static_cast<int>(query[i]) << " ";
	}
	cout << "...";
      }

      NGT::SearchQuery		sc(query);
      NGT::ObjectDistances	objects;
      sc.setResults(&objects);
      sc.setSize(10);
      sc.setEpsilon(0.1);

      index.search(sc);
      cout << endl << "Rank\tID\tDistance: Object" << std::showbase << endl;
      for (size_t i = 0; i < objects.size(); i++) {
	cout << i + 1 << "\t" << objects[i].id << "\t" << objects[i].distance << "\t: ";
	NGT::ObjectSpace &objectSpace = index.getObjectSpace();
	uint8_t *object = static_cast<uint8_t*>(objectSpace.getObject(objects[i].id));
	for (size_t idx = 0; idx < 5; idx++) {
	  cout << static_cast<int>(object[idx]) << " ";
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
#endif
  return 0;
}


