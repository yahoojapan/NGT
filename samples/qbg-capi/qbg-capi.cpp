
#include	"NGT/Index.h"
#include	"NGT/NGTQ/Capi.h"
int
main(int argc, char **argv)
{
#if !defined(NGT_SHARED_MEMORY_ALLOCATOR)
  std::string indexPath  = "qbg-index";
  std::string objectFile = "sift-128-euclidean.tsv";
  std::string queryFile  = "query.tsv";

  std::cerr << "run the following commands to prepare data for this sample program." << std::endl;
  std::cerr << "  curl -L -O https://github.com/yahoojapan/NGT/raw/main/tests/datasets/ann-benchmarks/sift-128-euclidean.tsv" << std::endl;
  std::cerr << "  curl -L -O https://github.com/yahoojapan/NGT/raw/main/tests/datasets/ann-benchmarks/sift-128-euclidean_query.tsv" << std::endl;
  std::cerr << "  head -1 sift-128-euclidean_query.tsv > query.tsv" << std::endl;
  std::cerr << std::endl;
  std::cerr << "index path=" << indexPath << std::endl;
  std::cerr << "object file=" << objectFile << std::endl;
  std::cerr << "query file=" << queryFile << std::endl;
  std::cerr << std::endl;

  {
    std::cerr << "remove the existing index. " << indexPath << std::endl;
    const std::string com = "rm -rf " + indexPath;
    if (system(com.c_str()) == -1) {
      std::cerr << "Cannot exec. " << com << std::endl;
    }
  }

  size_t dimension = 128;
  NGTError err = ngt_create_error_object();

  std::cerr << "create an empty index..." << std::endl;
  QBGConstructionParameters constructionParameters;
  qbg_initialize_construction_parameters(&constructionParameters);
  constructionParameters.dimension = dimension;
  constructionParameters.number_of_subvectors = 64;
  constructionParameters.number_of_blobs = 0;
  if (!qbg_create(indexPath.c_str(), &constructionParameters, err)) {
    std::cerr << "Cannot create" << std::endl;
    std::cerr << ngt_get_error_string(err) << std::endl;
    return 1;
  }

  std::cerr << "append objects..." << std::endl;
  auto index = qbg_open_index(indexPath.c_str(), false, err);
  if (index == 0) {
    std::cerr << "Cannot open" << std::endl;
    std::cerr << ngt_get_error_string(err) << std::endl;
    return 1;
  }

  try {
    std::ifstream is(objectFile);
    std::string	line;
    while (getline(is, line)) {
      std::vector<float> obj;
      std::stringstream	linestream(line);
      while (!linestream.eof()) {
	float value;
	linestream >> value;
	if (linestream.fail()) {
	  obj.clear();
	  break;
	}
	obj.push_back(value);
      }
      if (obj.empty()) {
	std::cerr << "An empty line or invalid value: " << line << std::endl;
	return 1;
      }
      if (qbg_append_object(index, obj.data(), dimension, err) == 0) {
	std::cerr << ngt_get_error_string(err) << std::endl;
	return 1;
      }
    }
  } catch (...) {
    std::cerr << "Error" << std::endl;
    return 1;
  }

  qbg_save_index(index, err);
  qbg_close_index(index);

  std::cerr << "building the index..." << std::endl;
  QBGBuildParameters buildParameters;
  qbg_initialize_build_parameters(&buildParameters);
  buildParameters.number_of_objects = 500;		
  auto status = qbg_build_index(indexPath.c_str(), &buildParameters, err);
  if (!status) {
    std::cerr << "Cannot build. " << ngt_get_error_string(err) << std::endl;
    return 1;
  }

  index = qbg_open_index(indexPath.c_str(), true, err);
  if (index == 0) {
    std::cerr << "Cannot open. " << ngt_get_error_string(err) << std::endl;
    return 1;
  }

  std::ifstream	is(queryFile);
  if (!is) {
    std::cerr << "Cannot open the specified file. " << queryFile << std::endl;
    return 1;
  }

  std::string line;
  float queryVector[dimension];
  if (getline(is, line)) {
    std::vector<double> queryObject;
    {
      std::vector<std::string> tokens;
      NGT::Common::tokenize(line, tokens, " \t");
      if (tokens.size() != dimension) {
	std::cerr << "dimension of the query is invalid. dimesion=" << tokens.size() << ":" << dimension << std::endl;
	return 1;
      }
      for (std::vector<std::string>::iterator ti = tokens.begin(); ti != tokens.end(); ++ti) {
	queryVector[distance(tokens.begin(), ti)] = NGT::Common::strtod(*ti);
      }
    }
    QBGObjectDistances result = ngt_create_empty_results(err);
    QBGQuery query;
    qbg_initialize_query(&query);
    query.query = &queryVector[0];
    std::cerr << "search the index for the specified query..." << std::endl;
    auto status = qbg_search_index(index, query, result, err);
    if (!status) {
      std::cerr << "Cannot search. " << ngt_get_error_string(err) << std::endl;
      return 1;
    }
    auto rsize = qbg_get_result_size(result, err);
    std::cout << "Rank\tID\tDistance" << std::endl;
    for (size_t i = 0; i < rsize; i++) {
      NGTObjectDistance object = qbg_get_result(result, i, err);
      std::cout << i + 1 << "\t" << object.id << "\t" << object.distance << std::endl;
    }

    qbg_destroy_results(result);
  }

  qbg_close_index(index);
  ngt_destroy_error_object(err);
#endif
  return 0;
}
