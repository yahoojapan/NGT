
#include	"NGT/Index.h"
#include	"NGT/NGTQ/Capi.h"
int
main(int argc, char **argv)
{
#if !defined(NGT_SHARED_MEMORY_ALLOCATOR)
  std::string indexPath  = "qg-index";
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

  NGTError err = ngt_create_error_object();
  NGTProperty prop = ngt_create_property(err);
  if (prop == NULL) {
    std::cerr << ngt_get_error_string(err) << std::endl;
    return 1;
  }
  size_t dimension = 128;
  ngt_set_property_dimension(prop, dimension, err);

  std::cerr << "create an empty index..." << std::endl;
  NGTIndex index = ngt_create_graph_and_tree(indexPath.c_str(), prop, err);
  if (index == NULL) {
    std::cerr << ngt_get_error_string(err) << std::endl;
    return 1;
  }

  std::cerr << "insert objects..." << std::endl;
  try {
    std::ifstream is(objectFile);
    std::string	line;
    while (getline(is, line)) {
      std::vector<double> obj;
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
      if (ngt_insert_index(index, obj.data(), dimension, err) == 0) {
	std::cerr << ngt_get_error_string(err) << std::endl;
	return 1;
      }
    }
  } catch (NGT::Exception &err) {
    std::cerr << "Error " << err.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Error" << std::endl;
    return 1;
  }

  std::cerr << "build the index..." << std::endl;
  if (ngt_create_index(index, 100, err) == false) {
    std::cerr << "Error:" << ngt_get_error_string(err) << std::endl;
    return 1;
  }

  std::cerr << "save the index..." << std::endl;
  if (ngt_save_index(index, indexPath.c_str(), err) == false) {
    std::cerr << ngt_get_error_string(err) << std::endl;
    return 1;
  }

  std::cerr << "close the index..." << std::endl;
  ngt_close_index(index);

  NGTQGQuantizationParameters quantizationParameters;
  ngtqg_initialize_quantization_parameters(&quantizationParameters);

  std::cerr << "quantize the index..." << std::endl;
  if (ngtqg_quantize(indexPath.c_str(), quantizationParameters, err) == false) {
    std::cerr << ngt_get_error_string(err) << std::endl;
    return 1;
  }

  std::cerr << "open the quantized index..." << std::endl;
  index = ngtqg_open_index(indexPath.c_str(), err);
  if (index == NULL) {
    std::cerr << ngt_get_error_string(err) << std::endl;
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
    NGTObjectDistances result = ngt_create_empty_results(err);
    NGTQGQuery query;
    ngtqg_initialize_query(&query);
    query.query = queryVector;
    query.size = 10;
    query.result_expansion = 100;
    query.epsilon = 0.1;
    std::cerr << "search the index for the specified query..." << std::endl;
    ngtqg_search_index(index, query, result, err);

    auto rsize = ngt_get_result_size(result, err);
    std::cout << "Rank\tID\tDistance" << std::endl;
    for (size_t i = 0; i < rsize; i++) {
      NGTObjectDistance object = ngt_get_result(result, i, err);
      std::cout << i + 1 << "\t" << object.id << "\t" << object.distance << std::endl;
    }

    ngt_destroy_results(result);
  }

  std::cerr << "close the quantized index" << std::endl;
  ngtqg_close_index(index);
  ngt_destroy_error_object(err);
#endif
  return 0;
}
