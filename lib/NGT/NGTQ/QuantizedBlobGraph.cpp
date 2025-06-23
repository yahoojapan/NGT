//
// Copyright (C) 2021 Yahoo Japan Corporation
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

#include "NGT/NGTQ/QuantizedBlobGraph.h"

void 
QBG::Index::append(const std::string &indexName, // index file
                   const std::string &data,      // data file
                   size_t dataSize,          // data size
                   bool verbose) {
  NGT::StdOstreamRedirector redirector(!verbose);
  redirector.begin();
  QBG::Index index(indexName);
  auto &quantizer = index.getQuantizer();
  istream *is;
  if (data == "-") {
    is = &cin;
  } else {
    ifstream *ifs = new ifstream;
    ifs->ifstream::open(data);
    if (!(*ifs)) {
      std::stringstream msg;
      msg << "Cannot open the specified file. " << data;
      NGTThrowException(msg);
    }
    is = ifs;
  }
  string line;
  size_t idx   = quantizer.objectList.size() == 0 ? 0 : quantizer.objectList.size() - 1;
  size_t count = 0;
  // extract objects from the file and insert them to the object list.
  while (getline(*is, line)) {
    idx++;
    count++;
    std::vector<float> object;
    NGT::Common::extractVector(line, " ,\t", object);
    if (object.empty()) {
      cerr << "Empty line or invalid value: " << line << endl;
      continue;
    }
    if ((quantizer.property.distanceType == NGTQ::DistanceType::DistanceTypeInnerProduct) &&
        (object.size() + 1 == quantizer.objectList.genuineDimension)) {
      object.emplace_back(0);
    }
    if (object.size() != quantizer.objectList.genuineDimension) {
      std::stringstream msg;
      msg << "The dimension of the specified object is inconsistent with the dimension of the index. "
          << object.size() << ":" << quantizer.objectList.genuineDimension;
      NGTThrowException(msg);
    }
    index.insert(idx, object);

    if (count % 100000 == 0) {
      std::cerr << "appended " << static_cast<float>(count) / 1000000.0 << "M objects.";
      if (count != idx) {
        std::cerr << " # of the total objects=" << static_cast<float>(idx) / 1000000.0 << "M";
      }
      cerr << " virtual memory(kbyte)=" << NGT::Common::getProcessVmSize() << std::endl;
    }
  }
  if (data != "-") {
    delete is;
  }

  index.save();
  index.close();
  redirector.end();
}

void 
QBG::Index::append(const std::string &indexName, // index file
                   NGT::ObjectSpace &objectSpace, // object space including objects
                   bool verbose) {
  NGT::StdOstreamRedirector redirector(!verbose);
  redirector.begin();
  QBG::Index index(indexName);
  auto &quantizer = index.getQuantizer();
  auto &repo      = objectSpace.getRepository();

  size_t count = 0;
  // extract objects from the file and insert them to the object list.
  for (size_t idx = 1; idx < repo.size(); idx++) {
    std::vector<float> object;
    object.clear();
    if (repo.isEmpty(idx)) {
      continue;
    }
    objectSpace.getObject(idx, object);
    count++;
    if (object.size() != quantizer.objectList.genuineDimension) {
      std::stringstream msg;
      msg << "The dimension of the specified object is inconsistent with the dimension of the index. "
          << object.size() << ":" << quantizer.objectList.genuineDimension;
      NGTThrowException(msg);
    }
    index.insert(idx, object);

    if (count % 100000 == 0) {
      std::cerr << "appended " << static_cast<float>(count) / 1000000.0 << "M objects.";
      if (count != idx) {
        std::cerr << " # of the total objects=" << static_cast<float>(idx) / 1000000.0 << "M";
      }
      cerr << " virtual memory(kbyte)=" << NGT::Common::getProcessVmSize() << std::endl;
    }
  }

  redirector.end();
}

void
QBG::Index::appendBinary(const std::string &indexName, // index file
                         const std::string &data,      // data file
                         size_t dataSize,          // data size
                         bool verbose) {
  NGT::StdOstreamRedirector redirector(!verbose);
  redirector.begin();
  QBG::Index index(indexName);
  std::vector<std::string> tokens;
  NGT::Common::tokenize(data, tokens, ".");
  if (tokens.size() < 2) {
    std::stringstream msg;
    msg << "Invalid file name format. " << data;
    NGTThrowException(msg);
  }
  auto &quantizer = index.getQuantizer();
  StaticObjectFileLoader loader(data, tokens[tokens.size() - 1]);
  size_t idx   = quantizer.objectList.size() == 0 ? 0 : quantizer.objectList.size() - 1;
  size_t count = 0;
  while (!loader.isEmpty()) {
    idx++;
    count++;
    if (dataSize > 0 && idx > dataSize) {
      break;
    }
    auto object = loader.getObject();
    if ((quantizer.property.distanceType == NGTQ::DistanceType::DistanceTypeInnerProduct) &&
        (object.size() + 1 == quantizer.objectList.genuineDimension)) {
      object.emplace_back(0);
    }
    if (object.size() != quantizer.objectList.genuineDimension) {
      std::stringstream msg;
      msg << "The dimension of the specified object is inconsistent with the dimension of the index. "
          << object.size() << ":" << quantizer.objectList.genuineDimension;
      NGTThrowException(msg);
    }
    index.insert(idx, object);
    if (count % 1000000 == 0) {
      std::cerr << "appended " << static_cast<float>(count) / 1000000.0 << "M objects.";
      if (count != idx) {
        std::cerr << " # of the total objects=" << static_cast<float>(idx) / 1000000.0 << "M";
      }
      cerr << " virtual memory(kbyte)=" << NGT::Common::getProcessVmSize() << std::endl;
    }
  }
  index.save();
  index.close();
  redirector.end();
}


void
QBG::Index::preprocessingForNGT(std::string &indexPath, std::string &objectPath, bool verbose) {
  NGT::Property prop;
  {
    if (verbose) {
      std::cerr << "Opening the NGT index to get the property..." << std::endl;
    }
    NGT::Index index(indexPath);
    index.getProperty(prop);
  }
  if (verbose) {
    std::cerr << "object type:" << prop.objectType << std::endl;
    std::cerr << "dimension:" << prop.dimension << std::endl;
    std::cerr << "distance:" << prop.distanceType << std::endl;
  }
  if (prop.objectType != NGT::ObjectSpace::Qint4) {
    std::stringstream msg;
    msg << "Not qint4 object type. " << prop.objectType << std::endl;
    NGTThrowException(msg);
  }
  if (prop.distanceType != NGT::ObjectSpace::DistanceTypeL2) {
    std::stringstream msg;
    msg << "Not L2 distance type. " << prop.distanceType << std::endl;
    NGTThrowException(msg);
  }

  QBG::BuildParameters buildParameters;
  buildParameters.creation.localClusterDataType = NGTQ::ClusterDataTypePQ4;
  buildParameters.creation.genuineDimension = prop.dimension;
  buildParameters.creation.dimension = ((buildParameters.creation.genuineDimension + 15) / 16) * 16;
  buildParameters.creation.numOfSubvectors = prop.dimension;
  buildParameters.creation.distanceType = NGTQ::DistanceType::DistanceTypeL2;
  buildParameters.creation.genuineDataType = ObjectFile::DataTypeFloat;
  buildParameters.creation.globalObjectType = NGT::ObjectSpace::ObjectType::Float;

  std::string qbgIndexPath = indexPath + "/" + NGT::Quantizer::getQbgIndex();
  if (verbose) {
    std::cerr << "qbg: creating.." << std::endl;
  }
  QBG::Index::create(qbgIndexPath, buildParameters, 0, objectPath);

  if (verbose) {
    std::cerr << "qbg: appending..." << std::endl;
  }
  size_t dataSize = 0;
  std::string mode = "";
  QBG::Index::append(qbgIndexPath, objectPath, dataSize, verbose);


  QBG::Optimizer optimizer;

  optimizer.unifiedPQ     = true;
  optimizer.rotation      = false;
  optimizer.repositioning = false;
  optimizer.globalType = QBG::Optimizer::GlobalTypeZero;

  if (verbose) {
    std::cerr << "qbg: optimizing..." << std::endl;
  }
  optimizer.optimize(qbgIndexPath);

  if (verbose) {
    std::cerr << "qbg: building..." << std::endl;
  }
  QBG::Index::buildNGTQ(qbgIndexPath, verbose);
}
