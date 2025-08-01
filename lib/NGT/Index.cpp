//
// Copyright (C) 2015 Yahoo Japan Corporation
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

#include <algorithm>

#include "NGT/defines.h"
#include "NGT/Common.h"
#include "NGT/ObjectSpaceRepository.h"
#include "NGT/Index.h"
#include "NGT/Thread.h"
#include "NGT/GraphReconstructor.h"
#include "NGT/Version.h"
#include "NGT/NGTQ/ObjectFile.h"

using namespace std;
using namespace NGT;

void Index::version(ostream &os) {
  os << "libngt:" << endl;
  Version::get(os);
}

string Index::getVersion() { return Version::getVersion(); }

size_t Index::getDimension() { return static_cast<NGT::GraphIndex &>(getIndex()).getProperty().dimension; }

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
NGT::Index::Index(NGT::Property &prop, const string &database) : redirect(false) {
  if (prop.dimension == 0) {
    NGTThrowException("Index::Index. Dimension is not specified.");
  }
  Index *idx = 0;
  mkdir(database);
  if (prop.indexType == NGT::Index::Property::GraphAndTree) {
    idx = new NGT::GraphAndTreeIndex(database, prop);
  } else if (prop.indexType == NGT::Index::Property::Graph) {
    idx = new NGT::GraphIndex(database, prop);
  } else {
    NGTThrowException("Index::Index: Not found IndexType in property file.");
  }
  if (idx == 0) {
    stringstream msg;
    msg << "Index::Index: Cannot construct. ";
    NGTThrowException(msg);
  }
  index = idx;
  path  = "";
}
#else
NGT::Index::Index(NGT::Property &prop) : redirect(false) {
  if (prop.dimension == 0) {
    NGTThrowException("Index::Index. Dimension is not specified.");
  }
  Index *idx = 0;
  if (prop.indexType == NGT::Index::Property::GraphAndTree) {
    idx = new NGT::GraphAndTreeIndex(prop);
  } else if (prop.indexType == NGT::Index::Property::Graph) {
    idx = new NGT::GraphIndex(prop);
  } else {
    NGTThrowException("Index::Index: Not found IndexType in property file.");
  }
  if (idx == 0) {
    stringstream msg;
    msg << "Index::Index: Cannot construct. ";
    NGTThrowException(msg);
  }
  index = idx;
  path  = "";
}
#endif

float NGT::Index::getEpsilonFromExpectedAccuracy(double accuracy) {
  return static_cast<NGT::GraphIndex &>(getIndex()).getEpsilonFromExpectedAccuracy(accuracy);
}

void NGT::Index::open(const string &database, bool rdOnly, NGT::Index::OpenType openType) {
  NGT::Property prop;
  prop.load(database);
  Index *idx = 0;
  if ((prop.indexType == NGT::Index::Property::GraphAndTree) &&
      ((openType & NGT::Index::OpenTypeGraphDisabled) == 0)) {
    idx = new NGT::GraphAndTreeIndex(database, rdOnly);
  } else if ((prop.indexType == NGT::Index::Property::Graph) ||
             ((openType & NGT::Index::OpenTypeGraphDisabled) != 0)) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    idx = new NGT::GraphIndex(database, rdOnly);
#else
    idx = new NGT::GraphIndex(database, rdOnly, openType);
#endif
  } else {
    NGTThrowException("Index::Open: Not found IndexType in property file.");
  }
  if (idx == 0) {
    stringstream msg;
    msg << "Index::open: Cannot open. " << database;
    NGTThrowException(msg);
  }
  index = idx;
  path  = database;
}

void NGT::Index::createGraphAndTree(const string &database, NGT::Property &prop, const string &dataFile,
                                    size_t dataSize, bool redirect) {
  if (prop.dimension == 0) {
    NGTThrowException("Index::createGraphAndTree. Dimension is not specified.");
  }
  prop.indexType = NGT::Index::Property::IndexType::GraphAndTree;
  Index *idx     = 0;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  mkdir(database);
  idx = new NGT::GraphAndTreeIndex(database, prop);
#else
  idx = new NGT::GraphAndTreeIndex(prop);
#endif
  assert(idx != 0);
  StdOstreamRedirector redirector(redirect);
  redirector.begin();
  try {
    if (idx->getObjectSpace().isQintObjectType()) {
      idx->saveIndex(database);
      idx->close();
      auto append     = true;
      auto refinement = false;
      if (!dataFile.empty()) {
        appendFromTextObjectFile(database, dataFile, dataSize, append, refinement, prop.threadPoolSize);
      }
    } else {
      loadAndCreateIndex(*idx, database, dataFile, prop.threadPoolSize, dataSize);
    }
  } catch (Exception &err) {
    delete idx;
    redirector.end();
    throw err;
  }
  delete idx;
  redirector.end();
}

void NGT::Index::createGraph(const string &database, NGT::Property &prop, const string &dataFile,
                             size_t dataSize, bool redirect) {
  if (prop.dimension == 0) {
    NGTThrowException("Index::createGraphAndTree. Dimension is not specified.");
  }
  prop.indexType = NGT::Index::Property::IndexType::Graph;
  Index *idx     = 0;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  mkdir(database);
  idx = new NGT::GraphIndex(database, prop);
#else
  idx = new NGT::GraphIndex(prop);
#endif
  assert(idx != 0);
  StdOstreamRedirector redirector(redirect);
  redirector.begin();
  try {
    if (idx->getObjectSpace().isQintObjectType()) {
      idx->saveIndex(database);
      idx->close();
      auto append     = true;
      auto refinement = false;
      if (!dataFile.empty()) {
        appendFromTextObjectFile(database, dataFile, dataSize, append, refinement, prop.threadPoolSize);
      }
    } else {
      loadAndCreateIndex(*idx, database, dataFile, prop.threadPoolSize, dataSize);
    }
  } catch (Exception &err) {
    delete idx;
    redirector.end();
    throw err;
  }
  delete idx;
  redirector.end();
}

void NGT::Index::loadAndCreateIndex(Index &index, const string &database, const string &dataFile,
                                    size_t threadSize, size_t dataSize) {
  NGT::Timer timer;
  timer.start();
  if (dataFile.size() != 0) {
    index.load(dataFile, dataSize);
  } else {
    index.saveIndex(database);
    return;
  }
  timer.stop();
  cerr << "loadAndCreateIndex: Data loading time=" << timer.time << " (sec) " << timer.time * 1000.0
       << " (msec)" << endl;
  if (index.getObjectRepositorySize() == 0) {
    NGTThrowException("Index::create: Data file is empty.");
  }
  cerr << "# of objects=" << index.getObjectRepositorySize() - 1 << endl;
  timer.reset();
  timer.start();
  index.createIndex(threadSize);
  timer.stop();
  index.saveIndex(database);
  cerr << "Index creation time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
}

void NGT::Index::append(const string &database, const string &dataFile, size_t threadSize, size_t dataSize) {
  NGT::Index index(database);
  NGT::Timer timer;
  timer.start();
  if (dataFile.size() != 0) {
    index.append(dataFile, dataSize);
  }
  timer.stop();
  cerr << "append: Data loading time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << index.getObjectRepositorySize() - 1 << endl;
  timer.reset();
  timer.start();
  size_t nOfObjects            = index.getObjectSpace().getRepository().size();
  size_t endOfAppendedObjectID = (nOfObjects == 0 ? 1 : nOfObjects) + dataSize;
  index.createIndex(threadSize, endOfAppendedObjectID);
  timer.stop();
  index.saveIndex(database);
  cerr << "Index creation time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  return;
}

void NGT::Index::append(const string &database, const float *data, size_t dataSize, size_t threadSize) {
  NGT::Index index(database);
  NGT::Timer timer;
  timer.start();
  if (data != 0 && dataSize != 0) {
    index.append(data, dataSize);
  } else {
    NGTThrowException("Index::append: No data.");
  }
  timer.stop();
  cerr << "Data loading time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << index.getObjectRepositorySize() - 1 << endl;
  timer.reset();
  timer.start();
  index.createIndex(threadSize);
  timer.stop();
  index.saveIndex(database);
  cerr << "Index creation time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  return;
}

void NGT::Index::appendFromRefinementObjectFile(const std::string &indexPath, size_t threadSize) {
  NGT::Index index(indexPath);
  index.appendFromRefinementObjectFile();
  index.createIndex(threadSize);
  index.save();
  index.close();
}

void NGT::Index::appendFromRefinementObjectFile() {
  NGT::Property prop;
  getProperty(prop);
  float maxMag    = prop.maxMagnitude;
  bool maxMagSkip = false;
  if (maxMag > 0.0) maxMagSkip = true;
  auto &ros     = getRefinementObjectSpace();
  auto &rrepo   = ros.getRepository();
  size_t dim    = getDimension();
  auto dataSize = rrepo.size();
  std::vector<float> addedElement(dataSize);
  if (prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
    NGT::Timer timer;
    timer.start();
    for (size_t idx = 1; idx < rrepo.size(); idx++) {
      if (rrepo[idx] == 0) {
        continue;
      }
      std::vector<float> object;
      ros.getObject(idx, object);
      if (object.size() != dim) {
        if (object.size() == dim + 1) {
          object.resize(dim);
        } else {
          std::stringstream msg;
          msg << "Fatal inner error! iInvalid dimension. " << dim << ":" << object.size();
          ;
          NGTThrowException(msg);
        }
      }
      if (prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
        double mag = 0.0;
        for (auto &v : object) {
          mag += v * v;
        }
        if (!maxMagSkip && mag > maxMag) {
          maxMag = mag;
        }
        addedElement[idx] = mag;
      }
      if (idx % 2000000 == 0) {
        timer.stop();
        std::cerr << "processed " << static_cast<float>(idx) / 1000000.0 << "M objects."
                  << " maxMag=" << maxMag << " time=" << timer << std::endl;
        timer.restart();
      }
    }
    timer.stop();
    std::cerr << "time=" << timer << std::endl;
    std::cerr << "maxMag=" << maxMag << std::endl;
    std::cerr << "dataSize=" << dataSize << std::endl;
    if (prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
      if (static_cast<NGT::GraphIndex &>(getIndex()).property.maxMagnitude <= 0.0 && maxMag > 0.0) {
        static_cast<NGT::GraphIndex &>(getIndex()).property.maxMagnitude = maxMag;
      }
    }
  }
  if (getObjectSpace().isQintObjectType() && prop.clippingRate >= 0.0) {
    std::priority_queue<float> min;
    std::priority_queue<float, vector<float>, std::greater<float>> max;
    {
      NGT::Timer timer;
      timer.start();
      auto clippingSize = static_cast<float>(dataSize * dim) * prop.clippingRate;
      clippingSize      = clippingSize == 0 ? 1 : clippingSize;
      size_t counter    = 0;
      for (size_t idx = 1; idx < rrepo.size(); idx++) {
        if (rrepo[idx] == 0) continue;
        std::vector<float> object;
        ros.getObject(idx, object);
        if (object.size() != dim) object.resize(dim);
        if (getObjectSpace().isNormalizedDistance()) {
          ObjectSpace::normalize(object);
        }
        if (prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
          float v = maxMag - addedElement[idx];
          object.emplace_back(sqrt(v >= 0.0 ? v : 0.0));
        }
        for (auto &v : object) {
          if (max.size() < clippingSize) {
            max.push(v);
          } else if (max.top() <= v) {
            max.push(v);
            max.pop();
          }
          if (min.size() < clippingSize) {
            min.push(v);
          } else if (min.top() >= v) {
            min.push(v);
            min.pop();
          }
        }
        counter++;
      }
      std::cerr << "time=" << timer << std::endl;
      if (counter != 0) {
        std::cerr << "max:min=" << max.top() << ":" << min.top() << std::endl;
        setQuantizationFromMaxMin(max.top(), min.top());
      }
    }
  }
  {

    for (size_t idx = 1; idx < rrepo.size(); idx++) {
      if (rrepo[idx] == 0) continue;
      std::vector<float> object;
      ros.getObject(idx, object);
      if (object.size() != dim) object.resize(dim);
      if (prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
        object.emplace_back(sqrt(maxMag - addedElement[idx]));
      }
      append(object);
      if (idx + 1 != getObjectRepositorySize()) {
        std::stringstream msg;
        msg << "The object repository and refinement repository are inconsistent. " << idx + 1 << ":"
            << getObjectRepositorySize();
        NGTThrowException(msg);
      }
    }
  }
}

void NGT::Index::insertFromRefinementObjectFile() {
  NGT::Property prop;
  getProperty(prop);
  float maxMag = prop.maxMagnitude;
  if (prop.maxMagnitude <= 0.0) {
    std::stringstream msg;
    msg << "Max magnitude is not set yet. " << maxMag;
    NGTThrowException(msg);
  }
  auto &ros     = getRefinementObjectSpace();
  auto &rrepo   = ros.getRepository();
  auto &repo    = getObjectSpace().getRepository();
  size_t dim    = getDimension();
  auto dataSize = rrepo.size();
  std::vector<float> addedElement(dataSize);

  for (size_t idx = 1; idx < rrepo.size(); idx++) {
    if (rrepo[idx] == 0) continue;
    if (repo.size() > idx && repo[idx] != 0) continue;
    std::vector<float> object;
    ros.getObject(idx, object);
    if (object.size() != dim) {
      if (object.size() == dim + 1) {
        object.resize(dim);
      } else {
        std::stringstream msg;
        msg << "Fatal inner error! iInvalid dimension. " << dim << ":" << object.size();
        ;
        NGTThrowException(msg);
      }
    }
    if (prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
      double mag = 0.0;
      for (auto &v : object) {
        mag += v * v;
      }
      if (mag > maxMag) {
        maxMag = mag;
      }
      object.emplace_back(sqrt(maxMag - mag));
    }
    try {
      insert(idx, object);
    } catch (NGT::Exception &err) {
      std::stringstream msg;
      msg << "Cannot insert. " << idx << " " << err.what();
      NGTThrowException(msg);
    }
    if (idx + 1 > getObjectRepositorySize()) {
      std::stringstream msg;
      msg << "The object repository and refinement repository are inconsistent. " << idx + 1 << ":"
          << getObjectRepositorySize();
      NGTThrowException(msg);
    }
  }
}

void NGT::Index::appendFromTextObjectFile(const std::string &indexPath, const std::string &data,
                                          size_t dataSize, bool append, bool refinement, size_t threadSize) {
  //#define APPEND_TEST

  NGT::Index index(indexPath);
  index.appendFromTextObjectFile(data, dataSize, append, refinement);
  index.createIndex(threadSize);
  index.save();
  index.close();
}

void NGT::Index::appendFromTextObjectFile(const std::string &data, size_t dataSize, bool append,
                                          bool refinement) {
  StdOstreamRedirector redirector(redirect);
  redirector.begin();
  NGT::Property prop;
  getProperty(prop);
  float maxMag    = prop.maxMagnitude;
  bool maxMagSkip = false;
  if (maxMag > 0.0) maxMagSkip = true;
  std::vector<float> addedElement;
  size_t dim = prop.dimension;
  bool warn  = false;
  if (append && prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
    NGT::Timer timer;
    timer.start();
    ifstream is(data);
    if (!is) {
      std::stringstream msg;
      msg << "Cannot open the specified data file. " << data;
      NGTThrowException(msg);
    }
    std::string line;
    size_t counter = 0;
    while (getline(is, line)) {
      if (is.eof()) break;
      if (dataSize > 0 && counter > dataSize) break;
      vector<float> object;
      vector<string> tokens;
      NGT::Common::tokenize(line, tokens, "\t, ");
      if (tokens.back() == "") tokens.pop_back();
      if (tokens.size() > dim && warn == false) {
        warn = true;
        std::cerr << "Warning! Invalid dimension of the specified data. The specified data is "
                  << tokens.size() << ". The index is " << dim << "." << std::endl;
        std::cerr << "Cut the tail." << std::endl;
      }
      if (tokens.size() < dim) {
        std::stringstream msg;
        msg << "The dimensions are not inconsist. " << counter << ":" << dim << "x" << tokens.size() << data;
        NGTThrowException(msg);
      }
      if (tokens.size() > dim) {
        tokens.resize(dim);
      }
      if (prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
        double mag = 0.0;
        for (auto &vstr : tokens) {
          auto v = NGT::Common::strtof(vstr);
          mag += v * v;
        }
        if (!maxMagSkip && mag > maxMag) {
          maxMag = mag;
        }
        addedElement.emplace_back(mag);
      }
      counter++;
      if (counter % 2000000 == 0) {
        timer.stop();
        std::cerr << "processed " << static_cast<float>(counter) / 1000000.0 << "M objects."
                  << " maxMag=" << maxMag << " time=" << timer << std::endl;
        timer.restart();
      }
    }
    timer.stop();
    dataSize = counter;
    std::cerr << "time=" << timer << std::endl;
    std::cerr << "maxMag=" << maxMag << std::endl;
    std::cerr << "dataSize=" << dataSize << std::endl;
    if (prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
      if (static_cast<NGT::GraphIndex &>(getIndex()).property.maxMagnitude <= 0.0 && maxMag > 0.0) {
        static_cast<NGT::GraphIndex &>(getIndex()).property.maxMagnitude = maxMag;
      }
    }
  }
  if (append && getObjectSpace().isQintObjectType() && prop.clippingRate >= 0.0) {
    std::priority_queue<float> min;
    std::priority_queue<float, vector<float>, std::greater<float>> max;
    {
      NGT::Timer timer;
      timer.start();
      ifstream is(data);
      if (!is) {
        std::stringstream msg;
        msg << "Cannot open the specified data file. " << data;
        NGTThrowException(msg);
      }
      auto clippingSize = static_cast<float>(dataSize * dim) * prop.clippingRate;
      clippingSize      = clippingSize == 0 ? 1 : clippingSize;
      std::string line;
      size_t counter = 0;
      while (getline(is, line)) {
        if (is.eof()) break;
        if (dataSize > 0 && counter > dataSize) break;
        vector<float> object;
        vector<string> tokens;
        NGT::Common::tokenize(line, tokens, "\t, ");
        if (tokens.back() == "") tokens.pop_back();
        if (tokens.size() > dim && warn == false) {
          warn = true;
          std::cerr << "Warning! Invalid dimension of the specified data. The specified data is "
                    << tokens.size() << ". The index is " << dim << "." << std::endl;
          std::cerr << "Cut the tail." << std::endl;
        }
        for (auto &vstr : tokens) {
          auto v = NGT::Common::strtof(vstr);
          object.emplace_back(v);
        }
        if (object.size() > dim) {
          object.resize(dim);
        }
        if (getObjectSpace().isNormalizedDistance()) {
          ObjectSpace::normalize(object);
        }
        if (prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
          float v = maxMag - addedElement[counter];
          object.emplace_back(sqrt(v >= 0.0 ? v : 0.0));
        }
        for (auto &v : object) {
          if (max.size() < clippingSize) {
            max.push(v);
          } else if (max.top() <= v) {
            max.push(v);
            max.pop();
          }
          if (min.size() < clippingSize) {
            min.push(v);
          } else if (min.top() >= v) {
            min.push(v);
            min.pop();
          }
        }
        counter++;
      }
      timer.stop();
      std::cerr << "time=" << timer << std::endl;
      if (counter != 0) {
        std::cerr << "max:min=" << max.top() << ":" << min.top() << std::endl;
        setQuantizationFromMaxMin(max.top(), min.top());
      }
    }
  }
  if (append || refinement) {

    ifstream is(data);
    if (!is) {
      std::stringstream msg;
      msg << "Cannot open the specified data file. " << data;
      NGTThrowException(msg);
    }
    std::string line;
    size_t counter = 0;
    while (getline(is, line)) {
      if (is.eof()) break;
      if (dataSize > 0 && counter > dataSize) break;
      vector<float> object;
      vector<string> tokens;
      NGT::Common::tokenize(line, tokens, "\t, ");
      if (tokens.back() == "") tokens.pop_back();
      if (tokens.size() > dim && warn == false) {
        warn = true;
        std::cerr << "Warning! Invalid dimension of the specified data. The specified data is "
                  << tokens.size() << ". The index is " << dim << "." << std::endl;
        std::cerr << "Cut the tail." << std::endl;
      }
      for (auto &vstr : tokens) {
        auto v = NGT::Common::strtof(vstr);
        object.emplace_back(v);
      }
      if (object.size() > dim) {
        object.resize(dim);
      }
#ifdef NGT_REFINEMENT
      if (refinement) {
        appendToRefinement(object);
      }
#endif
      if (append) {
        if (prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct && maxMag > 0.0) {
          float v = maxMag - addedElement[counter];
          object.emplace_back(sqrt(v >= 0.0 ? v : 0.0));
        }
        NGT::Index::append(object);
      }
      counter++;
    }
  }
  redirector.end();
}

void NGT::Index::appendFromBinaryObjectFile(const std::string &indexPath, const std::string &data,
                                            size_t dataSize, bool append, bool refinement,
                                            size_t threadSize) {
  NGT::Index index(indexPath);
  index.appendFromBinaryObjectFile(data, dataSize, append, refinement);
  index.createIndex(threadSize);
  index.save();
  index.close();
}

void NGT::Index::appendFromBinaryObjectFile(const std::string &data, size_t dataSize, bool append,
                                            bool refinement) {
  NGT::Property prop;
  getProperty(prop);
  float maxMag    = prop.maxMagnitude;
  bool maxMagSkip = false;
  if (maxMag > 0.0) maxMagSkip = true;
  std::vector<float> addedElement;
  size_t dim = 0;
  if (append && prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
    NGT::Timer timer;
    timer.start();
    StaticObjectFileLoader loader(data);
    size_t counter = 0;
    while (!loader.isEmpty()) {
      if (dataSize > 0 && counter > dataSize) break;
      auto object = loader.getObject();
      if (dim == 0) {
        dim = object.size();
      } else if (dim != object.size()) {
        std::stringstream msg;
        msg << "The dimensions are not inconsist. " << counter << ":" << dim << "x" << object.size() << data;
        NGTThrowException(msg);
      }
      if (prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
        double mag = 0.0;
        for (auto &v : object) {
          mag += v * v;
        }
        if (!maxMagSkip && mag > maxMag) {
          maxMag = mag;
        }
        addedElement.emplace_back(mag);
      }
      counter++;
      if (counter % 2000000 == 0) {
        timer.stop();
        std::cerr << "processed " << static_cast<float>(counter) / 1000000.0 << "M objects."
                  << " maxMag=" << maxMag << " time=" << timer << std::endl;
        timer.restart();
      }
    }
    timer.stop();
    dataSize = counter;
    std::cerr << "time=" << timer << std::endl;
    std::cerr << "maxMag=" << maxMag << std::endl;
    std::cerr << "dataSize=" << dataSize << std::endl;
    if (prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
      if (static_cast<NGT::GraphIndex &>(getIndex()).property.maxMagnitude <= 0.0 && maxMag > 0.0) {
        static_cast<NGT::GraphIndex &>(getIndex()).property.maxMagnitude = maxMag;
      }
    }
  }
  if (append && getObjectSpace().isQintObjectType() && prop.clippingRate >= 0.0) {
    std::priority_queue<float> min;
    std::priority_queue<float, vector<float>, std::greater<float>> max;
    {
      NGT::Timer timer;
      timer.start();
      auto clippingSize = static_cast<float>(dataSize * dim) * prop.clippingRate;
      clippingSize      = clippingSize == 0 ? 1 : clippingSize;
      StaticObjectFileLoader loader(data);
      size_t counter = 0;
      while (!loader.isEmpty()) {
        if (dataSize > 0 && counter > dataSize) break;
        auto object = loader.getObject();
        if (prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
          float v = maxMag - addedElement[counter];
          object.emplace_back(sqrt(v >= 0.0 ? v : 0.0));
        }
        if (getObjectSpace().isNormalizedDistance()) {
          ObjectSpace::normalize(object);
        }
        for (auto &v : object) {
          if (max.size() < clippingSize) {
            max.push(v);
          } else if (max.top() <= v) {
            max.push(v);
            max.pop();
          }
          if (min.size() < clippingSize) {
            min.push(v);
          } else if (min.top() >= v) {
            min.push(v);
            min.pop();
          }
        }
        counter++;
      }
      timer.stop();
      std::cerr << "time=" << timer << std::endl;
      if (counter != 0) {
        std::cerr << "max:min=" << max.top() << ":" << min.top() << std::endl;
        setQuantizationFromMaxMin(max.top(), min.top());
      }
    }
  }
  if (append || refinement) {
    StaticObjectFileLoader loader(data);
    size_t counter = 0;
    while (!loader.isEmpty()) {
      if (dataSize > 0 && counter > dataSize) break;
      auto object = loader.getObject();
#ifdef NGT_REFINEMENT
      if (refinement) {
        appendToRefinement(object);
      }
#endif
      if (append) {
        if (prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
          object.emplace_back(sqrt(maxMag - addedElement[counter]));
        }
        NGT::Index::append(object);
      }
      counter++;
    }
  }
}

void NGT::Index::remove(const string &database, vector<ObjectID> &objects, bool force) {
  NGT::Index index(database);
  NGT::Timer timer;
  timer.start();
  for (vector<ObjectID>::iterator i = objects.begin(); i != objects.end(); i++) {
    try {
      index.remove(*i, force);
    } catch (Exception &err) {
      cerr << "Warning: Cannot remove the node. ID=" << *i << " : " << err.what() << endl;
      continue;
    }
  }
  timer.stop();
  cerr << "Data removing time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << index.getObjectRepositorySize() - 1 << endl;
  index.saveIndex(database);
  return;
}

void NGT::Index::importIndex(const string &database, const string &file) {
  Index *idx = 0;
  NGT::Property property;
  property.importProperty(file);
  NGT::Timer timer;
  timer.start();
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  property.databaseType = NGT::Index::Property::DatabaseType::MemoryMappedFile;
  mkdir(database);
#else
  property.databaseType = NGT::Index::Property::DatabaseType::Memory;
#endif
  if (property.indexType == NGT::Index::Property::IndexType::GraphAndTree) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    idx = new NGT::GraphAndTreeIndex(database, property);
#else
    idx = new NGT::GraphAndTreeIndex(property);
#endif
    assert(idx != 0);
  } else if (property.indexType == NGT::Index::Property::IndexType::Graph) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    idx = new NGT::GraphIndex(database, property);
#else
    idx = new NGT::GraphIndex(property);
#endif
    assert(idx != 0);
  } else {
    NGTThrowException("Index::Open: Not found IndexType in property file.");
  }
  idx->importIndex(file);
  timer.stop();
  cerr << "Data importing time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << idx->getObjectRepositorySize() - 1 << endl;
  idx->saveIndex(database);
  delete idx;
}

void NGT::Index::exportIndex(const string &database, const string &file) {
  NGT::Index idx(database);
  NGT::Timer timer;
  timer.start();
  idx.exportIndex(file);
  timer.stop();
  cerr << "Data exporting time=" << timer.time << " (sec) " << timer.time * 1000.0 << " (msec)" << endl;
  cerr << "# of objects=" << idx.getObjectRepositorySize() - 1 << endl;
}

void NGT::Index::searchUsingOnlyGraph(NGT::SearchContainer &sc) {
  static_cast<GraphIndex &>(getIndex()).search(sc);
}

void NGT::Index::searchUsingOnlyGraph(NGT::SearchQuery &searchQuery) {
  static_cast<GraphIndex &>(getIndex()).GraphIndex::search(searchQuery);
}

std::vector<float> NGT::Index::makeSparseObject(std::vector<uint32_t> &object) {
  if (static_cast<NGT::GraphIndex &>(getIndex()).getProperty().distanceType !=
      NGT::ObjectSpace::DistanceType::DistanceTypeSparseJaccard) {
    NGTThrowException("NGT::Index::makeSparseObject: Not sparse jaccard.");
  }
  size_t dimension = getObjectSpace().getDimension();
  if (object.size() + 1 > dimension) {
    std::stringstream msg;
    dimension = object.size() + 1;
  }
  std::vector<float> obj(dimension, 0.0);
  for (size_t i = 0; i < object.size(); i++) {
    float fv = *reinterpret_cast<float *>(&object[i]);
    obj[i]   = fv;
  }
  return obj;
}

void NGT::Index::setQuantizationFromMaxMin(float max, float min) {

  float offset;
  float scale;
  if (getObjectSpace().getObjectType() == typeid(NGT::qsint8)) {
    offset = 0.0;
    scale  = std::max(fabs(max), fabs(min));
  } else {
    offset = min;
    scale  = max - offset;
  }
  setQuantization(scale, offset);
}

void NGT::Index::setQuantization(float scale, float offset) {
  static_cast<NGT::GraphIndex &>(getIndex()).property.quantizationScale  = scale;
  static_cast<NGT::GraphIndex &>(getIndex()).property.quantizationOffset = offset;
  getObjectSpace().setQuantization(scale, offset);
}

void NGT::Index::extractInsertionOrder(InsertionOrder &insertionOrder) {
  static_cast<NGT::GraphIndex &>(getIndex()).extractInsertionOrder(insertionOrder);
}

void NGT::Index::createIndex(size_t threadNumber, size_t sizeOfRepository) {
  StdOstreamRedirector redirector(redirect);
  redirector.begin();
  try {
    InsertionOrder insertionOrder;
    NGT::Property prop;
    getProperty(prop);
    if (prop.objectType == NGT::ObjectSpace::ObjectType::Qsuint8) {
      auto &ros = getRefinementObjectSpace();
      auto &os  = getObjectSpace();
      if (&ros != 0 && ros.getRepository().size() > os.getRepository().size()) {
        if (os.getRepository().size() <= 1) {
          if (ros.getRepository().size() < 100) {
            std::cerr << "Warning! # of refinement objects is too small. " << ros.getRepository().size()
                      << std::endl;
          }
          appendFromRefinementObjectFile();
        } else {
          if (prop.quantizationScale <= 0.0) {
            stringstream msg;
            msg << "Fatal inner error! Scalar quantization parameters are not set yet. "
                << prop.quantizationScale << ":" << prop.quantizationOffset;
            NGTThrowException(msg);
          }
          insertFromRefinementObjectFile();
        }
      }
    } else {
      if (prop.distanceType == ObjectSpace::DistanceTypeInnerProduct) {
        size_t beginId                        = 1;
        NGT::GraphRepository &graphRepository = static_cast<NGT::GraphIndex &>(getIndex()).repository;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
        auto &graphNodes       = static_cast<PersistentRepository<GraphNode> &>(graphRepository);
        auto &graphNodeVectors = reinterpret_cast<PersistentRepository<void> &>(graphNodes);
#else
        auto &graphNodes       = static_cast<Repository<GraphNode> &>(graphRepository);
        auto &graphNodeVectors = reinterpret_cast<Repository<void> &>(graphNodes);
#endif
        if (prop.maxMagnitude <= 0.0) {
          getObjectSpace().setMagnitude(prop.maxMagnitude, graphNodeVectors, beginId);
        } else {
          auto maxMag = getObjectSpace().computeMaxMagnitude(beginId);
          static_cast<NGT::GraphIndex &>(getIndex()).property.maxMagnitude = maxMag;
          getObjectSpace().setMagnitude(maxMag, graphNodeVectors, beginId);
        }
      }
    }
    if (prop.nOfNeighborsForInsertionOrder != 0) {
      insertionOrder.nOfNeighboringNodes = prop.nOfNeighborsForInsertionOrder;
      insertionOrder.epsilon             = prop.epsilonForInsertionOrder;
      extractInsertionOrder(insertionOrder);
    }
    createIndexWithInsertionOrder(insertionOrder, threadNumber, sizeOfRepository);
  } catch (Exception &err) {
    redirector.end();
    throw err;
  }
  redirector.end();
}

void NGT::Index::Property::set(NGT::Property &prop) {
  if (prop.dimension != -1) dimension = prop.dimension;
  if (prop.threadPoolSize != -1) threadPoolSize = prop.threadPoolSize;
  if (prop.objectType != ObjectSpace::ObjectTypeNone) objectType = prop.objectType;
#ifdef NGT_REFINEMENT
  if (prop.refinementObjectType != ObjectSpace::ObjectTypeNone)
    refinementObjectType = prop.refinementObjectType;
#endif
  if (prop.distanceType != DistanceType::DistanceTypeNone) distanceType = prop.distanceType;
  if (prop.indexType != IndexTypeNone) indexType = prop.indexType;
  if (prop.databaseType != DatabaseTypeNone) databaseType = prop.databaseType;
  if (prop.objectAlignment != ObjectAlignmentNone) objectAlignment = prop.objectAlignment;
  if (prop.pathAdjustmentInterval != -1) pathAdjustmentInterval = prop.pathAdjustmentInterval;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  if (prop.graphSharedMemorySize != -1) graphSharedMemorySize = prop.graphSharedMemorySize;
  if (prop.treeSharedMemorySize != -1) treeSharedMemorySize = prop.treeSharedMemorySize;
  if (prop.objectSharedMemorySize != -1) objectSharedMemorySize = prop.objectSharedMemorySize;
#endif
  if (prop.prefetchOffset != -1) prefetchOffset = prop.prefetchOffset;
  if (prop.prefetchSize != -1) prefetchSize = prop.prefetchSize;
  if (prop.accuracyTable != "") accuracyTable = prop.accuracyTable;
  if (prop.maxMagnitude != -1) maxMagnitude = prop.maxMagnitude;
  if (prop.quantizationScale != -1) quantizationScale = prop.quantizationScale;
  if (prop.quantizationOffset != -1) quantizationOffset = prop.quantizationOffset;
  if (prop.clippingRate != -1) clippingRate = prop.clippingRate;
  if (prop.nOfNeighborsForInsertionOrder != -1)
    nOfNeighborsForInsertionOrder = prop.nOfNeighborsForInsertionOrder;
  if (prop.epsilonForInsertionOrder != -1) epsilonForInsertionOrder = prop.epsilonForInsertionOrder;
}

void NGT::Index::Property::get(NGT::Property &prop) {
  prop.dimension      = dimension;
  prop.threadPoolSize = threadPoolSize;
  prop.objectType     = objectType;
#ifdef NGT_REFINEMENT
  prop.refinementObjectType = refinementObjectType;
#endif
  prop.distanceType           = distanceType;
  prop.indexType              = indexType;
  prop.databaseType           = databaseType;
  prop.pathAdjustmentInterval = pathAdjustmentInterval;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
  prop.graphSharedMemorySize  = graphSharedMemorySize;
  prop.treeSharedMemorySize   = treeSharedMemorySize;
  prop.objectSharedMemorySize = objectSharedMemorySize;
#endif
  prop.prefetchOffset                = prefetchOffset;
  prop.prefetchSize                  = prefetchSize;
  prop.accuracyTable                 = accuracyTable;
  prop.maxMagnitude                  = maxMagnitude;
  prop.quantizationScale             = quantizationScale;
  prop.quantizationOffset            = quantizationOffset;
  prop.clippingRate                  = clippingRate;
  prop.nOfNeighborsForInsertionOrder = nOfNeighborsForInsertionOrder;
  prop.epsilonForInsertionOrder      = epsilonForInsertionOrder;
}

class CreateIndexJob {
 public:
  CreateIndexJob() {}
  CreateIndexJob &operator=(const CreateIndexJob &d) {
    id       = d.id;
    results  = d.results;
    object   = d.object;
    batchIdx = d.batchIdx;
    return *this;
  }
  friend bool operator<(const CreateIndexJob &ja, const CreateIndexJob &jb) {
    return ja.batchIdx < jb.batchIdx;
  }
  NGT::ObjectID id;
  NGT::Object *object; // this will be a node of the graph later.
  NGT::ObjectDistances *results;
  size_t batchIdx;
};

class CreateIndexSharedData {
 public:
  CreateIndexSharedData(NGT::GraphIndex &nngt) : graphIndex(nngt) {}
  NGT::GraphIndex &graphIndex;
};

class CreateIndexThread : public NGT::Thread {
 public:
  CreateIndexThread() {}
  virtual ~CreateIndexThread() {}
  virtual int run();
};

typedef NGT::ThreadPool<CreateIndexJob, CreateIndexSharedData *, CreateIndexThread> CreateIndexThreadPool;

int CreateIndexThread::run() {

  NGT::ThreadPool<CreateIndexJob, CreateIndexSharedData *, CreateIndexThread>::Thread &poolThread =
      (NGT::ThreadPool<CreateIndexJob, CreateIndexSharedData *, CreateIndexThread>::Thread &)*this;

  CreateIndexSharedData &sd   = *poolThread.getSharedData();
  NGT::GraphIndex &graphIndex = sd.graphIndex;

  for (;;) {
    CreateIndexJob job;
    try {
      poolThread.getInputJobQueue().popFront(job);
    } catch (NGT::ThreadTerminationException &err) {
      break;
    } catch (NGT::Exception &err) {
      cerr << "CreateIndex::search:Error! popFront " << err.what() << endl;
      break;
    }
    ObjectDistances *rs = new ObjectDistances;
    Object &obj         = *job.object;
    try {
      if (graphIndex.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeKNNG) {
        graphIndex.searchForKNNGInsertion(obj, job.id, *rs); // linear search
      } else {
        graphIndex.searchForNNGInsertion(obj, *rs);
      }
    } catch (NGT::Exception &err) {
      stringstream msg;
      msg << "CreateIndex::search:Fatal error! ID=" << job.id << " " << err.what();
      NGTThrowException(msg);
    }
    job.results = rs;
    poolThread.getOutputJobQueue().pushBack(job);
  }

  return 0;
}

class BuildTimeController {
 public:
  BuildTimeController(GraphIndex &graph, NeighborhoodGraph::Property &prop) : property(prop) {
    noOfInsertedObjects            = graph.objectSpace->getRepository().size() - graph.repository.size();
    interval                       = 10000;
    count                          = interval;
    edgeSizeSave                   = property.edgeSizeForCreation;
    insertionRadiusCoefficientSave = property.insertionRadiusCoefficient;
    buildTimeLimit                 = property.buildTimeLimit;
    time                           = 0.0;
    timer.start();
  }
  ~BuildTimeController() {
    property.edgeSizeForCreation        = edgeSizeSave;
    property.insertionRadiusCoefficient = insertionRadiusCoefficientSave;
  }
  void adjustEdgeSize(size_t c) {
    if (buildTimeLimit > 0.0 && count <= c) {
      timer.stop();
      double estimatedTime = time + timer.time / interval * (noOfInsertedObjects - count);
      estimatedTime /= 60 * 60; // hour
      const size_t edgeInterval  = 5;
      const int minimumEdge      = 5;
      const float radiusInterval = 0.02;
      if (estimatedTime > buildTimeLimit) {
        if (property.insertionRadiusCoefficient - radiusInterval >= 1.0) {
          property.insertionRadiusCoefficient -= radiusInterval;
        } else {
          property.edgeSizeForCreation -= edgeInterval;
          if (property.edgeSizeForCreation < minimumEdge) {
            property.edgeSizeForCreation = minimumEdge;
          }
        }
      }
      time += timer.time;
      count += interval;
      timer.start();
    }
  }

  size_t noOfInsertedObjects;
  size_t interval;
  size_t count;
  size_t edgeSizeSave;
  double insertionRadiusCoefficientSave;
  Timer timer;
  double time;
  double buildTimeLimit;
  NeighborhoodGraph::Property &property;
};

void NGT::GraphIndex::constructObjectSpace(NGT::Property &prop) {
  assert(prop.dimension != 0);
  size_t dimension = prop.dimension;
  if (prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeSparseJaccard ||
      prop.distanceType == NGT::ObjectSpace::DistanceType::DistanceTypeInnerProduct) {
    dimension++;
  }

  switch (prop.objectType) {
  case NGT::ObjectSpace::ObjectType::Float:
    objectSpace = new ObjectSpaceRepository<float, double>(dimension, typeid(float), prop.distanceType);
    break;
  case NGT::ObjectSpace::ObjectType::Uint8:
    objectSpace =
        new ObjectSpaceRepository<unsigned char, int>(dimension, typeid(uint8_t), prop.distanceType);
    break;
#ifdef NGT_HALF_FLOAT
  case NGT::ObjectSpace::ObjectType::Float16:
    objectSpace = new ObjectSpaceRepository<float16, float>(dimension, typeid(float16), prop.distanceType);
    break;
#endif
  case NGT::ObjectSpace::ObjectType::Qsuint8:
    objectSpace = new ObjectSpaceRepository<qsint8, float>(dimension, typeid(qsint8), prop.distanceType);
    break;
#ifdef NGT_PQ4
  case NGT::ObjectSpace::ObjectType::Qint4:
    objectSpace = new ObjectSpaceRepository<qint4, float>(dimension, typeid(qint4), prop.distanceType);
    break;
#endif
  default:
    stringstream msg;
    msg << "Invalid Object Type in the property. " << prop.objectType;
    NGTThrowException(msg);
  }
  objectSpace->setQuantization(prop.quantizationScale, prop.quantizationOffset);
#ifdef NGT_REFINEMENT
  auto dtype = prop.distanceType;
  dtype      = dtype == ObjectSpace::DistanceTypeInnerProduct ? ObjectSpace::DistanceTypeDotProduct
                                                              : prop.distanceType;
  switch (prop.refinementObjectType) {
  case NGT::ObjectSpace::ObjectType::Float:
    refinementObjectSpace = new ObjectSpaceRepository<float, double>(dimension, typeid(float), dtype);
    break;
  case NGT::ObjectSpace::ObjectType::Uint8:
    refinementObjectSpace = new ObjectSpaceRepository<unsigned char, int>(dimension, typeid(uint8_t), dtype);
    break;
#ifdef NGT_HALF_FLOAT
  case NGT::ObjectSpace::ObjectType::Float16:
    refinementObjectSpace = new ObjectSpaceRepository<float16, float>(dimension, typeid(float16), dtype);
    break;
#endif
#ifdef NGT_BFLOAT
  case NGT::ObjectSpace::ObjectType::Bfloat16:
    refinementObjectSpace = new ObjectSpaceRepository<bfloat16, float>(dimension, typeid(bfloat16), dtype);
    break;
#endif
  default:
    stringstream msg;
    msg << "Invalid Refinement Object Type in the property. " << prop.refinementObjectType;
    NGTThrowException(msg);
  }
#endif
}

void NGT::GraphIndex::loadGraph(const string &ifile, NGT::GraphRepository &graph) {
  ifstream isg(ifile + "/grp");
  graph.deserialize(isg);
}

void NGT::GraphIndex::loadIndex(const string &ifile, bool readOnly, NGT::Index::OpenType openType) {
  if ((openType & NGT::Index::OpenTypeObjectDisabled) == 0) {
    objectSpace->deserialize(ifile + "/obj");
  }
#ifdef NGT_PQ4
  if (objectSpace != 0) {
    objectSpace->openQuantizer(ifile);
  }
#endif
#ifdef NGT_REFINEMENT
  try {
    refinementObjectSpace->deserialize(ifile + "/robj");
  } catch (Exception &err) {
    std::cerr << "Warning. Cannot open the refinment objects. " << err.what() << std::endl;
  }
#endif
  if ((openType & NGT::Index::OpenTypeGraphDisabled) == 0) {
    try {
#ifdef NGT_GRAPH_READ_ONLY_GRAPH
      if (readOnly) {
        GraphIndex::NeighborhoodGraph::loadSearchGraph(ifile);
      } else {
        loadGraph(ifile, repository);
        checkEdgeLengths(1000);
      }
#else
      loadGraph(ifile, repository);
      checkEdgeLengths(1000);
#endif
    } catch (Exception &err) {
      std::stringstream msg;
      msg << "Fatal error! Cannot load the graph. :" << err.what();
      NGTThrowException(msg);
    }
  }
}

void NGT::GraphIndex::saveProperty(const std::string &file) { NGT::Property::save(*this, file); }

void NGT::GraphIndex::exportProperty(const std::string &file) { NGT::Property::exportProperty(*this, file); }

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
NGT::GraphIndex::GraphIndex(const string &allocator, bool rdonly) : readOnly(rdonly) {
  NGT::Property prop;
  prop.load(allocator);
  if (prop.databaseType != NGT::Index::Property::DatabaseType::MemoryMappedFile) {
    NGTThrowException("GraphIndex: Cannot open. Not memory mapped file type.");
  }
  initialize(allocator, prop);
#ifdef NGT_GRAPH_READ_ONLY_GRAPH
  searchUnupdatableGraph = NeighborhoodGraph::Search::getMethod(prop.distanceType, prop.objectType,
                                                                objectSpace->getRepository().size());
#endif
}

NGT::GraphAndTreeIndex::GraphAndTreeIndex(const string &allocator, NGT::Property &prop)
    : GraphIndex(allocator, prop) {
  initialize(allocator, prop.treeSharedMemorySize);
}

void GraphAndTreeIndex::createTreeIndex() {
  ObjectRepository &fr = GraphIndex::objectSpace->getRepository();
  for (size_t id = 0; id < fr.size(); id++) {
    if (id % 100000 == 0) {
      cerr << " Processed id=" << id << endl;
    }
    if (fr.isEmpty(id)) {
      continue;
    }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    Object *f = GraphIndex::objectSpace->allocateObject(*fr[id]);
    DVPTree::InsertContainer tiobj(*f, id);
#else
    DVPTree::InsertContainer tiobj(*fr[id], id);
#endif
    try {
      DVPTree::insert(tiobj);
    } catch (Exception &err) {
      cerr << "GraphAndTreeIndex::createTreeIndex: Warning. ID=" << id << ":";
      cerr << err.what() << " continue.." << endl;
    }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    GraphIndex::objectSpace->deleteObject(f);
#endif
  }
}

void NGT::GraphIndex::initialize(const string &allocator, NGT::Property &prop) {
  constructObjectSpace(prop);
  repository.open(allocator + "/grp", prop.graphSharedMemorySize);
  objectSpace->open(allocator + "/obj", prop.objectSharedMemorySize);
#ifdef NGT_REFINEMENT
  refinementObjectSpace->open(allocator + "/robj", prop.objectSharedMemorySize);
#endif
  setProperty(prop);
}
#else // NGT_SHARED_MEMORY_ALLOCATOR
NGT::GraphIndex::GraphIndex(const string &database, bool rdOnly, NGT::Index::OpenType openType)
    : readOnly(rdOnly) {
  NGT::Property prop;
  prop.load(database);
  if (prop.databaseType != NGT::Index::Property::DatabaseType::Memory) {
    NGTThrowException("GraphIndex: Cannot open. Not memory type.");
  }
  assert(prop.dimension != 0);
  initialize(prop);
  loadIndex(database, readOnly, openType);
#ifdef NGT_GRAPH_READ_ONLY_GRAPH
  if (prop.searchType == "Large") {
    searchUnupdatableGraph =
        NeighborhoodGraph::Search::getMethod(prop.distanceType, prop.objectType, 10000000);
  } else if (prop.searchType == "Small") {
    searchUnupdatableGraph = NeighborhoodGraph::Search::getMethod(prop.distanceType, prop.objectType, 0);
  } else {
    searchUnupdatableGraph = NeighborhoodGraph::Search::getMethod(prop.distanceType, prop.objectType,
                                                                  objectSpace->getRepository().size());
  }
#endif
}
#endif

void GraphIndex::extractSparseness(InsertionOrder &insertionOrder) {
  if (getNumberOfIndexedObjects() == 0) {
    NGTThrowException("extractInsertionOrder: No indexed objects.");
  }
  auto nOfThreads =
      insertionOrder.nOfThreads == 0 ? std::thread::hardware_concurrency() : insertionOrder.nOfThreads;
  NGT::Timer timer;
  timer.start();
  std::cerr << "extractInsertionOrder" << std::endl;
  std::cerr << "VM size=" << NGT::Common::getProcessVmSizeStr() << std::endl;
  std::cerr << "Peak VM size=" << NGT::Common::getProcessVmPeakStr() << std::endl;
  std::cerr << "searching..." << std::endl;

  if (getObjectRepositorySize() != getGraphRepositorySize()) {
    std::stringstream msg;
    msg << "extractInsertionOrder: # of objects and # of indexed objects are not consistent. "
        << getObjectRepositorySize() << ":" << getGraphRepositorySize();
    NGTThrowException(msg);
  }

  omp_set_num_threads(nOfThreads);

  std::cerr << "search size=" << insertionOrder.nOfNeighboringNodes << std::endl;
  std::vector<uint32_t> counter(nOfThreads);
  std::vector<uint32_t> indegrees[nOfThreads];
  for (size_t tidx = 0; tidx < nOfThreads; tidx++) {
    indegrees[tidx].resize(getGraphRepositorySize());
  }
  std::vector<std::pair<float, uint32_t>> length;
  length.resize(getObjectRepositorySize());
#pragma omp parallel for
  for (NGT::ObjectID query = 1; query < getObjectRepositorySize(); query++) {
    auto thdID = omp_get_thread_num();
    counter[thdID]++;
    if (query % 100000 == 0) {
      size_t n = 0;
      for (auto &c : counter)
        n += c;
      timer.stop();
      std::cerr << "# of the processed objects=" << n << " VM size=" << NGT::Common::getProcessVmSizeStr()
                << " Peak VM size=" << NGT::Common::getProcessVmPeakStr() << " Time=" << timer << std::endl;
      timer.restart();
    }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    NGT::Object *object = getObjectSpace().allocateObject(*getObjectSpace().getRepository().get(query));
#else
    NGT::Object *object = getObjectSpace().getRepository().get(query);
#endif
    {
      NGT::SearchContainer sc(*object);
      NGT::ObjectDistances objects;
      sc.setResults(&objects);
      sc.setSize(insertionOrder.nOfNeighboringNodes);
      sc.setEpsilon(insertionOrder.epsilon);
      sc.setEdgeSize(-2);
      NGT::Timer timer;
      try {
        timer.start();
        search(sc);
        timer.stop();
      } catch (NGT::Exception &err) {
        std::cerr << "extractSparseness: Warning! " << err.what() << ":" << query << std::endl;
      }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      getObjectSpace().deleteObject(object);
#endif
      float len    = 0.0;
      size_t count = 0;
      if (objects.size() == 0) {
        std::stringstream msg;
        msg << "extractInsertionOrder: Error! # of the searched objects is zero. " << query << ":"
            << getPath() << std::endl;
        NGTThrowException(msg);
      }
      for (size_t i = 0; i < objects.size(); i++) {
        if (query == objects[i].id) continue;
        len += objects[i].distance;
        count++;
        if (objects[i].id >= indegrees[thdID].size()) {
          std::cerr << "too large. " << objects[i].id << ":" << indegrees[thdID].size() << std::endl;
          exit(1);
        }
        indegrees[thdID][objects[i].id]++;
      }
      length[query].first  = len / count;
      length[query].second = query;
    }
  }

  std::sort(length.begin(), length.end());

  size_t max = 0;
  for (NGT::ObjectID id = 1; id < getObjectRepositorySize(); id++) {
    for (size_t tidx = 1; tidx < nOfThreads; tidx++) {
      indegrees[0][id] += indegrees[tidx][id];
    }
    if (indegrees[0][id] > max) max = indegrees[0][id];
  }
  std::cerr << "max=" << max << std::endl;
  if (insertionOrder.indegreeOrder) {
    std::vector<std::vector<uint32_t>> sortedIndegrees;
    sortedIndegrees.resize(max + 1);
    for (uint32_t oid = 1; oid < indegrees[0].size(); oid++) {

      sortedIndegrees[indegrees[0][oid]].push_back(oid);
    }
    {
      size_t c = 0;
      insertionOrder.reserve(getObjectRepositorySize());
      for (uint32_t ind = 0; ind < sortedIndegrees.size(); ind++) {
        c += ind * sortedIndegrees[ind].size();
        if (sortedIndegrees[ind].size() != 0) {
          for (auto &id : sortedIndegrees[ind]) {
            insertionOrder.push_back(id);
          }
        }
      }
      std::cerr << "total number of the incoming edges=" << c << ":"
                << (insertionOrder.nOfThreads - 1) * (getObjectRepositorySize() - 1) << std::endl;
    }
  } else {
    insertionOrder.reserve(getObjectRepositorySize());
    for (NGT::ObjectID id = getObjectRepositorySize() - 1; id != 0; id--) {
      insertionOrder.push_back(length[id].second);
    }
  }
}

void GraphIndex::extractInsertionOrder(InsertionOrder &insertionOrder) {
#ifndef NGT_SHARED_MEMORY_ALLOCATOR
  if (getNumberOfObjects() == 0) {
    NGTThrowException("extractInsertionOrder: No objects.");
  }
  auto edgeSizeBackup                             = NeighborhoodGraph::property.edgeSizeForCreation;
  NeighborhoodGraph::property.edgeSizeForCreation = 10;
  auto nOfThreads =
      insertionOrder.nOfThreads == 0 ? std::thread::hardware_concurrency() : insertionOrder.nOfThreads;

  try {
    InsertionOrder io;
    GraphIndex::createIndexWithInsertionOrder(io, nOfThreads);
  } catch (Exception &err) {
    NeighborhoodGraph::property.edgeSizeForCreation = edgeSizeBackup;
    throw err;
  }
  NeighborhoodGraph::property.edgeSizeForCreation = edgeSizeBackup;

  extractSparseness(insertionOrder);

  NeighborhoodGraph::repository.initialize();
#endif
}

void GraphIndex::createIndexWithSingleThread() {
  GraphRepository &anngRepo = repository;
  ObjectRepository &fr      = objectSpace->getRepository();
  size_t pathAdjustCount    = property.pathAdjustmentInterval;
  NGT::ObjectID id          = 1;
  size_t count              = 0;
  BuildTimeController buildTimeController(*this, NeighborhoodGraph::property);
  for (; id < fr.size(); id++) {
    if (id < anngRepo.size() && anngRepo[id] != 0) {
      continue;
    }
    insert(id);
    buildTimeController.adjustEdgeSize(++count);
    if (pathAdjustCount > 0 && pathAdjustCount <= id) {
      GraphReconstructor::adjustPathsEffectively(static_cast<GraphIndex &>(*this));
      pathAdjustCount += property.pathAdjustmentInterval;
    }
  }
}

static size_t searchMultipleQueryForCreation(GraphIndex &neighborhoodGraph, NGT::ObjectID &id,
                                             CreateIndexJob &job, CreateIndexThreadPool &threads,
                                             size_t sizeOfRepository, Index::InsertionOrder &insertionOrder) {
  ObjectRepository &repo    = neighborhoodGraph.objectSpace->getRepository();
  GraphRepository &anngRepo = neighborhoodGraph.repository;
  size_t cnt                = 0;
  for (; id < repo.size(); id++) {
    if (sizeOfRepository > 0 && id >= sizeOfRepository) {
      break;
    }
    auto oid = insertionOrder.empty() ? id : insertionOrder.getID(id);
    if (repo[oid] == 0) {
      continue;
    }
    if (neighborhoodGraph.NeighborhoodGraph::property.graphType != NeighborhoodGraph::GraphTypeBKNNG) {
      if (oid < anngRepo.size() && anngRepo[oid] != 0) {
        continue;
      }
    }
    job.id = oid;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    job.object = neighborhoodGraph.objectSpace->allocateObject(*repo[oid]);
#else
    job.object = repo[oid];
#endif
    job.batchIdx = cnt;
    threads.pushInputQueue(job);
    cnt++;
    if (cnt >= (size_t)neighborhoodGraph.NeighborhoodGraph::property.batchSizeForCreation) {
      id++;
      break;
    }
  } // for
  return cnt;
}

static void insertMultipleSearchResults(GraphIndex &neighborhoodGraph,
                                        CreateIndexThreadPool::OutputJobQueue &output, ObjectID id,
                                        size_t dataSize) {
  if (neighborhoodGraph.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeRANNG ||
      neighborhoodGraph.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeRIANNG) {
#pragma omp parallel for
    for (size_t i = 0; i < dataSize; i++) {
      neighborhoodGraph.deleteShortcutEdges(*output[i].results);
    }
  }
  // compute distances among all of the resultant objects
  if (neighborhoodGraph.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeANNG ||
      neighborhoodGraph.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeIANNG ||
      neighborhoodGraph.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeONNG ||
      neighborhoodGraph.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeDNNG ||
      neighborhoodGraph.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeRANNG ||
      neighborhoodGraph.NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeRIANNG) {
    // This processing occupies about 30% of total indexing time when batch size is 200.
    // Only initial batch objects should be connected for each other.
    // The number of nodes in the graph is checked to know whether the batch is initial.
    //size_t size = NeighborhoodGraph::property.edgeSizeForCreation;
    size_t size = neighborhoodGraph.NeighborhoodGraph::property.edgeSizeForCreation;
    // add distances from a current object to subsequence objects to imitate of sequential insertion.

    sort(output.begin(), output.end()); // sort by batchIdx

    for (size_t idxi = 0; idxi < dataSize; idxi++) {
      // add distances
      ObjectDistances &objs = *output[idxi].results;
      for (size_t idxj = 0; idxj < idxi; idxj++) {
        ObjectDistance r;
        r.distance =
            neighborhoodGraph.objectSpace->getComparator()(*output[idxi].object, *output[idxj].object);
        r.id = output[idxj].id;
        objs.push_back(r);
      }
      // sort and cut excess edges
      std::sort(objs.begin(), objs.end());
      if (objs.size() > size) {
        objs.resize(size);
      }
    } // for (size_t idxi ....
  } // if (neighborhoodGraph.graphType == NeighborhoodGraph::GraphTypeUDNNG)
  // insert resultant objects into the graph as edges
  for (size_t i = 0; i < dataSize; i++) {
    CreateIndexJob &gr = output[i];
    if ((*gr.results).size() == 0) {
    }
    auto targetID = id == 0 ? gr.id : (id + i);
    if (static_cast<int>(targetID) > neighborhoodGraph.NeighborhoodGraph::property.edgeSizeForCreation &&
        static_cast<int>(gr.results->size()) <
            neighborhoodGraph.NeighborhoodGraph::property.edgeSizeForCreation) {
      if (neighborhoodGraph.NeighborhoodGraph::property.graphType != NeighborhoodGraph::GraphTypeRANNG &&
          neighborhoodGraph.NeighborhoodGraph::property.graphType != NeighborhoodGraph::GraphTypeRIANNG) {
        cerr << "createIndex: Warning. The specified number of edges could not be acquired, because the "
                "pruned parameter [-S] might be set."
             << endl;
        cerr << "  The node id=" << gr.id << ":" << id + i << ":" << targetID << endl;
        cerr << "  The number of edges for creation="
             << neighborhoodGraph.NeighborhoodGraph::property.edgeSizeForCreation << endl;
        cerr << "  The number of edges for the node=" << gr.results->size() << endl;
        cerr << "  The pruned parameter (edgeSizeForSearch [-S])="
             << neighborhoodGraph.NeighborhoodGraph::property.edgeSizeForSearch << endl;
      }
    }
    try {
      neighborhoodGraph.insertNode(gr.id, *gr.results);
    } catch (NGT::Exception &err) {
      std::stringstream msg;
      msg << " Cannot insert the node. " << gr.id << ". " << err.what();
      NGTThrowException(msg);
    }
  }
}

void GraphIndex::createIndexWithInsertionOrder(InsertionOrder &insertionOrder, size_t threadPoolSize,
                                               size_t sizeOfRepository) {
  if (NeighborhoodGraph::property.edgeSizeForCreation == 0) {
    return;
  }
  if (!insertionOrder.empty()) {
    if (objectSpace->getRepository().size() - 1 != insertionOrder.size()) {
      stringstream msg;
      msg << "Index::createIndex: The insertion order size is invalid. "
          << (objectSpace->getRepository().size() - 1) << ":" << insertionOrder.size();
      NGTThrowException(msg);
    }
  }
  threadPoolSize = threadPoolSize == 0 ? std::thread::hardware_concurrency() : threadPoolSize;
  threadPoolSize = threadPoolSize == 0 ? 8 : threadPoolSize;
  if (threadPoolSize <= 1) {
    createIndexWithSingleThread();
  } else {
    Timer timer;
    size_t timerInterval = 100000;
    size_t timerCount    = timerInterval;
    size_t count         = 0;
    timer.start();

    size_t pathAdjustCount = property.pathAdjustmentInterval;
    CreateIndexThreadPool threads(threadPoolSize);
    CreateIndexSharedData sd(*this);

    threads.setSharedData(&sd);
    threads.create();
    CreateIndexThreadPool::OutputJobQueue &output = threads.getOutputJobQueue();

    BuildTimeController buildTimeController(*this, NeighborhoodGraph::property);

    try {
      CreateIndexJob job;
      NGT::ObjectID id = 1;
      for (;;) {
        // search for the nearest neighbors
        size_t cnt =
            searchMultipleQueryForCreation(*this, id, job, threads, sizeOfRepository, insertionOrder);
        if (cnt == 0) {
          break;
        }
        // wait for the completion of the search
        threads.waitForFinish();
        if (output.size() != cnt) {
          cerr << "NNTGIndex::insertGraphIndexByThread: Warning!! Thread response size is wrong." << endl;
          cnt = output.size();
        }
        // insertion
        insertMultipleSearchResults(*this, output, insertionOrder.empty() ? 0 : (id - cnt), cnt);

        while (!output.empty()) {
          delete output.front().results;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
          GraphIndex::objectSpace->deleteObject(output.front().object);
#endif
          output.pop_front();
        }

        count += cnt;
        if (timerCount <= count) {
          timer.stop();
          cerr << "Processed " << timerCount << " time= " << timer << endl;
          timerCount += timerInterval;
          timer.start();
        }
        buildTimeController.adjustEdgeSize(count);
        if (pathAdjustCount > 0 && pathAdjustCount <= count) {
          GraphReconstructor::adjustPathsEffectively(static_cast<GraphIndex &>(*this));
          pathAdjustCount += property.pathAdjustmentInterval;
        }
      }
    } catch (Exception &err) {
      threads.terminate();
      throw err;
    }
    threads.terminate();
  }
}

void GraphIndex::setupPrefetch(NGT::Property &prop) {
  assert(GraphIndex::objectSpace != 0);
  prop.prefetchOffset = GraphIndex::objectSpace->setPrefetchOffset(prop.prefetchOffset);
  prop.prefetchSize   = GraphIndex::objectSpace->setPrefetchSize(prop.prefetchSize);
}

bool NGT::GraphIndex::showStatisticsOfGraph(NGT::GraphIndex &outGraph, char mode, size_t edgeSize) {
  long double distance             = 0.0;
  size_t numberOfNodes             = 0;
  size_t numberOfOutdegree         = 0;
  size_t numberOfNodesWithoutEdges = 0;
  size_t maxNumberOfOutdegree      = 0;
  size_t minNumberOfOutdegree      = SIZE_MAX;
  std::vector<int64_t> indegreeCount;
  std::vector<size_t> outdegreeHistogram;
  std::vector<size_t> indegreeHistogram;
  std::vector<std::vector<float>> indegree;
  NGT::GraphRepository &graph = outGraph.repository;
  NGT::ObjectRepository &repo = outGraph.objectSpace->getRepository();
#ifdef NGT_REFINEMENT
  auto &rrepo = outGraph.refinementObjectSpace->getRepository();
#endif
  indegreeCount.resize(graph.size(), 0);
  indegree.resize(graph.size());
  size_t removedObjectCount = 0;
  bool valid                = true;
  auto &comparator          = outGraph.objectSpace->getComparator();
  size_t nOfEdges           = 0;
  size_t nOfDifferentEdges  = 0;
  for (size_t id = 1; id < graph.size(); id++) {
    if (repo[id] == 0) {
      removedObjectCount++;
      continue;
    }
    NGT::GraphNode *node = 0;
    try {
      node = outGraph.getNode(id);
    } catch (NGT::Exception &err) {
      std::cerr << "ngt info: Error. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
      valid = false;
      continue;
    }
    numberOfNodes++;
    if (numberOfNodes % 1000000 == 0) {
      std::cerr << "Processed " << numberOfNodes << std::endl;
    }
    size_t esize = node->size() > edgeSize ? edgeSize : node->size();
    if (esize == 0) {
      numberOfNodesWithoutEdges++;
    }
    if (esize > maxNumberOfOutdegree) {
      maxNumberOfOutdegree = esize;
    }
    if (esize < minNumberOfOutdegree) {
      minNumberOfOutdegree = esize;
    }
    if (outdegreeHistogram.size() <= esize) {
      outdegreeHistogram.resize(esize + 1);
    }
    outdegreeHistogram[esize]++;
    if (mode == 'e') {
      std::cout << id << "," << esize << ": ";
    }
    NGT::PersistentObject *obj = 0;
    if (mode == 'd' || mode == 'D') {
      if (!repo.isEmpty(id)) {
        obj = repo.get(id);
      }
    }
    for (size_t i = 0; i < esize; i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      NGT::ObjectDistance &n = (*node).at(i, graph.allocator);
#else
      NGT::ObjectDistance &n = (*node)[i];
#endif
      if (std::isnan(n.distance)) {
        stringstream msg;
        msg << "Index::showStatisticsOfGraph: Fatal inner error! The graph has a node with nan distance. "
            << id << ":" << n.id << ":" << n.distance;
        NGTThrowException(msg);
      }
      if (n.id == 0) {
        std::cerr << "ngt info: Warning. id is zero." << std::endl;
        valid = false;
        continue;
      }
      if (mode == 'd' || mode == 'D') {
        if (!repo.isEmpty(id) && !repo.isEmpty(n.id)) {
          nOfEdges++;
          float d = comparator(*obj, *repo.get(n.id));
          if (d != n.distance) {
            nOfDifferentEdges++;
            if (mode == 'D') {
              std::cerr << "The current edge length is different from the indexed length. "
                        << std::setprecision(15) << d << ":" << n.distance << std::setprecision(6) << " "
                        << nOfDifferentEdges << "/" << nOfEdges << std::endl;
            }
          }
        }
      }
      indegreeCount[n.id]++;
      indegree[n.id].push_back(n.distance);
      numberOfOutdegree++;
      double d = n.distance;
      if (mode == 'e') {
        std::cout << n.id << ":" << d << " ";
      }
      distance += d;
    }
    if (mode == 'e') {
      std::cout << std::endl;
    }
  }

  if (mode == 'a') {
    size_t count = 0;
    for (size_t id = 1; id < graph.size(); id++) {
      if (repo[id] == 0) {
        continue;
      }
      NGT::GraphNode *n = 0;
      try {
        n = outGraph.getNode(id);
      } catch (NGT::Exception &err) {
        continue;
      }
      NGT::GraphNode &node = *n;
      for (size_t i = 0; i < node.size(); i++) {
        NGT::GraphNode *nn = 0;
        try {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
          nn = outGraph.getNode(node.at(i, graph.allocator).id);
#else
          nn = outGraph.getNode(node[i].id);
#endif
        } catch (NGT::Exception &err) {
          count++;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
          std::cerr << "Directed edge! " << id << "->" << node.at(i, graph.allocator).id << " no object. "
                    << node.at(i, graph.allocator).id << std::endl;
#else
          std::cerr << "Directed edge! " << id << "->" << node[i].id << " no object. " << node[i].id
                    << std::endl;
#endif
          continue;
        }
        NGT::GraphNode &nnode = *nn;
        bool found            = false;
        for (size_t i = 0; i < nnode.size(); i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
          if (nnode.at(i, graph.allocator).id == id) {
#else
          if (nnode[i].id == id) {
#endif
            found = true;
            break;
          }
        }
        if (!found) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
          std::cerr << "Directed edge! " << id << "->" << node.at(i, graph.allocator).id << " no edge. "
                    << node.at(i, graph.allocator).id << "->" << id << std::endl;
#else
          std::cerr << "Directed edge! " << id << "->" << node[i].id << " no edge. " << node[i].id << "->"
                    << id << std::endl;
#endif
          count++;
        }
      }
    }
    std::cerr << "The number of directed edges=" << count << std::endl;
  }

  // calculate outdegree distance 10
  size_t d10count        = 0;
  long double distance10 = 0.0;
  size_t d10SkipCount    = 0;
  const size_t dcsize    = 10;
  for (size_t id = 1; id < graph.size(); id++) {
    if (repo[id] == 0) {
      continue;
    }
    if (graph.isEmpty(id)) {
      continue;
    }
    NGT::GraphNode *n = 0;
    try {
      n = outGraph.getNode(id);
    } catch (NGT::Exception &err) {
      std::cerr << "ngt info: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
      continue;
    }
    NGT::GraphNode &node = *n;
    if (node.size() < dcsize) {
      d10SkipCount++;
      continue;
    }
    for (size_t i = 0; i < dcsize; i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      distance10 += node.at(i, graph.allocator).distance;
#else
      distance10 += node[i].distance;
#endif
      d10count++;
    }
  }
  if (d10count != 0) {
    distance10 /= (long double)d10count;
  }

  // calculate indegree distance 10
  size_t ind10count              = 0;
  long double indegreeDistance10 = 0.0;
  size_t ind10SkipCount          = 0;
  for (size_t id = 1; id < indegree.size(); id++) {
    std::vector<float> &node = indegree[id];
    if (node.size() < dcsize) {
      ind10SkipCount++;
      continue;
    }
    std::sort(node.begin(), node.end());
    for (size_t i = 0; i < dcsize; i++) {
      if (i > 0 && node[i - 1] > node[i]) {
        stringstream msg;
        msg << "Index::showStatisticsOfGraph: Fatal inner error! Wrong distance order " << node[i - 1] << ":"
            << node[i];
        NGTThrowException(msg);
      }
      indegreeDistance10 += node[i];
      ind10count++;
    }
  }
  if (ind10count != 0) {
    indegreeDistance10 /= (long double)ind10count;
  }

  // calculate variance
  double averageNumberOfOutdegree = (double)numberOfOutdegree / (double)numberOfNodes;
  double sumOfSquareOfOutdegree   = 0;
  double sumOfSquareOfIndegree    = 0;
  for (size_t id = 1; id < graph.size(); id++) {
    if (repo[id] == 0) {
      continue;
    }
    NGT::GraphNode *node = 0;
    try {
      node = outGraph.getNode(id);
    } catch (NGT::Exception &err) {
      std::cerr << "ngt info: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
      continue;
    }
    size_t esize = node->size();
    sumOfSquareOfOutdegree +=
        ((double)esize - averageNumberOfOutdegree) * ((double)esize - averageNumberOfOutdegree);
    sumOfSquareOfIndegree += ((double)indegreeCount[id] - averageNumberOfOutdegree) *
                             ((double)indegreeCount[id] - averageNumberOfOutdegree);
  }

  size_t numberOfNodesWithoutIndegree = 0;
  size_t maxNumberOfIndegree          = 0;
  size_t minNumberOfIndegree          = INT64_MAX;
  for (size_t id = 1; id < graph.size(); id++) {
    if (graph[id] == 0) {
      continue;
    }
    if (indegreeCount[id] == 0) {
      numberOfNodesWithoutIndegree++;
      std::cerr << "Warning! The node without incoming edges. " << id << std::endl;
      valid = false;
    }
    if (indegreeCount[id] > static_cast<int>(maxNumberOfIndegree)) {
      maxNumberOfIndegree = indegreeCount[id];
    }
    if (indegreeCount[id] < static_cast<int64_t>(minNumberOfIndegree)) {
      minNumberOfIndegree = indegreeCount[id];
    }
    if (static_cast<int>(indegreeHistogram.size()) <= indegreeCount[id]) {
      indegreeHistogram.resize(indegreeCount[id] + 1);
    }
    indegreeHistogram[indegreeCount[id]]++;
  }

  size_t count         = 0;
  int medianOutdegree  = -1;
  size_t modeOutdegree = 0;
  size_t max           = 0;
  double c95           = 0.0;
  double c99           = 0.0;
  for (size_t i = 0; i < outdegreeHistogram.size(); i++) {
    count += outdegreeHistogram[i];
    if (medianOutdegree == -1 && count >= numberOfNodes / 2) {
      medianOutdegree = i;
    }
    if (max < outdegreeHistogram[i]) {
      max           = outdegreeHistogram[i];
      modeOutdegree = i;
    }
    if (count > numberOfNodes * 0.95) {
      if (c95 == 0.0) {
        c95 += i * (count - numberOfNodes * 0.95);
      } else {
        c95 += i * outdegreeHistogram[i];
      }
    }
    if (count > numberOfNodes * 0.99) {
      if (c99 == 0.0) {
        c99 += i * (count - numberOfNodes * 0.99);
      } else {
        c99 += i * outdegreeHistogram[i];
      }
    }
  }
  c95 /= (double)numberOfNodes * 0.05;
  c99 /= (double)numberOfNodes * 0.01;

  count               = 0;
  int medianIndegree  = -1;
  size_t modeIndegree = 0;
  max                 = 0;
  double c5           = 0.0;
  double c1           = 0.0;
  for (size_t i = 0; i < indegreeHistogram.size(); i++) {
    if (count < numberOfNodes * 0.05) {
      if (count + indegreeHistogram[i] >= numberOfNodes * 0.05) {
        c5 += i * (numberOfNodes * 0.05 - count);
      } else {
        c5 += i * indegreeHistogram[i];
      }
    }
    if (count < numberOfNodes * 0.01) {
      if (count + indegreeHistogram[i] >= numberOfNodes * 0.01) {
        c1 += i * (numberOfNodes * 0.01 - count);
      } else {
        c1 += i * indegreeHistogram[i];
      }
    }
    count += indegreeHistogram[i];
    if (medianIndegree == -1 && count >= numberOfNodes / 2) {
      medianIndegree = i;
    }
    if (max < indegreeHistogram[i]) {
      max          = indegreeHistogram[i];
      modeIndegree = i;
    }
  }
  c5 /= (double)numberOfNodes * 0.05;
  c1 /= (double)numberOfNodes * 0.01;

  std::cerr << "The number of the objects:\t" << outGraph.getNumberOfObjects() << std::endl;
  std::cerr << "The number of the indexed objects:\t" << outGraph.getNumberOfIndexedObjects() << std::endl;
  std::cerr << "The size of the object repository (not the number of the objects):\t"
            << (repo.size() == 0 ? 0 : repo.size() - 1) << std::endl;
#ifdef NGT_REFINEMENT
  std::cerr << "The size of the refinement object repository (not the number of the objects):\t"
            << (rrepo.size() == 0 ? 0 : rrepo.size() - 1) << std::endl;
#endif
  std::cerr << "The number of the removed objects:\t" << removedObjectCount << "/"
            << (repo.size() == 0 ? 0 : repo.size() - 1) << std::endl;
  std::cerr << "The number of the nodes:\t" << numberOfNodes << std::endl;
  std::cerr << "The number of the edges:\t" << numberOfOutdegree << std::endl;
  std::cerr << "The mean of the edge lengths:\t" << std::setprecision(10)
            << (numberOfOutdegree != 0.0 ? distance / (double)numberOfOutdegree : 0) << std::endl;
  std::cerr << "The mean of the number of the edges per node:\t"
            << (numberOfNodes != 0.0 ? (double)numberOfOutdegree / (double)numberOfNodes : 0) << std::endl;
  std::cerr << "The number of the nodes without edges:\t" << numberOfNodesWithoutEdges << std::endl;
  std::cerr << "The maximum of the outdegrees:\t" << maxNumberOfOutdegree << std::endl;
  if (minNumberOfOutdegree == SIZE_MAX) {
    std::cerr << "The minimum of the outdegrees:\t-NA-" << std::endl;
  } else {
    std::cerr << "The minimum of the outdegrees:\t" << minNumberOfOutdegree << std::endl;
  }
  std::cerr << "The number of the nodes where indegree is 0:\t" << numberOfNodesWithoutIndegree << std::endl;
  std::cerr << "The maximum of the indegrees:\t" << maxNumberOfIndegree << std::endl;
  if (minNumberOfIndegree == INT64_MAX) {
    std::cerr << "The minimum of the indegrees:\t-NA-" << std::endl;
  } else {
    std::cerr << "The minimum of the indegrees:\t" << minNumberOfIndegree << std::endl;
  }
  std::cerr << "The mean of the edge lengths for 10 edges:\t" << std::setprecision(10) << distance10 << "/"
            << d10count << std::endl;
  if (mode == 'd' || mode == 'D') {
    std::cerr << "The number of the different length edges:\t" << nOfDifferentEdges << "/" << nOfEdges
              << std::endl;
  }
  std::cerr
      << "#-nodes,#-edges,#-no-indegree,avg-edges,avg-dist,max-out,min-out,v-out,max-in,min-in,v-in,med-out,"
         "med-in,mode-out,mode-in,c95,c5,o-distance(10),o-skip,i-distance(10),i-skip:"
      << numberOfNodes << ":" << numberOfOutdegree << ":" << numberOfNodesWithoutIndegree << ":"
      << std::setprecision(10) << (double)numberOfOutdegree / (double)numberOfNodes << ":"
      << distance / (double)numberOfOutdegree << ":" << maxNumberOfOutdegree << ":" << minNumberOfOutdegree
      << ":" << sumOfSquareOfOutdegree / (double)numberOfOutdegree << ":" << maxNumberOfIndegree << ":"
      << minNumberOfIndegree << ":" << sumOfSquareOfIndegree / (double)numberOfOutdegree << ":"
      << medianOutdegree << ":" << medianIndegree << ":" << modeOutdegree << ":" << modeIndegree << ":" << c95
      << ":" << c5 << ":" << c99 << ":" << c1 << ":" << distance10 << ":" << d10SkipCount << ":"
      << indegreeDistance10 << ":" << ind10SkipCount << std::endl;
  if (mode == 'h') {
    std::cerr << "#\tout\tin" << std::endl;
    for (size_t i = 0; i < outdegreeHistogram.size() || i < indegreeHistogram.size(); i++) {
      size_t out = outdegreeHistogram.size() <= i ? 0 : outdegreeHistogram[i];
      size_t in  = indegreeHistogram.size() <= i ? 0 : indegreeHistogram[i];
      std::cerr << i << "\t" << out << "\t" << in << std::endl;
    }
  } else if (mode == 'p') {
    std::cerr << "ID\toutdegree\tindegree" << std::endl;
    for (size_t id = 1; id < graph.size(); id++) {
      std::cerr << id << "\t" << outGraph.getNode(id)->size() << "\t" << indegreeCount[id] << std::endl;
    }
  }
  return valid;
}

NGT::GraphIndex::GraphStatistics NGT::GraphIndex::getGraphStatistics(NGT::GraphIndex &outGraph, char mode,
                                                                     size_t edgeSize) {
  NGT::GraphIndex::GraphStatistics stats = {};
  long double distance                   = 0.0;
  size_t numberOfNodes                   = 0;
  size_t numberOfOutdegree               = 0;
  size_t numberOfNodesWithoutEdges       = 0;
  size_t maxNumberOfOutdegree            = 0;
  size_t minNumberOfOutdegree            = SIZE_MAX;
  std::vector<int64_t> indegreeCount;
  std::vector<size_t> outdegreeHistogram;
  std::vector<size_t> indegreeHistogram;
  std::vector<std::vector<float>> indegree;
  NGT::GraphRepository &graph = outGraph.repository;
  NGT::ObjectRepository &repo = outGraph.objectSpace->getRepository();

#ifdef NGT_REFINEMENT
  auto &rrepo = outGraph.refinementObjectSpace->getRepository();
#endif

  indegreeCount.resize(graph.size(), 0);
  indegree.resize(graph.size());
  size_t removedObjectCount = 0;
  bool valid                = true;

  for (size_t id = 1; id < graph.size(); id++) {
    if (repo[id] == 0) {
      removedObjectCount++;
      continue;
    }
    NGT::GraphNode *node = nullptr;
    try {
      node = outGraph.getNode(id);
    } catch (NGT::Exception &err) {
      std::cerr << "ngt info: Error. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
      valid = false;
      continue;
    }
    if (node == nullptr) continue;
    numberOfNodes++;
    size_t esize = std::min(node->size(), edgeSize); // edge size limitation by using edgeSize argument.
    if (esize == 0) {
      numberOfNodesWithoutEdges++;
    }
    maxNumberOfOutdegree = std::max(maxNumberOfOutdegree, esize);
    minNumberOfOutdegree = std::min(minNumberOfOutdegree, esize);
    if (outdegreeHistogram.size() <= esize) {
      outdegreeHistogram.resize(esize + 1);
    }
    outdegreeHistogram[esize]++;
    for (size_t i = 0; i < esize; i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      NGT::ObjectDistance &n = (*node).at(i, graph.allocator);
#else
      NGT::ObjectDistance &n = (*node)[i];
#endif
      if (std::isnan(n.distance)) {
        std::stringstream msg;
        msg << "NGT::GraphIndex::getGraphStatistics: Fatal inner error! The graph has a node with nan "
               "distance. "
            << id << ":" << n.id << ":" << n.distance;
        NGTThrowException(msg);
      }
      if (n.id == 0) {
        std::cerr << "ngt info: Warning. id is zero." << std::endl;
        valid = false;
      }
      indegreeCount[n.id]++;
      indegree[n.id].push_back(n.distance);
      numberOfOutdegree++;
      distance += n.distance;
    }
  }

  // if mode is 'a' process additional edge checking
  if (mode == 'a') {
    size_t count = 0;
    for (size_t id = 1; id < graph.size(); id++) {
      if (repo[id] == 0) continue;
      NGT::GraphNode *n = nullptr;
      try {
        n = outGraph.getNode(id);
      } catch (NGT::Exception &err) {
        continue;
      }
      if (n == nullptr) continue;
      NGT::GraphNode &node = *n;
      for (size_t i = 0; i < node.size(); i++) {
        NGT::GraphNode *nn = nullptr;
        try {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
          nn = outGraph.getNode(node.at(i, graph.allocator).id);
#else
          nn = outGraph.getNode(node[i].id);
#endif
        } catch (NGT::Exception &err) {
          count++;
          continue;
        }
        NGT::GraphNode &nnode = *nn;
        bool found            = false;
        for (size_t i = 0; i < nnode.size(); i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
          if (nnode.at(i, graph.allocator).id == id) {
#else
          if (nnode[i].id == id) {
#endif
            found = true;
            break;
          }
        }
        if (!found) count++;
      }
    }
    std::cerr << "The number of directed edges=" << count << std::endl;
  }

  // Calculate outdegree distance for the first 10 edges
  size_t d10count        = 0;
  long double distance10 = 0.0;
  size_t d10SkipCount    = 0;
  const size_t dcsize    = 10;
  for (size_t id = 1; id < graph.size(); id++) {
    if (repo[id] == 0) continue;
    NGT::GraphNode *n = nullptr;
    try {
      n = outGraph.getNode(id);
    } catch (NGT::Exception &err) {
      std::cerr << "ngt info: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
      continue;
    }
    if (n == nullptr) continue;
    NGT::GraphNode &node = *n;
    if (node.size() < dcsize) {
      d10SkipCount++;
      continue;
    }
    for (size_t i = 0; i < dcsize; i++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
      distance10 += node.at(i, graph.allocator).distance;
#else
      distance10 += node[i].distance;
#endif
      d10count++;
    }
  }
  if (d10count != 0) {
    distance10 /= static_cast<long double>(d10count);
  }

  // Calculate indegree distance for the first 10 edges
  size_t ind10count              = 0;
  long double indegreeDistance10 = 0.0;
  size_t ind10SkipCount          = 0;
  for (size_t id = 1; id < indegree.size(); id++) {
    std::vector<float> &node = indegree[id];
    if (node.size() < dcsize) {
      ind10SkipCount++;
      continue;
    }
    std::sort(node.begin(), node.end());
    for (size_t i = 0; i < dcsize; i++) {
      if (i > 0 && node[i - 1] > node[i]) {
        stringstream msg;
        msg << "NGT::GraphIndex::getGraphStatistics: Fatal inner error! Wrong distance order " << node[i - 1]
            << ":" << node[i];
        NGTThrowException(msg);
      }
      indegreeDistance10 += node[i];
      ind10count++;
    }
  }
  if (ind10count != 0) {
    indegreeDistance10 /= static_cast<long double>(ind10count);
  }

  // Calculate variance
  double averageNumberOfOutdegree = static_cast<double>(numberOfOutdegree) / numberOfNodes;
  double sumOfSquareOfOutdegree   = 0;
  double sumOfSquareOfIndegree    = 0;
  for (size_t id = 1; id < graph.size(); id++) {
    if (repo[id] == 0) continue;
    NGT::GraphNode *node = nullptr;
    try {
      node = outGraph.getNode(id);
    } catch (NGT::Exception &err) {
      std::cerr << "ngt info: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
      continue;
    }
    size_t esize = node->size();
    sumOfSquareOfOutdegree += (static_cast<double>(esize) - averageNumberOfOutdegree) *
                              (static_cast<double>(esize) - averageNumberOfOutdegree);
    sumOfSquareOfIndegree += (static_cast<double>(indegreeCount[id]) - averageNumberOfOutdegree) *
                             (static_cast<double>(indegreeCount[id]) - averageNumberOfOutdegree);
  }

  size_t numberOfNodesWithoutIndegree = 0;
  size_t maxNumberOfIndegree          = 0;
  size_t minNumberOfIndegree          = INT64_MAX;
  for (size_t id = 1; id < graph.size(); id++) {
    if (graph[id] == 0) continue;
    if (indegreeCount[id] == 0) {
      numberOfNodesWithoutIndegree++;
      std::cerr << "Warning! The node without incoming edges. " << id << std::endl;
      valid = false;
    }
    maxNumberOfIndegree = std::max(maxNumberOfIndegree, static_cast<size_t>(indegreeCount[id]));
    minNumberOfIndegree = std::min(minNumberOfIndegree, static_cast<size_t>(indegreeCount[id]));
    if (indegreeHistogram.size() <= static_cast<size_t>(indegreeCount[id])) {
      indegreeHistogram.resize(indegreeCount[id] + 1);
    }
    indegreeHistogram[indegreeCount[id]]++;
  }

  size_t count         = 0;
  int medianOutdegree  = -1;
  size_t modeOutdegree = 0;
  size_t max           = 0;
  double c95           = 0.0;
  double c99           = 0.0;
  for (size_t i = 0; i < outdegreeHistogram.size(); i++) {
    count += outdegreeHistogram[i];
    if (medianOutdegree == -1 && count >= numberOfNodes / 2) {
      medianOutdegree = i;
    }
    if (max < outdegreeHistogram[i]) {
      max           = outdegreeHistogram[i];
      modeOutdegree = i;
    }
    if (count > numberOfNodes * 0.95) {
      if (c95 == 0.0) {
        c95 += i * (count - numberOfNodes * 0.95);
      } else {
        c95 += i * outdegreeHistogram[i];
      }
    }
    if (count > numberOfNodes * 0.99) {
      if (c99 == 0.0) {
        c99 += i * (count - numberOfNodes * 0.99);
      } else {
        c99 += i * outdegreeHistogram[i];
      }
    }
  }
  c95 /= static_cast<double>(numberOfNodes) * 0.05;
  c99 /= static_cast<double>(numberOfNodes) * 0.01;

  count               = 0;
  int medianIndegree  = -1;
  size_t modeIndegree = 0;
  max                 = 0;
  double c5           = 0.0;
  double c1           = 0.0;
  for (size_t i = 0; i < indegreeHistogram.size(); i++) {
    if (count < numberOfNodes * 0.05) {
      if (count + indegreeHistogram[i] >= numberOfNodes * 0.05) {
        c5 += i * (numberOfNodes * 0.05 - count);
      } else {
        c5 += i * indegreeHistogram[i];
      }
    }
    if (count < numberOfNodes * 0.01) {
      if (count + indegreeHistogram[i] >= numberOfNodes * 0.01) {
        c1 += i * (numberOfNodes * 0.01 - count);
      } else {
        c1 += i * indegreeHistogram[i];
      }
    }
    count += indegreeHistogram[i];
    if (medianIndegree == -1 && count >= numberOfNodes / 2) {
      medianIndegree = i;
    }
    if (max < indegreeHistogram[i]) {
      max          = indegreeHistogram[i];
      modeIndegree = i;
    }
  }
  c5 /= static_cast<double>(numberOfNodes) * 0.05;
  c1 /= static_cast<double>(numberOfNodes) * 0.01;

  stats.setNumberOfObjects(outGraph.getNumberOfObjects());
  stats.setNumberOfIndexedObjects(outGraph.getNumberOfIndexedObjects());
  stats.setSizeOfObjectRepository(repo.size() == 0 ? 0 : repo.size() - 1);
#ifdef NGT_REFINEMENT
  stats.setSizeOfRefinementObjectRepository(rrepo.size() == 0 ? 0 : rrepo.size() - 1);
#endif
  stats.setNumberOfRemovedObjects(removedObjectCount);
  stats.setNumberOfNodes(numberOfNodes);
  stats.setNumberOfEdges(numberOfOutdegree);
  stats.setMeanEdgeLength(numberOfOutdegree != 0.0 ? distance / static_cast<double>(numberOfOutdegree) : 0.0);
  stats.setMeanNumberOfEdgesPerNode(numberOfNodes != 0.0 ? static_cast<double>(numberOfOutdegree) /
                                                               static_cast<double>(numberOfNodes)
                                                         : 0.0);
  stats.setNumberOfNodesWithoutEdges(numberOfNodesWithoutEdges);
  stats.setMaxNumberOfOutdegree(maxNumberOfOutdegree);
  stats.setMinNumberOfOutdegree(minNumberOfOutdegree == SIZE_MAX ? static_cast<size_t>(-1)
                                                                 : minNumberOfOutdegree);
  stats.setNumberOfNodesWithoutIndegree(numberOfNodesWithoutIndegree);
  stats.setMaxNumberOfIndegree(maxNumberOfIndegree);
  stats.setMinNumberOfIndegree(minNumberOfIndegree == INT64_MAX ? static_cast<size_t>(-1)
                                                                : minNumberOfIndegree);
  stats.setMeanEdgeLengthFor10Edges(distance10);
  stats.setNodesSkippedFor10Edges(d10SkipCount);
  stats.setMeanIndegreeDistanceFor10Edges(indegreeDistance10);
  stats.setNodesSkippedForIndegreeDistance(ind10SkipCount);
  stats.setVarianceOfOutdegree(sumOfSquareOfOutdegree / static_cast<double>(numberOfOutdegree));
  stats.setVarianceOfIndegree(sumOfSquareOfIndegree / static_cast<double>(numberOfOutdegree));
  stats.setMedianOutdegree(medianOutdegree);
  stats.setModeOutdegree(modeOutdegree);
  stats.setC95Outdegree(c95);
  stats.setC99Outdegree(c99);
  stats.setMedianIndegree(medianIndegree);
  stats.setModeIndegree(modeIndegree);
  stats.setC5Indegree(c5);
  stats.setC1Indegree(c1);
  stats.setIndegreeCount(std::move(indegreeCount));
  stats.setOutdegreeHistogram(std::move(outdegreeHistogram));
  stats.setIndegreeHistogram(std::move(indegreeHistogram));
  stats.setValid(valid);

  return stats;
}

void GraphAndTreeIndex::createIndexWithInsertionOrder(InsertionOrder &insertionOrder, size_t threadPoolSize,
                                                      size_t sizeOfRepository) {
  if (NeighborhoodGraph::property.edgeSizeForCreation == 0) {
    return;
  }
  if (!insertionOrder.empty()) {
    if (GraphIndex::objectSpace->getRepository().size() - 1 != insertionOrder.size()) {
      stringstream msg;
      msg << "Index::createIndex: The insertion order size is invalid. "
          << (GraphIndex::objectSpace->getRepository().size() - 1) << ":" << insertionOrder.size();
      NGTThrowException(msg);
    }
  }
  threadPoolSize = threadPoolSize == 0 ? std::thread::hardware_concurrency() : threadPoolSize;
  threadPoolSize = threadPoolSize == 0 ? 8 : threadPoolSize;
  Timer timer;
  size_t timerInterval = 100000;
  size_t timerCount    = timerInterval;
  size_t count         = 0;
  timer.start();
  size_t pathAdjustCount = property.pathAdjustmentInterval;
  CreateIndexThreadPool threads(threadPoolSize);

  CreateIndexSharedData sd(*this);

  threads.setSharedData(&sd);
  threads.create();
  CreateIndexThreadPool::OutputJobQueue &output = threads.getOutputJobQueue();

  BuildTimeController buildTimeController(*this, NeighborhoodGraph::property);

  try {
    CreateIndexJob job;
    NGT::ObjectID id = 1;
    for (;;) {
      size_t cnt = searchMultipleQueryForCreation(*this, id, job, threads, sizeOfRepository, insertionOrder);
      if (cnt == 0) {
        break;
      }
      threads.waitForFinish();

      if (output.size() != cnt) {
        cerr << "NNTGIndex::insertGraphIndexByThread: Warning!! Thread response size is wrong." << endl;
        cnt = output.size();
      }

      insertMultipleSearchResults(*this, output, insertionOrder.empty() ? 0 : (id - cnt), cnt);

      for (size_t i = 0; i < cnt; i++) {
        CreateIndexJob &job = output[i];
        if (GraphIndex::objectSpace->isNormalizedDistance()) {
          if (job.results->size() > 0) {
            auto *o                    = GraphIndex::getObjectRepository().get((*job.results)[0].id);
            (*job.results)[0].distance = GraphIndex::objectSpace->compareWithL1(*job.object, *o);
          }
        }
        if (((job.results->size() > 0) && ((*job.results)[0].distance != 0.0)) ||
            (job.results->size() == 0)) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
          Object *f = GraphIndex::objectSpace->allocateObject(*job.object);
          DVPTree::InsertContainer tiobj(*f, job.id);
#else
          DVPTree::InsertContainer tiobj(*job.object, job.id);
#endif
          try {
            DVPTree::insert(tiobj);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
            GraphIndex::objectSpace->deleteObject(f);
#endif
          } catch (Exception &err) {
            cerr << "NGT::createIndex: Fatal error. ID=" << job.id << ":";
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
            GraphIndex::objectSpace->deleteObject(f);
#endif
            if (NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeKNNG) {
              cerr << err.what() << " continue.." << endl;
            } else {
              throw err;
            }
          }
        }
      } // for

      while (!output.empty()) {
        delete output.front().results;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
        GraphIndex::objectSpace->deleteObject(output.front().object);
#endif
        output.pop_front();
      }

      count += cnt;
      if (timerCount <= count) {
        timer.stop();
        cerr << "Processed " << timerCount << " objects. time= " << timer
             << " vm size=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr()
             << endl;
        timerCount += timerInterval;
        timer.restart();
      }
      buildTimeController.adjustEdgeSize(count);
      if (pathAdjustCount > 0 && pathAdjustCount <= count) {
        GraphReconstructor::adjustPathsEffectively(static_cast<GraphIndex &>(*this));
        pathAdjustCount += property.pathAdjustmentInterval;
      }
    }
  } catch (Exception &err) {
    threads.terminate();
    throw err;
  }
  threads.terminate();
}

void GraphAndTreeIndex::createIndex(const vector<pair<NGT::Object *, size_t>> &objects,
                                    vector<InsertionResult> &ids, float range, size_t threadPoolSize) {
  Timer timer;
  size_t timerInterval = 100000;
  size_t timerCount    = timerInterval;
  size_t count         = 0;
  timer.start();
  if (threadPoolSize <= 0) {
    cerr << "Not implemented!!" << endl;
    abort();
  } else {
    CreateIndexThreadPool threads(threadPoolSize);
    CreateIndexSharedData sd(*this);
    threads.setSharedData(&sd);
    threads.create();
    CreateIndexThreadPool::OutputJobQueue &output = threads.getOutputJobQueue();
    try {
      CreateIndexJob job;
      size_t idx = 0;
      for (;;) {
        size_t cnt = 0;
        {
          for (; idx < objects.size(); idx++) {
            if (objects[idx].first == 0) {
              ids.push_back(InsertionResult());
              continue;
            }
            job.id       = 0;
            job.results  = 0;
            job.object   = objects[idx].first;
            job.batchIdx = ids.size();
            // insert an empty entry to prepare.
            ids.push_back(InsertionResult(job.id, false, 0.0));
            threads.pushInputQueue(job);
            cnt++;
            if (cnt >= (size_t)NeighborhoodGraph::property.batchSizeForCreation) {
              idx++;
              break;
            }
          }
        }
        if (cnt == 0) {
          break;
        }
        threads.waitForFinish();
        if (output.size() != cnt) {
          cerr << "NNTGIndex::insertGraphIndexByThread: Warning!! Thread response size is wrong." << endl;
          cnt = output.size();
        }
        {
          size_t size = NeighborhoodGraph::property.edgeSizeForCreation;
          sort(output.begin(), output.end());
          for (size_t idxi = 0; idxi < cnt; idxi++) {
            // add distances
            ObjectDistances &objs = *output[idxi].results;
            for (size_t idxj = 0; idxj < idxi; idxj++) {
              if (output[idxi].batchIdx == output[idxj].batchIdx) {
                continue;
              }
              if (output[idxj].id == 0) {
                continue;
              }
              ObjectDistance r;
              r.distance =
                  GraphIndex::objectSpace->getComparator()(*output[idxi].object, *output[idxj].object);
              r.id = output[idxj].id;
              objs.emplace_back(r);
            }
            std::sort(objs.begin(), objs.end());
            if (objs.size() > size) {
              objs.resize(size);
            }
            if ((objs.size() > 0) && (range >= 0.0) && (objs[0].distance <= range)) {
              // The line below was replaced by the line above to consider EPSILON for float comparison. 170702
              // if ((objs.size() > 0) && (range < 0.0 || (objs[0].distance <= range))) {
              // An identical or similar object already exits
              ids[output[idxi].batchIdx].identical = true;
              ids[output[idxi].batchIdx].id        = objs[0].id;
              ids[output[idxi].batchIdx].distance  = objs[0].distance;
              output[idxi].id                      = 0;
            } else {
              assert(output[idxi].id == 0);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
              PersistentObject *obj = GraphIndex::objectSpace->allocatePersistentObject(*output[idxi].object);
              output[idxi].id       = GraphIndex::objectSpace->insert(obj);
#else
              output[idxi].id = GraphIndex::objectSpace->insert(output[idxi].object);
#endif
              ids[output[idxi].batchIdx].id = output[idxi].id;
            }
          }
        }
        // insert resultant objects into the graph as edges
        for (size_t i = 0; i < cnt; i++) {
          CreateIndexJob &job = output.front();
          if (job.id != 0) {
            if (property.indexType == NGT::Property::GraphAndTree) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
              Object *f = GraphIndex::objectSpace->allocateObject(*job.object);
              DVPTree::InsertContainer tiobj(*f, job.id);
#else
              DVPTree::InsertContainer tiobj(*job.object, job.id);
#endif
              try {
                DVPTree::insert(tiobj);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
                GraphIndex::objectSpace->deleteObject(f);
#endif
              } catch (Exception &err) {
                cerr << "NGT::createIndex: Fatal error. ID=" << job.id << ":" << err.what();
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
                GraphIndex::objectSpace->deleteObject(f);
#endif
                if (NeighborhoodGraph::property.graphType == NeighborhoodGraph::GraphTypeKNNG) {
                  cerr << err.what() << " continue.." << endl;
                } else {
                  throw err;
                }
              }
            }
            if (((*job.results).size() == 0) && (job.id != 1)) {
              cerr << "insert warning!! No searched nodes!. If the first time, no problem. " << job.id
                   << endl;
            }
            try {
              GraphIndex::insertNode(job.id, *job.results);
            } catch (NGT::Exception &err) {
              std::stringstream msg;
              msg << " Cannot insert the node. " << job.id << ". " << err.what();
              NGTThrowException(msg);
            }
          }
          if (job.results != 0) {
            delete job.results;
          }
          output.pop_front();
        }

        count += cnt;
        if (timerCount <= count) {
          timer.stop();
          cerr << "Processed " << timerCount << " time= " << timer << endl;
          timerCount += timerInterval;
          timer.start();
        }
      }
    } catch (Exception &err) {
      cerr << "thread terminate!" << endl;
      threads.terminate();
      throw err;
    }
    threads.terminate();
  }
}

static bool findPathAmongIdenticalObjects(GraphAndTreeIndex &graph, size_t srcid, size_t dstid) {
  stack<size_t> nodes;
  unordered_set<size_t> done;
  nodes.push(srcid);
  while (!nodes.empty()) {
    auto tid = nodes.top();
    nodes.pop();
    done.insert(tid);
    GraphNode &node = *graph.GraphIndex::getNode(tid);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    for (auto i = node.begin(graph.GraphIndex::repository.allocator);
         i != node.end(graph.GraphIndex::repository.allocator); ++i) {
#else
    for (auto i = node.begin(); i != node.end(); ++i) {
#endif
      if ((*i).distance != 0.0) {
        break;
      }
      if ((*i).id == dstid) {
        return true;
      }
      if (done.count((*i).id) == 0) {
        nodes.push((*i).id);
      }
    }
  }
  return false;
}

bool GraphAndTreeIndex::verify(vector<uint8_t> &status, bool info, char mode) {
  bool valid = GraphIndex::verify(status, info);
  if (!valid) {
    cerr << "The graph or object is invalid!" << endl;
  }
  bool treeValid = DVPTree::verify(GraphIndex::objectSpace->getRepository().size(), status);
  if (!treeValid) {
    cerr << "The tree is invalid" << endl;
  }
  valid = valid && treeValid;
  // status: tree|graph|object
  cerr << "Started checking consistency..." << endl;
  for (size_t id = 1; id < status.size(); id++) {
    if (id % 100000 == 0) {
      cerr << "The number of processed objects=" << id << endl;
    }
    if (status[id] != 0x00 && status[id] != 0x07) {
      if (status[id] == 0x03) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
        NGT::Object *po = GraphIndex::objectSpace->allocateObject(*GraphIndex::getObjectRepository().get(id));
        NGT::SearchContainer sc(*po);
#else
        NGT::SearchContainer sc(*GraphIndex::getObjectRepository().get(id));
#endif
        NGT::ObjectDistances objects;
        sc.setResults(&objects);
        sc.id                     = 0;
        sc.radius                 = 0.0;
        sc.explorationCoefficient = 1.1;
        sc.edgeSize               = 0;
        ObjectDistances seeds;
        seeds.push_back(ObjectDistance(id, 0.0));
        objects.clear();
        try {
          GraphIndex::search(sc, seeds);
        } catch (Exception &err) {
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
          GraphIndex::objectSpace->deleteObject(po);
#endif
          cerr << "Fatal Error!: Cannot search! " << err.what() << endl;
          objects.clear();
        }
        size_t n                       = 0;
        bool registeredIdenticalObject = false;
        for (; n < objects.size(); n++) {
          if (objects[n].id != id && status[objects[n].id] == 0x07) {
            registeredIdenticalObject = true;
            break;
          }
        }
        if (!registeredIdenticalObject) {
          if (info) {
            cerr << "info: not found the registered same objects. id=" << id << " size=" << objects.size()
                 << endl;
          }
          sc.id                     = 0;
          sc.radius                 = FLT_MAX;
          sc.explorationCoefficient = 1.2;
          sc.edgeSize               = 0;
          sc.size                   = objects.size() < 100 ? 100 : objects.size() * 2;
          ObjectDistances seeds;
          seeds.push_back(ObjectDistance(id, 0.0));
          objects.clear();
          try {
            GraphIndex::search(sc, seeds);
          } catch (Exception &err) {
            cerr << "Fatal Error!: Cannot search! " << err.what() << endl;
            objects.clear();
          }
          registeredIdenticalObject = false;
          for (n = 0; n < objects.size(); n++) {
            if (objects[n].distance != 0.0) break;
            if (objects[n].id != id && status[objects[n].id] == 0x07) {
              registeredIdenticalObject = true;
              if (info) {
                cerr << "info: found by using mode accurate search. " << objects[n].id << endl;
              }
              break;
            }
          }
        }
        if (!registeredIdenticalObject && mode != 's') {
          if (info) {
            cerr << "info: not found by using more accurate search." << endl;
          }
          sc.id                     = 0;
          sc.radius                 = 0.0;
          sc.explorationCoefficient = 1.1;
          sc.edgeSize               = 0;
          sc.size                   = SIZE_MAX;
          objects.clear();
          linearSearch(sc);
          n                         = 0;
          registeredIdenticalObject = false;
          for (; n < objects.size(); n++) {
            if (objects[n].distance != 0.0) break;
            if (objects[n].id != id && status[objects[n].id] == 0x07) {
              registeredIdenticalObject = true;
              if (info) {
                cerr << "info: found by using linear search. " << objects[n].id << endl;
              }
              break;
            }
          }
        }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
        GraphIndex::objectSpace->deleteObject(po);
#endif
        if (registeredIdenticalObject) {
          if (info) {
            cerr << "Info ID=" << id << ":" << static_cast<int>(status[id]) << endl;
            cerr << "  found the valid same objects. " << objects[n].id << endl;
          }
          GraphNode &fromNode = *GraphIndex::getNode(id);
          bool fromFound      = false;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
          for (auto i = fromNode.begin(GraphIndex::repository.allocator);
               i != fromNode.end(GraphIndex::repository.allocator); ++i) {
#else
          for (auto i = fromNode.begin(); i != fromNode.end(); ++i) {
#endif
            if ((*i).id == objects[n].id) {
              fromFound = true;
            }
          }
          GraphNode &toNode = *GraphIndex::getNode(objects[n].id);
          bool toFound      = false;
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
          for (auto i = toNode.begin(GraphIndex::repository.allocator);
               i != toNode.end(GraphIndex::repository.allocator); ++i) {
#else
          for (auto i = toNode.begin(); i != toNode.end(); ++i) {
#endif
            if ((*i).id == id) {
              toFound = true;
            }
          }
          if (!fromFound || !toFound) {
            if (info) {
              if (!fromFound && !toFound) {
                cerr << "Warning no undirected edge between " << id << "(" << fromNode.size() << ") and "
                     << objects[n].id << "(" << toNode.size() << ")." << endl;
              } else if (!fromFound) {
                cerr << "Warning no directed edge from " << id << "(" << fromNode.size() << ") to "
                     << objects[n].id << "(" << toNode.size() << ")." << endl;
              } else if (!toFound) {
                cerr << "Warning no reverse directed edge from " << id << "(" << fromNode.size() << ") to "
                     << objects[n].id << "(" << toNode.size() << ")." << endl;
              }
            }
            if (!findPathAmongIdenticalObjects(*this, id, objects[n].id)) {
              cerr << "Warning no path from " << id << " to " << objects[n].id << endl;
            }
            if (!findPathAmongIdenticalObjects(*this, objects[n].id, id)) {
              cerr << "Warning no reverse path from " << id << " to " << objects[n].id << endl;
            }
          }
        } else {
          if (mode == 's') {
            cerr << "Warning: not found the valid same object, but not try to use linear search." << endl;
            cerr << "Error! ID=" << id << ":" << static_cast<int>(status[id]) << endl;
          } else {
            cerr << "Warning: not found the valid same object even by using linear search." << endl;
            cerr << "Error! ID=" << id << ":" << static_cast<int>(status[id]) << endl;
            valid = false;
          }
        }
      } else if (status[id] == 0x01) {
        if (info) {
          cerr << "Warning! ID=" << id << ":" << static_cast<int>(status[id]) << endl;
          cerr << "  not inserted into the indexes" << endl;
        }
      } else {
        cerr << "Error! ID=" << id << ":" << static_cast<int>(status[id]) << endl;
        valid = false;
      }
    }
  }
  return valid;
}
