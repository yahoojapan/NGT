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

#include "NGT/defines.h"
#include "NGT/Common.h"
#include "NGT/ObjectSpace.h"
#include "NGT/ObjectRepository.h"

NGT::Distance NGT::ObjectSpace::compareWithL1(NGT::Object &o1, NGT::Object &o2) {
  auto dim = getPaddedDimension();
  NGT::Distance d;
  if (getObjectType() == typeid(uint8_t) || getObjectType() == typeid(quint8) ||
      getObjectType() == typeid(qsint8)) {
    d = PrimitiveComparator::compareL1(reinterpret_cast<uint8_t *>(o1.getPointer()),
                                       reinterpret_cast<uint8_t *>(o2.getPointer()), dim);
#ifdef NGT_HALF_FLOAT
  } else if (getObjectType() == typeid(float16)) {
    d = PrimitiveComparator::compareL1(reinterpret_cast<float16 *>(o1.getPointer()),
                                       reinterpret_cast<float16 *>(o2.getPointer()), dim);
#endif
  } else if (getObjectType() == typeid(float)) {
    d = PrimitiveComparator::compareL1(reinterpret_cast<float *>(o1.getPointer()),
                                       reinterpret_cast<float *>(o2.getPointer()), dim);
  } else {
    std::stringstream msg;
    msg << "ObjectSpace::compareWithL1: Fatal Inner Error! Unexpected object type. "
        << getObjectType().name();
    NGTThrowException(msg);
  }
  return d;
}

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
NGT::Distance NGT::ObjectSpace::compareWithL1(NGT::Object &o1, NGT::PersistentObject &o2) {
  auto dim = getPaddedDimension();
  NGT::Distance d;
  if (getObjectType() == typeid(uint8_t)) {
    d = PrimitiveComparator::compareL1(
        reinterpret_cast<uint8_t *>(o1.getPointer()),
        reinterpret_cast<uint8_t *>(o2.getPointer(getRepository().getAllocator())), dim);
#ifdef NGT_HALF_FLOAT
  } else if (getObjectType() == typeid(float16)) {
    d = PrimitiveComparator::compareL1(
        reinterpret_cast<float16 *>(o1.getPointer()),
        reinterpret_cast<float16 *>(o2.getPointer(getRepository().getAllocator())), dim);
#endif
  } else if (getObjectType() == typeid(float)) {
    d = PrimitiveComparator::compareL1(
        reinterpret_cast<float *>(o1.getPointer()),
        reinterpret_cast<float *>(o2.getPointer(getRepository().getAllocator())), dim);
  } else {
    std::stringstream msg;
    msg << "ObjectSpace::compareWithL1: Fatal Inner Error! Unexpected object type.";
    NGTThrowException(msg);
  }
  return d;
}

#endif
