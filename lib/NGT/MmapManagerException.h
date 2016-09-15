
//
// Copyright (C) 2015-2016 Yahoo! JAPAN Research
//
// This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
// To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
//

#pragma once

#include <iostream>
#include <exception>
#include <stdexcept>

namespace MemoryManager{
  class MmapManagerException : public std::domain_error{
  public:
  MmapManagerException(const std::string &msg) : std::domain_error(msg){}
  };
}
