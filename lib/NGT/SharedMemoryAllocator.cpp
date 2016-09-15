
//
// Copyright (C) 2015-2016 Yahoo! JAPAN Research
//
// This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
// To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
//

#include	"NGT/SharedMemoryAllocator.h"



void* operator 
new(size_t size, SharedMemoryAllocator &allocator) 
{
  void *addr = allocator.allocate(size);
#ifdef MEMORY_ALLOCATOR_INFO
  std::cerr << "new:" << size << " " << addr << " " << allocator.getTotalSize() << std::endl;
#endif
  return addr;
}

void* operator 
new[](size_t size, SharedMemoryAllocator &allocator) 
{

  void *addr = allocator.allocate(size);
#ifdef MEMORY_ALLOCATOR_INFO
  std::cerr << "new[]:" << size << " " << addr << " " << allocator.getTotalSize() << std::endl;
#endif
  return addr;
}
