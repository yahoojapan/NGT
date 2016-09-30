
//
// Copyright (C) 2015-2016 Yahoo! JAPAN Research
//
// This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
// To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
//

#pragma once

#include <sys/types.h>
#include <stdint.h>
#include <string>
#include <memory>

#define ABS_ADDR(x, y) (void *)(x + (char *)y);

#define USE_MMAP_MANAGER

namespace MemoryManager{
  
  typedef enum _option_reuse_t{
    REUSE_DATA_CLASSIFY,
    REUSE_DATA_QUEUE,
    REUSE_DATA_QUEUE_PLUS,    
  }option_reuse_t;

  typedef enum _reuse_state_t{
    REUSE_STATE_OK, 
    REUSE_STATE_FALSE, 
    REUSE_STATE_ALLOC, 
  }reuse_state_t;

  typedef enum _check_statistics_t{
    CHECK_STATS_USE_SIZE,  
    CHECK_STATS_USE_NUM,   
    CHECK_STATS_FREE_SIZE, 
    CHECK_STATS_FREE_NUM,  
  }check_statistics_t;

  typedef struct _init_option_st{
    bool use_expand;  
    option_reuse_t reuse_type; 
  }init_option_st;


  class MmapManager{
  public:
    MmapManager();
    ~MmapManager();
        
    bool init(const std::string &filePath, size_t size, const init_option_st *optionst = NULL) const;
    bool openMemory(const std::string &filePath);
    void closeMemory(const bool force = false);
    off_t alloc(const size_t size, const bool not_reuse_flag = false);
    void free(const off_t p);
    off_t reuse(const size_t size, reuse_state_t &reuse_state);
    void *getAbsAddr(off_t p) const;
    off_t getRelAddr(const void *p) const;

    size_t getTotalSize() const; 
    size_t getUseSize() const;   
    uint64_t getUseNum() const;  
    size_t getFreeSize() const;  
    uint64_t getFreeNum() const; 
    uint16_t getUnitNum() const; 
    size_t getQueueCapacity() const; 
    uint64_t getQueueNum() const;
    uint64_t getLargeListNum() const;

    void dumpHeap() const;

    bool isOpen() const;
    void *getEntryHook() const;
    void setEntryHook(const void *entry_p);

    // static method --- 
    static void setDefaultOptionValue(init_option_st &optionst);
    static size_t getAlignSize(size_t size);

  private:
    //    static const uint64_t VERSION = 5;
    
    class Impl;
    std::unique_ptr<Impl> _impl;
  };

  std::string getErrorStr(int32_t err_num);
}
