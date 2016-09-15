
//
// Copyright (C) 2015-2016 Yahoo! JAPAN Research
//
// This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
// To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
//

#pragma once

#include <iostream>
#include <unistd.h>
#include <stdint.h>
#include <cstring>

#define ABS_ADDR(x, y) (void *)(x + (char *)y);

#define USE_MMAP_MANAGER

namespace MemoryManager{
  static const bool MMAP_DEFAULT_ALLOW_EXPAND = false; 
  static const uint64_t MMAP_CNTL_FILE_RANGE = 16; 
  static const size_t MMAP_CNTL_FILE_SIZE = MMAP_CNTL_FILE_RANGE * sysconf(_SC_PAGESIZE); 
  static const uint64_t MMAP_MAX_FILE_NAME_LENGTH = 1024; 
  static const std::string MMAP_CNTL_FILE_SUFFIX = "c"; 
  
  static const size_t MMAP_LOWER_SIZE = 1; 
  static const size_t MMAP_MEMORY_ALIGN = 8;   
  static const size_t MMAP_MEMORY_ALIGN_EXP = 3; 

#ifndef MMANAGER_TEST_MODE
  static const uint64_t MMAP_MAX_UNIT_NUM = 1024; 
#else
  static const uint64_t MMAP_MAX_UNIT_NUM = 8; 
#endif

  
  static const uint64_t MMAP_FREE_LIST_NUM = 64;  
  
  
  static const uint64_t MMAP_FREE_QUEUE_SIZE = 1024; 

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


  typedef struct _boot_st{
    uint32_t version; 
    uint64_t reserve; 
    size_t size;      
  }boot_st;

  typedef struct _head_st{
    off_t break_p;    
    uint64_t chunk_num; 
    uint64_t reserve; 
  }head_st;

  
  typedef struct _free_list_st{
    off_t free_p;
    off_t free_last_p;
  }free_list_st;

  
  typedef struct _free_st{
    free_list_st large_list;
    free_list_st free_lists[MMAP_FREE_LIST_NUM];
  }free_st;

  
  typedef struct _free_queue_st{
    off_t data;
    size_t capacity;
    uint64_t tail;
  }free_queue_st;

  
  
  typedef struct _control_st{
    bool use_expand;   
    uint16_t unit_num;   
    uint16_t active_unit; 
    uint64_t reserve;   
    size_t base_size;   
    off_t entry_p;      
    option_reuse_t reuse_type;  
    free_st free_data;   
    free_queue_st free_queue; 
    head_st data_headers[MMAP_MAX_UNIT_NUM]; 
  }control_st;

  typedef struct _chunk_head_st{
    bool delete_flg; 
    uint16_t unit_id; 
    off_t free_next;  
    size_t size; 
  }chunk_head_st;

  typedef struct _init_option_st{
    bool use_expand;  
    option_reuse_t reuse_type; 
  }init_option_st;


  class MmapManager{
  public:
    //    MmapManager():_isOpen(false), _mmapCntlAddr(NULL), _mmapCntlHead(NULL);
    MmapManager();
    ~MmapManager(){}
    
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

    void dumpHeap();    


    bool isOpen() const {
      return _isOpen;
    }
    
    void *getEntryHook() const {
      return getAbsAddr(_mmapCntlHead->entry_p);
    }
    
    void setEntryHook(const void *entry_p){
      _mmapCntlHead->entry_p = getRelAddr(entry_p);
    }

    // static method --- 
    static void setDefaultOptionValue(init_option_st &optionst)
    {
      optionst.use_expand = MMAP_DEFAULT_ALLOW_EXPAND;
      optionst.reuse_type = REUSE_DATA_CLASSIFY;
    }

    static size_t getAlignSize(size_t size){
      if((size % MMAP_MEMORY_ALIGN) == 0){
	
	return size;
      }else{
	
	return ( (size >> MMAP_MEMORY_ALIGN_EXP ) + 1 ) * MMAP_MEMORY_ALIGN;	
      }
    }

  private:
    
    static const uint64_t VERSION = 5; 

    bool _isOpen;
    void *_mmapCntlAddr;
    control_st *_mmapCntlHead;
    std::string _filePath;         
    void *_mmapDataAddr[MMAP_MAX_UNIT_NUM];


    void _initBootStruct(boot_st &bst, size_t size) const;
    void _initFreeStruct(free_st &fst) const;
    void _initFreeQueue(free_queue_st &fqst) const;
    void _initControlStruct(control_st &cntlst, size_t size) const;

    void _setupChunkHead(chunk_head_st *chunk_head, const bool delete_flg, const uint16_t unit_id, const off_t free_next, const size_t size) const;
    bool _expandMemory();
    int32_t _formatFile(const std::string &filePath, size_t size) const;
    void _clearChunk(const off_t chunk_off) const;

    void _free_data_classify(const off_t p, const bool force_large_list = false) const;    
    off_t _reuse_data_classify(const size_t size, reuse_state_t &reuse_state, const bool force_large_list = false) const;
    void _free_data_queue(const off_t p);
    off_t _reuse_data_queue(const size_t size, reuse_state_t &reuse_state);
    void _free_data_queue_plus(const off_t p);
    off_t _reuse_data_queue_plus(const size_t size, reuse_state_t &reuse_state);
    
    bool _scanAllData(void *target, const check_statistics_t stats_type) const;

    void _upHeap(free_queue_st *free_queue, uint64_t index) const;
    void _downHeap(free_queue_st *free_queue)const;
    bool _insertHeap(free_queue_st *free_queue, const off_t p) const;
    bool _getHeap(free_queue_st *free_queue, off_t *p) const;
    size_t _getMaxHeapValue(free_queue_st *free_queue) const;
    
    bool _jointChunk(const off_t chunk_offset);
    void _divChunk(const off_t chunk_offset, const size_t size);
  };

  std::string getErrorStr(int32_t err_num);
}
