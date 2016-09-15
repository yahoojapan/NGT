
//
// Copyright (C) 2015-2016 Yahoo! JAPAN Research
//
// This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
// To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
//

#include <cstdio>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <sstream>
#include <cerrno>
#include <cassert>

#include "MmapManager.h"
#include "MmapManagerException.h"

namespace MemoryManager{

  MmapManager::MmapManager():_isOpen(false), _mmapCntlAddr(NULL), _mmapCntlHead(NULL)
  {
    for(uint64_t i = 0; i < MMAP_MAX_UNIT_NUM; ++i){
      _mmapDataAddr[i] = NULL;
    }
  }
  
  bool MmapManager::init(const std::string &filePath, size_t size, const init_option_st *optionst) const
  {
    try{
      const std::string controlFile = filePath + MMAP_CNTL_FILE_SUFFIX;
      
      struct stat st;
      if(stat(controlFile.c_str(), &st) == 0){
        return false;
      }
      if(filePath.length() > MMAP_MAX_FILE_NAME_LENGTH){
        std::cerr << "too long filepath" << std::endl;
        return false;
      }
      if((size % sysconf(_SC_PAGESIZE) != 0) || ( size < MMAP_LOWER_SIZE )){
        std::cerr << "input size error" << std::endl;
        return false;
      }
       
      int32_t fd = _formatFile(controlFile, MMAP_CNTL_FILE_SIZE);
      assert(fd >= 0);
      
      errno = 0;
      char *cntl_p = (char *)mmap(NULL, MMAP_CNTL_FILE_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
      if(cntl_p == MAP_FAILED){
        const std::string err_str = getErrorStr(errno);
        if(close(fd) == -1) std::cerr << controlFile << "[WARN] : filedescript cannot close" << std::endl;
        throw MmapManagerException(controlFile + " " + err_str);
      }
      if(close(fd) == -1) std::cerr << controlFile << "[WARN] : filedescript cannot close" << std::endl;
      
      try {
	fd = _formatFile(filePath, size);
      } catch (MmapManagerException &err) {
        if(munmap(cntl_p, MMAP_CNTL_FILE_SIZE) == -1) {
	  throw MmapManagerException("[ERR] : munmap error : " + getErrorStr(errno) +
				     " : Through the exception : " + err.what());
	}
        throw err;
      }
      if(close(fd) == -1) std::cerr << controlFile << "[WARN] : filedescript cannot close" << std::endl;
      
      boot_st bootStruct = {0};
      control_st controlStruct = {0};
      _initBootStruct(bootStruct, size); 
      _initControlStruct(controlStruct, size);     
      
      char *cntl_head = cntl_p;
      cntl_head += sizeof(boot_st);

      if(optionst != NULL){
        controlStruct.use_expand = optionst->use_expand;
        controlStruct.reuse_type = optionst->reuse_type;
      }
      
      memcpy(cntl_p, (char *)&bootStruct, sizeof(boot_st));
      memcpy(cntl_head, (char *)&controlStruct, sizeof(control_st));
      
      errno = 0;
      if(munmap(cntl_p, MMAP_CNTL_FILE_SIZE) == -1) throw MmapManagerException("munmap error : " + getErrorStr(errno));
      
      return true;
    }catch(MmapManagerException e){
      std::cerr << "init error. " << e.what() << std::endl;
      throw e;
    }
  }

  bool MmapManager::openMemory(const std::string &filePath)
  {
    try{
      if(_isOpen == true){
        std::string err_str = "[ERROR] : openMemory error (double open).";
        throw MmapManagerException(err_str);
      }
      
      const std::string controlFile = filePath + MMAP_CNTL_FILE_SUFFIX;
      _filePath = filePath; 
      
      int32_t fd;
      
      errno = 0;
      if((fd = open(controlFile.c_str(), O_RDWR, 0666)) == -1){
        const std::string err_str = getErrorStr(errno);
        throw MmapManagerException("file open error" + err_str);
      }
      
      errno = 0;
      boot_st *boot_p = (boot_st*)mmap(NULL, MMAP_CNTL_FILE_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
      if(boot_p == MAP_FAILED){
        const std::string err_str = getErrorStr(errno);
        if(close(fd) == -1) std::cerr << controlFile << "[WARN] : filedescript cannot close" << std::endl;
        throw MmapManagerException(controlFile + " " + err_str);
      }
      if(close(fd) == -1) std::cerr << controlFile << "[WARN] : filedescript cannot close" << std::endl;
      
      if(boot_p->version != VERSION){
        std::cerr << "[WARN] : version error" << std::endl;
        errno = 0;
        if(munmap(boot_p, MMAP_CNTL_FILE_SIZE) == -1) throw MmapManagerException("munmap error : " + getErrorStr(errno));
        throw MmapManagerException("MemoryManager version error");	
      }
      
      errno = 0;
      if((fd = open(filePath.c_str(), O_RDWR, 0666)) == -1){
        const std::string err_str = getErrorStr(errno);
        errno = 0;
        if(munmap(boot_p, MMAP_CNTL_FILE_SIZE) == -1) throw MmapManagerException("munmap error : " + getErrorStr(errno));
        throw MmapManagerException("file open error = " + std::string(filePath.c_str())  + err_str);
      }
      
      _mmapCntlHead = (control_st*)( (char *)boot_p + sizeof(boot_st)); 
      _mmapCntlAddr = (void *)boot_p;
      
      for(uint64_t i = 0; i < _mmapCntlHead->unit_num; i++){
        off_t offset = _mmapCntlHead->base_size * i;
        errno = 0;
        _mmapDataAddr[i] = mmap(NULL, _mmapCntlHead->base_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, offset);
        if(_mmapDataAddr[i] == MAP_FAILED){
	  if (errno == EINVAL) {
	    std::cerr << "MmapManager::openMemory: Fatal error. EINVAL" << std::endl
		      << "  If you use valgrind, this error might occur when the DB is created." << std::endl
		      << "  In the case of that, reduce bsize in SharedMemoryAllocator." << std::endl;
	    assert(errno != EINVAL);
	  }
          const std::string err_str = getErrorStr(errno);
          if(close(fd) == -1) std::cerr << controlFile << "[WARN] : filedescript cannot close" << std::endl;
          closeMemory(true); 
          throw MmapManagerException(err_str);
        }
      }
      if(close(fd) == -1) std::cerr << controlFile << "[WARN] : filedescript cannot close" << std::endl;
      
      _isOpen = true;
      return true;
    }catch(MmapManagerException e){
      std::cerr << "open error" << std::endl;
      throw e;
    }
  }

  void MmapManager::closeMemory(const bool force)
  {
    try{
      if(force || _isOpen){
        uint16_t count = 0;
        void *error_ids[MMAP_MAX_UNIT_NUM] = {0};
        for(uint16_t i = 0; i < _mmapCntlHead->unit_num; i++){
          if(_mmapDataAddr[i] != NULL){
            if(munmap(_mmapDataAddr[i], _mmapCntlHead->base_size) == -1){
              error_ids[i] = _mmapDataAddr[i];;
              count++;
            }
            _mmapDataAddr[i] = NULL;
          }
        }
        
        if(count > 0){
          std::string msg = "";
          
          for(uint16_t i = 0; i < count; i++){
            std::stringstream ss;
            ss <<  error_ids[i];
            msg += ss.str() + ", ";
          }
          throw MmapManagerException("unmap error : ids = " + msg);
        }
        
        if(_mmapCntlAddr != NULL){
          if(munmap(_mmapCntlAddr, MMAP_CNTL_FILE_SIZE) == -1) throw MmapManagerException("munmap error : " + getErrorStr(errno));
          _mmapCntlAddr = NULL;
        }
        _isOpen = false;
      }
    }catch(MmapManagerException e){
      std::cerr << "close error" << std::endl;
      throw e;
    }
  }

  off_t MmapManager::alloc(const size_t size, const bool not_reuse_flag)
  {
    try{
      if(!_isOpen){
        std::cerr << "not open this file" << std::endl;
        return -1;
      }
      if(size > _mmapCntlHead->base_size + sizeof(chunk_head_st)){
        std::cerr << "alloc size over. size=" << size << "." << std::endl;
        return -1;
      }

      size_t alloc_size = getAlignSize(size);

      if(!not_reuse_flag){
        if(  _mmapCntlHead->reuse_type == REUSE_DATA_CLASSIFY 
	    || _mmapCntlHead->reuse_type == REUSE_DATA_QUEUE
	    || _mmapCntlHead->reuse_type == REUSE_DATA_QUEUE_PLUS){	  
          off_t ret_offset;
          reuse_state_t reuse_state = REUSE_STATE_OK;
          ret_offset = reuse(alloc_size, reuse_state);
          if(reuse_state != REUSE_STATE_ALLOC){
            return ret_offset;
          }
        }
      }
      
      head_st *unit_header = &_mmapCntlHead->data_headers[_mmapCntlHead->active_unit];
      if((unit_header->break_p + sizeof(chunk_head_st) + alloc_size) >= _mmapCntlHead->base_size){
        if(_mmapCntlHead->use_expand == true){
          if(_expandMemory() == false){
            std::cerr << __func__ << ": cannot expand" << std::endl;
            return -1;
          }
          unit_header = &_mmapCntlHead->data_headers[_mmapCntlHead->active_unit];
        }else{
          std::cerr << __func__ << ": total size over" << std::endl;
          return -1;
        }
      }
      
      const off_t file_offset = _mmapCntlHead->active_unit * _mmapCntlHead->base_size;
      const off_t ret_p = file_offset + ( unit_header->break_p + sizeof(chunk_head_st) );
      
      chunk_head_st *chunk_head = (chunk_head_st*)(unit_header->break_p + (char *)_mmapDataAddr[_mmapCntlHead->active_unit]);
      _setupChunkHead(chunk_head, false, _mmapCntlHead->active_unit, -1, alloc_size);
      unit_header->break_p += alloc_size + sizeof(chunk_head_st);
      unit_header->chunk_num++;
      
      return ret_p;
    }catch(MmapManagerException e){
      std::cerr << "allocation error" << std::endl;
      throw e;
    }
  }

  void MmapManager::free(const off_t p)
  {
    switch(_mmapCntlHead->reuse_type){
    case REUSE_DATA_CLASSIFY:
      _free_data_classify(p);
      break;
    case REUSE_DATA_QUEUE:
      _free_data_queue(p);
      break;
    case REUSE_DATA_QUEUE_PLUS:
      _free_data_queue_plus(p);
      break;
    default:
      _free_data_classify(p);
      break;
    }
  }

  off_t MmapManager::reuse(const size_t size, reuse_state_t &reuse_state)
  {
    off_t ret_off;

    switch(_mmapCntlHead->reuse_type){
    case REUSE_DATA_CLASSIFY:
      ret_off = _reuse_data_classify(size, reuse_state);
      break;
    case REUSE_DATA_QUEUE:
      ret_off = _reuse_data_queue(size, reuse_state);
      break;
    case REUSE_DATA_QUEUE_PLUS:
      ret_off = _reuse_data_queue_plus(size, reuse_state);
      break;
    default:
      ret_off = _reuse_data_classify(size, reuse_state);
      break;
    }

    return ret_off;
  }

  void *MmapManager::getAbsAddr(off_t p) const 
  {
    if(p < 0){
      return NULL;
    }
    const uint16_t unit_id = p / _mmapCntlHead->base_size;
    const off_t file_offset = unit_id * _mmapCntlHead->base_size;
    const off_t ret_p = p - file_offset;

    return ABS_ADDR(ret_p, _mmapDataAddr[unit_id]);
  }
  
  off_t MmapManager::getRelAddr(const void *p) const 
  {
    const chunk_head_st *chunk_head = (chunk_head_st *)((char *)p - sizeof(chunk_head_st));
    const uint16_t unit_id = chunk_head->unit_id;
    const off_t file_offset = unit_id * _mmapCntlHead->base_size;
    off_t ret_p = (off_t)((char *)p - (char *)_mmapDataAddr[unit_id]);
    ret_p += file_offset;

    return ret_p;
  }

  std::string getErrorStr(int32_t err_num){
    char err_msg[256];
    strerror_r(err_num, err_msg, 256);
    return std::string(err_msg);
  }

  size_t MmapManager::getTotalSize() const
  {
    const uint16_t active_unit = _mmapCntlHead->active_unit;
    const size_t ret_size = ((_mmapCntlHead->unit_num - 1) * _mmapCntlHead->base_size) + _mmapCntlHead->data_headers[active_unit].break_p;

    return ret_size;
  }

  size_t MmapManager::getUseSize() const
  {
    size_t total_size = 0;
    void *ref_addr = (void *)&total_size;
    _scanAllData(ref_addr, CHECK_STATS_USE_SIZE);

    return total_size;
  }

  uint64_t MmapManager::getUseNum() const
  {
    uint64_t total_chunk_num = 0;
    void *ref_addr = (void *)&total_chunk_num;
    _scanAllData(ref_addr, CHECK_STATS_USE_NUM);

    return total_chunk_num;
  }
  
  size_t MmapManager::getFreeSize() const
  {
    size_t total_size = 0;
    void *ref_addr = (void *)&total_size;
    _scanAllData(ref_addr, CHECK_STATS_FREE_SIZE);

    return total_size;
  }

  uint64_t MmapManager::getFreeNum() const
  {
    uint64_t total_chunk_num = 0;
    void *ref_addr = (void *)&total_chunk_num;
    _scanAllData(ref_addr, CHECK_STATS_FREE_NUM);

    return total_chunk_num;
  }

  uint16_t MmapManager::getUnitNum() const
  {
    return _mmapCntlHead->unit_num;
  }

  size_t MmapManager::getQueueCapacity() const
  {
    free_queue_st *free_queue = &_mmapCntlHead->free_queue;
    return free_queue->capacity;
  }

  uint64_t MmapManager::getQueueNum() const
  {
    free_queue_st *free_queue = &_mmapCntlHead->free_queue;
    return free_queue->tail;
  }


  void MmapManager::_initBootStruct(boot_st &bst, size_t size) const 
  {
    bst.version = VERSION;
    bst.reserve = 0;
    bst.size = size;
  }
  
  void MmapManager::_initFreeStruct(free_st &fst) const
  {
    fst.large_list.free_p = -1;
    fst.large_list.free_last_p = -1;	
    for(uint32_t i = 0; i < MMAP_FREE_LIST_NUM; ++i){
      fst.free_lists[i].free_p = -1;
      fst.free_lists[i].free_last_p = -1;	
    }
  }
  
  void MmapManager::_initFreeQueue(free_queue_st &fqst) const
  {
    fqst.data = -1;
    fqst.capacity = MMAP_FREE_QUEUE_SIZE;
    fqst.tail = 1;
  }
  
  void MmapManager::_initControlStruct(control_st &cntlst, size_t size) const
  {
    cntlst.use_expand = MMAP_DEFAULT_ALLOW_EXPAND;
    cntlst.unit_num = 1;
    cntlst.active_unit = 0;
    cntlst.reserve = 0;
    cntlst.base_size = size;
    cntlst.entry_p = 0;
    cntlst.reuse_type = REUSE_DATA_CLASSIFY;
    _initFreeStruct(cntlst.free_data);
    _initFreeQueue(cntlst.free_queue);
    memset(cntlst.data_headers, 0, sizeof(head_st) * MMAP_MAX_UNIT_NUM);
  }
  
  void MmapManager::_setupChunkHead(chunk_head_st *chunk_head, const bool delete_flg, const uint16_t unit_id, const off_t free_next, const size_t size) const
  {
      chunk_head_st chunk_buffer;
      chunk_buffer.delete_flg = delete_flg;
      chunk_buffer.unit_id    = unit_id;
      chunk_buffer.free_next  = free_next;
      chunk_buffer.size       = size;
      
      memcpy(chunk_head, &chunk_buffer, sizeof(chunk_head_st));
  }

  bool MmapManager::_expandMemory()
  {
    const uint16_t new_unit_num = _mmapCntlHead->unit_num + 1;
    const size_t new_file_size = _mmapCntlHead->base_size * new_unit_num;
    const off_t old_file_size = _mmapCntlHead->base_size * _mmapCntlHead->unit_num;

    if(new_unit_num >= MMAP_MAX_UNIT_NUM){
      std::cerr << "over max unit num" << std::endl;
      return false;
    }

    int32_t fd = _formatFile(_filePath, new_file_size);
    assert(fd >= 0);

    const off_t offset = _mmapCntlHead->base_size * _mmapCntlHead->unit_num;
    errno = 0;
    void *new_area = mmap(NULL, _mmapCntlHead->base_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, offset);
    if(new_area == MAP_FAILED){
      const std::string err_str = getErrorStr(errno);
      
      errno = 0;
      if(ftruncate(fd, old_file_size) == -1){
        const std::string err_str = getErrorStr(errno);
        throw MmapManagerException("truncate error" + err_str);
      }
      
      if(close(fd) == -1) std::cerr << _filePath << "[WARN] : filedescript cannot close" << std::endl;
      throw MmapManagerException("mmap error" + err_str);
    }
    if(close(fd) == -1) std::cerr << _filePath << "[WARN] : filedescript cannot close" << std::endl;
    
    _mmapDataAddr[_mmapCntlHead->unit_num] = new_area;
    
    _mmapCntlHead->unit_num = new_unit_num;
    _mmapCntlHead->active_unit++;
    
    return true;
  }

  int32_t MmapManager::_formatFile(const std::string &filePath, size_t size) const
  {
    const char *c = "";
    int32_t fd;
    
    errno = 0;
    if((fd = open(filePath.c_str(), O_RDWR|O_CREAT, 0666)) == -1){
      std::stringstream ss;
      ss << "[ERR] Cannot open the file. " << filePath << " " << getErrorStr(errno);
      throw MmapManagerException(ss.str());
    }
    errno = 0;
    if(lseek(fd, (off_t)size-1, SEEK_SET) < 0){
      std::stringstream ss;
      ss << "[ERR] Cannot seek the file. " << filePath << " " << getErrorStr(errno);
      if(close(fd) == -1) std::cerr << filePath << "[WARN] : filedescript cannot close" << std::endl;
      throw MmapManagerException(ss.str());
    }
    errno = 0;
    if(write(fd, &c, sizeof(char)) == -1){
      std::stringstream ss;
      ss << "[ERR] Cannot write the file. Check the disk space. " << filePath << " " << getErrorStr(errno);
      if(close(fd) == -1) std::cerr << filePath << "[WARN] : filedescript cannot close" << std::endl;
      throw MmapManagerException(ss.str());
    }
    
    return fd;
  }

  void MmapManager::_clearChunk(const off_t chunk_off) const
  {
    chunk_head_st *chunk_head = (chunk_head_st *)getAbsAddr(chunk_off);
    const off_t payload_off = chunk_off + sizeof(chunk_head_st);
    
    chunk_head->delete_flg = false;
    chunk_head->free_next = -1;
    char *payload_addr = (char *)getAbsAddr(payload_off);
    memset(payload_addr, 0, chunk_head->size);
  }


  void MmapManager::_free_data_classify(const off_t p, const bool force_large_list) const
  {
    const off_t chunk_offset = p - sizeof(chunk_head_st);
    chunk_head_st *chunk_head = (chunk_head_st *)getAbsAddr(chunk_offset);
    const size_t p_size = chunk_head->size;

    
    
    const size_t border_size = MMAP_MEMORY_ALIGN * MMAP_FREE_LIST_NUM;

    free_list_st *free_list;
    if(p_size <= border_size && force_large_list == false){
      uint32_t index = (p_size / MMAP_MEMORY_ALIGN) - 1; 
      free_list = &_mmapCntlHead->free_data.free_lists[index];
    }else{
      free_list = &_mmapCntlHead->free_data.large_list;
    }
    
    if(free_list->free_p == -1){
      free_list->free_p = free_list->free_last_p = chunk_offset;
    }else{
      off_t last_off = free_list->free_last_p;
      chunk_head_st *tmp_chunk_head = (chunk_head_st *)getAbsAddr(last_off);
      free_list->free_last_p = tmp_chunk_head->free_next = chunk_offset;
    }
    chunk_head->delete_flg = true;
  }

  off_t MmapManager::_reuse_data_classify(const size_t size, reuse_state_t &reuse_state, const bool force_large_list) const
  {
    
    
    const size_t border_size = MMAP_MEMORY_ALIGN * MMAP_FREE_LIST_NUM;

    free_list_st *free_list;
    if(size <= border_size && force_large_list == false){
      uint32_t index = (size / MMAP_MEMORY_ALIGN) - 1; 
      free_list = &_mmapCntlHead->free_data.free_lists[index];
    }else{
      free_list = &_mmapCntlHead->free_data.large_list;
    }
    
    if(free_list->free_p == -1){
      reuse_state = REUSE_STATE_ALLOC;
      return -1;
    }

    off_t current_off = free_list->free_p;
    off_t ret_off = 0;
    chunk_head_st *current_chunk_head = (chunk_head_st *)getAbsAddr(current_off);
    chunk_head_st *ret_chunk_head = NULL;

    if( (size <= border_size) && (free_list->free_last_p == free_list->free_p) ){
      ret_off = current_off;
      ret_chunk_head = current_chunk_head;
      free_list->free_p = free_list->free_last_p = -1;
    }else{
      off_t ret_before_off = -1, before_off = -1;
      bool found_candidate_flag = false;

      
      while(current_chunk_head != NULL){
	if( current_chunk_head->size >= size ) found_candidate_flag = true;

	if(found_candidate_flag){
	  ret_off = current_off;
	  ret_chunk_head = current_chunk_head;
	  ret_before_off = before_off;
	  break;
	}
	before_off = current_off;
	current_off = current_chunk_head->free_next;
	current_chunk_head = (chunk_head_st *)getAbsAddr(current_off);
      }

      if(!found_candidate_flag){
	reuse_state = REUSE_STATE_ALLOC;
	return -1;
      }

      const off_t free_next = ret_chunk_head->free_next;
      if(free_list->free_p == ret_off){
	free_list->free_p = free_next;
      }else{
	chunk_head_st *before_chunk = (chunk_head_st *)getAbsAddr(ret_before_off);
	before_chunk->free_next = free_next;
      }

      if(free_list->free_last_p == ret_off){
	free_list->free_last_p = ret_before_off;
      }
    }

    _clearChunk(ret_off);

    ret_off = ret_off + sizeof(chunk_head_st);
    return ret_off;
  }

  void MmapManager::_free_data_queue(const off_t p)
  {
    free_queue_st *free_queue = &_mmapCntlHead->free_queue;
    if(free_queue->data == -1){
      
      const size_t queue_size = sizeof(off_t) * free_queue->capacity;
      const off_t alloc_offset = alloc(queue_size);
      if(alloc_offset == -1){
	
	return _free_data_classify(p, true);
      }
      free_queue->data = alloc_offset;
    }else if(free_queue->tail >= free_queue->capacity){
      
      const off_t tmp_old_queue = free_queue->data;
      const size_t old_size = sizeof(off_t) * free_queue->capacity;
      const size_t new_capacity = free_queue->capacity * 2;
      const size_t new_size = sizeof(off_t) * new_capacity;

      if(new_size > _mmapCntlHead->base_size){
	
	
	return _free_data_classify(p, true);
      }else{
	const off_t alloc_offset = alloc(new_size);	
	if(alloc_offset == -1){
	  
	  return _free_data_classify(p, true);
	}
	free_queue->data = alloc_offset;
	const off_t *old_data = (off_t *)getAbsAddr(tmp_old_queue);
	off_t *new_data = (off_t *)getAbsAddr(free_queue->data);
	memcpy(new_data, old_data, old_size);
      
	free_queue->capacity = new_capacity;
	this->free(tmp_old_queue);
      }
    }

    const off_t chunk_offset = p - sizeof(chunk_head_st);
    if(!_insertHeap(free_queue, chunk_offset)){
      
      return;
    }

    chunk_head_st *chunk_head = (chunk_head_st*)getAbsAddr(chunk_offset);
    chunk_head->delete_flg = 1;
    
    return;
  }
  
  off_t MmapManager::_reuse_data_queue(const size_t size, reuse_state_t &reuse_state)
  {
    free_queue_st *free_queue = &_mmapCntlHead->free_queue;    
    if(free_queue->data == -1){    
      
      reuse_state = REUSE_STATE_ALLOC;
      return -1;
    }

    if(_getMaxHeapValue(free_queue) < size){
      reuse_state = REUSE_STATE_ALLOC;
      return -1;
    }

    off_t ret_off;
    if(!_getHeap(free_queue, &ret_off)){
      
      reuse_state = REUSE_STATE_ALLOC;
      return -1;
    }

    
    reuse_state_t list_state = REUSE_STATE_OK;
    
    off_t candidate_off = _reuse_data_classify(MMAP_MEMORY_ALIGN, list_state, true);
    if(list_state == REUSE_STATE_OK){
      
      this->free(candidate_off);
    }
    
    const off_t c_ret_off = ret_off;
    _divChunk(c_ret_off, size);

    _clearChunk(ret_off);

    ret_off = ret_off + sizeof(chunk_head_st);

    return ret_off;
  }

  void MmapManager::_free_data_queue_plus(const off_t p)
  {
    const off_t chunk_offset = p - sizeof(chunk_head_st);
    chunk_head_st *chunk_head = (chunk_head_st *)getAbsAddr(chunk_offset);
    const size_t p_size = chunk_head->size;

    
    
    const size_t border_size = MMAP_MEMORY_ALIGN * MMAP_FREE_LIST_NUM;

    if(p_size <= border_size){
      _free_data_classify(p);
    }else{
      _free_data_queue(p);
    }
  }

  off_t MmapManager::_reuse_data_queue_plus(const size_t size, reuse_state_t &reuse_state)
  {
    
    
    const size_t border_size = MMAP_MEMORY_ALIGN * MMAP_FREE_LIST_NUM;

    off_t ret_off;
    if(size <= border_size){
      ret_off = _reuse_data_classify(size, reuse_state);
      if(reuse_state == REUSE_STATE_ALLOC){
	
	reuse_state = REUSE_STATE_OK;
	ret_off = _reuse_data_queue(size, reuse_state);	
      }
    }else{
      ret_off = _reuse_data_queue(size, reuse_state);
    }
    
    return ret_off;
  }
  

  bool MmapManager::_scanAllData(void *target, const check_statistics_t stats_type) const
  {
    const uint16_t unit_num = _mmapCntlHead->unit_num;
    size_t total_size = 0;
    uint64_t total_chunk_num = 0;

    for(int i = 0; i < unit_num; i++){
      const head_st *target_unit_head = &_mmapCntlHead->data_headers[i];
      const uint64_t chunk_num = target_unit_head->chunk_num;
      const off_t base_offset = i * _mmapCntlHead->base_size;
      off_t target_offset = base_offset;
      chunk_head_st *target_chunk;

      for(uint64_t j = 0; j < chunk_num; j++){
        target_chunk = (chunk_head_st*)getAbsAddr(target_offset);

        if(stats_type == CHECK_STATS_USE_SIZE){
          if(target_chunk->delete_flg == false){
            total_size += target_chunk->size;
          }
        }else if(stats_type == CHECK_STATS_USE_NUM){
          if(target_chunk->delete_flg == false){
            total_chunk_num++;
          }
        }else if(stats_type == CHECK_STATS_FREE_SIZE){
          if(target_chunk->delete_flg == true){
            total_size += target_chunk->size;
          }
        }else if(stats_type == CHECK_STATS_FREE_NUM){
          if(target_chunk->delete_flg == true){
            total_chunk_num++;
          }
        }

        const size_t chunk_size = sizeof(chunk_head_st) + target_chunk->size;
        target_offset += chunk_size;
      }
    }

    if(stats_type == CHECK_STATS_USE_SIZE || stats_type == CHECK_STATS_FREE_SIZE){
      size_t *tmp_size = (size_t *)target;
      *tmp_size = total_size;
    }else if(stats_type == CHECK_STATS_USE_NUM || stats_type == CHECK_STATS_FREE_NUM){
      uint64_t *tmp_chunk_num = (uint64_t *)target;
      *tmp_chunk_num = total_chunk_num;
    }
      
    return true;
  }

  void MmapManager::_upHeap(free_queue_st *free_queue, uint64_t index) const 
  {
    off_t *queue = (off_t *)getAbsAddr(free_queue->data);

    while(index > 1){
      uint64_t parent = index / 2;
      
      const off_t parent_chunk_offset = queue[parent];
      const off_t index_chunk_offset = queue[index];
      const chunk_head_st *parent_chunk_head = (chunk_head_st *)getAbsAddr(parent_chunk_offset);
      const chunk_head_st *index_chunk_head = (chunk_head_st *)getAbsAddr(index_chunk_offset);    
      
      if(parent_chunk_head->size < index_chunk_head->size){
	
	const off_t tmp = queue[parent];
	queue[parent] = queue[index];
	queue[index] = tmp;
      }
      index = parent;
    }
  }

  void MmapManager::_downHeap(free_queue_st *free_queue)const
  {
    off_t *queue = (off_t *)getAbsAddr(free_queue->data);
    uint64_t index = 1;    

    while(index * 2 <= free_queue->tail){
      uint64_t child = index * 2;

      const off_t index_chunk_offset = queue[index];
      const chunk_head_st *index_chunk_head = (chunk_head_st *)getAbsAddr(index_chunk_offset);          

      if(child + 1 < free_queue->tail){
	const off_t left_chunk_offset = queue[child];
	const off_t right_chunk_offset = queue[child+1];
	const chunk_head_st *left_chunk_head = (chunk_head_st *)getAbsAddr(left_chunk_offset);
	const chunk_head_st *right_chunk_head = (chunk_head_st *)getAbsAddr(right_chunk_offset);

	
	if(left_chunk_head->size < right_chunk_head->size){
	  child = child + 1;
	}
      }

      
      const off_t child_chunk_offset = queue[child];
      const chunk_head_st *child_chunk_head = (chunk_head_st *)getAbsAddr(child_chunk_offset);

      if(child_chunk_head->size > index_chunk_head->size){
	
	const off_t tmp = queue[child];
	queue[child] = queue[index];
	queue[index] = tmp;
	index = child;
      }else{
	break;
      }
    }
  }
  
  bool MmapManager::_insertHeap(free_queue_st *free_queue, const off_t p) const 
  {
    off_t *queue = (off_t *)getAbsAddr(free_queue->data);    
    uint64_t index;
    if(free_queue->capacity < free_queue->tail){
      return false;
    }

    index = free_queue->tail;
    queue[index] = p;
    free_queue->tail += 1;

    _upHeap(free_queue, index);
    
    return true;
  }

  bool MmapManager::_getHeap(free_queue_st *free_queue, off_t *p) const
  {
    
    if( (free_queue->tail - 1) <= 0){     
      return false;
    }

    off_t *queue = (off_t *)getAbsAddr(free_queue->data);        
    *p = queue[1];
    free_queue->tail -= 1;
    queue[1] = queue[free_queue->tail];
    _downHeap(free_queue);

    return true;
  }

  size_t MmapManager::_getMaxHeapValue(free_queue_st *free_queue) const
  {
    if(free_queue->data == -1){
      return 0;
    }
    const off_t *queue = (off_t *)getAbsAddr(free_queue->data);
    const chunk_head_st *chunk_head = (chunk_head_st *)getAbsAddr(queue[1]);
    
    return chunk_head->size;
  }

  void MmapManager::dumpHeap()
  {
    free_queue_st *free_queue = &_mmapCntlHead->free_queue;
    if(free_queue->data == -1){
      std::cout << "heap unused" << std::endl;
      return;
    }
    
    off_t *queue = (off_t *)getAbsAddr(free_queue->data);        
    for(uint32_t i = 1; i < free_queue->tail; ++i){
      const off_t chunk_offset = queue[i];
      const off_t payload_offset = chunk_offset + sizeof(chunk_head_st);
      const chunk_head_st *chunk_head = (chunk_head_st *)getAbsAddr(chunk_offset);
      const size_t size = chunk_head->size;
      std::cout << "[" << chunk_offset << "(" << payload_offset << "), " << size << "] ";
    }
    std::cout << std::endl;
  }
  


  
  
  

  
  bool MmapManager::_jointChunk(const off_t chunk_offset)
  {
    return true;
  }

  void MmapManager::_divChunk(const off_t chunk_offset, const size_t size)
  {
    if((_mmapCntlHead->reuse_type != REUSE_DATA_QUEUE)
       && (_mmapCntlHead->reuse_type != REUSE_DATA_QUEUE_PLUS)){
      return;
    }
    
    chunk_head_st *chunk_head = (chunk_head_st *)getAbsAddr(chunk_offset);
    const size_t border_size = sizeof(chunk_head_st) + MMAP_MEMORY_ALIGN;
    const size_t align_size = getAlignSize(size);
    const size_t rest_size = chunk_head->size - align_size;

    if(rest_size < border_size){
      return;
    }

    
    chunk_head->size = align_size;

    const off_t new_chunk_offset = chunk_offset + sizeof(chunk_head_st) + align_size;
    chunk_head_st *new_chunk_head = (chunk_head_st *)getAbsAddr(new_chunk_offset);
    const size_t new_size = rest_size - sizeof(chunk_head_st);
    _setupChunkHead(new_chunk_head, true, chunk_head->unit_id, -1, new_size);
    
    
    head_st *unit_header = &_mmapCntlHead->data_headers[_mmapCntlHead->active_unit];    
    unit_header->chunk_num++;    

    
    const off_t payload_offset = new_chunk_offset + sizeof(chunk_head_st);
    this->free(payload_offset);
    
    return;
  }
}
