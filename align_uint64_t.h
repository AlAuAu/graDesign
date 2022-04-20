#include<x86intrin.h>
#include<stdint.h>

#ifndef _ALIGN_UINT64_T_H_
#define _ALIGN_UINT64_T_H_
#define ALIGNED_OPERATOR_NEW \
  void* operator new (std::size_t count) { \
    void* original = ::operator new(count + 32); \
    void* aligned = reinterpret_cast<void*>((reinterpret_cast<size_t>(original) & ~size_t(32 - 1)) + 32); \
    *(reinterpret_cast<void**>(aligned) - 1) = original; \
    return aligned;\
  } \
  void operator delete (void* ptr) { \
    ::operator delete(*(reinterpret_cast<void**>(ptr) - 1)); \
  }

class align_uint64_t{
    public:
            uint64_t value;
            align_uint64_t(uint64_t _value){
                 value=_value;
            }
            align_uint64_t(){
             
            }
           
           ALIGNED_OPERATOR_NEW
          


};
#endif

