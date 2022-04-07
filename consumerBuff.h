#include<stdio.h>
#include<vector>
#include<stdint.h>
#include<string>
#include<cstring>
#include<iostream>
#ifndef _CONSUMERBUFF_H_
#define _CONSUMERBUFF_H_
class consumerBuff{
    public:
        int Rank;//消费者的序号
        char *strobemerBuff;//每次都生成的strobemer的缓冲区
        int buffsize;
        std::vector<uint64_t> kmer_hashes; //计算strobemer需要的vector
        std::vector<unsigned int> position_ofSeq;//计算strobemer需要的vector
        consumerBuff(int _Rank,int buffersize);
        ~consumerBuff();
        void reset(int size);
    
        
};
 
#endif
