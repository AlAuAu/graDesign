#include<stdio.h>
#include<vector>
#include<stdint.h>
#include<string>
#include<cstring>
#include<iostream>
#include"align_uint64_t.h"
#ifndef _CONSUMERBUFF_H_
#define _CONSUMERBUFF_H_
class consumerBuff{
    public:
        int Rank;//消费者的序号
        FILE * outputStream;
        //char *strobemerBuff;//每次都生成的strobemer的缓冲区
        int *positionSet;//生成的strobemers的strobe在kmer中的大小 行数为validlength 列数为nkmer;
        //int buffsize;
        int hashSetSize;//记录kmerhash数组的大小的
        int buffRowsize;
        int buffColsize;
        int validLength;
        align_uint64_t *kmer_hashes; //计算strobemer需要的vector
        std::vector<unsigned int> position_ofSeq;//计算strobemer需要的vector
        


        consumerBuff(int _Rank,int _buffRowsize, int _buffColsize,std::string &filename);
        ~consumerBuff();
        void reset(int size);
        void reserve(int _validlength);
        void write_to_file();
    
        
};
#endif
