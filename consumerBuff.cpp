#include"consumerBuff.h"
consumerBuff::consumerBuff(int _Rank,int _buffRowsize, int _buffColsize,std::string &filename){
    Rank=_Rank;
    //buffsize=_buffsize;
    //strobemerBuff=new char[_buffsize];
    hashSetSize=10000;
    kmer_hashes=new align_uint64_t[10000];
    buffRowsize=_buffRowsize;
    buffColsize=_buffColsize;
    //std::cout<<"执行buffer的构建"<<std::endl;
    validLength=buffRowsize;
    positionSet=new int[buffRowsize*buffColsize];
    outputStream=fopen(filename.c_str(),"a");

}
consumerBuff::~consumerBuff(){
    //std::cout<<"调用consumerBuff的析构函数"<<std::endl;
    fclose(outputStream);
    delete [] kmer_hashes;
    //delete []strobemerBuff;
    delete []positionSet;
}

void::consumerBuff::reset(int size){
    if(hashSetSize<size){
        hashSetSize=size;
        delete [] kmer_hashes;
        kmer_hashes=new align_uint64_t[hashSetSize];
    }
    position_ofSeq.clear();
    position_ofSeq.reserve(size);
    //memset(strobemerBuff,'\0',sizeof(strobemerBuff));
    
    
}

void::consumerBuff::reserve(int _validlength){
    validLength=_validlength;
    if(validLength>buffRowsize){
        buffRowsize=_validlength;
        delete []positionSet;
        positionSet=new int[buffRowsize*buffColsize];
    }
     
}

