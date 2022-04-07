#include"consumerBuff.h"
consumerBuff::consumerBuff(int _Rank,int _buffsize){
    Rank=_Rank;
    buffsize=_buffsize;
    strobemerBuff=new char[consumerBuff::buffsize];

}
consumerBuff::~consumerBuff(){
    //std::cout<<"调用consumerBuff的析构函数"<<std::endl;
    
    delete []strobemerBuff;
}

void::consumerBuff::reset(int size){
    kmer_hashes.clear();
    kmer_hashes.reserve(size);
    position_ofSeq.clear();
    position_ofSeq.reserve(size);
    memset(strobemerBuff,'\0',sizeof(strobemerBuff));
    
}