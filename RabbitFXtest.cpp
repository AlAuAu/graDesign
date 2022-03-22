#include<string>
#include"io/FastxStream.h"
#include<vector>
#include<thread>
#include"io/Formater.h"
#include"io/FastxChunk.h"
#include"io/DataQueue.h"
#include<iostream>
#include<pthread.h>

//临界区定义
pthread_mutex_t  output_mutex;

int producer_task(std::string file,rabbit::fa::FastaDataPool & fastapool,rabbit::core::TDataQueue<rabbit::fa::FastaChunk> &dq){
    rabbit::fa::FastaFileReader faFileReader(file,fastapool);
    rabbit::int64 n_chunks=0;

    while(true){
        rabbit::fa::FastaChunk* fachunk;
        fachunk = faFileReader.readNextChunkList();
        if (fachunk== NULL) break;
        n_chunks++;
        dq.Push(n_chunks,fachunk);
    }
    dq.SetCompleted();
    pthread_mutex_lock(&output_mutex);
    std::cout<<"file  "<<file<<"  has  "<<n_chunks<<"  chunks  "<<std::endl;
    pthread_mutex_unlock(&output_mutex);
    return 0;
}

int consumer_task(rabbit::fa::FastaDataPool& fastapool,rabbit::core::TDataQueue<rabbit::fa::FastaChunk> & dq){
    long line_sum=0;
    rabbit::int64 id=0;
    std::vector<Reference> data;
    rabbit::fa::FastaChunk *fachunk;
    data.resize(10000);
    while(dq.Pop(id,fachunk)){
       line_sum += rabbit::fa::chunkFormat(*fachunk,data);
       fastapool.Release(fachunk->chunk);
    }
    pthread_mutex_lock(&output_mutex);
    for(int i=0;i<data.size();i++){
        std::cout<<data[i].name;
    }
    std::cout<<"line_sum    "<<line_sum<<std::endl;
    pthread_mutex_unlock(&output_mutex);
    return 0;
}

int main(int argc,char** argv){
    //初始化临界区
    pthread_mutex_init(&output_mutex,NULL);
    std::string file="/home/old_home/ljj/graDesign/code/data/GCF_000003135.1_ASM313v1_genomic.fna";
    int th=3;
    rabbit::fa::FastaDataPool datapool(32,1<<22);
    rabbit::core::TDataQueue<rabbit::fa::FastaChunk> queue1(64,1);
    std::thread producer(producer_task,file,std::ref(datapool),std::ref(queue1));
    std::vector<std::thread> consumers;

    for(int i=0;i<th;i++){
        consumers.emplace_back(std::thread(consumer_task,std::ref(datapool),std::ref(queue1)));
    }

    producer.join();

    for(int t=0;t<th;t++){
        consumers[t].join();
    }
    pthread_mutex_destroy(&output_mutex);
    return 0;
}