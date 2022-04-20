#include <zlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <regex>
#include <omp.h>
#include<chrono>
#include<math.h>
#include<pthread.h>
#include<vector>
#include<stdio.h>
//#include<unistd.h>

#include"consumerBuff.h"
#include "originCode/kseq.h"
#include "strobemer.h"
#include"io/FastxChunk.h"
#include"io/FastxStream.h"
#include"io/DataQueue.h"
#include"io/Formater.h"

KSEQ_INIT(gzFile, gzread)
/*
1看看chopstrobemer和另外的方法能不能合并
2.stl容器替换
3.枚举的替换
*/

//打印命令行的用法
void print_usage() {
    std::cerr << "\n";
    std::cerr << "Make sequence into strobemers\n";
    std::cerr << "\n";
    std::cerr << "main [options] <filename.fasta>\n";
    std::cerr << "options:\n";
    std::cerr << "\t-n INT number of strobes [2]\n";
    std::cerr << "\t-k INT strobe length, limited to 32 [20]\n";
    std::cerr << "\t-v INT strobe w_min offset [k+1]\n";
    std::cerr << "\t-w INT strobe w_max offset [70]\n";
    std::cerr << "\t-t INT number of threads [3]\n";
    std::cerr << "\t-o name of output tsv-file [output.fasta]\n";
    std::cerr << "\t-c Choice of protocol to use; minstrobes, hybridstrobes, randstrobes [randstrobes]. \n";
    std::cerr <<"\t-C Choice of method to use;0 means choping by Kmers and 1 means not[0]. \n";
    
//    std::cerr << "\t-u Produce NAMs only from unique strobemers (w.r.t. reference sequences). This provides faster mapping.\n";
}

//判断是否是以.fna或者.fasta结尾
bool check(char *s){
    int L;
    L = strlen(s);
    if (strcmp(s+L-3,"fna")==0 ) return true;
    if (strcmp(s+L-5,"fasta")==0) return true;
    if (strcmp(s+L-2,"fa")==0)  return true;
    return false;
}

//临界区定义
//pthread_mutex_t  output_mutex;
//int file_count;

//定义需要计算的时间
std::chrono::duration<double, std::micro> elapsed;
std::chrono::duration<double, std::micro> elapsed1;


int producer_task(std::string file,rabbit::fa::FastaDataPool & fastapool,rabbit::core::TDataQueue<rabbit::fa::FastaChunk> &dq){
    
    rabbit::int64 n_chunks=0;
    
    //std::cout<<"数据的path是"<<file<<std::endl;
    rabbit::fa::FastaFileReader faFileReader(file,fastapool);
    std::cout<<"数据的path是"<<file<<std::endl;
     while(true){
            rabbit::fa::FastaChunk* fachunk;
            fachunk = faFileReader.readNextChunk();
            if (fachunk== NULL) break;
            n_chunks++;
            
            dq.Push(n_chunks,fachunk);
    }
    
    faFileReader.Close();
    std::cout<<"生产者已经完成所有生产"<<std::endl;
    dq.SetCompleted();
    
    return 0;
}

int consumer_task(std::string output_file,rabbit::fa::FastaDataPool& fastapool,rabbit::core::TDataQueue<rabbit::fa::FastaChunk> & dq,consumerBuff &consumerbuff){
    std::thread::id thread_id=this_thread::get_id();
    std::stringstream sin;
    sin <<thread_id;
    

    long line_sum=0;
    int chunks_count=0;
    rabbit::int64 id=0;
    std::vector<Reference> data;
    rabbit::fa::FastaChunk *fachunk;
    data.resize(10000);
    
    while(dq.Pop(id,fachunk)){
       //std::cout<<"data的size为"<<data.size()<<std::endl;
       data.clear();
       line_sum += rabbit::fa::chunkFormat(*fachunk,data);
      //std::cout<<"line_sum是"<<line_sum<<std::endl;
       fastapool.Release(fachunk->chunk);
       chunks_count++;
       for(int i=0;i<data.size();i++){
            Reference current=data[i];
            
            int number =current.length-strobemer::strobmer_span()+1; 
            if(number<0) continue;
            
            auto start1 = std::chrono::high_resolution_clock::now();
            strobemer::chop_minstrobemers(current,consumerbuff,false);
            auto finish1 = std::chrono::high_resolution_clock::now();
            elapsed += finish1 - start1;

            // auto start2 = std::chrono::high_resolution_clock::now();
            // strobemer::chop_randstrobemers(current,consumerbuff);
            // auto finish2 = std::chrono::high_resolution_clock::now();
            // elapsed1 += finish2 - start2;
            // for(int i=0;i<consumerbuff.strobemers.size();i++){
            //     std::cout<<consumerbuff.strobemers[i]<<std::endl;
            // }
           
            
                  
        }

        // if(chunks_count%10==0){
           
        //     std::cout<<"线程"<<sin.str()<<"完成"<<chunks_count<<"chunks"<<std::endl;
            
        // }
    }
    
    std::vector<Reference>().swap(data);
    return 0;
}




//看看chop里的length是不是序列长度 看看seq的l是不是长度
int main(int argc, char *argv[])
{   

    // Record  start time
    auto start = std::chrono::high_resolution_clock::now(); 
    if (argc < 2) {
        print_usage();
        return 0;
    }

    // Default parameters
    std::string choice = "randstrobes";

    int n = 2;
    int k = 20;
    //std::string output_file_name = "./output/output.fasta";
    std::string output_path="/home/old_home/ljj/graDesign/code/output/";
    int w_min = k+1;
    int w_max = 70;
   
    int n_threads = 1;
    int opn = 1;

    int chopMethod=0;

    while (opn < argc) {
        bool flag = false;
        if (argv[opn][0] == '-') {
            if (argv[opn][1] == 'n') {
                n = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'k') {
                k = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'o') {
                output_path= argv[opn + 1];
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'v') {
                w_min = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'w') {
                w_max = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'c') {
                choice = argv[opn + 1];
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 't') {
                n_threads = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if(argv[opn][1]=='C'){
                chopMethod=std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            }else {
				std::cout<<"Please check you parameter"<<std::endl;
                print_usage();
				return 0;
            }
        }
        if (!flag)
            break;
    }

    //omp_set_num_threads(n_threads); // set number of threads in "parallel" blocks
    std::cout << "Using" << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "k: " << k << std::endl;
    std::cout << "w_min: " << w_min << std::endl;
    std::cout << "w_max: " << w_max << std::endl;
    std::cout << "t: " << n_threads << std::endl;
    std::cout << "output_path:"<<output_path<<std::endl;
	std::cout<<"choice:"<<choice<<std::endl;


	
	 if (choice=="randstrobes")
	 {   
	     strobemer::init(n,k,w_min,w_max,strobemer_type::randstrobe);  	 
	 }else if (choice=="minstrobes")
	 {  
	    strobemer::init(n,k,w_min,w_max,strobemer_type::minstrobe);  
	 }else if (choice=="hybridstrobes")
	 {   
	     strobemer::init(n,k,w_min,w_max,strobemer_type::hybridstrobe);  
	 }else{
		 std::cout<<"Please set the choice in {minstrobes,randstrobes,hybridestrobes}";
		 return 0;
	 }
	 

     
   
    std::string file=argv[opn];
    rabbit::fa::FastaDataPool datapool(256,1<<22);
    rabbit::core::TDataQueue<rabbit::fa::FastaChunk> queue1(128,1);

    std::vector<consumerBuff> consumerBuffSet;
    consumerBuffSet.reserve(n_threads);
    
    for(int i=0;i<n_threads;i++){
        std::string name=output_path+"output"+to_string(i)+".fna";
        consumerBuffSet.emplace_back(i,10000,strobemer::nkmer,name);
    }
       
        

    std::thread producer(producer_task,file,std::ref(datapool),std::ref(queue1));
    std::vector<std::thread> consumers;
        
        

    for(int i=0;i<n_threads;i++){
            consumers.emplace_back(std::thread(consumer_task,output_path,std::ref(datapool),std::ref(queue1),std::ref(consumerBuffSet[i])));
    }

    producer.join();
        //std::this_thread::sleep_for(std::chrono::milliseconds(2000));

    for(int t=0;t<n_threads;t++){
         consumers[t].join();
    }
	std::cout<<"所有的消费者已经返回"<<std::endl;
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> totalElapsed = finish - start;
    std::cout << "用randstrobemers计算完所有的strobemer需要 " << elapsed.count() << "微秒" <<  std::endl;
    //std::cout << "用randstrobemers1计算完所有的strobemer需要 " << elapsed1.count() << "微秒" <<  std::endl;
    std::cout << "总耗时 " << totalElapsed.count() << "秒" <<  std::endl;




    
	return 0;
}
