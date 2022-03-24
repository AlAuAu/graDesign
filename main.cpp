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


#include "originCode/kseq.h"
#include "originCode/strobemer.h"
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
    
//    std::cerr << "\t-u Produce NAMs only from unique strobemers (w.r.t. reference sequences). This provides faster mapping.\n";
}

//判断是否是以.fna或者.fasta结尾
bool check(char *s){
    int L;
    L = strlen(s);
    if (strcmp(s+L-3,"fna")==0 ) return true;
    if (strcmp(s+L-3,"fasta")==0) return true;
    return false;
}

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
    //pthread_mutex_lock(&output_mutex);
    //std::cout<<"file  "<<file<<"  has  "<<n_chunks<<"  chunks  "<<std::endl;
    //pthread_mutex_unlock(&output_mutex);
    return 0;
}

int consumer_task(std::string output_file,rabbit::fa::FastaDataPool& fastapool,rabbit::core::TDataQueue<rabbit::fa::FastaChunk> & dq,int  n_threads){

    long line_sum=0;
    rabbit::int64 id=0;
    std::vector<Reference> data;
    rabbit::fa::FastaChunk *fachunk;
    data.resize(10000);
    while(dq.Pop(id,fachunk)){
       line_sum += rabbit::fa::chunkFormat(*fachunk,data);
       fastapool.Release(fachunk->chunk);
    }

    

    
    // for(int j=0;j<data.size();j++){
    //     if(data[j].length!=0)
    //     std:cout<<data[j].length<<std::endl;
    // }
    #pragma omp parallel for num_threads(n_threads)
    for(int i=0;i<data.size();i++){
       
        Reference current=data[i];
        int number =current.length-strobemer::strobmer_span()+1; 

        if(number<0) continue;

        //std::cout<<"seq是："<<current.seq<<std::endl;
        // std::cout<<"length的值是："<<current.length<<std::endl;
        // std::cout<<"number的值是："<<number<<std::endl;
        // std::cout<<"span的值是："<<strobemer::strobmer_span()<<std::endl;
        int validLength=0;
        //std::cout<<"程序断点1"<<std::endl;
        strobemer * buff = new strobemer[number];
        //std::cout<<"程序断点2"<<std::endl;
        strobemer::chop_randstrobe_byKmer(current.seq.c_str(),current.length,buff,validLength);

        pthread_mutex_lock(&output_mutex);
        std::ofstream outputStream;
        outputStream.open(output_file);//这是多线程共享的吗 是不是需要临界区
        //std::cout<<"程序断点3"<<std::endl;
        outputStream<<">"<<current.name<<'\n';
        for(int j=0;j<=validLength;j++){
            //std::cout<<"validLength为"<<validLength<<std::endl;
            //if(buff[i].valid){
                 outputStream<<buff[j].to_string()<<'\n';
                 //std::cout<<buff[i].to_string()<<"\n";
            //}
           
        }
        //std::cout<<"程序断点4"<<std::endl;
        outputStream.close();
        pthread_mutex_unlock(&output_mutex);
        delete []buff;
        
        
    }
    
    //std::cout<<"line_sum    "<<line_sum<<std::endl;
    std::vector<Reference>().swap(data);//释放vector
    return 0;
}




//看看chop里的length是不是序列长度 看看seq的l是不是长度
int main(int argc, char *argv[])
{   
    if (argc < 2) {
        print_usage();
        return 0;
    }

    // Default parameters
    std::string choice = "randstrobes";

    int n = 2;
    int k = 20;
    //std::string output_file_name = "./output/output.fasta";
    std::string output_path="../output";
    int w_min = k+1;
    int w_max = 70;
    int n_threads = 3;
    int opn = 1;
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
            } else {
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
	 

     //初始化临界区
    pthread_mutex_init(&output_mutex,NULL);
   
	// Record  start time
    auto start = std::chrono::high_resolution_clock::now(); 
    
    DIR *dirp;
    struct dirent *direntp;
    dirp=opendir(argv[opn]);
   
    
	
    // gzFile fp;
    // kseq_t *seq;
    // int l;
     
   
    std::string path=argv[opn];
    path=path+'/';
    output_path=output_path+'/';
    
    
    long int count=0;
    long int total=500;
   
    // while (((direntp=readdir(dirp))!=NULL))
    // {
    //     total++;
    //     std::cout<<total<<std::endl;
    // }

    // rewinddir(dirp);

    
   
    
    //#pragma omp parallel for num_threads(n_threads)
    
    for (int i = 1; i <= total; i++){   
    //while (((direntp=readdir(dirp))!=NULL)){    
        direntp=readdir(dirp);
        //尝试读取fasta文件中的序列
        //std::cout<<direntp->d_name<<std::endl;
        if(!check(direntp->d_name)) continue;
        
        std::string d_name=direntp->d_name;

        rabbit::fa::FastaDataPool datapool(32,1<<22);
        rabbit::core::TDataQueue<rabbit::fa::FastaChunk> queue1(64,1);

        std::thread producer(producer_task,path+d_name,std::ref(datapool),std::ref(queue1));
        std::vector<std::thread> consumers;

        for(int i=0;i<n_threads;i++){
            consumers.emplace_back(std::thread(consumer_task,output_path+d_name,std::ref(datapool),std::ref(queue1),n_threads));
        }

        producer.join();

        for(int t=0;t<n_threads;t++){
            consumers[t].join();
            
        }
       
        
        // gzFile fp;
        // kseq_t *seq;
        // int l;
         
        
        //std::cout<<"执行到这,文件路径是："<<path+d_name<<std::endl;
        
        // std::ofstream output_file;
        // output_file.open(output_path+d_name);
        //std::cout<<"执行到这,文件路径是："<<output_path+d_name<<std::endl;

        // fp=gzopen((path+d_name).c_str(),"r");
        // seq= kseq_init(fp);
       
        // while ((l = kseq_read(seq)) >= 0) {
    
        //     int number = seq->seq.l-strobemer::strobmer_span()+1; 
        //     strobemer * buff = new strobemer[number];
        //     strobemer::chop_strobemer(seq->seq.s,seq->seq.l,buff);
        //     output_file<<">"<<seq->name.s<<'\n';
        //     //std::cout<<"序列名："<<seq->name.s<<std::endl;
        //     for(int i = 0 ; i< number ; i++ ){
        //         if(buff[i].valid)
        //             output_file<<buff[i].to_string()<<'\n'; // or do whatever you want ...
        //             //std::cout<<buff[i].to_string()<<'\n';;
        //     }
            
        //     delete [] buff;
	    // }
        // output_file.close();
        // kseq_destroy(seq);
	    // gzclose(fp);
        

        //int rank=omp_get_thread_num();
        //std::cout<<"线程"<<rank<<"已完成第"<<i<<"个文件"<<std::endl;
        if(i%100==0) {
            std::cout<<"已完成："<<i<<"/"<<total<<std::endl;
        }
   
    }
    closedir(dirp);
    pthread_mutex_destroy(&output_mutex);
    
    
    // while (((direntp=readdir(dirp))!=NULL)){
        
        
    //     //尝试读取fasta文件中的序列
    //     if(!check(direntp->d_name)) continue;
        
        
    //     std::string d_name=direntp->d_name;
    //     //std::cout<<"执行到这,文件路径是："<<path+d_name<<std::endl;
        
    //     std::ofstream output_file;
    //     output_file.open(output_path+d_name);
    //     //std::cout<<"执行到这,文件路径是："<<output_path+d_name<<std::endl;

    //     fp=gzopen((path+d_name).c_str(),"r");
    //     seq= kseq_init(fp);
    //     //test
    //     while ((l = kseq_read(seq)) >= 0) {
    //         // printf("name: %s\n", seq->name.s);//序列的序列号
    //         // if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
    //         // printf("seq: %s\n", seq->seq.s);//基因序列
    //         // if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
            
    //         int number = seq->seq.l-strobemer::strobmer_span()+1; 
    //         strobemer * buff = new strobemer[number];
    //         strobemer::chop_strobemer(seq->seq.s,seq->seq.l,buff);
    //         output_file<<">"<<seq->name.s<<'\n';
    //         //std::cout<<"序列名："<<seq->name.s<<std::endl;
    //         for(int i = 0 ; i< number ; i++ ){
    //             if(buff[i].valid)
    //                 output_file<<buff[i].to_string()<<'\n'; // or do whatever you want ...
    //                 //std::cout<<buff[i].to_string()<<'\n';;
    //         }
            
    //         delete [] buff;
	//     }
    //     output_file.close();	

    //     count++;
    //     if(count%100==0) {
    //         std::cout<<"已完成："<<count<<"/"<<total<<std::endl;
    //     }
    //     if (count==500) break;

        
        
    //}
    //std::cout<<"count:"<<count<<std::endl;
    
	//output_file.close();
	// kseq_destroy(seq);
	// gzclose(fp);

	auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_mers = finish - start;
    float rounded = truncf(elapsed_mers.count() * 10) / 10;
    std::cout << "总耗时 " << rounded << " s\n" <<  std::endl;



    
	return 0;
}