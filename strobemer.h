#ifndef __STROBEMER_H_
#define __STROBEMER_H_

#include <stdint.h>
#include <string>
#include<vector>
#include"io/Formater.h"
#include"consumerBuff.h"

// type of strobemer
enum strobemer_type{
    unknow=0,
    minstrobe=1,
    randstrobe=2,
    hybridstrobe=3,
};

class strobemer {
    public:
        static void init(int n, int k, int w_min, int w_max , strobemer_type t) ;
        static void chop_strobemer(const char * seq,int len, strobemer * buff);
        //static void chop_minstrobe_byKmer(const char * seq,int len, strobemer * buff,int &validLength);
        static void chop_strobemer_byKmer(const char *seq,int len,strobemer *buff,int& validLength);
        static void chop_randstrobemers(Reference &data,consumerBuff &consumerbuff);
        static void chop_minstrobemers(Reference &data,consumerBuff &consumerbuff,bool write_flag);
        //static void chop_randstrobemers1(Reference &data,consumerBuff &consumerbuff);
        
        //static void chop_randstrobe_byKmer(const char * seq,int len,strobemer *buff,int &validLength);
        //static void chop_hybridstrobe_byKmer(const char * seq,int len, strobemer * buff,int &validLength);
        static int strobmer_span() { return span ; }
        static int nkmer;
        static int ksize;
        static uint64_t kmask;
    private:
        
        static int wsize; // w_max here
        static int wmin;  // w_min
        static int span ; //W=k+2w=7+2*12=31
        static int kspan; //kspan=w-k+1; define a minimizer of how many kmers.
        static strobemer_type type;
        //kmask=(1ULL<<2*k)-1
    public:
        bool valid ;
        strobemer(); 
        ~strobemer();
        strobemer(const strobemer & other);
        strobemer & operator =(const strobemer & other);
        void print() const;
        std::string to_string() const ;
    private:
        char * kmer_forward ;//[nkmer][ksize];
    private:
        static void chop_minstrobe(const char * seq,int len, strobemer * buff);
        static void chop_minstrobe_byKmer(const char * seq,int len, strobemer * buff,int &validLength);
        static void chop_randstrobe(const char * seq,int len, strobemer * buff);
        static void chop_randstrobe_byKmer(const char * seq,int len,strobemer *buff,int &validLength);
        static void chop_hybridstrobe(const char * seq,int len, strobemer * buff);
        static void chop_hybridstrobe_byKmer(const char * seq,int len, strobemer * buff,int &validLength);
}; // class strobemer


#endif //__MINSTROBE_H_
