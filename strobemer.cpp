#include "strobemer.h"
#include <cassert>
#include <cstring>
#include <iostream>
#include<vector>



/**********************************************************
 *
 * hash kmer into uint64
 *
 * *******************************************************/
// copy from minimap2:sketch.c :
static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}//hash64

/**********************************************************
 *
 * char tables
 *
 * *******************************************************/
char complementary_base[256] = {
    'T', 'G', 'C', 'A',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G',  'N', 'N', 'N', 'C',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'A', 'A', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G',  'N', 'N', 'N', 'C',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'A', 'A', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
}; //complementary_base

static unsigned char seq_nt4_table[256] = {
    0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
}; //seq_nt4_table

/**********************************************************
 *
 * to get reverse-complementary sequences
 *
 * *******************************************************/
int inline rc_index(int index, int len){
    return len - 1 - index ;
}

void reverse_complete( const char * seq , int len, char * buff ){
    for( int i = 0 ; i < len ; i++ ){
        buff[i] = complementary_base[seq[rc_index(i,len)]];
    }
}


int inline check_index(__mmask8 checkseq){
  //用此函数返回mask从左到右第一位为1的位置（即局部最小值得位置）
  int index=0;
//   __mmask8 checkFlag=_cvtu32_mask8(1);
  unsigned char c = checkseq;
  c=(c&(-c))-1;
  index+=(c>>0)&1;
  index+=(c>>1)&1;
  index+=(c>>2)&1;
  index+=(c>>3)&1;
  index+=(c>>4)&1;
  index+=(c>>5)&1;
  index+=(c>>6)&1;
  index+=(c>>7)&1;
 
  return index;
};
/**********************************************************
 *
 * to get hashed-kmers
 *
 * *******************************************************/
class binary_kmer {
    public:
        static int      ksize;
        static uint64_t kmask;
    public:
        static void InitK(int k) {
            assert(k>2 && k<33);
            ksize = k;
            kmask=(1ULL<<2*ksize) - 1;
        } // InitK

    public:
        binary_kmer() {
            kvalue=0; // forward kmer
            kspan = 0;
        } // Kmer
        void inline AddChar(char nucleotide){
            int c = seq_nt4_table[(uint8_t)nucleotide];
            if( c>3 ) { // hit invalid nucleotide
                kspan = 0;
                kvalue=0;
                return ;
            }
            kvalue = (kvalue << 2 | c) & kmask;
            kspan ++ ;
            if( kspan > ksize)  kspan = ksize ;;
        } // AddChar
        inline bool valid() const {
            return kspan>=ksize ;      /*check kmer length*/
        } // valid
        inline uint64_t khash64() const {
            assert(valid());
            return hash64(kvalue, kmask);
        } // khash64
    private:
        uint64_t kvalue;
        int      kspan;
}; // class binary_kmer

int      binary_kmer::ksize;
uint64_t binary_kmer::kmask;

struct hash_kmer {
    uint64_t hash_forward;
    bool valid;
    void init(const binary_kmer & bkmer){
        valid = bkmer.valid();
        if( ! valid ) return ;
        hash_forward = bkmer.khash64() ;
    }
    //
    // @Brief : chop sequence if length l into l-k+1 kmers.
    // @In    : seq and len
    // @Out   : buff
    // @Warn  : user must ensure the buff is large enough for l-k+1 kmers.
    //
    static void chop_kmers(const char * seq , int len, hash_kmer * buffer){
        binary_kmer bkmer;
        for( int i = 0 ; i<len ; i++ ){
            bkmer.AddChar(seq[i]);
            int start_pos = i-bkmer.ksize+1 ;
            if(start_pos >= 0){
                buffer[start_pos].init(bkmer);
            }
        }
    } //chop_kmers
}; // class hash_kmer



int strobemer::nkmer;
int strobemer::ksize;
int strobemer::wsize;
int strobemer::wmin;
int strobemer::span ;
int strobemer::kspan;
strobemer_type strobemer::type ;
uint64_t strobemer::kmask;

void strobemer::init(int n, int k , int w_min , int w_max, strobemer_type t){
    assert(n>1);
    assert(w_min>=k);
    assert(w_max>w_min);
    
    //assert(t == strobemer_type::minstrobe || t == strobemer_type::randstrobe||t==strobemer_type::hybridstrobe );
    nkmer = n;
    ksize = k;
    wsize = w_max;
    wmin = w_min ;
    span=k+(wmin-k)+(n-1)*wsize;
    kspan=wsize - wmin +1;
    type = t ;
    kmask=(1ULL<<2*ksize) - 1;
    
    binary_kmer::InitK(k);

    std::cerr<<"INFO: init strobemer with"
        <<" n="<<n
        <<" k="<<k
        <<" wmin="<<wmin
        <<" wmax="<<wsize;
    if( t == strobemer_type::minstrobe ) 
        std::cerr<<" t="<<"minstrobe";
    else if ( t == strobemer_type::randstrobe ) 
        std::cerr<<" t="<<"randstrobe";
    else 
        std::cerr<<" t="<<"hybridstrobe";
    std::cerr<<" span="<<span
        <<" kspan="<<kspan
        <<std::endl;
}

void strobemer::chop_strobemer(const char * seq,int len, strobemer * buff){
    if( type == strobemer_type::minstrobe ) {
        chop_minstrobe(seq,len,buff);
    } else if ( type == strobemer_type::randstrobe) {
        chop_randstrobe(seq,len,buff);
    } else if (  type == strobemer_type::hybridstrobe ){
        chop_hybridstrobe(seq,len,buff);
    } else {
        assert(0);
    }
}

void strobemer::chop_strobemer_byKmer(const char *seq,int len,strobemer *buff,int &validLength){
    if( type == strobemer_type::minstrobe ) {
        chop_minstrobe_byKmer(seq,len,buff,validLength);
    } else if ( type == strobemer_type::randstrobe) {
        chop_randstrobe_byKmer(seq,len,buff,validLength);
    } else if (  type == strobemer_type::hybridstrobe ){
        chop_hybridstrobe_byKmer(seq,len,buff,validLength);
    } else {
        assert(0);
    }
}

static inline void make_seq_to_kmer(const char * seq, int len, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, int k, uint64_t kmask) {

    //unsigned int hash_count = 0;
    int l;
    int i;
    uint64_t x = 0;
    for (int i = l = 0; i < len; i++) {
        int c = seq_nt4_table[(uint8_t) seq[i]];
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & kmask;                 
            if (++l >= k) { 
                uint64_t hash_k = hash64(x, kmask);
                string_hashes.emplace_back(hash_k);
                pos_to_seq_choord.emplace_back( i - k + 1);
                //hash_count ++;
            }
        } else {
            l = 0, x = 0; // if there is an "N", restart
        }
    }
};

static inline void make_seq_to_kmer(const char * seq, int len, align_uint64_t *string_hashes, std::vector<unsigned int> &pos_to_seq_choord, int k, uint64_t kmask,int &size) {

    //unsigned int hash_count = 0;
    int l;
    int i;
    uint64_t x = 0;
    size=0;
    for (int i = l = 0; i < len; i++) {
        int c = seq_nt4_table[(uint8_t) seq[i]];
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & kmask;                 
            if (++l >= k) { 
                uint64_t hash_k = hash64(x, kmask);
                string_hashes[size++].value=hash_k;
                //string_hashes.emplace_back(hash_k);
                pos_to_seq_choord.emplace_back( i - k + 1);
                //hash_count ++;
            }
        } else {
            l = 0, x = 0; // if there is an "N", restart
        }
    }
  
};


void strobemer::chop_minstrobe(const char * seq,int len, strobemer * buff){
    // sanity check ...
    assert(seq!=nullptr);
    assert(len>=strobemer::span);
    
    // alloc buffers
    char * rc_seq = new char[len];
    hash_kmer* hkmer_buffer = new hash_kmer[len-binary_kmer::ksize+1];
    // sanity check ...
    assert(rc_seq!=nullptr);
    assert(hkmer_buffer!=nullptr);
    // fill buffers
    reverse_complete(seq,len,rc_seq);
    hash_kmer::chop_kmers(seq,len,hkmer_buffer);
    // clean buffer
    for( int i=0 ; i<=len-strobemer::span; i++ ){
        buff[i].valid = false;
    }

    for( int i=0 ; i<=len-strobemer::span; i++ ){
        // [i,w1_s) [w1_s , w1_e) [w1_e,w2_e)
        bool valid = true;
        // check the validation of span area
        for( int j = i ; j < i+strobemer::kspan; j++ ){
            if( ! hkmer_buffer[j].valid ) {
                i = j ; // iterator jump all infected area.
                valid =false ;
                break;
            }
        }
        if( ! valid ) continue ;
        buff[i].valid = true;
        //
        // -----------------------------------------------------------------------
        //                                     k            w       ....
        // 1. construct forward strobemer in [i,w1_s) [w1_s , w1_e) ....
        int p0 = i ;
        strncpy(buff[i].kmer_forward, &seq[p0],strobemer::ksize);
        // 1.1 find the ki minimizer in [wi_s , wi_e)
        int wi_s = i+strobemer::wmin;
        int wi_e = wi_s+strobemer::kspan ;
        int k_shift = strobemer::ksize;
        for( int ki=1; ki<strobemer::nkmer ; ki++){
            int p_next = -1;
            uint64_t temp_min ;
            uint64_t h_now;
            for( int j = wi_s ; j < wi_e ; j++ ){
                h_now = hkmer_buffer[j].hash_forward ; // minstrobe calculate hash value independantly
                if ( p_next == -1 || h_now < temp_min ) {
                    temp_min = h_now;
                    p_next = j ;
                }
            }
            strncpy(&buff[i].kmer_forward[k_shift], &seq[p_next],strobemer::ksize);
            wi_s=wi_s+strobemer::wsize;
            wi_e=wi_s+strobemer::kspan;
            k_shift += strobemer::ksize;
        }
    }
    // free memory
    delete [] hkmer_buffer;
    delete [] rc_seq;
}


void strobemer::chop_randstrobe(const char * seq,int len, strobemer * buff){
    // sanity check ...
    assert(seq!=nullptr);
    assert(len>=strobemer::span);


    // alloc buffers
    char * rc_seq = new char[len];
    hash_kmer* hkmer_buffer = new hash_kmer[len-binary_kmer::ksize+1];
    // sanity check ...
    assert(rc_seq!=nullptr);
    assert(hkmer_buffer!=nullptr);
    // fill buffers
    reverse_complete(seq,len,rc_seq);
    hash_kmer::chop_kmers(seq,len,hkmer_buffer);
    // clean buffer
    

    for( int i=0 ; i<=len-strobemer::span; i++ ){
        // [i,w1_s) [w1_s , w1_e) [w1_e,w2_e)
        bool valid = true;
        // check the validation of span area
        for( int j = i ; j < i+strobemer::kspan; j++ ){
            if( ! hkmer_buffer[j].valid ) {
                i = j ; // iterator jump all infected area.
                valid =false ;
                break;
            }
        }
        if( ! valid ) continue ;
        buff[i].valid = true;
        //
        // -----------------------------------------------------------------------
        //                                     k            w       ....
        // 1. construct forward strobemer in [i,w1_s) [w1_s , w1_e) ....
        int p0 = i ;
        strncpy(buff[i].kmer_forward, &seq[p0],strobemer::ksize);
        uint64_t h_prev = hkmer_buffer[i].hash_forward ;
        // 1.1 find the ki minimizer in [wi_s , wi_e)
        int wi_s = i+strobemer::wmin;
        int wi_e = wi_s+strobemer::kspan ;
        int k_shift = strobemer::ksize;
        for( int ki=1; ki<strobemer::nkmer ; ki++){
            int p_next = -1;
            uint64_t temp_min ;
            uint64_t h_now;
            for( int j = wi_s ; j < wi_e ; j++ ){
                h_now = h_prev ^hkmer_buffer[j].hash_forward ;
                if ( p_next == -1 || h_now < temp_min ) {
                    temp_min = h_now;
                    p_next = j ;
                }
            }
            strncpy(&buff[i].kmer_forward[k_shift], &seq[p_next],strobemer::ksize);
            wi_s=wi_s+strobemer::wsize;
            wi_e=wi_s+strobemer::kspan;
            k_shift += strobemer::ksize;
            h_prev = h_now ;
        }
    }
    // free memory
    delete [] hkmer_buffer;
    delete [] rc_seq;
}

void strobemer::chop_hybridstrobe(const char * seq,int len, strobemer * buff){
    // sanity check ...
    assert(seq!=nullptr);
    assert(len>=strobemer::span);
    // alloc buffers
    char * rc_seq = new char[len];
    hash_kmer* hkmer_buffer = new hash_kmer[len-binary_kmer::ksize+1];
    // sanity check ...
    assert(rc_seq!=nullptr);
    assert(hkmer_buffer!=nullptr);
    // fill buffers
    reverse_complete(seq,len,rc_seq);
    hash_kmer::chop_kmers(seq,len,hkmer_buffer);
    // clean buffer
    for( int i=0 ; i<=len-strobemer::span; i++ ){
        buff[i].valid = false;
    }

    for( int i=0 ; i<=len-strobemer::span; i++ ){
        // [i,w1_s) [w1_s , w1_e) [w1_e,w2_e)
        bool valid = true;
        // check the validation of span area
        for( int j = i ; j < i+strobemer::kspan; j++ ){
            if( ! hkmer_buffer[j].valid ) {
                i = j ; // iterator jump all infected area.
                valid =false ;
                break;
            }
        }
        if( ! valid ) continue ;
        buff[i].valid = true;
        //
        // -----------------------------------------------------------------------
        //                                     k            w       ....
        // 1. construct forward strobemer in [i,w1_s) [w1_s , w1_e) ....
        int p0 = i ;
        strncpy(buff[i].kmer_forward, &seq[p0],strobemer::ksize);
        uint64_t h_prev = hkmer_buffer[i].hash_forward ;
        // 1.1 find the ki minimizer in [wi_s , wi_e)
        int wi_s = i+strobemer::wmin;
        int wi_e = wi_s+strobemer::kspan ;
        int k_shift = strobemer::ksize;
        for( int ki=1; ki<strobemer::nkmer ; ki++){
            int p_next = -1;
            uint64_t temp_min ;
            uint64_t h_now;
            int small_window_size = strobemer::kspan /3 ;
            int current_window = h_prev % 3 ;
            for( int j = wi_s + current_window * small_window_size; j < wi_s + current_window * small_window_size + small_window_size ; j++ ){
                h_now = hkmer_buffer[j].hash_forward ;
                if ( p_next == -1 || h_now < temp_min ) {
                    temp_min = h_now;
                    p_next = j ;
                }
            }
            strncpy(&buff[i].kmer_forward[k_shift], &seq[p_next],strobemer::ksize);
            wi_s=wi_s+strobemer::wsize;
            wi_e=wi_s+strobemer::kspan;
            k_shift += strobemer::ksize;
            h_prev = h_now ;
        }
    }
    // free memory
    delete [] hkmer_buffer;
    delete [] rc_seq;
}

void strobemer::chop_minstrobe_byKmer(const char * seq,int len, strobemer * buff,int &validLength){
    // sanity check ...
    assert(seq!=nullptr);
    assert(len>=strobemer::span);
    //std::cout<<"调用的是minstrobes"<<std::endl;
    std::vector<uint64_t> kmer_hashes;
    std::vector<unsigned int> position_ofSeq;
    make_seq_to_kmer(seq,len,kmer_hashes,position_ofSeq,strobemer::ksize,strobemer::kmask);
    validLength=kmer_hashes.size()-strobemer::span;

    for(int i=0;i<=validLength;i++){
        // 1. construct forward strobemer in [i,w1_s) [w1_s , w1_e) ....
        unsigned int p0 = position_ofSeq[i] ;
        strncpy(buff[i].kmer_forward, &seq[p0],strobemer::ksize);
        // 1.1 find the ki minimizer in [wi_s , wi_e)
        int wi_s = i+strobemer::wmin;
        int wi_e = wi_s+strobemer::kspan ;
        int k_shift = strobemer::ksize;
        for( int ki=1; ki<strobemer::nkmer ; ki++){
            unsigned int p_next = -1;
            uint64_t temp_min ;
            uint64_t h_now;
            for( int j = wi_s ; j < wi_e ; j++ ){
                h_now = kmer_hashes[j] ; // minstrobe calculate hash value independantly
                if ( p_next == -1 || h_now < temp_min ) {
                    temp_min = h_now;
                    p_next = position_ofSeq[j] ;
                }
            }
            strncpy(&buff[i].kmer_forward[k_shift], &seq[p_next],strobemer::ksize);
            wi_s=wi_s+strobemer::wsize;
            wi_e=wi_s+strobemer::kspan;
            k_shift += strobemer::ksize;
        }

    }

}

void strobemer::chop_randstrobe_byKmer(const char * seq,int len, strobemer * buff,int &validLength){
    // sanity check ...
    assert(seq!=nullptr);
    assert(len>=strobemer::span);
    //std::cout<<"调用的是randstrobes"<<std::endl;
    std::vector<uint64_t> kmer_hashes;//kmer对应的hash值
    std::vector<unsigned int> positon_ofSeq;//每个kmer对应在seq中的位置
    make_seq_to_kmer(seq,len,kmer_hashes,positon_ofSeq,strobemer::ksize,strobemer::kmask);
    validLength=kmer_hashes.size()-strobemer::span;
    
    for(int i=0;i<=validLength;i++){
        //      k            w       ....
        // 1. construct forward strobemer in [i,w1_s) [w1_s , w1_e) ....
        unsigned int p0 =positon_ofSeq[i] ;
        strncpy(buff[i].kmer_forward, &seq[p0],strobemer::ksize);
        uint64_t h_prev = kmer_hashes[i];
        // 1.1 find the ki minimizer in [wi_s , wi_e)
        int wi_s = i+strobemer::wmin;
        int wi_e = wi_s+strobemer::kspan ;
        int k_shift = strobemer::ksize;
        for( int ki=1; ki<strobemer::nkmer ; ki++){
            unsigned int p_next = -1;
            uint64_t temp_min ;
            uint64_t h_now;
            for( int j = wi_s ; j < wi_e ; j++ ){
                h_now = h_prev ^ kmer_hashes[j];
                if ( p_next == -1 || h_now < temp_min ) {
                    temp_min = h_now;
                    p_next = positon_ofSeq[j] ;
                }
            }
            strncpy(&buff[i].kmer_forward[k_shift], &seq[p_next],strobemer::ksize);
            wi_s=wi_s+strobemer::wsize;
            wi_e=wi_s+strobemer::kspan;
            k_shift += strobemer::ksize;
            h_prev = h_now ;
        
        }
        
    }
}


void strobemer::chop_hybridstrobe_byKmer(const char * seq,int len, strobemer * buff,int &validLength){
    // sanity check ...
    assert(seq!=nullptr);
    assert(len>=strobemer::span);
    //std::cout<<"调用的是hybridstrobes"<<std::endl;
    std::vector<uint64_t> kmer_hashes;
    std::vector<unsigned int> position_ofSeq;
    make_seq_to_kmer(seq,len,kmer_hashes,position_ofSeq,strobemer::ksize,strobemer::kmask);
    validLength=kmer_hashes.size()-strobemer::span;

    for(int i=0;i<=validLength;i++){
        // 1. construct forward strobemer in [i,w1_s) [w1_s , w1_e) ....
        unsigned int p0 = i ;
        strncpy(buff[i].kmer_forward, &seq[p0],strobemer::ksize);
        uint64_t h_prev = kmer_hashes[i];
        // 1.1 find the ki minimizer in [wi_s , wi_e)
        int wi_s = i+strobemer::wmin;
        int wi_e = wi_s+strobemer::kspan ;
        int k_shift = strobemer::ksize;
        for( int ki=1; ki<strobemer::nkmer ; ki++){
            unsigned int p_next = -1;
            uint64_t temp_min ;
            uint64_t h_now;
            int small_window_size = strobemer::kspan /3 ;
            int current_window = h_prev % 3 ;
            for( int j = wi_s + current_window * small_window_size; j < wi_s + current_window * small_window_size + small_window_size ; j++ ){
                h_now = kmer_hashes[j] ;
                if ( p_next == -1 || h_now < temp_min ) {
                    temp_min = h_now;
                    p_next = position_ofSeq[j] ;
                }
            }
            strncpy(&buff[i].kmer_forward[k_shift], &seq[p_next],strobemer::ksize);
            wi_s=wi_s+strobemer::wsize;
            wi_e=wi_s+strobemer::kspan;
            k_shift += strobemer::ksize;
            h_prev = h_now ;
        }

    }

};

// void strobemer::chop_randstrobemers(Reference &data,consumerBuff &consumerbuff){
//     // sanity check ...
//     const char *seq=data.seq.c_str();
//     int len=data.length;

//     assert(seq!=nullptr);
//     assert(len>=strobemer::span);
//     consumerbuff.reset(len);
//     make_seq_to_kmer(seq,len,consumerbuff.kmer_hashes,consumerbuff.position_ofSeq,strobemer::ksize,strobemer::kmask);
//     int validLength=consumerbuff.kmer_hashes.size()-strobemer::span;
    
    
//     for(int i=0;i<=validLength;i++){
       
//         unsigned int p0 =consumerbuff.position_ofSeq[i] ;
//         strncpy(consumerbuff.strobemerBuff, &seq[p0],strobemer::ksize);
        
//         uint64_t h_prev = consumerbuff.kmer_hashes[i];
//         // 1.1 find the ki minimizer in [wi_s , wi_e)
//         int wi_s = i+strobemer::wmin;
//         int wi_e = wi_s+strobemer::kspan ;
//         int k_shift = strobemer::ksize;
        
//         for( int ki=1; ki<strobemer::nkmer ; ki++){
//             unsigned int p_next = -1;
//             uint64_t temp_min ;
//             uint64_t h_now;
//             for( int j = wi_s ; j < wi_e ; j++ ){
//                 h_now = h_prev ^ consumerbuff.kmer_hashes[j];
//                 if ( p_next == -1 || h_now < temp_min ) {
//                     temp_min = h_now;
//                     p_next = consumerbuff.position_ofSeq[j] ;
//                 }
//             }
            
            
//             strncpy(&consumerbuff.strobemerBuff[k_shift], &seq[p_next],strobemer::ksize);
//             //std::cout<<strobemerBuff<<std::endl;
//             wi_s=wi_s+strobemer::wsize;
//             wi_e=wi_s+strobemer::kspan;
//             k_shift += strobemer::ksize;
//             h_prev = h_now ;
        
//         }
           
        
//     }
    
// }

void strobemer::chop_randstrobemers(Reference &data,consumerBuff &consumerbuff){
    // sanity check ...
    const char *seq=data.seq.c_str();
    int len=data.length;

    assert(seq!=nullptr);
    assert(len>=strobemer::span);
    consumerbuff.reset(len);
    int size;
    make_seq_to_kmer(seq,len,consumerbuff.kmer_hashes,consumerbuff.position_ofSeq,strobemer::ksize,strobemer::kmask,size);
    int validLength=size-strobemer::span;
    consumerbuff.reserve(validLength);
    
    
    for(int i=0;i<=validLength;i++){
        
        int index=i*consumerbuff.buffColsize;
        unsigned int p0 =consumerbuff.position_ofSeq[i] ;
        //consumerbuff.positionSet[index++]=p0;
        //index++;
        //strncpy(consumerbuff.strobemerBuff, &seq[p0],strobemer::ksize);
        
        uint64_t h_prev = consumerbuff.kmer_hashes[i].value;
        // 1.1 find the ki minimizer in [wi_s , wi_e)
        int wi_s = i+strobemer::wmin;
        int wi_e = wi_s+strobemer::kspan ;
        int k_shift = strobemer::ksize;
        
        for( int ki=1; ki<strobemer::nkmer ; ki++){
            unsigned int p_next = -1;
            uint64_t temp_min ;
            uint64_t h_now;
            for( int j = wi_s ; j < wi_e ; j++ ){
                h_now = h_prev ^ consumerbuff.kmer_hashes[j].value;
                if ( p_next == -1 || h_now < temp_min ) {
                    temp_min = h_now;
                    p_next = consumerbuff.position_ofSeq[j] ;
                }
            }
            
            //consumerbuff.positionSet[index++]=p_next;
            
            //strncpy(&consumerbuff.strobemerBuff[k_shift], &seq[p_next],strobemer::ksize);
            //std::cout<<strobemerBuff<<std::endl;
            wi_s=wi_s+strobemer::wsize;
            wi_e=wi_s+strobemer::kspan;
            k_shift += strobemer::ksize;
            h_prev = h_now ;
        
        }
           
        
    }

}

void strobemer::chop_minstrobemers(Reference &data,consumerBuff &consumerbuff,bool write_flag){
    const char *seq=data.seq.c_str();
    int len=data.length;
    // sanity check ...
    assert(seq!=nullptr);
    assert(len>=strobemer::span);
    char *strobemerbuff;

    if(write_flag){
        std::string write_name=">"+data.name;
        const char *write_name_char=write_name.c_str();
        fwrite(write_name_char,strlen(write_name_char),1,consumerbuff.outputStream);
        fwrite("\r\n",2,1,consumerbuff.outputStream);
        strobemerbuff=new char[strobemer::ksize*strobemer::nkmer];
         
    }
    
    //std::cout<<"调用的是minstrobes"<<std::endl;

    consumerbuff.reset(len);
    int size;
    make_seq_to_kmer(seq,len,consumerbuff.kmer_hashes,consumerbuff.position_ofSeq,strobemer::ksize,strobemer::kmask,size);
    int validLength=size-strobemer::span;
    consumerbuff.reserve(validLength+1);

    for(int i=0;i<=validLength;i++){
        int Setindex=i*consumerbuff.buffColsize;
        unsigned int p0 =consumerbuff.position_ofSeq[i] ;
        consumerbuff.positionSet[Setindex++]=p0;
        if(write_flag){
            strncpy(strobemerbuff, &seq[p0],strobemer::ksize);
        }
        
        // 1.1 find the ki minimizer in [wi_s , wi_e)
        int wi_s = i+strobemer::wmin;
        int wi_e = wi_s+strobemer::kspan ;
        int k_shift = strobemer::ksize;
        for( int ki=1; ki<strobemer::nkmer ; ki++){
            unsigned int p_next = -1;
            uint64_t temp_min=UINT64_MAX ;
            //uint64_t h_now;
            for(int j = wi_s ; j < wi_e-(strobemer::kspan%8) ; j+=8){
              
                //std::cout<<"kmer_hashes[j]"<<kmer_hashes[j].value<<std::endl;
                __m512i data=_mm512_loadu_si512((__m512i *) &(consumerbuff.kmer_hashes[j].value));
               // _mm512_print_epi64(data);
                uint64_t minNum=_mm512_reduce_min_epu64(data);
                __m512i dMin=_mm512_set1_epi64(minNum);
                 //_mm512_print_epi64(dMin);
                __mmask8 cmpResutl=_mm512_cmpeq_epi64_mask(data,dMin);
                //std::cout<<"__mmask8为"<<_cvtmask8_u32(cmpResutl)<<std::endl;
                int index=check_index(cmpResutl)+j;
                //std::cout<<"index为"<<index<<std::endl;
                if(minNum<temp_min){
                    temp_min=minNum;
                    p_next=consumerbuff.position_ofSeq[index];
                }

            }
            for(int j= wi_e-(strobemer::kspan%8);j<wi_e;j++){
                if(consumerbuff.kmer_hashes[j].value<temp_min){
                    temp_min=consumerbuff.kmer_hashes[j].value;
                    p_next=consumerbuff.position_ofSeq[j];
                }

            }
            //resultposition[index++]=p_next;
            //strncpy(&buff[i].kmer_forward[k_shift], &seq[p_next],strobemer::ksize);
            if(write_flag){
                strncpy(&strobemerbuff[k_shift], &seq[p_next],strobemer::ksize);  
            }
            consumerbuff.positionSet[Setindex++]=p_next;
            wi_s=wi_s+strobemer::wsize;
            wi_e=wi_s+strobemer::kspan;
            k_shift += strobemer::ksize;
            
            //std::cout<<"p_next是"<<p_next<<std::endl;
        }
        if(write_flag){
            fwrite(strobemerbuff,strobemer::nkmer*strobemer::ksize,1,consumerbuff.outputStream);
            fwrite("\r\n",2,1,consumerbuff.outputStream);
        }

        //buff.emplace_back(s);
        //delete []s;

    }

    if(write_flag){
        delete [] strobemerbuff;
    }

}



strobemer::strobemer() {
    kmer_forward = new char[nkmer*ksize];
    valid = false ;
}
strobemer::~strobemer(){
    delete [] kmer_forward;
}
strobemer::strobemer(const strobemer & other){
    kmer_forward = new char[nkmer*ksize];
    strncpy(kmer_forward,other.kmer_forward,nkmer*ksize);
}
strobemer & strobemer::operator =(const strobemer & other){
    if( &other == this) return *this;
    assert(kmer_forward!=nullptr);
    strncpy(kmer_forward,other.kmer_forward,nkmer*ksize);
    return *this;
}
void strobemer::print() const{
    for(int i = 0 ; i <nkmer*ksize; i++)
        std::cout<<kmer_forward[i];
    std::cout<<'\n';
}

std::string strobemer::to_string() const {
    return std::string(kmer_forward,0,nkmer*ksize);
}
