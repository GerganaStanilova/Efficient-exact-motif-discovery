#include <iostream>
#include <algorithm> 
#include <cstdlib>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/file.h>
#include <seqan/modifier.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>

using namespace std;
using namespace seqan;

unsigned h(){
    unsigned hx;
    return hx;
}

int projection(unsigned k, unsigned s, unsigned m, String<DnaString> sequences)
{
    unsigned l;
    //for i = 1 to m (we iterate m times; we have m trials)
        //chose k random positions and save them in an array called kPositions
        //first generate all possible positions of the motif
        String<int> lPositions;
        for (int i=0; i<l; ++i)
        {
            lPositions += i;
        } 
    
        srand(1); //set seed
        //shuffle the array lPositions
        for (int i=(l-1); i>0; --i){ //--i or i--???????????????????????????????????????????
            int j = rand()%i;         //chose a random position to swap
            int temp = lPositions[i]; //temporary variable
            lPositions[i] = lPositions[j];
            lPositions[j] = temp;
        }
    
        //get the first k elements of lPositions into kPositions
        String<int> kPositions;
        for (int i=0; i<k; ++i){
            kPositions += lPositions[i];
        } 
    
        //create a map (dictionary) called buckets with 4^k key-value pairs
        //the key will be the result hx from the hash function h
        //the value will be the k-mer x which was given as input in the hash function
        //h(x) = hx

                            /*RANDOM PROJECTIONS*/

        //for each l-mer in each of the t sequences 


        DnaString lmer;
        DnaString x;
        for (unsigned t = 0; t < length(sequences); ++t)
        {
            for(unsigned j = 0; j < length(sequences[t])-l; ++j)
            {
                
                lmer = infix(sequences[t], j, j+k);
                
                for (unsigned c = 0; c < k; ++c)
                {
                    x += lmer[kPositions[c]];
                    std::cout << x << std::endl;
                }
            } 
        }
            //create the k-mer x which resulted from taking the chars
            //at the k random positions and concatinating them

            //compute h(x)

            //store x in the map as value with key hx

                            /*MOTIF REFINEMENT*/
        //s  is the bucket treshold
        //y is the size of a bucket
        //for each bucket with y=>s elements
            //refine bucket with EM algorithm


        //choose consensus pattern of best bucket 


    return 0;
}

int main(){

    DnaString motif;
    motif = "AAAAAA";
    unsigned l = length(motif);
    unsigned d = 2; //allowed mutations
    unsigned k = l-d-1; //the number of positions to be projected

    //create example data
    //t = 30 sequences of length n = 13, where the mutated motif m' occurs
    //save them in an array called sequences

    String<DnaString> sequences;
    DnaString s0 = "GGAAAAAAGTCTG";
    DnaString s1 = "GGGAAAAAAGTTG";
    appendValue(sequences, s0);
    appendValue(sequences, s1);
    //choose a value for s (bucket treshold)

    //calculate the optimal number of trials m

    //working only with the motif for now
    DnaString x;
    x = motif;
    
    //projection(unsigned k, unsigned s, unsigned m, String<DnaString> sequences)
    projection(k, 1, 1, sequences);

    
    
    
    //concatinate the chars from the motif at the positions kPositions


    for (int i=0; i<k; i++){
        //std::cout << kPositions[i] << std::endl;
    }
    
}
