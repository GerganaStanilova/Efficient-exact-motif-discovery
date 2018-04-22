#include <iostream>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/file.h>
#include <seqan/modifier.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
using namespace seqan;

unsigned h(){
    unsigned hx;
    return hx;
}

int projection(){

    //for i = 1 to m (we iterate m times; we have m trials)
        //chose k random positions and save them in an array called kPositions

        //create a map (dictionary) called buckets with 4^k key-value pairs
        //the key will be the result hx from the hash function h
        //the value will be the k-mer x which was given as input in the hash function
        //h(x) = hx

                            /*RANDOM PROJECTIONS*/

        //for each l-mer in each of the t sequences 
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
    String<Dna5> motif;
    motif = "AAAAAA";
    unsigned l = length(motif);
    unsigned d = 2; //allowed mutations
    unsigned k = l-d-1; //the number of positions to be projected

    //create example data
    //t = 30 sequences of length n = 13, where the mutated motif m' occurs
    //save them in an array called sequences

    //choose a value for s (bucket treshold)

    //calculate the optimal number of trials m

    //working only with the motif for now
    String<Dna5> x;
    x = motif;
    
    //projection(unsigned k, unsigned s, unsigned m, array sequences)

    
    String<int> kPositions;
    srand (1); //set seed
    for (unsigned i=0; i < k; i++){
        kPositions +=  rand() % (l-1);
    }

    for (int i=0; i<length(kPositions); i++){
        std::cout << kPositions[i] << std::endl;
    }
    
}
