#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/stream.h>
#include<bits/stdc++.h>

//#include <algorithm> 
//#include <cstdlib>
//#include <seqan/basic.h>
//#include <seqan/file.h>
//#include <seqan/modifier.h>
//#include <seqan/arg_parse.h>


using namespace std;
using namespace seqan;

pair<DnaString,int> get_best_conseq(vector<pair<DnaString,int>> conseqs) {
    pair<DnaString,int> best_conseq;
    int cur_min_score = -1;
    for(auto conseq : conseqs) {
        // lowest possible score, return right away
        if(conseq.second == 0) 
            return conseq;
        // for the first run, set cur_min_score to first conseq score
        if(cur_min_score == -1) cur_min_score = conseq.second;

        // look for the smallest conseq score
        if(conseq.second < cur_min_score) {
            cur_min_score = conseq.second;
            best_conseq = conseq;
        }
    }
    return best_conseq;
}

char getOneOrZero(int number) { //int to char function for the bitmap
    if(number == 0) return '0';
    if(number == 1) return '1';
    return -1;
}


// from https://www.geeksforgeeks.org/hamming-distance-two-strings/
int hammingDist(DnaString str1, DnaString str2) {
    int i = 0, count = 0;
    for(int x = 0; x < length(str1); x++) {
        if (str1[i] != str2[i])
            count++;
        i++;
    }
    return count;
}


int projection(int l, int k, int s, int m, int d, StringSet<DnaString> sequences){

    bool print_bucket_contents = false;

    int motif_length = l;
    //we iterate m times; we have m trials
    int seq_amount = length(sequences);
    int seq_length = length(sequences[0]);
    vector<pair<DnaString,int>> final_conseqs;
    for(int tr = 0; tr < m; tr++){ // tr stands for trials
                                /*RANDOM PROJECTIONS*/
        //for each l-mer in each of the t sequences 
        //to chose k random positions 
        //we first generate an array of length l with k ones and l-k zeros
        String<int> bitmapNumbers;
        String<char> bitmap;

        for (int i=0; i<k; ++i){
            bitmapNumbers += 1;
        } 
        for (int i=k; i<motif_length; ++i){
            bitmapNumbers += 0;
        } 
    
        srand(1); //set seed
        //shuffle the string of ints bitmapNumbers
        for (int i=(motif_length-1); i>0; --i){ 
            int j = rand()%i;         //chose a random position to swap
            int temp = bitmapNumbers[i]; //temporary variable
            bitmapNumbers[i] = bitmapNumbers[j];
            bitmapNumbers[j] = temp;
        }

        //we need the map in the form of a string of char
        for (int i=0; i<motif_length; i++){
            append(bitmap, getOneOrZero(bitmapNumbers[i]));
        }
        //create a map (dictionary) called buckets
        //the key will be the hash value of the hashed k-mer
        //the value will be a vector of pairs
        //in each pair we save the number of the sequence in the string set sequences 
        //and the position of the l-mer  
        map<int, vector<pair<int,int>>> buckets;
        
        //iterator for the string set
        typedef Iterator<StringSet<DnaString> >::Type TStringSetIterator;   

        //for each sequence in sequences
        for (TStringSetIterator it = begin(sequences); it != end(sequences); ++it){
            //use the bitmap on a generic shape
            Index<DnaString, IndexQGram<GenericShape> > index(*it); 
            stringToShape(indexShape(index), bitmap);

            //turn the iterator into an int so we can use it to save the sequence number
            int seq = distance(begin(sequences),it); 
            
            for (int i = 0; i < length(*it) - motif_length + 1; ++i){
                int hashV = seqan::hash(indexShape(index), begin(*it) + i);
                if(buckets.count(hashV) > 0) {
                    buckets[hashV].push_back(pair<int,int>(seq,i));
                } else {
                    vector<pair<int,int>> tmp_vec;
                    tmp_vec.push_back(pair<int,int>(seq,i));
                    buckets.insert( pair<int, vector<pair<int,int>>>(hashV, tmp_vec));
                }
            }
        }
        //print the content of buckets (seqNr and positon in seq) and its size
        if(print_bucket_contents) {
            for(auto elem : buckets) {
                cout << "buckets[" << elem.first << "] (size: "<< elem.second.size()<< ") => {" << endl;
                for(auto p : elem.second) {
                    //cout << "sequence: " << value(sequences,p.first) << endl;
                    cout << "    " << p.first << "," << p.second << "\n";
                    //cout << "starting pos of lmer: " << value(sequences,p.first)[p.second] << endl;
                } 
                cout << "}" << endl << endl;
            }
        }

                            /*MOTIF REFINEMENT*/

        vector <pair<DnaString,int>> final_kmer_conseqs;  

        //s is the bucket treshold
        //y is the size of a bucket
        //refine each bucket h with y=>s elements with EM algorithm
        for(auto elem : buckets){
            int y = elem.second.size(); //number of pairs = number of k-mers in the bucket
            if(y >= s){
            
                float P[4] = {0.25, 0.25, 0.25, 0.25}; //background probability distribution
                //initialize weight matrix Wh for prob of a base in the motif with 0-s
                float Wh[4][motif_length]; 
                std::fill(Wh[0], Wh[0] + 4 * motif_length, 0);
    
                for(auto p : elem.second) { //for each k-mer
                    int i = p.first; // seqNr
                    int j = p.second; // position of l-mer in sequence
                    //calculate the frequencies of the bases at every position u in the l-mers
                    //we access the lmer with value(sequences,i)
                    for (int u = 0; u < motif_length; u++)
                        Wh[ordValue(value(sequences,i)[u+j])][u]++;   
                } 

                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < motif_length; j++)
                        Wh[i][j] = Wh[i][j]/y + P[i]; //in percent + Laplace correction

                //cout << "for bucket with hashValue: " << elem.first << endl;
                /*for (int i = 0; i < 4; i++){
                    cout << "i is " << i << "\n";
                    for (int j = 0; j < motif_length; j++){
                        cout << "j is " << j << endl;
                        cout << "value in matrix is " << Wh[i][j] << endl;
                    }     
                }*/
                //-------------Wh is ready

	/*Wh[0][0] =  0.1; Wh[0][1] =  0.5; Wh[0][2] =  0.2;
	Wh[1][0] =  0.3; Wh[1][1] =  0.2; Wh[1][2] =  0.1;
	Wh[2][0] =  0.3; Wh[2][1] =  0.1; Wh[2][2] =  0.4;
	Wh[3][0] =  0.3; Wh[3][1] =  0.2; Wh[3][2] =  0.3;*/



                StringSet<DnaString> T; //multiset of l-mers
                //use the Wh* matrix on every sequence and save the l-mer with the highest 
                //likelihood ratio in T 
                int score_of_T = 0;

                //initialize a position matrix posM which shows the probability that the motif starts at
                //a certain position in the sequence (using the initial matrix Wh)
                float posM[seq_amount][seq_length-motif_length+1];
                /*
                for (TStringSetIterator it = begin(sequences); it != end(sequences); ++it){
                    int seq = distance(begin(sequences),it); //seqNr
                    //Index<DnaString, IndexEsa<> > index(*it);
                    float denominator = 0; //Nenner
                    for (int i = 0; i < length(*it) - motif_length + 1; ++i){
                        float numerator = 1;
                        for (int u = 0; u < motif_length; u++){
                            numerator *= Wh[ordValue((*it)[u+i])][u];    
                        }
                        
                        posM[seq][i] = numerator; //Zaehler
                        denominator += numerator;
                    }
                    for(int i = 0; i < length(*it) - motif_length + 1; ++i){
                        posM[seq][i] = posM[seq][i]/denominator;
                        cout << posM[seq][i] << "\t";
                    }
                    cout << endl;
                }
                cout << endl << endl;
                */
                
                //Refine weight matrix Wh until it converges
                float ref_Wh_tmp[4][motif_length];
                float ref_Wh[4][motif_length];
                float sums[motif_length];
                
                for(int refine_iter = 0; refine_iter < 100; refine_iter++){

                    for (TStringSetIterator it = begin(sequences); it != end(sequences); ++it){
                        int seq = distance(begin(sequences),it); //seqNr
                        //Index<DnaString, IndexEsa<> > index(*it);
                        float denominator = 0; //Nenner
                        for (int i = 0; i < length(*it) - motif_length + 1; ++i){
                            float numerator = 1;
                            for (int u = 0; u < motif_length; u++){
                                numerator *= Wh[ordValue((*it)[u+i])][u];    
                            }
                            
                            posM[seq][i] = numerator; //Zaehler
                            denominator += numerator;
                        }
                        for(int i = 0; i < length(*it) - motif_length + 1; ++i){
                            posM[seq][i] = posM[seq][i]/denominator;
                            //cout << posM[seq][i] << "\t";
                        }
                        //cout << endl;
                    }
                    //cout << endl << endl;

                    // make Wh the current ref_Wh_tmp:                    
                    for(int i = 0; i<4; i++)
                        for(int j = 0; j<motif_length; j++)
                            ref_Wh_tmp[i][j] = Wh[i][j];
                            
                    std::fill(ref_Wh[0], ref_Wh[0] + 4 * motif_length, 0);
                    std::fill(sums, sums+4, 0);

                    /*
                    start is the position in the motif whose probability we are researching P’ ordValue(), start
                    i ist the number of the sequence
                    pos is the position in the window
                    j is the position in the sequence
                    */
                    for(int start = 0; start < motif_length; start++)
                        for(int i=0; i < seq_amount; i++)
                            for (int j=start, pos=0; (pos < seq_length - motif_length + 1) ; j++, pos++)
                                ref_Wh_tmp[ordValue(sequences[i][j])][start] += posM[i][pos]; //nominator


                    for(int i = 0; i<4; i++) //iterating through Wh
                        for(int j = 0; j < motif_length; j++) {
                            //ref_Wh_tmp[i][j] += Wh[i][j];
                            sums[j] += ref_Wh_tmp[i][j]; //for the denominator
                            ref_Wh[i][j] = ref_Wh_tmp[i][j]; //for better overview
                        }

                    for(int i = 0; i<4; i++) {
                        for(int j = 0; j < motif_length; j++) {
                            Wh[i][j] = ref_Wh[i][j] / sums[j];
                            //cout << "  \t" << Wh[i][j];
                        }
                        //cout << endl;
                    }    
                    //cout << endl;  
                } // End of refinement
                
                /*for(int i = 0; i<4; i++) {
                    for(int j = 0; j < motif_length; j++) {
                        cout << "  \t" << Wh[i][j];
                    }
                    cout << endl;
                } */   


                for(int i = 0; i < seq_amount; i++) {
                    DnaString current_lmer;
                    for(int j = 0; j < seq_length - motif_length + 1; j++) {
                        if(posM[i][j] == 1) {
                            for(int k = 0; k<motif_length; k++)
                                current_lmer += sequences[i][j+k];
                            break;
                        }
                    }
                    appendValue(T, current_lmer);
                    //cout << "current_lmer: " << current_lmer << endl;
                } // T is filled with this buckets lmers now


                DnaString Ct;
                for(int j=0; j<motif_length;j++) {
                    int scores[4] = {0};
                    int max_score = 0;
                    char character = 'A';
                    for(int i=0; i<length(T); i++) {
                        scores[ordValue(T[i][j])]++;
                        if(max_score < scores[ordValue(T[i][j])]) {
                            max_score = scores[ordValue(T[i][j])];
                            character = T[i][j];
                        }
                    }
                    Ct += character;
                }
                //cout << "Ct: " << Ct << endl;

                for(auto lmer : T)
                    if(hammingDist(Ct, lmer) > d) score_of_T++;

                final_kmer_conseqs.push_back(pair<DnaString,int>(Ct,score_of_T));
                //cout << "score: " << score_of_T << endl;

                //Ct is the consensus of T
                //s(T) is the number of elements of T whose Hamming distance to Ct exceeds d

            } // End of (for each bucket in threshold) (this is where stuff happens)
        } // End of (for each bucket) (nothing really happens between above and here)
        
        pair<DnaString,int> best_conseq_of_kmer = get_best_conseq(final_kmer_conseqs);
        // best_conseq_of_kmer (or an euqally good one) is found and saved
        final_conseqs.push_back(best_conseq_of_kmer);

    } // End of trials
    //choose consensus pattern of best bucket 
    //the best bucket is the one with the smallest s(T)

    
    pair<DnaString,int> best_conseq = get_best_conseq(final_conseqs);
    cout << "searched consensus sequence: [" << best_conseq.first << "] with a score of " << best_conseq.second << endl;

    return 0;
}

int main(){

    //iterator for the string set
    typedef Iterator<StringSet<DnaString> >::Type TStringSetIterator;
    int l = 8;
    int motif_length = l;
    int d = 2; //allowed mutations
    int k = 5; //the number of positions to be projected

    //choose a value for s (bucket treshold)
    int s = 2;
    //calculate the optimal number of trials m
    int m = 10;
    
    //create example data
    //t = 3 sequences of length n = 12, where the mutated motif m' occurs
    //save them in an a string set called sequences
    
    StringSet<DnaString> sequences;
    DnaString str0 = "ACAGTTGCACA";
    appendValue(sequences, str0);
    DnaString str1 = "AGGCAGTGGCA";
    appendValue(sequences, str1);
    DnaString str2 = "TCAATTGCATC";
    appendValue(sequences, str2);

    //call the function
    projection(l, k, s, m, d, sequences);

    
    //std::cout << typeid(a).name()  << std::endl;
    
}
