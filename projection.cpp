#include <algorithm>
#include <bits/stdc++.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/modifier.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>


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
        if(cur_min_score == -1) {
            cur_min_score = conseq.second;
            best_conseq = conseq;
        } 

        // look for the smallest conseq score
        if(conseq.second < cur_min_score) {
            cur_min_score = conseq.second;
            best_conseq = conseq;
        }
    }
    return best_conseq;
}

double rwp(double num, int precision) { //round with precision = rwp
    double factor = pow(10.0,(double) precision+1);
    return (int)(((num * factor)+5)/10)/(factor/10);
} 

char getOneOrZero(int number) { //1 if it's positive (for k), 0 if it's negative (for l-k)
    return (number > 0) ? '1' : '0';
}


// from https://www.geeksforgeeks.org/hamming-distance-two-strings/
int hammingDist(DnaString& str1, DnaString& str2) {
    int i = 0, count = 0;
    for(int x = 0; x < length(str1); x++) {
        if (str1[i] != str2[i])
            count++;
        i++;
    }
    return count;
}


int projection(int l, int k, int s, int m, int d, StringSet<DnaString> sequences){

    bool print_trial_number = false;
    bool print_bucket_avg_size = false;
    
    bool print_bucket_size = false;
    bool print_bitmap = true;
    bool print_bucket_contents = false;
    bool print_current_posM = false;
    bool print_first_posM = false;
    bool print_refined_posM = false;
    bool print_initial_Wh = false;
    bool print_first_refined_Wh = false;
    bool print_current_lmer_added_to_T = false;
    bool print_conseqs_in_bucket = false;
    bool print_current_hamm = false;
    bool print_final_kmer_conseqs = false;
    bool print_best_conseq_of_kmer = false;
    bool print_final_conseqs = false;
    bool print_best_conseq = false;

    int motif_length = l;
    int projection_length = k;
    int threshold = s;
    int number_of_trials = m;
    int number_of_mutations = d;
    int seq_quantity = length(sequences);
    int seq_length = length(sequences[0]);
    vector<pair<DnaString,int>> final_conseqs;
    //we iterate m times; we have m trials
    for(int tr = 1; tr <= m; tr++){ // tr stands for trials
        if(print_trial_number) cout << "trial number: " << tr << endl;
                                
                               
                               
                                /*RANDOM PROJECTIONS*/
        //for each l-mer in each of the t sequences 
        //randomly and uniformly distribute k 1s and l-k 0s
        String<char> bitmap;

        random_device rd;
        mt19937 mt(rd());
        uniform_real_distribution<double> distribution(0.0,1.0);

        for (int i=0, j=projection_length, k=projection_length - motif_length; i<motif_length; i++)
            append(bitmap, (distribution(mt) > 0.5 && j != 0 || k == 0) ? getOneOrZero(j--) : getOneOrZero(k++));

        // Bitmap generated 

        //%debug
        if(print_bitmap){
            cout << "the bitmap is: ";
            for (int i=0; i<motif_length; i++){
                cout << bitmap [i];
             }
            cout << endl;
        } 
        
        //create a map (dictionary) called buckets
        //the key will be the hash value of the hashed k-mer
        //the value will be a vector of pairs of ints
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
            int cur_seq = distance(begin(sequences),it); 
            
            for (int i = 0; i < length(*it) - motif_length + 1; ++i){
                int hashV = seqan::hash(indexShape(index), begin(*it) + i);
                if(buckets.count(hashV) > 0) { //if it already exist
                    buckets[hashV].push_back(pair<int,int>(cur_seq,i));
                } else { //if it doesn't
                    vector<pair<int,int>> tmp_vec;
                    tmp_vec.push_back(pair<int,int>(cur_seq,i));
                    buckets.insert( pair<int, vector<pair<int,int>>>(hashV, tmp_vec));
                }
            }
        }


        //%debug
        //print the content of buckets (seqNr and positon in seq) and its size
        if(print_bucket_contents) {
            for(auto elem : buckets) {
                cout << "buckets[" << elem.first << "] (size: "<< elem.second.size()<< ") => {" << endl;
                for(auto p : elem.second) {
                    cout << "sequence: " << value(sequences,p.first) << endl;
                    cout << "    " << p.first << "," << p.second << "\n";
                    cout << "starting pos of lmer: " << value(sequences,p.first)[p.second] << endl;
                } 
                cout << "}" << endl << endl;
            }
        }
        //%debug
        if(print_bucket_size) {
            for(auto elem : buckets) {
                cout << "size: "<< elem.second.size()<< endl;
            }
        }
        //%debug
        if(print_bucket_avg_size) {
            int num_of_buckets = 0;
            int bucket_avg_size = 0;
            int size_of_largest_bucket = 0;
            for(auto elem : buckets) {
                
                num_of_buckets++;
                bucket_avg_size += elem.second.size();
                int size_of_bucket = elem.second.size();
                if(size_of_bucket > size_of_largest_bucket){
                    size_of_largest_bucket = size_of_bucket;
                }

            }
            bucket_avg_size =  bucket_avg_size/num_of_buckets;
            cout << "the avg size of the buckets is: " << bucket_avg_size << endl;
            cout << "the size of the largest bucket is: " << size_of_largest_bucket << endl;
        }




                            /*MOTIF REFINEMENT*/

        vector <pair<DnaString,int>> final_kmer_conseqs;  

        //s is the bucket treshold
        //y is the size of a bucket
        //refine each bucket h with y=>s elements with EM algorithm
        for(auto elem : buckets){
            int y = elem.second.size(); //number of pairs = number of elements in the bucket
            if(y >= s){
                cout << "for this bucket" << endl;
                float P[4] = {0.25, 0.25, 0.25, 0.25}; //background probability distribution
                //initialize weight matrix Wh for prob of a base in the motif with 0-s
                float Wh[4][motif_length]; 
                std::fill(Wh[0], Wh[0] + 4 * motif_length, 0);
    
                for(auto p : elem.second) { //for each element in the bucket
                    int i = p.first; // seqNr
                    int j = p.second; // position of l-mer in sequence
                    //calculate the frequencies of the bases at every position u in the l-mers
                    //ordValue turns a nt into an int
                    for (int u = 0; u < motif_length; u++)
                        Wh[ordValue(sequences[i][u+j])][u]++;   
                } 

       
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < motif_length; j++)
                        Wh[i][j] = Wh[i][j]/y + P[i]; //in percent + Laplace correction
               
                //-------------Wh is ready

                //%debug
                if(print_initial_Wh) {
                    cout << "\ninitial Wh:\n";
                    for (int i = 0; i < 4; i++) {
                        for (int j = 0; j < motif_length; j++) {
                            cout << "\t" << Wh[i][j];
                        }
                        cout << endl;
                    }
                }


                StringSet<DnaString> T; //multiset of l-mers
                //the score of T is the number of l-mers whose hamming distance to the consensus sequence is larger than d
                int score_of_T = 0;

                //initialize a position matrix posM which shows the probability that the motif starts at
                //a certain position in the sequence (using the initial matrix Wh)
                float posM[seq_quantity][seq_length-motif_length+1];
                std::fill(posM[0], posM[0] + seq_quantity * (seq_length-motif_length+1), 0);          
                
                
                float ref_Wh_tmp[4][motif_length];
                float ref_Wh[4][motif_length];
                float sums[motif_length];
                //Refine weight matrix W and posM until convergence
                for(int refine_iter = 0; refine_iter < 200; refine_iter++){
                    if(print_current_posM || (print_first_posM && refine_iter == 0))
                        cout << "\ncurrent posM:\n";
                    for (TStringSetIterator it = begin(sequences); it != end(sequences); ++it){
                        int seq = distance(begin(sequences),it); //seqNr
                        //Index<DnaString, IndexEsa<> > index(*it);
                        float denominator = 0; //denominator
                        for (int i = 0; i < seq_length - motif_length + 1; ++i){
                            float numerator = 1;
                            for (int u = 0; u < motif_length; u++){
                                numerator *= Wh[ordValue((*it)[u+i])][u];    
                            }
                            
                            posM[seq][i] = numerator; //numerator
                            denominator += numerator;
                        }
                        for(int i = 0; i < seq_length - motif_length + 1; ++i){
                            posM[seq][i] = posM[seq][i]/denominator;
                            if(print_current_posM || (print_first_posM && refine_iter == 0)) {
                                cout << posM[seq][i] << "\t";
                            } 
                        }
                        if(print_current_posM || (print_first_posM && refine_iter == 0)) cout << endl;
                    }
                    if(print_current_posM || (print_first_posM && refine_iter == 0)) cout << endl << endl;

                    // make Wh the current ref_Wh_tmp:                    
                    for(int i = 0; i<4; i++)
                        for(int j = 0; j<motif_length; j++)
                            ref_Wh_tmp[i][j] = Wh[i][j];
                            
                    std::fill(ref_Wh[0], ref_Wh[0] + 4 * motif_length, 0);
                    std::fill(sums, sums+motif_length, 0);

                    /*
                    start is the position in the motif whose probability we are researching P’ ordValue(), start
                    i ist the number of the sequence
                    pos is the position in the window
                    j is the position in the sequence
                    */
                    for(int start = 0; start < motif_length; start++)
                        for(int i=0; i < seq_quantity; i++)
                            for (int j=start, pos=0; (pos < seq_length - motif_length + 1) ; j++, pos++)
                                ref_Wh_tmp[ordValue(sequences[i][j])][start] += posM[i][pos]; //numerator


                    for(int i = 0; i<4; i++) //iterating through Wh
                        for(int j = 0; j < motif_length; j++) {
                            //ref_Wh_tmp[i][j] += Wh[i][j];
                            sums[j] += ref_Wh_tmp[i][j]; //for the denominator
                            ref_Wh[i][j] = ref_Wh_tmp[i][j]; //for better overview of numerator and denominator
                        }

                    for(int i = 0; i<4; i++) {
                        for(int j = 0; j < motif_length; j++) {
                            Wh[i][j] = ref_Wh[i][j] / sums[j];
                            if(isnan(Wh[i][j])) {
                                cout << "\nref_Wh["<<i<<"]["<<j<<"] = "<< ref_Wh[i][j] << " and sums["<<j<<"] is "<< sums[j] <<endl;
                                abort();
                            }
                            if(sums[j] == 0) abort();
                        }
                    }    
                    
                    if(refine_iter == 0 && print_first_refined_Wh) {
                        cout << "\nfirst refined Wh:\n";
                        for(int i = 0; i<4; i++) {
                            for(int j = 0; j < motif_length; j++) {
                                cout << "  \t" << Wh[i][j];
                            }
                            cout << endl;
                        }    
                    }
                } // End of refinement

                if(print_refined_posM) {
                    cout << "\nfinal posM:"<< seq_length <<"\n";
                    for(int i = 0; i < seq_quantity; i++) {
                        for(int j = 0; j < seq_length - motif_length + 1; j++) {
                            cout << "\t" << rwp(posM[i][j],2);
                        }
                        cout << endl;
                    }
                    cout << endl;
                }



                            /*CONSENSUS SEQUENCE*/

                //create the stringSet T of the l-mers
                //take an l-mer from each sequence using the posM
                for(int i = 0; i < seq_quantity; i++) {
                    DnaString current_lmer;
                    int row = 0;
                    int col = 0;
                    int local_max = 0;
                    for(int j = 0; j < seq_length - motif_length + 1; j++) {
                        if(posM[i][j] > local_max) {
                            local_max = posM[i][j];
                            row = i;
                            col = j;
                        }
                    }
                    for(int k = 0; k<motif_length; k++)
                        current_lmer += sequences[row][col+k];
                    appendValue(T, current_lmer);
                    if(print_current_lmer_added_to_T) cout << "current_lmer: " << current_lmer << endl;
                } // here T is filled with this bucket's lmers

                DnaString Ct; //consensus sequnce
                for(int j=0; j<motif_length;j++) { ;
                    int scores[4] = {0};
                    int max_score = 0;
                    char character; // = 'A'; // Platzhalter
                    for(int i=0; i<length(T); i++) {
                        scores[ordValue(T[i][j])]++;
                        if(max_score < scores[ordValue(T[i][j])]) {
                            max_score = scores[ordValue(T[i][j])];
                            character = T[i][j];
                        }
                    }
                    Ct += character;
                }

                if(print_conseqs_in_bucket) cout << "\nThis buckets conseq is " << Ct << endl;

                for(auto lmer : T) {
                    int hamm = hammingDist(Ct, lmer);
                    if(hamm > d) score_of_T++;
                    if(print_current_hamm) cout << "\nCurrent hamm between "<<Ct<<" and "<<lmer<< " is " << hamm << endl; 
                }

                final_kmer_conseqs.push_back(pair<DnaString,int>(Ct,score_of_T));
                
                if(print_final_kmer_conseqs) {
                    cout << "\nFinal kmer conseqs:\n";
                    for(auto conseq : final_kmer_conseqs) {
                        cout << "kmer: " << conseq.first << "\t|\tscore: " <<conseq.second<<endl;
                    }
                    cout << endl;
                }

                //Ct is the consensus of T
                //s(T) is the number of elements of T whose Hamming distance to Ct exceeds d

            } // End of (for each bucket in threshold)
        } // End of (for each bucket)
        
        pair<DnaString,int> best_conseq_of_kmer = get_best_conseq(final_kmer_conseqs);

        if(print_best_conseq_of_kmer) {
            cout << "\nbest conseq of kmer:\n";
            cout << "kmer: " << best_conseq_of_kmer.first << "\t|\tscore: " <<best_conseq_of_kmer.second<<endl;
        }
        // best_conseq_of_kmer (or an euqally good one) is found and saved
        final_conseqs.push_back(best_conseq_of_kmer);

    } // End of trials
    //choose consensus sequence of best bucket 
    //the best bucket is the one with the smallest s(T)

    if(print_final_conseqs) {
        cout << "\nFinal conseqs:\n";
        for(auto conseq : final_conseqs) {
            cout << "kmer: " << conseq.first << "\t|\tscore: " <<conseq.second<<endl;
        }
        cout << endl;
    }
    
    pair<DnaString,int> best_conseq = get_best_conseq(final_conseqs);
    if(print_best_conseq) {
        cout << "\nbest conseq:\n";
        cout << "kmer: " << best_conseq.first << "\t|\tscore: " <<best_conseq.second<<endl;
    }

    cout << "searched consensus sequence: [" << best_conseq.first << "] with a score of " << best_conseq.second << endl;

    return 0;
}

int main(int argc, char const ** argv){
    typedef Iterator<StringSet<DnaString> >::Type TStringSetIterator;    
    StringSet<DnaString> sequences;

   //get the fasta file via the console

    if (argc != 2)
    {
        std::cerr << "USAGE: build_fai FILE.fa\n";
        return 0;
    }

    FaiIndex faiIndex;
    if (!build(faiIndex, argv[1]))
    {
        std::cerr << "ERROR: Could not build FAI index for file " << argv[1] << ".\n";
        return 0;
    }

    CharString faiFilename = argv[1];
    append(faiFilename, ".fai");

    if (!save(faiIndex, toCString(faiFilename)))
    {
        std::cerr << "ERROR: Could not write the index to file!\n";
        return 0;
    }

    std::cout << "Index file " << faiFilename << " was successfully created.\n";

    //use the fai file to get the content of the fasta file
    CharString pathToFile = argv[1];

    if (!open(faiIndex, toCString(pathToFile)))
        std::cout << "ERROR: Could not load FAI index " << pathToFile << ".fai\n";

    
    
    unsigned num_of_seqs = numSeqs(faiIndex);
    cout << "Num of seqs: " << num_of_seqs << endl;
    for(unsigned idx = 0; idx < num_of_seqs; idx++){
        DnaString seq_in_file;
        readSequence(seq_in_file, faiIndex, idx);
        appendValue(sequences, seq_in_file); //save each sequence in the stringSet sequences
    }

    int l = 15;
    int motif_length = l;
    int d = 4; //allowed mutations
    int k = 7; //the number of positions to be projected

    //choose a value for s (bucket threshold)
    int s = 50;
    //calculate the optimal number of trials m
    int m = 1;
    
    //call the function
    projection(l, k, s, m, d, sequences);

}
