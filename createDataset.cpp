#include <iostream>
#include <algorithm> 
#include <cstdlib>
#include <vector> 
#include <cmath>
#include <random>
#include <fstream>
#include <string>
#include <sstream>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/file.h>
#include <seqan/modifier.h>

using namespace std;
using namespace seqan;

int main()
{
    

    //t sequences each of length n
    int t = 20;
    int n = 600;
    
    //(15,4)
    //int motif_length = 15; //l
    //int d = 4;

    //(14,4)
    //int motif_length = 14; //l
    //int d = 4;
    
    //(16,5) 
    //int motif_length = 16; //l
    //int d = 5;

    //(18,6)
    //int motif_length = 18; //l
    //int d = 6;
    
    //smaller dataset
    int motif_length = 11; //l
    int d = 2;

    //for the sake of reproducibility
    std::mt19937 rng;

    DnaString planted_motif;
    for (int i = 0; i < motif_length; ++i)
            appendValue(planted_motif, Dna(rng() % 4));
    //randomly chose a planted motif
    cout << "planted_motif " << planted_motif << endl;

    //save the planted motif in a file
    ofstream planted_motif_File;
    ostringstream ossP;
    ossP << "planted_motif_" << motif_length << "_" << d << ".fasta";
    string planted_motif_FileName = ossP.str();
    planted_motif_File.open (planted_motif_FileName);
    planted_motif_File << planted_motif;
    planted_motif_File.close();

    //save the sequences in a file
    ofstream fastaFile;
    ostringstream oss;
    oss << "exampleDataset_" << motif_length << "_" << d << ".fasta";
    string fastaFileName = oss.str();
    fastaFile.open (fastaFileName);

    //from https://stackoverflow.com/questions/7560114/random-number-c-in-some-range
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator

    //randomly chosing at how many positions the planted motif will have mutations
    std::uniform_int_distribution<> mut_num_distr(0, d); //how many mutations
    std::uniform_int_distribution<> mut_pos_distr(0, motif_length-1); //at which position
    std::uniform_int_distribution<> motif_pos_distr(0, n - motif_length); //at which position
    
    //for each of the t sequences make a permutation of the planted motif
    for (int seqNum = 0; seqNum < t; seqNum++){
        DnaString planted_motif_in_seq;
        planted_motif_in_seq = planted_motif;
        //randomly chose the number of positions at which the planted motif will have mutations
        int number_of_mut = mut_num_distr(eng);
        if(number_of_mut > 0){
            //for each mutation
            for(int i=0; i<number_of_mut; ++i){
                //randomly chose the positions at which the mutation will occur
                int mutation_pos = mut_pos_distr(eng);
                planted_motif_in_seq[mutation_pos] = Dna(rng() % 4); //the mutation (problem: could replace with the same letter)
            }
        }
        //randomly chose the positions at which the motif will be in the sequence
        int pos_of_motif = motif_pos_distr(eng);
        DnaString sequence;
        for (int j = 0; j <= n; j++){
            if(j == pos_of_motif){
                sequence += planted_motif_in_seq;
            } else if(j > pos_of_motif &&(j < (pos_of_motif+length(planted_motif)))){
                continue;
            } else if (j < pos_of_motif || (j > (pos_of_motif+length(planted_motif)))){
                appendValue(sequence, Dna(rng() % 4));
            }
        }
        
        //saving the sequences in a fasta file
        ostringstream oss2;
        oss2 << ">seq" << seqNum << "\n";
        string seq_identifier = oss2.str();
        fastaFile << seq_identifier;

        fastaFile << sequence;

        ostringstream oss3;
        oss3 << "\n";
        string new_line = oss3.str();
        fastaFile << new_line;

        
    }
    
    //fastaFile << "Writing this to a file.\n";
    fastaFile.close();
    return 0;
}
