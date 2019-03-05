#include <seqan/index.h>
#include <seqan/seq_io.h>


using namespace std;
using namespace seqan;

bool freePos(int& i, vector<int>& p) {
    for(int a : p)
        if(i == a)
            return false;
    return true;
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

map<string,string> args;
vector<string> keyless_args;
vector<string> valueless_args;

int main(int argc, char *argv[])
{
	for(int i=0; i < argc; i++) {
		if(argv[i][0] != '-') { keyless_args.push_back(string(argv[i])); continue; }
		string arg(argv[i]);
		size_t ftpos = arg.find("=");
		if(ftpos == -1) { valueless_args.push_back(arg.substr(1,ftpos-1)); continue; }
		string key = arg.substr(1,ftpos-1);
		string val = arg.substr(ftpos+1);
		args.insert(pair<string,string>(key,val));
	}

	bool help = false;

	for(string s : valueless_args) {
		cout << s << endl;
		if(s == "help")
			help = true;
	}

	if(help) {
		cout
			<< "Syntax"<<endl
			<< "./createDataset -key=value" << endl
			<< endl
			<< "Keys:" << endl
			<< "help                 // Displays this overview" << endl
			<< "t (default 20)       // The number of sequences" << endl
			<< "n (default 600)      // The length of a sequence" << endl
			<< "l (default 11)       // The length of the planted motif" << endl
			<< "d (default 2)        // The amount of mutations in each planted motif" << endl
			<< "r (default 10)       // The number of datasets to be created" << endl
			<< "prefix (default syn) // Prefix for the files" << endl
			<< endl
			<< "If a key is left out, the default value is taken." << endl;
		return 0;
	}

	//for(pair<string,string> p : args)	cout << p.first << ":" << p.second << endl;

    //t sequences each of length n
    int t = (args.find("t") == args.end()) ? 20 : stoi(args["t"]);
    int n = (args.find("n") == args.end()) ? 600 : stoi(args["n"]);
    //smaller dataset
    int l = (args.find("l") == args.end()) ? 11 : stoi(args["l"]);
    int d = (args.find("d") == args.end()) ? 2 : stoi(args["d"]);

    int rounds = (args.find("r") == args.end()) ? 10 : stoi(args["r"]);

    int motif_length = l;

    string setPrefix = (args.find("prefix") == args.end()) ? "syn" : args["prefix"];
    string setSuffix = to_string(motif_length)+"_"+to_string(d);

    cout << "t is " << t << endl << "n is " << n << endl << "l is " << l << endl << "d is " << d << endl << "prefix is " << setPrefix << endl;
    cout << "r is " << rounds << endl;

    //return -1;
    //for the sake of reproducibility
    std::mt19937 rng;

    //from https://stackoverflow.com/questions/7560114/random-number-c-in-some-range
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator

    //randomly chosing at how many positions the planted motif will have mutations
    std::uniform_int_distribution<> mut_num_distr(0, d); //how many mutations
    std::uniform_int_distribution<> mut_pos_distr(0, motif_length-1); //at which position
    std::uniform_int_distribution<> motif_pos_distr(0, n - motif_length); //at which position
    std::uniform_int_distribution<> dna(0, 3); //Dna nucleotide

    for(int i=1; i<=rounds; i++) {
		DnaString planted_motif;
		for (int i = 0; i < motif_length; ++i)
			appendValue(planted_motif, Dna(dna(eng)));
		//randomly chose a planted motif
		cout << "planted_motif " << planted_motif << endl;

		//save the planted motif in a file
		ofstream planted_motif_File;
		ostringstream ossP;
		ossP << setPrefix << "_planted_motif_" << setSuffix << "_" << i << ".csv";
		string planted_motif_FileName = ossP.str();
		planted_motif_File.open (planted_motif_FileName);
		planted_motif_File << "seq\tpos\tmotif\tmutations\tmut_pos" << endl;
		planted_motif_File << "X\tX\t" << planted_motif << "\t0\t";


		//save the sequences in a file
		ofstream fastaFile;
		ostringstream oss;
		oss << setPrefix << "_synthetic_data_" << setSuffix << "_" << i << ".fasta";
		string fastaFileName = oss.str();
		fastaFile.open (fastaFileName);

		//for each of the t sequences make a permutation of the planted motif
		for (int seqNum = 0; seqNum < t; seqNum++){
			vector<int> mutation_positions;
			DnaString planted_motif_in_seq;
			planted_motif_in_seq = planted_motif;
			//randomly chose the number of positions at which the planted motif will have mutations
			//int number_of_mut = mut_num_distr(eng);
			int number_of_mut = d;
			if(number_of_mut > 0){
				//for each mutation
				for(int i=0; i<number_of_mut; ++i){
					//randomly chose the positions at which the mutation will occur
					int mutation_pos = mut_pos_distr(eng);
					while(!freePos(mutation_pos, mutation_positions)) {
						mutation_pos = mut_pos_distr(eng);
					};
					mutation_positions.push_back(mutation_pos);
					DnaString r = Dna(dna(eng));
					while(planted_motif_in_seq[mutation_pos] == r)
						r = Dna(dna(eng));
					planted_motif_in_seq[mutation_pos] = r[0]; //the mutation
				}
			}
			//randomly chose the positions at which the motif will be in the sequence
			int pos_of_motif = motif_pos_distr(eng);
			string sequence;
			for (int j = 0; j < n; j++){
				if(j!=0 && j%70==0) sequence += "\n";
				if(j >= pos_of_motif &&  j<=pos_of_motif+length(planted_motif)-1){
					sequence += planted_motif_in_seq[j-pos_of_motif];
				} else {
					appendValue(sequence, Dna(dna(eng)));
				}
			}

			//saving the sequences in a fasta file
			fastaFile << ">seq" << seqNum << "\n" << sequence << endl;

			planted_motif_File << endl;
			planted_motif_File << seqNum << "\t" << pos_of_motif << "\t";
			planted_motif_File << planted_motif_in_seq << "\t" <<  hammingDist(planted_motif, planted_motif_in_seq) << "\t";
			for(int i : mutation_positions)
				planted_motif_File << "|" << i;
		}

		//fastaFile << "Writing this to a file.\n";
		fastaFile.close();
		planted_motif_File.close();
    }
    return 0;
}

