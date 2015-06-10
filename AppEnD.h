/*******     File:    AppEnD.h
 *******   Author:    Joshua Welch
 *******     Date:    6/9/15
 *******  Purpose:    This program takes indexed BAM input and identifies
 *******  untemplated 3' RNA additions.
 *******/

#include <iostream>
#include <iomanip>
#include <set>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
using namespace BamTools;
using namespace std;

class Parameters
{
public:
	bool read1_start;			//Look for tails at beginning of read 1?
	bool read1_end;				//Look for tails at end of read 1?	
	bool read2_start;			//Look for tails at beginning of read 2?
	bool read2_end;				//Look for tails at end of read 2?
	bool single_start;			//Look for tails at beginning of non-paired-end read?
	bool single_end;			//Look for tails at end of non-paired-end read?
	bool exact_linker_match;	//Require exact linker sequence?
	bool output_reads;			//Output read sequences?
	bool output_read_ids;		//Output read IDs?
	bool output_tail_seqs;		//Output tail sequences?
	bool output_tail_comp;		//Output tail nucleotide compositions?
	bool output_linker;			//Output identified linker sequence?
	int min_linker_match;		//Number of nucleotides of linker required (used only if exact linker match required)
	int min_tail_len;			//Minimum tail length (>= 0)
	int max_tail_len;			//Maximum tail length (0 <= max_tail_len <= read length)
	string linker;				//Sequence of linker marking tail position
	string input_type;			//Sequencing file format ("SAM", "BAM", or "FASTQ")
	string input_file;			//Sequencing input file
	string gene_file;			//Gene/transcript annotation file in GTF format
	Parameters();
	Parameters(string param_file);
};

class Tail
{
public:
	string read;
	string read_id;
	int quality;
	string sequence;
	string linker;
	string chr;
	int start;
	int length;
	double A_pct;
	double C_pct;
	double G_pct;
	double T_pct;
	bool misprime;
	vector<int> mRNA_pos;
	Tail();
	Tail(const Parameters & params, const BamAlignment & al, map<int,string> & chr_ids);
	friend ostream & operator << (ostream & out, const Tail & tail);
private:
	Parameters params;
	void get_tail(const BamAlignment & al);
	string get_clip(const BamAlignment & al);
	void get_tail_start(const BamAlignment & al, const string & clip);
	bool look_at_start(const BamAlignment & al) const;
	int needleman_wunsch(const string & linker, const string & clip, string & aligned_linker) const;
	string reverse_complement(const string & s) const;
	void compute_tail_comp();
};

void read_chr_ids (BamReader & bam_file, map<int,string> & chr_ids);
vector<string> split(const string& s, const string& delim, const bool keep_empty = false);
void read_regions(const string& regions_file, vector<pair<string,pair<int,int> > > & regions); 
bool no_indels(const BamAlignment & al);
void write_header(Parameters params, ostream & out);
bool read_of_interest(const BamAlignment & al, const Parameters & params);
void write_bedgraph(map<string, map<int, int> > & tail_counts, string file);
bool atob(const string & a);
void convert_to_indexed_bam(Parameters & params);
