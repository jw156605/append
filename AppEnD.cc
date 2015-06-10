/*******     File:    AppEnD.cc
 *******   Author:    Joshua Welch
 *******     Date:    11/25/14
 *******  Purpose:    This program takes indexed BAM input and identifies
 *******  untemplated 3' RNA additions.
 *******/

#include "AppEnD.h"

int main (int argc, char* argv[])
{
	string param_file(argv[1]);

	cout << "Reading parameters..." << endl;
	Parameters params(param_file);

	cout << "Opening BAM file..." << endl;
	BamReader reader;
	if ( !reader.Open(params.input_file) ) {
	   cerr << "Could not open input BAM file." << endl;
	   cout << params.input_file << endl;
	   return EXIT_SUCCESS;
	}

	map<int,string> chr_ids;
	read_chr_ids(reader,chr_ids);
	// retrieve metadata from BAM files--required by BamWriter 
	const SamHeader header = reader.GetHeader();
	const RefVector references = reader.GetReferenceData();

	// attempt to open BamWriter
	BamWriter writer;
	if ( ! writer.Open("Aligned.out.tails.bam", header, references) ) {
		cout << "Could not open output BAM file" << endl;
		return EXIT_SUCCESS;  
	}

	ofstream fout ("tails.txt");
	ofstream misprime ("misprimes.txt");

	int tot_reads = 0;
	int num_reads = 0;
	int tot_misprime = 0;
	BamAlignment al;
	map <string, map<int, int> > tail_counts;

	cout << "Finding tails..." << endl;
	write_header(params,fout);
	while ( reader.GetNextAlignment(al) ) {
		++tot_reads;
		if (tot_reads % 1000000 == 0)
		{
			cout << tot_reads << " reads processed" << endl;
		}
		
		if (read_of_interest(al,params))
		  {
			  writer.SaveAlignment(al);
			  Tail tail(params,al,chr_ids);
			  fout << tail;
			  if (tail.chr != "")
			  {
					++num_reads;
					if (tail_counts.find(tail.chr) == tail_counts.end())
						tail_counts[tail.chr] = map<int, int>();
					if (tail_counts[tail.chr].find(tail.start) == tail_counts[tail.chr].end())
						tail_counts[tail.chr][tail.start] = 0;
					++tail_counts[tail.chr][tail.start];
			  }
			  if (tail.misprime)
			  {
			    ++tot_misprime;
			    misprime << tail.read_id << "\t" << tail.sequence << "\t" << chr_ids[al.RefID] << "\t" <<  al.Position << endl;
			  }
		  }
	}

	// close the reader & writer
	reader.Close();
	writer.Close();
	write_bedgraph(tail_counts,"tails.bedgraph");
    cout << "Found " << num_reads << " tails and " << tot_misprime <<  " misprimes in " << tot_reads << " reads."  << endl;
    return EXIT_SUCCESS;
}

void write_header(Parameters params, ostream & out)
{
	if (params.output_read_ids)
		out << "Read ID\t";
	if (params.output_reads)
		out << "Read\t";
	if (params.output_tail_seqs)
	{
		out << "Tail\t";
		out << "Tail Length\t";
	}
	if (params.output_linker)
	{
		out << "Linker\t";
	}
	if (params.output_tail_comp)
	{
		out	<< "% A\t% C\t% G\t% T\t";
	}
	out	<< "Read Alignment Quality\t"
		<< "Tail Chromosome\t"
		<< "Tail Start" << endl;
}

//Returns true if a given read alignment shows the correct user-specified soft clipping pattern
bool read_of_interest(const BamAlignment & al, const Parameters & params)
{
	if (!al.IsMapped()) return false;
	
	bool read1_start = params.read1_start;
	bool read1_end = params.read1_end;
	bool read2_start = params.read2_start;
	bool read2_end = params.read2_end;
	bool single_start = params.single_start;
	bool single_end = params.single_end;
	
	bool first = al.IsFirstMate();
	bool second = al.IsSecondMate();
	bool single = (!first && !second);

	if (!single && !al.IsProperPair()) //Use only reads where both reads in the pair align together
	  return false;

	vector<int> sizes, rpos, gpos;
	vector<CigarOp> cigar_ops;
	bool has_soft_clip = al.GetSoftClips(sizes,rpos,gpos);
	cigar_ops = al.CigarData;
	bool clip_at_start = (!al.IsReverseStrand() && cigar_ops[0].Type == 'S') || (al.IsReverseStrand() && cigar_ops[cigar_ops.size()-1].Type == 'S');
	bool clip_at_end = (al.IsReverseStrand() && cigar_ops[0].Type == 'S') || (!al.IsReverseStrand() && cigar_ops[cigar_ops.size()-1].Type == 'S');

	if (read1_start && first && has_soft_clip && clip_at_start)
	{
		return true;
	}
	if (read1_end && first && has_soft_clip && clip_at_end)
	{
		return true;
	}
	if (read2_start && second && has_soft_clip && clip_at_start)
	{
		return true;
	}
	if (read2_end && second && has_soft_clip && clip_at_end)
	{
		return true;
	}
	if (single_start && single && has_soft_clip && clip_at_start)
	{
		return true;
	}
	if (single_end && single && has_soft_clip && clip_at_end)
	{
		return true;
	}
	return false;
}

void write_bedgraph(map<string, map<int, int> > & tail_counts, string file)
{
	ofstream fout(file.c_str());
	fout << "track type=bedGraph name=\"Tail Counts\" graphType=bar\n";
	for (map<string, map<int, int> >::iterator i = tail_counts.begin(); i != tail_counts.end(); ++i)
	{
		for (map<int, int>::iterator j = tail_counts[i->first].begin(); j != tail_counts[i->first].end(); ++j)
		{
			fout << i->first << "\t" << (j->first)-1 << "\t" << j->first << "\t" << j->second << endl;
		}
	}
}

bool atob (const string & a)
{
	return (a == "true" || a == "1");
}

/******************************************
****** Parameters member functions ********
*******************************************/

Parameters::Parameters()
{
	read1_start = true; 
	read1_end = false; 
	read2_start = false;
	read2_end = false; 
	single_start = false;
	single_end = false;
	exact_linker_match = false;
	output_reads = true;
	output_read_ids = true;
	output_tail_seqs = true;
	output_tail_comp = true;
	min_linker_match = 0;
	min_tail_len = 0;
	max_tail_len = 0;
	linker = "ATCACCGACTGCCCATAGAGAGGTG";
	input_type = "BAM";
	input_file = "Aligned.out.sort.bam";
}

Parameters::Parameters(string param_file)
{
	ifstream fin(param_file.c_str());
	if (fin.fail())
	{
		cout << "Parameter file failed to open." << endl;
		return;
	}
	string line;
	vector<string> line_parts;

	getline(fin, line);
	int tab_pos = line.find("\t");
	read1_start = atob(line.substr(tab_pos+1,line.length()-tab_pos-1));

	getline(fin, line);
	tab_pos = line.find("\t");
	read1_end = atob(line.substr(tab_pos+1,line.length()-tab_pos-1));

	getline(fin, line);
	tab_pos = line.find("\t");
	read2_start = atob(line.substr(tab_pos+1,line.length()-tab_pos-1));
	
	getline(fin, line);
	tab_pos = line.find("\t");
	read2_end = atob(line.substr(tab_pos+1,line.length()-tab_pos-1));

	getline(fin, line);
	tab_pos = line.find("\t");
	single_start = atob(line.substr(tab_pos+1,line.length()-tab_pos-1));

	getline(fin, line);
	tab_pos = line.find("\t");
	single_end = atob(line.substr(tab_pos+1,line.length()-tab_pos-1));

	getline(fin, line);
	tab_pos = line.find("\t");
	exact_linker_match = atob(line.substr(tab_pos+1,line.length()-tab_pos-1));

	getline(fin, line);
	tab_pos = line.find("\t");
	output_reads = atob(line.substr(tab_pos+1,line.length()-tab_pos-1));

	getline(fin, line);
	tab_pos = line.find("\t");
	output_read_ids = atob(line.substr(tab_pos+1,line.length()-tab_pos-1));

	getline(fin, line);
	tab_pos = line.find("\t");
	output_tail_seqs = atob(line.substr(tab_pos+1,line.length()-tab_pos-1));

	getline(fin, line);
	tab_pos = line.find("\t");
	output_tail_comp = atob(line.substr(tab_pos+1,line.length()-tab_pos-1));

	getline(fin, line);
	tab_pos = line.find("\t");
	output_linker = atob(line.substr(tab_pos+1,line.length()-tab_pos-1));

	getline(fin, line);
	tab_pos = line.find("\t");
	min_linker_match = atoi(line.substr(tab_pos+1,line.length()-tab_pos-1).c_str());

	getline(fin, line);
	tab_pos = line.find("\t");
	min_tail_len = atoi(line.substr(tab_pos+1,line.length()-tab_pos-1).c_str());

	getline(fin, line);
	tab_pos = line.find("\t");
	max_tail_len = atoi(line.substr(tab_pos+1,line.length()-tab_pos-1).c_str());

	getline(fin, line);
	tab_pos = line.find("\t");
	linker = line.substr(tab_pos+1,line.length()-tab_pos-1);

	getline(fin, line);
	tab_pos = line.find("\t");
	input_type = line.substr(tab_pos+1,line.length()-tab_pos-1);

	getline(fin, line);
	tab_pos = line.find("\t");
	input_file = line.substr(tab_pos+1,line.length()-tab_pos-1);
}

/******************************************
********* Tail member functions ***********
*******************************************/

Tail::Tail()
{
}

Tail::Tail(const Parameters & params, const BamAlignment & al, map<int,string> & chr_ids)
{
	this->params = params;
	read = al.QueryBases;
	if (params.output_read_ids)
		read_id = al.Name;
	quality = al.MapQuality;
	chr = chr_ids[al.RefID];
	misprime = false;
	get_tail(al);
	if (!params.output_tail_seqs)
		sequence = "";
	if (params.output_reads)
	{
		read = reverse_complement(read);
	}
	else
	{
		read = "";
	}
	length = sequence.length();
	if (length < params.min_tail_len || length > params.max_tail_len)
		chr = "";
	if (params.output_tail_comp)
	{
		compute_tail_comp();
	}
	else
	{
		A_pct = C_pct = G_pct = T_pct = -1;
	}
}

void Tail::get_tail(const BamAlignment & al)
{
	string clip = get_clip(al);
	get_tail_start(al, clip);
	string linker = reverse_complement(params.linker);
	length = 0;
	
	if (!params.exact_linker_match)
	{
	  string aligned_linker;
		int linker_score = needleman_wunsch(linker,clip,aligned_linker);
		for (int i = aligned_linker.length()-1; i > 0; i--)
		{
			if (aligned_linker[i] == '-')
			{
				length++;
				linker_score++;
			}
			else
			{
				break;
			}
		}
	}
	else
	{
	  if (params.min_linker_match > 0){
		string exact_linker_stub = linker.substr(linker.length()-params.min_linker_match,params.min_linker_match);
		int linker_pos = clip.rfind(exact_linker_stub);
		if (linker_pos == -1) 
		{
			chr = "";
			if (read.rfind(exact_linker_stub)) //Couldn't find linker in clip, but found it in aligned portion of read
			  {
			    misprime = true;
			  }
			return;
		}
		//Example:
		//GGTGATAA
		//012345678
		//linker_pos = 2, clip.length() = 8, params.min_linker_match = 4, length = 2
		length = clip.length() - linker_pos - params.min_linker_match;
		sequence = clip.substr(0,length);
	  }
	  else
	  {
	    length = clip.length();
	    sequence = clip;
	  }
	}
	if (length == 0)
	{
		sequence = "";
	}
	else
	{
		sequence = reverse_complement(clip.substr(clip.length()-length, length));
	}
	this->linker = reverse_complement(clip.substr(0,clip.length()-length));
}

string Tail::get_clip(const BamAlignment & al)
{
	//Example (tail in lowercase, linker uppercase):
	//transcript+linker: 5' ...uuATC 3'
	//cDNA 1st   strand: 5' ...GATaa... 3'
	//cDNA 2nd   strand: 3' ...CTAtt... 5'
	//Tail at     start: 5' ...GATaa... 3'
	//Tail at       end: 5' ...ttATC... 3'
	//Bottom line: reverse complementing one gives the other

	vector<CigarOp> cigar_ops;
	cigar_ops = al.CigarData;
	read = al.QueryBases;
	bool look_start = look_at_start(al);
	int clip_length;

	//Clip is at beginning of read if we're looking for tails 
	//at the beginning and alignment is forward strand
	//or we're looking for tails at the end and alignment
	//is reverse strand.
	if ( look_start ) 
	{
	  clip_length = cigar_ops[0].Length;     
	}

	//If we're looking for tails at the beginning and alignment 
	//is reverse strand or we're looking for tails at the end 
	//and alignment is forward strand, take the complement
	// so that the soft clip is in the correct position.
	if ( !look_start )
	{
	  clip_length = cigar_ops[cigar_ops.size()-1].Length;
	  read = reverse_complement(read);
	}
	return read.substr(0,clip_length);
}

void Tail::get_tail_start(const BamAlignment & al, const string & clip)
{
	bool look_start = look_at_start(al);
	start = al.Position; //Note: BamTools uses 0-based coords

	//Need to take indels into account when computing tail start pos in this case:
	if ( !look_start )	
	{
		vector<CigarOp> cigar_ops;
		cigar_ops = al.CigarData;
		for (int i = 0; i < cigar_ops.size(); ++i)
		{
			if (cigar_ops[i].Type == 'D')
			{
				start += cigar_ops[i].Length;
			}
			if (cigar_ops[i].Type == 'I')
			{
				start -= cigar_ops[i].Length;
			}
		}
		start += read.length() - clip.length();
	}
	else
	{
	  ++start;
	}
}

bool Tail::look_at_start(const BamAlignment & al) const
{
	if (!al.IsMapped()) return false;
	
	bool read1_start = params.read1_start;
	bool read1_end = params.read1_end;
	bool read2_start = params.read2_start;
	bool read2_end = params.read2_end;
	bool single_start = params.single_start;
	bool single_end = params.single_end;
	
	bool first = al.IsFirstMate();
	bool second = al.IsSecondMate();
	bool single = (!first && !second);

	vector<int> sizes, rpos, gpos;
	vector<CigarOp> cigar_ops;
	bool has_soft_clip = al.GetSoftClips(sizes,rpos,gpos);
	cigar_ops = al.CigarData;
	bool clip_at_start = (!al.IsReverseStrand() && cigar_ops[0].Type == 'S') || (al.IsReverseStrand() && cigar_ops[cigar_ops.size()-1].Type == 'S');
	bool clip_at_end = (al.IsReverseStrand() && cigar_ops[0].Type == 'S') || (!al.IsReverseStrand() && cigar_ops[cigar_ops.size()-1].Type == 'S');

	if (read1_start && first && has_soft_clip && clip_at_start && !al.IsReverseStrand())
	{
		return true;
	}
	else if (read2_start && second && has_soft_clip && clip_at_start && !al.IsReverseStrand())
	{
		return true;
	}
	else if (single_start && single && has_soft_clip && clip_at_start && !al.IsReverseStrand())
	{
		return true;
	}
	else if (read1_end && first && has_soft_clip && clip_at_end && al.IsReverseStrand())
	{
		return true;
	}
	else if (read2_end && second && has_soft_clip && clip_at_end && al.IsReverseStrand())
	{
		return true;
	}
	else if (single_end && single && has_soft_clip && clip_at_end && al.IsReverseStrand())
	{
		return true;
	}
	else
	{
		return false;
	}
}

int Tail::needleman_wunsch(const string & linker, const string & clip, string & aligned_linker) const
{
	int match=2;
    int mismatch=-1;
    int gep=-2;
    string res0 = linker;
    string res1 = clip;
    int len0 = res0.length();
    int len1 = res1.length();
	vector< vector<int> > smat (len0+1, vector<int> (len1+1));
    vector< vector<int> > tb (len0+1, vector<int> (len1+1));
	string aln0;
	string aln1;

    int i;
	int j;
   
    for (i = 0; i <= len0; ++i)
	{
		smat[i][0]= i * gep;
		tb[i][0]= 1;
	}
    for (j = 0; j <= len1; ++j)
	{
		smat[0][j] = j*gep;
		tb[0][j] = -1;
	}

    //Compute the optimal solution value and backtracking values bottom-up
    for (i = 1; i <= len0; ++i)
    {
		for (j = 1; j <= len1; j++)
		{
			//calculate alignment score
			int s;
			if (res0[i-1] == res1[j-1])
			{
				s = match;
			}
			else 
			{
				s = mismatch;
			}
	    
			int sub = smat[i-1][j-1] + s;
			int del = smat[i][j-1] + gep;
			int ins = smat[i-1][j] + gep;
	    
			if (sub > del && sub > ins)
			{
				smat[i][j] = sub;
				tb[i][j]=0;
			}
			else if (del > ins)
			{
				smat[i][j] = del;
				tb[i][j] = -1;
			}
			else 
			{
				smat[i][j] = ins;
				tb[i][j] = 1;
			}
		}
    }

    i = len0;
    j = len1;
    int aln_len=0;
	int optimal_score = smat[i][j];

	//Reconstruct optimal solution
    while (!(i==0 && j==0))
    {
		if (tb[i][j] == 0)
		{
			aln0.push_back(res0[--i]);
			aln1.push_back(res1[--j]);
		}
		else if (tb[i][j] == -1)
		{
			aln0.push_back('-');
			aln1.push_back(res1[--j]);
		}
		else if (tb[i][j] == 1)
		{
			aln0.push_back(res0[--i]);
			aln1.push_back('-');
		}
		++aln_len;
    }

	aligned_linker.clear();
	for (i = 0; i <= aln_len; ++i)
	{
		aligned_linker.push_back(aln0[aln_len-i]);
	}
    return optimal_score;
}

string Tail::reverse_complement(const string & s) const
{
	string rc = s;

	for (int i = 0; i < s.length(); ++i)
	{
		switch(s[i])
		{
		case 'A':
			rc[s.length()-i-1] = 'T';
			break;
		case 'C':
			rc[s.length()-i-1] = 'G';
			break;
		case 'G':
			rc[s.length()-i-1] = 'C';
			break;
		case 'T':
		case 'U':
			rc[s.length()-i-1] = 'A';
			break;
		}
	}

	return rc;
}

void Tail::compute_tail_comp()
{
	A_pct = 0;
	C_pct = 0;
	G_pct = 0;
	T_pct = 0;
	if (sequence == "")
		return;
	for (int i = 0; i < length; ++i)
	{
		switch(sequence[i])
		{
		case 'A':
			A_pct++;
			break;
		case 'C':
			C_pct++;
			break;
		case 'G':
			G_pct++;
			break;
		case 'T':
			T_pct++;
			break;
		}
	}
	A_pct /= length;
	C_pct /= length;
	G_pct /= length;
	T_pct /= length;
	A_pct *= 100;
	C_pct *= 100;
	G_pct *= 100;
	T_pct *= 100;
}

ostream & operator << (ostream & out, const Tail & tail)
{
	Parameters params = tail.params;
	if (tail.chr == "")	//No tail position called; do not output
		return out;
	int num_lines = 1;
	for (int i = 0; i < num_lines; ++i)
	{
		out << tail.read_id << "\t"
			<< tail.read << "\t";
		if (params.output_tail_seqs)
		{
			out << tail.sequence << "\t";
			out << tail.length << "\t";
		}
		if (params.output_linker)
		{
			out << tail.linker << "\t";
		}
		if (params.output_tail_comp)
		{
			out	<< tail.A_pct << "\t"
				<< tail.C_pct << "\t"
				<< tail.G_pct << "\t"
				<< tail.T_pct << "\t";
		}
		out	<< tail.quality << "\t"
			<< tail.chr << "\t"
			<< tail.start << "\n";
	}
	return out;
}

vector<string> split(const string& s, const string& delim, const bool keep_empty) {
    vector<string> result;
    if (delim.empty()) {
        result.push_back(s);
        return result;
    }
    string::const_iterator substart = s.begin(), subend;
    while (true) {
        subend = search(substart, s.end(), delim.begin(), delim.end());
        string temp(substart, subend);
        if (keep_empty || !temp.empty()) {
            result.push_back(temp);
        }
        if (subend == s.end()) {
            break;
        }
        substart = subend + delim.size();
    }
	result.pop_back();
    return result;
}

void read_chr_ids (BamReader & bam_file, map<int,string> & chr_ids)
{
	RefVector chrs = bam_file.GetReferenceData();
	for (int i = 0; i < bam_file.GetReferenceCount(); ++i)
	{
		chr_ids[i] = chrs[i].RefName;
	}
}

void read_regions(const string& regions_file, vector<pair<string,pair<int,int> > > & regions)
{
	ifstream fin;
	string temp;
	string chr; 
	int start, end;
	int num_regions = 0;
	fin.open(regions_file.c_str());
	if (fin.fail())
	{
		cout << "Regions file failed to open" << endl;
		return;
	}
	string line;
	getline(fin, line);
	while (!fin.eof())
		{
			vector<string> coords = split(line,",");
			chr = coords[1];
			start = atoi(coords[2].c_str());
			end = atoi(coords[3].c_str());
			pair<string,pair<int,int> > temp(chr,pair<int,int>(start,end));
			regions.push_back(temp);
			++num_regions;
			getline(fin, line);
		}
	cout << num_regions << " regions of interest read" << endl;
}

bool no_indels(const BamAlignment & al)
{
	for(int i = 0; i < al.CigarData.size(); ++i)
	{
		if (al.CigarData[i].Type == 'I' || al.CigarData[i].Type == 'D')
			return false;
	}
	return true;
}
