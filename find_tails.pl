#!/usr/bin/perl -w
######################################################################
# Program: find_tails.pl
#
# Author: Joshua Welch
#
# Date: 5/28/2014
#
# Description: Finds non-templated tails added to mRNA transcripts. 
# Takes either FASTQ or SAM input. If the user specifies FASTQ input,
# the program aligns the reads to the specified reference genome,
# finds tailed reads, and outputs the tails overlapping regions of
# interest. The user also has the option of specifying a SAM file
# containing only tailed reads; in this case, the tails overlapping
# regions of interest are output.
######################################################################

use POSIX;
use strict;
use IO::File;
use Getopt::Long;
use Data::Dumper;

sub get_csv_column_pairs
{
    my $filename = $_[0];
    my $column1 = $_[1];
    my $column2 = $_[2];
    my %column_pairs;
    open(FILE, "$filename") or die("Unable to open file $filename");
    my @file = <FILE>;
    foreach my $line (@file)
    {
      my @columns = split(/,/, $line);
      $columns[$column1]=~s/^\s+//;
      $columns[$column1]=~s/\s+$//;
      $columns[$column2]=~s/^\s+//;
      $columns[$column2]=~s/\s+$//;
      my $col1 = $columns[$column1];
      my $col2 = $columns[$column2];
      $column_pairs{$col1} = $col2;
     }
    return \%column_pairs;
}

sub run_aligner
{
    my ($aligner, $index_dir, $fastq1, $fastq2, $chr_dir) = @_;
    if ($aligner eq 'star')
    {
	system("/proj/prins_lab/Josh/STAR_2.3.1x/STAR --genomeDir $index_dir --readFilesIn $fastq1 $fastq2 --runThreadN 8 --genomeLoad NoSharedMemory --alignIntronMax 1 --alignIntronMin 2 --outFilterMismatchNmax 4");
    }
    if ($aligner eq 'mapsplice')
    {
	print "python /nas02/apps/mapsplice-2.1.4/src/MapSplice-v2.1.4/mapsplice.py -x $index_dir -c $chr_dir -1 $fastq1 -2 $fastq2 -p 8 --min-map-len 25\n";
	if ($index_dir eq '')
	{
	    system("python /nas02/apps/mapsplice-2.1.4/src/MapSplice-v2.1.4/mapsplice.py -c $chr_dir -1 $fastq1 -2 $fastq2 -p 8 --min-map-len 25");    
	}
	else
	{
	    system("python /nas02/apps/mapsplice-2.1.4/src/MapSplice-v2.1.4/mapsplice.py -x $index_dir -c $chr_dir -1 $fastq1 -2 $fastq2 -p 8 --min-map-len 40");    
	}
	system("mv ./mapsplice_out/alignments.sam Aligned.out.sam");
    }
    if ($aligner eq 'bowtie2')
    {
	system("bowtie2 --local -x $index_dir -1 $fastq1 -2 $fastq2 -S Aligned.out.sam -p 8"); #-U $fastq1
    }
    if ($aligner eq 'bowtie')
    {
	print "bowtie -a -1 $fastq1 -2 $fastq2 -S -p 8 $index_dir > Aligned.out.sam\n";
	system("bowtie -a -1 $fastq1 -2 $fastq2 -S -p 8 $index_dir > Aligned.out.sam");
    }
}

sub read_genome
{
    my $genome_dir = $_[0];
    opendir(DIR, $genome_dir) or die "cannot find genome directory $genome_dir";
    my @fastas = grep(/\.fa$/,readdir(DIR));
    my %genome;
    foreach my $file (@fastas) {
        $file =~ /(.*)\.fa/;
        my $chr_name = $1;
        my $seq_obj = Bio::SeqIO->new(-file => "$genome_dir/$file", '-format' => 'fasta');
        $genome{$chr_name} = $seq_obj;
    }
    return \%genome;
}

sub get_num_alignments
{
    my $sam_line = $_[0];
    my $tag = 'NH';
    $sam_line =~ /(NH|IH|XN|XM):i:(\d+)/; #STAR indicates number of alignments with the NH field, but MapSplice uses IH
    if ($1 eq 'XN' or $1 eq 'XM')
    {
	return 1;
    }
    return $2;
}

sub has_soft_clip
{
    my $sam_line = $_[0];
    my @line_parts = split(/\t/, $sam_line);
    my $cigar_string = $line_parts[5];
    return ($cigar_string =~ /[0-9]+S/);
}

sub get_clip_length
{
    my $sam_line = $_[0];
    my @line_parts = split(/\t/, $sam_line);
    my $cigar_string = $line_parts[5];
    $cigar_string =~ /([0-9]+)S/;
    return $1;
}

sub get_tailed_reads
{
    my $sam_in = $_[0];
    my $sam_out = $_[1];
    my %multimapped;
    open(SAM_IN, $sam_in);
    open(OUT_FILE, ">$sam_out");
    my $num_tailed = 0;
    while (my $sam_line = <SAM_IN>)
    {
	if ($sam_line =~ /^@.+/) #Header line
	{
	    print OUT_FILE $sam_line;
	    next;
	}
	if (get_flag($sam_line,2)) #Skip unaligned reads
	{
	    next;
	}
        my $second = get_flag($sam_line,7); #Only first read in pair contains tail
        my $num_alignments = get_num_alignments($sam_line);
	if (has_soft_clip($sam_line) and !$second and $num_alignments <= 1)
	{
	    print OUT_FILE $sam_line;
	    $num_tailed++;
	}
        if (has_soft_clip($sam_line) and !$second and $num_alignments > 1)
        {
            my @line_parts = split(/\s+/,$sam_line);
            my $read_id = $line_parts[0];
            if (!exists $multimapped{$read_id})
            {
                my @temp;
                $multimapped{$read_id} = [@temp];
            }
            push(@{$multimapped{$read_id}},$sam_line);
        }
    }
    foreach my $read (keys %multimapped)
    {
        my $num_mappings = @{$multimapped{$read}};
        my $rand_ind = int(rand($num_mappings));
        print OUT_FILE $multimapped{$read}[$rand_ind];
    }
    close SAM_IN;
    close OUT_FILE;
}

sub get_flag
{
    my $sam_line = $_[0];
    my $flag_num = $_[1];
    my @line_parts = split(/\t/, $sam_line);
    my $value = $line_parts[1];
    return ($value & (2**$flag_num));
}

sub reverse_complement
{
    my $dna = $_[0];
    # reverse the DNA sequence
    my $revcomp = reverse($dna);
    # complement the reversed DNA sequence
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

sub is_rev_comp
{
    my $sam_line = $_[0];
    my $read_rev = get_flag($sam_line,4);
    my $first = get_flag($sam_line,6);
    my $mate = get_flag($sam_line,7);
    my $mate_rev = get_flag($sam_line,5);
    return (($first and $read_rev) or ($mate and $mate_rev) or (!$first and !$mate and $read_rev));
}

sub get_tail
{
    my $clip = $_[0];
    my $adapter = $_[1];
    my $aligned_adapter = needleman_wunsch($adapter,$clip);
    my $tail_length = 0;
    for (my $i = @{$aligned_adapter}-1; $i > 0; $i--)
    {
	if ($aligned_adapter->[$i] eq '-')
	{
	    $tail_length++;
	}
	else
	{
	    last;
	}
    }
    if ($tail_length == 0)
    {
	return '';
    }
    return reverse_complement(substr($clip, length($clip)-$tail_length, $tail_length));
}

sub get_adapter
{
    my $clip = $_[0];
    my $adapter = $_[1];
    my $aligned_adapter = needleman_wunsch($adapter,$clip);
    my $tail_length = 0;
    for (my $i = @{$aligned_adapter}-1; $i > 0; $i--)
    {
	if ($aligned_adapter->[$i] eq '-')
	{
	    $tail_length++;
	}
	else
	{
	    last;
	}
    }
    if ($tail_length == 0)
    {
	return $clip;
    }
    return substr($clip, 0, length($clip)-$tail_length+1);
}

sub get_read_start
{
    my $sam_line = $_[0];
    my $first = get_flag($sam_line,6);
    my $mate = get_flag($sam_line,7);
    my @line_parts = split(/\t/, $sam_line);
    if ($first)
    {
	return $line_parts[3]-1;    
    }
    elsif ($mate)
    {
	return $line_parts[7]-1;
    }
    else
    {
        return $line_parts[3]-1;
    }
}

sub needleman_wunsch
{
    my $match=10;
    my $mismatch=-10;
    my $gop=-10;
    my $gep=-10;
    my $al0 = $_[0];
    my $al1 = $_[1];
    my @name_list0;
    my @name_list1;
    my @seq_list0;
    my @seq_list1;
    my @res0;
    my @res1;
    my @smat;
    my @tb;
    my @aln0;
    my @aln1;

    my $i;
    my $j;
    #Split the sequences into character arrays
    for ($i=0; $i < length($al0); $i++)
    {
        $res0[$i] = substr($al0, $i, 1);
    }
    for ($i=0; $i < length($al1); $i++)
    {
        $res1[$i] = substr($al1, $i, 1);
    }

    my $len0=@res0;
    my $len1=@res1;
    
    #Initialize dynamic programming matrix
    for ($i=0; $i<=$len0; $i++){$smat[$i][0]=$i*$gep;$tb[$i][0 ]= 1;}
    for ($j=0; $j<=$len1; $j++){$smat[0][$j]=$j*$gep;$tb[0 ][$j]=-1;}
	    
    #Compute the optimal solution value and backtracking values bottom-up
    for ($i=1; $i<=$len0; $i++)
    {
	for ($j=1; $j<=$len1; $j++)
	{
	    #calculate alignment score
	    my $s;
	    if ($res0[$i-1] eq $res1[$j-1]){$s=$match;}
	    else {$s=$mismatch;}
	    
	    my $sub=$smat[$i-1][$j-1]+$s;
	    my $del=$smat[$i  ][$j-1]+$gep;
	    my $ins=$smat[$i-1][$j  ]+$gep;
	    
	    if   ($sub>$del && $sub>$ins){$smat[$i][$j]=$sub;$tb[$i][$j]=0;}
	    elsif($del>$ins){$smat[$i][$j]=$del;$tb[$i][$j]=-1;}
	    else {$smat[$i][$j]=$ins;$tb[$i][$j]=1;}
	}
    }
    
    $i=$len0;
    $j=$len1;
    my $aln_len=0;
    #Reconstruct optimal solution
    while (!($i==0 && $j==0))
    {
	if ($tb[$i][$j]==0)
		{
		$aln0[$aln_len]=$res0[--$i];
		$aln1[$aln_len]=$res1[--$j];
		}
	elsif ($tb[$i][$j]==-1)
		{
		$aln0[$aln_len]='-';
		$aln1[$aln_len]=$res1[--$j];
		}
	elsif ($tb[$i][$j]==1)
		{
		$aln0[$aln_len]=$res0[--$i];
		$aln1[$aln_len]='-';
		}
	$aln_len++;
    }
    #Output solution and return aligned sequence 1
    #print "Adapter:\n";
    #for ($i=$aln_len-1; $i>=0; $i--){print $aln0[$i];}
    #print "\n";
    #print "Clipped sequence:\n";
    #for ($j=$aln_len-1; $j>=0; $j--){print $aln1[$j];}
    #print "\n";
    my @forward_aln = reverse(@aln0);
    return \@forward_aln;
}

sub output_tails
{
    my ($sam_file,$adapter,$chrs,$starts,$coding_starts,$coding_ends,$ends,$strands,$short,$long,$out) = @_;
    $adapter = reverse_complement($adapter);
    open(SAM_IN, $sam_file);
    open(OUT, ">$out");
    print OUT "Read ID\tRead\tGene\tmRNA Position\tTail Sequence\tTail Length\tA Count\tC Count\tG Count\tT Count\tQuality\tTail Chromosome\tTail Start\n";
    my $num_tails = 0;
    while (my $sam_line = <SAM_IN>)
    {
	if ($sam_line =~ /^@.+/) #Header line
	{
	    next;
	}
	my @line_parts = split(/\t/, $sam_line);
	my $read_id = $line_parts[0];
        my $quality = $line_parts[4];
	my $read = $line_parts[9];
	my $chr = $line_parts[2];
	my $start = get_read_start($sam_line);

	#Get length of soft clip
	my $clip_length = get_clip_length($sam_line);
	#Reverse if necessary
	if (is_rev_comp($sam_line))
	{
            $read = reverse_complement($read);
	}
	#Get clipped sequence and trim adapter; tail is what's left
	my $clip = substr($read,0,$clip_length);
	my $adapter_seq = get_adapter($clip, $adapter);
        my $tail = get_tail($clip, $adapter);
	my $tail_len = length($tail);

	#Check for case in which first nt of adapter matches genome
	if ($tail_len == 0 and 
	    substr($clip,$clip_length-1,1) ne substr($adapter,-1,1) and
	    substr($adapter,-1,1) eq substr($read,$clip_length,1))
	{
	    ++$clip_length;
	    if (!is_rev_comp($sam_line))
	    {
		++$start;
	    }
	    #print "----\nClip: $clip\nAdapter: $adapter\nRead: $read\n";
	}
	#else
	#{
	#    if ($tail_len == 0)
	#    {
	#	print "----\nClip: $clip\nAdapter: $adapter\nRead: $read\n";
	#	print "". substr($clip,$clip_length-1,1) ."\n";
	#	print "". substr($adapter,-1,1) ."\n";
	#	print "". substr($adapter,-1,1) ."\n";
	#	print "". substr($read,$clip_length,1) ."\n----\n";
	#    }
	#}

	if ($short <= $tail_len and $tail_len <= $long)
	{
	    if (is_rev_comp($sam_line))
            {
                my $cigar = $line_parts[5];
		my @cigar_ops = split(/([MIDNSHPX=])/,$cigar);
		
		for (my $i = 0; $i < @cigar_ops; $i++)
		{
		    if ($cigar_ops[$i+1] eq 'D')
		    {
			$start += $cigar_ops[$i];
		    }
		    if ($cigar_ops[$i+1] eq 'I')
		    {
			$start -= $cigar_ops[$i];
		    }
		    $i++;
		}
		$start = $start + length($read) - $clip_length;
            }
            my $numAs = ($tail =~ tr/A//);
            my $numCs = ($tail =~ tr/C//);
            my $numGs = ($tail =~ tr/G//);
            my $numTs = ($tail =~ tr/T//);
	    foreach my $gene_id (keys %{$chrs})
            {
                if ($chr eq $chrs->{$gene_id} and $starts->{$gene_id} <= $start and $start <= $ends->{$gene_id})
                {
                    #print OUT "Read ID\tRead\tGene\tmRNA Position\tTail Sequence\tTail Length\tA Count\tC Count\tG Count\tT Count\tQuality\tTail Chromosome\tTail Start\n";
		    #my $mRNA_pos = $start-$coding_ends->{$gene_id};
		    my $mRNA_pos = $ends->{$gene_id}-$start;
		    if ($strands->{$gene_id} eq '-')
		    {
			#$mRNA_pos = $coding_starts->{$gene_id}-$start;
			$mRNA_pos = $start-$starts->{$gene_id};
		    }
                    $read = reverse_complement($read);
	            print OUT "$read_id\t$read\t$gene_id\t$mRNA_pos\t$tail\t" . length($tail) . "\t$numAs\t$numCs\t$numGs\t$numTs\t$quality\t$chr\t$start\n";
                }
            }
	}
    }
    close SAM_IN;
    close OUT;
}

sub draw_graphs
{
    my ($chrs,$starts,$ends,$short,$long,$tails_in,$sample_id) = @_;
    my %num_tails;
    open (TAILS, $tails_in);
    my $line = <TAILS>;
    while ($line = <TAILS>)
    {
	my @line_parts = split(/\t/, $line);
	my $chr = $line_parts[11];
	my $start = $line_parts[12];
	$chr =~s/^\s+//;
	$chr =~s/\s+$//;
	$start =~s/^\s+//;
	$start =~s/\s+$//;
        if (!exists $num_tails{$chr}{$start})
	{
	    $num_tails{$chr}{$start} = 0;
	}
	$num_tails{$chr}{$start}++;
    }
    open (OUT, ">$sample_id.bedgraph");
    print OUT "track type=bedGraph name=\"$sample_id\" graphType=bar\n";
    foreach my $chr (sort (keys %num_tails))
    {
	foreach my $position (sort {$a <=> $b} (keys %{$num_tails{$chr}}))
	{
	    if ($chr ne 'chrM_rCRS')
	    {
		print OUT "$chr $position $position ".$num_tails{$chr}{$position}."\n";
	    }    
	}
    }
    close OUT;
}

sub graph_each_gene
{
    my ($chrs,$starts,$ends,$short,$long,$start,$end,$tails_in,$sample_id) = @_;
    my %num_tails;
    my %total_per_gene;
    open (TAILS, $tails_in);
    my $line = <TAILS>;
    while ($line = <TAILS>)
    {
	my @line_parts = split(/\t/, $line);
	my $gene = $line_parts[2];
	my $mRNA_pos = $line_parts[3];
	my $tail_length = $line_parts[5];
	if (!exists $total_per_gene{$gene})
	{
	    $total_per_gene{$gene} = 0;
	}
	$total_per_gene{$gene}++;
	if ($tail_length < $short or $tail_length > $long or $mRNA_pos < $start or $mRNA_pos > $end)
	{
	    next;
	}
	$gene =~ s/^\s+//;
	$gene =~ s/\s+$//;
	$mRNA_pos =~ s/^\s+//;
	$mRNA_pos =~ s/\s+$//;
        if (!exists $num_tails{$gene}{$mRNA_pos})
	{
	    $num_tails{$gene}{$mRNA_pos} = 0;
	}
	
	$num_tails{$gene}{$mRNA_pos}++;
    }
    mkdir "graphs";
    
    foreach my $gene (keys %num_tails)
    {
	#foreach my $pos (keys %{$num_tails{$gene}})
	#{
	#    $num_tails{$gene}{$pos} /= $total_per_gene{$gene};
	#    $num_tails{$gene}{$pos} *= 100;
	#}
	my @positions = sort { $a <=> $b } keys %{$num_tails{$gene}};
	my $min_pos = $positions[0];
	my $max_pos = $positions[@positions-1];
	@positions = reverse($min_pos..$max_pos);
	
	my @counts = map($num_tails{$gene}{$_},@positions); #/$total_per_gene{$gene}*100
	my @data = (\@positions, \@counts);
	my $graph = GD::Graph::bars->new(800, 600);
	$graph->set( 
	    x_label           => 'Nucleotides from 3\' End',
	    y_label           => 'Number of tails mapped to Gene',
	    title             => "$short-$long nt Tails for $gene",
	    x_label_skip      => 10,
	    dclrs       => ['red']
	    
	) or die "Error plotting graph: " . $graph->error;
	my $plot;
	eval{$plot = $graph->plot(\@data) or die "Error plotting graph for $gene $short-$long: " . $graph->error;};
	if ($@)
	{
	    print $@;
	}
	else
	{
	    open (IMG, ">graphs/$gene"."_$short"."_$long"."_$start"."_$end".".png");
	    binmode IMG;
	    my $gene_length = $ends->{$gene}-$starts->{$gene}+1;
	    if (scalar(@positions) > 0 and scalar(@counts) > 0)
	    {
		#print "Gene: $gene Gene length: $gene_length Tail lengths: $short-$long Positions: " .scalar(@positions). " counts: " .scalar(@counts) ."\n";
		#print Dumper(@data);
		print IMG $plot->png;    
	    }
	    close IMG;    
	}
    }
}

sub find_tailed_loci
{
    my ($chrs,$starts,$ends,$bg_in,$sample_id) = @_;
    my %num_tails;
    my $total_tails = 0;
    open (TAILS, $bg_in);
    my $line = <TAILS>;
    while ($line = <TAILS>)
    {
	my @line_parts = split(/ /, $line);
	my $chr = $line_parts[0];
	my $start = $line_parts[1];
	my $num = $line_parts[3];
	$total_tails += $num;
	foreach my $ID (%{$chrs})
	{
	    if (!exists $chrs->{$ID} or !exists $starts->{$ID} or !exists $ends->{$ID})
	    {
		next;
	    }
	    if ($chr eq $chrs->{$ID} and $starts->{$ID} <= $start and $start <= $ends->{$ID})
	    {
		if (!exists $num_tails{$ID})
		{
		    $num_tails{$ID} = 0;
		}
		$num_tails{$ID} += $num;
	    }
	}
    }
    open (OUT, ">$sample_id"."_loci.csv");
    print OUT "all,,,,$total_tails\n";
    foreach my $ID (sort keys %{$chrs})
    {
	if (!exists $num_tails{$ID})
	{
	    $num_tails{$ID} = 0;
	}
	print OUT "$ID,$chrs->{$ID},$starts->{$ID},$ends->{$ID},$num_tails{$ID}\n";
    }
}

sub plot_tail_composition
{
    my ($tails_in,$short,$long,$start,$end,$sample_id) = @_;
    
    my %pwm;
    $pwm{'A'}=[0 x ($short-$long+1)];
    $pwm{'C'}=[0 x ($short-$long+1)];
    $pwm{'G'}=[0 x ($short-$long+1)];
    $pwm{'T'}=[0 x ($short-$long+1)];
    my %di_nt_freqs;
    foreach my $nt (keys %pwm)
    {
	for (my $i = $short; $i <= $long; $i++)
	{
	    $pwm{$nt}[$i-$short]=0;
	}
    }
    open (TAILS, $tails_in);
    my $line = <TAILS>;
    while ($line = <TAILS>)
    {
	my @line_parts = split(/\t/, $line);
	my $tail_pos = $line_parts[3];
	my $tail_seq = $line_parts[4];
	my @tail = split("",$tail_seq);
	my $tail_length = $line_parts[5];
	if ($tail_pos < $start or $tail_pos > $end or $tail_length < $short or $tail_length > $long) #Skip tails outside specified lengths and mRNA positions
	{
	    next;
	}
	#print Dumper(@tail);
	if ($tail_length == 2)
	{
	    if (!exists $di_nt_freqs{$tail_seq})
	    {
		$di_nt_freqs{$tail_seq} = 0;
	    }
	    $di_nt_freqs{$tail_seq}++;
	}
	for (my $i = $short; $i <= @tail and $i <= $long; $i++)
	{
	    if (exists $pwm{$tail[$i-1]}) #Count ACGT only
	    {
		$pwm{$tail[$i-1]}[$i-$short]++;
	    }
	}
    }
    
    for (my $i = $short; $i <= $long; $i++)
    {
	my $col_sum = 0;
	foreach my $nt (keys %pwm)
	{
	    $col_sum += $pwm{$nt}[$i-$short];
	}
	foreach my $nt (keys %pwm)
	{
	    if ($col_sum != 0){
		$pwm{$nt}[$i-$short] /= $col_sum;
	    }
	}
    }
    my $pwm_file = "$sample_id"."_pwm.csv";
    open (PWM_OUT, ">$pwm_file");
    foreach my $nt (sort keys %pwm)
    {
	#print $nt;
	for (my $i = $short; $i < $long; $i++)
	{
	    my $rounded = sprintf("%.3f",$pwm{$nt}[$i-$short]);
	    print PWM_OUT "$rounded,";
	    #print " $rounded";
	}
	my $rounded = sprintf("%.3f",$pwm{$nt}[$long-$short]);
	print PWM_OUT "$rounded\n";
	#print " $rounded\n";
    }
    #system ("Rscript /proj/prins_lab/Josh/code/berryplot.R $pwm_file $sample_id"."_berryplot.png");
    close PWM_OUT;
    open (OUT, ">$sample_id"."_di_nt.csv");
    foreach my $di_nt (keys %di_nt_freqs)
    {
	print OUT "$di_nt,".($di_nt_freqs{$di_nt})."\n";
    }
}

my ($fastq1, $fastq2, $index_dir, $region_file, $adapter, $sam_file, $graph, $short, $long, $aligner_tool, $chr_dir);
if (!GetOptions( '-1=s' => \$fastq1,
		 '-2=s' => \$fastq2,
                 '-i=s' => \$index_dir,
                 '-r=s' => \$region_file,
                 '-a=s' => \$adapter,
                 '-b=s' => \$sam_file,
                 '-s=s' => \$short,
                 '-l=s' => \$long,
                 'graph' => \$graph,
		 '-t=s' => \$aligner_tool,
		 '-c=s' => \$chr_dir
		 ))
{
print <<END_USAGE;
Usage:
  -1 <FASTQ file> [-2 <FASTQ file for second pair>] -i <genome index directory> -r <region file> -a <adapter sequence> -c <chromosome directory> [-t <RNA-seq alignment tool>] [-b <SAM file with tailed reads>] [--graph] [-s <shortest allowed tail>] [-l <longest allowed tail>]
END_USAGE
  exit(1);
}
my $sample_id = 'tail_run';
#Determine whether FASTQ or SAM input is specified and run aligner if necessary
if (defined $fastq1 and defined $sam_file)
{
    print "Please specify either FASTQ or SAM input, but not both.\n";
    print <<END_USAGE;
Usage:
  -1 <FASTQ file> [-2 <FASTQ file for second pair>] -i <STAR index directory> -r <region file> -a <adapter sequence> [-b <SAM file with tailed reads>] [-g] [-s <shortest allowed tail>] [-l <longest allowed tail>]
END_USAGE
  exit(1);
}
elsif (!defined $fastq1 and !defined $sam_file)
{
    print "Please specify either FASTQ or SAM input.\n";
    print <<END_USAGE;
Usage:
  -1 <FASTQ file> [-2 <FASTQ file for second pair>] -i <STAR index directory> -r <region file> -a <adapter sequence> [-b <SAM file with tailed reads>] [-g] [-s <shortest allowed tail>] [-l <longest allowed tail>]
END_USAGE
exit(1);
}
elsif (defined $fastq1)
{
    if (!defined $aligner_tool)
    {
	$aligner_tool = 'bowtie2';
    }
    print "FASTQ input specified. Aligning reads...\n";
    run_aligner($aligner_tool,$index_dir, $fastq1, $fastq2,$chr_dir);
    print "Finding tailed reads...\n";
    get_tailed_reads('Aligned.out.sam','Aligned.out.tails.sam');
    $sam_file = 'Aligned.out.tails.sam';
    $fastq1 =~ /.*\/(.*)\.fastq$/;
    $sample_id = $1;
}
if (!defined $short)
{
    $short = 1;
}
if (!defined $long)
{
    $long = 15;
}

print "Scanning SAM file for tails overlapping specified regions...\n";

#Load region information
my $chrs = get_csv_column_pairs($region_file, 0, 1);
my $starts = get_csv_column_pairs($region_file, 0, 2);
my $ends = get_csv_column_pairs($region_file, 0, 3);
my $strands = get_csv_column_pairs($region_file, 0, 4);
my $coding_starts = get_csv_column_pairs($region_file, 0, 5);
my $coding_ends = get_csv_column_pairs($region_file, 0, 6);

output_tails($sam_file,$adapter,$chrs,$starts,$coding_starts,$coding_ends,$ends,$strands,$short,$long,"$sample_id"."_tails.txt");
plot_tail_composition("$sample_id"."_tails.txt",1,$long,0,200,$sample_id);
print "Identifying annotated loci with tails...\n";
find_tailed_loci($chrs,$starts,$ends,"$sample_id".".bedgraph",$sample_id);
print "Done.\n";
