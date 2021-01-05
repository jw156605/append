#!/usr/bin/perl -w
######################################################################
# Program: graph_tails.pl
#
# Author: Joshua Welch
#
# Date: 5/28/2014
#
# Description: Given a tab-delimited text file containing tails
# resulting from the EnD-seq protocol, this script graphs the
# positional distribution of tails.
######################################################################

use strict;
use Getopt::Long;
use Data::Dumper;
use GD::Graph::bars;
use GD::Graph::lines;

sub gnuplot_stacked_bar
{
    my ($data_file, $gene, $out_file,$long,$start,$end) = @_;
    open GNUPLOT, "| gnuplot";
    print GNUPLOT <<gnuplot_commands ;
    set boxwidth 0.75 absolute
    set style fill   solid 1.00 border lt -1
    set key inside left top vertical Left reverse noenhanced autotitles columnhead nobox
    set key invert samplen 4 spacing 1 width 0 height 0 
    set style histogram rowstacked title  offset character 0,0,0
    set datafile missing '-'
    set style data histograms
    set xtics border in scale 0,0 nomirror
    set xtics  norangelimit font ",8"
    set title "Tail counts for $gene" 
    set format y "%.0f"
    set xlabel "Nucleotides from 3\' end"
    set ylabel "Count"
    i = 5
    set style line 1 lt rgbcolor "blue"
    set style line 2 lt rgbcolor "green"
    set style line 3 lt rgbcolor "orange"
    set style line 4 lt rgbcolor "red"
    set style increment user
    set term postscript eps color solid "Arial" 8
    set output "$out_file"
gnuplot_commands
if ($long == 1)
{
    if ($end == 30)
    {
	print GNUPLOT<<gnuplot_commands
	set xtics ("30" 1, "20" 11, "10" 21, "0" 31)
	plot "$data_file" using 2, for [i=3:5] '' using i
gnuplot_commands
    }
    else
    {
	print GNUPLOT<<gnuplot_commands
	set xtics ("150" 1, "140" 11, "130" 21, "120" 31, "110" 41, "100" 51, "90" 61, "80" 71, "70" 81, "60" 91, "50" 101, "40" 111, "30" 121, "20" 131, "10" 141, "0" 151) 
	plot "$data_file" using 2, for [i=3:5] '' using i
gnuplot_commands
    }
}
else
{
    print GNUPLOT<<gnuplot_commands
    plot "$data_file" using 2:xtic(1), for [i=3:5] '' using i
gnuplot_commands
}
}

sub gnuplot_single_tail_length
{
    my ($data_file, $gene, $option, $out_file,$long) = @_;
    my $column;
    my $color;
    if ($option == 0)
    {
        $column = 2;
        $color = "blue";
    }
    if ($option == 1)
    {
        $column = 3;
        $color = "green";
    }
    if ($option == 2)
    {
        $column = 4;
        $color = "orange";
    }
    if ($option == 3)
    {
        $column = 5;
        $color = "red";
    }
    open GNUPLOT, "| gnuplot";
    print GNUPLOT <<gnuplot_commands ;
    set boxwidth 0.75 absolute
    set style fill   solid 1.00 border lt -1
    set key inside right top vertical Left reverse noenhanced autotitles columnhead nobox
    set key invert samplen 4 spacing 1 width 0 height 0 
    set style histogram rowstacked title  offset character 0,0,0
    set datafile missing '-'
    set style data histograms
    set xtics border in scale 0,0 nomirror
    set xtics  norangelimit font ",8"
    set title "Tail counts for $gene" 
    set format y "%.0f"
    set xlabel "Nucleotides from 3\' end"
    set ylabel "Count"
    i = 5
    set term postscript eps color solid "Arial" 8
    set output "$out_file"
gnuplot_commands
if ($long == 1)
{
    #set xtics ("150" 1, "140" 11, "130" 21, "120" 31, "110" 41, "100" 51, "90" 61, "80" 71, "70" 81, "60" 91, "50" 101, "40" 111, "30" 121, "20" 131, "10" 141, "0" 151)
    print GNUPLOT<<gnuplot_commands
    set xtics ("30" 1, "20" 11, "10" 21, "0" 31)
    plot "$data_file" using $column title columnheader lt rgbcolor "$color"
gnuplot_commands
}
else
{
    print GNUPLOT<<gnuplot_commands
    plot "$data_file" using $column:xtic(1) title columnheader lt rgbcolor "$color"
gnuplot_commands
}
}

sub graph_transcript_ends
{
    my ($tails_in,$plot_type,$start,$end) = @_;
    my @shorts = (0,1,2,3);
    my @longs = (0,1,2,15);
    my @positions = reverse($start..$end);
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
	for (my $i = 0; $i < @shorts; $i++)
	{
	    if ($tail_length >= $shorts[$i] and $tail_length <= $longs[$i] and $mRNA_pos >= $start and $mRNA_pos <= $end)
	    {
		$gene =~ s/^\s+//;
		$gene =~ s/\s+$//;
		$mRNA_pos =~ s/^\s+//;
		$mRNA_pos =~ s/\s+$//;
		if (!exists $num_tails{$gene})
		{
		    $num_tails{$gene}= [[0,0,0,0,0,0,0,0,0,0,0],
					[0,0,0,0,0,0,0,0,0,0,0],
					[0,0,0,0,0,0,0,0,0,0,0],
					[0,0,0,0,0,0,0,0,0,0,0]];
		}
		
		$num_tails{$gene}->[$i]->[$mRNA_pos]++;
	    }
	}
    }
    mkdir "graphs";
    foreach my $gene (sort keys %num_tails)
    {
	open (OUT, ">graphs/$gene.dat");
	print OUT "Nucleotide\tNo_Tails\t1-nt_Tails\t2-nt_Tails\t>2-nt_Tails\n";
	my $i = 0;
	my @counts = @{$num_tails{$gene}};
	my @data = (\@positions, $counts[0], $counts[1], $counts[2], $counts[3]);
	print OUT "$gene\n";
	foreach my $position (@positions)
	{
	    print OUT "$position\t";
	    if (defined $counts[0]->[$position])
	    {
		print OUT $counts[0]->[$position]."\t"; 
	    }
	    else 
	    {
		print OUT "0\t";
	    }
	    if (defined $counts[1]->[$position])
	    {
		print OUT $counts[1]->[$position]."\t"; 
	    }
	    else 
	    {
		print OUT "0\t";
	    }
	    if (defined $counts[2]->[$position])
	    {
		print OUT $counts[2]->[$position]."\t"; 
	    }
	    else 
	    {
		print OUT "0\t";
	    }
	    if (defined $counts[3]->[$position])
	    {
		print OUT $counts[3]->[$position]."\n"; 
	    }
	    else 
	    {
		print OUT "0\n";
	    }
	}
	close OUT;
        if ($plot_type == 4) #Plot 0, 1, 2, 3-15-nt tails on single plot
        {
	    my $long = 0;
	    if ($end-$start+1 > 20)
	    {
		$long = 1;
	    }
            gnuplot_stacked_bar("graphs/$gene.dat",$gene,"graphs/$gene"."_$start"."_$end.eps",$long,$start,$end);    
        }
        else #Single tail length only
        {
	    my $long = 0;
	    if ($end-$start+1 > 20)
	    {
		$long = 1;
	    }
            gnuplot_single_tail_length("graphs/$gene.dat",$gene,$plot_type,"graphs/$gene"."_$start"."_$end.eps",$long)
        }
	
	#my $graph = GD::Graph::bars->new(800, 600);
	#$graph->set_legend_font("GD::gdFontLarge");
        #$graph->set_legend('No Tails', '1-nt Tails', '2-nt Tails','> 2-nt Tails');
	#$graph->set( 
	#    x_label           => 'Nucleotides from 3\' End',
	#    y_label           => 'Tails mapped to Gene',
	#    title             => "Tail Counts for $gene",
	#    x_label_skip      => 10,
	#    dclrs             => ['lred','lyellow','blue','lyellow'],
	#    overwrite => 2
	    
	#) or die "Error plotting graph: " . $graph->error;
	#my $plot;
	#eval{$plot = $graph->plot(\@data) or die "Error plotting graph for $gene end: " . $graph->error;};
	#if ($@)
	#{
	#    print $@;
	#}
	#else
	#{
        #    mkdir "graphs";
	#    open (IMG, ">graphs/$gene"."_end.png");
	#    binmode IMG;
	#    print IMG $plot->png;    
	#    close IMG;    
	#}
    }
}

my ($tails_in,$plot_type, $start, $end);
if (!GetOptions( '-t=s' => \$tails_in,
                 '-p=s' => \$plot_type,
                 '-s=s' => \$start,
                 '-e=s' => \$end
		 ))
{
print <<END_USAGE;
Usage:
  -t <tail file> -s <sample ID> -p <plot type = 0,1,2,3, or 4> -s <mRNA start position> -e <mRNA end position>
END_USAGE
  exit(1);
}
if (!defined $start)
{
    $start = 0;
}
if (!defined $end)
{
    $end = 10;
}
graph_transcript_ends($tails_in,$plot_type,$start,$end);
