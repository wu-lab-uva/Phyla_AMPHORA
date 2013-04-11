# Phyla_AMPHORA A Phylum-Specific Automated Phylogenomic Inference Application for Bacterial Sequences. 
# Copyright 2012 by Martin Wu
 
# This file is part of Phyla_AMPHORA.

# Phyla_AMPHORA is free software: you may redistribute it and/or modify its under the terms of the 
# GNU General Public License as published by the Free Software Foundation; either version 2 of
# the License, or any later version.

# Phyla_AMPHORA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details (http://www.gnu.org/licenses/).
 

# For any other inquiries send an Email to Martin Wu
#       mw4yv@virginia.edu
 
# When publishing work that is based on the results from AMPHORA please cite:
# Wang Z and Wu M: A Phylum-level Bacterial Phylogenetic Marker Database. Mol. Biol. Evol. Advance Access publication March 21, 2013. doi:10.1093/molbev/mst059

#! /usr/bin/perl
use strict;
use Bio::SeqIO;
use Getopt::Long;
use Bio::Root::Root;
$Bio::Root::Root::DEBUG = -1;

my $Phyla_AMPHORA_home = $ENV{'Phyla_AMPHORA_home'};
my (%markerlist, %seq) = ();
my $help = undef;
my $evalue = 1e-7;
my $ref_dir = "$Phyla_AMPHORA_home/Marker";
my $is_dna = undef;
my $phylum_index = 0;


my $usage = qq~
This tool will search for bacterial and archaeal phylogenetic markers for a given fasta sequence file.
The output of this tool is a collection of marker protein sequences in fasta format. For example, rpoB.pep, dnaG.pep.

When DNA sequences are used, this program first identifies ORFs longer than 100 bp in all six reading frames, then scans the translated peptide sequences for the phylogenetic markers.

Usage: $0 <options> sequence-file

Options:
	-Phylum: 	0. All (Default)
			1. Alphaproteobacteria
			2. Betaproteobacteria
			3. Gammaproteobacteria
			4. Deltaproteobacteria
			5. Epsilonproteobacteria
			6. Acidobacteria
			7. Actinobacteria
			8. Aquificae
			9. Bacteroidetes
			10. Chlamydiae/Verrucomicrobia
			11. Chlorobi
			12. Chloroflexi
			13. Cyanobacteria
			14. Deinococcus/Thermus
			15. Firmicutes
			16. Fusobacteria
			17. Planctomycetes
			18. Spirochaetes
			19. Tenericutes
			20. Thermotogae
	-DNA: input sequences are DNA. Default: no
	-Evalue: HMMER evalue cutoff. Default: 1e-7 
	-ReferenceDirectory: the file directory that contain the reference alignments, hmms and masks. Default: $Phyla_AMPHORA_home/Marker
	-Help: print help;
~;



GetOptions (	'DNA'=>\$is_dna,
		'Evalue=f'=>\$evalue,
		'Phylum=i'=>\$phylum_index,
		'ReferenceDirectory=s'=>\$ref_dir,
		'Help'=>\$help) || die "Invalid command line options\n";

die $usage if $help;
die $usage unless $ARGV[0];

my %phylum = ('0', '*', '1', 'Alpha', '2', 'Beta', '3', 'Gamma', '4', 'Delta', '5', 'Epsilon', '6', 'Acidobacteria', '7', 'Actinobacteria', '8', 'Aquificae', '9', 'Bacteroidetes', '10', 'ChlamydiaeVerrucomicrobia', '11', 'Chlorobi', '12', 'Chloroflexi', '13', 'Cyanobacteria', '14', 'DeinococcusThermus', '15', 'Firmicutes', '16', 'Fusobacteria', '17', 'Planctomycetes', '18', 'Spirochaetes', '19', 'Tenericutes', '20', 'Thermotogae');

my $input_seq = $ARGV[0];

if ($is_dna) {
	system ("getorf -sequence $ARGV[0] -outseq $ARGV[0].orf -table 1 -minsize 100");
	$input_seq = "$ARGV[0].orf";
}

# Stage HMMs
system ("cat $ref_dir/$phylum{$phylum_index}*hmm > $$.hmm");

# Get candidates
my $candidate = get_candidates();

# Get orthologs
get_orthologs($candidate);

# clean up
system ("rm $$.*");

####################################################################################################################
	
sub get_candidates {
	my (%score, %candidates, %query_length, %hmm_length, %query_match, %hmm_match, %hits) =();
	system ("hmmsearch -Z 5000 -E $evalue --domE $evalue --domtblout $$.hmmsearch -o /dev/null $$.hmm $input_seq");		# fix the number of sequences in the database for E-value calculation

	open (IN, "$$.hmmsearch") || die "Can't open $$.hmmsearch";
	while (<IN>) {
		chop;
		next if /^#/;
		my ($query, $query_accession, $qlength, $hmm, $hmm_accession, $hmm_length, $evalue, $score, $bias, $domain, $domain_number, $dom_evalue, $ievalue, $dom_score, $dom_bias, $hmm_start, $hmm_stop, $query_start, $query_stop, $rest) = split /\s+/;
		if ((! exists $score{$query}) or ($score{$query} < $score)) {
			$score{$query} = $score;
			$candidates{$query} = $hmm;
		}
		$query_length{$query} = $qlength;
		$hmm_length{$hmm} = $hmm_length;
		$query_match{$hmm}{$query} += ($query_stop - $query_start);
		$hmm_match{$hmm}{$query} += ($hmm_stop - $hmm_start);
	}
	close IN;
	
	while ( my ($query, $hmm) = each (%candidates) ) {
		next if ( ($query_match{$hmm}{$query}/$query_length{$query} < 0.7) and ($hmm_match{$hmm}{$query}/$hmm_length{$hmm} < 0.7) ); 	#ignore the hit if the match is partial
		$hits{$hmm}{$query} = 1;
	}

	my $seqin = new Bio::SeqIO('-file'=>$input_seq);
	while (my $seq = $seqin->next_seq) {
		$seq{$seq->id} = $seq;
	}
	
	for my $marker (keys %hits) {
		my $seqout = new Bio::SeqIO('-file'=>">>$$.$marker.candidate.pep",'-format'=>'fasta');
		for my $seqid (keys %{$hits{$marker}}) {
			$seqout->write_seq($seq{$seqid});
		}
	}
	return \%hits;
}

sub get_orthologs {
	my $list = shift;

	MARKER:for my $marker (keys %{$list}) {
		unless (-e "$ref_dir/Others/$marker.pep") {
			system("cp $$.$marker.candidate.pep $marker.pep");
			next MARKER;
		}
		my (%score, %tophit, @seqs, %marker) =();
		
		open (IN, "$ref_dir/$marker.stock") || die "Cannot open $ref_dir/$marker.stock";
		while (<IN>) {
			if (/^REF-(\S+)/) {
				$marker{$1} = 1;
			}
		}
		close IN;

		system ("phmmer -Z 5000 -E $evalue --noali --tblout $$.$marker.phmmer -o /dev/null $$.$marker.candidate.pep $ref_dir/Others/$marker.pep");

		open (IN, "$$.$marker.phmmer") || die "Can't open $$.$marker.phmmer";

		while (<IN>) {
			chop;
			next if /^#/;
			my ($target, $target_accession, $query, $query_accession, $evalue, $score, $rest) = split /\s+/;

			if ((! exists $score{$query}) or ($score{$query} < $score)) {
				$score{$query} = $score;
				$tophit{$query} = $target;
			}
		}
		close IN;

		for my $query (keys %score) {
			if (exists $marker{$tophit{$query}}) {
				push @seqs, $seq{$query};
			}
		}
			
		if (@seqs) {
			my $seqout = new Bio::SeqIO('-file'=>">$marker.pep",'-format'=>'fasta');
			for (@seqs) {
				$seqout->write_seq($_);
			}
		}
	}
}

	
