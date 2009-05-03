#!/usr/bin/perl 
use strict;
#####################################################################################################
#PhosCalc is a free script, it is provided with the hope that you will enjoy, you may freely redistribute it at will. We would be greatful if you would keep these acknowledgements with it. 
#
# Put together from info in Olsen and Mann 2004 and Olsen et al 2006 (see provided documentation for further information) 
# by Dan MacLean
# dan.maclean@sainsbury-laboratory.ac.uk
# 
# This program is free academic software; academic and non-profit
# users may redistribute it freely.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
#
#
#LAST EDITED: 23rd April 2009
#
use lib ".";
use Math::Combinatorics;
#####################################################################################################
# Math::Combinatorics is a perl module produced by Allen Day <allenday@ucla.edu>, 
# with algorithmic contributions from Christopher Eltschka and Tye.
# Copyright (c) 2004-2005 Allen Day. All rights reserved. Math::Combinatorics is free software; you can redistribute it and/or modify it under the same terms as Perl itself.
#For more info: http://search.cpan.org/~allenday/Math-Combinatorics-0.08/lib/Math/Combinatorics.pm#AUTHOR
#####################################################################################################
use POSIX qw(log10);
use Getopt::Std;
use Data::Dumper;
##### do brief checks on input:

my %options=();
getopts("p:x:e:o:h:t:m:MAb:D",\%options);

my $vers = "1.2 (Nov 2008)";
my $usage = "\n\n-----PhosCalc $vers -----\nusage = $0 <options> \nFlags:\n -p <peptide spectrum list>\n -t <spectrum type [dta or mgf] (default = dta)>\n -m <mgf file>\n -x <exp type MS2 or MS3 (default = MS2)>\n -e <error margin (default = 0.4)>\n -D <use dehydration>\n -M <use alkylated cysteine (carboxmethylation) = cysteine (103.00919) + 58.00548>\n -A <use alkylated cysteine (carboxamidomethylation) = cysteine (103.00919) + 57.021464>\n -o <output_file>\n -b <best guess range (default = 1000)>\n -h Print this help message\n\n\n";


my $listfile = $options{p} || die "\nERROR: no list file provided\n$usage";
my $experiment_type = 'MS2';
$experiment_type = $options{x} if defined $options{x};
my $spectrum_type = 'dta';
$spectrum_type = $options{t} if defined $options{t};
my $mgf_file = $options{m};
die "\nERROR: unknown input file type $spectrum_type .  Please use \'dta\' or \'mgf\'\n$usage" unless $spectrum_type =~ m/^dta$/i or $spectrum_type =~ m/^mgf$/i;
die "\nERROR: mgf file set as spectrum type but no mgf file provided\n$usage" if $spectrum_type =~ m/^mgf$/i and !$mgf_file;

die "\nERROR: two types of cysteine modification selected. Please select one of:\n carboxmethylation = cysteine (103.00919) + 58.00548 = '-M'\ncarboxamidomethylation = cysteine (103.00919) + 57.021464 = '-A'\n$usage" if $options{M} and $options{A};
my $cysteine = '103.00919';
$cysteine += 58.00548 if $options{M};
$cysteine += 57.021464 if $options{A};

my $window_size = '0.4';
$window_size = $options{e} if defined $options{e};




my $guess_range = '1000';
die "\nERROR: best guess range must be a number\n$usage" if defined $options{b} and $options{b} !~ /^\d+$/; 
$guess_range = $options{b} if defined $options{b};

my $outfile = $options{p} . '_phoscalc_output.tab';
$outfile = $options{o} if defined $options{o};

die $usage if $options{h};

my %peptides;
my $id =0;
my $linenum =0;
open LISTFILE, "<$listfile"  or die "Cant open $listfile\n";
#my $experiment_type = $ARGV[1];
$experiment_type = uc($experiment_type);
#my $window_size = $ARGV[2];
open OUTFILE, ">$outfile" or die "Outfile cannot be opened!\n";

while ( my $line = <LISTFILE> )
{ ++$linenum;
	chomp $line;
	
	my @tmp = split(/\t/,$line);
	$tmp[0] =~ s/\d//g;
	$tmp[0] =~ s/\s//g;
	$tmp[0] = uc($tmp[0]);
	


	if ( check_peptide_sequence($tmp[0]) ) 
	{
	die_because_of_input_problem("error with peptide sequence on line $linenum") ;
	}
	
	#count potential phosphorylation sites
	my $possible_sites = 0;	
	my $number_of_phos_sites = 0;
	for (my $i = 0; $i <= length($tmp[0]); $i++)
	{
		if (substr( $tmp[0], $i, 1) eq '#' or  substr( $tmp[0], $i, 1) eq '@' or substr( $tmp[0], $i, 1) eq '^')
		
		{
			++$number_of_phos_sites;
		}
		if (substr( $tmp[0], $i, 1) eq 'S' or  substr( $tmp[0], $i, 1) eq 'T' or substr( $tmp[0], $i, 1) eq 'Y')
		
		{	
			++$possible_sites;
		}
	}

	unless ($number_of_phos_sites > 0)
	{
		$number_of_phos_sites = 1
	}
	
	if ($possible_sites == 1){

		warn "\n$tmp[0] : .. only one phosphorylation site\n";

	}
	

	my $seq_num = "$tmp[0]$number_of_phos_sites";

	my $concat = "$seq_num:$tmp[1]";
	$peptides{$id} = $concat;
	++$id;	

}
close LISTFILE;

###################################################
#start calculations:                                                                                                           #
###################################################

##set up hash structure of putative ions
my %putative_ions;

##Calculate all possible b and y ions, save in nested hash structure: keys = peptide => ion type => ion => calculated mass
foreach my $peptide (keys %peptides)
{	
	my $concat = $peptides{$peptide};
	
	my @peps = (split/:/,$concat);	
	my $sequence = $peps[0];
	my $spectra_file = $peps[1];
	
	my @variants = get_phospho_site_variants($sequence);

	die "no variants error\n" unless $variants[0];
	foreach my $variant (@variants)
		{
			
		my @b_ions = get_b_ions($variant);
		foreach my $b_ion (@b_ions)
		{
		
			my $mass = calculate_mass($b_ion, 'b', $experiment_type);
			$putative_ions{$peptide}{$variant}{b_ions}{$b_ion} = $mass;
		
		}
	
		my @y_ions = get_y_ions($variant);
		foreach my $y_ion (@y_ions)
		{
		
			my $mass = calculate_mass($y_ion, 'y', $experiment_type); 
			$putative_ions{$peptide}{$variant}{y_ions}{$y_ion} = $mass;
		
		}
	}
}



foreach my $peptide_id (keys %peptides)
{
	
	my @variants = (keys %{$putative_ions{$peptide_id} });

	my $concat = $peptides{$peptide_id};
	
	my @peps = (split/:/,$concat);	
	
	my $sequence = $peps[0];
	chop $sequence;
	my $spectra_file = "$peps[1]";
	
my $ptm_scores = '';
foreach my $variant (@variants)
{
	my %b_ions = %{$putative_ions{$peptide_id}{$variant}{b_ions} };
	my %y_ions = %{$putative_ions{$peptide_id}{$variant}{y_ions} };
	
	
	my %most_intense_peaks = get_spectra_peaks($spectra_file, '1', $spectrum_type, $mgf_file);
	my $charge_state = $most_intense_peaks{charge_state}; #warn "In main $charge_state\n";
	delete $most_intense_peaks{charge_state};
	##add the extra variant ions if the charge state is greater than 2 
	if ($charge_state >= 2 ){
		#warn "Using +2 ions because charge state = $charge_state\n";
		%b_ions = add_two_plus_ions(\%b_ions);
		%y_ions = add_two_plus_ions(\%y_ions);
	}
	if($options{D}){

		%b_ions = add_dehydro_ions(\%b_ions);
		%y_ions = add_dehydro_ions(\%y_ions);

	}
 


	# filter out the putative b and y ions that dont fit in the mass range of the spectra
	#(mass range is stored in %most_intense_peaks keys are mass windows)	
		
	

	my @mass_intervals = keys (%most_intense_peaks);
	@mass_intervals = sort {$a <=> $b} @mass_intervals; 
	my $min_mass = $mass_intervals[0]; 
	my $max_mass = $mass_intervals[-1] + 100;
	
	foreach my $ion (keys %b_ions)
	{
		if ($b_ions{$ion} < $min_mass or $b_ions{$ion} > $max_mass)
		{
			delete $b_ions{$ion}; 
		}
	}
	
	foreach my $ion (keys %y_ions)
	{
		if ($y_ions{$ion} < $min_mass or $y_ions{$ion} >$max_mass)
		{
			delete $y_ions{$ion};
		}
	}
	
	die "no ions within range $min_mass .. $max_mass for file $spectra_file\n " unless keys %b_ions and keys %y_ions;
	# calculate the number of matches between putative ions and actual measured masses
	
	my @matches = identify_matches(\%b_ions, \%y_ions, \%most_intense_peaks, $charge_state, $window_size);
	

	#reformat the sequences for people reading
	$variant =~ s/1/[S]/g;
	$variant =~ s/2/[T]/g;
	$variant =~ s/3/[Y]/g;
	
	#$variant =~ s/W//gi;
	$variant =~ s/[\#@\*\^]//gi;
	my $sequence = $peptides{$peptide_id};

	my @seq = split(/\:/,$sequence);

	#chop $sequence;	
	$seq[0] =~ s/\W//gi; 
	$seq[0] =~ s/[\#@\*\^\d]//gi;
	#unless (exists $failed_ids{$peptide_id} )


	#my @seq = split(/\:/, $sequence); 
	#my @no_b_ions = get_b_ions($seq[0]);
	#my @no_y_ions = get_y_ions($seq[0]);
	#my %summ_ions;
	#foreach (@no_b_ions){ $summ_ions{$_} = 1;}
	#foreach (@no_y_ions){ $summ_ions{$_} = 1;}
	#my $total_no_ions = keys %summ_ions;
	my $total_no_ions = scalar(keys %b_ions) + scalar(keys %y_ions);
	my @probabilities = calculate_scores($total_no_ions,$matches[1]);
	


	######OUTPUT#########
	{
	##prepare for best guess... 
	if ($ptm_scores =~ m/\w/){
		$ptm_scores = $ptm_scores . $variant . ':' . (sprintf '%.2e',"$probabilities[0]") . ';';

	}
	else {
	      $ptm_scores = $variant . ':' . (sprintf '%.2e',"$probabilities[0]") . ';';

	} 


	print OUTFILE "$variant\t$total_no_ions\t$matches[1]\t$probabilities[0]\t$probabilities[1]";

	###special case where number matched > than ions: report case to user

	if ($matches[1] >= $total_no_ions){

	my @special_probs = calculate_scores($total_no_ions, $total_no_ions);
	print OUTFILE "\tSpecial Case: Matched Ions > Total Ions\t Maximum probability for this sequence ($total_no_ions / $total_no_ions) = \t$special_probs[0]\t Maximum score for this sequence ($total_no_ions / $total_no_ions) = \t$special_probs[1]";
	}

	print OUTFILE "\n";	
	
	}
	
}
	#if ( exists $options{b} ){
	#warn $guess_range, "\n";
	my $best_guess = make_format_string($ptm_scores, $guess_range);
	my $peptide_info = $peptides{$peptide_id};
	my @info = split(/\:/,$peptide_info);

	#chop $sequence;	
	$info[0] =~ s/[\W\d]//gi;
	$info[0] =~ s/[\#@\*\^]//gi;
	
	print OUTFILE "BEST GUESS:\t$info[0]\t$info[1]\t$best_guess\n\\\\\n";
	#}
}

close OUTFILE;
print "written to $outfile\n\n";
##################################################
#sub-routines                                                                                                                
##################################################


sub add_dehydro_ions {
my %ions = %{$_[0]};

		foreach my $ion (keys %ions){
			
			if($ion =~ m/_2plus/){

				my $mass = $ions{$ion} - (18.01057/2);
				my $ion_dehydro = $ion . '_dehydro';
				$ions{$ion_dehydro} = $mass;

			}
			else{
				
				my $mass = $ions{$ion} - 18.01057;
				my $ion_dehydro = $ion . '_dehydro';
				$ions{$ion_dehydro} = $mass;

			}

		} 

return %ions;


}

sub add_two_plus_ions{
my %ions = %{$_[0]};

		foreach my $ion (keys %ions){
			
			my $ion_2_plus = $ion . '_2plus';
			my $mass = ($ions{$ion} + (length($ion) *  1.0078) ) / 2;
			$ions{$ion_2_plus} = $mass;

		} 

return %ions;
}
sub calculate_scores{

my $n = shift; 
my $k = shift; 

my $n_fac = factorial($n); 
my $k_fac = factorial($k); 

my $n_k = $n - $k; 

my $n_k_fac = factorial($n_k); 

my $c = 0.96 ** ($n_k); 
my $b = 0.04 ** $k; 
my $a = $n_fac / ( $k_fac *  $n_k_fac ); 
my $p = $a * $b * $c;

my $score = -10 * ( log10($p) );

return ($p, $score);
}

sub factorial{
my $counter  = shift;
my $factorial = 1;
while ( $counter >=1 ) 
{
	$factorial *= $counter--; 
}
return $factorial;
}


##################################################

sub identify_matches{

my %b_ions = %{$_[0]};
my %y_ions = %{$_[1]};
my %most_intense_peaks = %{$_[2]};
my $charge_state = $_[3];
my $error_margin = "$_[4]";  
unless (defined $error_margin) {$error_margin = '0.4';}

my @total_b_ions = keys %b_ions;
my @total_y_ions = keys %y_ions;
my $total_b = @total_b_ions;
my $total_y = @total_y_ions;
my $total_ions = $total_b + $total_y; 

my $matched_ions = 0;

foreach my $b_ion (keys %b_ions)
{
	foreach my $window (keys %most_intense_peaks)
	{
		my @mass = split(/:/, $most_intense_peaks{$window});
		foreach my $intensity_mass (@{$most_intense_peaks{$window} })
		{
			my @mass = split(/:/, $intensity_mass);
			if ( ($mass[1] + $error_margin)  >= $b_ions{$b_ion}  and  ( $mass[1] - $error_margin ) <=$b_ions{$b_ion})
			{
			++$matched_ions; #warn "matched an ion!\n";
			}
			else { }#warn "didnt match -> ", "\t", $mass[0], "\t", $mass[1], "\n";} 
		}
	}
}



foreach my $y_ion (keys %y_ions)
{
	foreach my $window (keys %most_intense_peaks)
	{
		my @mass = split(/:/, $most_intense_peaks{$window});
		foreach my $intensity_mass (@{$most_intense_peaks{$window} })
		{
			my @mass = split(/:/, $intensity_mass);
			if ( ($mass[1] + $error_margin)  >= $y_ions{$y_ion}  and  ( $mass[1] - $error_margin ) <=$y_ions{$y_ion})
			{
			++$matched_ions;
			}
		}
	}
}
return ($total_ions, $matched_ions);

}   
##################################################

sub get_phospho_site_variants{

my $unedited_sequence = shift;  
my $num_phos_sites = chop($unedited_sequence); 
my @phos_site_indices;
my @chars = split(undef, $unedited_sequence);


my @variants;	
	my $i = 0;
	foreach my $char (@chars)
	{	
	if ($char eq 'S' or $char eq 'T' or $char eq 'Y')
	{
		push @phos_site_indices, $i;
	}
	++$i;
	}
	
my @combinations = combine($num_phos_sites,@phos_site_indices);

	foreach my $combination (@combinations)
	{ 
	my @one_combination = @$combination;
	my $new_seq = $unedited_sequence;	
	foreach my $residue (@one_combination)
	{
	
	if (substr($unedited_sequence, $residue, 1) eq 'S' ) 
		{
		my @chars = split(undef, $new_seq);
		$chars[$residue] = '1';
		$new_seq =join('',@chars);
		}

		if (substr($unedited_sequence, $residue, 1) eq 'T' ) 
		{
		my @chars = split(undef, $new_seq);
		$chars[$residue] = '2';
		$new_seq =join('',@chars);
		}
		if (substr($unedited_sequence, $residue, 1) eq 'Y' ) 
		{
		my @chars = split(undef, $new_seq);
		$chars[$residue] = '3';
		$new_seq =join('',@chars);
		}
	
	}	
	push(@variants, $new_seq); 
	}

	
	
	
return @variants;
}

##################################################
sub get_spectra_peaks{

#$spectra_file, $filter, $spectrum_type, $mgf_file

my $spectrum = shift; #"phos_tool_uploaded_data/$upload_file";
my $filter = shift;
my $spectrum_type = shift;
my $mgf_file = shift || 'null';
my $charge_state = 0;

my %spectra_info;
my $i = 0;





if($spectrum_type =~ m/^dta$/i){
	#warn "Using a .dta file $spectrum\n";
	open FILE, "<$spectrum" or die "Couldn't open $spectrum .. !\n";
	while (my $line = <FILE>)
	{	
	#ignore first line
	if ($i == '0')
	{
	chomp $line;
	my @tmp = split(/\s/,$line);
	$charge_state = $tmp[1];
	}
	else
	{
	chomp $line;
	my @tmp = split(/\s/,$line);
	$spectra_info{$tmp[0]} = $tmp[1];
	}
	++$i;
}


close FILE;
}
elsif($spectrum_type =~ m/^mgf$/i){
	my $found = 0;
	local $/ = 'END IONS';
	#warn "Using a .mgf file $mgf_file\n";
	open FILE, "<$mgf_file" or die "Couldn't open $mgf_file .. !\n";

	#warn "spec is ", $spectrum, "\n";
	while (my $record = <FILE>){

	if ($record =~ m/TITLE=$spectrum/){
		$found = 1;
		$record =~ m/CHARGE=(\d)/; 
		$charge_state = $1; #warn "In sub $charge_state\n";

		while ($record =~ m/(\d+\.\d+\s\d+\.\d+)/g){

#			warn $1, "\n";
			my @tmp = split(/\s/, $1);
			#warn $tmp[0], "\n";
			$spectra_info{$tmp[0]} = $tmp[1];
			
		}


		

	}

	}
	close FILE;
	
	
	die "Couldn't find spectrum $spectrum in file $mgf_file .. \n" unless $found;

}



my @sorted_masses = sort {$a <=> $b} keys %spectra_info; 		
my $min_reported_mass = $sorted_masses[0];
my $max_reported_mass = $sorted_masses[-1];

my %windows;

for (my $i = $min_reported_mass;  $i <= ($max_reported_mass + 101); $i = $i + 100)
{
	foreach my $mass (keys %spectra_info)
	{
		if ($mass >= $i and $mass <= $i + 100)
		{	my $combined_mass_intensity = "$spectra_info{$mass}:$mass"; 
			push @{$windows{$i} }, $combined_mass_intensity;
		}
	}
}

####sort out the four most intense values for each window 
if ($filter){
foreach my $window (keys %windows)
{	
	my %intensity_mass;
	my $i = 0;
	my @combined = @{ $windows{$window} };
	foreach my $el (@combined)
	{
		my @tmp = split(/:/, $el);
		$intensity_mass{$el} = $tmp[0];
	}
	my @most_intense;
	my @sorted_intensities = sort {$b <=> $a} values %intensity_mass;
	@sorted_intensities = @sorted_intensities[0 .. 3 ];
	foreach my $el (@sorted_intensities)
	{
		foreach my $key (keys %intensity_mass)
		{
			my @tmp = split(/:/, $key);
			if ($el eq $tmp[0])
			{
			push (@most_intense, $key);
			}
		}
	}

	 @{ $windows{$window} } = @most_intense[0 .. 3];
	
}
}
else{

}
$windows{charge_state} = $charge_state;
return %windows;

}

##################################################


sub calculate_mass{
####args to this subroutine are 0 = sequence, 1 = type of ion, either b,a,y,c or i


my $ion = shift;
my $type_of_ion = shift;
my $experiment_type = shift; ### get the experiment type from cgi form
my $mass = 0;
my %residue_masses;
if ($experiment_type eq 'MS2'){

%residue_masses =
(
       
    G    =>    '57.02147'    ,
    A    =>    '71.03712'    ,
    S    =>    '87.03203'    ,
    P    =>    '97.05277'    ,
    V    =>    '99.06842'    ,
    T    =>    '101.04768'    ,
    C    =>    $cysteine    ,
    L    =>    '113.08407'    ,
    I    =>    '113.08407'    ,
    N    =>    '114.04293'    ,
    D    =>    '115.02695'    ,
    Q    =>    '128.05858'    ,
    K    =>    '128.09497'    ,
    E    =>    '129.04260'    ,
    M    =>    '131.04049'    ,
    H    =>    '137.05891'    ,
    F    =>    '147.06842'    ,
    R    =>    '156.10112'    ,
    Y    =>    '163.06333'    ,
    W    =>    '186.07932'    ,
    #extra symbols representing phosphorylated residues from ms2
    1    =>   '167.03203'      ,  # S + +79.9663
    2    =>   '181.04768'    ,   # T + +79.9663
    3    =>    '243.02963'    ,    # Y + +79.9663
    4      =>    '147.03541'        # M* + 15.9949
 );

for (my $i = 0; $i < length($ion); $i++)
{

	$mass = $mass + ( $residue_masses{substr($ion, $i, 1)} );
	if (substr($ion, $i, 2) eq 'M*' ){$mass = $mass + 15.9949;}
}



if ($type_of_ion eq 'b') {$mass = $mass + 1.0078;}
elsif ($type_of_ion eq 'y'){$mass = $mass + 19.0184;}
}




elsif ($experiment_type eq 'MS3'){
%residue_masses =
(
       
    G    =>    '57.02147'    ,
    A    =>    '71.03712'    ,
    S    =>    '87.03203'    ,
    P    =>    '97.05277'    ,
    V    =>    '99.06842'    ,
    T    =>    '101.04768'    ,
    C    =>    $cysteine    ,
    L    =>    '113.08407'    ,
    I    =>    '113.08407'    ,
    N    =>    '114.04293'    ,
    D    =>    '115.02695'    ,
    Q    =>    '128.05858'    ,
    K    =>    '128.09497'    ,
    E    =>    '129.04260'    ,
    M    =>    '131.04049'    ,
    H    =>    '137.05891'    ,
    F    =>    '147.06842'    ,
    R    =>    '156.10112'    ,
    Y    =>    '163.06333'    ,
    W    =>    '186.07932'    ,
    #extra symbols representing phosphorylated residues
    1    =>   '69.01703'      ,  # S - 18.0105
    2    =>   '83.03268'    ,   # T - 18.0105
    3    =>    '163.06333',        # Y
    4      =>    '147.03541'        # M* + 15.9949

 );


for (my $i = 0; $i < length($ion); $i++)
{

	$mass = $mass + ( $residue_masses{substr($ion, $i, 1)} );
	if (substr($ion, $i, 2) eq 'M*' ){$mass = $mass + 15.9949;}

}
if ($type_of_ion eq 'b') {$mass = $mass + 1.0079;}
elsif ($type_of_ion eq 'y'){$mass = $mass + 19.0184;}



}

return $mass

}

##################################################
sub get_sequence {

	my $peptide_id = shift;
	my @y_ions;
	my $sequence = ##get sequence from cgi form##
	my $number_of_phos_sites;

	for (my $i = 0; $i <= length($sequence); $i++)
	{
		if (substr( $sequence, $i, 1) eq '#' or  substr( $sequence, $i, 1) eq '@' or substr( $sequence, $i, 1) eq '^')
		
		{
			++$number_of_phos_sites;
		}
	}
		
my $concat = "$sequence$number_of_phos_sites";
return $concat;

}



##################################################
sub get_b_ions {
	my $sequence = shift;
	my @b_ions;
		for (my $i = 1; $i < length($sequence); $i++)
		{
			push(@b_ions, substr($sequence, 0, $i))
		}

	return @b_ions;

}

##################################################
sub get_y_ions {
	my $sequence = shift;
	my @y_ions;	
		for (my $i = 1; $i < length($sequence); $i++)
		{
			push(@y_ions, substr($sequence, (length($sequence) - $i), length($sequence)));
		}

		

		
	@y_ions = reverse @y_ions;
	return @y_ions;

}

##################################################
sub die_because_of_input_problem {

my $reason=shift;
print "Input error: $reason\n";
exit;
}

###################################################
sub check_input_file{
my $found_something_wrong = 0;
my $filename = shift;

open FILE, "<$filename" || die "can't open $filename\n";
while (my $line = <FILE>)
{
	chomp $line;
	my @tmp = split(/\s/,$line);
	if($tmp[0] =~ m/[^\d\.]/g or $tmp[1] =~ m/[^\d\.]/g)
	{
		$found_something_wrong = 1;
	}
}
close FILE;
return $found_something_wrong;
}
###################################################
sub check_peptide_sequence{
my $se = shift;
my $something_wrong = 0;
$se =~ s/\d//g;
$se =~ s/\s//g;
if ($se =~ m/[^ACDEFGHIKLMNPQRSTVWY@\#^*]/gi)
{ 
	
	$something_wrong = 1;

}
return $something_wrong;
}
#####################################################
sub make_format_string{
	my %scores;
	my $lowest;
	my $passed;
	my $ptm_scores = shift;
	my $range = shift;
	$ptm_scores =~ s/\s+//g;
	my @sc = split(/\;/, $ptm_scores);
	foreach my $s (@sc){

		my @i = split(/\:/, $s);

		my $dec = expand($i[1]);
		#warn $dec, "\n"; sleep 1;
		$scores{$i[0]} = $dec;
		
	}

	foreach my $value (sort {$scores{$a} cmp $scores{$b} }
           keys %scores){
		if (!defined $lowest){
			$lowest = $value;
			$passed = $value;
			#warn "$lowest is highest\n$scores{$lowest}\n";
			my $min = $scores{$lowest} * $range;
			#warn "min is ", $min, "\n";

		}
		elsif(defined $lowest and ($scores{$value} <= $scores{$lowest} * $range)){

			$passed = $passed . ':' . $value;			

		}

	

	}

	my @passes = split(/\:/, $passed);

	if ($passes[1]){
		
		my %posns;
		foreach my $pass (@passes){
		#warn $pass, "\n"; sleep 2;
		$pass =~ s/[\#@\*\^]//g;
		while($pass =~ m/(\[)/){
			my $match = $`; #`
			my $pos = length($match);
			#warn $pos, " is pos\n";
			$posns{$pos} = 1;
			$pass =~ s/\[//;
			$pass =~ s/\]//;
			#warn "added $pos\n"; 
			}
		}

		#foreach my $key (keys %posns){

		#	warn $key, "=> is key\n";

		#}

		$passes[0] =~ s/[\[\]]//;;

		#warn $passes[0], "\n";

		
		my @res = split(/|/, $passes[0]);
		my $format_string = '';
		for (my $i = 0; $i < length($passes[0]); $i++){

		if (exists $posns{$i}){

			$format_string = $format_string . '[' . lc($res[$i]) . ']';
			
		}
		else {

			$format_string = $format_string . $res[$i];

		}


		}

		return $format_string;

	}
	else{
		my %posns; 
		$passes[0] =~ s/[\#@\*\^]//g;
		$passes[0] =~ m/(\[)/;
		my $match = $`; #`
		my $pos = length($match);
		$posns{$pos} = 1;
		$passes[0] =~ s/[\[\]]//g;
		my @res = split(/|/, $passes[0]);
		my $format_string = '';
		for (my $i = 0; $i < length $passes[0]; $i++){

		if (exists $posns{$i}){

			$format_string = $format_string . '[p' . $res[$i] . ']' ;

		}
		else {

			$format_string = $format_string . $res[$i];

		}


		}

		return $format_string;

	}

	
}
####################################################################

sub expand {
        my $n = shift;
        return $n unless $n =~ /^(.*)e([-+]?)(.*)$/;
        my ($num, $sign, $exp) = ($1, $2, $3);
        my $sig = $sign eq '-' ? "." . ($exp - 1 + length $num) : '';
        return sprintf "%${sig}f", $n;
}
