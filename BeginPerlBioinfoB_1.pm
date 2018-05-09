#poly_sites_syn
#a subroutine to return an array of positions of synonymous polymorphic sites
sub poly_sites_syn
{  my $gap=0;
   my @align =@_;
    length ($align[0]);
    my($position);
    my($comp);
    my($count) = 0;
   my @polysites=();


 for ($position=0; $position<(length($align[0])-2); $position+=3)
{ 

 $gap=0;
 
	for (my $seq=0; $seq<scalar (@align); ++$seq)
{
for ($j=0; $j<3; ++$j)
{

if ((substr ($align[$seq], ($position+$j), 1) ne 'A') and (substr ($align[$seq], ($position+$j), 1) ne 'G') and (substr ($align[$seq], ($position+$j), 1) ne 'T') and (substr ($align[$seq], ($position+$j), 1) ne 'C'))
{$seq=scalar(@align);
 $gap=1;}
}}
 if ($gap ne 1)
{
 for ($w=0; $w<3; ++$w)
      {  for ($comp=1; $comp<(scalar(@align)); ++$comp)
{

    $base1=substr($align[0],($position+$w),1);
    $base2=substr($align[$comp],($position+$w),1);


if ($base1 ne $base2)
{
$codon1=substr($align[0], $position, 3);


if ($w eq 0)
{
$codon2=substr($align[$comp],$position,1).substr($align[0],$position+1,1).substr($align[0],$position+2,1);
}
elsif ($w eq 1)
{$codon2=substr($align[0],$position,1).substr($align[$comp],$position+1,1).substr($align[0],$position+2,1);
}
elsif ($w eq 2)
{$codon2=substr($align[0],$position,1).substr($align[0],$position+1,1).substr($align[$comp],$position+2,1);}


$amino1=codon2aa($codon1);
$amino2=codon2aa($codon2);

#print "$position $comp $codon1 $codon2 $amino1 $amino2\n";


if ($amino1 eq $amino2)
{


    $polysites[$count]=$position+$w;
$count+=1;
 

$w+=1;     
$comp=scalar[@align];

}

}}
}
}
}
return @polysites;

}


#poly_sites_non
#a subroutine to return an array of positions of nonsynonymous polymorphic sites
sub poly_sites_non
{  my $gap=0;
   my @align =@_;
    length ($align[0]);
    my($position);
    my($comp);
    my($count) = 0;
   my @polysites=();


 for ($position=0; $position<(length($align[0])-2); $position+=3)
{ 

 $gap=0;
 
	for (my $seq=0; $seq<scalar (@align); ++$seq)
{
for ($j=0; $j<3; ++$j)
{

if ((substr ($align[$seq], ($position+$j), 1) ne 'A') and (substr ($align[$seq], ($position+$j), 1) ne 'G') and (substr ($align[$seq], ($position+$j), 1) ne 'T') and (substr ($align[$seq], ($position+$j), 1) ne 'C'))
{$seq=scalar(@align);
 $gap=1;}
}}
 if ($gap ne 1)
{for ($w=0; $w<3; ++$w)

      { for ($comp=1; $comp<(scalar(@align)); ++$comp)

{

    $base1=substr($align[0],($position+$w),1);
    $base2=substr($align[$comp],($position+$w),1);


if ($base1 ne $base2)
{
$codon1=substr($align[0], $position, 3);


if ($w eq 0)
{
$codon2=substr($align[$comp],$position,1).substr($align[0],$position+1,1).substr($align[0],$position+2,1);
}
elsif ($w eq 1)
{$codon2=substr($align[0],$position,1).substr($align[$comp],$position+1,1).substr($align[0],$position+2,1);
}
elsif ($w eq 2)
{$codon2=substr($align[0],$position,1).substr($align[0],$position+1,1).substr($align[$comp],$position+2,1);}


$amino1=codon2aa($codon1);
$amino2=codon2aa($codon2);

#print "$position $comp $codon1 $codon2 $amino1 $amino2\n";


if ($amino1 ne $amino2)
{

$comp=scalar(@align);
    $polysites[$count]=$position+$w;
$count+=1;
 

$w+=1;     
}

}}
}
}
}
return @polysites;

}




sub read_dir{    $dir=$_[0];    $pattern=$_[1];    opendir (DIR, "$dir") || die " can not read directory '$dir' $!";    @files = sort grep(/$pattern/, readdir(DIR));    closedir DIR;    return @files;}
# From Chapter 8

#
# codon2aa
#
# A subroutine to translate a DNA 3-character codon to an amino acid
#   Version 3, using hash lookup

sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '_',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{

            print  "Bad codon \"$codon\"!!\n";
            
    }
}


# From Chapter 8

#
# codonprefs
#
# A subroutine to test whether a codon is preferred ('P'), unpreferred #('U'), or not defined ('N')

sub codonprefs {
    my($codonpref) = @_;

    $codonpref = uc $codonpref;
 
    my(%codon_prefs) = (
    
    'TCA' => 'U',    # Serine
    'TCC' => 'P',    # Serine
    'TCG' => 'U',    # Serine
    'TCT' => 'P',    # Serine
    'TTC' => 'N',    # Phenylalanine
    'TTT' => 'N',    # Phenylalanine
    'TTA' => 'N',    # Leucine
    'TTG' => 'N',    # Leucine
    'TAC' => 'N',    # Tyrosine
    'TAT' => 'N',    # Tyrosine
    'TAA' => 'N',    # Stop
    'TAG' => 'N',    # Stop
    'TGC' => 'N',    # Cysteine
    'TGT' => 'N',    # Cysteine
    'TGA' => 'N',    # Stop
    'TGG' => 'N',    # Tryptophan
    'CTA' => 'U',    # Leucine
    'CTC' => 'P',    # Leucine
    'CTG' => 'U',    # Leucine
    'CTT' => 'P',    # Leucine
    'CCA' => 'N',    # Proline
    'CCC' => 'N',    # Proline
    'CCG' => 'N',    # Proline
    'CCT' => 'N',    # Proline
    'CAC' => 'N',    # Histidine
    'CAT' => 'N',    # Histidine
    'CAA' => 'N',    # Glutamine
    'CAG' => 'N',    # Glutamine
    'CGA' => 'N',    # Arginine
    'CGC' => 'N',    # Arginine
    'CGG' => 'N',    # Arginine
    'CGT' => 'N',    # Arginine
    'ATA' => 'N',    # Isoleucine
    'ATC' => 'N',    # Isoleucine
    'ATT' => 'N',    # Isoleucine
    'ATG' => 'N',    # Methionine
    'ACA' => 'U',    # Threonine
    'ACC' => 'P',    # Threonine
    'ACG' => 'U',    # Threonine
    'ACT' => 'P',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'N',    # Lysine
    'AAG' => 'N',    # Lysine
    'AGC' => 'N',    # Serine
    'AGT' => 'N',    # Serine
    'AGA' => 'N',    # Arginine
    'AGG' => 'N',    # Arginine
    'GTA' => 'U',    # Valine
    'GTC' => 'P',    # Valine
    'GTG' => 'U',    # Valine
    'GTT' => 'P',    # Valine
    'GCA' => 'N',    # Alanine
    'GCC' => 'N',    # Alanine
    'GCG' => 'N',    # Alanine
    'GCT' => 'N',    # Alanine
    'GAC' => 'N',    # Aspartic Acid
    'GAT' => 'N',    # Aspartic Acid
    'GAA' => 'N',    # Glutamic Acid
    'GAG' => 'N',    # Glutamic Acid
    'GGA' => 'U',    # Glycine
    'GGC' => 'P',    # Glycine
    'GGG' => 'U',    # Glycine
    'GGT' => 'P',    # Glycine
    );

    if(exists $codon_prefs{$codonpref}) {
        return $codon_prefs{$codonpref};
    }else{

            print  "Bad codon \"$codonpref\"!!\n";
            
    }
}


# From Chapter 8

#
# codonprefs_ATneut
#
# A subroutine to test whether a codon is preferred ('P'), unpreferred #('U'), or not defined ('N')

sub codonprefs_ATneut {
    my($codonpref) = @_;

    $codonpref = uc $codonpref;
 
    my(%codon_prefs) = (
    
    'TCA' => 'U',    # Serine
    'TCC' => 'P',    # Serine
    'TCG' => 'U',    # Serine
    'TCT' => 'P',    # Serine
    'TTC' => 'N',    # Phenylalanine
    'TTT' => 'N',    # Phenylalanine
    'TTA' => 'N',    # Leucine
    'TTG' => 'N',    # Leucine
    'TAC' => 'N',    # Tyrosine
    'TAT' => 'N',    # Tyrosine
    'TAA' => 'N',    # Stop
    'TAG' => 'N',    # Stop
    'TGC' => 'N',    # Cysteine
    'TGT' => 'N',    # Cysteine
    'TGA' => 'N',    # Stop
    'TGG' => 'N',    # Tryptophan
    'CTA' => 'U',    # Leucine
    'CTC' => 'P',    # Leucine
    'CTG' => 'U',    # Leucine
    'CTT' => 'P',    # Leucine
    'CCA' => 'N',    # Proline
    'CCC' => 'N',    # Proline
    'CCG' => 'N',    # Proline
    'CCT' => 'N',    # Proline
    'CAC' => 'N',    # Histidine
    'CAT' => 'N',    # Histidine
    'CAA' => 'N',    # Glutamine
    'CAG' => 'N',    # Glutamine
    'CGA' => 'N',    # Arginine
    'CGC' => 'N',    # Arginine
    'CGG' => 'N',    # Arginine
    'CGT' => 'N',    # Arginine
    'ATA' => 'N',    # Isoleucine
    'ATC' => 'N',    # Isoleucine
    'ATT' => 'N',    # Isoleucine
    'ATG' => 'N',    # Methionine
    'ACA' => 'U',    # Threonine
    'ACC' => 'P',    # Threonine
    'ACG' => 'U',    # Threonine
    'ACT' => 'P',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'N',    # Lysine
    'AAG' => 'N',    # Lysine
    'AGC' => 'N',    # Serine
    'AGT' => 'N',    # Serine
    'AGA' => 'N',    # Arginine
    'AGG' => 'N',    # Arginine
    'GTA' => 'U',    # Valine
    'GTC' => 'P',    # Valine
    'GTG' => 'U',    # Valine
    'GTT' => 'P',    # Valine
    'GCA' => 'N',    # Alanine
    'GCC' => 'N',    # Alanine
    'GCG' => 'N',    # Alanine
    'GCT' => 'N',    # Alanine
    'GAC' => 'N',    # Aspartic Acid
    'GAT' => 'N',    # Aspartic Acid
    'GAA' => 'N',    # Glutamic Acid
    'GAG' => 'N',    # Glutamic Acid
    'GGA' => 'U',    # Glycine
    'GGC' => 'P',    # Glycine
    'GGG' => 'U',    # Glycine
    'GGT' => 'P',    # Glycine
    );

    if(exists $codon_prefs{$codonpref}) {
        return $codon_prefs{$codonpref};
    }else{

            print  "Bad codon \"$codonpref\"!!\n";
            
    }
}

#calculate_synonymous_GC_changes
#subroutine to calculate, for a pair of sequences, the number of times sequence 1 has a synonymous A/T and sequence 2 has a synonymous G/C
sub calculate_synonymous_GC_changes
{ $pol=0;
  @site=();
  $site[0]=0;
  $site[1]=0;
    my (@alignment) = @_;
    for ($i=0; $i<(length($alignment[0])-2); $i+=3)
{ $codon1 = substr($alignment[0],$i,3);
  $codon2 = substr($alignment[1],$i,3);

  $test=0;
  for ($w=0; $w<2; ++$w)
{for($x=0; $x<3; ++$x)
{if ((substr($codon1, $x, 1) ne "A") and (substr($codon1, $x, 1) ne "T") and (substr($codon1, $x, 1) ne "C") and (substr($codon1, $x, 1) ne "G"))
{ $test+=1;} 
}}

  for ($w=0; $w<2; ++$w)
{for($x=0; $x<3; ++$x)
{if ((substr($codon2, $x, 1) ne "A") and (substr($codon2, $x, 1) ne "T") and (substr($codon2, $x, 1) ne "C") and (substr($codon2, $x, 1) ne "G"))
{ $test+=1;} 
}}





if ($test eq 0)
{
if (($codon1 ne $codon2) and (codon2aa($codon1) eq codon2aa($codon2)))
{ if ((substr ($codon1, 2, 1) ne "T") and (substr($codon1, 2, 1) ne "A") and (substr($codon2, 2, 1) ne "G") and (substr($codon2,2,1) ne "C"))

{ $site[0]+=1}

elsif ((substr ($codon2, 2, 1) ne "T") and (substr($codon2, 2, 1) ne "A") and (substr($codon1, 2, 1) ne "G") and (substr($codon1,2,1) ne "C"))
{$site[1]+=1;}


}


}}

return @site;

}

#count_nucs
#a subroutine to return an array of nucleotide counts of a DNA sequnece
sub count_nucs
{my ($seq)=@_;


$nucs[0]=($seq =~ tr/Aa//);
$nucs[1]=($seq =~ tr/Tt//);
$nucs[2]=($seq =~ tr/Cc//);
$nucs[3]=($seq =~ tr/Gg//);
return @nucs;
}
#calculate_codonpref_diffs
#subroutine to calculate, for a pair of sequences, the number of times sequence 1 has a preferred codon and sequence 2 has an unpreferred codon, and vice versa
sub calculate_codonpref_diffs
{ $pol=0;
  @site=();
  $site[0]=0;
  $site[1]=0;
    my (@alignment) = @_;
    for ($i=0; $i<(length($alignment[0])-2); $i+=3)
{ $codon1 = substr($alignment[0],$i,3);
  $codon2 = substr($alignment[1],$i,3);

  $test=0;
  for ($w=0; $w<2; ++$w)
{for($x=0; $x<3; ++$x)
{if ((substr($codon1, $x, 1) ne "A") and (substr($codon1, $x, 1) ne "T") and (substr($codon1, $x, 1) ne "C") and (substr($codon1, $x, 1) ne "G"))
{ $test+=1;} 
}}

  for ($w=0; $w<2; ++$w)
{for($x=0; $x<3; ++$x)
{if ((substr($codon2, $x, 1) ne "A") and (substr($codon2, $x, 1) ne "T") and (substr($codon2, $x, 1) ne "C") and (substr($codon2, $x, 1) ne "G"))
{ $test+=1;} 
}}





if ($test eq 0)
{
if (($codon1 ne $codon2) and (codon2aa($codon1) eq codon2aa($codon2)))
{ if ((codonprefs_ATneut($codon1) eq "P") and (codonprefs_ATneut($codon2) eq "U"))

{ $site[0]+=1}

elsif ((codonprefs_ATneut($codon1) eq "U") and (codonprefs_ATneut($codon2) eq "P"))
{$site[1]+=1;}


}


}}

return @site;

}














#get_file_data
#subroutine to get data from a file given its filename
sub get_file_data {
    my($filename) = @_;
    use strict;
   
    my @filedata = ();
    unless ( open(GET_FILE_DATA, $filename) ) {
	print STDERR "Cannot open file \"$filename\"\n\n";
	exit;
    }
    @filedata = <GET_FILE_DATA>;
    close GET_FILE_DATA;
    return @filedata;
}


#extract_names_from_fasta_data
#subroutine that stores the names of the sequences of a fasta file
sub extract_names_from_fasta_data
{my (@fasta_file_data) =@_;
 use strict;
 my $i=-1;
 my @names;
foreach my $line (@fasta_file_data)
{ if ($line =~ /^>/)
  {$i=$i+1;
 chomp($line);
  $names[$i]= $line;
   
   
 }
}
   return @names;
}

#sizeof_fasta
#subroutine to determine number of sequences in fasta file
sub sizeof_fasta
{ my(@fasta_file_data) = @_;
    use strict;
    my $i=0;
    my @sequence;
    foreach my $line (@fasta_file_data) {
	#discard blank lines
       
	
	if ($line =~ /^>/){
	    $i +=1;
	  
	    }
    }
	return $i;
}



#extract_sequence_from_fasta_data
# A subroutine to extract FASTA sequence data from an array
sub extract_sequence_from_fasta_data {
    my(@fasta_file_data) = @_;
    use strict;
    my $i=-1;
    my @sequence;
    foreach my $line (@fasta_file_data) {
	#discard blank lines
	if ($line =~ /^\s*$/) {
	    next;
	}
	elsif ($line =~ /^>/){
	    $i +=1;
	    next;
	    


}

	else { $sequence[$i] .= $line;
	   }
    }
for (my $j=0; $j<($i+1); ++$j)
{   
    $sequence[$j] =~ s/\s//g;}
    return @sequence;
}



#extract_fullsequence_from_fasta_data
# A subroutine to extract FASTA sequence data from an array, including names
sub extract_fullsequence_from_fasta_data {
    my(@fasta_file_data) = @_;
    use strict;
    my $i=-1;
    my @sequence;
    foreach my $line (@fasta_file_data) {
       
       if ($line =~ /^>/){
	    $i +=1;
$sequence[$i] .= $line;	    

next;
	    


}

	else { $sequence[$i] .= $line;
	   }
    }

    return @sequence;
}

#extract_qualsequence_from_fasta_data
# A subroutine to extract FASTA quality data from an array-doesn't remove spaces!
sub extract_qualsequence_from_fasta_data {
    my(@fasta_file_data) = @_;
    use strict;
    my $i=-1;
    my @sequence;
    foreach my $line (@fasta_file_data) {
	#discard blank lines
	if ($line =~ /^\s*$/) {
	    next;
	}
	elsif ($line =~ /^>/){
	    $i +=1;
	    next;
	    }

	else { 
	    chomp ($line);
$sequence[$i] .= $line;
	   }
    }

    return @sequence;

}

#print_sequence
# A subroutine to format and print sequence data

sub print_sequence {
    my($sequence, $length) = @_;

    use strict;
   

    #print sequence in lines of $length

    for ( my $pos = 0; $pos < length($sequence); $pos +=$length) {
	print substr($sequence, $pos, $length), "\n";
    }
}


#gapped_sites
#a program to record positions nucleotide sites that have gaps or nonstandard 
#bases
sub gapped_sites
{
    @align =@_;
    my($length) = length ($align[0]);
    my($position);
    my($comp);
    my($count) = 0;
    my @gapsites;
 for ($position=0; $position < $length ; ++$position)
    {

 for ($comp=0; $comp<(scalar(@align)); ++$comp)
      { 
if ((substr ($align[$comp], $position, 1) ne "a") and (substr ($align[$comp], $position, 1) ne "g") and (substr ($align[$comp], $position, 1) ne "t") and (substr ($align[$comp], $position, 1) ne "c"))
{
    $gapsites[$count]=$position;
    $count +=1;
$comp=scalar(@align);}


}}

return @gapsites;

}



#de_gap
#a subroutine to remove positions with gaps from an alignment
sub de_gap
{ my (@align) = @_;
  for (my $position=0; $position< length ($align[0]); ++$position)
  {  for ( my $seq=0; $seq<scalar (@align); ++$seq)
     { if ((substr($align[$seq], $position, 1) ne 'a') and (substr($align[$seq], $position, 1) ne 't') and (substr($align[$seq], $position, 1) ne 'g') and (substr($align[$seq], $position, 1) ne 'c'))
       { $seq = scalar (@align);
	 for (my $al=0; $al<scalar (@align); ++$al)
	 { substr ($align[$al], $position, 1) = "";
	  
       }
 $position=$position-1;   

 }
 }
 }

  return @align;

}












#count_SNPs
#a subroutine to check for SNPs and count the total- 
#excludes indels, gaps, and#segregating sites in indels 
#multiple substitutions per site are counted as 1 SNP
sub count_snps
{   my $gap=0;
   my @align =@_;
    my($length) = length ($align[0]);
    my($position);
    my($comp);
    my($count) = 0;
   
 for ($position=0; $position < $length ; ++$position)
 { $gap=0;
 
	for (my $seq=0; $seq<scalar (@align); ++$seq)
{
if ((substr ($align[$seq], $position, 1) ne 'a') and (substr ($align[$seq], $position, 1) ne 'g') and (substr ($align[$seq], $position, 1) ne 't') and (substr ($align[$seq], $position, 1) ne 'c'))
{$seq=scalar(@align);
 $gap=1;}
}
 if ($gap ne 1)
{
 for ($comp=1; $comp<(scalar(@align)); ++$comp)
      { 
if (substr($align[0], $position, 1) ne substr($align[$comp], $position, 1))
{++$count;
 $comp=scalar(@align);

}}
}
}

return $count;

}








#count_SNPs_missdata
#a subroutine to check for SNPs and count the total- 
#includes sites with missing data
sub count_snps_missdata
{   my $gap=0;
   my @align =@_;
    my($length) = length ($align[0]);
    my($position);
    my($comp);
    my($count) = 0;

 for ($position=0; $position < $length ; ++$position)
 { $gap=0;
 
	for (my $seq=0; $seq<scalar (@align); ++$seq)
{
if ((substr ($align[$seq], $position, 1) ne 'A') and (substr ($align[$seq], $position, 1) ne 'G') and (substr ($align[$seq], $position, 1) ne 'T') and (substr ($align[$seq], $position, 1) ne 'C'))
{
    $gap=1;
}

else {$compseq=$seq;
      $seq=scalar(@align);
      $gap=0; 
      
}
}



 if ($gap ne 1)
{
 for ($comp=0; $comp<(scalar(@align)); ++$comp)
      { 
if (substr($align[$compseq], $position, 1) ne substr($align[$comp], $position, 1))
{


++$countit;
 $comp=scalar(@align);


}}
}
}

return $countit;

}




#poly_sites
#a subroutine to return an array of positions of polymorphic sites
sub poly_sites
{  my $gap=0;
   my @align =@_;
    length ($align[0]);
    my($position);
    my($comp);
    my($count) = 0;
   
 for ($position=0; $position < length ($align[0]) ; ++$position)
 { $gap=0;
 
	for (my $seq=0; $seq<scalar (@align); ++$seq)
{
if ((substr ($align[$seq], $position, 1) ne 'A') and (substr ($align[$seq], $position, 1) ne 'G') and (substr ($align[$seq], $position, 1) ne 'T') and (substr ($align[$seq], $position, 1) ne 'C'))
{$seq=scalar(@align);
 $gap=1;}
}
 if ($gap ne 1)
{
 for ($comp=1; $comp<(scalar(@align)); ++$comp)
      { 
if (substr($align[0], $position, 1) ne substr($align[$comp], $position, 1))
{
    $polysites[$count]=$position;
$count+=1;
 $comp=scalar(@align);

}}
}
}

return @polysites;

}
#a subroutine to return an array of positions of polymorphic sites
sub poly_sites_aminoacid
{  my $gap=0;
   my @align =@_;
    length ($align[0]);
    my($position);
    my($comp);
    my($count) = 0;
   
 for ($position=0; $position < length ($align[0]) ; ++$position)
 { $gap=0;
 
	for (my $seq=0; $seq<scalar (@align); ++$seq)
{
if (substr ($align[$seq], $position, 1) eq '?')
{$seq=scalar(@align);
 $gap=1;}
}
 if ($gap ne 1)
{
 for ($comp=1; $comp<(scalar(@align)); ++$comp)
      { 
if (substr($align[0], $position, 1) ne substr($align[$comp], $position, 1))
{
    $polysites[$count]=$position;
$count+=1;
 $comp=scalar(@align);

}}
}
}

return @polysites;

}

#total_sites
#a subroutine to return the total number of sites containing standard nucleotides (no gaps, n's, etc)
sub total_sites
{
    @align =@_;
    my($length) = length ($align[0]);
    my($position);
    my($comp);
    my($count) = 0;
    my @polysites;
 for ($position=0; $position < $length ; ++$position)
    {

 for ($comp=0; $comp<(scalar(@align)); ++$comp)
      { 
if ((substr ($align[$comp], $position, 1) ne "a") and (substr ($align[$comp], $position, 1) ne "g") and (substr ($align[$comp], $position, 1) ne "t") and (substr ($align[$comp], $position, 1) ne "c"))
{
$comp=scalar(@align);}

else 
{if ($comp == (scalar(@align-1)))
{
$count +=1;
} 

}
}}

return $count;

}


#degap_position
#a subroutine to determine the position of a nucleotide site 
#once gaps have been removed from the sequence
sub degap_sequence
{    my ($seq, $pos) = @_;
my $gapnum = 0;
for (my $i=0; $i<$pos; ++$i)
 { if (substr($seq, $i, 1) eq "-")
         {++$gapnum;}  
}

my $newpos = $pos-$gapnum;


return $newpos;
}

#degap_polsites
#reasigns polymorphic site positions from a fasta alignment to their positions #without gaps
sub degap_polsites
{ my (@align, @polsites) = @_;
}




#calculate_pi
sub calculate_pi
{
    $diffs=0;
    $comps=0;
 my @align = @_;
  for (my $i=0; $i<scalar (@align); ++$i)
{ for ($t=($i+1); $t< scalar (@align); ++$t)
{$comps+=1;
  { for ($j=0; $j< length($align[0]); ++$j)
    { 
	
if (substr ($align[$i], $j, 1) ne substr ($align[$t], $j, 1))
      {$diffs+=1;}
  

}
}
}
}
if ($comps>0){    $pi= $diffs/$comps;}
else {$pi=0;}
    return $pi;
}

sub count_haplotypes
{

    my @align=@_;
   
my $hapnum=1;
    my @haptype;
    $haptype[0]= $align[0];
for (my $r=1; $r<scalar(@align); ++$r)
{   for (my $s=0; $s<$hapnum; ++$s)
{if ($align[$r] ne $haptype[$s])
{$test=$s+1;
if ($test==$hapnum)
{$haptype[$hapnum]=$align[$r];
 $hapnum+=1;
 $s=$hapnum;


}
}
else {$s=$hapnum;}

}


}

  




return $hapnum;
}

#count_Rm
#carries out the algorithm of Hudson and Kaplan (1985) to count the minimum number of recombination events

sub count_Rm
{my @align=@_;
 my @type=();
 my $pairnum=0;
 my @pair = ();
 for (my $w=0; $w<length($align[0]); ++$w)
{
 
    for (my $r=($w+1); $r<length($align[0]); ++$r)
{
    $type[0]=substr($align[0], $w, 1). substr($align[0],$r,1);
  my  $typenum=1;

 for (my $s=1; $s<scalar (@align); ++$s)
 { for (my $q=0; $q<$typenum; ++$q)
{ 


if ((substr($align[$s], $w, 1). substr($align[$s],$r,1) eq $type[$q]))
{$q=$typenum;

}

else { if ($q==($typenum-1))
{$type[$typenum]=substr($align[$s], $w, 1). substr($align[$s],$r,1);
 $typenum+=1;
 $q=$typenum;
 


}


if ($typenum==4)
{$pair[$pairnum] = "$w $r\n";
 print "$pair[$pairnum]";
$pairnum+=1;
 

$q=$typenum;
 $s = scalar(@align);
}



}


}
}




}
}

for ($z=0; $z<$pairnum; ++$z)

{for ($s=($z+1); $s<$pairnum; ++$s)
{
    
 my @arraypair1 = split(/\s/,$pair[$z]);
   my  @arraypair2 = split(/\s/,$pair[$s]);
    $M=$arraypair1[0];
    $I=$arraypair2[0];
    $N=$arraypair1[1];
  $J=$arraypair2[1];



if (($M<=$I) && ($J<=$N) && ($z ne $s))
{    splice (@pair,$z,1);
     
     $z=0;
     $s=0;
     print "reject: $M $N over $I $J ($z $s)\n";     
    
}


elsif (($M==$I) && ($J>$N))
{    splice (@pair,$s,1);
     print "Reject: $I $J ($z $s)\n";

     $z=0;
     $s=0;

    
    }

$pairnum = scalar (@pair);


}}










for ($z=0; $z<$pairnum; ++$z)
{ for ($s=($z+1); $s<$pairnum; ++$s)
{ my @arraypair1 = split(/\s/,$pair[$z]);
   my  @arraypair2 = split(/\s/,$pair[$s]);
    $I1=$arraypair1[0];
    $M=$arraypair2[0];
    $J1=$arraypair1[1];
  $N=$arraypair2[1];

if (($M>$I1) && ($M<$J1))
{    splice (@pair,$s,1);
     $s=$s-1;
     print "REJect: $M $N\n";
    
}


$pairnum=scalar (@pair);




}}

print "\n";
for ($try=0; $try<$pairnum; ++$try)
{ 

    print "$pair[$try]\n";
}
return $pairnum;

}



sub compute_D {

    my ($N) = $_[0];
    my ($Theta) = $_[1];
    my ($Pi)=$_[2];
    my ($S) =$_[3];
    my ($D);
      
    
   my ($a1) =0; #denom of theta
   for ($i=1; $i<$N; $i++) {
      $a1=$a1 + 1.0/$i;
   }

    my ($a2) =0; #denom of theta
    for ($i=1; $i<$N; $i++) {
       $a2=$a2 + 1.0/($i*$i);
    }
		
    $b1=($N+1.0)/(3.0*($N-1.0));
    $b2=(2.0*($N*$N+$N+3.0))/(9.0*$N*($N-1.0));
	 
    $c1=$b1-(1.0/$a1);
    $e1=$c1/$a1;
 
    $c2 = $b2-($N+2.0)/($a1*$N) + $a2/($a1*$a1);
    $e2 = $c2/(($a1*$a1)+$a2);

    unless ($S==0) {$D = ($Pi-$Theta)/(sqrt($e1*$S+$e2*$S*($S-1.0)));}
    if ($S==0) {$D=0;}
    
    return $D;
}



1;













