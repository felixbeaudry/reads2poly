##Perl subroutines for Hudson data (from ms), edited from S. Wright


#calculate_pi
sub calculate_pi{
  $diffs=0;
  $comps=0;
  my @align = @_;
  if ($align[0]){
    for (my $i=0; $i<scalar(@align); ++$i){
      for ($t=($i+1); $t< scalar (@align); ++$t){
        $comps+=1;
        for ($j=0; $j< length($align[0]); ++$j){ 
	         if (substr ($align[$i], $j, 1) ne substr ($align[$t], $j, 1)){$diffs+=1;}
    } } }
    if ($comps>0){$pi= $diffs/$comps;}
    else {$pi=0;}
  }
  else {$pi=0;}
  return $pi;
}

#extract_sequence_from_hudson_data
# A subroutine to extract hudson sequence data from an array
sub extract_sequence_from_hudson_data {
  my(@hudson_file_data) = @_;
  use strict;
  my $i=0;
  my @sequence;
  foreach my $line (@hudson_file_data) {
    #discard blank lines
    if ($line =~ /^\s*$/) {next;}
    #discard header lines
    elsif ($line =~ /^segsites/){next;}
    elsif ($line =~ /^rho/){next;}
    elsif ($line =~ /^positions/) {next;}
    elsif ($line =~ /^\/\//){next;}
    elsif ($line =~ /ms/){next;}
    else { 
      $sequence[$i] .= $line; 
      $i+=1;
    }
  }
  return @sequence;
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

    
}


elsif (($M==$I) && ($J>$N))
{    splice (@pair,$s,1);


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

    
}


$pairnum=scalar (@pair);




}}



return $pairnum;

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





#get_rho
# a subroutine to extract the value of rho from hudsonsubs
sub get_rho
{   my(@hudson_file_data) = @_;
    
    my $i=0;
    my @sequence;
    foreach my $line (@hudson_file_data) {
if ($line =~ /^rho/)
{my @dat=split (' ', $line);
 my $rho=$dat[1]; 
 print "$rho\n";

}


}

return $rho;
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

#count_SNPs
#a subroutine to check for SNPs and count the total- 
#excludes indels, gaps, and#segregating sites in indels 
#multiple substitutions per site are counted as 1 SNP
sub count_snps
{
    @align =@_;
    my($length) = length ($align[0]);
    my($position);
    my($comp);
    my($count) = 0;
   
 for ($position=0; $position < $length ; ++$position)
    {
 if (substr ($align[0], $position, 1) ne "-")
{
 for ($comp=1; $comp<(scalar(@align)); ++$comp)
      { 
if (substr ($align[$comp], $position, 1) eq "-")
{
$comp=scalar(@align);}

elsif (substr($align[0], $position, 1) ne substr($align[$comp], $position, 1))
{++$count;
 $comp=scalar(@align);

}}
}
}



return $count;

}

#poly_sites
#a subroutine to return an array of positions of polymorphic sites
sub poly_sites
{
    @align =@_;
    my($length) = length ($align[0]);
    my($position);
    my($comp);
    my($count) = 0;
    my @polysites;
 for ($position=0; $position < $length ; ++$position)
    {
 if (substr ($align[0], $position, 1) ne "-")
{
 for ($comp=1; $comp<(scalar(@align)); ++$comp)
      { 
if (substr ($align[$comp], $position, 1) eq "-")
{
$comp=scalar(@align);}

elsif (substr($align[0], $position, 1) ne substr($align[$comp], $position, 1))
{$polysites[$count]=($position);
++$count;
 $comp = scalar(@align);

}}
}
}

return @polysites;

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







1;













