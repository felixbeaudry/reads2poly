#!perl 
##Perl script for extracting Hudson data (from ms), edited from S. Wright


use hudsonsubsrho_fb;
my $dirfile=$ARGV[0];
my $popSplit=int($ARGV[1]);

@file_data = get_file_data ($dirfile);
@data = extract_sequence_from_hudson_data (@file_data);

#$reps=$ARGV[1];

my $i=0;
my @sequence;
my @sequence_one;
my @sequence_two;
my $pi_tot = 0;
my $pi_one = 0;
my $pi_two = 0;
my $fst = 0;
my $rep = 0;
print "rep\tpi_tot\tfst\n"; 
foreach my $line (@file_data) {
	if ($line =~ /ms/){next;}
	elsif ($line =~ /^\s*$/) {next;}
	elsif ($line =~ /\/\//){
		if ($rep == 0){
			$rep +=1;
		}
		else{
			$rep +=1;
			$pi_tot = calculate_pi(@sequence);
			$pi_one = calculate_pi(@sequence_one);
			$pi_two = calculate_pi(@sequence_two);
			$fst = ( $pi_tot - ( ($pi_one + $pi_two) / 2 )) / $pi_tot;
			if ($fst < 0 ){$fst=0;}
			print $rep-1, "\t", $pi_tot / 1000, "\t", $fst, "\n";


			#reset
			my $pi_tot = 0;
			my $pi_one = 0;
			my $pi_two = 0;
			my $fst = 0;
	    	$i=0;
	    	my @sequence;
	    	my @sequence_one;
			my @sequence_two;
	    }
  	}
	elsif ($line =~ /segsites/){next;}
	#elsif ($line =~ /^rho/){next;}
	elsif ($line =~ /positions/) {next;}
	else { 
	  $sequence[$i] = $line; 
	  
	  if ($i < $popSplit){
	  	$sequence_one[$i] = $line;

	  }  
	  else{
	  	$sequence_two[$i - $popSplit] = $line;

	  }
	  $i+=1;
	}
}  







