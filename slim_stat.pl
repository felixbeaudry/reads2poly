

use strict;

my $count = 0;
my $sum = 0;
my $first;
my $last;

open(BLOC, "$ARGV[0]") or die $!;
while(<BLOC>) {
	        my $line = <BLOC>;
	        $sum = $sum + $line;
	        $count = $count + 1;
	    }
close(BLOC);
my $mean = $sum / $count;
my $lastCount = $count;
my $total_dev = 0;
my $count = 0;

open(BLOC, "$ARGV[0]") or die $!;
while(<BLOC>) {
	        my $line = <BLOC>;
	       	if ($count == 0){
	        	$first = $line;
	        }
	        elsif($count == $lastCount){
	        	$last = $line;
	        }
	        else{}
	        my $dev = ( $line - $mean )**2;
	        $total_dev = $total_dev + $dev;
	        $count = $count + 1;

	    }
close(BLOC);
my $sd =  sqrt($total_dev / ($count - 1));
my $se = $sd / sqrt($count);
my $slope = ($last - $first) / $count;

print $mean, " ", $se, " ", $slope, "\n";

######


