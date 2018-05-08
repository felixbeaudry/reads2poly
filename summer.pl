use strict;

my $sum = 0;
my $count = 0;

open(BLOC, "$ARGV[0]") or die $!;
while(<BLOC>) {
	        my $line = <BLOC>;
	        $sum = $sum + $line;
	        $count = $count + 1;

	    }
close(BLOC);
my $mean = $sum / $count;
print $mean, "/n";


