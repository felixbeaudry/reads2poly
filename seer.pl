


use strict;


open(mean_file, "$ARGV[1]") or die $!;
my $mean = <mean_file>
close(mean_file);
print $mean


my $total_dev = 0;
my $count = 0;

open(BLOC, "$ARGV[0]") or die $!;
while(<BLOC>) {
	        my $line = <BLOC>;
	        my $dev = ( $line - $mean )**2;
	        $total_dev = $total_dev + $dev;
	        $count = $count + 1;

	    }
close(BLOC);
my $sd =  sqrt($total_dev / ($count - 1));
my $se = $sd / sqrt($count);


print $se, "\n";



