


my $count;
my $sum;
my $first;
my $last;

open(BLOC, "$ARGV[0]") or die $!;
while(<BLOC>) {
	        my $line = <BLOC>;
	        $line = sprintf("%.3f", $line);
	        $sum = $sum + $line;
	        print $line,"\n";
	        print $sum,"\n";
	        $count = $count + 1;
	    }
close(BLOC);
print $sum,"\n";
print $count,"\n";

my $mean = $sum / $count;
my $lastCount = $count;
my $total_dev;
my $count;

open(BLOC, "$ARGV[0]") or die $!;
while(<BLOC>) {
	        my $line = <BLOC>;
	        $line = sprintf("%.3f", $line);
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


