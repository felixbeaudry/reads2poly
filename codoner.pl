
use strict;
my $file =  "$ARGV[0]";
open my $info, $file or die "Could not open $file: $!";
my $dash_count = 0;
while(my $line = <$info>){
	#print "start\n";
	chomp $line;
	#print $line;
	
	if($line =~ />/) {
		print $line, "\n";
	} 
	elsif(length($line) == 60){
		my $count = () = $line =~ /-/g;
		$dash_count =  $count + $dash_count;
		print $line, "\n";
	}
	else{
		print $line;
		my $extra_ncl = (length($line) - $dash_count)  % 3;
		my $x = 0;
		if ($extra_ncl > 0){
		for ($x=0; $x<(3-$extra_ncl);++$x){ 
			print "N";
		}}
		#print " ", length($line), " ", $extra_ncl, " ", $dash_count;	
		print "\n";
	}
}	
