use Bio::SeqIO;
use Getopt::Long qw(:config no_ignore_case);

my %opts;
GetOptions(\%opts, "i:s");
my $usage = <<USAGE;

    Program: Sort
    Contact: Felix Beaudry(felix.beaudry\@utoronto.ca)

    Usage:	$0 -i intake_directory

    			-i input


 	Example $0 -i test_files


USAGE

die $usage unless ($opts{i});


open my $file, '<', $opts{i} or die "Could not open $file: $!";


$in  = Bio::SeqIO->new(-file => "$opts{i}",
                       -format => 'Fasta');
$out = Bio::SeqIO->new(-file => ">$opts{i}.sort",
                       -format => 'Fasta');
while ( my $seq = $in->next_seq() ) {$out->write_seq($seq); }