#! /bin/env perl

use strict;

my %delete_sites;
my %select_sites;
open (VCF, "allsites.vcf") or die $!;
my $line = <VCF>;

while (<VCF>) {
    next if (/^#/);
    next if (/\.\/\./);
    next if (/lowQUAL\t/);
    chomp;
    $line = $_;
    my ($scf, $pis, $id, $ref, $alt, $qual,$filter, $info, $format, @gts) = split(/\t/, $line);
    next if ($qual < 60);
    next if ($alt eq '.');
    my @bases = split(/,/, $alt);
    unshift @bases, $ref;
    my $indel_flag = 0;
    foreach my $tmpstr (@bases) {
        if (length($tmpstr)>1) {
            $indel_flag++;
        }
    }
    if ($indel_flag>0) {
        my $insert_length = 0;
        if (length($bases[0]) > 1) {
            $insert_length = length($bases[0]) -1;
        }
        for ($pis-10 .. $pis+10+$insert_length) {
            my $identi = $scf . "__" . $_;
            $delete_sites{$identi} = 0;
        }
        next;
    }
    my $ident = $scf . "__" . $pis;

    my @gts = map { s/:.+//; $_ } @gts;
	if ( $gts[0] eq $gts[4] && $gts[0] eq $gts[5] && $gts[0] eq $gts[6] && $gts[0] eq $gts[7] && $gts[0] eq $gts[8] && $gts[0] eq $gts[12] && $gts[0] eq $gts[13] && $gts[0] eq $gts[14] && $gts[0] eq $gts[15] && $gts[0] eq $gts[16] && $gts[0] eq $gts[17] && $gts[1] eq $gts[2] && $gts[1] eq $gts[3] && $gts[1] eq $gts[9] && $gts[1] eq $gts[10] && $gts[1] eq $gts[11] && $gts[1] eq $gts[18] && $gts[1] eq $gts[19] && $gts[1] eq $gts[20] && $gts[1] eq $gts[21] && $gts[1] eq $gts[22] && $gts[1] eq $gts[23] && $gts[1] eq $gts[24] && $gts[1] eq $gts[25] && $gts[1] eq $gts[27] && $gts[1] eq $gts[28] && $gts[1] eq $gts[29] && $gts[1] eq $gts[30] && $gts[1] eq $gts[35] && $gts[1] eq $gts[36] && $gts[1] eq $gts[37] && $gts[1] eq $gts[38] && $gts[1] eq $gts[40] && $gts[1] eq $gts[41] && $gts[1] eq $gts[42] && $gts[1] eq $gts[43] && $gts[1] eq $gts[44] && $gts[1] eq $gts[45] && $gts[1] eq $gts[46] && $gts[0] ne $gts[1]
	){
	$select_sites{$ident} = $line;
    }
}

foreach (keys %delete_sites) {
    delete($select_sites{$_});
}

foreach (sort keys %select_sites) {
    print $select_sites{$_},"\n";
}
