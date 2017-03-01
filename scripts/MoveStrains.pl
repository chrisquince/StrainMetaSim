#!/usr/bin/perl

use strict;
use warnings;

my $speciesSelect = $ARGV[0];

mkdir "Strain_${speciesSelect}";

while(my $line = <STDIN>){
    chomp($line);

    my @tokens = split(/\t/,$line);

    my $species = $tokens[6];

    if($species eq $speciesSelect){
        system("cp -r $tokens[0] Strain_${speciesSelect}/$tokens[0]");
    }

}
