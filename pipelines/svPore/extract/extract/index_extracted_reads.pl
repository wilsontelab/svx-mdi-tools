use strict;
use warnings;

# create an index for rapid retrieval of SV-containing reads

# constants
use constant {
    molType => 0, # extract nodes fields
    qName => 1,
    qSeq => 2,
    qQual => 3
};

my $offset = 0;
while(my $line = <STDIN>){
    my $nChar = length($line);    
    my @line = split("\t", $line);
    print join("\t", $line[qName], $offset, $nChar), "\n";
    $offset += $nChar;
}
