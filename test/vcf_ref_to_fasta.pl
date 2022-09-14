use strict;
use warnings;
my $next_pos=1;
while(<>){
    my @fields=split;
    my $line_pos=$fields[0];
    for(;$next_pos<$line_pos;$next_pos++){
        print "N" ;
    }
    print $fields[1];
    $next_pos++;
}
print "\n";