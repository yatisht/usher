use warnings;
use strict;

my $POS_IDX=1;
my $REF_IDX=3;
my $ALT_IDX=4;
my $SAMPLE_START_IDX=9;

open(my $file1,"<",$ARGV[0]);
open(my $file2,"<",$ARGV[1]);

##parse header
my $headerLine1="";
until($headerLine1=~/CHROM/){
    $headerLine1=<$file1>;
}
my $headerLine2="";
until ($headerLine2=~/CHROM/){
    $headerLine2=<$file2>;
}

my %sample2;
my @headerFields2=split("\t",$headerLine2);
my $index=0;
foreach(@headerFields2){
    chomp $_;
    $sample2{$_}=$index;
    $index++;
}

my @sampleMap;
my @headerFields1=split("\t",$headerLine1);
foreach(@headerFields1){
    chomp $_;
    if(not defined $sample2{$_}){
        print("$_ not found\n");
    }
    push @sampleMap, [$_,$sample2{$_}];
}

while(<$file1>){
    my @file1_fields=split(/\s/,$_);
    my @file2_fields=split(/\s/,<$file2>);
    if(!($file1_fields[$POS_IDX] eq $file2_fields[$POS_IDX])){
        die("Different Position, file1: $file1_fields[$POS_IDX] file2:$file2_fields[$POS_IDX]\n");
    }
    my @file1Nuc=($file1_fields[$REF_IDX]);
    push @file1Nuc, split(",",$file1_fields[$ALT_IDX]);

    my @file2Nuc=($file2_fields[$REF_IDX]);
    push @file2Nuc, split(",",$file2_fields[$ALT_IDX]);

    for(my $i=$SAMPLE_START_IDX;$i<(scalar @file1_fields);$i++){
        my $sample2Idx=$sampleMap[$i]->[1];
        if(not defined $sample2Idx){
            print("a");
        }
        $file1_fields[$i]=~/(\d+)/;
        my $sample1Nuc=$file1Nuc[$1];
        $file2_fields[$sample2Idx]=~/(\d+)/;
        my $sample2Nuc=$file2Nuc[$1];
        if($sample1Nuc ne $sample2Nuc){
            print ("At $file1_fields[$POS_IDX] of $sampleMap[$i]->[0], $ARGV[0] is $sample1Nuc, but $ARGV[1] is $sample2Nuc \n");
        }
    }
}
unless(eof($file2)){
    my @file2_fields=split("\t",<$file2>);
    print ("File 2 have extra lines, position starts at $file2_fields[1] \n");
}