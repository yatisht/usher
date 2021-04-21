use warnings;
use strict;

my $POS_IDX=1;
my $REF_IDX=3;
my $ALT_IDX=4;
my $SAMPLE_START_IDX=9;

open(my $file1,"<",$ARGV[0]);
open(my $file2,"<",$ARGV[1]);

##parse header
my $headerLine=<$file1>;
$headerLine=<$file2>;

my %sample2;
$headerLine=<$file2>;
my @headerFields2=split("\t",$headerLine);
my $index=0;
foreach(@headerFields2){
    $sample2{$_}=$index;
    $index++;
}

my @sampleMap;
$headerLine=<$file1>;
my @headerFields1=split("\t",$headerLine);
foreach(@headerFields1){
    push @sampleMap, [$_,$sample2{$_}];
}


while(<$file1>){
    my @file1_fields=split("\t",$_);
    my @file2_fields=split("\t",<$file2>);
    if(!($file1_fields[$POS_IDX] eq $file2_fields[$POS_IDX])){
        die("Different Position, file1: $file1_fields[$POS_IDX] file2:$file2_fields[$POS_IDX]\n");
    }
    my @file1Nuc=($file1_fields[$REF_IDX]);
    push @file1Nuc, split(",",$file1_fields[$ALT_IDX]);

    my @file2Nuc=($file2_fields[$REF_IDX]);
    push @file2Nuc, split(",",$file2_fields[$ALT_IDX]);

    for(my $i=$SAMPLE_START_IDX;$i<(scalar @file1_fields);$i++){
        my $sample2Idx=$sampleMap[$i]->[1];
        my $sample1Nuc=$file1Nuc[$file1_fields[$i]];
        my $sample2Nuc=$file2Nuc[$file2_fields[$sample2Idx]];
        if($sample1Nuc ne $sample2Nuc){
            print ("At $file1_fields[$POS_IDX] of $sampleMap[$i]->[0], $ARGV[0] is $sample1Nuc, but $ARGV[1] is $sample2Nuc \n");
        }
    }
}