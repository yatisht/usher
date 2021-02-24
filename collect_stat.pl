use warnings;
use strict;
sub get_stat{
	my ($path)=@_;
	my ($radius, $test_case)=$path=~/.*r(\d+)\/(.*)\/log/;
	open(my $FH,"<",$path);
	my $before;
	my $after;
	my $n_nodes;
	my $realtime;
	my $usertime;
	while(<$FH>){
		if($_=~/Before refinement: (\d+)/ and (not defined $before)){
			$before=$1;
		}
		if($_=~/After refinement: (\d+)/){
			$after=$1;
		}
		if($_=~/real	(.*)/){
			$realtime=$1;
		}
		if($_=~/user	(.*)/){
			$usertime=$1;
		}
		if($_=~/Current tree size \(#nodes\): (\d*)/){
			$n_nodes=$1;
		}
	}
	if(defined $realtime){
		print("$radius,$test_case,$n_nodes,$before,$after,$realtime,$usertime\n");
	}
}

my @logs=glob("testout/r*/*/log");
print("Radius,Test Case,Numer of Nodes,Before,After,Real Time,CPU Time\n");
foreach(@logs){
	get_stat($_);
}

