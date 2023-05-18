#!/usr/bin/env perl
#########################################################################
# File Name: random_selection.pl
# Author: Rui
# mail: dongruipicb@gmail.com
# Created Time: Mon 19 Jan 2015 01:45:45 PM CST
#########################################################################

use strict;
use warnings;
use List::Util qw(shuffle);
use Getopt::Long;

sub usage {
	print <<"END_USAGE";
Usage: perl $0
	--file
	--num
END_USAGE
	exit;
}


my ($file,$num);
GetOptions (
	'file=s'=>\$file,
	'num=s'=>\$num,
) or usage();
usage if (!$file or !$num);

open my $file_in,"$file";
my @all_info;
while (<$file_in>){
	chomp;
	push @all_info,$_;
}
my @shuffle_all_info=shuffle(@all_info);
for my $num_count(0..$num-1){
	print "$shuffle_all_info[$num_count]\n";
}
