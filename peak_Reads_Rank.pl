#!/usr/bin/perl -w
## Here, the absolute count of reads is used (no need for normalization in a sample).
## This code is taken from "ratio_PeakFlank.pl".
## 该脚本对read数目进行按照从大到小排列。当相同reads数目情况下，排列对应和排列跳跃的算法;
## 结合array 和hash同时实现;
use strict;
use vars qw/%peaks @rank/;

if (@ARGV!=3)
{
	die " [Usage]  perl  $0 <peak.bed>  <input.wig>  <output.txt> \n";
	exit 0;
}

open PEAK, "$ARGV[0]";
while(<PEAK>)
{
	chomp;
	my ($chr, $start, $end) = split/\t/;
	$peaks{$chr}->{$start}->{End} = $end;
}
close PEAK;

my $count = {};
open WIG, "$ARGV[1]";
$/ = "variableStep";
<WIG>;
while(<WIG>)
{
	chomp;
	my ($chr) = ($_=~/chrom\=(\S+)/x);

	my @each = split/\n/;
	foreach my $num (1..$#each)
	{
		my ($site, $density) = split/\t/, $each[$num];
		$count->{$chr}->{$site} = $density;
	}
}
close WIG;

foreach my $chrm (sort keys %peaks)
{
	if ( $count->{$chrm} )
	{
		foreach my $coord (sort keys %{$peaks{$chrm}} )
		{
			if ( $count->{$chrm}->{$coord} )
			{
				my @peak = ();
				foreach my $site ($coord..$peaks{$chrm}->{$coord}->{End} )
				{
					if ($count->{$chrm}->{$site} )
					{
						push @peak, $count->{$chrm}->{$site};
					}
				}
				my $mean_peak = mean(\@peak);         ## average read count in peak region;
				push @rank, $mean_peak;
				$peaks{$chrm}->{$coord}->{Rank} = $mean_peak;
			}
		}
	}
}

my %rank_read = ();
my $k = 1;
my @ranked = sort { $b <=> $a} @rank;        ## descending order;
my $init = { $ranked[0] => 1 };
$rank_read{$ranked[0]} = $k;        ## 定义起始状态;

 foreach my $m (1..$#ranked)
{
	if ( $init->{$ranked[$m]} )        ## 每次与上一次比较，判断是否存在;
	{
		$rank_read{$ranked[$m]}= $k;
		$init->{$ranked[$m]} = 1;
	}

	else
	{
		$k = $m + 1;          ## 对于不存在的情况，保持rank一致性;
		$rank_read{$ranked[$m]}= $m+1;
		$init->{$ranked[$m]} = 1;
	}
}
open OUT, ">$ARGV[2]";
foreach my $chr (sort keys %peaks)
{
	foreach my $begin (sort keys %{$peaks{$chr}} )
	{
		print OUT $chr, "\t", $begin, "\t", $peaks{$chr}->{$begin}->{End}, "\t";
		print OUT $peaks{$chr}->{$begin}->{Rank}, "\t";
		print OUT $rank_read{$peaks{$chr}->{$begin}->{Rank}}, "\n";
	}
}
close OUT;
exit 0;


#####################
sub mean
{
	my $value = shift;
	my $sum = 0;
	my $size = 0;
	my $mean ;
	foreach (0..$#$value)
	{
		$sum += $value->[$_];
		$size++;
	}
	if ($size == 0)    ## null array;
	{
		$mean = 1;         ## 人为定义，因为区域内如果没有reads，那么就将其当作背景1来处理;
	}
	else
	{
		$mean = sprintf ("%.3f", $sum/$size);
	}

	return $mean;
}
__END__
