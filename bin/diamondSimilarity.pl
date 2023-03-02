#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(first_index);
my ($fasta,$result,$output, $minLen, $minPercent, $minPval);
&GetOptions("fasta=s"=> \$fasta,
            "result=s"=> \$result,
            "output=s"=> \$output,
            "minLen=i"=> \$minLen,
            "minPercent=i"=> \$minPercent,
            "minPval=f"=> \$minPval);
$minPval = $minPval ? $minPval : 1e-5;
$minLen = $minLen ? $minLen : 10;
$minPercent = $minPercent ? $minPercent : 20; 
my %subjectCountHash = &getSubjectCount($fasta,$result,$minLen,$minPercent,$minPval);
my @deflines = &getSubjectIds($fasta);
open(my $data, '<', $result) || die "Could not open file $result: $!";
open(OUT,">$output");
open(LOG,">diamondSimilarity.log");
my $previousQseqId;
my $previousSseqId = "hello";
my $counter = 0;
while (my $line = <$data>) {
    chomp $line;
    my ($qseqid,$qlen,$sseqid,$slen,$qstart,$qend,$sstart,$send,$evalue,$bitscore,$length,$nident,$pident,$positive,$qframe,$qstrand) = split(/\t/, $line);
    my $subjectCount = $subjectCountHash{">$qseqid"};
    if ($counter == 0) {
	&addSubjectsNoCount($previousQseqId,$qseqid,1,@deflines);
	$subjectCount == 0 ? print LOG "No HSPS passed requirements for $qseqid\n" : print OUT ">$qseqid ($subjectCount subjects)\n";
	$previousQseqId = $qseqid;
    }
    if ($qseqid eq $previousQseqId) {
	if ($previousSseqId eq $sseqid) {
	    print LOG "Multiple HSPS for $qseqid, $sseqid\n";
	    die;
	}
	if ($length >= $minLen && $pident >= $minPercent && $evalue <= $minPval) {
            print OUT "  Sum: $sseqid:$bitscore:$evalue:$sstart:$send:$qstart:$qend:1:$length:$nident:$positive:$qstrand:$qframe:$qlen:$slen:100\n";
            print OUT "   HSP1: $sseqid:$nident:$positive:$length:$bitscore:$evalue:$sstart:$send:$qstart:$qend:$qstrand:$qframe\n";
	}
    } else {
	&addSubjectsNoCount($previousQseqId,$qseqid,0,@deflines);
	$subjectCount == 0 ? print LOG "No HSPS passed requirements for $qseqid\n" : print OUT ">$qseqid ($subjectCount subjects)\n";
	if ($length >= $minLen && $pident >= $minPercent && $evalue <= $minPval) {
          print OUT "  Sum: $sseqid:$bitscore:$evalue:$sstart:$send:$qstart:$qend:1:$length:$nident:$positive:$qstrand:$qframe:$qlen:$slen:100\n";
          print OUT "   HSP1: $sseqid:$nident:$positive:$length:$bitscore:$evalue:$sstart:$send:$qstart:$qend:$qstrand:$qframe\n";
	}
    }
    $previousQseqId = $qseqid;
    $previousSseqId = $sseqid; 
    $counter += 1;
}
close($data);
close OUT;
close LOG;
sub getSubjectIds{
    my ($fasta) = @_;
    my @deflines = `grep "^>" $fasta | awk 'sub(/\\slength=[0-9]+/, "")'`;
    chomp(@deflines);
    return @deflines;
}
sub getSubjectCount{
    my ($fasta,$outFile,$minLen,$minPercent,$minPval) = @_;
    my @deflines = &getSubjectIds($fasta);
    my %hash = map { $_ => 0 } @deflines;
    open(my $data, '<', $outFile) || die "Could not open file $outFile: $!";
    while (my $line = <$data>) {
        chomp $line;
        my ($qseqid,$qlen,$sseqid,$slen,$qstart,$qend,$sstart,$send,$evalue,$bitscore,$length,$nident,$pident,$positive,$qframe,$qstrand) = split(/\t/, $line);
	if ($length >= $minLen && $pident >= $minPercent && $evalue <= $minPval) {
	    $hash{">$qseqid"} += 1;
	}
    }
    close($data);
    return %hash;
}
sub addSubjectsNoCount {
    my ($prevQSeqId,$qSeqId,$isFirstQSeq,@subjectIds) = @_;
    my $counter = 0;
    if ($isFirstQSeq) {
        until ($subjectIds[$counter] eq ">$qSeqId") {
            print OUT "$subjectIds[$counter] (0 subjects)\n";
	    $counter += 1;
	}
    }
    else {
        my $firstIndex = first_index { $_ eq ">$prevQSeqId" } @subjectIds;
	$firstIndex += 1;
	my $currentIndex = first_index { $_ eq ">$qSeqId" } @subjectIds;
        until ($firstIndex == $currentIndex) {
	    print OUT "$subjectIds[$firstIndex] (0 subjects)\n";
	    $firstIndex += 1;
        }
    }
}
