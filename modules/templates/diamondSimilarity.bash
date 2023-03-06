#!/usr/bin/env bash

set -euo pipefail
 
diamond blastp \
	-d $database \
	-q $fasta \
	-o out.txt \
	-f 6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length nident pident positive qframe qstrand gaps \
	--comp-based-stats 0 \
	$blastArgs

perl /usr/bin/diamondSimilarity.pl \
     --fasta $fasta \
     --result out.txt \
     --output diamondSimilarity.out \
     --minLen $lengthCutoff \
     --minPercent $percentCutoff \
     --minPval $pValCutoff
