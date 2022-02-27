#!/bin/bash

echo "grouping read pairs to unique consensus source molecules"

slurp -s 100M $NAME_BAM_FILE |
samtools view -F 4 - | # strip SAM header; suppress unmapped reads (orphaned reads persist)

perl $ACTION_DIR/group/parse_read_pairs.pl | # extract signatures of all unique input source molecules


sort --parallel=$N_CPU -T $TMP_DIR_WRK -S $MAX_SORT_RAM --compress-program=pigz -k1,1 -k2,2n | # sort by mol_key/pos1

perl $ACTION_DIR/group/group_reads.pl | # apply fuzzy position logic to establish molecule groups
perl $ACTION_DIR/consensus/make_consensus.pl # merge reads and build strand/molecule consensuses
checkPipe

echo "done"
echo
