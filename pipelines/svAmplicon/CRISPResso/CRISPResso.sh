# do the work
CRISPResso \
--fastq_r1 ${FASTQ_R1} \
--fastq_r2 ${FASTQ_R2} \
--amplicon_seq ${AMPLICON_SEQ} \
--output_folder ${TASK_DIR} \
--file_prefix ${DATA_NAME} \
--name ${DATA_NAME} \
--n_processes ${N_CPU}
checkPipe


# options:
#   -h, --help            show this help message and exit
#   --version             show program's version number and exit
#   -r1 FASTQ_R1, --fastq_r1 FASTQ_R1
#                         First fastq file (default: )
#   -r2 FASTQ_R2, --fastq_r2 FASTQ_R2
#                         Second fastq file for paired end reads (default: )
#   -a AMPLICON_SEQ, --amplicon_seq AMPLICON_SEQ
#                         Amplicon Sequence (can be comma-separated list of multiple sequences) (default: None)
#   -an AMPLICON_NAME, --amplicon_name AMPLICON_NAME
#                         Amplicon Name (can be comma-separated list of multiple names, corresponding to amplicon sequences given in --amplicon_seq (default:
#                         Reference)
#   -amas AMPLICON_MIN_ALIGNMENT_SCORE, --amplicon_min_alignment_score AMPLICON_MIN_ALIGNMENT_SCORE
#                         Amplicon Minimum Alignment Score; score between 0 and 100; sequences must have at least this homology score with the amplicon to be
#                         aligned (can be comma-separated list of multiple scores, corresponding to amplicon sequences given in --amplicon_seq) (default: )
#   --default_min_aln_score DEFAULT_MIN_ALN_SCORE, --min_identity_score DEFAULT_MIN_ALN_SCORE
#                         Default minimum homology score for a read to align to a reference amplicon (default: 60)
#   --expand_ambiguous_alignments
#                         If more than one reference amplicon is given, reads that align to multiple reference amplicons will count equally toward each amplicon.
#                         Default behavior is to exclude ambiguous alignments. (default: False)
#   --assign_ambiguous_alignments_to_first_reference
#                         If more than one reference amplicon is given, ambiguous reads that align with the same score to multiple amplicons will be assigned to the
#                         first amplicon. Default behavior is to exclude ambiguous alignments. (default: False)
#   -g GUIDE_SEQ, --guide_seq GUIDE_SEQ, --sgRNA GUIDE_SEQ
#                         sgRNA sequence, if more than one, please separate by commas. Note that the sgRNA needs to be input as the guide RNA sequence (usually 20
#                         nt) immediately adjacent to but not including the PAM sequence (5' of NGG for SpCas9). If the PAM is found on the opposite strand with
#                         respect to the Amplicon Sequence, ensure the sgRNA sequence is also found on the opposite strand. The CRISPResso convention is to depict
#                         the expected cleavage position using the value of the parameter '--quantification_window_center' nucleotides from the 3' end of the guide.
#                         In addition, the use of alternate nucleases besides SpCas9 is supported. For example, if using the Cpf1 system, enter the sequence
#                         (usually 20 nt) immediately 3' of the PAM sequence and explicitly set the '--cleavage_offset' parameter to 1, since the default setting of
#                         -3 is suitable only for SpCas9. (default: )
#   -gn GUIDE_NAME, --guide_name GUIDE_NAME
#                         sgRNA names, if more than one, please separate by commas. (default: )
#   -fg FLEXIGUIDE_SEQ, --flexiguide_seq FLEXIGUIDE_SEQ
#                         sgRNA sequence (flexible) (can be comma-separated list of multiple flexiguides). The flexiguide sequence will be aligned to the amplicon
#                         sequence(s), as long as the guide sequence has homology as set by --flexiguide_homology. (default: None)
#   -fh FLEXIGUIDE_HOMOLOGY, --flexiguide_homology FLEXIGUIDE_HOMOLOGY
#                         flexiguides will yield guides in amplicons with at least this homology to the flexiguide sequence. (default: 80)
#   -fgn FLEXIGUIDE_NAME, --flexiguide_name FLEXIGUIDE_NAME
#                         flexiguide name (default: )
#   --discard_guide_positions_overhanging_amplicon_edge
#                         If set, for guides that align to multiple positions, guide positions will be discarded if plotting around those regions would included bp
#                         that extend beyond the end of the amplicon. (default: False)
#   -e EXPECTED_HDR_AMPLICON_SEQ, --expected_hdr_amplicon_seq EXPECTED_HDR_AMPLICON_SEQ
#                         Amplicon sequence expected after HDR (default: )
#   -c CODING_SEQ, --coding_seq CODING_SEQ
#                         Subsequence/s of the amplicon sequence covering one or more coding sequences for frameshift analysis. If more than one (for example, split
#                         by intron/s), please separate by commas. (default: )
#   -q MIN_AVERAGE_READ_QUALITY, --min_average_read_quality MIN_AVERAGE_READ_QUALITY
#                         Minimum average quality score (phred33) to keep a read (default: 0)
#   -s MIN_SINGLE_BP_QUALITY, --min_single_bp_quality MIN_SINGLE_BP_QUALITY
#                         Minimum single bp score (phred33) to keep a read (default: 0)
#   --min_bp_quality_or_N MIN_BP_QUALITY_OR_N
#                         Bases with a quality score (phred33) less than this value will be set to "N" (default: 0)
#   --file_prefix FILE_PREFIX
#                         File prefix for output plots and tables (default: )
#   -n NAME, --name NAME  Output name of the report (default: the name is obtained from the filename of the fastq file/s used in input) (default: )
#   -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
#                         Output folder to use for the analysis (default: current folder) (default: )
#   -v VERBOSITY, --verbosity VERBOSITY
#                         Verbosity level of output to the console (1-4), 4 is the most verbose (default: 3)
#   --split_interleaved_input, --split_paired_end
#                         Splits a single fastq file containing paired end reads into two files before running CRISPResso (default: False)
#   --trim_sequences      Enable the trimming of Illumina adapters with Trimmomatic (default: False)
#   --trimmomatic_command TRIMMOMATIC_COMMAND
#                         Command to run trimmomatic (default: trimmomatic)
#   --trimmomatic_options_string TRIMMOMATIC_OPTIONS_STRING
#                         Override options for Trimmomatic, e.g. "ILLUMINACLIP:/data/NexteraPE-PE.fa:0:90:10:0:true" (default: )
#   --flash_command FLASH_COMMAND
#                         Command to run flash (default: flash)
#   --min_paired_end_reads_overlap MIN_PAIRED_END_READS_OVERLAP
#                         Parameter for the FLASH read merging step. Minimum required overlap length between two reads to provide a confident overlap. (default: 10)
#   --max_paired_end_reads_overlap MAX_PAIRED_END_READS_OVERLAP
#                         Parameter for the FLASH merging step. Maximum overlap length expected in approximately 90% of read pairs. Please see the FLASH manual for
#                         more information. (default: 100)
#   --stringent_flash_merging
#                         Use stringent parameters for flash merging. In the case where flash could merge R1 and R2 reads ambiguously, the expected overlap is
#                         calculated as 2*average_read_length - amplicon_length. The flash parameters for --min-overlap and --max-overlap will be set to prefer
#                         merged reads with length within 10bp of the expected overlap. These values override the --min_paired_end_reads_overlap or
#                         --max_paired_end_reads_overlap CRISPResso parameters. (default: False)
#   -w QUANTIFICATION_WINDOW_SIZE, --quantification_window_size QUANTIFICATION_WINDOW_SIZE, --window_around_sgrna QUANTIFICATION_WINDOW_SIZE
#                         Defines the size (in bp) of the quantification window extending from the position specified by the "--cleavage_offset" or "--
#                         quantification_window_center" parameter in relation to the provided guide RNA sequence(s) (--sgRNA). Mutations within this number of bp
#                         from the quantification window center are used in classifying reads as modified or unmodified. A value of 0 disables this window and
#                         indels in the entire amplicon are considered. Default is 1, 1bp on each side of the cleavage position for a total length of 2bp. Multiple
#                         quantification window sizes (corresponding to each guide specified by --guide_seq) can be specified with a comma-separated list. (default:
#                         1)
#   -wc QUANTIFICATION_WINDOW_CENTER, --quantification_window_center QUANTIFICATION_WINDOW_CENTER, --cleavage_offset QUANTIFICATION_WINDOW_CENTER
#                         Center of quantification window to use within respect to the 3' end of the provided sgRNA sequence. Remember that the sgRNA sequence must
#                         be entered without the PAM. For cleaving nucleases, this is the predicted cleavage position. The default is -3 and is suitable for the
#                         Cas9 system. For alternate nucleases, other cleavage offsets may be appropriate, for example, if using Cpf1 this parameter would be set to
#                         1. For base editors, this could be set to -17 to only include mutations near the 5' end of the sgRNA. Multiple quantification window
#                         centers (corresponding to each guide specified by --guide_seq) can be specified with a comma-separated list. (default: -3)
#   --exclude_bp_from_left EXCLUDE_BP_FROM_LEFT
#                         Exclude bp from the left side of the amplicon sequence for the quantification of the indels (default: 15)
#   --exclude_bp_from_right EXCLUDE_BP_FROM_RIGHT
#                         Exclude bp from the right side of the amplicon sequence for the quantification of the indels (default: 15)
#   --use_legacy_insertion_quantification
#                         If set, the legacy insertion quantification method will be used (i.e. with a 1bp quantification window, indels at the cut site and 1bp
#                         away from the cut site would be quantified). By default (if this parameter is not set) with a 1bp quantification window, only insertions
#                         at the cut site will be quantified. (default: False)
#   --ignore_substitutions
#                         Ignore substitutions events for the quantification and visualization (default: False)
#   --ignore_insertions   Ignore insertions events for the quantification and visualization (default: False)
#   --ignore_deletions    Ignore deletions events for the quantification and visualization (default: False)
#   --discard_indel_reads
#                         Discard reads with indels in the quantification window from analysis (default: False)
#   --needleman_wunsch_gap_open NEEDLEMAN_WUNSCH_GAP_OPEN
#                         Gap open option for Needleman-Wunsch alignment (default: -20)
#   --needleman_wunsch_gap_extend NEEDLEMAN_WUNSCH_GAP_EXTEND
#                         Gap extend option for Needleman-Wunsch alignment (default: -2)
#   --needleman_wunsch_gap_incentive NEEDLEMAN_WUNSCH_GAP_INCENTIVE
#                         Gap incentive value for inserting indels at cut sites (default: 1)
#   --needleman_wunsch_aln_matrix_loc NEEDLEMAN_WUNSCH_ALN_MATRIX_LOC
#                         Location of the matrix specifying substitution scores in the NCBI format (see ftp://ftp.ncbi.nih.gov/blast/matrices/) (default: EDNAFULL)
#   --plot_histogram_outliers
#                         If set, all values will be shown on histograms. By default (if unset), histogram ranges are limited to plotting data within the 99
#                         percentile. (default: False)
#   --plot_window_size PLOT_WINDOW_SIZE, --offset_around_cut_to_plot PLOT_WINDOW_SIZE
#                         Defines the size of the window extending from the quantification window center to plot. Nucleotides within plot_window_size of the
#                         quantification_window_center for each guide are plotted. (default: 20)
#   --min_frequency_alleles_around_cut_to_plot MIN_FREQUENCY_ALLELES_AROUND_CUT_TO_PLOT
#                         Minimum % reads required to report an allele in the alleles table plot. (default: 0.2)
#   --expand_allele_plots_by_quantification
#                         If set, alleles with different modifications in the quantification window (but not necessarily in the plotting window (e.g. for another
#                         sgRNA)) are plotted on separate lines, even though they may have the same apparent sequence. To force the allele plot and the allele table
#                         to be the same, set this parameter. If unset, all alleles with the same sequence will be collapsed into one row. (default: False)
#   --allele_plot_pcts_only_for_assigned_reference
#                         If set, in the allele plots, the percentages will show the percentage as a percent of reads aligned to the assigned reference. Default
#                         behavior is to show percentage as a percent of all reads. (default: False)
#   -qwc QUANTIFICATION_WINDOW_COORDINATES, --quantification_window_coordinates QUANTIFICATION_WINDOW_COORDINATES
#                         Bp positions in the amplicon sequence specifying the quantification window. This parameter overrides values of the "--
#                         quantification_window_center", "--cleavage_offset", "--window_around_sgrna" or "--window_around_sgrna" values. Any indels/substitutions
#                         outside this window are excluded. Indexes are 0-based, meaning that the first nucleotide is position 0. Ranges are separted by the dash
#                         sign (e.g. "start-stop"), and multiple ranges can be separated by the underscore (_). A value of 0 disables this filter. (can be comma-
#                         separated list of values, corresponding to amplicon sequences given in --amplicon_seq e.g. 5-10,5-10_20-30 would specify the 5th-10th bp
#                         in the first reference and the 5th-10th and 20th-30th bp in the second reference) (default: None)
#   --annotate_wildtype_allele ANNOTATE_WILDTYPE_ALLELE
#                         Wildtype alleles in the allele table plots will be marked with this string (e.g. **). (default: )
#   --keep_intermediate   Keep all the intermediate files (default: False)
#   --dump                Dump numpy arrays and pandas dataframes to file for debugging purposes (default: False)
#   --write_detailed_allele_table
#                         If set, a detailed allele table will be written including alignment scores for each read sequence. (default: False)
#   --fastq_output        If set, a fastq file with annotations for each read will be produced. (default: False)
#   --bam_output          If set, a bam file with alignments for each read will be produced. (default: False)
#   -x BOWTIE2_INDEX, --bowtie2_index BOWTIE2_INDEX
#                         Basename of Bowtie2 index for the reference genome (default: )
#   --zip_output          If set, the output will be placed in a zip folder. (default: False)
#   --max_rows_alleles_around_cut_to_plot MAX_ROWS_ALLELES_AROUND_CUT_TO_PLOT
#                         Maximum number of rows to report in the alleles table plot. (default: 50)
#   --suppress_report     Suppress output report (default: False)
#   --place_report_in_output_folder
#                         If true, report will be written inside the CRISPResso output folder. By default, the report will be written one directory up from the
#                         report output. (default: False)
#   --suppress_plots      Suppress output plots (default: False)
#   --base_editor_output  Outputs plots and tables to aid in analysis of base editor studies. (default: False)
#   --conversion_nuc_from CONVERSION_NUC_FROM
#                         For base editor plots, this is the nucleotide targeted by the base editor (default: C)
#   --conversion_nuc_to CONVERSION_NUC_TO
#                         For base editor plots, this is the nucleotide produced by the base editor (default: T)
#   --prime_editing_pegRNA_spacer_seq PRIME_EDITING_PEGRNA_SPACER_SEQ
#                         pegRNA spacer sgRNA sequence used in prime editing. The spacer should not include the PAM sequence. The sequence should be given in the
#                         RNA 5'->3' order, so for Cas9, the PAM would be on the right side of the given sequence. (default: )
#   --prime_editing_pegRNA_extension_seq PRIME_EDITING_PEGRNA_EXTENSION_SEQ
#                         Extension sequence used in prime editing. The sequence should be given in the RNA 5'->3' order, such that the sequence starts with the RT
#                         template including the edit, followed by the Primer-binding site (PBS). (default: )
#   --prime_editing_pegRNA_extension_quantification_window_size PRIME_EDITING_PEGRNA_EXTENSION_QUANTIFICATION_WINDOW_SIZE
#                         Quantification window size (in bp) at flap site for measuring modifications anchored at the right side of the extension sequence. Similar
#                         to the --quantification_window parameter, the total length of the quantification window will be 2x this parameter. Default: 5bp (10bp
#                         total window size) (default: 5)
#   --prime_editing_pegRNA_scaffold_seq PRIME_EDITING_PEGRNA_SCAFFOLD_SEQ
#                         If given, reads containing any of this scaffold sequence before extension sequence (provided by --prime_editing_extension_seq) will be
#                         classified as 'Scaffold-incorporated'. The sequence should be given in the 5'->3' order such that the RT template directly follows this
#                         sequence. A common value is 'GGCACCGAGUCGGUGC'. (default: )
#   --prime_editing_pegRNA_scaffold_min_match_length PRIME_EDITING_PEGRNA_SCAFFOLD_MIN_MATCH_LENGTH
#                         Minimum number of bases matching scaffold sequence for the read to be counted as 'Scaffold-incorporated'. If the scaffold sequence matches
#                         the reference sequence at the incorporation site, the minimum number of bases to match will be minimally increased (beyond this parameter)
#                         to disambiguate between prime-edited and scaffold-incorporated sequences. (default: 1)
#   --prime_editing_nicking_guide_seq PRIME_EDITING_NICKING_GUIDE_SEQ
#                         Nicking sgRNA sequence used in prime editing. The sgRNA should not include the PAM sequence. The sequence should be given in the RNA
#                         5'->3' order, so for Cas9, the PAM would be on the right side of the sequence (default: )
#   --prime_editing_override_prime_edited_ref_seq PRIME_EDITING_OVERRIDE_PRIME_EDITED_REF_SEQ
#                         If given, this sequence will be used as the prime-edited reference sequence. This may be useful if the prime-edited reference sequence has
#                         large indels or the algorithm cannot otherwise infer the correct reference sequence. (default: )
#   --prime_editing_override_sequence_checks
#                         If set, checks to assert that the prime editing guides and extension sequence are in the proper orientation are not performed. This may be
#                         useful if the checks are failing inappropriately, but the user is confident that the sequences are correct. (default: False)
#   --crispresso1_mode    Parameter usage as in CRISPResso 1 (default: False)
#   --dsODN DSODN         Label reads with the dsODN sequence provided (default: )
#   --auto                Infer amplicon sequence from most common reads (default: False)
#   --debug               Show debug messages (default: False)
#   --no_rerun            Don't rerun CRISPResso2 if a run using the same parameters has already been finished. (default: False)
#   -p N_PROCESSES, --n_processes N_PROCESSES
#                         Specify the number of processes to use for analysis. Please use with caution since increasing this parameter will significantly increase
#                         the memory required to run CRISPResso. Can be set to 'max'. (default: 1)
#   --bam_input BAM_INPUT
#                         Aligned reads for processing in bam format (default: )
#   --bam_chr_loc BAM_CHR_LOC
#                         Chromosome location in bam for reads to process. For example: "chr1:50-100" or "chrX". (default: )
