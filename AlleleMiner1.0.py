import sys
import subprocess
import os.path
import pathlib
import re
import glob
import shutil
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from os.path import basename
from collections import defaultdict
import csv
from datasketch import MinHash

# parser settings
parser = argparse.ArgumentParser(
    prog='AlleleMiner',
    usage='Software for automatically phasing gene sequences',
    add_help=True,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

# Setting arguments
# FASTA of amino acid sequences of reference organisms
parser.add_argument('--pro',
                    help='Specify the FASTA of amino acid sequences of organism that become reference.',
                    required=False,
                    metavar=""
                    )
# GFF file of the reference organism
parser.add_argument('--gff',
                    help='Specify the GFF3 file of organism that become reference.',
                    required=True,
                    metavar=""
                    )
# Whole genome FASTA of the reference organism's base sequence
parser.add_argument('--base',
                    help='Specify the FASTA of whole genome sequence of organism that become reference.',
                    required=True,
                    metavar=""
                    )
# FASTQ files of leading organisms
parser.add_argument('--read',
                    help='Specify the read of organism. FASTQ is required as file type.',
                    required=True,
                    metavar=""
                    )

# Specifying the number of threads
parser.add_argument('-t', '--threads',
                    help='Number of CPU threads used by AlleleMiner',
                    default=2,
                    metavar="")

# Specifying the name of the organism to be analyzed
parser.add_argument('-N', '--name',
                    help='Specify the name of the organism to be analyzed.',
                    default='organism',
                    metavar=""
                    )

# Specifying the ploidy of the target organism to be analyzed
parser.add_argument('-P', '--ploidy',
                    help='Specify the ploidy of the target organism.',
                    required=True,
                    metavar=""
                    )

# Whether to use BUSCO
parser.add_argument('-b', '--busco',
                    help='Use busco to handle only genes common to land plants. '
                         '[ex) -b embryophyta_odb10] '
                         'Please refer to the BUSCO database to see which orthologs to analyze.',
                    type=str,
                    required=False,
                    metavar=""
                    )

# Setting the mapping quality for Minimap2 for the Read-to-Reference Mapping
parser.add_argument('--minimap_Q',
                    help='Specification of the minimum level of mapping quality in Minimap2:'
                         'The higher the value, the higher the quality of the mapping results; however,'
                         ' a value greater than 60 may result in no mapping results.',
                    default=60,
                    metavar="")

# Setting the quality of the Read-to-Reference Mapping
parser.add_argument('--map_q1',
                    help='Mapping quality:'
                         'This is used in Read-to-Reference Mapping. '
                         'Higher values indicate higher quality, but may decrease the number of outputs.',
                    default=0.95,
                    metavar=""
                    )

# Setting the quality of the Phased_sequence-to-Reference Mapping
parser.add_argument('--map_q2',
                    help='Mapping quality:'
                         'It will be used when modifying the PAF file output from the Phased_sequence-to-Reference Mapping. '
                         'Decreasing the value increases the likelihood of obtaining appropriate results '
                         'but may reduce the number of outputs. Increasing the value has the opposite effect.',
                    default=20,
                    metavar=""
                    )

# Read-to-Reference Mapping minimap2 k-mer
parser.add_argument('--k_mini_1',
                    help='Specify the k-mer of Minimap2 for Read-to-Reference Mapping.',
                    default=15,
                    metavar="")

# Phased_sequence-to-Reference Mapping minimap2 k-mer
parser.add_argument('--k_mini_2',
                    help='Specify the k-mer of Minimap2 for Phased_seqence-to-Reference Mapping.',
                    default=15,
                    metavar="")

# Hifiasm k-mer
parser.add_argument('--k_hifiasm',
                    help='Specify the k-mer of Hifiasm.',
                    default=51,
                    metavar="")

# Use only flye when assembling
parser.add_argument('-f', '--only_flye',
                    help='Only flye will be used in assembly.',
                    action='store_true'
                    )

# Using "--keep-haplotypes" with Flye
parser.add_argument('--hap_solved',
                    help='When you want to use the "--keep-haplotypes" option in Flye.',
                    action='store_true'
                    )

# Specify the input sequencing type in Flye
parser.add_argument('--flye_input',
                    help='You can specify the type of input long reads.'
                         'Please refer to the Flye manual. When specifying type of input long reads, '
                         'remove the initial --.',
                    default='pacbio-hifi',
                    metavar=""
                    )

# Phased sequences extraction based on user settings
parser.add_argument('-x', '--x_cut_out',
                    help='Cut out user-defined upstream and downstream regions of a gene from the phased sequences.'
                         'The upstream area is cut out according to the first value specified, '
                         'and the downstream area is cut out according to the second value specified.'
                         '[ex) -x 200 100] In this case, 200 bp upstream is cut out and 100 bp downstream',
                    type=int,
                    nargs=2,
                    metavar="")

# Extract phased sequences without adjacent genes (only possible when using GFF files containing all gene information)
parser.add_argument('-n', '--no_neighbor',
                    help='Extract phased sequences without adjacent genes. '
                         'This option is available only if the input file is a GFF file '
                         'that contains information on all genes of an organism. '
                         'in this case, it can also be used with the -b option.',
                    action='store_true'
                    )

# Perform alignment
parser.add_argument('-a', '--align',
                    help='When this option is specified, alignment will be performed. When using this option,'
                         'three arguments must be provided.'
                         'The first argument specifies the alignment tool, which must be '
                         'one of "Muscle", "PRANK", or "Clustal".'
                         'The second argument determines whether to align the gene region +2Kbp sequence. '
                         'If "Y" is specified, alignment is performed; if "N" is specified, alignment is not performed.'
                         'The third argument determines whether to align the gene region +user-defined sequence.'
                         'If "xY" is specified, alignment is performed; if "xN" is specified, alignment is not performed.'
                         'For example, if the alignment tool is "Muscle", '
                         'alignment of the gene region +2Kbp sequence is not performed, '
                         'and alignment of the gene region +user-specified sequence is performed, '
                         'the command should be written as \n[ex) -a Muscle N xY].',
                    type=str,
                    nargs=3,
                    metavar="")

# When muscle is selected in the "--align" option
parser.add_argument('--muscle_setting',
                    help='When the "-a" option is used with Muscle, it becomes a required option. '
                         'Two arguments must be specified when using this option.'
                         'The first argument determines whether to use "align" or "super5" in Muscle.'
                         'The second argument specifies an alternative alignment tool in case Muscle fails to run. '
                         'You can specify either "PRANK" or "Clustal".'
                         '[ex) --muscle_settings align PRANK].',
                    type=str,
                    nargs=2,
                    metavar="")

# K-mer for MinHash
parser.add_argument('-K', '--Kmer',
                    help='k-mer used to convert phased sequences to MinHash values.'
                         'The larger the k-mer, the more sensitive it is to differences between sequences.',
                    default=15,
                    metavar="")

# num_perm for MinHash
parser.add_argument('--num_perm',
                    help='The number of hash functions used for MinHash calculation. '
                         'A larger value improves the accuracy of Jaccard similarity estimation but '
                         'increases computation time and memory usage. Recommended values include 128, 256, 512, etc.',
                    default=128,
                    metavar="")


args = parser.parse_args()

# Assign arguments to variables
protein_fasta = args.pro
GFF = args.gff
base_fasta = args.base
read = args.read
busco_db = str(args.busco)
organism = args.name
ploidy = int(args.ploidy)
map_q1 = float(args.map_q1)
map_q2 = float(args.map_q2)
minimap_Q = float(args.minimap_Q)
threads = int(args.threads)
k_m_1 = int(args.k_mini_1)
k_m_2 = int(args.k_mini_2)
k_h = int(args.k_hifiasm)
flye_input = str(args.flye_input)
K = int(args.Kmer)
num_perm = int(args.num_perm)


# Get the analysis execution directory
pwd_dir = os.getcwd()
os.mkdir(f'{pwd_dir}/{organism}')
running_dir = f'{pwd_dir}/{organism}'

# Conditions for stopping analysis.
# If the FASTA amino acid sequences is not specified when BUSCO is specified, the analysis will be stopped.
if args.busco:
    if not args.pro:
        print('AlleleMiner Error: If you use -b option, you must also specify the --pro option.')
        sys.exit()

# Eliminate conflicts related to alignment
if args.align:
    tool, align_2kb, x_option_align = args.align
    if tool == 'Muscle' or tool == 'PRANK' or tool == 'Clustal':
        pass
    else:
        print('AlleleMiner Error: The first argument of the -a option must be either “Muscle” or “PRANK” or "Clustal".')
        sys.exit()
    if tool == 'Muscle':
        if not args.muscle_setting:
            print('AlleleMiner Error: If you chose "Muscle" for the -a option, you must also specify the --muscle_setting option.')
            sys.exit()
    if align_2kb == 'Y' or align_2kb == 'N':
        pass
    else:
        print('AlleleMiner Error: The second argument of the -a option must be either “Y” or “N”.')
        sys.exit()
    if x_option_align == 'xY' or x_option_align == 'xN':
        if x_option_align == 'xY':
            if not args.x_cut_out:
                print('AlleleMiner Error: If you chose "xY" for the -a option, you must also specify the -x option.')
                sys.exit()
    else:
        print('AlleleMiner Error: The third argument of the -a option must be either “xY” or “xN”.')
        sys.exit()

if args.muscle_setting:
    Command, sub = args.muscle_setting
    if not args.align:
        print('AlleleMiner Error: If you use --muscle_setting option, you must also specify the -a option')
        sys.exit()
    if Command == 'align' or Command == 'super5':
        pass
    else:
        print('AlleleMiner Error: The second argument of the --muscle_setting option must be either “align” or “super5”.')
        sys.exit()
    if sub == 'PRANK' or sub == 'Clustal' or sub == 'NONE':
        pass
    else:
        print('AlleleMiner Error: The fifth argument of the --muscle_setting option must be either "PRANK" or "Clustal" or "NONE".')
        sys.exit()

# When the "-b" option is specified.
# Run BUSCO.
if args.busco:
    output_dir_busco = f'busco_result'
    cmd_busco = f'busco -i {protein_fasta} -l {busco_db} -o {output_dir_busco} -m protein'
    subprocess.run(cmd_busco, shell=True)

    # Combine busco_singlecopy_seq into a single file
    busco_singlecopy_seq = f'{running_dir}/busco_result/run_{busco_db}/busco_sequences/single_copy_busco_sequences'
    print('make busco_singlecopy fasta')
    os.mkdir(f'{running_dir}/busco_pro_fa')
    busco_singlecopy_seq_paths = []
    busco_singlecopy_seq_files = glob.glob(f'{busco_singlecopy_seq}/*')
    for file in busco_singlecopy_seq_files:
        busco_singlecopy_seq_paths.append(file)
        # Put all paths in busco_singlecopy_seq_files into the list of busco_singlecopy_seq_paths
    busco_singlecopy_file = pathlib.Path(f'{running_dir}/busco_pro_fa/busco_singlecopy_pro.fa')
    busco_singlecopy_file.touch()  # Creating an empty busco_singlecopy.fa file
    with open(busco_singlecopy_file, 'w') as busco_singlecopy_file_in:
        for path in busco_singlecopy_seq_paths:
            with open(path) as path_in:
                for line in path_in:
                    busco_singlecopy_file_in.write(line)  # BUSCO's single copy fa combined into one
    print('DONE')

    # Extract only the ID from busco_singlecopy.fa
    print('make busco_singlecopy ID txt')
    busco_ID_file = pathlib.Path(f'{running_dir}/busco_pro_fa/busco_ID.txt')
    busco_ID_file.touch()  # Creating an empty busco_only_ID.txt file
    busco_ID_list = []
    cmd_seqkit = f'seqkit seq {busco_singlecopy_file} -n -i -o {busco_ID_file}'
    subprocess.run(cmd_seqkit, shell=True)  # Extract only the ID
    with open(busco_ID_file) as busco_ID_file_in:
        for line in busco_ID_file_in:
            Line = line.rstrip('\n')
            busco_ID_list.append(Line)  # Add ID to busco_ID_list
    print('DONE')

# Convert GFF files to GFF containing only genes or mRNA
print('make only_gene.gff')
os.mkdir(f'{running_dir}/all_output_gff')
only_gene_backup_gff = pathlib.Path(f'{running_dir}/all_output_gff/only_gene_backup.gff')
only_gene_backup_gff.touch()
pattern_gene = 'gene'
pattern_mrna = 'mRNA'
gff_path = f'{GFF}'
switch = 0
with open(gff_path, 'r') as gff:
    with open(only_gene_backup_gff, 'w') as only_gene_gff_in:
        for line in gff:
            fields = line.split()
            if len(fields) > 2:
                if fields[2] == pattern_gene:
                    only_gene_gff_in.write(line)
                    switch += 1
                elif switch == 0 and fields[2] == pattern_mrna:
                    only_gene_gff_in.write(line)
print('DONE')

# If it is written as “mRNA” instead of “gene,” correct it
print('make only_gene_backup.gff and only_gene.gff')
only_gene_gff = pathlib.Path(f'{running_dir}/all_output_gff/only_gene.gff')
only_gene_gff.touch()
shutil.copy(only_gene_backup_gff, only_gene_gff)
with open(only_gene_gff, encoding='cp932') as single_copy_gene_gff_in:
    gff_lines = single_copy_gene_gff_in.read()
gff_lines = gff_lines.replace('mRNA', 'gene')
with open(only_gene_gff, mode='w', encoding='cp932') as only_gene_gff_in:
    only_gene_gff_in.write(gff_lines)
only_gene_gff = f'{running_dir}/all_output_gff/only_gene.gff'
print('DONE')

# Creating single_copy_gene.gff
if args.busco:
    busco_ID_file = f'{running_dir}/busco_pro_fa/busco_ID.txt'
    busco_ID_list = []
    with open(busco_ID_file) as busco_ID_file_in:
        for line in busco_ID_file_in:
            Line = line.rstrip('\n')
            busco_ID_list.append(Line)  # Add ID to busco_ID_list
    print('make single_copy_gene.gff')
    only_gene_gff = f'{running_dir}/all_output_gff/only_gene.gff'
    single_copy_gene_gff = pathlib.Path(f'{running_dir}/all_output_gff/single_copy_gene.gff')
    single_copy_gene_gff.touch()  # Creating an empty single_copy_gene.gff
    cnt = 0
    while cnt < len(busco_ID_list):
        pattern_singlegene = busco_ID_list[cnt]
        with open(only_gene_gff, 'r') as only_gene_gff_in:
            with open(single_copy_gene_gff, 'a') as single_copy_gene_gff_in:
                for line in only_gene_gff_in:
                    if re.search(pattern_singlegene, line):  # Compare busco_ID_list and only_gene_gff
                        single_copy_gene_gff_in.write(line)
        cnt = cnt + 1
    print('DONE')

# Extract only the headers from the whole genome FASTA files
headers = []
for record in SeqIO.parse(f'{base_fasta}', 'fasta'):
    headers.append(f'{record.id}\n')
with open(f'{running_dir}/all_output_gff/fa_list.txt', 'w') as f:
    f.writelines(headers)

# Compare gff and whole genome FASTA.

# If there are chromosome names in gff that are not in whole genome FASTA,
# delete those lines and store them in a different file.
# This is to prevent errors in "gff3toolkit".

headers = []
with open(f'{running_dir}/all_output_gff/fa_list.txt', 'r') as f:
    for line_txt in f:
        headers.append(line_txt.strip())  # List the headers again
# Perform gff selection
use_gff_list = []
no_use_gff_list = []
use_gff_list_single = []
no_use_gff_list_single = []

# Processing only_gene_gff
gff_file = f'{running_dir}/all_output_gff/only_gene.gff'
with open(gff_file) as gff_file_in:
    for line in gff_file_in:
        fields = line.split()
        if fields[0] in headers:
            use_gff_list.append(line)
        else:
            no_use_gff_list.append(line)
with open(f'{running_dir}/all_output_gff/only_gene_analysis.gff', 'w') as new_gff:
    new_gff.writelines(use_gff_list)
with open(f'{running_dir}/all_output_gff/only_gene_del.gff', 'w') as new_gff_del:
    new_gff_del.writelines(no_use_gff_list)

# When the "-b" option is specified
if args.busco:
    gff_file_single = f'{running_dir}/all_output_gff/single_copy_gene.gff'
    with open(gff_file_single) as gff_file_in:
        for line in gff_file_in:
            fields = line.split()
            if fields[0] in headers:
                use_gff_list_single.append(line)
            else:
                no_use_gff_list_single.append(line)
if args.busco:
    with open(f'{running_dir}/all_output_gff/single_copy_gene_analysis.gff', 'w') as new_gff_single:
        new_gff_single.writelines(use_gff_list_single)
    with open(f'{running_dir}/all_output_gff/single_copy_gene_del.gff', 'w') as new_gff_del_single:
        new_gff_del_single.writelines(no_use_gff_list_single)

# Extract gene sequences from the FASTA base sequence of the reference organism (creation of multi fasta).
if args.busco:
    new_gene_gff = f'{running_dir}/all_output_gff/single_copy_gene_analysis.gff'
    print('make ref_single_copy_gene_multi FASTA')
    os.mkdir(f'{running_dir}/ref_single_copy_gene_multi_fa')
    output_file_gff3toolkit = f'{running_dir}/ref_single_copy_gene_multi_fa/pre_single_copy_ref_multi'
else:
    new_gene_gff = f'{running_dir}/all_output_gff/only_gene_analysis.gff'
    print('make ref_gene_multi fasta')
    os.mkdir(f'{running_dir}/ref_gene_multi_fa')
    output_file_gff3toolkit = f'{running_dir}/ref_gene_multi_fa/pre_ref_multi'
cmd_gff3toolkit = \
    f'gff3_to_fasta -g {new_gene_gff} -f {base_fasta} -st gene -d complete -o {output_file_gff3toolkit}'
subprocess.run(cmd_gff3toolkit, shell=True)

# Creating a reference sequence.
# The header is “>chromosome_number;gene_ID.”
if args.busco:
    output_file_gff3toolkit = f'{running_dir}/ref_single_copy_gene_multi_fa/pre_single_copy_ref_multi_gene.fa'
else:
    output_file_gff3toolkit = f'{running_dir}/ref_gene_multi_fa/pre_ref_multi_gene.fa'

count_gff3toolkit_file = 0
analysis_file_gff3toolkit_list = []
with open(f'{output_file_gff3toolkit}', 'r') as gff3toolkit_file:
    for line in gff3toolkit_file:
        count_gff3toolkit_file += 1
with open(f'{output_file_gff3toolkit}', 'r') as gff3toolkit_file:
    gff3toolkit_file_line = gff3toolkit_file.read().splitlines()
    for num in range(0, count_gff3toolkit_file):
        if num % 2 == 0:
            gff3toolkit_gene_name = [gff3toolkit_file_line[num]]
            gff3toolkit_gene_name_new = gff3toolkit_gene_name[0].split(':')
            gff3toolkit_gene_name_new2 = gff3toolkit_gene_name_new[2].split('|')
            # print(gff3toolkit_gene_name_new)
            analysis_file_gff3toolkit_list.append(f'{gff3toolkit_gene_name_new[0]};{gff3toolkit_gene_name_new2[2]}\n')
        if num % 2 == 1:
            analysis_file_gff3toolkit_list.append(f'{gff3toolkit_file_line[num]}\n')
if args.busco:
    with open(f'{running_dir}/ref_single_copy_gene_multi_fa/single_copy_ref_multi_gene.fa', 'w') as f:
        f.writelines(analysis_file_gff3toolkit_list)
else:
    with open(f'{running_dir}/ref_gene_multi_fa/ref_multi_gene.fa', 'w') as f:
        f.writelines(analysis_file_gff3toolkit_list)
print('DONE')

# Creating a single FASTA (required for final output)
print('make single_fasta')
if args.busco:
    os.mkdir(f'{running_dir}/ref_single_copy_gene_multi_fa/ref_gene_single_fa')
    fa_dir = f'{running_dir}/ref_single_copy_gene_multi_fa/ref_gene_single_fa'
    multi_fasta = f'{running_dir}/ref_single_copy_gene_multi_fa/single_copy_ref_multi_gene.fa'
else:
    os.mkdir(f'{running_dir}/ref_gene_multi_fa/ref_gene_single_fa')
    fa_dir = f'{running_dir}/ref_gene_multi_fa/ref_gene_single_fa'
    multi_fasta = f'{running_dir}/ref_gene_multi_fa/ref_multi_gene.fa'
count = 0
with open(multi_fasta, 'r') as multi_fasta_in:
    for line in multi_fasta_in:
        count = count + 1
with open(multi_fasta, 'r') as multi_fasta_in:
    multi_fasta_lines = multi_fasta_in.read().splitlines()
for i in range(0, count - 1):
    if i % 2 == 0:
        multi_fasta_target_lines = multi_fasta_lines[i:i + 2]
        ref_gene_ID = multi_fasta_target_lines[0].split(';')[1].split('=')[1]
        for line in multi_fasta_target_lines:
            with open(f'{fa_dir}/ref_gene_{ref_gene_ID}.fa', 'a') as f:
                f.write(line + '\n')
    if i % 2 == 1:
        pass
print('DONE')

# Run Minimap2 (Read-to-Reference Mapping)
os.mkdir(f'{running_dir}/mapping_result')
if args.busco:
    ref_fa = f'{running_dir}/ref_single_copy_gene_multi_fa/single_copy_ref_multi_gene.fa'
else:
    ref_fa = f'{running_dir}/ref_gene_multi_fa/ref_multi_gene.fa'
output_minimap2 = f'{running_dir}/mapping_result/mapping_result.paf'
cmd_minimap2 = f'minimap2 -t {threads} -k {k_m_1} -K 50M --secondary=no {ref_fa} {read} > {output_minimap2}'  # 変数の設定
subprocess.run(cmd_minimap2, shell=True)
print(f'finished mapping')

# Divide the mapping results (PAF) (divide them in stages, first by chromosome and then by gene).
print('splitting mapping_result.paf')
print('phase1')
os.mkdir(f'{running_dir}/mapping_result_split')
os.mkdir(f'{running_dir}/all_output_gff/ch_gff')
if args.busco:
    gff_file_map_distinguish = f'{running_dir}/all_output_gff/single_copy_gene_analysis.gff'
else:
    gff_file_map_distinguish = f'{running_dir}/all_output_gff/only_gene_analysis.gff'
ch_list = []
dir_name_list = []
with open(f'{gff_file_map_distinguish}', 'r') as f:
    for line in f:
        fields = line.split()
        ch_list.append(fields[0])
ch_list = list(dict.fromkeys(ch_list))  # Create a list of chromosomes
ch_list = sorted(ch_list)
print(f'chromosome_list\n{ch_list}')

# Create folders for each chromosome
for ch in ch_list:
    os.mkdir(f'{running_dir}/mapping_result_split/{ch}')
    os.mkdir(f'{running_dir}/all_output_gff/ch_gff/{ch}')

# Create GFF files for each chromosome
with open(f'{gff_file_map_distinguish}', 'r') as f:
    for line in f:
        fields = line.split()
        for ch in ch_list:
            path = f'{running_dir}/mapping_result_split/{ch}'
            dir_name = basename(path)
            if fields[0] == dir_name:
                with open(f'{running_dir}/mapping_result_split/{ch}/gene_{ch}.gff', 'a') as file:
                    file.write(line)
# Create GFF (all genes) files for each chromosome in the all_output_gff folder
with open(f'{running_dir}/all_output_gff/only_gene_analysis.gff') as f:
    for line in f:
        fields = line.split()
        for ch in ch_list:
            path = f'{running_dir}/mapping_result_split/{ch}'
            dir_name = basename(path)
            if fields[0] == dir_name:
                with open(f'{running_dir}/all_output_gff/ch_gff/{ch}/gene_{ch}.gff', 'a') as file:
                    file.write(line)

# Split PAF files by chromosome
for ch in ch_list:
    ch_line = []
    path = f'{running_dir}/mapping_result_split/{ch}'
    dir_name = basename(path)
    with open(f'{running_dir}/mapping_result/mapping_result.paf', 'r') as f:
        for line_paf in f:
            fields_paf = line_paf.split()
            fields_paf_5_split = fields_paf[5].split(';')
            fields_paf_5_ch = fields_paf_5_split[0]
            if fields_paf_5_ch == dir_name:
                ch_line.append(line_paf)
    with open(f'{running_dir}/mapping_result_split/{ch}/mapping_result_{ch}.paf', 'w') as file:
        file.writelines(ch_line)

# Obtain gene IDs from GFF files for each chromosome and list them.
for ch in ch_list:
    GENE_ID_list = []
    with open(f'{running_dir}/mapping_result_split/{ch}/gene_{ch}.gff', 'r') as f:
        for line in f:
            fields = line.split()
            if str(';') in fields[8]:
                field_split = fields[8].split(';')[0]
                GENE_ID = field_split.split('=')[1]
                GENE_ID_list.append(GENE_ID)
            else:
                GENE_ID = field_split.split('=')[1]
                GENE_ID_list.append(GENE_ID)
    # print(GENE_ID_list)

    # Create folders for each gene and create GFF files for each gene.
    for num in range(0, len(GENE_ID_list)):
        os.mkdir(f'{running_dir}/mapping_result_split/{ch}/{GENE_ID_list[num]}')
        with open(f'{running_dir}/mapping_result_split/{ch}/{GENE_ID_list[num]}/{GENE_ID_list[num]}.gff', 'w') as file:
            with open(f'{running_dir}/mapping_result_split/{ch}/gene_{ch}.gff', 'r') as f:
                line = f.readlines()
                file.write(line[num])
print('phase1 finished')

# Split mapping results by gene.
# The name of the reference in PAF should be “chromosome;gene_ID”.
# Genes for which no mapping results were obtained should be written to a text file.
no_map = []
print('phase2')
for ch in ch_list:
    print(f'{ch} split')
    gene_dir_list = []
    ch_path = f'{running_dir}/mapping_result_split/{ch}'
    for List in os.listdir(ch_path):
        if os.path.isdir(os.path.join(ch_path, List)):
            gene_dir_list.append(List)
    gene_dir_list = sorted(gene_dir_list)

    for gene in gene_dir_list:
        gene_line = []
        print(f'splitting mapping_result.paf:{gene}')
        if os.path.exists(f'{running_dir}/mapping_result_split/{ch}/mapping_result_{ch}.paf'):
            with open(f'{running_dir}/mapping_result_split/{ch}/mapping_result_{ch}.paf', 'r') as f:
                for line_paf in f:
                    fields_paf = line_paf.split()
                    fields_paf_5_split = fields_paf[5].split(';')
                    fields_paf_5_ID = fields_paf_5_split[1]
                    fields_paf_5_gene = fields_paf_5_ID.split('=')[1]

                    if gene == fields_paf_5_gene:
                        gene_line.append(line_paf)
        with open(f'{running_dir}/mapping_result_split/{ch}/{gene}/{gene}.paf', 'w') as paf_file:
            paf_file.writelines(gene_line)
        if os.path.exists(f'{running_dir}/mapping_result_split/{ch}/{gene}/{gene}.paf'):
            if os.stat(f'{running_dir}/mapping_result_split/{ch}/{gene}/{gene}.paf').st_size == 0:
                no_map.append(f'{gene}\n')
                os.remove(f'{running_dir}/mapping_result_split/{ch}/{gene}/{gene}.paf')
with open(f'{running_dir}/mapping_result_split/no_mapping_list.txt', 'w') as no_txt:
    no_txt.writelines(no_map)
print('phase2 finished')

# Extract only high-quality mapping results
print('make only high quality mapping paf')
os.mkdir(f'{running_dir}/mapping_high_q')
low_q_list = []
for ch in ch_list:
    os.mkdir(f'{running_dir}/mapping_high_q/{ch}')
    gene_dir_list = []
    ch_path = f'{running_dir}/mapping_result_split/{ch}'
    for List in os.listdir(ch_path):
        if os.path.isdir(os.path.join(ch_path, List)):
            gene_dir_list.append(List)
    gene_dir_list = sorted(gene_dir_list)
    for gene in gene_dir_list:
        print(f'picking out {ch};{gene} high quality mapping')
        high_q_list = []
        if os.path.exists(f'{running_dir}/mapping_result_split/{ch}/{gene}/{gene}.paf'):
            with open(f'{running_dir}/mapping_result_split/{ch}/{gene}/{gene}.paf', 'r') as paf_file:
                for line in paf_file:
                    fields = line.split()
                    if float(fields[11]) >= minimap_Q and float(fields[8]) - float(fields[7]) > float(
                            fields[6]) * map_q1:
                        high_q_list.append(line)
                    with open(f'{running_dir}/mapping_high_q/{ch}/{gene}_hq.paf', 'w') as hq_paf_file:
                        hq_paf_file.writelines(high_q_list)
        if os.path.exists(f'{running_dir}/mapping_high_q/{ch}/{gene}_hq.paf'):
            if os.stat(f'{running_dir}/mapping_high_q/{ch}/{gene}_hq.paf').st_size == 0:
                low_q_list.append(f'{gene}\n')
                os.remove(f'{running_dir}/mapping_high_q/{ch}/{gene}_hq.paf')
with open(f'{running_dir}/mapping_high_q/no_hq_mapping_list.txt', 'w') as no_hq_txt:
    no_hq_txt.writelines(low_q_list)
print('DONE')

# Extract reads that straddled references
print('make straddled mapping paf')
os.mkdir(f'{running_dir}/mapping_straddled')
no_straddled_list = []
for ch in ch_list:
    os.mkdir(f'{running_dir}/mapping_straddled/{ch}')
    gene_dir_list = []
    ch_path = f'{running_dir}/mapping_result_split/{ch}'
    for List in os.listdir(ch_path):
        if os.path.isdir(os.path.join(ch_path, List)):
            gene_dir_list.append(List)
    gene_dir_list = sorted(gene_dir_list)
    for gene in gene_dir_list:
        print(f'picking out {ch};{gene} straddled mapping')
        straddled_list = []
        if os.path.exists(f'{running_dir}/mapping_high_q/{ch}/{gene}_hq.paf'):
            with open(f'{running_dir}/mapping_high_q/{ch}/{gene}_hq.paf', 'r') as paf_file:
                for line in paf_file:
                    fields = line.split()
                    if fields[4] == '+':  # When strand is "+"
                        if int(fields[1]) > int(fields[6]) and int(fields[2]) > int(fields[7]) and (
                                int(fields[1]) - int(fields[2])) > (int(fields[6]) - int(fields[7])):
                            straddled_list.append(line)
                    if fields[4] == '-':  # When strand is "-"
                        new_query_start = (int(fields[1]) - int(fields[3])) - 1
                        if int(fields[1]) > int(fields[6]) and int(new_query_start) > int(fields[7]) and (
                                int(fields[3]) + 1) > (int(fields[6]) - int(fields[7])):
                            straddled_list.append(line)
                    with open(f'{running_dir}/mapping_straddled/{ch}/{gene}_straddled.paf', 'w') as straddled_file:
                        straddled_file.writelines(straddled_list)
        if os.path.exists(f'{running_dir}/mapping_straddled/{ch}/{gene}_straddled.paf'):
            if os.stat(f'{running_dir}/mapping_straddled/{ch}/{gene}_straddled.paf').st_size == 0:
                no_straddled_list.append(f'{gene}\n')
                os.remove(f'{running_dir}/mapping_straddled/{ch}/{gene}_straddled.paf')
with open(f'{running_dir}/mapping_straddled/no_straddled_mapping_list.txt', 'w') as no_hq_txt:
    no_hq_txt.writelines(no_straddled_list)
print('DONE')


# Convert FASTQ reads to FASTA. Then, based on the PAF information,
# extract the read sequences from the converted FASTA in FASTA format for each gene.

# First, convert FASTQ to FASTA.
fa_file_name = os.path.basename(read).split('.', 1)[0]
print('converting read file from fastq to fasta')
read_fasta = f'{running_dir}/{fa_file_name}.fa'
cmd_seqkit2 = f'seqkit fq2fa {read} -o {read_fasta}'
subprocess.run(cmd_seqkit2, shell=True)
print('DONE')

# Extract only lead and gene IDs from PAF
print('extract ID')
os.mkdir(f'{running_dir}/read_extract')
os.mkdir(f'{running_dir}/read_extract/read_ID')
os.mkdir(f'{running_dir}/read_extract/ID')
for ch in ch_list:
    os.mkdir(f'{running_dir}/read_extract/read_ID/{ch}')
    os.mkdir(f'{running_dir}/read_extract/ID/{ch}')
    gene_dir_list = []
    ch_path = f'{running_dir}/mapping_result_split/{ch}'
    for List in os.listdir(ch_path):
        if os.path.isdir(os.path.join(ch_path, List)):
            gene_dir_list.append(List)
    gene_dir_list = sorted(gene_dir_list)
    for gene in gene_dir_list:
        print(f'extract {ch};{gene} ID')
        ID_list = []
        ID_ref_list = []
        if os.path.exists(f'{running_dir}/mapping_straddled/{ch}/{gene}_straddled.paf'):
            with open(f'{running_dir}/mapping_straddled/{ch}/{gene}_straddled.paf', 'r') as f:
                for line in f:
                    fields = line.split()
                    ID_list.append(f'{fields[0]}\n')
                    ID_ref_list.append(f'read:{fields[0]}#reference:{fields[5]}\n')
            with open(f'{running_dir}/read_extract/read_ID/{ch}/{gene}_read_ID.txt', 'w') as file_1:
                file_1.writelines(ID_list)
            with open(f'{running_dir}/read_extract/ID/{ch}/{gene}_ID.txt', 'w') as file_2:
                file_2.writelines(ID_ref_list)

# Extract read sequences for each gene from all read FASTA files in FASTA format based on read ID information.
print(f'extract read genome from {fa_file_name}.fa')
os.mkdir(f'{running_dir}/read_extract/read_fasta')
for ch in ch_list:
    os.mkdir(f'{running_dir}/read_extract/read_fasta/{ch}')
    gene_dir_list = []
    ch_path = f'{running_dir}/mapping_result_split/{ch}'
    for List in os.listdir(ch_path):
        if os.path.isdir(os.path.join(ch_path, List)):
            gene_dir_list.append(List)
    gene_dir_list = sorted(gene_dir_list)
    for gene in gene_dir_list:
        read_sequence_list = []
        if os.path.exists(f'{running_dir}/read_extract/read_ID/{ch}/{gene}_read_ID.txt'):
            print(f'this is {ch};{gene} analysis')
            input_file = f'{running_dir}/{fa_file_name}.fa'
            read_ID_list = f'{running_dir}/read_extract/read_ID/{ch}/{gene}_read_ID.txt'
            output_file = f'{running_dir}/read_extract/read_fasta/{ch}/read_{gene}.fa'
            cmd_seqkit3 = f'seqkit grep -w 0 -f {read_ID_list} {input_file} > {output_file}'
            subprocess.run(cmd_seqkit3, shell=True)

# Start assembly.
# When the "-f" command is used
if args.only_flye:
    print('assembly start')
    os.mkdir(f'{running_dir}/assembly_result')
# Without the "-f" command.
# Assemble the read with hifiasm.
else:
    print('Hybrid assembly start')
    print('Hifiasm start')
    os.mkdir(f'{running_dir}/assembly_result')
    os.mkdir(f'{running_dir}/assembly_result/hifiasm_result')
    os.mkdir(f'{running_dir}/assembly_result/hifiasm_result_hetero')
    os.mkdir(f'{running_dir}/assembly_result/hifiasm_result_homo')
    os.mkdir(f'{running_dir}/assembly_result/hifiasm_result_failed')
    for ch in ch_list:
        gene_dir_list = []
        os.mkdir(f'{running_dir}/assembly_result/hifiasm_result/{ch}')
        os.mkdir(f'{running_dir}/assembly_result/hifiasm_result_hetero/{ch}')
        os.mkdir(f'{running_dir}/assembly_result/hifiasm_result_homo/{ch}')
        os.mkdir(f'{running_dir}/assembly_result/hifiasm_result_failed/{ch}')
        ch_path = f'{running_dir}/mapping_result_split/{ch}'
        for List in os.listdir(ch_path):
            if os.path.isdir(os.path.join(ch_path, List)):
                gene_dir_list.append(List)
        gene_dir_list = sorted(gene_dir_list)
        for gene in gene_dir_list:
            if os.path.exists(f'{running_dir}/read_extract/read_fasta/{ch}/read_{gene}.fa'):
                os.mkdir(f'{running_dir}/assembly_result/hifiasm_result/{ch}/hifiasm_result_{gene}')
                input_file = f'{running_dir}/read_extract/read_fasta/{ch}/read_{gene}.fa'
                hifiasm_cmd = f'hifiasm -o {running_dir}/assembly_result/hifiasm_result/{ch}/hifiasm_result_{gene}/{gene} -t {threads} -k {k_h} -f 0 -l 3 -i --n-hap {ploidy} -n 0 {input_file}'
                subprocess.run(hifiasm_cmd, shell=True)
                # If the hap1 and hap2 files are not present, consider the assembly a failure and move the folder.
                # If the contents of hap1 or hap2 do not exist,
                # the assembly is considered a failure and the folder is moved.
                # Distinguishing between homo and hetero.
                if os.path.exists(f'{running_dir}/assembly_result/hifiasm_result/{ch}/hifiasm_result_{gene}/{gene}.bp.hap1.p_ctg.gfa') and \
                    os.path.exists(f'{running_dir}/assembly_result/hifiasm_result/{ch}/hifiasm_result_{gene}/{gene}.bp.hap2.p_ctg.gfa'):
                    if os.stat(f'{running_dir}/assembly_result/hifiasm_result/{ch}/hifiasm_result_{gene}/{gene}.bp.hap1.p_ctg.gfa').st_size == 0 or \
                            os.stat(f'{running_dir}/assembly_result/hifiasm_result/{ch}/hifiasm_result_{gene}/{gene}.bp.hap2.p_ctg.gfa').st_size == 0:
                        shutil.move(f'{running_dir}/assembly_result/hifiasm_result/{ch}/hifiasm_result_{gene}',
                                    f'{running_dir}/assembly_result/hifiasm_result_failed/{ch}')
                        print('assembly Hifiasm failed!' + '\n')
                    else:
                        with open(
                                f'{running_dir}/assembly_result/hifiasm_result/{ch}/hifiasm_result_{gene}/{gene}.bp.hap1.p_ctg.gfa') as f1:
                            line1_1 = f1.readline().rstrip()
                            fields_hap1 = line1_1.split()
                        with open(
                                f'{running_dir}/assembly_result/hifiasm_result/{ch}/hifiasm_result_{gene}/{gene}.bp.hap2.p_ctg.gfa') as f2:
                            line2_1 = f2.readline().rstrip()
                            fields_hap2 = line2_1.split()
                        if fields_hap1[2] == fields_hap2[2]:
                            shutil.move(f'{running_dir}/assembly_result/hifiasm_result/{ch}/hifiasm_result_{gene}',
                                        f'{running_dir}/assembly_result/hifiasm_result_homo/{ch}')
                        elif fields_hap1[2] != fields_hap2[2]:
                            shutil.move(f'{running_dir}/assembly_result/hifiasm_result/{ch}/hifiasm_result_{gene}',
                                        f'{running_dir}/assembly_result/hifiasm_result_hetero/{ch}')
                        print('assembly Hifiasm OK!' + '\n')
                else:
                    shutil.move(f'{running_dir}/assembly_result/hifiasm_result/{ch}/hifiasm_result_{gene}',
                                f'{running_dir}/assembly_result/hifiasm_result_failed/{ch}')
                    print('assembly Hifiasm failed!' + '\n')
        os.rmdir(f'{running_dir}/assembly_result/hifiasm_result/{ch}')
    os.rmdir(f'{running_dir}/assembly_result/hifiasm_result')
    print('hifism end')

# If the "-f" command is specified, assembly will start from here.
# If the -f command is not specified, reads that failed in Hifiasm will be assembled in Flye.
print('flye start')
os.mkdir(f'{running_dir}/assembly_result/flye_result')
os.mkdir(f'{running_dir}/assembly_result/flye_result_hetero')
os.mkdir(f'{running_dir}/assembly_result/flye_result_homo')
os.mkdir(f'{running_dir}/assembly_result/flye_result_failed')
for ch in ch_list:
    gene_dir_list = []
    os.mkdir(f'{running_dir}/assembly_result/flye_result/{ch}')
    os.mkdir(f'{running_dir}/assembly_result/flye_result_hetero/{ch}')
    os.mkdir(f'{running_dir}/assembly_result/flye_result_homo/{ch}')
    os.mkdir(f'{running_dir}/assembly_result/flye_result_failed/{ch}')
    ch_path = f'{running_dir}/mapping_result_split/{ch}'
    for List in os.listdir(ch_path):
        if os.path.isdir(os.path.join(ch_path, List)):
            gene_dir_list.append(List)
    gene_dir_list = sorted(gene_dir_list)
    for gene in gene_dir_list:
        if args.only_flye:  # If the "-f" command is specified
            if args.hap_solved:
                if os.path.exists(f'{running_dir}/read_extract/read_fasta/{ch}/read_{gene}.fa'):
                    input_file = f'{running_dir}/read_extract/read_fasta/{ch}/read_{gene}.fa'
                    flye_cmd = f'flye --{flye_input} {input_file} --out-dir {running_dir}/assembly_result/flye_result/{ch}/flye_result_{gene} -t {threads} --keep-haplotypes'
                    subprocess.run(flye_cmd, shell=True)
            else:
                if os.path.exists(f'{running_dir}/read_extract/read_fasta/{ch}/read_{gene}.fa'):
                    input_file = f'{running_dir}/read_extract/read_fasta/{ch}/read_{gene}.fa'
                    flye_cmd = f'flye --{flye_input} {input_file} --out-dir {running_dir}/assembly_result/flye_result/{ch}/flye_result_{gene} -t {threads}'
                    subprocess.run(flye_cmd, shell=True)
        else:
            if args.hap_solved:
                if os.path.exists(f'{running_dir}/assembly_result/hifiasm_result_failed/{ch}/hifiasm_result_{gene}'):
                    if os.path.exists(f'{running_dir}/read_extract/read_fasta/{ch}/read_{gene}.fa'):
                        input_file = f'{running_dir}/read_extract/read_fasta/{ch}/read_{gene}.fa'
                        flye_cmd = f'flye --{flye_input} {input_file} --out-dir {running_dir}/assembly_result/flye_result/{ch}/flye_result_{gene} -t {threads} --keep-haplotypes'
                        subprocess.run(flye_cmd, shell=True)
            else:
                if os.path.exists(f'{running_dir}/assembly_result/hifiasm_result_failed/{ch}/hifiasm_result_{gene}'):
                    if os.path.exists(f'{running_dir}/read_extract/read_fasta/{ch}/read_{gene}.fa'):
                        input_file = f'{running_dir}/read_extract/read_fasta/{ch}/read_{gene}.fa'
                        flye_cmd = f'flye --{flye_input} {input_file} --out-dir {running_dir}/assembly_result/flye_result/{ch}/flye_result_{gene} -t {threads}'
                        subprocess.run(flye_cmd, shell=True)
        # If assembly_graph.gfa does not exist, assembly will fail.
        # If the number of contigs exceeds {ploidy}, the assembly is considered a failure.
        if os.path.exists(f'{running_dir}/assembly_result/flye_result/{ch}/flye_result_{gene}/assembly_graph.gfa'):
            cnt = 0
            with open(f'{running_dir}/assembly_result/flye_result/{ch}/flye_result_{gene}/assembly_graph.gfa') as f:
                for line in f:
                    cnt += 1
            cnt = cnt - 1
            contig = cnt // 2
            if 1 < contig <= ploidy:
                shutil.move(f'{running_dir}/assembly_result/flye_result/{ch}/flye_result_{gene}',
                            f'{running_dir}/assembly_result/flye_result_hetero/{ch}')
                print('assembly Flye OK!' + '\n')
            elif contig == 1:
                shutil.move(f'{running_dir}/assembly_result/flye_result/{ch}/flye_result_{gene}',
                            f'{running_dir}/assembly_result/flye_result_homo/{ch}')
                print('assembly Flye OK!' + '\n')
            else:
                shutil.move(f'{running_dir}/assembly_result/flye_result/{ch}/flye_result_{gene}',
                            f'{running_dir}/assembly_result/flye_result_failed/{ch}')
                print('assembly Flye failed!' + '\n')
        else:
            if os.path.exists(f'{running_dir}/assembly_result/flye_result/{ch}/flye_result_{gene}'):
                shutil.move(f'{running_dir}/assembly_result/flye_result/{ch}/flye_result_{gene}',
                            f'{running_dir}/assembly_result/flye_result_failed/{ch}')
                print('assembly Flye failed!')
    os.rmdir(f'{running_dir}/assembly_result/flye_result/{ch}')
os.rmdir(f'{running_dir}/assembly_result/flye_result')
print('flye end')
if args.only_flye:
    print('assembly DONE')
else:
    print('hybrid assembly DONE')


# Count the number of files obtained in the assembly.
if os.path.exists(f'{running_dir}/assembly_result/hifiasm_result_hetero'):
    file_cnt = 0
    for ch in ch_list:
        elements = glob.glob(f'{running_dir}/assembly_result/hifiasm_result_hetero/{ch}/*')
        for element in elements:
            file_cnt += 1
    with open(f'{running_dir}/assembly_result/file_count.txt', 'a') as file:
        file.write(f'hifiasm_result_hetero = {file_cnt}\n')
if os.path.exists(f'{running_dir}/assembly_result/hifiasm_result_homo'):
    file_cnt = 0
    for ch in ch_list:
        elements = glob.glob(f'{running_dir}/assembly_result/hifiasm_result_homo/{ch}/*')
        for element in elements:
            file_cnt += 1
    with open(f'{running_dir}/assembly_result/file_count.txt', 'a') as file:
        file.write(f'hifiasm_result_homo = {file_cnt}\n')
if os.path.exists(f'{running_dir}/assembly_result/hifiasm_result_failed'):
    file_cnt = 0
    for ch in ch_list:
        elements = glob.glob(f'{running_dir}/assembly_result/hifiasm_result_failed/{ch}/*')
        for element in elements:
            file_cnt += 1
    with open(f'{running_dir}/assembly_result/file_count.txt', 'a') as file:
        file.write(f'hifiasm_result_failed = {file_cnt}\n')
if os.path.exists(f'{running_dir}/assembly_result/flye_result_hetero'):
    file_cnt = 0
    for ch in ch_list:
        elements = glob.glob(f'{running_dir}/assembly_result/flye_result_hetero/{ch}/*')
        for element in elements:
            file_cnt += 1
    with open(f'{running_dir}/assembly_result/file_count.txt', 'a') as file:
        file.write(f'flye_result_hetero = {file_cnt}\n')
if os.path.exists(f'{running_dir}/assembly_result/flye_result_homo'):
    file_cnt = 0
    for ch in ch_list:
        elements = glob.glob(f'{running_dir}/assembly_result/flye_result_homo/{ch}/*')
        for element in elements:
            file_cnt += 1
    with open(f'{running_dir}/assembly_result/file_count.txt', 'a') as file:
        file.write(f'flye_result_homo = {file_cnt}\n')
if os.path.exists(f'{running_dir}/assembly_result/flye_result_failed'):
    file_cnt = 0
    for ch in ch_list:
        elements = glob.glob(f'{running_dir}/assembly_result/flye_result_failed/{ch}/*')
        for element in elements:
            file_cnt += 1
    with open(f'{running_dir}/assembly_result/file_count.txt', 'a') as file:
        file.write(f'flye_result_failed = {file_cnt}\n')


# Final outputs.
# Extract the contig and arrange it in a phased sequence.
print('make phased sequence')
os.mkdir(f'{running_dir}/phased_seq')
os.mkdir(f'{running_dir}/phased_seq_sub')
os.mkdir(f'{running_dir}/phased_seq/FASTA')
os.mkdir(f'{running_dir}/phased_seq_sub/FASTA')
for ch in ch_list:
    gene_dir_list = []
    os.mkdir(f'{running_dir}/phased_seq/FASTA/{ch}')
    os.mkdir(f'{running_dir}/phased_seq_sub/FASTA/{ch}')
    ch_path = f'{running_dir}/mapping_result_split/{ch}'
    for List in os.listdir(ch_path):
        if os.path.isdir(os.path.join(ch_path, List)):
            gene_dir_list.append(List)
    gene_dir_list = sorted(gene_dir_list)
    for gene in gene_dir_list:
        if os.path.exists(f'/{running_dir}/assembly_result/flye_result_failed/{ch}/flye_result_{gene}/assembly_graph.gfa'):
            seq_cnt = 0
            with open(f'/{running_dir}/assembly_result/flye_result_failed/{ch}/flye_result_{gene}/assembly_graph.gfa', 'r') as file:
                for line in file:
                    if line.startswith('S'):
                        seq_cnt +=1
            with open(f'/{running_dir}/phased_seq_sub/FASTA/{ch}/phased_{gene}.fa', 'a') as f:
                hap = 1
                for i in range(seq_cnt):
                    alle_f = f'>{organism}#{ch}#{gene}#raw#hap{hap}#flye'
                    f.write(alle_f + '\n')
                    List = []
                    with open(
                            f'/{running_dir}/assembly_result/flye_result_failed/{ch}/flye_result_{gene}/assembly_graph.gfa',
                            'r') as file:
                        for line in file:
                            split_line = line.split()
                            List.append(split_line)
                        seq_f = List[hap][2]
                    f.write(seq_f + '\n')
                    hap += 1
        else:
            with open(f'/{running_dir}/phased_seq/FASTA/{ch}/phased_{gene}.fa', 'a') as f:
                if os.path.exists(
                        f'/{running_dir}/assembly_result/hifiasm_result_hetero/{ch}/hifiasm_result_{gene}'):
                    hap = 1
                    for i in range(2):
                        alle_h = f'>{organism}#{ch}#{gene}#raw#hap{hap}#hifiasm'
                        f.write(alle_h + '\n')
                        List = []
                        with open(
                                f'/{running_dir}/assembly_result/hifiasm_result_hetero/{ch}/hifiasm_result_{gene}/{gene}.bp.hap{hap}.p_ctg.gfa',
                                'r') as file:
                            for line in file:
                                split_line = line.split()
                                List.append(split_line)
                            seq_h = List[0][2]
                        f.write(seq_h + '\n')
                        hap += 1
                if os.path.exists(
                        f'/{running_dir}/assembly_result/hifiasm_result_homo/{ch}/hifiasm_result_{gene}'):
                    alle_h = f'>{organism}#{ch}#{gene}#raw#homo#hifiasm'
                    f.write(alle_h + '\n')
                    List = []
                    with open(
                            f'/{running_dir}/assembly_result/hifiasm_result_homo/{ch}/hifiasm_result_{gene}/{gene}.bp.hap1.p_ctg.gfa',
                            'r') as file:
                        for line in file:
                            split_line = line.split()
                            List.append(split_line)
                        seq_h = List[0][2]
                    f.write(seq_h + '\n')
                if os.path.exists(
                        f'/{running_dir}/assembly_result/flye_result_hetero/{ch}/flye_result_{gene}'):
                    hap = 1
                    for i in range(ploidy):
                        alle_f = f'>{organism}#{ch}#{gene}#raw#hap{hap}#flye'
                        f.write(alle_f + '\n')
                        List = []
                        with open(
                                f'/{running_dir}/assembly_result/flye_result_hetero/{ch}/flye_result_{gene}/assembly_graph.gfa') as file:
                            for line in file:
                                split_line = line.split()
                                List.append(split_line)
                            seq_f = List[hap][2]
                        f.write(seq_f + '\n')
                        hap += 1
                if os.path.exists(
                        f'/{running_dir}/assembly_result/flye_result_homo/{ch}/flye_result_{gene}'):
                    alle_f = f'>{organism}#{ch}#{gene}#raw#homo#flye'
                    f.write(alle_f + '\n')
                    List = []
                    with open(
                            f'/{running_dir}/assembly_result/flye_result_homo/{ch}/flye_result_{gene}/assembly_graph.gfa') as file:
                        for line in file:
                            split_line = line.split()
                            List.append(split_line)
                        seq_f = List[1][2]
                    f.write(seq_f + '\n')

        # Remove if the file is empty
        if os.path.exists(f'/{running_dir}/phased_seq/FASTA/{ch}/phased_{gene}.fa'):
            if os.stat(f'/{running_dir}/phased_seq/FASTA/{ch}/phased_{gene}.fa').st_size == 0:
                os.remove(f'{running_dir}/phased_seq/FASTA/{ch}/phased_{gene}.fa')

# Run Minimap2 again for final output (Phased_sequence-to-Reference Mapping).
# phased_seq
os.mkdir(f'{running_dir}/phased_seq/phased_seq_mapping_result')
for ch in ch_list:
    gene_dir_list = []
    os.mkdir(f'{running_dir}/phased_seq/phased_seq_mapping_result/{ch}')
    ch_path = f'{running_dir}/mapping_result_split/{ch}'
    for List in os.listdir(ch_path):
        if os.path.isdir(os.path.join(ch_path, List)):
            gene_dir_list.append(List)
    gene_dir_list = sorted(gene_dir_list)
    for gene in gene_dir_list:
        if os.path.exists(f'{running_dir}/phased_seq/FASTA/{ch}/phased_{gene}.fa'):
            if args.busco:
                ref_fa = f'{running_dir}/ref_single_copy_gene_multi_fa/ref_gene_single_fa/ref_gene_{gene}.fa'
            else:
                ref_fa = f'{running_dir}/ref_gene_multi_fa/ref_gene_single_fa/ref_gene_{gene}.fa'
            read = f'{running_dir}/phased_seq/FASTA/{ch}/phased_{gene}.fa'
            output_pre_ref = f'{running_dir}/phased_seq/phased_seq_mapping_result/{ch}/phased_seq_mapping_result_{gene}.paf'
            cmd_minimap2 = f'minimap2 -t {threads} -k {k_m_2} -K 50M --secondary=no {ref_fa} {read} > {output_pre_ref}'  # 変数の設定
            subprocess.run(cmd_minimap2, shell=True)
# phased_seq_sub
os.mkdir(f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result')
for ch in ch_list:
    gene_dir_list = []
    os.mkdir(f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/{ch}')
    ch_path = f'{running_dir}/mapping_result_split/{ch}'
    for List in os.listdir(ch_path):
        if os.path.isdir(os.path.join(ch_path, List)):
            gene_dir_list.append(List)
    gene_dir_list = sorted(gene_dir_list)
    for gene in gene_dir_list:
        if os.path.exists(f'{running_dir}/phased_seq_sub/FASTA/{ch}/phased_{gene}.fa'):
            if args.busco:
                ref_fa = f'{running_dir}/ref_single_copy_gene_multi_fa/ref_gene_single_fa/ref_gene_{gene}.fa'
            else:
                ref_fa = f'{running_dir}/ref_gene_multi_fa/ref_gene_single_fa/ref_gene_{gene}.fa'
            read = f'{running_dir}/phased_seq_sub/FASTA/{ch}/phased_{gene}.fa'
            output_pre_ref = f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/{ch}/phased_seq_sub_mapping_result_{gene}.paf'
            cmd_minimap2 = f'minimap2 -t {threads} -k {k_m_2} -K 50M --secondary=no {ref_fa} {read} > {output_pre_ref}'
            subprocess.run(cmd_minimap2, shell=True)

# Correct any mapping results that are incorrect.
# phased_seq
os.mkdir(f'{running_dir}/phased_seq/phased_seq_mapping_result/mapping_error_modifiable')
os.mkdir(f'{running_dir}/phased_seq/phased_seq_mapping_result/mapping_error_trash')
for ch in ch_list:
    gene_dir_list = []
    ch_path = f'{running_dir}/mapping_result_split/{ch}'
    for List in os.listdir(ch_path):
        if os.path.isdir(os.path.join(ch_path, List)):
            gene_dir_list.append(List)
    gene_dir_list = sorted(gene_dir_list)
    for gene in gene_dir_list:
        if os.path.exists(f'{running_dir}/phased_seq/FASTA/{ch}/phased_{gene}.fa' and
                          f'{running_dir}/phased_seq/phased_seq_mapping_result/{ch}/phased_seq_mapping_result_{gene}.paf'):
            if os.stat(
                    f'{running_dir}/phased_seq/phased_seq_mapping_result/{ch}/phased_seq_mapping_result_{gene}.paf').st_size != 0:
                con_fa_count = 0
                with open(f'{running_dir}/phased_seq/FASTA/{ch}/phased_{gene}.fa') as con_fa:
                    for line in con_fa:
                        con_fa_count += 1
                    con_fa_count = con_fa_count // 2
                con_paf_count = 0
                with open(
                        f'{running_dir}/phased_seq/phased_seq_mapping_result/{ch}/phased_seq_mapping_result_{gene}.paf') as con_paf:
                    for line in con_paf:
                        con_paf_count += 1
                if con_fa_count != con_paf_count:
                    print(f'phased_seq_mapping_result_{gene} error! modifying PAF')
                    modi_PAF = shutil.move(
                        f'{running_dir}/phased_seq/phased_seq_mapping_result/{ch}/phased_seq_mapping_result_{gene}.paf',
                        f'{running_dir}/phased_seq/phased_seq_mapping_result/mapping_error_modifiable')
                    good_map_list = []
                    with open(f'{modi_PAF}', 'r') as modi_PAF_in:
                        for line in modi_PAF_in:
                            fields = line.split()
                            if (int(fields[8]) - int(fields[7])) >= int(fields[6]) - map_q2:
                                good_map_list.append(line)
                    if len(good_map_list) == 0:
                        print(f'phased_seq_mapping_result_{gene} could not modify')
                        trash_PAF = shutil.move(
                            f'{running_dir}/phased_seq/phased_seq_mapping_result/mapping_error_modifiable/phased_seq_mapping_result_{gene}.paf',
                            f'{running_dir}/phased_seq/phased_seq_mapping_result/mapping_error_trash')
                    else:
                        with open(
                                f'{running_dir}/phased_seq/phased_seq_mapping_result/{ch}/phased_seq_mapping_result_{gene}.paf',
                                'w') as new_paf:
                            new_paf.writelines(good_map_list)
            else:
                print(f'phased_seq_mapping_result_{gene} no data')
                with open(f'{running_dir}/phased_seq/phased_seq_mapping_result/mapping_error_trash/trash_list',
                          'a') as file:
                    trash_PAF = shutil.move(
                        f'{running_dir}/phased_seq/phased_seq_mapping_result/{ch}/phased_seq_mapping_result_{gene}.paf',
                        f'{running_dir}/phased_seq/phased_seq_mapping_result/mapping_error_trash')
        # Tagging mapping results
        input_file = f'{running_dir}/phased_seq/phased_seq_mapping_result/{ch}/phased_seq_mapping_result_{gene}.paf'
        temp_file = f'{running_dir}/phased_seq/phased_seq_mapping_result/{ch}/phased_seq_mapping_result_{gene}_temp.paf'
        counter = defaultdict(int)
        if os.path.exists(input_file):
            with open(input_file) as fin, open(temp_file, 'w') as fout:
                for line in fin:
                    fields = line.strip().split('\t')
                    query = fields[0]
                    counter[query] += 1
                    unique_query = f"{query}*{counter[query]}"
                    fields[0] = unique_query
                    fout.write('\t'.join(fields) + '\n')
            os.replace(temp_file, input_file)  # Replace input_file with the contents of temp_file.

# phased_seq_sub
os.mkdir(f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/mapping_error_modifiable')
os.mkdir(f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/mapping_error_trash')
for ch in ch_list:
    gene_dir_list = []
    ch_path = f'{running_dir}/mapping_result_split/{ch}'
    for List in os.listdir(ch_path):
        if os.path.isdir(os.path.join(ch_path, List)):
            gene_dir_list.append(List)
    gene_dir_list = sorted(gene_dir_list)
    for gene in gene_dir_list:
        if os.path.exists(f'{running_dir}/phased_seq_sub/FASTA/{ch}/phased_{gene}.fa' and
                          f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/{ch}/phased_seq_sub_mapping_result_{gene}.paf'):
            if os.stat(
                    f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/{ch}/phased_seq_sub_mapping_result_{gene}.paf').st_size != 0:
                con_fa_count = 0
                with open(f'{running_dir}/phased_seq_sub/FASTA/{ch}/phased_{gene}.fa') as con_fa:
                    for line in con_fa:
                        con_fa_count += 1
                    con_fa_count = con_fa_count // 2
                con_paf_count = 0
                with open(
                        f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/{ch}/phased_seq_sub_mapping_result_{gene}.paf') as con_paf:
                    for line in con_paf:
                        con_paf_count += 1
                if con_fa_count != con_paf_count:
                    print(f'phased_seq_sub_mapping_result_{gene} error! modifying PAF')
                    modi_PAF = shutil.move(
                        f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/{ch}/phased_seq_sub_mapping_result_{gene}.paf',
                        f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/mapping_error_modifiable')
                    good_map_list = []
                    with open(f'{modi_PAF}', 'r') as modi_PAF_in:
                        for line in modi_PAF_in:
                            fields = line.split()
                            if (int(fields[8]) - int(fields[7])) >= int(fields[6]) - map_q2:
                                good_map_list.append(line)
                    if len(good_map_list) == 0:
                        print(f'phased_seq_sub_mapping_result_{gene} could not modify')
                        trash_PAF = shutil.move(
                            f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/mapping_error_modifiable/phased_seq_sub_mapping_result_{gene}.paf',
                            f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/mapping_error_trash')
                    else:
                        with open(
                                f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/{ch}/phased_seq_sub_mapping_result_{gene}.paf',
                                'w') as new_paf:
                            new_paf.writelines(good_map_list)
            else:
                print(f'phased_seq_sub_mapping_result_{gene} no data')
                with open(f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/mapping_error_trash/trash_list',
                          'a') as file:
                    trash_PAF = shutil.move(
                        f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/{ch}/phased_seq_sub_mapping_result_{gene}.paf',
                        f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/mapping_error_trash')
        # Tagging mapping results
        input_file = f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/{ch}/phased_seq_sub_mapping_result_{gene}.paf'
        temp_file = f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/{ch}/phased_seq_sub_mapping_result_{gene}_temp.paf'
        counter = defaultdict(int)
        if os.path.exists(input_file):
            with open(input_file) as fin, open(temp_file, 'w') as fout:
                for line in fin:
                    fields = line.strip().split('\t')
                    query = fields[0]
                    counter[query] += 1
                    unique_query = f"{query}*{counter[query]}"
                    fields[0] = unique_query
                    fout.write('\t'.join(fields) + '\n')
            os.replace(temp_file, input_file)  # Replace input_file with the contents of temp_file.

# Output the extracted sequences
phased_dir_list = ['phased_seq', 'phased_seq_sub']
for d in phased_dir_list:
    os.mkdir(f'{running_dir}/{d}/gene_region_FASTA')
    os.mkdir(f'{running_dir}/{d}/gene_+2Kb_FASTA')
    if args.x_cut_out:
        os.mkdir(f'{running_dir}/{d}/gene_+user_setting_FASTA')
    if args.no_neighbor:
        os.mkdir(f'{running_dir}/{d}/gene_except_neighbor_FASTA')
    for ch in ch_list:
        gene_dir_list = []
        os.mkdir(f'{running_dir}/{d}/gene_region_FASTA/{ch}')
        os.mkdir(f'{running_dir}/{d}/gene_+2Kb_FASTA/{ch}')
        if args.x_cut_out:
            os.mkdir(f'{running_dir}/{d}/gene_+user_setting_FASTA/{ch}')
        if args.no_neighbor:
            os.mkdir(f'{running_dir}/{d}/gene_except_neighbor_FASTA/{ch}')
        ch_path = f'{running_dir}/mapping_result_split/{ch}'
        for List in os.listdir(ch_path):
            if os.path.isdir(os.path.join(ch_path, List)):
                gene_dir_list.append(List)
        gene_dir_list = sorted(gene_dir_list)

        # List information in gff files (required when using the "-n" option)
        ch_gff_gene_list = []
        ch_gff_start_list = []
        ch_gff_end_list = []
        gff_strand_list = []
        with open(f'{running_dir}/all_output_gff/ch_gff/{ch}/gene_{ch}.gff', 'r') as ch_gff:
            for line_gff in ch_gff:
                fields_gff = line_gff.split()
                gene_name_gff = fields_gff[8].split(';')[0].split('=')[1]
                ch_gff_gene_list.append(gene_name_gff)
                ch_gff_start_list.append(int(fields_gff[3]))
                ch_gff_end_list.append(int(fields_gff[4]))
                gff_strand_list.append(fields_gff[6])

        for gene in gene_dir_list:
            if os.path.exists(f'{running_dir}/{d}/FASTA/{ch}/phased_{gene}.fa' and
                              f'{running_dir}/{d}/{d}_mapping_result/{ch}/{d}_mapping_result_{gene}.paf'):
                to_other_output_con_path = f'{running_dir}/{d}/FASTA/{ch}/phased_{gene}.fa'
                # Create a dictionary of {header: array}
                phased_seq_dict = {}
                with open(to_other_output_con_path, 'r') as seq_fasta:
                    lines = seq_fasta.readlines()
                    for line in range(0, len(lines), 2):
                        header = lines[line].strip()
                        seq = lines[line + 1].strip()
                        id_ = header.lstrip('>')
                        phased_seq_dict[id_] = seq  # Dictionary creation
                # Variable settings
                with open(
                        f'{running_dir}/{d}/{d}_mapping_result/{ch}/{d}_mapping_result_{gene}.paf') as con_paf:
                    for line in con_paf:
                        fields_con_paf = line.split()
                        consensus_name = fields_con_paf[0]

                        # Redefining header names
                        header_name_split = consensus_name.split('#')
                        header_name_gene = header_name_split[:]
                        header_name_2000 = header_name_split[:]
                        header_name_user = header_name_split[:]
                        header_name_neighbor = header_name_split[:]
                        header_name_gene[3] ='gene'
                        header_name_2000[3] = '2000'
                        header_name_user[3] = 'user'
                        header_name_neighbor[3] = 'neighbor'
                        header_rename_gene = '#'.join(header_name_gene)
                        header_rename_2000 = '#'.join(header_name_2000)
                        header_rename_user = '#'.join(header_name_user)
                        header_rename_neighbor = '#'.join(header_name_neighbor)

                        consensus_name_check = consensus_name.split('*')[0]
                        strand = fields_con_paf[4]
                        if consensus_name_check in phased_seq_dict:
                            edited_seq = phased_seq_dict[consensus_name_check]

                        # Extract only the gene region from the phased sequences.
                        if strand == '+':
                            gene_region_start = int(fields_con_paf[2]) - int(fields_con_paf[7])
                            # If gene_region_start becomes negative, set gene_region_start to 0 and list it.
                            if gene_region_start < 0:
                                gene_region_start = 0
                                with open(f'{running_dir}/{d}/no_gene_region_list', 'a') as no_region_list:
                                    no_region_list.write(f'{fields_con_paf[0]}#upstream_no_gene_region\n')
                            gene_region_end = int(fields_con_paf[3]) + (
                                    int(fields_con_paf[6]) - int(fields_con_paf[8]) - 1) + 1
                            if gene_region_end > int(fields_con_paf[1]) - 1:
                                gene_region_end = int(fields_con_paf[1]) - 1
                                with open(f'{running_dir}/{d}/no_gene_region_list', 'a') as no_region_list:
                                    no_region_list.write(f'{fields_con_paf[0]}#downstream_no_gene_region\n')
                        if strand == '-':
                            gene_region_start_minus = int(fields_con_paf[2]) - (
                                    int(fields_con_paf[6]) - int(fields_con_paf[8]))
                            # As with strand being "+"
                            if gene_region_start_minus < 0:
                                gene_region_start_minus = 0
                                with open(f'{running_dir}/{d}/no_gene_region_list', 'a') as no_region_list:
                                    no_region_list.write(f'{fields_con_paf[0]}#downstream_no_gene_region\n')
                            gene_region_end_minus = int(fields_con_paf[3]) + int(fields_con_paf[7])
                            if gene_region_end_minus > int(fields_con_paf[1]) - 1:
                                gene_region_end_minus = int(fields_con_paf[1]) - 1
                                with open(f'{running_dir}/{d}/no_gene_region_list', 'a') as no_region_list:
                                    no_region_list.write(f'{fields_con_paf[0]}#upstream_no_gene_region\n')
                        with open(
                                f'{running_dir}/{d}/gene_region_FASTA/{ch}/phased_gene_region_{gene}.fa',
                                'a') as new_file_gene:
                            new_file_gene.write('>' + header_rename_gene + '\n')
                            if strand == '+':
                                new_file_gene.write(
                                    edited_seq[gene_region_start:gene_region_end] + '\n')
                            # When strand is "-", make it reverse complementary.
                            if strand == '-':
                                seq = edited_seq[gene_region_start_minus:gene_region_end_minus]
                                seq_c = str(Seq(seq).reverse_complement())
                                new_file_gene.write(seq_c + '\n')

                        # Extract sequences containing 2kbp upstream and
                        # downstream from the gene region from the phased sequence.

                        # If upstream 2Kb or downstream 2KbP cannot be extracted,
                        # extract the sequence from the beginning or end of the phased sequence.

                        # Compile a list of genes
                        # for which upstream 2Kb or downstream 2Kb could not be extracted into a file (no_2Kb_list).

                        if strand == '+':
                            gene_from2Kb_start = (int(fields_con_paf[2]) - int(fields_con_paf[7])) - 2000
                            if gene_from2Kb_start < 0:
                                gene_from2Kb_start = 0
                                with open(f'{running_dir}/{d}/no_+2Kb_list', 'a') as no_2Kb_list:
                                    no_2Kb_list.write(f'{fields_con_paf[0]}#upstream_no_+2Kb\n')
                            gene_from2Kb_end = (int(fields_con_paf[3]) + (
                                    int(fields_con_paf[6]) - int(fields_con_paf[8]) - 1) + 1) + 2000
                            if gene_from2Kb_end > int(fields_con_paf[1]) - 1:
                                gene_from2Kb_end = int(fields_con_paf[1]) - 1
                                with open(f'{running_dir}/{d}/no_+2Kb_list', 'a') as no_2Kb_list:
                                    no_2Kb_list.write(f'{fields_con_paf[0]}#downstream_no_+2Kb\n')
                        if strand == '-':
                            gene_from2Kb_start_minus = (int(fields_con_paf[2]) - (
                                    int(fields_con_paf[6]) - int(fields_con_paf[8]))) - 2000
                            if gene_from2Kb_start_minus < 0:
                                gene_from2Kb_start_minus = 0
                                with open(f'{running_dir}/{d}/no_+2Kb_list', 'a') as no_2Kb_list:
                                    no_2Kb_list.write(f'{fields_con_paf[0]}#downstream_no_+2Kb\n')
                            gene_from2Kb_end_minus = (int(fields_con_paf[3]) + int(fields_con_paf[7])) + 2000
                            if gene_from2Kb_end_minus > int(fields_con_paf[1]) - 1:
                                gene_from2Kb_end_minus = int(fields_con_paf[1]) - 1
                                with open(f'{running_dir}/{d}/no_+2Kb_list', 'a') as no_2Kb_list:
                                    no_2Kb_list.write(f'{fields_con_paf[0]}#upstream_no_+2Kb\n')
                        with open(
                                f'{running_dir}/{d}/gene_+2Kb_FASTA/{ch}/phased_gene_+2Kb_{gene}.fa',
                                'a') as new_file_gene:
                            new_file_gene.write('>' + header_rename_2000 + '\n')
                            if strand == '+':
                                new_file_gene.write(
                                    edited_seq[gene_from2Kb_start:gene_from2Kb_end] + '\n')
                            if strand == '-':
                                seq = edited_seq[gene_from2Kb_start_minus:gene_from2Kb_end_minus]
                                seq_c = str(Seq(seq).reverse_complement())
                                new_file_gene.write(seq_c + '\n')

                        # Extract user-specified ranges from phased sequences.
                        if args.x_cut_out:
                            x_top, x_bottom = args.x_cut_out
                            if strand == '+':
                                gene_user_start = (int(fields_con_paf[2]) - int(fields_con_paf[7])) - x_top
                                if gene_user_start < 0:
                                    gene_user_start = 0
                                    with open(f'{running_dir}/{d}/no_+user_setting_list', 'a') as no_user_list:
                                        no_user_list.write(f'{fields_con_paf[0]}#upstream_no_+{x_top}\n')
                                gene_user_end = (int(fields_con_paf[3]) + (
                                        int(fields_con_paf[6]) - int(fields_con_paf[8]) - 1) + 1) + x_bottom
                                if gene_user_end > int(fields_con_paf[1]) - 1:
                                    gene_user_end = int(fields_con_paf[1]) - 1
                                    with open(f'{running_dir}/{d}/no_+user_setting_list', 'a') as no_user_list:
                                        no_user_list.write(f'{fields_con_paf[0]}#downstream_no_+{x_bottom}\n')
                            if strand == '-':
                                gene_user_start_minus = (int(fields_con_paf[2]) - (
                                        int(fields_con_paf[6]) - int(fields_con_paf[8]))) - x_bottom
                                if gene_user_start_minus < 0:
                                    gene_user_start_minus = 0
                                    with open(f'{running_dir}/{d}/no_+user_setting_list', 'a') as no_user_list:
                                        no_user_list.write(f'{fields_con_paf[0]}#downstream_no_{x_bottom}\n')
                                gene_user_end_minus = (int(fields_con_paf[3]) + int(fields_con_paf[7])) + x_top
                                if gene_user_end_minus > int(fields_con_paf[1]) - 1:
                                    gene_user_end_minus = int(fields_con_paf[1]) - 1
                                    with open(f'{running_dir}/{d}/no_+user_setting_list', 'a') as no_user_list:
                                        no_user_list.write(f'{fields_con_paf[0]}#upstream_no_{x_top}\n')
                            with open(
                                    f'{running_dir}/{d}/gene_+user_setting_FASTA/{ch}/phased_gene_+user_setting_{gene}.fa',
                                    'a') as new_file_gene:
                                new_file_gene.write('>' + header_rename_user + '\n')
                                if strand == '+':
                                    new_file_gene.write(
                                        edited_seq[gene_user_start:gene_user_end] + '\n')
                                if strand == '-':
                                    seq = edited_seq[gene_user_start_minus:gene_user_end_minus]
                                    seq_c = str(Seq(seq).reverse_complement())
                                    new_file_gene.write(seq_c + '\n')

                        # Extract sequences that do not contain adjacent genes from phased sequences
                        # (optional specification, only available for whole genome or BUSCO).
                        if args.no_neighbor:
                            if f'{gene}' in ch_gff_gene_list:
                                index_gene_list = ch_gff_gene_list.index(gene)
                                # print(index_gene_list)
                                len_gene_list = len(ch_gff_gene_list) - 1
                                # print(gene)
                                ch_gff_start = ch_gff_start_list[index_gene_list]
                                ch_gff_end = ch_gff_end_list[index_gene_list]
                                gff_strand = gff_strand_list[index_gene_list]
                                if index_gene_list == 0:
                                    pass
                                else:
                                    prev_ch_gff_end = ch_gff_end_list[index_gene_list - 1]
                                if index_gene_list == len_gene_list:
                                    pass
                                else:
                                    next_ch_gff_start = ch_gff_start_list[index_gene_list + 1]
                                if index_gene_list == 0:
                                    next_from_now = int(next_ch_gff_start) - int(ch_gff_end) - 1
                                    now_from_prev = 0
                                elif index_gene_list == len_gene_list:
                                    next_from_now = 0
                                    now_from_prev = int(ch_gff_start) - int(prev_ch_gff_end) - 1
                                else:
                                    next_from_now = int(next_ch_gff_start) - int(ch_gff_end) - 1
                                    now_from_prev = int(ch_gff_start) - int(prev_ch_gff_end) - 1

                                if strand == '+':
                                    gene_start = (int(fields_con_paf[2]) - int(fields_con_paf[7]))
                                    gene_end = (int(fields_con_paf[3]) + (
                                            int(fields_con_paf[6]) - int(fields_con_paf[8]) - 1) + 1)
                                    if gff_strand == '+':
                                        gene_no_neighbor_start = gene_start - now_from_prev
                                        if gene_no_neighbor_start < 0:
                                            gene_no_neighbor_start = 0
                                            with open(f'{running_dir}/{d}/no_except_neighbor_list',
                                                      'a') as no_no_neighbor_list:
                                                no_no_neighbor_list.write(
                                                    f'{fields_con_paf[0]}#upstream_no_except_neighbor\n')
                                        gene_no_neighbor_end = gene_end + next_from_now
                                        if gene_no_neighbor_end > int(fields_con_paf[1]) - 1:
                                            gene_no_neighbor_end = int(fields_con_paf[1]) - 1
                                            with open(f'{running_dir}/{d}/no_except_neighbor_list',
                                                      'a') as no_no_neighbor_list:
                                                no_no_neighbor_list.write(
                                                    f'{fields_con_paf[0]}#downstream_no_except_neighbor\n')
                                    if gff_strand == '-':
                                        gene_no_neighbor_start = gene_start - next_from_now
                                        if gene_no_neighbor_start < 0:
                                            gene_no_neighbor_start = 0
                                            with open(f'{running_dir}/{d}/no_except_neighbor_list',
                                                      'a') as no_no_neighbor_list:
                                                no_no_neighbor_list.write(
                                                    f'{fields_con_paf[0]}#upstream_no_except_neighbor\n')
                                        gene_no_neighbor_end = gene_end + now_from_prev
                                        if gene_no_neighbor_end > int(fields_con_paf[1]) - 1:
                                            gene_no_neighbor_end = int(fields_con_paf[1]) - 1
                                            with open(f'{running_dir}/{d}/no_except_neighbor_list',
                                                      'a') as no_no_neighbor_list:
                                                no_no_neighbor_list.write(
                                                    f'{fields_con_paf[0]}#downstream_no_except_neighbor\n')

                                if strand == '-':
                                    gene_start_minus = (int(fields_con_paf[2]) - (
                                            int(fields_con_paf[6]) - int(fields_con_paf[8])))
                                    gene_end_minus = (int(fields_con_paf[3]) + int(fields_con_paf[7]))
                                    if gff_strand == '+':
                                        gene_no_neighbor_start_minus = gene_start_minus - next_from_now
                                        if gene_no_neighbor_start_minus < 0:
                                            gene_no_neighbor_start_minus = 0
                                            with open(f'{running_dir}/{d}/no_except_neighbor_list',
                                                      'a') as no_no_neighbor_list:
                                                no_no_neighbor_list.write(
                                                    f'{fields_con_paf[0]}#downstream_no_except_neighbor\n')
                                        gene_no_neighbor_end_minus = gene_end_minus + now_from_prev
                                        if gene_no_neighbor_end_minus > int(fields_con_paf[1]) - 1:
                                            gene_no_neighbor_end_minus = int(fields_con_paf[1]) - 1
                                            with open(f'{running_dir}/{d}/no_except_neighbor_list',
                                                      'a') as no_no_neighbor_list:
                                                no_no_neighbor_list.write(
                                                    f'{fields_con_paf[0]}#upstream_no_except_neighbor\n')
                                    if gff_strand == '-':
                                        gene_no_neighbor_start_minus = gene_start_minus - now_from_prev
                                        if gene_no_neighbor_start_minus < 0:
                                            gene_no_neighbor_start_minus = 0
                                            with open(f'{running_dir}/{d}/no_except_neighbor_list',
                                                      'a') as no_no_neighbor_list:
                                                no_no_neighbor_list.write(
                                                    f'{fields_con_paf[0]}#downstream_no_except_neighbor\n')
                                        gene_no_neighbor_end_minus = gene_end_minus + next_from_now
                                        if gene_no_neighbor_end_minus > int(fields_con_paf[1]) - 1:
                                            gene_no_neighbor_end_minus = int(fields_con_paf[1]) - 1
                                            with open(f'{running_dir}/{d}/no_except_neighbor_list',
                                                      'a') as no_no_neighbor_list:
                                                no_no_neighbor_list.write(
                                                    f'{fields_con_paf[0]}#upstream_no_except_neighbor\n')
                                    # print(f'gene_start_minus {gene_start_minus}')
                                    # print(f'gene_end_minus {gene_end_minus}')
                                    # print(f'next_from_now_minus {next_from_now}')
                                    # print(f'now_from_prev_minus {now_from_prev}')
                                    # print(f'gene_no_neighbor_start_minus {gene_no_neighbor_start_minus}')
                                    # print(f'gene_no_neighbor_end_minus {gene_no_neighbor_end_minus}')
                                with open(
                                        f'{running_dir}/{d}/gene_except_neighbor_FASTA/{ch}/phased_gene_except_neighbor_{gene}.fa',
                                        'a') as new_file_gene:
                                    # print(new_seq_list[con_paf_cnt])
                                    # print(new_seq_list[row][gene_region_start:gene_region_end])
                                    new_file_gene.write('>' + header_rename_neighbor + '\n')
                                    if strand == '+':
                                        new_file_gene.write(
                                            edited_seq[gene_no_neighbor_start:gene_no_neighbor_end] + '\n')
                                    if strand == '-':
                                        seq = edited_seq[gene_no_neighbor_start_minus:gene_no_neighbor_end_minus]
                                        seq_c = str(Seq(seq).reverse_complement())
                                        new_file_gene.write(seq_c + '\n')

            elif os.path.exists(f'{running_dir}/{d}/FASTA/{ch}/phased_{gene}.fa') and not os.path.exists(
                    f'{running_dir}/{d}/{d}_mapping_result/{ch}/{d}_mapping_result_{gene}.paf'):
                print(f'sorry! No other output of phased_{gene} could be obtained.'
                      f'This is because the PAF file was not output')
                with open(f'{running_dir}/{d}/only_raw_FASTA_list', 'a') as file:
                    file.write(f'phased_{gene}\n')

            # Move phased_seq gene_region_FASTA sequences with more than {ploidy} sequences to sub.
            if d == 'phased_seq_sub':
                raw_file = f'{running_dir}/phased_seq/FASTA/{ch}/phased_{gene}.fa'
                gene_file = f'{running_dir}/phased_seq/gene_region_FASTA/{ch}/phased_gene_region_{gene}.fa'
                two_file = f'{running_dir}/phased_seq/gene_+2Kb_FASTA/{ch}/phased_gene_+2Kb_{gene}.fa'
                user_file = f'{running_dir}/phased_seq/gene_+user_setting_FASTA/{ch}/phased_gene_+user_setting_{gene}.fa'
                neighbor_file = f'{running_dir}/phased_seq/gene_except_neighbor_FASTA/{ch}/phased_gene_except_neighbor_{gene}.fa'
                map_file = f'{running_dir}/phased_seq/phased_seq_mapping_result/{ch}/phased_seq_mapping_result_{gene}.paf'
                map_file_err = f'{running_dir}/phased_seq/phased_seq_mapping_result/mapping_error_modifiable/phased_seq_mapping_result_{gene}.paf'
                if os.path.exists(gene_file):
                    cnt = 0
                    with open(gene_file, 'r') as f:
                        for line in f:
                            cnt += 1
                    seq_cnt = cnt // 2
                    if seq_cnt > ploidy:
                        move_path = shutil.move(gene_file, f'{running_dir}/phased_seq_sub/gene_region_FASTA/{ch}')
                        if os.path.exists(raw_file):
                            move_path = shutil.move(raw_file, f'{running_dir}/phased_seq_sub/FASTA/{ch}')
                        if os.path.exists(two_file):
                            move_path = shutil.move(two_file, f'{running_dir}/phased_seq_sub/gene_+2Kb_FASTA/{ch}')
                        if os.path.exists(user_file):
                            move_path = shutil.move(user_file, f'{running_dir}/phased_seq_sub/gene_+user_setting_FASTA/{ch}')
                        if os.path.exists(neighbor_file):
                            move_path = shutil.move(neighbor_file, f'{running_dir}/phased_seq_sub/gene_except_neighbor_FASTA/{ch}')
                        if os.path.exists(map_file):
                            move_path = shutil.move(map_file, f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/{ch}')
                        if os.path.exists(map_file_err):
                            move_path = shutil.move(map_file_err, f'{running_dir}/phased_seq_sub/phased_seq_sub_mapping_result/mapping_error_modifiable')

# Divide the output into homo and hetero.
phased_dir_list = ['phased_seq', 'phased_seq_sub']
for d in phased_dir_list:
    os.mkdir(f'{running_dir}/{d}/raw_FASTA')
    for ch in ch_list:
        gene_dir_list = []
        os.mkdir(f'{running_dir}/{d}/raw_FASTA/{ch}')
        os.mkdir(f'{running_dir}/{d}/raw_FASTA/{ch}/hetero')
        os.mkdir(f'{running_dir}/{d}/raw_FASTA/{ch}/homo')
        os.mkdir(f'{running_dir}/{d}/gene_region_FASTA/{ch}/hetero')
        os.mkdir(f'{running_dir}/{d}/gene_region_FASTA/{ch}/homo')
        os.mkdir(f'{running_dir}/{d}/gene_+2Kb_FASTA/{ch}/hetero')
        os.mkdir(f'{running_dir}/{d}/gene_+2Kb_FASTA/{ch}/homo')
        if args.x_cut_out:
            os.mkdir(f'{running_dir}/{d}/gene_+user_setting_FASTA/{ch}/hetero')
            os.mkdir(f'{running_dir}/{d}/gene_+user_setting_FASTA/{ch}/homo')
        if args.no_neighbor:
            os.mkdir(f'{running_dir}/{d}/gene_except_neighbor_FASTA/{ch}/hetero')
            os.mkdir(f'{running_dir}/{d}/gene_except_neighbor_FASTA/{ch}/homo')
        ch_path = f'{running_dir}/mapping_result_split/{ch}'
        for List in os.listdir(ch_path):
            if os.path.isdir(os.path.join(ch_path, List)):
                gene_dir_list.append(List)
        gene_dir_list = sorted(gene_dir_list)
        for gene in gene_dir_list:
            if os.path.exists(f'{running_dir}/{d}/FASTA/{ch}/phased_{gene}.fa'):
                with open(f'{running_dir}/{d}/FASTA/{ch}/phased_{gene}.fa', 'r') as f:
                    lines = f.readlines()
                    line_count = len(lines)
                if line_count >= 3:
                    consensus_hetero = shutil.move(
                        f'{running_dir}/{d}/FASTA/{ch}/phased_{gene}.fa',
                        f'{running_dir}/{d}/raw_FASTA/{ch}/hetero')
                else:
                    consensus_homo = shutil.move(
                        f'{running_dir}/{d}/FASTA/{ch}/phased_{gene}.fa',
                        f'{running_dir}/{d}/raw_FASTA/{ch}/homo')
            if os.path.exists(f'{running_dir}/{d}/gene_region_FASTA/{ch}/phased_gene_region_{gene}.fa'):
                with open(f'{running_dir}/{d}/gene_region_FASTA/{ch}/phased_gene_region_{gene}.fa', 'r') as f:
                    lines = f.readlines()
                    line_count = len(lines)
                if line_count >= 3:
                    consensus_hetero = shutil.move(
                        f'{running_dir}/{d}/gene_region_FASTA/{ch}/phased_gene_region_{gene}.fa',
                        f'{running_dir}/{d}/gene_region_FASTA/{ch}/hetero')
                else:
                    consensus_homo = shutil.move(
                        f'{running_dir}/{d}/gene_region_FASTA/{ch}/phased_gene_region_{gene}.fa',
                        f'{running_dir}/{d}/gene_region_FASTA/{ch}/homo')
            if os.path.exists(f'{running_dir}/{d}/gene_+2Kb_FASTA/{ch}/phased_gene_+2Kb_{gene}.fa'):
                with open(f'{running_dir}/{d}/gene_+2Kb_FASTA/{ch}/phased_gene_+2Kb_{gene}.fa', 'r') as f:
                    lines = f.readlines()
                    line_count = len(lines)
                if line_count >= 3:
                    consensus_hetero = shutil.move(
                        f'{running_dir}/{d}/gene_+2Kb_FASTA/{ch}/phased_gene_+2Kb_{gene}.fa',
                        f'{running_dir}/{d}/gene_+2Kb_FASTA/{ch}/hetero')
                else:
                    consensus_homo = shutil.move(
                        f'{running_dir}/{d}/gene_+2Kb_FASTA/{ch}/phased_gene_+2Kb_{gene}.fa',
                        f'{running_dir}/{d}/gene_+2Kb_FASTA/{ch}/homo')
            if os.path.exists(
                    f'{running_dir}/{d}/gene_+user_setting_FASTA/{ch}/phased_gene_+user_setting_{gene}.fa'):
                with open(f'{running_dir}/{d}/gene_+user_setting_FASTA/{ch}/phased_gene_+user_setting_{gene}.fa',
                          'r') as f:
                    lines = f.readlines()
                    line_count = len(lines)
                if line_count >= 3:
                    consensus_hetero = shutil.move(
                        f'{running_dir}/{d}/gene_+user_setting_FASTA/{ch}/phased_gene_+user_setting_{gene}.fa',
                        f'{running_dir}/{d}/gene_+user_setting_FASTA/{ch}/hetero')
                else:
                    consensus_homo = shutil.move(
                        f'{running_dir}/{d}/gene_+user_setting_FASTA/{ch}/phased_gene_+user_setting_{gene}.fa',
                        f'{running_dir}/{d}/gene_+user_setting_FASTA/{ch}/homo')
            if os.path.exists(
                    f'{running_dir}/{d}/gene_except_neighbor_FASTA/{ch}/phased_gene_except_neighbor_{gene}.fa'):
                with open(
                        f'{running_dir}/{d}/gene_except_neighbor_FASTA/{ch}/phased_gene_except_neighbor_{gene}.fa',
                        'r') as f:
                    lines = f.readlines()
                    line_count = len(lines)
                if line_count >= 3:
                    consensus_hetero = shutil.move(
                        f'{running_dir}/{d}/gene_except_neighbor_FASTA/{ch}/phased_gene_except_neighbor_{gene}.fa',
                        f'{running_dir}/{d}/gene_except_neighbor_FASTA/{ch}/hetero')
                else:
                    consensus_homo = shutil.move(
                        f'{running_dir}/{d}/gene_except_neighbor_FASTA/{ch}/phased_gene_except_neighbor_{gene}.fa',
                        f'{running_dir}/{d}/gene_except_neighbor_FASTA/{ch}/homo')
    shutil.rmtree(f'{running_dir}/{d}/FASTA')


# perform alignment
if args.align:
    print('Align genes that determined to be hetero.\n')
    align_dir_list = ['alignment', 'alignment_sub']
    phased_dir_list = ['phased_seq', 'phased_seq_sub']
    num = [0, 1]
    for i in num:
        da = align_dir_list[i]
        dp = phased_dir_list[i]
        os.mkdir(f'{running_dir}/{da}')
        tool, align_2kb, x_option_align,= args.align
        if tool == 'Muscle':
            Command, sub = args.muscle_setting
            for ch in ch_list:
                gene_dir_list = []
                ch_path = f'{running_dir}/mapping_result_split/{ch}'
                for List in os.listdir(ch_path):
                    if os.path.isdir(os.path.join(ch_path, List)):
                        gene_dir_list.append(List)
                gene_dir_list = sorted(gene_dir_list)
                if align_2kb == 'Y':
                    if not os.path.exists(f'{running_dir}/{da}/gene_+2Kb'):
                        os.mkdir(f'{running_dir}/{da}/gene_+2Kb')
                    os.mkdir(f'{running_dir}/{da}/gene_+2Kb/{ch}')
                if x_option_align == 'xY':
                    if not os.path.exists(f'{running_dir}/{da}/gene_+user_setting'):
                        os.mkdir(f'{running_dir}/{da}/gene_+user_setting')
                    os.mkdir(f'{running_dir}/{da}/gene_+user_setting/{ch}')
                for gene in gene_dir_list:
                    if align_2kb == 'Y':
                        if os.path.exists(
                            f'{running_dir}/{dp}/gene_+2Kb_FASTA/{ch}/hetero/phased_gene_+2Kb_{gene}.fa'):
                            cmd_muscle_1 = \
                                f'muscle -{Command} {running_dir}/{dp}/gene_+2Kb_FASTA/{ch}/hetero/phased_gene_+2Kb_{gene}.fa -output {running_dir}/{da}/gene_+2Kb/{ch}/align_gene_+2Kb_{gene}.aln -threads {threads} '
                            subprocess.run(cmd_muscle_1, shell=True)
                            if not os.path.exists(f'{running_dir}/{da}/gene_+2Kb/{ch}/align_gene_+2Kb_{gene}.aln') or os.path.getsize(f'{running_dir}/{da}/gene_+2Kb/{ch}/align_gene_+2Kb_{gene}.aln') == 0:
                                if os.path.exists(f'{running_dir}/{da}/gene_+2Kb/{ch}/align_gene_+2Kb_{gene}.aln'):
                                    os.remove(f'{running_dir}/{da}/gene_+2Kb/{ch}/align_gene_+2Kb_{gene}.aln')
                                if sub == 'PRANK':
                                    print('\nAlleleMiner rescue : Align with PRANK')
                                    cmd_prank_1 = \
                                        f'prank -d={running_dir}/{dp}/gene_+2Kb_FASTA/{ch}/hetero/phased_gene_+2Kb_{gene}.fa -o={running_dir}/{da}/gene_+2Kb/{ch}/align_gene_+2Kb_{gene}'
                                    subprocess.run(cmd_prank_1, shell=True)
                                    if os.path.exists(
                                            f'{running_dir}/{da}/gene_+2Kb/{ch}/align_gene_+2Kb_{gene}.best.fas'):
                                        os.rename(f'{running_dir}/{da}/gene_+2Kb/{ch}/align_gene_+2Kb_{gene}.best.fas',
                                                  f'{running_dir}/{da}/gene_+2Kb//{ch}align_gene_+2Kb_{gene}.prank.aln')
                                elif sub == 'Clustal':
                                    print('\nAlleleMiner rescue : Align with ClustalΩ')
                                    cmd_clustal_1 = \
                                        f'clustalo -i {running_dir}/{dp}/gene_+2Kb_FASTA/{ch}/hetero/phased_gene_+2Kb_{gene}.fa -t DNA -o {running_dir}/{da}/gene_+2Kb/{ch}/align_gene_+2Kb_{gene}.clustal.aln'
                                    print(f'Align {gene} by clustalΩ')
                                    subprocess.run(cmd_clustal_1, shell=True)
                                elif sub == 'NONE':
                                    pass
                    if x_option_align == 'xY':
                        if os.path.exists(
                            f'{running_dir}/{dp}/gene_+user_setting_FASTA/{ch}/hetero/phased_gene_+user_setting_{gene}.fa'):
                            cmd_muscle_2 = \
                                f'muscle -{Command} {running_dir}/{dp}/gene_+user_setting_FASTA/{ch}/hetero/phased_gene_+user_setting_{gene}.fa -output {running_dir}/{da}/gene_+user_setting/{ch}/align_gene_+user_setting_{gene}.aln -threads {threads}'
                            subprocess.run(cmd_muscle_2, shell=True)
                            if not os.path.exists(f'{running_dir}/{da}/gene_+user_setting/{ch}/align_gene_+user_setting_{gene}.aln') or os.path.getsize(f'{running_dir}/{da}/gene_+user_setting/{ch}/align_gene_+user_setting_{gene}.aln') == 0:
                                if os.path.exists(f'{running_dir}/{da}/gene_+user_setting/{ch}/align_gene_+user_setting_{gene}.aln'):
                                    os.remove(f'{running_dir}/{da}/gene_+user_setting/{ch}/align_gene_+user_setting_{gene}.aln')
                                if sub == 'PRANK':
                                    print('\nAlleleMiner rescue : Align with PRANK')
                                    cmd_prank_2 = \
                                        f'prank -d={running_dir}/{dp}/gene_+user_setting_FASTA/{ch}/hetero/phased_gene_+user_setting_{gene}.fa -o={running_dir}/{da}/gene_+user_setting/{ch}/align_gene_+user_setting_{gene}'
                                    subprocess.run(cmd_prank_2, shell=True)
                                    if os.path.exists(
                                            f'{running_dir}/{da}/gene_+user_setting/{ch}/align_gene_+user_setting_{gene}.best.fas'):
                                        os.rename(f'{running_dir}/{da}/gene_+user_setting/{ch}/align_gene_+user_setting_{gene}.best.fas',
                                                  f'{running_dir}/{da}/gene_+user_setting/{ch}/align_gene_+user_setting_{gene}.prank.aln')
                                elif sub == 'Clustal':
                                    print('\nAlleleMiner rescue : Align with ClustalΩ')
                                    cmd_clustal_2 = \
                                        f'clustalo -i {running_dir}/{dp}/gene_+user_setting_FASTA/{ch}/hetero/phased_gene_+user_setting_{gene}.fa -t DNA -o {running_dir}/{da}/gene_+user_setting/{ch}/align_gene_+user_setting_{gene}.clustal.aln'
                                    print(f'Align {gene} by clustalΩ')
                                    subprocess.run(cmd_clustal_2, shell=True)
                                elif sub == 'NONE':
                                    pass

        elif tool == 'PRANK':
            for ch in ch_list:
                gene_dir_list = []
                ch_path = f'{running_dir}/mapping_result_split/{ch}'
                for List in os.listdir(ch_path):
                    if os.path.isdir(os.path.join(ch_path, List)):
                        gene_dir_list.append(List)
                gene_dir_list = sorted(gene_dir_list)
                if align_2kb == 'Y':
                    if not os.path.exists(f'{running_dir}/{da}/gene_+2Kb'):
                        os.mkdir(f'{running_dir}/{da}/gene_+2Kb')
                    os.mkdir(f'{running_dir}/{da}/gene_+2Kb/{ch}')
                if x_option_align == 'xY':
                    if not os.path.exists(f'{running_dir}/{da}/gene_+user_setting'):
                        os.mkdir(f'{running_dir}/{da}/gene_+user_setting')
                    os.mkdir(f'{running_dir}/{da}/gene_+user_setting/{ch}')
                for gene in gene_dir_list:
                    if align_2kb == 'Y':
                        if os.path.exists(
                                f'{running_dir}/{dp}/gene_+2Kb_FASTA/{ch}/hetero/phased_gene_+2Kb_{gene}.fa'):
                            cmd_prank_a = \
                                f'prank -d={running_dir}/{dp}/gene_+2Kb_FASTA/{ch}/hetero/phased_gene_+2Kb_{gene}.fa -o={running_dir}/{da}/gene_+2Kb/{ch}/align_gene_+2Kb_{gene}'
                            subprocess.run(cmd_prank_a, shell=True)
                            if os.path.exists(
                                    f'{running_dir}/{da}/gene_+2Kb/{ch}/align_gene_+2Kb_{gene}.best.fas'):
                                os.rename(f'{running_dir}/{da}/gene_+2Kb/{ch}/align_gene_+2Kb_{gene}.best.fas',
                                          f'{running_dir}/{da}/gene_+2Kb/{ch}/align_gene_+2Kb_{gene}.prank.aln')
                    if x_option_align == 'xY':
                        if os.path.exists(
                                f'{running_dir}/{dp}/gene_+2Kb_FASTA/{ch}/hetero/phased_gene_+2Kb_{gene}.fa'):
                            cmd_prank_b = \
                                f'prank -d={running_dir}/{dp}/gene_+user_setting_FASTA/{ch}/hetero/phased_gene_+user_setting_{gene}.fa -o={running_dir}/{da}/gene_+user_setting/{ch}/align_gene_+user_setting_{gene}'
                            subprocess.run(cmd_prank_b, shell=True)
                            if os.path.exists(
                                    f'{running_dir}/{da}/gene_+user_setting/{ch}/align_gene_+user_setting_{gene}.best.fas'):
                                os.rename(f'{running_dir}/{da}/gene_+user_setting/{ch}/align_gene_+user_setting_{gene}.best.fas',
                                          f'{running_dir}/{da}/gene_+user_setting/{ch}/align_gene_+user_setting_{gene}.prank.aln')

        elif tool == 'Clustal':
            for ch in ch_list:
                gene_dir_list = []
                ch_path = f'{running_dir}/mapping_result_split/{ch}'
                for List in os.listdir(ch_path):
                    if os.path.isdir(os.path.join(ch_path, List)):
                        gene_dir_list.append(List)
                gene_dir_list = sorted(gene_dir_list)
                if align_2kb == 'Y':
                    if not os.path.exists(f'{running_dir}/{da}/gene_+2Kb'):
                        os.mkdir(f'{running_dir}/{da}/gene_+2Kb')
                    os.mkdir(f'{running_dir}/{da}/gene_+2Kb/{ch}')
                if x_option_align == 'xY':
                    if not os.path.exists(f'{running_dir}/{da}/gene_+user_setting'):
                        os.mkdir(f'{running_dir}/{da}/gene_+user_setting')
                    os.mkdir(f'{running_dir}/{da}/gene_+user_setting/{ch}')
                for gene in gene_dir_list:
                    if align_2kb == 'Y':
                        if os.path.exists(
                            f'{running_dir}/{dp}/gene_+2Kb_FASTA/{ch}/hetero/phased_gene_+2Kb_{gene}.fa'):
                            cmd_clustal_a = \
                                f'clustalo -i {running_dir}/{dp}/gene_+2Kb_FASTA/{ch}/hetero/phased_gene_+2Kb_{gene}.fa -t DNA -o {running_dir}/{da}/gene_+2Kb/{ch}/align_gene_+2Kb_{gene}.clustal.aln'
                            print(f'Align {gene} by clustalΩ')
                            subprocess.run(cmd_clustal_a, shell=True)
                    if x_option_align == 'xY':
                        if os.path.exists(
                            f'{running_dir}/{dp}/gene_+user_setting_FASTA/{ch}/hetero/phased_gene_+user_setting_{gene}.fa'):
                            cmd_clustal_b = \
                                f'clustalo -i {running_dir}/{dp}/gene_+user_setting_FASTA/{ch}/hetero/phased_gene_+user_setting_{gene}.fa -t DNA -o {running_dir}/{da}/gene_+user_setting/{ch}/align_gene_+user_setting_{gene}.clustal.aln'
                            print(f'Align {gene} by clustalΩ')
                            subprocess.run(cmd_clustal_b, shell=True)


# Divide outputs into final and intermediate
os.mkdir(f'{running_dir}/final')
os.mkdir(f'{running_dir}/intermediate')

if os.path.exists(f'{running_dir}/alignment'):
    final_dir1 = shutil.move(f'{running_dir}/alignment', f'{running_dir}/final')
if os.path.exists(f'{running_dir}/alignment_sub'):
    final_dir2 = shutil.move(f'{running_dir}/alignment_sub', f'{running_dir}/final')
final_dir3 = shutil.move(f'{running_dir}/phased_seq', f'{running_dir}/final')
final_dir4 = shutil.move(f'{running_dir}/phased_seq_sub', f'{running_dir}/final')
intermediate_dir1 = shutil.move(f'{running_dir}/all_output_gff', f'{running_dir}/intermediate')
intermediate_dir2 = shutil.move(f'{running_dir}/assembly_result', f'{running_dir}/intermediate')
if os.path.exists(f'{running_dir}/busco_downloads'):
    intermediate_dir3 = shutil.move(f'{running_dir}/busco_downloads', f'{running_dir}/intermediate')
if os.path.exists(f'{running_dir}/busco_pro_fa'):
    intermediate_dir4 = shutil.move(f'{running_dir}/busco_pro_fa', f'{running_dir}/intermediate')
if os.path.exists(f'{running_dir}/busco_result'):
    intermediate_dir5 = shutil.move(f'{running_dir}/busco_result', f'{running_dir}/intermediate')
intermediate_dir6 = shutil.move(f'{running_dir}/mapping_high_q', f'{running_dir}/intermediate')
intermediate_dir7 = shutil.move(f'{running_dir}/mapping_result', f'{running_dir}/intermediate')
intermediate_dir8 = shutil.move(f'{running_dir}/mapping_result_split', f'{running_dir}/intermediate')
intermediate_dir9 = shutil.move(f'{running_dir}/mapping_straddled', f'{running_dir}/intermediate')
intermediate_dir10 = shutil.move(f'{running_dir}/read_extract', f'{running_dir}/intermediate')
if os.path.exists(f'{running_dir}/ref_single_copy_gene_multi_fa'):
    intermediate_dir11 = shutil.move(f'{running_dir}/ref_single_copy_gene_multi_fa',
                                     f'{running_dir}/intermediate')
if os.path.exists(f'{running_dir}/ref_gene_multi_fa'):
    intermediate_dir12 = shutil.move(f'{running_dir}/ref_gene_multi_fa', f'{running_dir}/intermediate')

# Delete empty directories
for root, dirs, files in os.walk(running_dir, topdown=False):
        for d in dirs:
            full_path = os.path.join(root, d)
            try:
                # 空なら削除
                if not os.listdir(full_path):
                    os.rmdir(full_path)
                    print(f"delete: {full_path}")
            except Exception as e:
                print(f"delete failed: {full_path} ({e})")

# Make list
os.mkdir(f'{running_dir}/final/list')

# Create csv and FASTA files
print("make list")
analysis_lists = ["gene_region", "gene_+2Kb", "gene_+user_setting", "except_neighbor"]
for name in analysis_lists:
    if os.path.exists(f"{running_dir}/final/phased_seq/{name}_FASTA") or os.path.exists(f"{running_dir}/final/phased_seq_sub/{name}_FASTA"):
        with open(f'{running_dir}/final/list/{name}_summary.csv', 'w', newline='') as gene_list_csv, open(f'{running_dir}/final/list/{name}.fa', 'a') as gene_list_fa:
            writer = csv.writer(gene_list_csv)

        # Creating a header
            writer.writerow(["seq_name", "zygosity", "length", "MinHash_info", "hash_string", "jaccard_similarity"])

        # Preparing to write to csv
            for ch in ch_list:
                gene_dir_list = []
                ch_path = f'{running_dir}/intermediate/mapping_result_split/{ch}'
                for List in os.listdir(ch_path):
                    if os.path.isdir(os.path.join(ch_path, List)):
                        gene_dir_list.append(List)
                gene_dir_list = sorted(gene_dir_list)
                for gene in gene_dir_list:
                    pattern = f'{running_dir}/final/*/{name}_FASTA/{ch}/*/phased_{name}_{gene}.fa'
                    matched = glob.glob(pattern)
                    if matched:
                        # print(matched[0])
                        with open(matched[0], 'r') as phased_seq:
                            lines = phased_seq.readlines()
                            for line in lines:
                                gene_list_fa.write(line)
                            line_count = len(lines)
                            seq_count = len(lines) // 2
                            if seq_count <= 1:
                                zyogosity = 'homo'
                            else:
                                zyogosity = 'hetero'

                            seq_names = {}
                            min_dict = {}
                            hash_string_dict = {}
                            length_dict = {}
                            for i in range(line_count):
                                n = i // 2
                                if i % 2 == 0:
                                    seq_names[f"seq_name{n}"] = lines[i].lstrip('>')
                                else:
                                    seq = lines[i].strip()
                                    # k-merの抽出
                                    tokens = set(seq[ii:ii + K] for ii in range(len(seq) - K + 1))
                                    # 配列のスケッチの作成
                                    m = MinHash(num_perm=num_perm)
                                    for d in tokens:
                                        m.update(d.encode('utf8'))
                                    # MinHashスケッチのハッシュ値をビット列/文字列として表現
                                    # m1.hashvalues は num_perm個のuint64配列
                                    # これを文字列化するには、各要素を16進数などに変換して連結する
                                    hash_string_m = "".join(format(h, '016x') for h in m.hashvalues)

                                    min_dict[f"MinHash_info{n}"] = m
                                    hash_string_dict[f"hash_string{n}"] = hash_string_m
                                    length_dict[f"length{n}"] = str(len(seq))

                            jaccard_dict = {}
                            if seq_count <= 1:
                                jaccard_dict[f"jaccard_similarity0"] = "-"
                            else:
                                for i in range(seq_count):
                                    jaccard_list = []
                                    for j in range(seq_count):
                                        if i != j:
                                            jaccard = min_dict[f"MinHash_info{i}"].jaccard(min_dict[f"MinHash_info{j}"])
                                            jaccard_list.append(f"vs {seq_names[f'seq_name{j}']} = {jaccard}" )
                                            jaccard_dict[f"jaccard_similarity{i}"] = ", ".join(jaccard_list)
                            # Write to CSV
                            for i in range(seq_count):
                                writer.writerow([seq_names[f"seq_name{i}"],
                                                 zyogosity, length_dict[f"length{i}"],
                                                 min_dict[f"MinHash_info{i}"],
                                                 hash_string_dict[f"hash_string{i}"],
                                                 jaccard_dict[f"jaccard_similarity{i}"]])

# Count the number of gene_region_FASTA files
HE_T = 0
HO_T = 0
HE_S = 0
HO_S = 0
for ch in ch_list:
    if os.path.exists(f'{running_dir}/final/phased_seq/gene_region_FASTA/{ch}/hetero'):
        HE_T_count = len(os.listdir(f'{running_dir}/final/phased_seq/gene_region_FASTA/{ch}/hetero'))
        HE_T += HE_T_count
    if os.path.exists(f'{running_dir}/final/phased_seq/gene_region_FASTA/{ch}/homo'):
        HO_T_count = len(os.listdir(f'{running_dir}/final/phased_seq/gene_region_FASTA/{ch}/homo'))
        HO_T += HO_T_count
    if os.path.exists(f'{running_dir}/final/phased_seq_sub/gene_region_FASTA/{ch}/hetero'):
        HE_S_count = len(os.listdir(f'{running_dir}/final/phased_seq_sub/gene_region_FASTA/{ch}/hetero'))
        HE_S += HE_S_count
    if os.path.exists(f'{running_dir}/final/phased_seq_sub/gene_region_FASTA/{ch}/homo'):
        HO_S_count = len(os.listdir(f'{running_dir}/final/phased_seq_sub/gene_region_FASTA/{ch}/homo'))
        HO_S += HO_S_count

HE = HE_T + HE_S
HO = HO_T + HO_S
SUM = HE + HO

with open(f'{running_dir}/final/list/file_count_gene_region_FASTA.txt', 'w') as f:
    f.write(f"gene_region_FASTA outputs\n")
    f.write(f"Sum = {SUM}\n")
    f.write(f"hetero = {HE}\n")
    f.write(f"homo = {HO}\n")
    f.write(f"phased_seq hetero = {HE_T}\n")
    f.write(f"phased_seq homo = {HO_T}\n")
    f.write(f"phased_seq_sub hetero = {HE_S}\n")
    f.write(f"phased_seq_sub homo = {HO_S}\n")

print('AlleleMiner end')
