# AlleleMiner

A long-read pipeline for **gene-wise *de novo* allele phasing** in diploid genomes.

AlleleMiner reconstructs the two allele sequences of each target gene directly from PacBio HiFi
reads. The reference genome is used only to define the target locus and to recruit reads — the
sequence itself is rebuilt by *de novo* assembly at each locus, which keeps reference dependence
to a minimum.

📄 **Paper:** Kiryu *et al.* (2026) *DNA Research* **33**(2): dsag004 —
[https://doi.org/10.1093/dnares/dsag004](https://academic.oup.com/dnaresearch/article/33/2/dsag004/8503928)
(Open Access)

---

## Overview

AlleleMiner takes **three input files** and runs **six steps**:

| # | Step | Tool |
|---|---|---|
| 1 | Extract target gene region sequences (**TGRS**) from the reference genome using the GFF3 | GFF3toolkit |
| 2 | Map long reads to each TGRS | Minimap2 |
| 3 | Select high-quality reads that span the whole TGRS | — |
| 4 | Gene-wise *de novo* assembly of the selected reads | Hifiasm, Flye |
| 5 | Re-map phased contigs to the TGRS and cut out the allele sequences | Minimap2 |
| 6 | Summarise zygosity, length, hash value and Jaccard similarity | Datasketch |

A heterozygous locus yields **two** allele sequences; a homozygous locus yields **one**.
Flanking regions are extracted together with the gene body.

---

## Requirements

Linux / UNIX. Tested on Ubuntu 24.04.2 LTS.

### Core — required

| Tool | Version tested |
|---|---|
| Python | 3.10 or later |
| SeqKit | 2.2.0 |
| GFF3toolkit | 2.1.0 |
| Minimap2 | 2.29-r1283 |
| Hifiasm | 0.25.0-r726 |
| Flye | 2.9.5-b1801 |
| Biopython | 1.85 |
| Datasketch | 1.6.5 |

### Optional — not needed for a default run

Install these **only if you use the corresponding option**.

| Tool | Version tested | Needed for |
|---|---|---|
| BUSCO | 6.0.0 | `-b` — restrict the analysis to single-copy orthologues |
| HMMER | 3.4 | required by BUSCO |
| MUSCLE | 5.3 | `-a` — multiple alignment |
| Clustal Omega | 1.2.4 | `-a` — multiple alignment |
| PRANK | 170427 | `-a` — multiple alignment |

---

## Installation

The commands below assume [**Miniforge**](https://github.com/conda-forge/miniforge)
(conda with `conda-forge` as the default channel). They also work with Anaconda or Miniconda.

```bash
git clone https://github.com/YukinariKiryu/AlleleMiner.git
cd AlleleMiner
```

### Option A — use the bundled environment file (recommended)

```bash
conda env create -f yaml/AMpy311-complete.yml
conda activate AMpy311
pip install -r yaml/pip_requirements.txt
```

### Option B — build the environment yourself (for Mac User)

```bash
# 1. channels
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# 2. environment
conda create -n alleleminer python=3.11 \
    seqkit=2.2.0 \
    minimap2=2.29 \
    hifiasm=0.25.0 \
    flye=2.9.5 \
    biopython=1.85 \
    datasketch=1.6.5

conda activate alleleminer

# 3. GFF3toolkit — install last, with pip
pip install gff3tool==2.1.0
```

### Make `alleleminer` available on your PATH

```bash
mkdir -p ~/.local/bin
cp AlleleMiner1.0-r2.py ~/.local/bin/

cat > ~/.local/bin/alleleminer << 'EOF'
#!/bin/bash
python3 ~/.local/bin/AlleleMiner1.0-r2.py "$@"
EOF
chmod +x ~/.local/bin/alleleminer

export PATH="$HOME/.local/bin:$PATH"   # add this line to ~/.bashrc to make it permanent
```

Check the installation:

```bash
alleleminer -h
```

### Optional extras

BUSCO pulls in a large dependency tree, so install it in a **separate environment**.

```bash
# for the -b option
conda create -n alleleminer-busco busco=6.0.0 hmmer=3.4

# for the -a option (multiple alignment)
conda install -n alleleminer muscle=5.3 clustalo=1.2.4 prank=170427
```

---

## Test data

Two simulated datasets are provided in `tests/`, corresponding to the **high-variation** and
**low-variation** models used for validation in the paper. Use them to confirm that your
installation works before moving to your own data.

```
tests/
├── high_heterozigosity/
│   ├── query/     reads_h.fastq.tar.gz
│   └── reference/ genome_high_alleleA.zip
└── low_heterozigosity/
    ├── query/     reads_l.fastq.tar.gz
    └── reference/ genome_low_alleleA.zip
```

Each dataset contains a reference genome FASTA, its GFF3 annotation, and simulated HiFi reads.

Unpack and run:

```bash
conda activate AMpy311
cd tests/high_heterozigosity

tar -xzf query/reads_h.fastq.tar.gz -C query/
unzip    reference/genome_high_alleleA.zip -d reference/

alleleminer \
  --gff   reference/genome_high_alleleA/genome_high_alleleA.gff3 \
  --base  reference/genome_high_alleleA/genome_high_alleleA.fasta \
  --read  query/reads_h.fastq \
  -P 2 \
  -N test_high
```

`-P 2` specifies a diploid sample. The run takes **roughly fifteen minutes**.

If `test_high/final/list/gene_region.fa` contains allele sequences, your installation is working.
See [Output](#output) for what the other files contain.

The low-variation dataset is run the same way:

```bash
cd tests/low_heterozigosity

tar -xzf query/reads_l.fastq.tar.gz -C query/
unzip    reference/genome_low_alleleA.zip -d reference/

alleleminer \
  --gff   reference/genome_low_alleleA/genome_low_alleleA.gff3 \
  --base  reference/genome_low_alleleA/genome_low_alleleA.fasta \
  --read  query/reads_l.fastq \
  -P 2 \
  -N test_low
```

> These datasets are for checking the installation and for illustrating the effect of variation
> density. They are not intended to reproduce the full results in the paper.

---

## Quick start

```bash
conda activate AMpy311

alleleminer \
  --gff   reference.primaryTranscript.gff3 \
  --base  reference.genome.fasta \
  --read  sample.hifi.fastq.gz \
  -P 2 \
  -N my_sample
```

`--gff`, `--base`, `--read` and `-P` are required.
Results are written to a directory named after `-N` (default: `organism`).

> **Coverage.** About **30× HiFi** is recommended. Below that, one of the two alleles starts to
> drop out — and a dropped allele looks exactly like a homozygous locus.

---

## Input

| File | Option | Description |
|---|---|---|
| Reference genome FASTA | `--base` | Whole-genome sequence |
| GFF3 annotation | `--gff` | Positions of the genes you want to analyse |
| Long reads | `--read` | PacBio HiFi reads, FASTQ |

Make sure the sequence names in the FASTA and the GFF3 match exactly
(`chr1` and `Chr1` are different).

---

## Output

Two directories are created. **Only `final/` matters for normal use** —
`intermediate/` holds working files and is only needed for troubleshooting.

```
<sample>/
├── intermediate/        working files — normally ignore
└── final/               ← everything you need is here
    ├── list/
    │   ├── gene_region.fa                  all allele sequences, gene body
    │   ├── gene_region_summary.csv         summary for the above
    │   ├── gene_+2Kb.fa                    gene body ± 2 kb
    │   ├── gene_+2Kb_summary.csv           summary for the above
    │   ├── gene_+user_setting.fa           gene body ± user-defined region  (with -x)
    │   ├── gene_+user_setting_summary.csv  summary for the above            (with -x)
    │   └── file_count_gene_region_FASTA.txt   number of homo / hetero loci
    │
    ├── phased_seq/                         the same sequences, split per gene
    │   ├── gene_region_FASTA/<chr>/{hetero,homo}/phased_gene_region_<gene>.fa
    │   ├── gene_+2Kb_FASTA/<chr>/{hetero,homo}/phased_gene_+2Kb_<gene>.fa
    │   ├── gene_+user_setting_FASTA/ ...                                     (with -x)
    │   ├── raw_FASTA/<chr>/{hetero,homo}/phased_<gene>.fa    phased contigs, uncut
    │   └── phased_seq_mapping_result/      contig-to-TGRS mapping (PAF)
    │
    ├── phased_seq_sub/                     loci that produced more contigs than the ploidy
    └── alignment/                          multiple alignments               (with -a)
```

### Where to start

| You want | Look at |
|---|---|
| All allele sequences in one file | `final/list/gene_region.fa` |
| Allele sequences with promoter regions | `final/list/gene_+2Kb.fa` |
| One gene at a time | `final/phased_seq/gene_region_FASTA/` |
| Zygosity, length, and cross-sample comparison | `final/list/*_summary.csv` |

The summary CSV reports **zygosity, sequence length, hash value and Jaccard similarity** for
every allele. Identical hash values mean identical (or near-identical) sequences, so you can find
shared and sample-specific alleles across many samples **without running a full-length alignment**.

Loci in `phased_seq_sub/` produced more contigs than the specified ploidy. In a gene family this
often indicates paralogue co-amplification or copy-number variation, so treat them with care.

---

## Options

| Option | Description |
|---|---|
| `-h, --help` | Show help |
| `--gff` | GFF3 annotation of the reference **(required)** |
| `--base` | Reference whole-genome FASTA **(required)** |
| `--read` | Query reads, FASTQ **(required)** |
| `-P, --ploidy` | Ploidy of the sample; `2` for diploid **(required)** |
| `-N, --name` | Sample name; used for the output directory (default: `organism`) |
| `-t, --threads` | CPU threads (default: 2) |
| `--pro` | Amino-acid FASTA of the reference organism |
| `-b, --busco` | Restrict the analysis to BUSCO single-copy orthologues (e.g. `-b embryophyta_odb10`) |
| `-n, --no_neighbor` | Exclude neighbouring genes from the extracted sequence. Requires a GFF3 containing all genes; can be combined with `-b` |
| `-x, --x_cut_out` | Cut out user-defined flanking regions. Two arguments: upstream, downstream in bp (e.g. `-x 200 100`) |
| `--minimap_Q` | Minimum mapping quality for Minimap2 (default: 60; above 60 may return no hits) |
| `--map_q1` | Quality threshold for read-to-TGRS mapping, 0–1 (default: 0.95) |
| `--map_q2` | Quality threshold when correcting the contig-to-TGRS PAF, 0–1 (default: 0.95) |
| `--k_hifiasm` | k-mer size for Hifiasm (default: 51) |
| `-f, --only_flye` | Use Flye only for assembly |
| `--hap_solved` | Enable Flye's `--keep-haplotypes` |
| `--flye_input` | Read type passed to Flye, without the leading `--` (e.g. `--flye_input pacbio-hifi`) |
| `-a, --align` | Run multiple alignment. Three arguments: tool (`Muscle`/`PRANK`/`Clustal`), align gene ± 2 kb (`Y`/`N`), align gene ± user-defined region (`xY`/`xN`). Example: `-a Muscle N xY` |
| `--muscle_setting` | Required when `-a Muscle` is used. Two arguments: mode (`align`/`super5`), fallback tool (`PRANK`/`Clustal`). Example: `--muscle_setting align PRANK` |
| `--num_perm` | Number of hash functions for MinHash. Higher is more accurate but slower (128 / 256 / 512; default: 128) |
| `-K, --Kmer` | k-mer size for hashing the phased sequences (default: 15) |

---

## Repository layout

```
AlleleMiner/
├── AlleleMiner1.0-r2.py     the pipeline
├── tests/                   simulated test datasets
├── yaml/                    conda environment files and pip requirements
├── LICENSE
└── README.md
```

---

## Citation

If you use AlleleMiner, please cite:

> Kiryu, Y., Kawahara, Y., Endo, T., Horiike, T., Shirasawa, K., Isobe, S., Shimada, T. and
> Fujii, H. (2026) AlleleMiner: a long-read pipeline for gene-wise *de novo* allele phasing and
> variant detection in diploid citrus cultivars. *DNA Research* **33**(2): dsag004.
> https://doi.org/10.1093/dnares/dsag004

---

## Related resources

- **Mikan Genome Database (MiGD)** — citrus reference genomes and annotations
  <https://mikan.dna.naro.go.jp/migd2/>

## License

AlleleMiner is released under the MIT License. See [LICENSE](LICENSE).

## Contact

Yukinari Kiryu, Hiroshi Fujii (fujii.hiroshi@shizuoka.ac.jp) — Shizuoka University

## Funding

This work was supported by Japan Society for the Promotion of Science (JSPS) KAKENHI
Grant Numbers JP 23K26902 and JP 23K05224.
