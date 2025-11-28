# BIFS_Systems - RNA-Seq Analysis Pipeline

A comprehensive RNA-seq analysis pipeline that performs quality control, read trimming, expression quantification, and visualization.

## Features

- **Sequence Download**: Fetches sequences from GenBank/NCBI
- **FastQC-like Analysis**: Quality assessment of raw FASTQ files
  - Per-base quality scores
  - GC content distribution
  - Adapter contamination detection
- **Read Trimming**: Adapter removal and quality filtering (fastp-like)
  - Adapter sequence trimming
  - Low-quality base removal
  - Short read filtering
- **Expression Quantification**: Gene expression analysis (Salmon-like)
  - Raw read counts
  - TPM (Transcripts Per Million) normalization
  - CPM (Counts Per Million) normalization
- **Top Gene Identification**: Identifies top 10-20 expressed genes
- **Visualization**: Generates multiple publication-ready figures
  - Per-base quality plots
  - GC content distribution
  - Trimming statistics
  - Expression heatmaps
  - Top genes bar plots
  - Sample correlation heatmaps
- **Enrichment Analysis**: GO term enrichment analysis

## Installation

```bash
# Clone the repository
git clone https://github.com/ssamsuss/BIFS_Systems.git
cd BIFS_Systems

# Install dependencies
pip install -r requirements.txt
```

## Usage

```bash
# Run the complete pipeline
python src/rnaseq_pipeline.py
```

## Output Files

The pipeline generates the following files in the `results/` directory:

### Reports
- `fastqc_summary.txt` - FastQC-like quality assessment report
- `trimming_report.txt` - Read trimming statistics

### Count Tables
- `raw_counts.csv` - Raw read counts per gene
- `tpm_counts.csv` - TPM normalized counts
- `cpm_counts.csv` - CPM normalized counts
- `top_expressed_genes.csv` - Top 15 expressed genes

### Enrichment Analysis
- `enrichment_results.csv` - GO term enrichment results

### Visualizations (Figures)
1. `per_base_quality.png` - Per-base quality score plots
2. `gc_distribution.png` - GC content distribution
3. `trimming_stats.png` - Read trimming statistics
4. `expression_heatmap.png` - Heatmap of top expressed genes
5. `top_genes_barplot.png` - Bar plot of top expressed genes
6. `sample_correlation.png` - Sample correlation heatmap
7. `enrichment_dotplot.png` - GO term enrichment dot plot

## Pipeline Steps

1. **Download sequences from GenBank**
   - Fetches reference sequences using NCBI Entrez API
   - Uses human housekeeping genes (GAPDH, ACTB)

2. **Generate simulated FASTQ data**
   - Creates realistic RNA-seq reads
   - Includes sequencing errors and adapter contamination

3. **FastQC analysis**
   - Calculates per-base quality scores
   - Measures GC content distribution
   - Detects adapter contamination

4. **Read trimming and filtering**
   - Removes Illumina TruSeq adapters
   - Trims low-quality bases (Q<20)
   - Filters short reads (<30 bp)

5. **Expression quantification**
   - Maps reads to gene reference
   - Calculates raw counts, TPM, and CPM

6. **Top gene identification**
   - Ranks genes by mean expression
   - Identifies top 15 expressed genes

7. **Visualization**
   - Generates 6+ publication-ready figures

8. **Enrichment analysis**
   - Performs GO term enrichment
   - Calculates fold enrichment and FDR

## Dependencies

- Python 3.8+
- pandas
- numpy
- matplotlib
- seaborn
- biopython
- scipy

## License

MIT License