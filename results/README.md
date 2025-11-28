# Results Directory

This directory contains the output files from the RNA-seq analysis pipeline.

## Generated Files

### Reports
- `fastqc_summary.txt` - FastQC-like quality assessment report
- `trimming_report.txt` - Read trimming statistics

### Count Tables
- `raw_counts.csv` - Raw read counts per gene
- `tpm_counts.csv` - TPM (Transcripts Per Million) normalized counts
- `cpm_counts.csv` - CPM (Counts Per Million) normalized counts
- `top_expressed_genes.csv` - Top 15 expressed genes

### Enrichment Analysis
- `enrichment_results.csv` - GO term enrichment results

### Visualizations
- `per_base_quality.png` - Per-base quality score plots
- `gc_distribution.png` - GC content distribution
- `trimming_stats.png` - Read trimming statistics
- `expression_heatmap.png` - Heatmap of top expressed genes
- `top_genes_barplot.png` - Bar plot of top expressed genes
- `sample_correlation.png` - Sample correlation heatmap
- `enrichment_dotplot.png` - GO term enrichment dot plot
