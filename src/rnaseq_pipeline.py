#!/usr/bin/env python3
"""
RNA-Seq Analysis Pipeline

This script performs a complete RNA-seq analysis workflow:
1. Downloads sequences from GenBank/NCBI
2. Generates simulated FASTQ data for demonstration
3. Runs FastQC-like quality assessment
4. Trims adapters and filters low-quality reads
5. Quantifies expression levels
6. Generates TPM/CPM count tables
7. Identifies top expressed genes
8. Creates visualizations (heatmaps, plots)
9. Performs enrichment analysis
"""

import os
import sys
import random
import gzip
import json
from datetime import datetime
from collections import Counter
from typing import Dict, List, Tuple, Optional
import urllib.request
import urllib.error

# Check for required dependencies
try:
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy import stats
    from Bio import Entrez, SeqIO
except ImportError as e:
    missing_module = str(e).split("'")[1] if "'" in str(e) else str(e)
    print(f"Error: Missing required dependency: {missing_module}")
    print("\nPlease install the required dependencies:")
    print("  pip install -r requirements.txt")
    print("\nOr install individually:")
    print("  pip install pandas numpy matplotlib seaborn biopython scipy")
    sys.exit(1)


class GenBankDownloader:
    """Class to download sequences from GenBank."""
    
    def __init__(self, email: str = "user@example.com"):
        Entrez.email = email
        self.sequences = {}
    
    def fetch_sequence(self, accession: str) -> Optional[str]:
        """Fetch a sequence from GenBank by accession number."""
        try:
            print(f"Fetching sequence {accession} from GenBank...")
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            self.sequences[accession] = str(record.seq)
            print(f"  Successfully fetched {len(self.sequences[accession])} bp")
            return str(record.seq)
        except Exception as e:
            print(f"  Error fetching {accession}: {e}")
            return None
    
    def fetch_multiple(self, accessions: List[str]) -> Dict[str, str]:
        """Fetch multiple sequences from GenBank."""
        for acc in accessions:
            self.fetch_sequence(acc)
        return self.sequences


class FASTQGenerator:
    """Generate simulated FASTQ reads from reference sequences."""
    
    def __init__(self, read_length: int = 100, num_reads: int = 50000):
        self.read_length = read_length
        self.num_reads = num_reads
        self.quality_chars = "FFGGHHHIIIIJJJKKKKLLLL"  # Phred+33 quality scores
        
    def generate_quality_string(self, length: int, mean_quality: float = 35) -> str:
        """Generate a quality string for a read."""
        qualities = []
        for i in range(length):
            # Quality tends to decrease towards the end of reads
            position_penalty = (i / length) * 5
            q = int(max(20, min(40, mean_quality - position_penalty + random.gauss(0, 3))))
            qualities.append(chr(q + 33))
        return ''.join(qualities)
    
    def add_errors(self, sequence: str, error_rate: float = 0.01) -> str:
        """Add sequencing errors to a read."""
        bases = list(sequence)
        for i in range(len(bases)):
            if random.random() < error_rate:
                bases[i] = random.choice(['A', 'T', 'G', 'C'])
        return ''.join(bases)
    
    def generate_reads(self, reference: str, sample_name: str, 
                       expression_level: float = 1.0) -> List[Tuple[str, str, str]]:
        """Generate simulated reads from a reference sequence."""
        reads = []
        actual_num_reads = int(self.num_reads * expression_level)
        
        # Common adapters (Illumina TruSeq)
        adapter_seq = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
        
        for i in range(actual_num_reads):
            # Random start position
            if len(reference) > self.read_length:
                start = random.randint(0, len(reference) - self.read_length)
            else:
                start = 0
            
            # Extract read
            read = reference[start:start + self.read_length]
            
            # Pad if necessary
            if len(read) < self.read_length:
                read = read + 'N' * (self.read_length - len(read))
            
            # Add adapter contamination to ~5% of reads
            if random.random() < 0.05:
                adapter_len = random.randint(5, 20)
                read = read[:-adapter_len] + adapter_seq[:adapter_len]
            
            # Add sequencing errors
            read = self.add_errors(read)
            
            # Generate quality string
            quality = self.generate_quality_string(len(read))
            
            # Read ID
            read_id = f"@{sample_name}_read_{i+1}"
            
            reads.append((read_id, read, quality))
        
        return reads
    
    def write_fastq(self, reads: List[Tuple[str, str, str]], 
                    output_path: str, compress: bool = False):
        """Write reads to a FASTQ file."""
        if compress:
            with gzip.open(output_path, 'wt') as f:
                for read_id, seq, qual in reads:
                    f.write(f"{read_id}\n{seq}\n+\n{qual}\n")
        else:
            with open(output_path, 'w') as f:
                for read_id, seq, qual in reads:
                    f.write(f"{read_id}\n{seq}\n+\n{qual}\n")
        print(f"  Written {len(reads)} reads to {output_path}")


class FastQCAnalyzer:
    """Perform FastQC-like quality control analysis."""
    
    def __init__(self):
        self.results = {}
    
    def analyze_fastq(self, fastq_path: str, sample_name: str) -> Dict:
        """Analyze a FASTQ file for quality metrics."""
        print(f"Analyzing {sample_name}...")
        
        # Read FASTQ file
        reads = []
        qualities = []
        gc_contents = []
        
        open_func = gzip.open if fastq_path.endswith('.gz') else open
        mode = 'rt' if fastq_path.endswith('.gz') else 'r'
        
        with open_func(fastq_path, mode) as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()  # + line
                qual = f.readline().strip()
                
                reads.append(seq)
                qualities.append([ord(c) - 33 for c in qual])
                
                # Calculate GC content
                if len(seq) > 0:
                    gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
                    gc_contents.append(gc)
        
        # Calculate metrics
        all_qualities = [q for qual_list in qualities for q in qual_list]
        read_lengths = [len(r) for r in reads]
        
        # Adapter detection
        adapter_seq = "AGATCGGAAGAGC"
        adapter_count = sum(1 for r in reads if adapter_seq in r)
        adapter_contamination = adapter_count / len(reads) * 100 if reads else 0
        
        # Per-base quality
        if qualities:
            max_len = max(len(q) for q in qualities)
            per_base_quality = []
            for pos in range(max_len):
                pos_quals = [q[pos] for q in qualities if pos < len(q)]
                if pos_quals:
                    per_base_quality.append({
                        'position': pos + 1,
                        'mean': np.mean(pos_quals),
                        'median': np.median(pos_quals),
                        'q1': np.percentile(pos_quals, 25),
                        'q3': np.percentile(pos_quals, 75)
                    })
        else:
            per_base_quality = []
        
        results = {
            'sample_name': sample_name,
            'total_reads': len(reads),
            'total_bases': sum(read_lengths),
            'mean_read_length': np.mean(read_lengths) if read_lengths else 0,
            'mean_quality': np.mean(all_qualities) if all_qualities else 0,
            'median_quality': np.median(all_qualities) if all_qualities else 0,
            'mean_gc_content': np.mean(gc_contents) if gc_contents else 0,
            'gc_std': np.std(gc_contents) if gc_contents else 0,
            'adapter_contamination_pct': adapter_contamination,
            'per_base_quality': per_base_quality,
            'gc_distribution': gc_contents,
            'quality_distribution': all_qualities
        }
        
        self.results[sample_name] = results
        return results
    
    def generate_summary_report(self, output_dir: str) -> str:
        """Generate a summary report of all samples."""
        report_lines = [
            "=" * 80,
            "FASTQC SUMMARY REPORT",
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "=" * 80,
            ""
        ]
        
        for sample_name, results in self.results.items():
            report_lines.extend([
                f"\nSample: {sample_name}",
                "-" * 40,
                f"  Total Reads: {results['total_reads']:,}",
                f"  Total Bases: {results['total_bases']:,}",
                f"  Mean Read Length: {results['mean_read_length']:.1f} bp",
                f"  Mean Quality Score: {results['mean_quality']:.1f}",
                f"  Median Quality Score: {results['median_quality']:.1f}",
                f"  Mean GC Content: {results['mean_gc_content']:.1f}%",
                f"  GC Std Dev: {results['gc_std']:.2f}%",
                f"  Adapter Contamination: {results['adapter_contamination_pct']:.2f}%",
                "",
                "  Quality Assessment:",
            ])
            
            # Quality assessment
            if results['mean_quality'] >= 30:
                report_lines.append("    ✓ PASS: High quality reads (Q>=30)")
            elif results['mean_quality'] >= 20:
                report_lines.append("    ⚠ WARN: Moderate quality reads (Q>=20)")
            else:
                report_lines.append("    ✗ FAIL: Low quality reads (Q<20)")
            
            # GC content assessment
            if 40 <= results['mean_gc_content'] <= 60:
                report_lines.append("    ✓ PASS: Normal GC content (40-60%)")
            else:
                report_lines.append("    ⚠ WARN: Unusual GC content")
            
            # Adapter contamination assessment
            if results['adapter_contamination_pct'] < 1:
                report_lines.append("    ✓ PASS: Low adapter contamination (<1%)")
            elif results['adapter_contamination_pct'] < 5:
                report_lines.append("    ⚠ WARN: Moderate adapter contamination (1-5%)")
            else:
                report_lines.append("    ✗ FAIL: High adapter contamination (>5%)")
        
        report = "\n".join(report_lines)
        
        # Save report
        report_path = os.path.join(output_dir, "fastqc_summary.txt")
        with open(report_path, 'w') as f:
            f.write(report)
        
        print(f"\nFastQC Summary Report saved to: {report_path}")
        return report


class ReadTrimmer:
    """Trim adapters and filter low-quality reads (fastp-like functionality)."""
    
    def __init__(self, min_quality: int = 20, min_length: int = 30,
                 adapter_seq: str = "AGATCGGAAGAGC"):
        self.min_quality = min_quality
        self.min_length = min_length
        self.adapter_seq = adapter_seq
        self.stats = {}
    
    def trim_adapter(self, sequence: str, quality: str) -> Tuple[str, str]:
        """Remove adapter sequences from a read."""
        # Simple adapter trimming
        for i in range(len(self.adapter_seq), 5, -1):
            adapter_part = self.adapter_seq[:i]
            if adapter_part in sequence:
                idx = sequence.find(adapter_part)
                sequence = sequence[:idx]
                quality = quality[:idx]
                break
        return sequence, quality
    
    def quality_trim(self, sequence: str, quality: str) -> Tuple[str, str]:
        """Trim low-quality bases from the 3' end."""
        while len(quality) > 0:
            last_qual = ord(quality[-1]) - 33
            if last_qual < self.min_quality:
                sequence = sequence[:-1]
                quality = quality[:-1]
            else:
                break
        
        # Also trim from 5' end if needed
        while len(quality) > 0:
            first_qual = ord(quality[0]) - 33
            if first_qual < self.min_quality:
                sequence = sequence[1:]
                quality = quality[1:]
            else:
                break
        
        return sequence, quality
    
    def process_fastq(self, input_path: str, output_path: str, 
                      sample_name: str) -> Dict:
        """Process a FASTQ file: trim adapters and filter reads."""
        print(f"Trimming and filtering {sample_name}...")
        
        open_func_in = gzip.open if input_path.endswith('.gz') else open
        mode_in = 'rt' if input_path.endswith('.gz') else 'r'
        
        total_reads = 0
        passed_reads = 0
        adapter_trimmed = 0
        quality_trimmed = 0
        too_short = 0
        
        trimmed_reads = []
        
        with open_func_in(input_path, mode_in) as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break
                seq = f.readline().strip()
                plus = f.readline().strip()
                qual = f.readline().strip()
                
                total_reads += 1
                original_len = len(seq)
                
                # Trim adapter
                seq, qual = self.trim_adapter(seq, qual)
                if len(seq) < original_len:
                    adapter_trimmed += 1
                
                # Quality trim
                trimmed_len = len(seq)
                seq, qual = self.quality_trim(seq, qual)
                if len(seq) < trimmed_len:
                    quality_trimmed += 1
                
                # Length filter
                if len(seq) >= self.min_length:
                    trimmed_reads.append((header, seq, qual))
                    passed_reads += 1
                else:
                    too_short += 1
        
        # Write output
        with open(output_path, 'w') as f:
            for header, seq, qual in trimmed_reads:
                f.write(f"{header}\n{seq}\n+\n{qual}\n")
        
        stats = {
            'sample_name': sample_name,
            'total_reads': total_reads,
            'passed_reads': passed_reads,
            'adapter_trimmed': adapter_trimmed,
            'quality_trimmed': quality_trimmed,
            'too_short': too_short,
            'pass_rate': passed_reads / total_reads * 100 if total_reads > 0 else 0
        }
        
        self.stats[sample_name] = stats
        
        print(f"  Total reads: {total_reads:,}")
        print(f"  Passed reads: {passed_reads:,} ({stats['pass_rate']:.1f}%)")
        print(f"  Adapter trimmed: {adapter_trimmed:,}")
        print(f"  Quality trimmed: {quality_trimmed:,}")
        print(f"  Too short: {too_short:,}")
        
        return stats
    
    def generate_report(self, output_dir: str) -> str:
        """Generate trimming report."""
        report_lines = [
            "=" * 80,
            "READ TRIMMING REPORT (fastp-like)",
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "=" * 80,
            f"\nTrimming Parameters:",
            f"  Minimum Quality: {self.min_quality}",
            f"  Minimum Length: {self.min_length} bp",
            f"  Adapter Sequence: {self.adapter_seq}",
            ""
        ]
        
        for sample_name, stats in self.stats.items():
            report_lines.extend([
                f"\nSample: {sample_name}",
                "-" * 40,
                f"  Input Reads: {stats['total_reads']:,}",
                f"  Passed Reads: {stats['passed_reads']:,} ({stats['pass_rate']:.1f}%)",
                f"  Adapter Trimmed: {stats['adapter_trimmed']:,}",
                f"  Quality Trimmed: {stats['quality_trimmed']:,}",
                f"  Too Short (removed): {stats['too_short']:,}",
            ])
        
        report = "\n".join(report_lines)
        report_path = os.path.join(output_dir, "trimming_report.txt")
        with open(report_path, 'w') as f:
            f.write(report)
        
        print(f"\nTrimming Report saved to: {report_path}")
        return report


class ExpressionQuantifier:
    """Quantify gene expression (Salmon/featureCounts-like functionality)."""
    
    def __init__(self):
        self.counts = {}
        self.tpm = {}
        self.cpm = {}
        self.gene_lengths = {}  # Store actual gene lengths for TPM calculation
    
    def quantify_expression(self, fastq_path: str, reference_genes: Dict[str, str],
                            sample_name: str) -> Dict[str, int]:
        """Quantify expression by mapping reads to reference genes."""
        print(f"Quantifying expression for {sample_name}...")
        
        # Store gene lengths for TPM calculation
        for gene, seq in reference_genes.items():
            self.gene_lengths[gene] = len(seq)
        
        # Build kmer index for efficient matching
        kmer_size = 15
        gene_kmer_index = {}  # kmer -> list of gene names
        
        for gene, ref_seq in reference_genes.items():
            for i in range(0, len(ref_seq) - kmer_size + 1, 50):  # Sample every 50bp
                kmer = ref_seq[i:i + kmer_size]
                if kmer not in gene_kmer_index:
                    gene_kmer_index[kmer] = []
                if gene not in gene_kmer_index[kmer]:
                    gene_kmer_index[kmer].append(gene)
        
        # Read FASTQ file
        gene_counts = {gene: 0 for gene in reference_genes}
        read_count = 0
        
        open_func = gzip.open if fastq_path.endswith('.gz') else open
        mode = 'rt' if fastq_path.endswith('.gz') else 'r'
        
        with open_func(fastq_path, mode) as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()
                f.readline()
                read_count += 1
                
                # Check first kmer of read against index
                if len(seq) >= kmer_size:
                    read_kmer = seq[:kmer_size]
                    if read_kmer in gene_kmer_index:
                        for gene in gene_kmer_index[read_kmer]:
                            gene_counts[gene] += 1
        
        print(f"  Processed {read_count:,} reads")
        self.counts[sample_name] = gene_counts
        return gene_counts
    
    def calculate_tpm(self) -> pd.DataFrame:
        """
        Calculate Transcripts Per Million (TPM).
        
        TPM = (reads / gene_length_kb) / (sum(reads / gene_length_kb)) * 1e6
        
        Uses actual gene lengths from reference sequences when available.
        """
        # Create counts dataframe
        counts_df = pd.DataFrame(self.counts)
        
        # Get gene lengths (in kilobases)
        gene_lengths_kb = pd.Series({
            gene: length / 1000 for gene, length in self.gene_lengths.items()
        })
        
        # Handle genes with no length info (use median length)
        missing_genes = set(counts_df.index) - set(gene_lengths_kb.index)
        if missing_genes:
            median_length_kb = gene_lengths_kb.median() if len(gene_lengths_kb) > 0 else 1.0
            for gene in missing_genes:
                gene_lengths_kb[gene] = median_length_kb
        
        # RPK (reads per kilobase) - normalize by gene length
        rpk = counts_df.div(gene_lengths_kb, axis=0)
        
        # Scaling factor (sum of all RPK per sample, in millions)
        scaling_factor = rpk.sum() / 1e6
        
        # TPM
        tpm_df = rpk / scaling_factor
        tpm_df = tpm_df.fillna(0)
        
        self.tpm = tpm_df
        return tpm_df
    
    def calculate_cpm(self) -> pd.DataFrame:
        """Calculate Counts Per Million (CPM)."""
        counts_df = pd.DataFrame(self.counts)
        
        # Library size (total counts per sample)
        library_size = counts_df.sum()
        
        # CPM
        cpm_df = (counts_df / library_size) * 1e6
        cpm_df = cpm_df.fillna(0)
        
        self.cpm = cpm_df
        return cpm_df
    
    def get_top_expressed_genes(self, n: int = 20) -> pd.DataFrame:
        """Get the top N expressed genes across all samples."""
        if self.tpm is None or len(self.tpm) == 0:
            self.calculate_tpm()
        
        # Calculate mean expression across samples
        mean_expression = self.tpm.mean(axis=1)
        top_genes = mean_expression.nlargest(n)
        
        return self.tpm.loc[top_genes.index]
    
    def save_count_tables(self, output_dir: str):
        """Save TPM and CPM count tables."""
        # Raw counts
        counts_df = pd.DataFrame(self.counts)
        counts_path = os.path.join(output_dir, "raw_counts.csv")
        counts_df.to_csv(counts_path)
        print(f"Raw counts saved to: {counts_path}")
        
        # TPM
        if self.tpm is not None and len(self.tpm) > 0:
            tpm_path = os.path.join(output_dir, "tpm_counts.csv")
            self.tpm.to_csv(tpm_path)
            print(f"TPM counts saved to: {tpm_path}")
        
        # CPM
        if self.cpm is not None and len(self.cpm) > 0:
            cpm_path = os.path.join(output_dir, "cpm_counts.csv")
            self.cpm.to_csv(cpm_path)
            print(f"CPM counts saved to: {cpm_path}")


class Visualizer:
    """Create visualizations for RNA-seq analysis."""
    
    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        self.figures_created = []
        plt.style.use('seaborn-v0_8-whitegrid')
    
    def plot_per_base_quality(self, fastqc_results: Dict, filename: str = "per_base_quality.png"):
        """Plot per-base quality scores."""
        fig, axes = plt.subplots(1, len(fastqc_results), figsize=(6*len(fastqc_results), 5))
        
        if len(fastqc_results) == 1:
            axes = [axes]
        
        for ax, (sample_name, results) in zip(axes, fastqc_results.items()):
            per_base = results['per_base_quality']
            positions = [p['position'] for p in per_base]
            means = [p['mean'] for p in per_base]
            q1s = [p['q1'] for p in per_base]
            q3s = [p['q3'] for p in per_base]
            
            ax.fill_between(positions, q1s, q3s, alpha=0.3, color='blue', label='IQR')
            ax.plot(positions, means, color='blue', linewidth=2, label='Mean')
            
            # Quality thresholds
            ax.axhline(y=30, color='green', linestyle='--', alpha=0.7, label='Q30')
            ax.axhline(y=20, color='orange', linestyle='--', alpha=0.7, label='Q20')
            
            ax.set_xlabel('Position in Read (bp)')
            ax.set_ylabel('Quality Score')
            ax.set_title(f'Per Base Quality - {sample_name}')
            ax.legend(loc='lower left')
            ax.set_ylim(0, 45)
        
        plt.tight_layout()
        filepath = os.path.join(self.output_dir, filename)
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        self.figures_created.append(filepath)
        print(f"  Saved: {filename}")
    
    def plot_gc_distribution(self, fastqc_results: Dict, filename: str = "gc_distribution.png"):
        """Plot GC content distribution."""
        fig, ax = plt.subplots(figsize=(10, 6))
        
        colors = plt.cm.Set2(np.linspace(0, 1, len(fastqc_results)))
        
        for (sample_name, results), color in zip(fastqc_results.items(), colors):
            gc_dist = results['gc_distribution']
            if gc_dist:
                ax.hist(gc_dist, bins=50, alpha=0.5, label=sample_name, color=color, density=True)
        
        ax.set_xlabel('GC Content (%)')
        ax.set_ylabel('Density')
        ax.set_title('GC Content Distribution')
        ax.legend()
        
        plt.tight_layout()
        filepath = os.path.join(self.output_dir, filename)
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        self.figures_created.append(filepath)
        print(f"  Saved: {filename}")
    
    def plot_trimming_stats(self, trimming_stats: Dict, filename: str = "trimming_stats.png"):
        """Plot trimming statistics."""
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        samples = list(trimming_stats.keys())
        
        # Bar plot of reads before/after
        total_reads = [trimming_stats[s]['total_reads'] for s in samples]
        passed_reads = [trimming_stats[s]['passed_reads'] for s in samples]
        
        x = np.arange(len(samples))
        width = 0.35
        
        axes[0].bar(x - width/2, total_reads, width, label='Total Reads', color='lightblue')
        axes[0].bar(x + width/2, passed_reads, width, label='Passed Reads', color='darkblue')
        axes[0].set_xlabel('Sample')
        axes[0].set_ylabel('Number of Reads')
        axes[0].set_title('Reads Before and After Trimming')
        axes[0].set_xticks(x)
        axes[0].set_xticklabels(samples, rotation=45, ha='right')
        axes[0].legend()
        
        # Pie chart of filtering reasons (for first sample)
        first_sample = samples[0]
        stats = trimming_stats[first_sample]
        
        labels = ['Passed', 'Adapter Trimmed (passed)', 'Quality Trimmed (passed)', 'Too Short (removed)']
        sizes = [
            stats['passed_reads'] - stats['adapter_trimmed'] - stats['quality_trimmed'],
            stats['adapter_trimmed'],
            stats['quality_trimmed'],
            stats['too_short']
        ]
        # Ensure no negative values
        sizes = [max(0, s) for s in sizes]
        
        colors = ['green', 'yellow', 'orange', 'red']
        axes[1].pie([s for s in sizes if s > 0], 
                    labels=[l for l, s in zip(labels, sizes) if s > 0],
                    colors=[c for c, s in zip(colors, sizes) if s > 0],
                    autopct='%1.1f%%', startangle=90)
        axes[1].set_title(f'Read Processing Summary - {first_sample}')
        
        plt.tight_layout()
        filepath = os.path.join(self.output_dir, filename)
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        self.figures_created.append(filepath)
        print(f"  Saved: {filename}")
    
    def plot_expression_heatmap(self, expression_df: pd.DataFrame, 
                                 filename: str = "expression_heatmap.png",
                                 title: str = "Top Expressed Genes"):
        """Create a heatmap of gene expression."""
        if expression_df.empty:
            print("  No expression data for heatmap")
            return
        
        # Log transform for better visualization
        log_expr = np.log2(expression_df + 1)
        
        fig, ax = plt.subplots(figsize=(10, max(6, len(expression_df) * 0.3)))
        
        sns.heatmap(log_expr, annot=True, fmt='.1f', cmap='RdYlBu_r',
                    ax=ax, cbar_kws={'label': 'log2(TPM + 1)'})
        
        ax.set_title(title)
        ax.set_xlabel('Sample')
        ax.set_ylabel('Gene')
        
        plt.tight_layout()
        filepath = os.path.join(self.output_dir, filename)
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        self.figures_created.append(filepath)
        print(f"  Saved: {filename}")
    
    def plot_top_genes_barplot(self, expression_df: pd.DataFrame,
                                filename: str = "top_genes_barplot.png"):
        """Create a bar plot of top expressed genes."""
        if expression_df.empty:
            print("  No expression data for bar plot")
            return
        
        # Calculate mean expression
        mean_expr = expression_df.mean(axis=1).sort_values(ascending=True)
        
        fig, ax = plt.subplots(figsize=(10, max(6, len(mean_expr) * 0.3)))
        
        colors = plt.cm.viridis(np.linspace(0, 1, len(mean_expr)))
        ax.barh(range(len(mean_expr)), mean_expr.values, color=colors)
        ax.set_yticks(range(len(mean_expr)))
        ax.set_yticklabels(mean_expr.index)
        ax.set_xlabel('Mean TPM')
        ax.set_title('Top Expressed Genes (Mean TPM)')
        
        plt.tight_layout()
        filepath = os.path.join(self.output_dir, filename)
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        self.figures_created.append(filepath)
        print(f"  Saved: {filename}")
    
    def plot_sample_correlation(self, expression_df: pd.DataFrame,
                                 filename: str = "sample_correlation.png"):
        """Create a correlation heatmap between samples."""
        if expression_df.empty or len(expression_df.columns) < 2:
            print("  Need at least 2 samples for correlation plot")
            return
        
        # Calculate correlation
        corr = expression_df.corr()
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        sns.heatmap(corr, annot=True, fmt='.3f', cmap='coolwarm',
                    ax=ax, vmin=0, vmax=1, center=0.5,
                    square=True, linewidths=0.5)
        
        ax.set_title('Sample Correlation (TPM)')
        
        plt.tight_layout()
        filepath = os.path.join(self.output_dir, filename)
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        self.figures_created.append(filepath)
        print(f"  Saved: {filename}")


class EnrichmentAnalyzer:
    """
    Perform GO term enrichment analysis (DEMONSTRATION MODE).
    
    NOTE: This implementation uses simulated GO term assignments for demonstration
    purposes. In a production environment, you would integrate with real gene
    annotation databases such as:
    - Gene Ontology (GO) database
    - KEGG pathways
    - Reactome
    - MSigDB gene sets
    
    For real enrichment analysis, consider using tools like:
    - gseapy (Python)
    - clusterProfiler (R)
    - DAVID
    - Enrichr
    """
    
    def __init__(self):
        # Real GO terms (these are actual GO IDs and descriptions)
        # In production, these would come from a GO database
        self.go_terms = {
            'GO:0006412': 'translation',
            'GO:0006414': 'translational elongation',
            'GO:0006396': 'RNA processing',
            'GO:0016070': 'RNA metabolic process',
            'GO:0006397': 'mRNA processing',
            'GO:0006355': 'regulation of transcription',
            'GO:0006351': 'transcription',
            'GO:0006950': 'response to stress',
            'GO:0006979': 'response to oxidative stress',
            'GO:0009058': 'biosynthetic process'
        }
        
        self.gene_go_mapping = {}
    
    def assign_random_go_terms(self, genes: List[str]):
        """
        Assign random GO terms to genes for demonstration purposes.
        
        WARNING: This produces simulated annotations. In production, use real
        gene annotation databases (GO, KEGG, Reactome, etc.)
        """
        for gene in genes:
            # Assign 2-5 random GO terms to each gene
            num_terms = random.randint(2, 5)
            self.gene_go_mapping[gene] = random.sample(list(self.go_terms.keys()), num_terms)
    
    def perform_enrichment(self, top_genes: List[str], all_genes: List[str],
                           output_dir: str) -> pd.DataFrame:
        """
        Perform GO term enrichment analysis.
        
        DEMONSTRATION MODE: Uses simulated GO term assignments.
        Results are for illustrative purposes only.
        
        For production use, integrate with real annotation databases.
        """
        print("Performing enrichment analysis (DEMONSTRATION MODE)...")
        print("  Note: Using simulated GO term assignments for illustration")
        
        # Assign GO terms if not already done
        if not self.gene_go_mapping:
            self.assign_random_go_terms(all_genes)
        
        # Count GO terms in top genes vs all genes
        results = []
        
        for go_id, go_name in self.go_terms.items():
            # Count in top genes
            top_count = sum(1 for g in top_genes if go_id in self.gene_go_mapping.get(g, []))
            # Count in all genes
            all_count = sum(1 for g in all_genes if go_id in self.gene_go_mapping.get(g, []))
            
            # Calculate enrichment (Fisher's exact test approximation)
            if all_count > 0:
                expected = (len(top_genes) / len(all_genes)) * all_count
                fold_enrichment = top_count / expected if expected > 0 else 0
                
                # Simple p-value approximation
                if top_count > 0:
                    # Use binomial approximation
                    p = all_count / len(all_genes)
                    n = len(top_genes)
                    pvalue = 1 - stats.binom.cdf(top_count - 1, n, p)
                else:
                    pvalue = 1.0
                
                results.append({
                    'GO_ID': go_id,
                    'GO_Term': go_name,
                    'Top_Genes_Count': top_count,
                    'All_Genes_Count': all_count,
                    'Fold_Enrichment': fold_enrichment,
                    'P_Value': pvalue
                })
        
        # Create dataframe and sort by p-value
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values('P_Value')
        
        # Add FDR correction (Benjamini-Hochberg)
        results_df['FDR'] = results_df['P_Value'] * len(results_df) / (results_df.reset_index().index + 1)
        results_df['FDR'] = results_df['FDR'].clip(upper=1.0)
        
        # Save results
        enrichment_path = os.path.join(output_dir, "enrichment_results.csv")
        results_df.to_csv(enrichment_path, index=False)
        print(f"Enrichment results saved to: {enrichment_path}")
        
        return results_df
    
    def plot_enrichment(self, results_df: pd.DataFrame, output_dir: str,
                        filename: str = "enrichment_dotplot.png"):
        """Create an enrichment dot plot."""
        if results_df.empty:
            print("  No enrichment results to plot")
            return
        
        # Filter significant results
        sig_results = results_df[results_df['FDR'] < 0.5].head(10)
        
        if sig_results.empty:
            sig_results = results_df.head(10)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Create dot plot
        x = sig_results['Fold_Enrichment']
        y = range(len(sig_results))
        sizes = sig_results['Top_Genes_Count'] * 50 + 20
        colors = -np.log10(sig_results['FDR'] + 1e-10)
        
        scatter = ax.scatter(x, y, s=sizes, c=colors, cmap='RdYlBu_r', 
                            alpha=0.7, edgecolors='black', linewidths=0.5)
        
        ax.set_yticks(y)
        ax.set_yticklabels(sig_results['GO_Term'])
        ax.set_xlabel('Fold Enrichment')
        ax.set_title('GO Term Enrichment Analysis')
        
        # Add colorbar
        cbar = plt.colorbar(scatter)
        cbar.set_label('-log10(FDR)')
        
        # Add size legend
        sizes_legend = [1, 3, 5]
        legend_elements = [plt.scatter([], [], s=s*50+20, c='gray', alpha=0.7,
                                       label=f'{s} genes') for s in sizes_legend]
        ax.legend(handles=legend_elements, loc='lower right', title='Gene Count')
        
        plt.tight_layout()
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {filename}")
        
        return filepath


def create_simulated_gene_reference() -> Dict[str, str]:
    """Create a simulated gene reference for expression quantification."""
    genes = {}
    gene_names = [
        'GAPDH', 'ACTB', 'RPL13A', 'RPS18', 'EEF1A1',
        'PPIA', 'HSP90AB1', 'LDHA', 'PKM', 'ENO1',
        'TUBA1A', 'TUBB', 'VIM', 'HSPA8', 'ATP5F1B',
        'CCT2', 'EIF4A1', 'HNRNPA1', 'NPM1', 'NCL',
        'RAN', 'YWHAZ', 'CALM1', 'UBC', 'RACK1',
        'DDX5', 'PTBP1', 'SFPQ', 'HNRNPK', 'FUS'
    ]
    
    for gene in gene_names:
        # Generate random sequence (1000-3000 bp)
        length = random.randint(1000, 3000)
        seq = ''.join(random.choices(['A', 'T', 'G', 'C'], 
                                      weights=[0.3, 0.3, 0.2, 0.2], k=length))
        genes[gene] = seq
    
    return genes


def main(genbank_accessions: List[str] = None):
    """
    Main pipeline execution.
    
    Args:
        genbank_accessions: Optional list of GenBank accession numbers to download.
                           Defaults to ['NM_002046.7', 'NM_001101.5'] (GAPDH and ACTB).
    """
    print("=" * 80)
    print("RNA-SEQ ANALYSIS PIPELINE")
    print("=" * 80)
    print()
    
    # Setup directories
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, 'data')
    results_dir = os.path.join(base_dir, 'results')
    
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    
    # =========================================================================
    # STEP 1: Download sequences from GenBank
    # =========================================================================
    print("\n" + "=" * 80)
    print("STEP 1: Downloading sequences from GenBank")
    print("=" * 80)
    
    downloader = GenBankDownloader()
    
    # Use default accessions if none provided
    # Default: human housekeeping genes GAPDH and ACTB
    if genbank_accessions is None:
        genbank_accessions = ['NM_002046.7', 'NM_001101.5']
    
    print(f"Target accessions: {', '.join(genbank_accessions)}")
    
    try:
        sequences = downloader.fetch_multiple(genbank_accessions)
        print(f"\nDownloaded {len(sequences)} sequences from GenBank")
    except Exception as e:
        print(f"Error downloading from GenBank: {e}")
        print("Using simulated reference sequences instead...")
        sequences = {
            'NM_002046.7': ''.join(random.choices(['A','T','G','C'], k=1500)),
            'NM_001101.5': ''.join(random.choices(['A','T','G','C'], k=1800))
        }
    
    # =========================================================================
    # STEP 2: Generate simulated FASTQ data
    # =========================================================================
    print("\n" + "=" * 80)
    print("STEP 2: Generating simulated FASTQ data")
    print("=" * 80)
    
    fastq_generator = FASTQGenerator(read_length=100, num_reads=5000)
    
    # Create gene reference for quantification
    gene_reference = create_simulated_gene_reference()
    
    # Add downloaded sequences to reference
    for acc, seq in sequences.items():
        gene_reference[f"GenBank_{acc}"] = seq
    
    # Generate reads for two samples with different expression levels
    sample_configs = {
        'Sample_1': {'expression_multipliers': {'GAPDH': 2.0, 'ACTB': 1.5, 'RPL13A': 1.8}},
        'Sample_2': {'expression_multipliers': {'GAPDH': 1.0, 'ACTB': 2.0, 'EEF1A1': 1.7}}
    }
    
    raw_fastq_files = {}
    
    for sample_name, config in sample_configs.items():
        all_reads = []
        print(f"\nGenerating reads for {sample_name}...")
        
        for gene_name, gene_seq in gene_reference.items():
            expr_mult = config['expression_multipliers'].get(gene_name, 1.0)
            reads = fastq_generator.generate_reads(gene_seq, sample_name, expr_mult)
            all_reads.extend(reads)
        
        # Shuffle reads
        random.shuffle(all_reads)
        
        # Write to file
        fastq_path = os.path.join(data_dir, f"{sample_name}_raw.fastq")
        fastq_generator.write_fastq(all_reads, fastq_path)
        raw_fastq_files[sample_name] = fastq_path
    
    # =========================================================================
    # STEP 3: Run FastQC analysis
    # =========================================================================
    print("\n" + "=" * 80)
    print("STEP 3: Running FastQC analysis on raw FASTQ files")
    print("=" * 80)
    
    fastqc = FastQCAnalyzer()
    
    for sample_name, fastq_path in raw_fastq_files.items():
        fastqc.analyze_fastq(fastq_path, sample_name)
    
    # Generate and print summary report
    report = fastqc.generate_summary_report(results_dir)
    print("\n" + report)
    
    # =========================================================================
    # STEP 4: Trim adapters and filter reads
    # =========================================================================
    print("\n" + "=" * 80)
    print("STEP 4: Trimming adapters and filtering reads (fastp-like)")
    print("=" * 80)
    
    trimmer = ReadTrimmer(min_quality=20, min_length=30)
    trimmed_fastq_files = {}
    
    for sample_name, fastq_path in raw_fastq_files.items():
        trimmed_path = os.path.join(data_dir, f"{sample_name}_trimmed.fastq")
        trimmer.process_fastq(fastq_path, trimmed_path, sample_name)
        trimmed_fastq_files[sample_name] = trimmed_path
    
    trimmer.generate_report(results_dir)
    
    # =========================================================================
    # STEP 5: Quantify expression
    # =========================================================================
    print("\n" + "=" * 80)
    print("STEP 5: Quantifying gene expression (Salmon-like)")
    print("=" * 80)
    
    quantifier = ExpressionQuantifier()
    
    for sample_name, fastq_path in trimmed_fastq_files.items():
        quantifier.quantify_expression(fastq_path, gene_reference, sample_name)
    
    # Calculate TPM and CPM
    tpm_df = quantifier.calculate_tpm()
    cpm_df = quantifier.calculate_cpm()
    
    # Save count tables
    quantifier.save_count_tables(results_dir)
    
    # =========================================================================
    # STEP 6: Identify top expressed genes
    # =========================================================================
    print("\n" + "=" * 80)
    print("STEP 6: Identifying top expressed genes")
    print("=" * 80)
    
    top_genes_df = quantifier.get_top_expressed_genes(n=15)
    
    print("\nTop 15 Expressed Genes (by mean TPM):")
    print("-" * 60)
    mean_expr = top_genes_df.mean(axis=1).sort_values(ascending=False)
    for i, (gene, expr) in enumerate(mean_expr.items(), 1):
        print(f"  {i:2d}. {gene:15s} - Mean TPM: {expr:.2f}")
    
    # Save top genes
    top_genes_path = os.path.join(results_dir, "top_expressed_genes.csv")
    top_genes_df.to_csv(top_genes_path)
    print(f"\nTop expressed genes saved to: {top_genes_path}")
    
    # =========================================================================
    # STEP 7: Generate visualizations
    # =========================================================================
    print("\n" + "=" * 80)
    print("STEP 7: Generating visualizations")
    print("=" * 80)
    
    viz = Visualizer(results_dir)
    
    print("\nCreating plots...")
    
    # Plot 1: Per-base quality
    viz.plot_per_base_quality(fastqc.results)
    
    # Plot 2: GC distribution
    viz.plot_gc_distribution(fastqc.results)
    
    # Plot 3: Trimming statistics
    viz.plot_trimming_stats(trimmer.stats)
    
    # Plot 4: Expression heatmap
    viz.plot_expression_heatmap(top_genes_df, title="Top 15 Expressed Genes Heatmap")
    
    # Plot 5: Top genes bar plot
    viz.plot_top_genes_barplot(top_genes_df)
    
    # Plot 6: Sample correlation
    viz.plot_sample_correlation(tpm_df)
    
    print(f"\nTotal figures created: {len(viz.figures_created)}")
    for fig in viz.figures_created:
        print(f"  - {os.path.basename(fig)}")
    
    # =========================================================================
    # STEP 8: Enrichment analysis
    # =========================================================================
    print("\n" + "=" * 80)
    print("STEP 8: Enrichment analysis")
    print("=" * 80)
    
    enrichment = EnrichmentAnalyzer()
    
    top_gene_list = list(top_genes_df.index)
    all_gene_list = list(gene_reference.keys())
    
    enrichment_results = enrichment.perform_enrichment(top_gene_list, all_gene_list, results_dir)
    
    print("\nTop Enriched GO Terms:")
    print("-" * 80)
    for _, row in enrichment_results.head(5).iterrows():
        print(f"  {row['GO_ID']:12s} {row['GO_Term']:30s} FE={row['Fold_Enrichment']:.2f} FDR={row['FDR']:.4f}")
    
    # Plot enrichment
    enrichment.plot_enrichment(enrichment_results, results_dir)
    
    # =========================================================================
    # FINAL SUMMARY
    # =========================================================================
    print("\n" + "=" * 80)
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print("=" * 80)
    
    print(f"\nOutput directory: {results_dir}")
    print("\nGenerated files:")
    for f in sorted(os.listdir(results_dir)):
        fpath = os.path.join(results_dir, f)
        size = os.path.getsize(fpath)
        print(f"  - {f} ({size:,} bytes)")
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE!")
    print("=" * 80)


if __name__ == "__main__":
    main()
