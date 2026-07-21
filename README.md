# TFBS_TE_Pipeline
It's a standalone pipeline that builds a single JBrowse 2–ready GFF3 track combining gene structure, transposable elements (TEs), and predicted transcription-factor binding sites (TFBS) for a set of Cyp (cytochrome P450) genes flagged as TE-affected.

Pipeline stages:

Parse TE hits (build_tfbs_te_gff.py:258) — reads a RepeatMasker-derived table (<seqid> TAB <gene> TAB RepeatMasker .out fields) and de-duplicates rows (the format commonly repeats 2-3x).

Parse the real genome annotation GFF3 (build_tfbs_te_gff.py:422) — a 3-pass scan pulling gene → mRNA → exon records for just the Cyp genes named in the TE table, picking one representative transcript per gene and deriving intron features from the exon gaps. Handles multi-paralog merged loci (comma-joined IDs).

Build scan windows (build_tfbs_te_gff.py:606) — a promoter window around each gene's TSS (strand-aware) plus a flanked window around each TE hit, to search for motifs.

Extract sequence + run motif scanning:

Sequence comes from either a local FASTA (via a dependency-free, samtools-faidx-style random-access indexer, build_tfbs_te_gff.py:649) or NCBI E-utils by accession/coordinates.
Motifs come from the JASPAR CORE Insects collection (auto-downloaded/cached), plus a hand-built literature-derived PFM for the CncC:Maf-S antioxidant response element (build_tfbs_te_gff.py:841), since that TF pair (the master regulator of insect Cyp-mediated detox) isn't in JASPAR's insect set.
Runs MEME Suite's fimo (local binary or auto-fallback to a Docker container) to find motif hits, then remaps FIMO's window-local coordinates back to absolute genome coordinates.
Combine into one GFF3 (build_tfbs_te_gff.py:1236) — writes gene → mRNA → exon/intron, mobile_genetic_element (TEs), and TF_binding_site (TFBS) records, sorted and with ##sequence-region headers. Optionally also writes a standalone TE-only GFF3 for a separate JBrowse track/color.

Optional bgzip + tabix indexing (build_tfbs_te_gff.py:1332) for direct JBrowse 2 GFF3Tabix loading, again with local-binary-or-Docker auto-detection.

It's designed so a minimal invocation (--te-hits, --gff, --fasta) auto-detects everything else (FIMO engine, JBrowse indexing, output filename), while supporting a per-species INI config file (build_tfbs_te_gff.py:1462) so the same pipeline can be rerun across species without retyping every flag.
