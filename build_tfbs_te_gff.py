#!/usr/bin/env python3
"""
build_tfbs_te_gff.py

Builds a combined GFF3 annotation (genes, exons, transposable elements, and
predicted transcription-factor binding sites / TFBS) for a set of Cyp
(cytochrome P450) genes flagged as TE-affected.

------------------------------------------------------------------------
PIPELINE OVERVIEW
------------------------------------------------------------------------
1. Parse a RepeatMasker-derived "genes affected by TE" table (the same
   layout as D_suzukiiGenesAffectedByTE.txt: <seqid> TAB <gene> TAB
   <RepeatMasker .out fields...>) into TE records, de-duplicating rows
   (this input format commonly repeats rows 2-3x).

2. Parse a real genome annotation GFF3 (e.g. an NCBI RefSeq annotation)
   and pull out `gene`, `mRNA` and `exon` features for the Cyp genes named
   in the TE table.

3. For each gene, build search windows for motif scanning:
     - a promoter window around the TSS (--upstream / --downstream), and
     - each TE interval (+/- --te-flank bp of context).
   Sequence for these windows is extracted either from a local genome
   FASTA (indexed on the fly with a small built-in, dependency-free
   random-access reader - no samtools/pyfaidx required, verified against
   a real 140MB repeat-masked dm6 FASTA with plain accession headers like
   ">NC_004354.4") or fetched on-the-fly from NCBI E-utils by accession +
   coordinate range (no full genome download required if your TE table's
   seqids are NCBI RefSeq accessions, e.g. NW_/NT_/NC_*).

4. Download (or reuse a local copy of) the JASPAR CORE Insects
   non-redundant motif collection in MEME format, and run it through
   MEME Suite's `fimo` against the extracted sequences.

5. Remap FIMO's hit coordinates (which are local to each extracted
   window) back to absolute genome coordinates, and merge everything into
   one sorted, standard 9-column GFF3 file: `gene` -> `mRNA` (one
   representative transcript per gene) -> `exon` / `intron` (derived from
   the gaps between that transcript's exons), plus `mobile_genetic_element`
   (TE hits, always strand-aware) and `TF_binding_site` (JASPAR motif name
   in the Name= attribute, p-/q-value, matched sequence) records.

6. Optionally bgzip + tabix-index the output (--bgzip-index) so it can be
   loaded directly into JBrowse 2 as a GFF3Tabix track.

------------------------------------------------------------------------
REQUIRED SETUP (please install before running the full pipeline)
------------------------------------------------------------------------
Python packages: none required beyond the standard library for parsing,
sequence extraction (local FASTA or NCBI), and GFF3 output. `requests` is
optional/unused (urllib handles downloads); install it only if you plan
to extend the script.

MEME Suite (provides the `fimo` binary used for motif scanning):
    # Easiest - via conda/mamba (recommended):
    conda install -c bioconda meme

    # Debian/Ubuntu apt (may be an older version, check `fimo --version`):
    sudo apt-get update && sudo apt-get install meme

    # Or build from source (latest version, needed on some HPC systems):
    # https://meme-suite.org/meme/doc/install.html?man_type=web
    #   wget https://meme-suite.org/meme/meme-software/latest/meme-<ver>.tar.gz
    #   tar zxf meme-<ver>.tar.gz && cd meme-<ver>
    #   ./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt
    #   make && make test && make install
    #   export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-<ver>:$PATH

Verify installation:
    fimo --version

If `fimo` is not on PATH, pass its location explicitly with --fimo-path.

Optional, only needed for --bgzip-index (JBrowse 2 prep):
    conda install -c bioconda htslib     # provides bgzip + tabix
    # or: sudo apt-get install tabix

This script deliberately does NOT try to install MEME Suite for you -
compiling/pulling large bioinformatics toolchains inside an unattended
script is fragile across platforms, so please install it once using one
of the methods above.

------------------------------------------------------------------------
USAGE EXAMPLES
------------------------------------------------------------------------
# Simplest possible invocation - just the 3 input files. Everything else
# (TFBS-scan engine, JBrowse-2 indexing, output filename, JASPAR caching)
# is auto-detected/auto-derived, so this alone produces the fullest output
# your machine is capable of (--fasta/--gff/--te-file are short aliases
# for --genome-fasta/--annotation-gff/--te-hits):
python build_tfbs_te_gff.py \\
    --te-file D_suzukiiGenesAffectedByTE.txt \\
    --gff dsuzukii_annotation.gff3 \\
    --fasta dsuzukii_genome.fa

# Full pipeline, spelled out with the long flag names + explicit output:
python build_tfbs_te_gff.py \\
    --te-hits D_suzukiiGenesAffectedByTE.txt \\
    --annotation-gff dsuzukii_annotation.gff3 \\
    --genome-fasta dsuzukii_genome.fa \\
    --output combined_cyp_annotation.gff3

# No local genome FASTA on hand - fetch sequence windows from NCBI
# by accession/coordinates instead (seqids must be NCBI accessions):
python build_tfbs_te_gff.py \\
    --te-hits D_suzukiiGenesAffectedByTE.txt \\
    --annotation-gff dsuzukii_annotation.gff3 \\
    --sequence-source ncbi \\
    --output combined_cyp_annotation.gff3

# Skip TFBS scanning entirely and just merge genes/exons/TEs:
python build_tfbs_te_gff.py \\
    --te-hits D_suzukiiGenesAffectedByTE.txt \\
    --annotation-gff dsuzukii_annotation.gff3 \\
    --skip-tfbs \\
    --output genes_exons_te_only.gff3

# Running the SAME pipeline for a different species (same 3 input files -
# TE-hits table / annotation GFF3 / genome FASTA - just for that species):
# put its settings in a [section] of a config file instead of retyping
# every flag. See species_config.example.ini for a filled-in template.
python build_tfbs_te_gff.py --config species_config.ini --species suzukii
python build_tfbs_te_gff.py --config species_config.ini --species melanogaster
# (--species can be omitted if the config file only defines one section)
# Any flag given on the command line still overrides that species' config,
# e.g. to bump the FIMO threshold for just this run:
python build_tfbs_te_gff.py --config species_config.ini --species suzukii --fimo-thresh 1e-5
"""

import argparse
import configparser
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
import urllib.error
from collections import OrderedDict, defaultdict
from pathlib import Path
from types import SimpleNamespace

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

JASPAR_RELEASE_DEFAULT = "2026"
JASPAR_MEME_URL_TEMPLATE = (
    "https://jaspar.elixir.no/download/data/{release}/CORE/"
    "JASPAR{release}_CORE_insects_non-redundant_pfms_meme.txt"
)

NCBI_EFETCH_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    "&seq_start={start}&seq_stop={stop}"
)
NCBI_RATE_LIMIT_SECONDS = 0.35  # be polite to NCBI (~3 req/sec without an API key)

GENE_SYMBOL_ATTR_KEYS = ("ID", "gene", "gene_name", "Name", "locus_tag")

DEFAULT_UPSTREAM = 1000
DEFAULT_DOWNSTREAM = 200
DEFAULT_TE_FLANK = 50
DEFAULT_FIMO_THRESH = 1e-4
FIMO_DOCKER_IMAGE_DEFAULT = "memesuite/memesuite"
DOCKER_HTSLIB_IMAGE_DEFAULT = "quay.io/biocontainers/htslib:1.19.1--h81da01d_2"


# --------------------------------------------------------------------------
# Dependency checks (no auto-install - see setup instructions in docstring)
# --------------------------------------------------------------------------

def check_fimo(fimo_path):
    """Return (ok: bool, version_or_error: str)."""
    exe = shutil.which(fimo_path) or (fimo_path if Path(fimo_path).exists() else None)
    if exe is None:
        return False, f"'{fimo_path}' not found on PATH."
    try:
        out = subprocess.run([fimo_path, "--version"], capture_output=True, text=True, timeout=15)
        version = (out.stdout or out.stderr).strip()
        return True, version or "fimo found (version unknown)"
    except Exception as exc:  # noqa: BLE001
        return False, f"'{fimo_path}' found but failed to run: {exc}"


def check_docker():
    """Return (ok: bool, version_or_error: str) - used when --fimo-via-docker is set."""
    if shutil.which("docker") is None:
        return False, "'docker' not found on PATH (is Docker Desktop installed and running?)."
    try:
        out = subprocess.run(["docker", "--version"], capture_output=True, text=True, timeout=15)
        if out.returncode != 0:
            return False, (out.stderr or out.stdout).strip() or "docker found but returned an error."
        version = (out.stdout or out.stderr).strip()
        return True, version or "docker found (version unknown)"
    except Exception as exc:  # noqa: BLE001
        return False, f"'docker' found but failed to run: {exc}"


SETUP_HINT = """
Missing dependency detected. See the setup instructions at the top of
this script (`python build_tfbs_te_gff.py --help` also prints a summary),
or re-run with --skip-tfbs to build the gene/exon/TE-only GFF3 without
motif scanning.
""".strip()


def resolve_fimo_engine(args):
    """Decide how (or whether) to run fimo, so a plain `--te-hits/
    --annotation-gff/--genome-fasta` invocation with no other flags still
    gets the fullest output this machine can produce.

    Returns (engine, message) where engine is "docker", "fimo", or None
    (None means: skip TFBS scanning, message explains why).
      - args.fimo_via_docker is True  -> Docker required, no fallback.
      - args.fimo_via_docker is False -> local fimo required, no fallback
        (set via --no-fimo-docker).
      - args.fimo_via_docker is None  -> auto: try Docker first (most
        likely to "just work", especially on Windows), then a local fimo
        binary, then give up gracefully.
    """
    if args.fimo_via_docker is True:
        ok, msg = check_docker()
        return ("docker", msg) if ok else (None, msg)
    if args.fimo_via_docker is False:
        ok, msg = check_fimo(args.fimo_path)
        return ("fimo", msg) if ok else (None, msg)

    docker_ok, docker_msg = check_docker()
    if docker_ok:
        return "docker", docker_msg
    fimo_ok, fimo_msg = check_fimo(args.fimo_path)
    if fimo_ok:
        return "fimo", fimo_msg
    return None, f"no fimo engine available (Docker: {docker_msg}; local fimo: {fimo_msg})"


# --------------------------------------------------------------------------
# 1. TE hits table parsing
#    Format: <seqid> TAB <gene> TAB <RepeatMasker .out-style fields>
#    e.g.:
#    NW_023496800.1\tCyp4e3\t  5025    4.0  4.7  0.0  NW_023496800.1   \
#        6720897  6721034 (18868207) C rnd-1_family-204  RC/Helentron  \
#        (901)    672     535  7363
# --------------------------------------------------------------------------

def _strip_parens(token):
    """RepeatMasker .out files wrap some numeric fields in parentheses
    (e.g. "(18868207)" for the number of bases left in the TE consensus
    after the match) and sometimes suffix a '*' flag. Strip both so the
    caller gets the bare value."""
    return token.strip().lstrip("(").rstrip(")").rstrip("*")


def parse_te_hits(path):
    """Parse a TE-hits table into a de-duplicated list of TE record dicts.

    Each input line is <seqid> TAB <gene symbol> TAB <RepeatMasker .out
    fields, whitespace-separated>. The RepeatMasker fields (in order) are:
    Smith-Waterman score, %div, %del, %ins, query seqid (ignored - we
    already have seqid from column 1), query start, query end, bases left
    in query (ignored), strand ("C" for complement/minus, otherwise "+"),
    matching repeat name, repeat class/family, 3 more positional fields
    (str1/str2/str3, sometimes parenthesized), an optional RepeatMasker
    record ID, and an optional '*' flag meaning "overlaps a higher-scoring
    match". See https://www.repeatmasker.org/webrepeatmaskerhelp.html for
    the full .out format spec.

    Returns one dict per parsed line, later de-duplicated (this input
    format commonly repeats identical rows 2-3x)."""
    records = []
    with open(path, "r", newline="") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                # Not tab-delimited the way we expect - skip with a warning.
                print(f"[warn] te-hits line {lineno}: expected >=3 tab fields, got {len(parts)}; skipping", file=sys.stderr)
                continue
            seqid = parts[0].strip()
            gene = parts[1].strip()
            rest = "\t".join(parts[2:]).split()
            if len(rest) < 14:
                print(f"[warn] te-hits line {lineno}: only {len(rest)} RepeatMasker fields, expected >=14; skipping", file=sys.stderr)
                continue
            try:
                sw_score = float(rest[0])
                pct_div = float(rest[1])
                pct_del = float(rest[2])
                pct_ins = float(rest[3])
                q_start = int(rest[5])
                q_end = int(rest[6])
            except ValueError:
                print(f"[warn] te-hits line {lineno}: could not parse numeric fields; skipping", file=sys.stderr)
                continue
            strand_raw = rest[8]
            strand = "-" if strand_raw.upper() == "C" else "+"
            repeat_name = rest[9]
            repeat_class = rest[10]
            rep_field_a = _strip_parens(rest[11])
            rep_field_b = _strip_parens(rest[12])
            rep_field_c = _strip_parens(rest[13])
            rec_id = rest[14] if len(rest) > 14 else f"L{lineno}"
            overlap_flag = bool(len(rest) > 15 and rest[15].strip() == "*")

            start, end = (q_start, q_end) if q_start <= q_end else (q_end, q_start)

            records.append({
                "seqid": seqid,
                "gene": gene,
                "sw_score": sw_score,
                "pct_div": pct_div,
                "pct_del": pct_del,
                "pct_ins": pct_ins,
                "start": start,
                "end": end,
                "strand": strand,
                "repeat_name": repeat_name,
                "repeat_class": repeat_class,
                "rep_pos": (rep_field_a, rep_field_b, rep_field_c),
                "rm_id": rec_id,
                "overlaps_other": overlap_flag,
            })

    # De-duplicate: this input format commonly repeats identical rows.
    seen = OrderedDict()
    for r in records:
        key = (r["seqid"], r["gene"], r["start"], r["end"], r["repeat_name"], r["strand"])
        seen.setdefault(key, r)
    deduped = list(seen.values())
    print(f"[info] parsed {len(records)} TE-hit rows -> {len(deduped)} unique after de-duplication")
    return deduped


# --------------------------------------------------------------------------
# 2. Annotation GFF3 parsing (genes / mRNAs / exons for genes of interest)
# --------------------------------------------------------------------------

def _parse_gff3_attributes(attr_field):
    """Parse a GFF3 column-9 attributes string into a dict.

    GFF3 attributes are semicolon-separated `Key=value` pairs, e.g.
    "ID=gene1;Name=Cyp4e3;gene_biotype=protein_coding". Reserved keys
    (ID, Name, Parent, ...) are capitalized by convention but this parser
    doesn't enforce that - it just splits on '=' and ';' and keeps
    whatever key casing the source file used."""
    attrs = {}
    for chunk in attr_field.strip().split(";"):
        chunk = chunk.strip()
        if not chunk or "=" not in chunk:
            continue
        key, _, value = chunk.partition("=")
        attrs[key.strip()] = value.strip()
    return attrs


def _gff3_line_iter(path):
    # Some real-world GFF3 exports (e.g. ones that have been text-processed)
    # contain literal tab characters inside free-text attribute values
    # (model_evidence notes, Dbxref lists, etc.), which breaks a naive
    # split("\t"). GFF3's spec guarantees exactly 9 columns with attributes
    # as the last one, so split with maxsplit=8 and treat everything after
    # the 8th tab as the (possibly tab-containing) attributes field. Also
    # strip \r for files with CRLF line endings.
    with open(path, "r", newline="") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t", 8)
            if len(fields) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attr_field = fields
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            yield {
                "seqid": seqid,
                "source": source,
                "type": ftype,
                "start": start_i,
                "end": end_i,
                "score": score,
                "strand": strand,
                "phase": phase,
                "attrs": _parse_gff3_attributes(attr_field),
            }


def _attr_tokens(attrs, keys):
    """Yield every comma-separated token across the given attribute keys.

    Some annotations (e.g. merged/duplicated tandem gene models) pack
    multiple gene symbols into one attribute, comma-separated, such as
    ID=Cyp4e2,Cyp4e3,Cyp4e1 for a single locus record shared by several
    paralogs. Child mRNA/exon features then reuse that same comma-joined
    string as their Parent, so callers should match/register token by
    token rather than as one opaque string.
    """
    for key in keys:
        val = attrs.get(key)
        if not val:
            continue
        for token in val.split(","):
            token = token.strip()
            if token:
                yield token


def _symbol_matches(attrs, target_symbols_lower):
    """Return the list of raw tokens (across ID/gene/gene_name/Name/locus_tag)
    that match one of our target gene symbols, case-insensitively."""
    return [t for t in _attr_tokens(attrs, GENE_SYMBOL_ATTR_KEYS) if t.lower() in target_symbols_lower]


def parse_annotation_gff3(path, gene_symbols):
    """
    Three lightweight passes over the annotation GFF3:
      pass 1: find `gene` features matching our target Cyp gene symbols
      pass 2: find `mRNA` features whose Parent is one of those genes
      pass 3: find `exon` features whose Parent is one of those mRNAs, then
              pick one representative transcript per gene and derive
              `intron` features from the gaps between its exons
    Returns: dict gene_symbol -> {"gene": rec, "mrna": rec or None,
                                   "exons": [...], "introns": [...]}
    """
    target_lower = {g.lower() for g in gene_symbols}
    result = {g: {"gene": None, "mrna": None, "exons": [], "introns": []} for g in gene_symbols}

    # lowercase -> canonical-cased symbol, for fast lookup
    symbol_by_lower = {g.lower(): g for g in gene_symbols}

    gene_id_to_symbol = {}
    print("[info] annotation pass 1/3: scanning for gene features...")
    for rec in _gff3_line_iter(path):
        if rec["type"] != "gene":
            continue
        matched_tokens = _symbol_matches(rec["attrs"], target_lower)
        if not matched_tokens:
            continue
        matched_symbols = {symbol_by_lower[t.lower()] for t in matched_tokens}
        for canonical in matched_symbols:
            if result[canonical]["gene"] is None:
                result[canonical]["gene"] = rec
        # Register every comma-separated ID token that matches one of our
        # target genes, so downstream mRNA/exon Parent references (which
        # reuse this same raw ID string, split the same way) resolve
        # correctly - including merged multi-paralog loci.
        gid_raw = rec["attrs"].get("ID", "")
        for token in gid_raw.split(","):
            token = token.strip()
            if token.lower() in target_lower:
                gene_id_to_symbol[token] = symbol_by_lower[token.lower()]
        # Fallback: some annotations use a single opaque ID (not a gene
        # symbol) with the symbol only in Name/gene/gene_name - register
        # the whole raw ID pointing at whichever symbol(s) matched via
        # those other attributes, so Parent linkage still works.
        if gid_raw and gid_raw not in gene_id_to_symbol and gid_raw.lower() not in target_lower:
            gene_id_to_symbol[gid_raw] = next(iter(matched_symbols))

    found_genes = [g for g, v in result.items() if v["gene"]]
    print(f"[info] matched {len(found_genes)}/{len(gene_symbols)} gene symbols in annotation")

    if not gene_id_to_symbol:
        return result

    print("[info] annotation pass 2/3: scanning for mRNA features...")
    # mRNA ID -> set of symbols. A plain 1:1 dict would lose data when one
    # mRNA is shared by several merged paralogs (Parent=Cyp4e2,Cyp4e3,Cyp4e1
    # all resolving to the same physical transcript record) - every gene
    # symbol sharing that locus should also share its exons downstream.
    mrna_id_to_symbols = defaultdict(set)
    mrna_records = {}                    # mRNA ID -> record
    symbol_mrna_ids = defaultdict(list)  # symbol -> [mRNA ID, ...]
    for rec in _gff3_line_iter(path):
        if rec["type"] not in ("mRNA", "transcript", "primary_transcript"):
            continue
        parent = rec["attrs"].get("Parent", "")
        matched_symbols_here = set()
        for pid in parent.split(","):
            pid = pid.strip()
            if pid in gene_id_to_symbol:
                matched_symbols_here.add(gene_id_to_symbol[pid])
        if not matched_symbols_here:
            continue
        mid = rec["attrs"].get("ID")
        if not mid:
            continue
        mrna_records[mid] = rec
        mrna_id_to_symbols[mid] |= matched_symbols_here
        for symbol in matched_symbols_here:
            symbol_mrna_ids[symbol].append(mid)

    if not mrna_id_to_symbols:
        print("[warn] no mRNA/transcript features found under matched genes; exons/introns will be empty")
        return result

    print("[info] annotation pass 3/3: scanning for exon features...")
    mrna_exons = defaultdict(list)  # mRNA ID -> [exon rec, ...]
    for rec in _gff3_line_iter(path):
        if rec["type"] != "exon":
            continue
        parent = rec["attrs"].get("Parent", "")
        for pid in parent.split(","):
            pid = pid.strip()
            if pid in mrna_id_to_symbols:
                mrna_exons[pid].append(rec)

    # Pick one representative transcript per gene (most exons, tie-broken by
    # longest span) and derive intron features from the gaps between its
    # sorted, de-duplicated exons. A single flattened representative model
    # keeps the gene/exon/intron/strand structure unambiguous for display,
    # rather than merging exons across every isoform.
    for symbol in gene_symbols:
        mids = symbol_mrna_ids.get(symbol) or []
        if not mids:
            continue

        def _mrna_key(mid):
            n_exons = len(mrna_exons.get(mid, []))
            span = mrna_records[mid]["end"] - mrna_records[mid]["start"]
            return (n_exons, span)

        best_mid = max(mids, key=_mrna_key)
        result[symbol]["mrna"] = mrna_records[best_mid]

        seen = set()
        exons = []
        for e in sorted(mrna_exons.get(best_mid, []), key=lambda e: e["start"]):
            key = (e["seqid"], e["start"], e["end"], e["strand"])
            if key in seen:
                continue
            seen.add(key)
            exons.append(e)
        result[symbol]["exons"] = exons

        introns = []
        for i in range(len(exons) - 1):
            gap_start = exons[i]["end"] + 1
            gap_end = exons[i + 1]["start"] - 1
            if gap_end >= gap_start:
                introns.append({
                    "seqid": exons[i]["seqid"],
                    "start": gap_start,
                    "end": gap_end,
                    "strand": exons[i]["strand"],
                    "source": exons[i]["source"],
                    "score": ".",
                    "phase": ".",
                })
        result[symbol]["introns"] = introns

    n_with_transcript = sum(1 for g in gene_symbols if result[g]["mrna"])
    total_exons = sum(len(v["exons"]) for v in result.values())
    total_introns = sum(len(v["introns"]) for v in result.values())
    print(f"[info] selected representative transcripts for {n_with_transcript}/{len(gene_symbols)} genes: "
          f"{total_exons} exons, {total_introns} introns")
    return result


# --------------------------------------------------------------------------
# 3. Region computation (promoter windows + TE-flank windows)
# --------------------------------------------------------------------------

def compute_promoter_region(gene_rec, upstream, downstream):
    """Build a genome window around a gene's transcription start site (TSS)
    to scan for promoter-region TFBS motifs.

    GFF3 coordinates are always given left-to-right on the + strand
    regardless of gene orientation, so "the start of the gene record" is
    only the TSS for + strand genes. For - strand genes the TSS is
    biologically at the *right-hand* (higher-coordinate) end of the
    feature, i.e. gene_rec["end"], and "upstream" means increasing
    coordinates instead of decreasing ones. All coordinates here are
    1-based inclusive, matching GFF3's own convention; start is clamped
    to a minimum of 1 since there's no such thing as base 0 or negative
    coordinates in a GFF3/FASTA file."""
    if gene_rec["strand"] == "-":
        tss = gene_rec["end"]
        start = max(1, tss - downstream)
        end = tss + upstream
    else:
        tss = gene_rec["start"]
        start = max(1, tss - upstream)
        end = tss + downstream
    return {"seqid": gene_rec["seqid"], "start": start, "end": end, "strand": gene_rec["strand"]}


def compute_te_region(te_rec, flank):
    """Pad a TE hit's coordinates by `flank` bp on each side, so FIMO also
    sees a bit of surrounding sequence context rather than exactly the
    RepeatMasker-called boundaries (motifs can straddle the edge of a
    TE call). Same 1-based, clamped-to-1 coordinate convention as
    compute_promoter_region."""
    start = max(1, te_rec["start"] - flank)
    end = te_rec["end"] + flank
    return {"seqid": te_rec["seqid"], "start": start, "end": end, "strand": te_rec["strand"]}


def build_scan_regions(genes_info, te_records, upstream, downstream, te_flank):
    """
    Build a de-duplicated list of scan windows, each tagged with the gene(s)
    and origin (promoter / te) it belongs to, for sequence extraction.
    Returns: list of region dicts with a unique 'region_id'.
    """
    regions = OrderedDict()  # (seqid,start,end) -> region dict

    for gene, info in genes_info.items():
        if not info["gene"]:
            continue
        promo = compute_promoter_region(info["gene"], upstream, downstream)
        key = (promo["seqid"], promo["start"], promo["end"])
        if key not in regions:
            regions[key] = {**promo, "origin": "promoter", "genes": set(), "region_id": None}
        regions[key]["genes"].add(gene)

    for te in te_records:
        win = compute_te_region(te, te_flank)
        key = (win["seqid"], win["start"], win["end"])
        if key not in regions:
            regions[key] = {**win, "origin": "te_flank", "genes": set(), "region_id": None}
        regions[key]["genes"].add(te["gene"])

    out = []
    for i, (key, region) in enumerate(regions.items(), 1):
        region["region_id"] = f"region_{i}"
        out.append(region)
    print(f"[info] built {len(out)} sequence-scan windows ({sum(1 for r in out if r['origin']=='promoter')} promoter, "
          f"{sum(1 for r in out if r['origin']=='te_flank')} TE-flank)")
    return out


# --------------------------------------------------------------------------
# 4a. Sequence extraction - local FASTA (built-in indexer) or NCBI fallback
#
# A minimal, dependency-free re-implementation of samtools-style faidx
# random access: index each sequence's byte offset + line geometry once,
# then seek()+read() only the bytes needed for a given coordinate range.
# Validated byte-for-byte against a real 140MB repeat-masked D. melanogaster
# genome FASTA (GCF_000001215.4) with plain ">accession" headers.
# --------------------------------------------------------------------------

def build_fasta_index(fasta_path):
    """Scan a genome FASTA once and build an in-memory index for random
    access, the same idea as samtools' .fai index file (but kept in memory
    here instead of written to disk).

    A FASTA file is one or more records, each a ">seqid ..." header line
    followed by the sequence wrapped to a fixed line width (commonly 60,
    70, or 80 bases/line). Because every wrapped line for a given record
    has the *same* number of bases (except the last, shorter one), you can
    compute exactly which byte offset any base lives at from just 4
    numbers per sequence - no need to scan the whole multi-GB file for
    every region lookup. That's what this function precomputes and what
    fetch_fasta_region() then uses to seek() straight to the bytes it needs.

    Returns {seqid: (seq_start_byte_offset, seq_len_bases, line_bases, line_bytes)}
    where seq_start_byte_offset is the byte position right after the
    header's newline, seq_len_bases is the total ungapped sequence length,
    line_bases is how many bases are on each wrapped line, and line_bytes
    is how many bytes that line occupies on disk (line_bases plus 1-2
    bytes for the line ending, \\n or \\r\\n)."""
    index = {}
    with open(fasta_path, "rb") as fh:
        seqid = None
        seq_start_offset = None
        line_bases = None
        line_bytes = None
        seq_len = 0
        offset = 0
        for line in fh:
            line_len = len(line)
            if line.startswith(b">"):
                if seqid is not None:
                    index[seqid] = (seq_start_offset, seq_len, line_bases, line_bytes)
                seqid = line[1:].split()[0].decode()
                seq_len = 0
                seq_start_offset = offset + line_len
                line_bases = None
                line_bytes = None
            else:
                stripped = line.rstrip(b"\r\n")
                if line_bases is None and len(stripped) > 0:
                    line_bases = len(stripped)
                    line_bytes = line_len
                seq_len += len(stripped)
            offset += line_len
        if seqid is not None:
            index[seqid] = (seq_start_offset, seq_len, line_bases, line_bytes)
    return index


def fetch_fasta_region(fasta_path, index, seqid, start, end):
    """Fetch a 1-based inclusive [start, end] region for seqid using the index.

    Converts the requested base range into a byte offset (accounting for
    how many newline characters fall within wrapped lines before that
    point), seeks directly there, reads a generous over-estimate of bytes
    to cover any newlines in the read span, strips those newlines out, and
    trims to exactly the requested number of bases. This avoids reading
    the file sequentially from the start for every single region."""
    if seqid not in index:
        return None
    seq_start_offset, seq_len, line_bases, line_bytes = index[seqid]
    if line_bases is None or seq_len == 0:
        return ""
    start = max(1, start)
    end = min(seq_len, end)
    if start > end:
        return ""
    start0 = start - 1
    start_line = start0 // line_bases
    start_line_pos = start0 % line_bases
    byte_offset = seq_start_offset + start_line * line_bytes + start_line_pos
    n_bases = end - start0
    # Read generously (accounting for newline bytes on wrapped lines), then trim.
    length_to_read = n_bases + (n_bases // line_bases + 2) * (line_bytes - line_bases)
    with open(fasta_path, "rb") as fh:
        fh.seek(byte_offset)
        raw = fh.read(length_to_read)
    seq = raw.replace(b"\n", b"").replace(b"\r", b"")[:n_bases]
    return seq.decode()


def extract_regions_local_fasta(genome_fasta, regions):
    """Index the genome FASTA once, then pull the sequence for every scan
    window (promoter + TE-flank regions from build_scan_regions()) out of
    it. Returns {region_id: {"seq": str, "seqid": ..., "start": ..., "end": ...}},
    keyed the same way regions are, so remap_fimo_hits() can later map a
    FIMO hit's region-local coordinates back to genome coordinates."""
    print(f"[info] indexing genome FASTA {genome_fasta} ...")
    index = build_fasta_index(genome_fasta)
    print(f"[info] indexed {len(index)} sequences")
    seqs = {}
    missing = []
    for r in regions:
        seqid = r["seqid"]
        if seqid not in index:
            missing.append(seqid)
            continue
        seq_len = index[seqid][1]
        start = max(1, r["start"])
        end = min(seq_len, r["end"])
        if start > end:
            continue
        seq = fetch_fasta_region(genome_fasta, index, seqid, start, end)
        if seq:
            seqs[r["region_id"]] = {"seq": seq, "seqid": seqid, "start": start, "end": end}
    if missing:
        print(f"[warn] {len(set(missing))} seqid(s) from the TE/annotation files were not found in the genome FASTA "
              f"(e.g. {sorted(set(missing))[:5]})", file=sys.stderr)
    return seqs


def extract_regions_ncbi(regions, email=None, api_key=None):
    """Fetch each region's sequence directly from NCBI by accession + coordinates.
    Only works when region['seqid'] is a resolvable NCBI nucleotide accession."""
    seqs = {}
    for r in regions:
        url = NCBI_EFETCH_URL.format(accession=r["seqid"], start=r["start"], stop=r["end"])
        if api_key:
            url += f"&api_key={api_key}"
        if email:
            url += f"&email={email}"
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                text = resp.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, urllib.error.HTTPError) as exc:
            print(f"[warn] NCBI fetch failed for {r['seqid']}:{r['start']}-{r['end']}: {exc}", file=sys.stderr)
            time.sleep(NCBI_RATE_LIMIT_SECONDS)
            continue
        lines = text.splitlines()
        seq = "".join(l.strip() for l in lines if l and not l.startswith(">"))
        if seq:
            seqs[r["region_id"]] = {"seq": seq, "seqid": r["seqid"], "start": r["start"], "end": r["end"]}
        time.sleep(NCBI_RATE_LIMIT_SECONDS)
    return seqs


def write_regions_fasta(seqs, out_fasta_path):
    """Write the extracted scan-window sequences out as a FASTA file (one
    record per region, headered by its region_id) - this is the file
    that actually gets handed to `fimo` as its sequence input. Lines are
    wrapped to 70 bases, the conventional FASTA line width."""
    with open(out_fasta_path, "w") as fh:
        for region_id, info in seqs.items():
            fh.write(f">{region_id}\n")
            seq = info["seq"]
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


# --------------------------------------------------------------------------
# 4b. JASPAR motif retrieval
# --------------------------------------------------------------------------

def get_jaspar_motif_file(local_path, download_dir, release=JASPAR_RELEASE_DEFAULT, url_override=None):
    """Resolve the JASPAR motif file to use, downloading it only if needed.

    download_dir is meant to be a *persistent* cache directory (not a
    temp work_dir that gets deleted) - the same JASPAR release is valid
    for every species, so downloading it once and reusing it on every
    later run/species avoids repeat network calls (and, e.g., any local
    SSL/proxy trouble) entirely after the first successful download.
    """
    if local_path:
        p = Path(local_path)
        if not p.exists():
            raise FileNotFoundError(f"--jaspar-motif-file given but not found: {local_path}")
        return str(p)

    dest = Path(download_dir) / f"JASPAR{release}_CORE_insects_non-redundant_pfms_meme.txt"
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[info] reusing cached JASPAR motif file: {dest}")
        return str(dest)

    url = url_override or JASPAR_MEME_URL_TEMPLATE.format(release=release)
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[info] downloading JASPAR insect motif collection from:\n       {url}")
    try:
        urllib.request.urlretrieve(url, dest)
    except (urllib.error.URLError, urllib.error.HTTPError) as exc:
        raise RuntimeError(
            f"Could not download JASPAR motifs automatically ({exc}).\n"
            "Download manually from https://jaspar.elixir.no/downloads/ "
            "(CORE PFMs -> Insects -> MEME, non-redundant, single batch file) "
            "and pass the path via --jaspar-motif-file (or drop it straight into the "
            f"cache dir as {dest})."
        ) from exc
    print(f"[info] saved JASPAR motif file to {dest} (will be reused on future runs)")
    return str(dest)


# --------------------------------------------------------------------------
# 4b-bis. Cap-n-Collar (CncC) / Maf-S antioxidant/xenobiotic response
# element (ARE/MARE) - not in JASPAR's insect collection at all (checked:
# the 2026 CORE insects non-redundant release has 129 motifs, none of them
# CncC, Nrf2, cnc, or Maf-S), even though CncC:Maf-S is THE master
# regulator of CYP-mediated xenobiotic/insecticide detoxification in
# insects - directly relevant to TE-affected Cyp genes. So instead of
# relying on JASPAR, we build a small literature-derived PFM for its
# binding site from the published consensus and splice it into the same
# motif file FIMO scans, using the exact same promoter/TE windows already
# built for the JASPAR motifs - no separate pipeline needed.
#
# Consensus: 5'-TMAnnRTGAYnnnGCRwwww-3' (IUPAC), the CncC:Maf/Nrf2:Maf
# antioxidant response element (ARE) / Maf recognition element (MARE).
# First identified for Drosophila Cnc:Maf-S by Veraksa et al. 2000, and
# shown to be necessary and sufficient for xenobiotic-inducible Cyp6a2
# transcription by Misra JR, Horner MA, Lam G, Thummel CS. "Transcriptional
# regulation of xenobiotic detoxification in Drosophila." Genes Dev. 2011
# Sep 1;25(17):1796-806. doi: 10.1101/gad.17280911.
#
# This is a hand-built, consensus-based PFM (defined IUPAC positions get
# 0.85 probability, ambiguous positions split evenly among their allowed
# bases, N positions are uniform) - NOT an empirically-fit PWM from
# ChIP-seq/SELEX data the way the JASPAR motifs are. The real element is
# also known in the literature to tolerate substantial sequence
# variability, so treat hits as candidate/putative CncC:Maf-S sites
# worth follow-up validation, same as any consensus-based motif scan.
# --------------------------------------------------------------------------

CNCC_MAF_ARE_CONSENSUS = "TMAnnRTGAYnnnGCRwwww"
CNCC_MAF_ARE_MOTIF_ID = "CncC_Maf_ARE"
CNCC_MAF_ARE_MOTIF_NAME = "CncC::Maf-S_ARE_consensus_Veraksa2000_Misra2011"

_CNCC_MAF_ARE_IUPAC = {
    "A": (0.85, 0.05, 0.05, 0.05), "C": (0.05, 0.85, 0.05, 0.05),
    "G": (0.05, 0.05, 0.85, 0.05), "T": (0.05, 0.05, 0.05, 0.85),
    "M": (0.45, 0.45, 0.05, 0.05),  # A or C
    "R": (0.45, 0.05, 0.45, 0.05),  # A or G
    "Y": (0.05, 0.45, 0.05, 0.45),  # C or T
    "W": (0.45, 0.05, 0.05, 0.45),  # A or T
    "N": (0.25, 0.25, 0.25, 0.25),
}


def build_cncc_maf_are_meme_block():
    """Return a MEME-format MOTIF block (no header - meant to be appended
    to an existing MEME file, e.g. the downloaded JASPAR one) for the
    CncC:Maf-S ARE/MARE consensus. See module comment above for sourcing."""
    rows = []
    for base in CNCC_MAF_ARE_CONSENSUS.upper():
        a, c, g, t = _CNCC_MAF_ARE_IUPAC[base]
        rows.append(f" {a:.6f}  {c:.6f}  {g:.6f}  {t:.6f}")
    header = (
        f"MOTIF {CNCC_MAF_ARE_MOTIF_ID} {CNCC_MAF_ARE_MOTIF_NAME}\n"
        f"letter-probability matrix: alength= 4 w= {len(CNCC_MAF_ARE_CONSENSUS)} "
        f"nsites= 20 E= 0"
    )
    return header + "\n" + "\n".join(rows) + "\n"


def add_cncc_maf_are_motif(motif_file, work_dir):
    """Append the CncC:Maf-S ARE motif block to motif_file, writing the
    result into work_dir so the original (possibly cached) JASPAR file is
    never modified in place. Returns the path FIMO should actually scan."""
    combined_path = Path(work_dir) / "motifs_with_cncc_maf_are.meme"
    base_text = Path(motif_file).read_text()
    with open(combined_path, "w") as fh:
        fh.write(base_text)
        if not base_text.endswith("\n"):
            fh.write("\n")
        fh.write("\n")
        fh.write(build_cncc_maf_are_meme_block())
    print(f"[info] added CncC:Maf-S ARE consensus motif ({CNCC_MAF_ARE_MOTIF_ID}) to the motif set to scan")
    return str(combined_path)


# --------------------------------------------------------------------------
# 4c. Run FIMO and parse results
#
# FIMO ("Find Individual Motif Occurrences", part of the MEME Suite) scans
# a set of sequences against a set of motifs and reports every position
# where a motif matches well enough to beat the given p-value threshold.
# A "motif" here is a position weight matrix (PWM/PFM) - a per-position
# probability of A/C/G/T - not a fixed string, so a single motif matches a
# whole family of similar-but-not-identical sequences (this is how real
# transcription factors bind DNA: loosely, not to one exact sequence).
# --------------------------------------------------------------------------

def run_fimo(fimo_path, motif_file, fasta_file, out_dir, thresh):
    """Run the local `fimo` binary against motif_file (MEME-format motifs,
    e.g. the JASPAR download plus our CncC:Maf-S addition) and fasta_file
    (the scan-window sequences from write_regions_fasta()). FIMO writes
    several output files into out_dir (--oc = "output directory"); we only
    need fimo.tsv, its tab-separated hit table. Raises if fimo exits
    non-zero or doesn't produce that file."""
    cmd = [
        fimo_path,
        "--oc", str(out_dir),
        "--thresh", str(thresh),
        "--verbosity", "1",
        motif_file,
        fasta_file,
    ]
    print(f"[info] running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"fimo failed (exit {result.returncode}):\n{result.stderr}")
    tsv_path = Path(out_dir) / "fimo.tsv"
    if not tsv_path.exists():
        raise RuntimeError(f"fimo ran but no fimo.tsv found in {out_dir}")
    return tsv_path


def run_fimo_docker(image, motif_file, fasta_file, out_dir, thresh, work_dir):
    """Run fimo inside the official MEME Suite Docker container instead of a
    local binary - handy on Windows, where compiling MEME Suite natively is
    painful (MEME Suite's own docs recommend Docker or WSL there instead).

    The official image (memesuite/memesuite) expects the working directory
    to be bind-mounted at /home/meme, so everything fimo needs to read
    (motif file, sequence FASTA) and write (--oc output dir) must live
    under a single mounted directory - here, the script's work_dir.
    """
    work_dir = Path(work_dir).resolve()
    motif_path = Path(motif_file).resolve()
    fasta_path = Path(fasta_file).resolve()
    out_dir = Path(out_dir)

    # Everything Docker can see has to be under the one bind-mounted folder.
    if work_dir not in motif_path.parents:
        local_motif = work_dir / motif_path.name
        shutil.copy2(motif_path, local_motif)
        motif_path = local_motif
    if work_dir not in fasta_path.parents:
        local_fasta = work_dir / fasta_path.name
        shutil.copy2(fasta_path, local_fasta)
        fasta_path = local_fasta

    mount_point = "/home/meme"
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{work_dir}:{mount_point}",
        image,
        "fimo",
        "--oc", f"{mount_point}/{out_dir.name}",
        "--thresh", str(thresh),
        "--verbosity", "1",
        f"{mount_point}/{motif_path.name}",
        f"{mount_point}/{fasta_path.name}",
    ]
    print(f"[info] running via Docker (first run may take a while to pull the image):\n       {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"docker fimo failed (exit {result.returncode}):\n{result.stderr}")
    tsv_path = out_dir / "fimo.tsv"
    if not tsv_path.exists():
        raise RuntimeError(f"fimo ran via Docker but no fimo.tsv found in {out_dir}")
    return tsv_path


def parse_fimo_tsv(tsv_path):
    """Parse FIMO's standard TSV output into a list of hit dicts."""
    hits = []
    with open(tsv_path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                continue
            if len(fields) != len(header):
                # FIMO appends trailing comment/blank lines - stop there.
                break
            row = dict(zip(header, fields))
            hits.append(row)
    return hits


def remap_fimo_hits(hits, region_lookup):
    """Convert FIMO's window-local coordinates back to absolute genome coords."""
    remapped = []
    for h in hits:
        region_id = h.get("sequence_name")
        region = region_lookup.get(region_id)
        if region is None:
            continue
        try:
            local_start = int(h["start"])
            local_stop = int(h["stop"])
        except (KeyError, ValueError):
            continue
        genome_start = region["start"] + local_start - 1
        genome_end = region["start"] + local_stop - 1
        remapped.append({
            "seqid": region["seqid"],
            "start": genome_start,
            "end": genome_end,
            "strand": h.get("strand", "."),
            "motif_id": h.get("motif_id", "."),
            "motif_alt_id": h.get("motif_alt_id", "."),
            "score": h.get("score", "."),
            "pvalue": h.get("p-value", "."),
            "qvalue": h.get("q-value", "."),
            "matched_sequence": h.get("matched_sequence", "."),
            "genes": sorted(region.get("genes", [])),
        })
    return remapped


# --------------------------------------------------------------------------
# 5. Combine everything into one GFF3
# --------------------------------------------------------------------------

def _gff3_escape(value):
    """Percent-encode characters that are structurally significant in a
    GFF3 attributes field (';' separates key=value pairs, '=' separates
    key from value, ',' separates multi-value lists) so free-text values
    like gene descriptions can't be misparsed as extra attributes. Tabs
    are just replaced with a space since GFF3 is itself tab-delimited."""
    return str(value).replace(";", "%3B").replace("=", "%3D").replace(",", "%2C").replace("\t", " ")


# Additional descriptive fields to carry over from the source annotation's
# own gene record, beyond the ID/Name/gene_biotype already written. Real
# annotations vary in what they provide (Gnomon-predicted D. suzukii genes
# have Dbxref/gene_synonym but no description/cyt_map; RefSeq-curated
# D. melanogaster genes have all of them) and in attribute-key casing
# (lowercase "dbxref" vs GFF3's own reserved "Dbxref"), so every key is
# read case-insensitively and only written out if actually present.
def _extract_gene_annotation_attrs(attrs):
    """Return a list of 'Key=value' GFF3 attribute strings for whichever of
    Dbxref/Note(description)/locus_tag/gene_synonym/cyt_map exist in attrs.
    Multi-value fields (Dbxref, gene_synonym) are comma-joined per GFF3
    convention even when the source file space-separates them."""
    lower = {k.lower(): v for k, v in attrs.items()}
    parts = []

    dbxref = lower.get("dbxref")
    if dbxref:
        tokens = [t for t in re.split(r"[,\s]+", dbxref.strip()) if t]
        if tokens:
            parts.append(f"Dbxref={','.join(_gff3_escape(t) for t in tokens)}")

    description = lower.get("description")
    if description:
        parts.append(f"Note={_gff3_escape(description)}")

    locus_tag = lower.get("locus_tag")
    if locus_tag:
        parts.append(f"locus_tag={_gff3_escape(locus_tag)}")

    gene_synonym = lower.get("gene_synonym")
    if gene_synonym:
        tokens = [t for t in gene_synonym.strip().split() if t]
        if tokens:
            parts.append(f"gene_synonym={','.join(_gff3_escape(t) for t in tokens)}")

    cyt_map = lower.get("cyt_map")
    if cyt_map:
        parts.append(f"cyt_map={_gff3_escape(cyt_map)}")

    return parts


def gene_record_to_gff3(gene_rec, symbol):
    """Render a parsed `gene` feature (from parse_annotation_gff3) as one
    standard 9-column GFF3 line: seqid, source, type, start, end, score,
    strand, phase, attributes. `phase` is only meaningful for CDS features
    (it isn't used here) but the column must still be present - GFF3
    requires all 9 columns on every line, using "." for not-applicable."""
    attrs = gene_rec["attrs"]
    ident = attrs.get("ID", f"gene-{symbol}")
    attr_parts = [
        f"ID={ident}",
        f"Name={symbol}",
        f"gene_biotype={attrs.get('gene_biotype', 'protein_coding')}",
    ]
    attr_parts.extend(_extract_gene_annotation_attrs(attrs))
    fields = [
        gene_rec["seqid"], gene_rec["source"] or "annotation", "gene",
        str(gene_rec["start"]), str(gene_rec["end"]),
        gene_rec["score"] or ".", gene_rec["strand"] or ".", gene_rec["phase"] or ".",
        ";".join(attr_parts),
    ]
    return "\t".join(fields)


def mrna_record_to_gff3(mrna_rec, symbol, gene_id):
    """Returns (gff3_line, mrna_id) - mrna_id is used as the Parent for
    this transcript's exon/intron features."""
    attrs = mrna_rec["attrs"]
    ident = attrs.get("ID", f"mrna-{symbol}")
    fields = [
        mrna_rec["seqid"], mrna_rec["source"] or "annotation", "mRNA",
        str(mrna_rec["start"]), str(mrna_rec["end"]),
        mrna_rec["score"] or ".", mrna_rec["strand"] or ".", mrna_rec["phase"] or ".",
        f"ID={ident};Parent={gene_id};Name={symbol}-RA;gene={symbol}",
    ]
    return "\t".join(fields), ident


def exon_record_to_gff3(exon_rec, symbol, idx, parent_id):
    """Render one exon feature, linked to its transcript via Parent=
    (GFF3's mechanism for expressing the gene -> mRNA -> exon hierarchy
    that lets genome browsers draw introns as the gaps between exons)."""
    fields = [
        exon_rec["seqid"], exon_rec["source"] or "annotation", "exon",
        str(exon_rec["start"]), str(exon_rec["end"]),
        exon_rec["score"] or ".", exon_rec["strand"] or ".", exon_rec["phase"] or ".",
        f"ID=exon-{symbol}-{idx};Parent={parent_id};gene={symbol}",
    ]
    return "\t".join(fields)


def intron_record_to_gff3(intron_rec, symbol, idx, parent_id):
    """Render one intron feature (a gap between two consecutive exons of
    the representative transcript, computed in parse_annotation_gff3).
    Introns aren't a required GFF3 feature type, but writing them
    explicitly makes gene structure visually obvious in a browser without
    the viewer having to infer gaps from exon spacing itself."""
    fields = [
        intron_rec["seqid"], intron_rec.get("source") or "annotation", "intron",
        str(intron_rec["start"]), str(intron_rec["end"]),
        intron_rec.get("score") or ".", intron_rec.get("strand") or ".", intron_rec.get("phase") or ".",
        f"ID=intron-{symbol}-{idx};Parent={parent_id};gene={symbol}",
    ]
    return "\t".join(fields)


def te_record_to_gff3(te_rec, idx):
    # Matches the convention used for genuine annotated TEs in the real
    # D. melanogaster RefSeq GFF (type=mobile_genetic_element), extended
    # with extra attributes (repeat class, %div/%del/%ins) so the record
    # stays informative; strand is always populated (RepeatMasker gives it,
    # unlike the reference file's own injected TE-proximity rows).
    fields = [
        te_rec["seqid"], "RepeatMasker", "mobile_genetic_element",
        str(te_rec["start"]), str(te_rec["end"]),
        str(te_rec["sw_score"]), te_rec["strand"], ".",
        (
            f"ID=TE_{idx};Name={_gff3_escape(te_rec['repeat_name'])};"
            f"repeat_class={_gff3_escape(te_rec['repeat_class'])};"
            f"pct_div={te_rec['pct_div']};pct_del={te_rec['pct_del']};pct_ins={te_rec['pct_ins']};"
            f"Description=Within range of {te_rec['gene']}"
        ),
    ]
    return "\t".join(fields)


def tfbs_record_to_gff3(hit, idx):
    """Render one remapped FIMO hit (already in absolute genome
    coordinates, from remap_fimo_hits()) as a `TF_binding_site` GFF3
    feature. Hits from the literature-derived CncC:Maf-S motif are tagged
    with a different `source` column and Description text than genuine
    JASPAR hits, so the two are never confused when browsing the track -
    see the CncC/Maf-S module comment above for why that motif exists
    outside JASPAR at all."""
    genes_str = ",".join(hit["genes"]) if hit["genes"] else "NA"
    is_cncc_maf = hit["motif_id"] == CNCC_MAF_ARE_MOTIF_ID
    if is_cncc_maf:
        # Not a real JASPAR accession - flag it distinctly so it's never
        # mistaken for an empirically-derived JASPAR hit (see the
        # add_cncc_maf_are_motif()/module comment above for sourcing).
        source = "FIMO/literature-ARE"
        matrix_attr_key = "motif_source_id"
        desc = (
            f"Predicted CncC:Maf-S antioxidant/xenobiotic response element (ARE/MARE) "
            f"near {genes_str} - literature consensus motif (Veraksa 2000; Misra et al. "
            f"2011 Genes Dev), not an empirical JASPAR PWM; candidate site, not "
            f"experimentally validated at this locus"
        )
    else:
        source = "FIMO/JASPAR"
        matrix_attr_key = "jaspar_matrix_id"
        desc = (
            f"Predicted {_gff3_escape(hit['motif_alt_id'])} binding site "
            f"(JASPAR {hit['motif_id']}) near {genes_str}"
        )
    fields = [
        hit["seqid"], source, "TF_binding_site",
        str(hit["start"]), str(hit["end"]),
        hit["score"], hit["strand"], ".",
        (
            f"ID=TFBS_{idx};Name={_gff3_escape(hit['motif_alt_id'])};"
            f"{matrix_attr_key}={_gff3_escape(hit['motif_id'])};"
            f"pvalue={hit['pvalue']};qvalue={hit['qvalue']};"
            f"matched_sequence={hit['matched_sequence']};"
            f"Description={_gff3_escape(desc)}"
        ),
    ]
    return "\t".join(fields)


def write_combined_gff3(genes_info, te_records, tfbs_hits, out_path):
    """Assemble every feature type (gene/mRNA/exon/intron from the
    annotation, mobile_genetic_element from the TE table, TF_binding_site
    from FIMO) into one GFF3 file, sorted by (seqid, start) as required
    for tabix indexing later. Also emits one ##sequence-region header
    line per seqid, a GFF3 convention that declares the min/max
    coordinates used on that sequence (not the full chromosome length -
    genome browsers don't require the true length here)."""
    lines = ["##gff-version 3"]

    # ##sequence-region header lines (one per seqid touched by any feature)
    seqid_extents = defaultdict(lambda: [None, None])
    body_lines = []

    for symbol, info in genes_info.items():
        gene_rec = info["gene"]
        if not gene_rec:
            continue
        gene_id = gene_rec["attrs"].get("ID", f"gene-{symbol}")
        body_lines.append((gene_rec["seqid"], gene_rec["start"], gene_record_to_gff3(gene_rec, symbol)))

        mrna_rec = info.get("mrna")
        if mrna_rec:
            mrna_line, parent_id = mrna_record_to_gff3(mrna_rec, symbol, gene_id)
            body_lines.append((mrna_rec["seqid"], mrna_rec["start"], mrna_line))
        else:
            parent_id = gene_id  # no transcript found - attach exons (if any) straight to the gene

        for i, exon in enumerate(sorted(info["exons"], key=lambda e: e["start"]), 1):
            body_lines.append((exon["seqid"], exon["start"], exon_record_to_gff3(exon, symbol, i, parent_id)))
        for i, intron in enumerate(sorted(info.get("introns", []), key=lambda x: x["start"]), 1):
            body_lines.append((intron["seqid"], intron["start"], intron_record_to_gff3(intron, symbol, i, parent_id)))

    for i, te in enumerate(sorted(te_records, key=lambda t: (t["seqid"], t["start"])), 1):
        body_lines.append((te["seqid"], te["start"], te_record_to_gff3(te, i)))

    for i, hit in enumerate(sorted(tfbs_hits, key=lambda h: (h["seqid"], h["start"])), 1):
        body_lines.append((hit["seqid"], hit["start"], tfbs_record_to_gff3(hit, i)))

    for seqid, start, line in body_lines:
        fields = line.split("\t")
        s, e = int(fields[3]), int(fields[4])
        lo, hi = seqid_extents[seqid]
        seqid_extents[seqid] = [s if lo is None else min(lo, s), e if hi is None else max(hi, e)]

    for seqid, (lo, hi) in seqid_extents.items():
        lines.append(f"##sequence-region {seqid} {lo} {hi}")

    body_lines.sort(key=lambda t: (t[0], t[1]))
    lines.extend(l for _, _, l in body_lines)

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[info] wrote {len(body_lines)} feature records to {out_path}")


def write_te_only_gff3(te_records, out_path):
    """Write just the mobile_genetic_element (TE) records as their own
    standalone GFF3. Meant to be loaded as a SECOND JBrowse 2 track
    alongside the combined one - genes/exons/TFBS in one file/track and
    TEs all blend into the same rendering/color by default, so splitting
    TEs out lets you give their track its own distinct color instead."""
    lines = ["##gff-version 3"]
    seqid_extents = defaultdict(lambda: [None, None])
    body_lines = []

    for i, te in enumerate(sorted(te_records, key=lambda t: (t["seqid"], t["start"])), 1):
        body_lines.append((te["seqid"], te["start"], te_record_to_gff3(te, i)))

    for seqid, start, line in body_lines:
        fields = line.split("\t")
        s, e = int(fields[3]), int(fields[4])
        lo, hi = seqid_extents[seqid]
        seqid_extents[seqid] = [s if lo is None else min(lo, s), e if hi is None else max(hi, e)]

    for seqid, (lo, hi) in seqid_extents.items():
        lines.append(f"##sequence-region {seqid} {lo} {hi}")

    body_lines.sort(key=lambda t: (t[0], t[1]))
    lines.extend(l for _, _, l in body_lines)

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[info] wrote {len(body_lines)} TE-only feature records to {out_path} (separate JBrowse 2 track)")


def derive_te_track_path(output_path):
    """<output>.gff3 -> <output>.transposable_elements.gff3 (same folder,
    same extension), for the standalone TE track written alongside the
    main combined output."""
    p = Path(output_path)
    suffix = p.suffix or ".gff3"
    return str(p.with_name(p.stem + ".transposable_elements" + suffix))


# --------------------------------------------------------------------------
# 6. Optional JBrowse 2 prep: bgzip + tabix index
#
# JBrowse 2 loads GFF3 either unindexed (fine for small files) or as a
# bgzip+tabix "GFF3Tabix" track (recommended - required for large files and
# for random-access performance). Our output is already written sorted by
# (seqid, start), which is exactly the order tabix requires, so indexing is
# just bgzip + tabix -p gff.
# --------------------------------------------------------------------------

def resolve_bgzip_engine(args):
    """Decide how (or whether) to bgzip/tabix-index the output, mirroring
    resolve_fimo_engine()'s auto-detect philosophy.

    Returns (engine, message) where engine is "local", "docker", or None
    (None means: skip indexing, message explains why).
      - args.bgzip_via_docker is True  -> Docker required, no fallback.
      - args.bgzip_via_docker is False -> local bgzip/tabix required, no
        Docker fallback (set via --no-bgzip-docker).
      - args.bgzip_via_docker is None  -> auto: prefer local binaries if
        already on PATH (no image pull needed), else try Docker, then give
        up gracefully.
    """
    def _local():
        missing = [t for t in ("bgzip", "tabix") if shutil.which(t) is None]
        if missing:
            return False, f"{', '.join(missing)} not found on PATH"
        return True, "bgzip/tabix found on PATH"

    if args.bgzip_via_docker is True:
        ok, msg = check_docker()
        return ("docker", msg) if ok else (None, msg)
    if args.bgzip_via_docker is False:
        ok, msg = _local()
        return ("local", msg) if ok else (None, msg)

    local_ok, local_msg = _local()
    if local_ok:
        return "local", local_msg
    docker_ok, docker_msg = check_docker()
    if docker_ok:
        return "docker", docker_msg
    return None, f"no bgzip/tabix engine available (local: {local_msg}; Docker: {docker_msg})"


def bgzip_and_index_gff3_local(gff3_path):
    """bgzip + tabix-index a sorted GFF3 using local htslib binaries."""
    gz_path = f"{gff3_path}.gz"
    try:
        subprocess.run(["bgzip", "-f", "-k", str(gff3_path)], check=True, capture_output=True, text=True)
        subprocess.run(["tabix", "-f", "-p", "gff", gz_path], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"local bgzip/tabix failed: {exc.stderr}") from exc
    return gz_path


def bgzip_and_index_gff3_docker(gff3_path, image):
    """bgzip + tabix-index a sorted GFF3 inside a Docker container (a small
    biocontainers htslib image), for machines with no local htslib install -
    same bind-mount-the-parent-directory pattern as run_fimo_docker."""
    gff3_path = Path(gff3_path).resolve()
    mount_dir = gff3_path.parent
    mount_point = "/data"
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{mount_dir}:{mount_point}",
        image,
        "bash", "-c",
        f"bgzip -f -k {mount_point}/{gff3_path.name} && tabix -f -p gff {mount_point}/{gff3_path.name}.gz",
    ]
    print(f"[info] running via Docker:\n       {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"docker bgzip/tabix failed (exit {result.returncode}):\n{result.stderr}")
    gz_path = f"{gff3_path}.gz"
    if not Path(gz_path).exists():
        raise RuntimeError(f"bgzip/tabix ran via Docker but no {gz_path} was produced")
    return str(gz_path)


def index_gff3_for_jbrowse(gff3_path, args):
    """Resolve local-vs-Docker bgzip/tabix and run it, warning (without
    aborting the whole pipeline - the un-indexed GFF3 is still valid and
    already written) if neither engine is available."""
    engine, msg = resolve_bgzip_engine(args)
    if engine is None:
        print(
            f"[warn] {msg} - skipping JBrowse 2 indexing.\n"
            "Install htslib (e.g. `conda install -c bioconda htslib` or `sudo apt-get install "
            "tabix`) or make sure Docker Desktop is running, then run manually:\n"
            f"  bgzip -k {gff3_path}\n"
            f"  tabix -p gff {gff3_path}.gz\n"
            "then add the resulting .gff3.gz as a JBrowse 2 GFF3Tabix track "
            "(JBrowse 2 needs the matching .gff3.gz.tbi alongside it).",
            file=sys.stderr,
        )
        return None
    print(f"[info] indexing for JBrowse 2 via {engine} ({msg})")
    try:
        if engine == "docker":
            gz_path = bgzip_and_index_gff3_docker(gff3_path, args.htslib_docker_image)
        else:
            gz_path = bgzip_and_index_gff3_local(gff3_path)
    except RuntimeError as exc:
        print(f"[error] {exc}", file=sys.stderr)
        return None
    print(f"[info] JBrowse 2-ready files written: {gz_path} and {gz_path}.tbi")
    return gz_path


# --------------------------------------------------------------------------
# 7. Per-species config files (optional convenience layer over the CLI)
#
# Typing out every --te-hits/--annotation-gff/--genome-fasta/... flag by
# hand each time gets old fast once you're running this against more than
# one species. --config points at a simple INI file with one [section] per
# species; --species picks which section to use (skip --species if the
# file only defines one). CLI flags always win over the config file, so you
# can still override any single value (e.g. --upstream 2000) without
# editing the file. See species_config.example.ini for a ready-to-copy
# template with a filled-in D. suzukii section.
# --------------------------------------------------------------------------

# Every key a config section is allowed to set, and how to interpret its
# string value. Keys not listed here are passed through as plain strings.
CONFIG_PATH_KEYS = {"te_hits", "annotation_gff", "genome_fasta", "jaspar_motif_file", "output", "work_dir"}
CONFIG_BOOL_KEYS = {"fimo_via_docker", "skip_tfbs", "bgzip_index", "bgzip_via_docker", "keep_work_dir", "include_cncc_maf_motif", "split_te_track"}
CONFIG_INT_KEYS = {"upstream", "downstream", "te_flank"}
CONFIG_FLOAT_KEYS = {"fimo_thresh"}


def load_species_config(config_path, species=None):
    """Read one [species] section out of an INI config file and return it
    as a dict of CLI-dest-name -> typed value. Relative paths in the file
    (te_hits, annotation_gff, genome_fasta, jaspar_motif_file, output,
    work_dir) are resolved relative to the config file's own directory, so
    a species' data files can simply sit alongside its config."""
    cp = configparser.ConfigParser()
    if not cp.read(config_path):
        raise FileNotFoundError(f"--config file not found or unreadable: {config_path}")

    sections = cp.sections()
    if not sections:
        raise ValueError(f"{config_path} has no [species] sections defined.")

    if species is None:
        if len(sections) == 1:
            species = sections[0]
        else:
            raise ValueError(
                f"{config_path} defines multiple species sections {sections} - "
                "pass --species <name> to pick one."
            )
    elif species not in cp:
        raise ValueError(f"No [{species}] section in {config_path}. Available: {sections}")

    section = cp[species]
    config_dir = Path(config_path).resolve().parent
    values = {}
    for raw_key in section:
        key = raw_key.replace("-", "_")
        if key in CONFIG_BOOL_KEYS:
            values[key] = section.getboolean(raw_key)
        elif key in CONFIG_INT_KEYS:
            values[key] = section.getint(raw_key)
        elif key in CONFIG_FLOAT_KEYS:
            values[key] = section.getfloat(raw_key)
        elif key in CONFIG_PATH_KEYS:
            p = Path(section[raw_key])
            values[key] = str(p if p.is_absolute() else (config_dir / p))
        else:
            values[key] = section[raw_key]
    print(f"[info] loaded species '{species}' from {config_path} ({len(values)} setting(s))")
    return values


# Built-in fallback for every optional setting - used when neither the CLI
# nor a --config file supplies a value. Kept as one dict (rather than
# scattered argparse `default=`s) so the three-way CLI > config > built-in
# precedence in main() is a single, auditable merge step.
ARG_DEFAULTS = {
    # None here means "auto-derived at runtime" (see resolve_args/main), not
    # "off" - the goal is that giving just --te-hits/--annotation-gff/
    # --genome-fasta produces the fullest output the local machine can
    # manage, with no other flags required.
    "output": None,                 # auto: "<te-hits filename stem>_combined.gff3"
    "genome_fasta": None,
    "sequence_source": None,
    "ncbi_email": None,
    "ncbi_api_key": None,
    "upstream": DEFAULT_UPSTREAM,
    "downstream": DEFAULT_DOWNSTREAM,
    "te_flank": DEFAULT_TE_FLANK,
    "jaspar_motif_file": None,
    "jaspar_release": JASPAR_RELEASE_DEFAULT,
    "jaspar_url": None,
    "jaspar_cache_dir": None,       # auto: ./jaspar_cache (reused across species/runs)
    "include_cncc_maf_motif": True,  # on by default: scan for the CncC:Maf-S ARE too (JASPAR has no such motif)
    "fimo_path": "fimo",
    "fimo_thresh": DEFAULT_FIMO_THRESH,
    "fimo_via_docker": None,        # auto: try Docker first, then a local fimo binary
    "fimo_docker_image": FIMO_DOCKER_IMAGE_DEFAULT,
    "skip_tfbs": False,
    "bgzip_index": True,            # auto: index if bgzip/tabix are available, skip (with a warning) if not
    "bgzip_via_docker": None,       # auto: prefer local bgzip/tabix if present, else try Docker
    "htslib_docker_image": DOCKER_HTSLIB_IMAGE_DEFAULT,
    "split_te_track": True,         # on by default: also write a standalone TE-only GFF3 for its own JBrowse track
    "work_dir": None,
    "keep_work_dir": False,
}


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def build_arg_parser():
    """Define every CLI flag. Most flags default to `None` here rather
    than their "real" default (see ARG_DEFAULTS below) - that's what lets
    resolve_args() tell "user didn't pass this flag" apart from "user
    explicitly passed the same value as the default", which matters for
    the CLI > config-file > built-in-default precedence chain."""
    p = argparse.ArgumentParser(
        description="Merge TE hits, gene/exon annotation, and JASPAR/FIMO-predicted "
                     "TFBS motifs into a single combined GFF3 for TE-affected Cyp genes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    cfg_group = p.add_argument_group(
        "species config file (optional shortcut for everything below - see species_config.example.ini)")
    cfg_group.add_argument("--config", default=None,
                            help="INI file with one [species] section per organism (te_hits/annotation_gff/"
                                 "genome_fasta/etc.). Any flag given explicitly on the command line still "
                                 "overrides the value from this file.")
    cfg_group.add_argument("--species", default=None,
                            help="Which [section] of --config to use (only needed if the file has more than one)")

    # NOTE: te-hits/annotation-gff/genome-fasta are not `required=True` here
    # because they may instead come from --config; presence of te-hits/
    # annotation-gff is validated in main() after the config file (if any)
    # has been merged in. These three are the ONLY inputs you need to give -
    # everything else below has an auto-detected/auto-derived default so a
    # bare `--te-hits X --annotation-gff Y --genome-fasta Z` run produces the
    # fullest output this machine is capable of (TFBS scan if fimo/Docker is
    # available, JBrowse-2 index if bgzip/tabix are available, output
    # filename derived automatically).
    p.add_argument("--te-hits", "--te-file", dest="te_hits", default=None,
                    help="Path to TE-hits table (format of D_suzukiiGenesAffectedByTE.txt). "
                         "Required, either here or via --config.")
    p.add_argument("--annotation-gff", "--gff", dest="annotation_gff", default=None,
                    help="Path to a real genome annotation GFF3 (genes/mRNA/exon features). "
                         "Required, either here or via --config.")
    p.add_argument("--genome-fasta", "--fasta", dest="genome_fasta", default=None,
                    help="Local genome FASTA (indexed on the fly, no extra deps needed). "
                         "Optional only if --sequence-source ncbi is used instead.")
    p.add_argument("--output", default=None,
                    help="Output combined GFF3 path (default: auto-derived from the --te-hits filename, "
                         "e.g. D_suzukiiGenesAffectedByTE.txt -> D_suzukiiGenesAffectedByTE_combined.gff3)")

    seq_group = p.add_argument_group("sequence source (needed for TFBS scanning)")
    seq_group.add_argument("--sequence-source", choices=["fasta", "ncbi"], default=None,
                            help="Force sequence source; default: fasta if --genome-fasta given, else ncbi")
    seq_group.add_argument("--ncbi-email", default=None, help="Email to include in NCBI E-utils requests (recommended)")
    seq_group.add_argument("--ncbi-api-key", default=None, help="NCBI API key (raises the E-utils rate limit)")

    win_group = p.add_argument_group("scan window sizes")
    win_group.add_argument("--upstream", type=int, default=None, help=f"bp upstream of TSS to scan (default: {DEFAULT_UPSTREAM})")
    win_group.add_argument("--downstream", type=int, default=None, help=f"bp downstream of TSS to scan (default: {DEFAULT_DOWNSTREAM})")
    win_group.add_argument("--te-flank", type=int, default=None, help=f"bp of flanking context added around each TE (default: {DEFAULT_TE_FLANK})")

    jaspar_group = p.add_argument_group("JASPAR / FIMO (all auto-detected by default - flags below are overrides)")
    jaspar_group.add_argument("--jaspar-motif-file", default=None, help="Local JASPAR MEME-format motif file (skips download)")
    jaspar_group.add_argument("--jaspar-release", default=None, help=f"JASPAR release year to download (default: {JASPAR_RELEASE_DEFAULT})")
    jaspar_group.add_argument("--jaspar-url", default=None, help="Override the full JASPAR motif download URL")
    jaspar_group.add_argument("--jaspar-cache-dir", default=None,
                               help="Folder to cache the downloaded JASPAR motif file in, reused across every "
                                    "future run/species so it's only downloaded once (default: ./jaspar_cache)")
    jaspar_group.add_argument("--cncc-maf-motif", dest="include_cncc_maf_motif", action="store_true", default=None,
                               help="Add the literature-derived Cap-n-Collar (CncC)/Maf-S antioxidant/xenobiotic "
                                    "response element (ARE) consensus motif to the scan (default: on already - "
                                    "this flag is mainly useful to re-enable it after a --config file turned it off)")
    jaspar_group.add_argument("--no-cncc-maf-motif", dest="include_cncc_maf_motif", action="store_false",
                               help="Do not add the CncC/Maf-S ARE consensus motif to the scan (default: "
                                    "included, since JASPAR's insect collection has no CncC/Maf-S entry at all - "
                                    "see the README for the source/caveats)")
    jaspar_group.add_argument("--fimo-path", default=None, help="Path to the fimo executable (default: fimo)")
    jaspar_group.add_argument("--fimo-thresh", type=float, default=None, help=f"FIMO p-value threshold (default: {DEFAULT_FIMO_THRESH})")
    jaspar_group.add_argument("--fimo-via-docker", dest="fimo_via_docker", action="store_true", default=None,
                               help="Force fimo to run via the official MEME Suite Docker container "
                                    "(default: auto - tried first, before a local fimo binary)")
    jaspar_group.add_argument("--no-fimo-docker", dest="fimo_via_docker", action="store_false",
                               help="Force a local fimo binary only, even if Docker is also available")
    jaspar_group.add_argument("--fimo-docker-image", default=None, help=f"Docker image to use with Docker-based fimo (default: {FIMO_DOCKER_IMAGE_DEFAULT})")

    p.add_argument("--skip-tfbs", action="store_true", default=None, help="Skip motif scanning; output genes+exons+TEs only")
    p.add_argument("--bgzip-index", dest="bgzip_index", action="store_true", default=None,
                    help="bgzip + tabix-index the output for direct loading into JBrowse 2 "
                         "(default: on automatically whenever bgzip/tabix or Docker are available)")
    p.add_argument("--no-bgzip-index", dest="bgzip_index", action="store_false",
                    help="Skip bgzip/tabix-indexing entirely")
    p.add_argument("--bgzip-via-docker", dest="bgzip_via_docker", action="store_true", default=None,
                    help="Force JBrowse-2 indexing via Docker (a small htslib container) instead of local "
                         "bgzip/tabix binaries (default: auto - local binaries preferred if present, "
                         "Docker used otherwise)")
    p.add_argument("--no-bgzip-docker", dest="bgzip_via_docker", action="store_false",
                    help="Only use local bgzip/tabix binaries for indexing; do not fall back to Docker")
    p.add_argument("--htslib-docker-image", default=None,
                    help=f"Docker image to use for bgzip/tabix indexing (default: {DOCKER_HTSLIB_IMAGE_DEFAULT})")
    p.add_argument("--split-te-track", dest="split_te_track", action="store_true", default=None,
                    help="Also write a standalone <output>.transposable_elements.gff3 with just the TE "
                         "(mobile_genetic_element) features, so they can be loaded as a second JBrowse 2 "
                         "track and given a distinct color instead of blending into the combined track "
                         "(default: on)")
    p.add_argument("--no-split-te-track", dest="split_te_track", action="store_false",
                    help="Do not write the separate TE-only GFF3 - TEs stay only in the combined output")
    p.add_argument("--work-dir", default=None, help="Working directory for intermediate files (default: a temp dir)")
    p.add_argument("--keep-work-dir", action="store_true", default=None, help="Do not delete the temp working directory on exit")
    return p


def resolve_args(argv=None):
    """Parse the CLI, merge in an optional --config file, and apply built-in
    defaults last. Precedence: explicit CLI flag > --config file > ARG_DEFAULTS.
    Returns a SimpleNamespace with the same attribute names main() has always
    used, so nothing downstream needs to know config files exist at all."""
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    merged = dict(ARG_DEFAULTS)
    merged["te_hits"] = None
    merged["annotation_gff"] = None

    if args.config:
        try:
            merged.update(load_species_config(args.config, args.species))
        except (FileNotFoundError, ValueError) as exc:
            parser.error(str(exc))
    elif args.species:
        parser.error("--species was given without --config.")

    cli_values = vars(args)
    for key in list(merged.keys()) + ["te_hits", "annotation_gff"]:
        cli_val = cli_values.get(key)
        if cli_val is not None:
            merged[key] = cli_val

    if not merged["te_hits"] or not merged["annotation_gff"]:
        parser.error(
            "--te-hits and --annotation-gff are required, either as CLI flags or via "
            "--config/--species (see species_config.example.ini)."
        )

    if not merged["output"]:
        merged["output"] = f"{Path(merged['te_hits']).stem}_combined.gff3"

    return SimpleNamespace(**merged)


def main(argv=None):
    """Entry point - runs the full pipeline described in the module
    docstring at the top of this file: parse inputs (1-2), build scan
    windows and run TFBS motif scanning (3-4, skippable), write the
    combined GFF3 (5), and optionally bgzip/tabix-index it plus a
    standalone TE-only track (6). The numbered comments below correspond
    to those same pipeline stages."""
    args = resolve_args(argv)

    # Intermediate files (extracted scan-window FASTA, downloaded/combined
    # motif file, raw fimo.tsv) live here. A temp dir is used unless the
    # caller wants to inspect them afterwards via --work-dir/--keep-work-dir.
    work_dir = Path(args.work_dir) if args.work_dir else Path(tempfile.mkdtemp(prefix="tfbs_te_gff_"))
    work_dir.mkdir(parents=True, exist_ok=True)
    print(f"[info] working directory: {work_dir}")

    # 1. TE hits
    te_records = parse_te_hits(args.te_hits)
    gene_symbols = sorted({t["gene"] for t in te_records})
    print(f"[info] {len(gene_symbols)} unique Cyp gene symbols referenced in TE-hits file")

    # 2. Annotation (genes/exons) - only for the genes named in the TE table
    genes_info = parse_annotation_gff3(args.annotation_gff, gene_symbols)
    unmatched = [g for g, v in genes_info.items() if not v["gene"]]
    if unmatched:
        print(f"[warn] {len(unmatched)} gene symbol(s) from the TE-hits file were not found in the annotation GFF3: "
              f"{unmatched[:10]}{' ...' if len(unmatched) > 10 else ''}", file=sys.stderr)

    # 3-4. Build scan windows (promoter + TE-flank) and run TFBS motif
    # scanning against them, unless the user opted out or no scan engine
    # (local fimo / Docker) is available on this machine.
    tfbs_hits = []
    if args.skip_tfbs:
        print("[info] --skip-tfbs given; producing gene+exon+TE GFF3 only")
    else:
        engine, engine_msg = resolve_fimo_engine(args)
        if engine is None:
            print(f"[warn] {engine_msg}\n{SETUP_HINT}", file=sys.stderr)
            print("[warn] continuing without TFBS scanning (genes/exons/TEs will still be written)", file=sys.stderr)
        else:
            print(f"[info] using {engine} for motif scanning: {engine_msg}")
            seq_source = args.sequence_source or ("fasta" if args.genome_fasta else "ncbi")
            if seq_source == "fasta" and not args.genome_fasta:
                print("[warn] --sequence-source fasta given but no --genome-fasta path provided.", file=sys.stderr)
                seq_source = None

            if seq_source:
                regions = build_scan_regions(genes_info, te_records, args.upstream, args.downstream, args.te_flank)
                if seq_source == "fasta":
                    seqs = extract_regions_local_fasta(args.genome_fasta, regions)
                else:
                    print("[info] fetching sequence windows from NCBI E-utils (this can be slow for many regions)...")
                    seqs = extract_regions_ncbi(regions, email=args.ncbi_email, api_key=args.ncbi_api_key)
                print(f"[info] extracted sequence for {len(seqs)}/{len(regions)} scan windows")

                if seqs:
                    region_lookup = {r["region_id"]: r for r in regions}
                    scan_fasta = work_dir / "scan_regions.fa"
                    write_regions_fasta(seqs, scan_fasta)

                    jaspar_cache_dir = Path(args.jaspar_cache_dir) if args.jaspar_cache_dir else (Path.cwd() / "jaspar_cache")
                    try:
                        motif_file = get_jaspar_motif_file(
                            args.jaspar_motif_file, jaspar_cache_dir,
                            release=args.jaspar_release, url_override=args.jaspar_url,
                        )
                    except Exception as exc:  # noqa: BLE001
                        print(f"[error] {exc}", file=sys.stderr)
                        motif_file = None

                    if motif_file and args.include_cncc_maf_motif:
                        motif_file = add_cncc_maf_are_motif(motif_file, work_dir)

                    if motif_file:
                        fimo_out_dir = work_dir / "fimo_out"
                        try:
                            if engine == "docker":
                                tsv_path = run_fimo_docker(
                                    args.fimo_docker_image, motif_file, str(scan_fasta),
                                    fimo_out_dir, args.fimo_thresh, work_dir,
                                )
                            else:
                                tsv_path = run_fimo(args.fimo_path, motif_file, str(scan_fasta), fimo_out_dir, args.fimo_thresh)
                            raw_hits = parse_fimo_tsv(tsv_path)
                            tfbs_hits = remap_fimo_hits(raw_hits, region_lookup)
                            print(f"[info] {len(tfbs_hits)} TFBS motif hits (p <= {args.fimo_thresh})")
                        except Exception as exc:  # noqa: BLE001
                            print(f"[error] FIMO run failed: {exc}", file=sys.stderr)

    # 5. Write combined GFF3
    write_combined_gff3(genes_info, te_records, tfbs_hits, args.output)

    if args.bgzip_index:
        index_gff3_for_jbrowse(args.output, args)

    # 6. Optional standalone TE-only track (own file+color in JBrowse 2,
    # instead of TEs blending into the combined gene/exon/TFBS track)
    if args.split_te_track:
        te_track_path = derive_te_track_path(args.output)
        write_te_only_gff3(te_records, te_track_path)
        if args.bgzip_index:
            index_gff3_for_jbrowse(te_track_path, args)

    if not args.keep_work_dir and not args.work_dir:
        shutil.rmtree(work_dir, ignore_errors=True)
    else:
        print(f"[info] intermediate files kept in {work_dir}")


if __name__ == "__main__":
    main()
