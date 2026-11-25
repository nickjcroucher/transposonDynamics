#!/bin/env bash
# author: ph-u
# script: sangerHPC_ncbiMeta.sh
# desc: NCBI genome assembly metadata only; fetch data reports (datasets) and reformat to TSV (dataformat) via apptainer
# in: bsub < sangerHPC_ncbiMeta.sh
# out: ../data/ncbi_metadata/{assembly_data_report.jsonl,assembly_metadata.tsv}
# arg: 0
# date: 20261124

#BSUB -G ph17g
#BSUB -o ../work/pfam-%J-%I.o
#BSUB -e ../work/pfam-%J-%I.e
#BSUB -q week
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]"
#BSUB -J "ncbiMeta"

PATH="/software/isg/private/wrappers/apptainer/1.4.0:$PATH"

set -euo pipefail

# ---- settings (edit here, or export before bsub) ---------------------------------
ACC_FILE=${ACC_FILE:-../pfam-data/acc.txt}           # one accession per line
OUT_DIR=${OUT_DIR:-../pfam-data/ncbi_metadata}
IMAGE=${IMAGE:-staphb/ncbi-datasets:18.39.0}
SIF=${SIF:-../work/ncbi-datasets_18.39.0.sif}   # pulled on first run if missing
CHUNK_SIZE=${CHUNK_SIZE:-1000}                  # accessions per request
RETRIES=${RETRIES:-4}                           # attempts per batch before it is split
EXTRA_OPTS=${EXTRA_OPTS:-}                      # extra datasets flags, e.g. "--assembly-version all"
FIELDS=${FIELDS:-}                              # empty = dataformat's default columns; compact example:
# FIELDS=accession,organism-name,organism-tax-id,assminfo-name,assminfo-level,assminfo-release-date,assminfo-submitter,assminfo-bioproject,assminfo-biosample-accession,assminfo-biosample-strain,assminfo-biosample-collection-date,assminfo-biosample-geo-loc-name,assminfo-biosample-host,assminfo-biosample-isolation-source,assmstats-total-sequence-len,assmstats-number-of-contigs,assmstats-contig-n50,assmstats-gc-percent,checkm-completeness,checkm-contamination

# compute nodes reach Docker Hub and NCBI through the Sanger web proxy (existing settings win)
export HTTP_PROXY=${HTTP_PROXY:-http://wwwcache.sanger.ac.uk:3128}
export HTTPS_PROXY=${HTTPS_PROXY:-$HTTP_PROXY}
export http_proxy=${http_proxy:-$HTTP_PROXY} https_proxy=${https_proxy:-$HTTPS_PROXY}

log() { echo "[$(date '+%F %T')] $*" >&2; }
die() { log "ERROR: $*"; exit 1; }

# ---- container: pull once, then run both tools from the local image -------------
command -v apptainer >/dev/null || die "apptainer not found (check the wrapper PATH above)"
[[ -r $ACC_FILE ]] || die "cannot read $ACC_FILE from $PWD"
mkdir -p "$OUT_DIR" "$(dirname "$SIF")"
OUT_DIR=$(cd "$OUT_DIR" && pwd -P)              # absolute path, bound into the container
if [[ ! -s $SIF ]]; then
  log "pulling docker://$IMAGE -> $SIF"
  apptainer pull "$SIF" "docker://$IMAGE" || die "apptainer pull failed (no internet/proxy on this node?)"
fi
datasets()   { apptainer exec -B "$OUT_DIR" "$SIF" datasets "$@"; }
dataformat() { apptainer exec -B "$OUT_DIR" "$SIF" dataformat "$@"; }
log "$(datasets --version 2>&1 | tail -n 1) | dataformat $(dataformat version 2>&1 | tail -n 1)"

# ---- clean the input list (on the host) ------------------------------------------
WORK=$OUT_DIR/.ncbi_meta_work; CHUNKS=$WORK/chunks; mkdir -p "$CHUNKS"
JSONL=$OUT_DIR/assembly_data_report.jsonl
TSV=$OUT_DIR/assembly_metadata.tsv

ACC_RE='^(GC[AF]_[0-9]{9}(\.[0-9]+)?|PRJ[A-Z]{2}[0-9]+)$'      # assembly or BioProject
tr -d '\r' < "$ACC_FILE" | sed 's/#.*//' | awk 'NF { print $1 }' | sort -u > "$WORK/all.txt"
grep -E  "$ACC_RE" "$WORK/all.txt" > "$WORK/valid.txt"               || true
grep -vE "$ACC_RE" "$WORK/all.txt" > "$OUT_DIR/invalid_accessions.txt" || true
n_valid=$(wc -l < "$WORK/valid.txt")
log "$(wc -l < "$WORK/all.txt") unique IDs: $n_valid usable, $(wc -l < "$OUT_DIR/invalid_accessions.txt") skipped"
(( n_valid > 0 )) || die "no GCA_/GCF_/PRJ accessions (virus nucleotide IDs need 'datasets summary virus genome accession')"

# ---- fixed batches; a resubmitted job skips batches that already finished --------
sig="$(md5sum < "$WORK/valid.txt" | cut -d' ' -f1) $CHUNK_SIZE"
if [[ $(cat "$CHUNKS/.sig" 2>/dev/null || true) != "$sig" ]]; then
  rm -rf "$CHUNKS"; mkdir -p "$CHUNKS"
  split -l "$CHUNK_SIZE" -d -a 5 "$WORK/valid.txt" "$CHUNKS/chunk_"
  echo "$sig" > "$CHUNKS/.sig"
fi

retry() {  # retry <attempts> <command...>
  local max=$1 i; shift
  for (( i = 1; i <= max; i++ )); do
    if "$@"; then return 0; fi
    if (( i < max )); then log "  attempt $i/$max failed; retrying in $(( 30 * i ))s"; sleep $(( 30 * i )); fi
  done
  return 1
}

summary_once() {  # $1 = accession list -> $1.jsonl (rewritten on every attempt)
  # shellcheck disable=SC2086  # EXTRA_OPTS is word-split on purpose
  datasets summary genome accession --inputfile "$1" --as-json-lines $EXTRA_OPTS > "$1.jsonl"
}

fetch() {  # $1 = list, $2 = attempts; a failing batch is halved until bad IDs are isolated
  local list=$1 n half
  if retry "$2" summary_once "$list"; then return 0; fi
  rm -f "$list.jsonl"
  n=$(wc -l < "$list")
  if (( n == 1 )); then cat "$list" >> "$CUR_CHUNK.failed"; return 0; fi
  half=$(( (n + 1) / 2 ))
  log "  batch of $n failed; splitting it to isolate problem accessions"
  head -n "$half" "$list" > "$list.a"
  tail -n +"$(( half + 1 ))" "$list" > "$list.b"
  fetch "$list.a" 1
  fetch "$list.b" 1
}

# ---- download metadata --------------------------------------------------------------
chunks=("$CHUNKS"/chunk_?????)
for (( k = 0; k < ${#chunks[@]}; k++ )); do
  CUR_CHUNK=${chunks[k]}
  if [[ -e $CUR_CHUNK.done ]]; then continue; fi
  rm -f "$CUR_CHUNK".?*                       # leftovers from an interrupted run
  log "batch $(( k + 1 ))/${#chunks[@]} ($(wc -l < "$CUR_CHUNK") accessions)"
  fetch "$CUR_CHUNK" "$RETRIES"
  touch "$CUR_CHUNK.done"
done

# ---- merge and reformat -------------------------------------------------------------
find "$CHUNKS" -name 'chunk_*.jsonl' -print0 | sort -z | xargs -0 -r cat | awk '!seen[$0]++' > "$JSONL"
find "$CHUNKS" -name 'chunk_*.failed' -exec cat {} + | sort -u > "$OUT_DIR/failed_requests.txt"
[[ -s $JSONL ]] || die "no metadata returned; see the .e log"

if [[ -n $FIELDS ]]; then
  dataformat tsv genome --inputfile "$JSONL" --fields "$FIELDS" > "$TSV"
else
  dataformat tsv genome --inputfile "$JSONL" > "$TSV"
fi

# ---- QC: requested accessions that came back without metadata ----------------------
dataformat tsv genome --inputfile "$JSONL" --fields accession --elide-header | sort -u > "$WORK/returned.txt"
awk 'FILENAME == ARGV[1] { got[$1] = 1; b = $1; sub(/\.[0-9]+$/, "", b); base[b] = 1; next }
     /^PRJ/      { next }                            # BioProjects expand to many assemblies
     /\.[0-9]+$/ { if (!($1 in got)) print; next }   # versioned: needs that exact version
     !($1 in base)                                   # unversioned: any version counts
    ' "$WORK/returned.txt" "$WORK/valid.txt" > "$OUT_DIR/missing_accessions.txt"

rows=$(( $(wc -l < "$TSV") - 1 ))
n_asm=$(wc -l < "$WORK/returned.txt")
n_missing=$(wc -l < "$OUT_DIR/missing_accessions.txt")
n_failed=$(wc -l < "$OUT_DIR/failed_requests.txt")
log "wrote $TSV: $rows rows for $n_asm assemblies"
if (( rows > n_asm )); then
  log "NOTE: some assemblies span several rows; set FIELDS to single-valued columns for one row each"
fi
if (( n_missing + n_failed > 0 )); then
  log "WARNING: $n_missing accession(s) without metadata, $n_failed still failing after retries"
fi

for f in invalid_accessions.txt missing_accessions.txt failed_requests.txt; do
  [[ -s $OUT_DIR/$f ]] || rm -f "$OUT_DIR/$f"
done
[[ ${KEEP_WORK:-0} == 1 ]] || rm -rf "$WORK"
log "done"
