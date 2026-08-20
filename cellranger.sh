#!/usr/bin/env bash
set -Eeuo pipefail

PROGRAM=${0##*/}

RNA_MM10_DEFAULT="/home/quanyiz/genome/refdata-gex-mm10-2020-A_tdTlacZ"
RNA_MM39_DEFAULT="/home/quanyiz/genome/refdata-gex-GRCm39-2024-A_tdTlacZ"
RNA_HG38_DEFAULT="/home/quanyiz/genome/refdata-gex-GRCh38-2024-A"
ARC_MM10_DEFAULT="/home/quanyiz/genome/refdata-cellranger-arc-mm10-2020-A-2.0.0"
ARC_MM39_DEFAULT="/home/quanyiz/genome/refdata-cellranger-arc-GRCm39-2024-A"
ARC_HG38_DEFAULT="/home/quanyiz/genome/refdata-cellranger-arc-GRCh38-2024-A"

mode=""
genome=""
reference=""
threads=20
memory_gb=200
rna_version=10
aggregate=0
aggregate_csv=""
normalize="none"
secondary_mode="default"
run_id=""
sample=""
libraries_csv=""
gex_path=""
atac_path=""
output_dir=""
create_bam="true"
dry_run=0
strict_scan=0
continue_on_error=0

usage() {
    cat <<'EOF'
Usage:
  cellranger.sh -m <rna|atac> -g <mm10|mm39|hg38> [--id ID] --sample PREFIX [options] FASTQ_DIRS
  cellranger.sh -m <rna|atac> -g <mm10|mm39|hg38> [options] FASTQ_ROOT
  cellranger.sh -m multiome -g <mm10|mm39|hg38> [--id ID] --libraries CSV [options]
  cellranger.sh -m multiome -g <mm10|mm39|hg38> [--id ID] --gex_path DIR --atac_path DIR [options]
  cellranger.sh -m <rna|atac|multiome> -a [--id ID] [--csv CSV | COUNT_ROOT] [options]

FASTQ_DIRS and PREFIX may be comma-separated. All lanes selected by PREFIX are
passed to one count invocation; FASTQ files are never concatenated.

For rna/atac count, omitting --id/--run while providing --sample/--samples
generates a deterministic run ID. Omitting both enables automatic mode.
FASTQ_ROOT is scanned recursively, exact filename prefixes are grouped across
lanes, and one sample is run at a time. Different prefixes are never inferred
to be the same biological sample.

For multiome count, --gex_path and --atac_path can replace --libraries. The
wrapper discovers and pairs samples, writes libraries CSV files in the current
directory, and runs pairs sequentially. For aggr, omitting --csv scans direct
children of COUNT_ROOT, or the current directory when COUNT_ROOT is omitted,
and writes the generated aggregation CSV in the current directory.

When --output-dir is omitted, the current working directory is the output root.
Automatic-mode status events are appended to
OUTPUT_ROOT/cellranger_auto_MODE/status.log.

Required options:
  -m MODE             rna, atac, or multiome
  -g GENOME           mm10, mm39, or hg38; not needed with -x
  -x PATH             Custom compatible reference; overrides -g
  --id, --run ID      Optional pipestance ID; generated from count/aggr inputs
  --sample, --samples PREFIX
                       Sample filename prefix for rna/atac count
  --libraries CSV     Libraries CSV for multiome count
  --gex_path DIR      GEX FASTQ root for generated multiome libraries CSV
  --atac_path DIR     ATAC FASTQ root for generated multiome libraries CSV
  -a, --aggr          Run aggr instead of count
  -c, --csv CSV       Aggregation CSV; optional when scanning COUNT_ROOT

Resource and behavior options:
  -t N                Local cores (default: 20)
  -r GB               Local memory in GB (default: 200)
  -u VERSION          RNA Cell Ranger version: 7, 8, 9, or 10 (default: 10)
  --normalize MODE    Aggr normalization (default: none)
  -n                  Legacy alias for depth/mapped normalization
  -s                  Enable secondary analysis
  -S                  Disable secondary analysis
  --create-bam BOOL   RNA/multiome BAM creation: true or false (default: true)
  --output-dir PATH   Explicit output directory; automatic-mode output root
  --strict-scan       Abort automatic mode if any discovered sample is blocked
  --continue-on-error Continue automatic mode after a failed sample
  --dry               Explicit mode: Cell Ranger dry run; automatic mode: print only
  -h, --help          Show this help

RNA aggr is supported by the wrapper but is not used by the mouse atlas workflow.
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 2
}

require_positive_integer() {
    local name=$1 value=$2
    [[ $value =~ ^[1-9][0-9]*$ ]] || die "$name must be a positive integer: $value"
}

generate_count_id() {
    local assay=$1 sample_name=$2
    local safe_name digest suffix maximum_base_length
    safe_name=$(printf '%s' "$sample_name" | LC_ALL=C tr -cs 'A-Za-z0-9_-' '_')
    while [[ $safe_name == _* ]]; do safe_name=${safe_name#_}; done
    while [[ $safe_name == *_ ]]; do safe_name=${safe_name%_}; done
    [[ -n $safe_name ]] || safe_name="sample"
    if [[ $safe_name == "$sample_name" && ${#safe_name} -le 64 ]]; then
        printf '%s\n' "$safe_name"
        return 0
    fi
    digest=$(printf '%s\0%s' "$assay" "$sample_name" |
        sha256sum | awk '{print substr($1, 1, 8)}')
    suffix="__${digest}"
    maximum_base_length=$((64 - ${#suffix}))
    printf '%s%s\n' "${safe_name:0:maximum_base_length}" "$suffix"
}

generate_csv_id() {
    local assay=$1 operation=$2 csv_file=$3 normalization=$4 reference_key=$5
    local filename stem safe_name digest suffix base maximum_base_length
    filename=${csv_file##*/}
    stem=${filename%.*}
    safe_name=$(printf '%s' "$stem" | LC_ALL=C tr -cs 'A-Za-z0-9_-' '_')
    while [[ $safe_name == _* ]]; do safe_name=${safe_name#_}; done
    while [[ $safe_name == *_ ]]; do safe_name=${safe_name%_}; done
    [[ -n $safe_name ]] || safe_name="inputs"
    digest=$(
        {
            printf '%s\0%s\0%s\0%s\0' \
                "$assay" "$operation" "$normalization" "$reference_key"
            sha256sum "$csv_file"
        } | sha256sum | awk '{print substr($1, 1, 8)}'
    )
    suffix="__${digest}"
    base="${assay}_${operation}_${safe_name}"
    maximum_base_length=$((64 - ${#suffix}))
    printf '%s%s\n' "${base:0:maximum_base_length}" "$suffix"
}


generate_aggr_csv_from_counts() {
    local assay=$1 count_root=$2 destination=$3
    python3 - "$assay" "$count_root" "$destination" <<'PY'
import csv
import os
import sys
from pathlib import Path

assay, root_arg, destination_arg = sys.argv[1:]
root = Path(root_arg).resolve()
destination = Path(destination_arg)

layouts = {
    "rna": (
        ["sample_id", "molecule_h5"],
        ["molecule_info.h5"],
    ),
    "atac": (
        ["library_id", "fragments", "cells"],
        ["fragments.tsv.gz", "singlecell.csv"],
    ),
    "multiome": (
        [
            "library_id",
            "atac_fragments",
            "per_barcode_metrics",
            "gex_molecule_info",
        ],
        [
            "atac_fragments.tsv.gz",
            "per_barcode_metrics.csv",
            "gex_molecule_info.h5",
        ],
    ),
}
header, required_names = layouts[assay]
rows = []

for child in sorted(root.iterdir(), key=lambda path: (path.name.lower(), path.name)):
    generated_aggr_prefixes = ("aggr_", f"{assay}_aggr_")
    if (
        child.name.startswith(".")
        or child.name.startswith(generated_aggr_prefixes)
        or not child.is_dir()
    ):
        continue
    resolved_child = child.resolve()
    required_paths = [resolved_child / "outs" / name for name in required_names]
    if not all(path.is_file() and path.stat().st_size > 0 for path in required_paths):
        continue
    rows.append([child.name, *[str(path.resolve()) for path in required_paths]])

if not rows:
    required_text = ", ".join(f"outs/{name}" for name in required_names)
    raise SystemExit(
        f"No complete {assay} count directories were found directly under {root}; "
        f"required files: {required_text}"
    )

destination.parent.mkdir(parents=True, exist_ok=True)
temporary = destination.with_name(destination.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
os.replace(temporary, destination)
print(f"Generated aggregation CSV: {destination.resolve()} rows={len(rows)}")
PY
}


auto_discover_fastqs() {
    local scan_mode=$1
    local scan_root=$2
    local scan_output_root=$3
    local scan_manifest=$4

    python3 - "$scan_mode" "$scan_root" "$scan_output_root" "$scan_manifest" <<'PY'
import csv
import hashlib
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

mode, root_arg, output_arg, manifest_arg = sys.argv[1:]
root = Path(root_arg).resolve()
output_root = Path(output_arg).resolve()
manifest = Path(manifest_arg)
output_is_nested = output_root != root and root in output_root.parents

patterns = [
    re.compile(
        r"^(?P<prefix>.+)_S\d+_(?P<lane>L\d{3})_"
        r"(?P<read>I1|I2|R1|R2|R3)_(?P<chunk>\d{3})\.(?:fastq|fq)\.gz$",
        re.IGNORECASE,
    ),
    re.compile(
        r"^(?P<prefix>.+)_(?P<lane>L\d{3})_"
        r"(?P<read>I1|I2|R1|R2|R3)_(?P<chunk>\d{3})\.(?:fastq|fq)\.gz$",
        re.IGNORECASE,
    ),
]
required_reads = {"rna": {"R1", "R2"}, "atac": {"R1", "R2", "R3"}}[mode]
records = defaultdict(list)
seen_paths = set()

for current, directories, filenames in os.walk(root, topdown=True, followlinks=False):
    current_path = Path(current)
    kept_directories = []
    for name in directories:
        if name.startswith("."):
            continue
        candidate = (current_path / name).resolve()
        if output_is_nested and (
            candidate == output_root or output_root in candidate.parents
        ):
            continue
        kept_directories.append(name)
    directories[:] = kept_directories

    for filename in sorted(filenames):
        match = None
        for pattern in patterns:
            match = pattern.match(filename)
            if match:
                break
        if match is None:
            continue

        file_path = (current_path / filename).resolve()
        path_key = str(file_path)
        if path_key in seen_paths:
            continue
        seen_paths.add(path_key)
        try:
            size = file_path.stat().st_size
        except OSError:
            size = 0
        records[match.group("prefix")].append(
            {
                "path": path_key,
                "parent": str(file_path.parent),
                "lane": match.group("lane").upper(),
                "read": match.group("read").upper(),
                "chunk": match.group("chunk"),
                "size": size,
            }
        )

used_ids = {}
rows = []

for prefix in sorted(records, key=lambda value: (value.lower(), value)):
    sample_records = records[prefix]
    errors = []
    selected = []

    if any(character in prefix for character in (",", "\t", "\r", "\n")):
        errors.append("sample prefix contains a comma, tab, or newline")

    logical_groups = defaultdict(lambda: defaultdict(list))
    for record in sample_records:
        logical_key = (record["lane"], record["chunk"])
        logical_groups[logical_key][record["parent"]].append(record)

    for logical_key in sorted(logical_groups):
        lane, chunk = logical_key
        parent_groups = logical_groups[logical_key]
        complete_parents = []
        parent_summaries = []

        for parent in sorted(parent_groups):
            parent_records = parent_groups[parent]
            by_read = defaultdict(list)
            for record in parent_records:
                by_read[record["read"]].append(record)
            duplicate_reads = sorted(
                read_name for read_name, values in by_read.items() if len(values) > 1
            )
            available_reads = set(by_read)
            parent_summaries.append(
                f"{parent}={','.join(sorted(available_reads)) or 'none'}"
            )
            if duplicate_reads:
                errors.append(
                    f"duplicate read files for {lane}/{chunk} in {parent}: "
                    f"{','.join(duplicate_reads)}"
                )
                continue
            if required_reads.issubset(available_reads):
                complete_parents.append(parent)

        if len(complete_parents) == 1:
            selected.extend(parent_groups[complete_parents[0]])
        elif len(complete_parents) == 0:
            errors.append(
                f"no complete {lane}/{chunk} read set; "
                + "; ".join(parent_summaries)
            )
        else:
            errors.append(
                f"ambiguous complete copies for {lane}/{chunk}: "
                + ",".join(complete_parents)
            )

    if mode == "rna" and any(record["read"] == "R3" for record in selected):
        errors.append("R3 was found in an RNA read set; verify the assay type")

    selected_paths = sorted({record["path"] for record in selected})
    selected_dirs = sorted({record["parent"] for record in selected})
    if any("," in directory for directory in selected_dirs):
        errors.append("a FASTQ directory contains a comma")

    safe_id = re.sub(r"[^A-Za-z0-9_-]+", "_", prefix).strip("_-") or "sample"
    digest_source = prefix + "|" + "|".join(selected_paths)
    digest = hashlib.sha1(digest_source.encode("utf-8")).hexdigest()[:8]
    if len(safe_id) > 64:
        safe_id = f"{safe_id[:54]}__{digest}"
    if safe_id in used_ids and used_ids[safe_id] != prefix:
        safe_id = f"{safe_id[:54]}__{digest}"
    collision_index = 1
    base_id = safe_id
    while safe_id in used_ids and used_ids[safe_id] != prefix:
        suffix = f"__{collision_index:02d}"
        safe_id = base_id[: 64 - len(suffix)] + suffix
        collision_index += 1
    used_ids[safe_id] = prefix

    lane_count = len({(record["parent"], record["lane"]) for record in selected})
    selected_size = sum(record["size"] for record in selected)
    decision = "blocked" if errors else "ready"
    reason = " | ".join(dict.fromkeys(errors)) if errors else "complete read sets"
    display_prefix = (
        prefix.replace("\t", "\\t").replace("\r", "\\r").replace("\n", "\\n")
    )
    rows.append(
        {
            "mode": mode,
            "sample_prefix": display_prefix,
            "run_id": safe_id,
            "fastq_dirs": ",".join(selected_dirs),
            "lane_count": lane_count,
            "fastq_file_count": len(selected_paths),
            "fastq_bytes": selected_size,
            "decision": decision,
            "reason": reason,
        }
    )

manifest.parent.mkdir(parents=True, exist_ok=True)
temporary = manifest.with_name(manifest.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "mode",
            "sample_prefix",
            "run_id",
            "fastq_dirs",
            "lane_count",
            "fastq_file_count",
            "fastq_bytes",
            "decision",
            "reason",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
temporary.replace(manifest)
PY
}

auto_output_is_complete() {
    local assay=$1
    local run_directory=$2
    case "$assay" in
        rna)
            [[ -s "$run_directory/outs/filtered_feature_bc_matrix.h5" ]]
            ;;
        atac)
            [[ -s "$run_directory/outs/fragments.tsv.gz" &&
               -s "$run_directory/outs/singlecell.csv" ]]
            ;;
        multiome)
            [[ -s "$run_directory/outs/atac_fragments.tsv.gz" &&
               -s "$run_directory/outs/per_barcode_metrics.csv" &&
               -s "$run_directory/outs/gex_molecule_info.h5" ]]
            ;;
        *)
            return 1
            ;;
    esac
}

append_auto_status() {
    local status_file=$1
    local assay=$2
    local prefix=$3
    local child_id=$4
    local state=$5
    local exit_code=$6
    local output_path=$7
    local log_path=$8
    local note=$9
    local timestamp
    timestamp=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$timestamp" "$assay" "$prefix" "$child_id" "$state" "$exit_code" \
        "$output_path" "$log_path" "$note" >> "$status_file"
}

run_automatic_count() {
    [[ $mode == rna || $mode == atac ]] ||
        die "Automatic FASTQ discovery is supported only for rna and atac count"
    [[ ${#positionals[@]} -eq 1 ]] ||
        die "Automatic mode requires exactly one FASTQ_ROOT argument"

    local scan_root
    scan_root=${positionals[0]}
    [[ -d $scan_root ]] || die "FASTQ root does not exist: $scan_root"
    scan_root=$(realpath "$scan_root")

    local auto_output_root
    if [[ -n $output_dir ]]; then
        mkdir -p "$output_dir"
        auto_output_root=$(realpath "$output_dir")
    else
        auto_output_root=$(pwd -P)
    fi

    local state_dir="$auto_output_root/cellranger_auto_$mode"
    local discovery_file="$state_dir/discovery.tsv"
    local status_file="$state_dir/status.log"
    local logs_dir="$state_dir/logs"
    mkdir -p "$logs_dir"

    command -v flock >/dev/null 2>&1 ||
        die "flock is required for automatic mode"
    local legacy_lock="$auto_output_root/.cellranger_auto_$mode/launcher.lock"
    if [[ -f $legacy_lock ]] &&
        ! flock -n "$legacy_lock" -c true 2>/dev/null; then
        die "An older automatic $mode launcher holds $legacy_lock"
    fi
    local lock_file="$state_dir/launcher.lock"
    local lock_fd
    exec {lock_fd}>"$lock_file"
    flock -n "$lock_fd" ||
        die "Another automatic $mode launcher holds $lock_file"

    auto_discover_fastqs "$mode" "$scan_root" "$auto_output_root" "$discovery_file"

    local discovered ready blocked
    discovered=$(awk -F '\t' 'NR > 1 {count++} END {print count + 0}' "$discovery_file")
    ready=$(awk -F '\t' 'NR > 1 && $8 == "ready" {count++} END {print count + 0}' "$discovery_file")
    blocked=$(awk -F '\t' 'NR > 1 && $8 == "blocked" {count++} END {print count + 0}' "$discovery_file")

    echo "Automatic FASTQ discovery: mode=$mode discovered=$discovered ready=$ready blocked=$blocked"
    echo "Discovery manifest: $discovery_file"

    if [[ $blocked -gt 0 ]]; then
        echo "Blocked samples are recorded in the discovery manifest." >&2
        awk -F '\t' 'BEGIN {OFS="\t"} NR == 1 || ($8 == "blocked" && shown++ < 20) {
            print $2, $3, $8, $9
        }' "$discovery_file" >&2
        [[ $strict_scan -eq 0 ]] ||
            die "Automatic discovery found $blocked blocked sample(s)"
    fi
    [[ $ready -gt 0 ]] || die "No runnable samples were discovered"

    if [[ ! -s $status_file ]]; then
        printf 'timestamp\tmode\tsample_prefix\trun_id\tstatus\texit_code\toutput_path\tlog_path\tnote\n' \
            > "$status_file"
    fi

    local scan_timestamp
    scan_timestamp=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    awk -F '\t' -v OFS='\t' -v stamp="$scan_timestamp" \
        -v output_root="$auto_output_root" -v logs="$logs_dir" \
        'NR > 1 && $8 == "blocked" {
            print stamp, $1, $2, $3, "blocked", "", output_root "/" $3,
                  logs "/" $3 ".log", $9
        }' "$discovery_file" >> "$status_file"

    local script_path
    script_path=$(realpath "$0")
    local failures=0
    local row_mode prefix child_id fastq_dirs lane_count fastq_count fastq_bytes
    local decision reason run_directory log_path child_status
    local -a child

    while IFS=$'\t' read -r row_mode prefix child_id fastq_dirs lane_count \
        fastq_count fastq_bytes decision reason; do
        [[ $decision == ready ]] || continue

        if [[ -e "$state_dir/STOP" ]]; then
            echo "STOP detected; no additional samples will be started."
            append_auto_status "$status_file" "$mode" "$prefix" "$child_id" \
                "stopped" "" "" "" "STOP control file detected"
            break
        fi
        while [[ -e "$state_dir/PAUSE" ]]; do
            echo "PAUSE detected; waiting before $child_id"
            sleep 30
            if [[ -e "$state_dir/STOP" ]]; then
                break
            fi
        done
        if [[ -e "$state_dir/STOP" ]]; then
            echo "STOP detected; no additional samples will be started."
            append_auto_status "$status_file" "$mode" "$prefix" "$child_id" \
                "stopped" "" "" "" "STOP control file detected"
            break
        fi

        run_directory="$auto_output_root/$child_id"
        log_path="$logs_dir/$child_id.log"

        if auto_output_is_complete "$mode" "$run_directory"; then
            echo "Skipping complete sample: $child_id"
            append_auto_status "$status_file" "$mode" "$prefix" "$child_id" \
                "complete" "0" "$run_directory" "$log_path" \
                "validated output already exists"
            continue
        fi

        child=(
            "$script_path"
            -m "$mode"
            -x "$reference"
            -t "$threads"
            -r "$memory_gb"
        )
        if [[ $mode == rna ]]; then
            child+=(-u "$rna_version" --create-bam "$create_bam")
        fi
        case "$secondary_mode" in
            on) child+=(-s) ;;
            off) child+=(-S) ;;
        esac
        child+=(
            --id "$child_id"
            --sample "$prefix"
            --output-dir "$run_directory"
            "$fastq_dirs"
        )

        if [[ $dry_run -eq 1 ]]; then
            printf 'Would run:'
            printf ' %q' "${child[@]}"
            printf '\n'
            append_auto_status "$status_file" "$mode" "$prefix" "$child_id" \
                "dry_run" "0" "$run_directory" "$log_path" \
                "command printed without execution"
            continue
        fi

        {
            printf '%s Starting %s sample=%q lanes=%s files=%s bytes=%s\n' \
                "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$child_id" "$prefix" \
                "$lane_count" "$fastq_count" "$fastq_bytes"
            printf 'Command:'
            printf ' %q' "${child[@]}"
            printf '\n'
        } | tee -a "$log_path"
        append_auto_status "$status_file" "$mode" "$prefix" "$child_id" \
            "running" "" "$run_directory" "$log_path" "count started"

        if (
            cd "$auto_output_root"
            "${child[@]}"
        ) 2>&1 | tee -a "$log_path"; then
            child_status=0
        else
            child_status=${PIPESTATUS[0]}
        fi

        if [[ $child_status -eq 0 ]] &&
            auto_output_is_complete "$mode" "$run_directory"; then
            append_auto_status "$status_file" "$mode" "$prefix" "$child_id" \
                "complete" "0" "$run_directory" "$log_path" \
                "validated required output files"
            echo "Completed sample: $child_id"
        else
            if [[ $child_status -eq 0 ]]; then
                child_status=3
                echo "ERROR: Required output files are missing for $child_id" |
                    tee -a "$log_path" >&2
            fi
            append_auto_status "$status_file" "$mode" "$prefix" "$child_id" \
                "failed" "$child_status" "$run_directory" "$log_path" \
                "count failed or required outputs are missing"
            failures=$((failures + 1))
            echo "ERROR: Sample $child_id failed with status $child_status" >&2
            if [[ $continue_on_error -eq 0 ]]; then
                echo "Stopping after the first failure. Re-run the same command to resume." >&2
                echo "Status log: $status_file" >&2
                return "$child_status"
            fi
        fi
    done < "$discovery_file"

    echo "Automatic count pass finished. Status log: $status_file"
    if [[ $failures -gt 0 ]]; then
        return 1
    fi
    return 0
}


run_multiome_paths() {
    local resolved_gex resolved_atac
    resolved_gex=$(realpath "$gex_path")
    resolved_atac=$(realpath "$atac_path")

    local auto_output_root
    if [[ -n $output_dir ]]; then
        mkdir -p "$output_dir"
        auto_output_root=$(realpath "$output_dir")
    else
        auto_output_root=$(pwd -P)
    fi

    local state_dir="$auto_output_root/cellranger_auto_multiome"
    local gex_manifest="$state_dir/gex_discovery.tsv"
    local atac_manifest="$state_dir/atac_discovery.tsv"
    local pairs_file="$state_dir/pairs.tsv"
    local libraries_dir
    libraries_dir=$(pwd -P)
    local status_file="$state_dir/status.log"
    local logs_dir="$state_dir/logs"
    mkdir -p "$logs_dir"

    command -v flock >/dev/null 2>&1 ||
        die "flock is required for automatic multiome path mode"
    local legacy_lock="$auto_output_root/.cellranger_auto_multiome/launcher.lock"
    if [[ -f $legacy_lock ]] &&
        ! flock -n "$legacy_lock" -c true 2>/dev/null; then
        die "An older automatic multiome launcher holds $legacy_lock"
    fi
    local lock_file="$state_dir/launcher.lock"
    local lock_fd
    exec {lock_fd}>"$lock_file"
    flock -n "$lock_fd" ||
        die "Another automatic multiome launcher holds $lock_file"

    auto_discover_fastqs "rna" "$resolved_gex" "$auto_output_root" "$gex_manifest"
    auto_discover_fastqs "atac" "$resolved_atac" "$auto_output_root" "$atac_manifest"

    python3 - "$gex_manifest" "$atac_manifest" "$pairs_file" \
        "$libraries_dir" "$reference" "$run_id" <<'PY'
import csv
import hashlib
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

gex_arg, atac_arg, pairs_arg, libraries_arg, reference, explicit_id = sys.argv[1:]
gex_path = Path(gex_arg)
atac_path = Path(atac_arg)
pairs_path = Path(pairs_arg)
libraries_dir = Path(libraries_arg)


def read_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pairing_key(prefix):
    value = prefix.lower()
    match = re.match(r"^(.+?[-_](?:s)?\d+)(?:[-_].*)?$", value)
    if match:
        value = match.group(1)
    parts = [part for part in re.split(r"[-_]+", value) if part]
    assay_tokens = {
        "gex", "rna", "geneexpression", "atac", "dna", "chromatin",
        "accessibility", "arc", "multiome",
    }
    parts = [part for part in parts if part not in assay_tokens]
    return "_".join(parts) or re.sub(
        r"[^a-z0-9]+", "_", prefix.lower()
    ).strip("_")


gex_groups = defaultdict(list)
atac_groups = defaultdict(list)
for row in read_rows(gex_path):
    gex_groups[pairing_key(row["sample_prefix"])].append(row)
for row in read_rows(atac_path):
    atac_groups[pairing_key(row["sample_prefix"])].append(row)

keys = sorted(set(gex_groups) | set(atac_groups))
candidate_keys = []
evaluated = {}

for key in keys:
    gex_rows = gex_groups.get(key, [])
    atac_rows = atac_groups.get(key, [])
    errors = []
    blocked_gex = [row for row in gex_rows if row["decision"] != "ready"]
    blocked_atac = [row for row in atac_rows if row["decision"] != "ready"]
    ready_gex = [row for row in gex_rows if row["decision"] == "ready"]
    ready_atac = [row for row in atac_rows if row["decision"] == "ready"]

    for row in blocked_gex:
        errors.append(f"GEX {row['sample_prefix']}: {row['reason']}")
    for row in blocked_atac:
        errors.append(f"ATAC {row['sample_prefix']}: {row['reason']}")
    if len(ready_gex) == 0:
        errors.append("no ready GEX prefix")
    elif len(ready_gex) > 1:
        errors.append(
            "ambiguous GEX prefixes: "
            + ",".join(row["sample_prefix"] for row in ready_gex)
        )
    if len(ready_atac) == 0:
        errors.append("no ready ATAC prefix")
    elif len(ready_atac) > 1:
        errors.append(
            "ambiguous ATAC prefixes: "
            + ",".join(row["sample_prefix"] for row in ready_atac)
        )

    if not errors:
        candidate_keys.append(key)
    evaluated[key] = (ready_gex, ready_atac, errors)

if explicit_id and len(candidate_keys) != 1:
    raise SystemExit(
        "--id/--run in multiome path mode requires exactly one unambiguous pair; "
        f"found {len(candidate_keys)}"
    )

used_ids = {}
output_rows = []
libraries_dir.mkdir(parents=True, exist_ok=True)

for key in keys:
    ready_gex, ready_atac, errors = evaluated[key]
    gex_row = ready_gex[0] if len(ready_gex) == 1 else {}
    atac_row = ready_atac[0] if len(ready_atac) == 1 else {}
    identity = "|".join(
        [
            "multiome", "count", key,
            gex_row.get("sample_prefix", ""),
            gex_row.get("fastq_dirs", ""),
            atac_row.get("sample_prefix", ""),
            atac_row.get("fastq_dirs", ""),
            reference,
        ]
    )
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:8]
    safe_key = re.sub(r"[^A-Za-z0-9_-]+", "_", key).strip("_-") or "pair"
    suffix = f"__{digest}"
    base = f"multiome_count_{safe_key}"
    generated_id = base[: 64 - len(suffix)] + suffix
    child_id = explicit_id if explicit_id and not errors else generated_id
    if child_id in used_ids and used_ids[child_id] != key:
        errors.append(f"generated run ID collision with pairing key {used_ids[child_id]}")
    used_ids[child_id] = key

    libraries_csv = ""
    if not errors:
        libraries_path = libraries_dir / f"multiome_libraries_{child_id}.csv"
        temporary = libraries_path.with_name(libraries_path.name + ".tmp")
        with temporary.open("w", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["fastqs", "sample", "library_type"])
            writer.writerow(
                [gex_row["fastq_dirs"], gex_row["sample_prefix"], "Gene Expression"]
            )
            writer.writerow(
                [
                    atac_row["fastq_dirs"],
                    atac_row["sample_prefix"],
                    "Chromatin Accessibility",
                ]
            )
        os.replace(temporary, libraries_path)
        libraries_csv = str(libraries_path.resolve())

    output_rows.append(
        {
            "pair_key": key,
            "gex_prefix": gex_row.get("sample_prefix", ""),
            "atac_prefix": atac_row.get("sample_prefix", ""),
            "run_id": child_id,
            "gex_fastq_dirs": gex_row.get("fastq_dirs", ""),
            "atac_fastq_dirs": atac_row.get("fastq_dirs", ""),
            "libraries_csv": libraries_csv,
            "decision": "blocked" if errors else "ready",
            "reason": " | ".join(dict.fromkeys(errors))
            if errors
            else "unambiguous pair",
        }
    )

pairs_path.parent.mkdir(parents=True, exist_ok=True)
temporary = pairs_path.with_name(pairs_path.name + ".tmp")
with temporary.open("w", newline="") as handle:
    fields = [
        "pair_key", "gex_prefix", "atac_prefix", "run_id",
        "gex_fastq_dirs", "atac_fastq_dirs", "libraries_csv",
        "decision", "reason",
    ]
    writer = csv.DictWriter(
        handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(output_rows)
os.replace(temporary, pairs_path)
PY

    local discovered ready blocked
    discovered=$(awk -F '\t' 'NR > 1 {count++} END {print count + 0}' "$pairs_file")
    ready=$(awk -F '\t' 'NR > 1 && $8 == "ready" {count++} END {print count + 0}' "$pairs_file")
    blocked=$(awk -F '\t' 'NR > 1 && $8 == "blocked" {count++} END {print count + 0}' "$pairs_file")

    echo "Multiome pair discovery: discovered=$discovered ready=$ready blocked=$blocked"
    echo "GEX discovery: $gex_manifest"
    echo "ATAC discovery: $atac_manifest"
    echo "Pair manifest: $pairs_file"

    if [[ $blocked -gt 0 ]]; then
        echo "Blocked multiome pairs are recorded in the pair manifest." >&2
        awk -F '\t' 'BEGIN {OFS="\t"} NR == 1 || ($8 == "blocked" && shown++ < 20) {
            print $1, $2, $3, $8, $9
        }' "$pairs_file" >&2
        [[ $strict_scan -eq 0 ]] ||
            die "Multiome discovery found $blocked blocked pair(s)"
    fi
    [[ $ready -gt 0 ]] || die "No runnable multiome pairs were discovered"

    if [[ ! -s $status_file ]]; then
        printf 'timestamp\tmode\tsample_prefix\trun_id\tstatus\texit_code\toutput_path\tlog_path\tnote\n' \
            > "$status_file"
    fi

    local scan_timestamp
    scan_timestamp=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    awk -F '\t' -v OFS='\t' -v stamp="$scan_timestamp" \
        -v output_root="$auto_output_root" -v logs="$logs_dir" \
        'NR > 1 && $8 == "blocked" {
            print stamp, "multiome", $1, $4, "blocked", "", output_root "/" $4,
                  logs "/" $4 ".log", $9
        }' "$pairs_file" >> "$status_file"

    local script_path
    script_path=$(realpath "$0")
    local failures=0
    local pair_key gex_prefix atac_prefix child_id gex_dirs atac_dirs libraries
    local decision reason pair_label run_directory log_path child_status
    local -a child

    while IFS=$'\t' read -r pair_key gex_prefix atac_prefix child_id gex_dirs \
        atac_dirs libraries decision reason; do
        [[ $decision == ready ]] || continue
        pair_label="$pair_key:$gex_prefix+$atac_prefix"

        if [[ -e "$state_dir/STOP" ]]; then
            echo "STOP detected; no additional multiome pairs will be started."
            append_auto_status "$status_file" "multiome" "$pair_label" "$child_id" \
                "stopped" "" "" "" "STOP control file detected"
            break
        fi
        while [[ -e "$state_dir/PAUSE" ]]; do
            echo "PAUSE detected; waiting before $child_id"
            sleep 30
            [[ -e "$state_dir/STOP" ]] && break
        done
        if [[ -e "$state_dir/STOP" ]]; then
            echo "STOP detected; no additional multiome pairs will be started."
            append_auto_status "$status_file" "multiome" "$pair_label" "$child_id" \
                "stopped" "" "" "" "STOP control file detected"
            break
        fi

        run_directory="$auto_output_root/$child_id"
        log_path="$logs_dir/$child_id.log"
        if auto_output_is_complete "multiome" "$run_directory"; then
            echo "Skipping complete multiome pair: $child_id"
            append_auto_status "$status_file" "multiome" "$pair_label" "$child_id" \
                "complete" "0" "$run_directory" "$log_path" \
                "validated output already exists"
            continue
        fi

        child=(
            "$script_path" -m multiome -x "$reference"
            -t "$threads" -r "$memory_gb"
            --create-bam "$create_bam"
        )
        case "$secondary_mode" in
            on) child+=(-s) ;;
            off) child+=(-S) ;;
        esac
        child+=(
            --id "$child_id"
            --libraries "$libraries"
            --output-dir "$run_directory"
        )

        if [[ $dry_run -eq 1 ]]; then
            printf 'Would run:'
            printf ' %q' "${child[@]}"
            printf '\n'
            append_auto_status "$status_file" "multiome" "$pair_label" "$child_id" \
                "dry_run" "0" "$run_directory" "$log_path" \
                "command printed without execution"
            continue
        fi

        {
            printf '%s Starting %s pair=%q\n' \
                "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$child_id" "$pair_label"
            printf 'Command:'
            printf ' %q' "${child[@]}"
            printf '\n'
        } | tee -a "$log_path"
        append_auto_status "$status_file" "multiome" "$pair_label" "$child_id" \
            "running" "" "$run_directory" "$log_path" "multiome count started"

        if (
            cd "$auto_output_root"
            "${child[@]}"
        ) 2>&1 | tee -a "$log_path"; then
            child_status=0
        else
            child_status=${PIPESTATUS[0]}
        fi

        if [[ $child_status -eq 0 ]] &&
            auto_output_is_complete "multiome" "$run_directory"; then
            append_auto_status "$status_file" "multiome" "$pair_label" "$child_id" \
                "complete" "0" "$run_directory" "$log_path" \
                "validated required output files"
            echo "Completed multiome pair: $child_id"
        else
            if [[ $child_status -eq 0 ]]; then
                child_status=3
                echo "ERROR: Required multiome outputs are missing for $child_id" |
                    tee -a "$log_path" >&2
            fi
            append_auto_status "$status_file" "multiome" "$pair_label" "$child_id" \
                "failed" "$child_status" "$run_directory" "$log_path" \
                "count failed or required outputs are missing"
            failures=$((failures + 1))
            echo "ERROR: Multiome pair $child_id failed with status $child_status" >&2
            if [[ $continue_on_error -eq 0 ]]; then
                echo "Stopping after the first failure. Re-run the same command to resume." >&2
                echo "Status log: $status_file" >&2
                return "$child_status"
            fi
        fi
    done < "$pairs_file"

    echo "Automatic multiome count pass finished. Status log: $status_file"
    if [[ $failures -gt 0 ]]; then
        return 1
    fi
    return 0
}

parsed=$(getopt \
    -o g:m:x:t:r:u:ac:nsSh \
    --long help,aggr,csv:,id:,run:,sample:,samples:,libraries:,gex_path:,gex-path:,atac_path:,atac-path:,normalize:,create-bam:,output-dir:,dry,strict-scan,continue-on-error \
    -n "$PROGRAM" -- "$@") || exit 2
eval set -- "$parsed"

while true; do
    case "$1" in
        -g) genome=$2; shift 2 ;;
        -m) mode=$2; shift 2 ;;
        -x) reference=$2; shift 2 ;;
        -t) threads=$2; shift 2 ;;
        -r) memory_gb=$2; shift 2 ;;
        -u) rna_version=$2; shift 2 ;;
        -a|--aggr) aggregate=1; shift ;;
        -c|--csv) aggregate_csv=$2; shift 2 ;;
        -n)
            normalize="legacy"
            shift
            ;;
        -s) secondary_mode="on"; shift ;;
        -S) secondary_mode="off"; shift ;;
        --id|--run) run_id=$2; shift 2 ;;
        --sample|--samples) sample=$2; shift 2 ;;
        --libraries) libraries_csv=$2; shift 2 ;;
        --gex_path|--gex-path) gex_path=$2; shift 2 ;;
        --atac_path|--atac-path) atac_path=$2; shift 2 ;;
        --normalize) normalize=$2; shift 2 ;;
        --create-bam) create_bam=$2; shift 2 ;;
        --output-dir) output_dir=$2; shift 2 ;;
        --dry) dry_run=1; shift ;;
        --strict-scan) strict_scan=1; shift ;;
        --continue-on-error) continue_on_error=1; shift ;;
        -h|--help) usage; exit 0 ;;
        --) shift; break ;;
        *) die "Internal option parser error at $1" ;;
    esac
done

positionals=("$@")
[[ -n $mode ]] || die "-m is required"
case "$mode" in
    rna|atac|multiome) ;;
    *) die "Unsupported mode: $mode" ;;
esac
require_positive_integer "threads" "$threads"
require_positive_integer "memory" "$memory_gb"

auto_count=0
multiome_path_mode=0
if [[ $aggregate -eq 1 ]]; then
    [[ -z $gex_path && -z $atac_path && -z $libraries_csv ]] ||
        die "--gex_path, --atac_path, and --libraries are count-only options"
    if [[ -n $aggregate_csv ]]; then
        [[ ${#positionals[@]} -eq 0 ]] ||
            die "Aggr with --csv does not accept COUNT_ROOT"
        [[ -f $aggregate_csv ]] ||
            die "Aggregation CSV does not exist: $aggregate_csv"
        aggregate_csv=$(realpath "$aggregate_csv")
    else
        [[ ${#positionals[@]} -le 1 ]] ||
            die "Aggr accepts at most one COUNT_ROOT when --csv is omitted"
        aggr_scan_root=${positionals[0]:-$(pwd -P)}
        [[ -d $aggr_scan_root ]] ||
            die "Aggregation count root does not exist: $aggr_scan_root"
        aggr_scan_root=$(realpath "$aggr_scan_root")
        generated_inputs_dir="$(pwd -P)"
        aggr_root_digest=$(printf '%s' "$aggr_scan_root" |
            sha256sum | awk '{print substr($1, 1, 8)}')
        aggregate_csv="$generated_inputs_dir/${mode}_aggregation_${aggr_root_digest}.csv"
        generate_aggr_csv_from_counts "$mode" "$aggr_scan_root" "$aggregate_csv"
        aggregate_csv=$(realpath "$aggregate_csv")
    fi
    if [[ -z $run_id ]]; then
        run_id=$(generate_csv_id "$mode" "aggr" "$aggregate_csv" "$normalize" "${reference:-$genome}")
        echo "Generated aggregation ID: $run_id"
    fi
elif [[ $mode == rna || $mode == atac ]]; then
    [[ -z $gex_path && -z $atac_path && -z $libraries_csv ]] ||
        die "--gex_path, --atac_path, and --libraries are multiome-only options"
    if [[ -z $run_id && -z $sample ]]; then
        auto_count=1
    else
        [[ -n $sample ]] ||
            die "--sample/--samples is required when --id/--run is provided"
        [[ ${#positionals[@]} -eq 1 ]] ||
            die "Explicit count mode requires exactly one FASTQ_DIRS argument"
        if [[ -z $run_id ]]; then
            run_id=$(generate_count_id "$mode" "$sample")
            echo "Generated count ID: $run_id"
        fi
    fi
else
    [[ ${#positionals[@]} -eq 0 ]] ||
        die "Multiome count does not accept a FASTQ positional argument"
    [[ -z $sample ]] || die "--sample/--samples is not used for multiome count"
    if [[ -n $libraries_csv ]]; then
        [[ -z $gex_path && -z $atac_path ]] ||
            die "Use either --libraries or --gex_path/--atac_path, not both"
        [[ -f $libraries_csv ]] ||
            die "Libraries CSV does not exist: $libraries_csv"
        libraries_csv=$(realpath "$libraries_csv")
        if [[ -z $run_id ]]; then
            run_id=$(generate_csv_id "$mode" "count" "$libraries_csv" "none" "${reference:-$genome}")
            echo "Generated count ID: $run_id"
        fi
    else
        [[ -n $gex_path && -n $atac_path ]] ||
            die "Multiome count requires --libraries or both --gex_path and --atac_path"
        [[ -d $gex_path ]] || die "GEX FASTQ root does not exist: $gex_path"
        [[ -d $atac_path ]] || die "ATAC FASTQ root does not exist: $atac_path"
        gex_path=$(realpath "$gex_path")
        atac_path=$(realpath "$atac_path")
        multiome_path_mode=1
    fi
fi

if [[ -n $run_id ]]; then
    [[ $run_id =~ ^[A-Za-z0-9_-]+$ ]] ||
        die "--id must contain only letters, digits, underscores, or hyphens"
    if [[ $mode != rna && ${#run_id} -gt 64 ]]; then
        die "--id exceeds the 64-character limit for $mode"
    fi
fi
case "$create_bam" in true|false) ;; *) die "--create-bam must be true or false" ;; esac

if [[ -z $reference ]]; then
    [[ -n $genome ]] || die "-g or -x is required"
    case "$mode:$genome" in
        rna:mm10) reference=$RNA_MM10_DEFAULT ;;
        rna:mm39) reference=$RNA_MM39_DEFAULT ;;
        rna:hg38) reference=$RNA_HG38_DEFAULT ;;
        atac:mm10|multiome:mm10) reference=$ARC_MM10_DEFAULT ;;
        atac:mm39|multiome:mm39) reference=$ARC_MM39_DEFAULT ;;
        atac:hg38|multiome:hg38) reference=$ARC_HG38_DEFAULT ;;
        *) die "Unsupported genome for $mode: $genome" ;;
    esac
fi
[[ -d $reference ]] || die "Reference directory does not exist: $reference"
reference=$(realpath "$reference")

case "$mode" in
    rna)
        case "$rna_version" in
            10) binary="cellranger" ;;
            9) binary="cellranger9" ;;
            8) binary="cellranger8" ;;
            7) binary="cellranger7" ;;
            *) die "RNA version must be 7, 8, 9, or 10" ;;
        esac
        ;;
    atac) binary="cellranger-atac" ;;
    multiome) binary="cellranger-arc" ;;
esac
command -v "$binary" >/dev/null 2>&1 || die "$binary is not available in PATH"

if [[ $normalize == legacy ]]; then
    if [[ $mode == rna ]]; then
        normalize="mapped"
    else
        normalize="depth"
    fi
    echo "WARNING: -n enables normalization; use --normalize none to retain full depth" >&2
fi
if [[ $mode == rna ]]; then
    [[ $normalize == none || $normalize == mapped ]] || die "RNA aggr normalization must be none or mapped"
else
    [[ $normalize == none || $normalize == depth ]] || die "$mode aggr normalization must be none or depth"
fi

if [[ $auto_count -eq 1 ]]; then
    if run_automatic_count; then
        exit 0
    else
        exit $?
    fi
fi
if [[ $multiome_path_mode -eq 1 ]]; then
    if run_multiome_paths; then
        exit 0
    else
        exit $?
    fi
fi

common_args=(--id "$run_id" --localcores "$threads" --localmem "$memory_gb" --disable-ui)
[[ -n $output_dir ]] && common_args+=(--output-dir "$output_dir")
[[ $dry_run -eq 1 ]] && common_args+=(--dry)

if [[ $secondary_mode == off || ($secondary_mode == default && $aggregate -eq 1) ]]; then
    common_args+=(--nosecondary)
fi

if [[ $aggregate -eq 1 ]]; then
    [[ -n $aggregate_csv ]] || die "--csv is required for aggr"
    [[ -f $aggregate_csv ]] || die "Aggregation CSV does not exist: $aggregate_csv"
    aggregate_csv=$(realpath "$aggregate_csv")
    cmd=("$binary" aggr "${common_args[@]}" --csv "$aggregate_csv" --normalize "$normalize")
    if [[ $mode != rna ]]; then
        cmd+=(--reference "$reference")
    fi
else
    case "$mode" in
        rna|atac)
            [[ ${#positionals[@]} -eq 1 ]] || die "Exactly one FASTQ_DIRS argument is required"
            [[ -n $sample ]] || die "--sample is required; one invocation must represent one biological sample"
            fastqs=${positionals[0]}
            IFS=',' read -r -a fastq_dirs <<< "$fastqs"
            resolved_dirs=()
            for directory in "${fastq_dirs[@]}"; do
                [[ -d $directory ]] || die "FASTQ directory does not exist: $directory"
                resolved_dirs+=("$(realpath "$directory")")
            done
            fastqs=$(IFS=,; echo "${resolved_dirs[*]}")
            cmd=("$binary" count "${common_args[@]}" --fastqs "$fastqs" --sample "$sample")
            if [[ $mode == rna ]]; then
                cmd+=(--transcriptome "$reference")
                if [[ $rna_version -ne 7 ]]; then
                    cmd+=(--create-bam "$create_bam")
                fi
            else
                cmd+=(--reference "$reference")
            fi
            ;;
        multiome)
            [[ ${#positionals[@]} -eq 0 ]] || die "Multiome count uses --libraries, not a FASTQ positional argument"
            [[ -n $libraries_csv ]] || die "--libraries is required for multiome count"
            [[ -f $libraries_csv ]] || die "Libraries CSV does not exist: $libraries_csv"
            cmd=("$binary" count "${common_args[@]}" --libraries "$(realpath "$libraries_csv")" --reference "$reference" --create-bam "$create_bam")
            ;;
    esac
fi

printf 'Running:'
printf ' %q' "${cmd[@]}"
printf '\n'
"${cmd[@]}"
status=$?
if [[ $status -eq 0 ]]; then
    echo "Cell Ranger command completed successfully."
else
    echo "Cell Ranger command failed with status $status." >&2
fi
exit "$status"
