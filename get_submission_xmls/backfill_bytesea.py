#!/usr/bin/env python3
"""Backfill BYTESea with genomes already submitted to the ENA.

For every assembly of the ENA umbrella study (latest version only), find the
NGL-BI analysis of project DKC by ToLID, upload its EAR report, generate a manifest whose
SAMPLE is the ENA sample of the assembly, and write the BYTESea JSON to the
pending directory. Genomes without an EAR report are skipped with a warning.
"""
import argparse
import csv
import gzip
import io
import os
import subprocess
import sys
import tempfile

import requests

from submit_genome import BYTESEA_PENDING_DIRECTORY, write_pending_bytesea_json

ENA_FILEREPORT_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
NGL_ANALYSES_URL = "http://ngl-bi.genoscope.cns.fr/api/analyses"
NGL_PROJECT = "DKC"
DEFAULT_STUDY = "PRJEB64126"


def warn(message):
    print(f"WARNING: {message}", file=sys.stderr)


def fetch_ena_assemblies(study):
    """Return {tolid: (assembly_name, ENA sample accession)} for the latest
    version of each primary assembly."""
    params = {
        "accession": study,
        "result": "assembly",
        "format": "tsv",
        "fields": "assembly_name,version,scientific_name,secondary_sample_accession",
    }
    res = requests.get(ENA_FILEREPORT_URL, params=params, timeout=120)
    res.raise_for_status()
    latest = {}
    for row in csv.DictReader(io.StringIO(res.text), delimiter="\t"):
        accession = row["accession"]
        if accession not in latest or int(row["version"]) > int(latest[accession]["version"]):
            latest[accession] = row

    assemblies = {}
    for row in latest.values():
        name = row["assembly_name"]
        # Alternate haplotypes share the primary's ToLID and sample; ENA names
        # them inconsistently ("alternate haplotype", "_AlternateHaplotype", ...).
        if "alternat" in name.lower():
            continue
        # assembly_name is "<tolid>.<version>", e.g. uyTelAmph1.1
        tolid = name.rsplit(".", 1)[0]
        assemblies[tolid] = (name, row["secondary_sample_accession"])
    return assemblies


def fetch_ngl_analyses():
    """Return {tolid: [(ba_code, ear_report_path or None), ...]} for the DKC analyses."""
    params = [
        ("projectCodes", NGL_PROJECT),
        ("includes", "code"),
        ("includes", "properties.tolid"),
        ("includes", "treatments.reviewing.pairs.earReport"),
    ]
    res = requests.get(
        NGL_ANALYSES_URL, params=params, headers={"User-Agent": "bot"}, timeout=600
    )
    res.raise_for_status()
    analyses = {}
    for analysis in res.json():
        tolid = ((analysis.get("properties") or {}).get("tolid") or {}).get("value")
        if not tolid:
            continue
        reviewing = (analysis.get("treatments") or {}).get("reviewing") or {}
        ear_path = ((reviewing.get("pairs") or {}).get("earReport") or {}).get("value")
        analyses.setdefault(tolid, []).append((analysis["code"], ear_path))
    return analyses


def set_manifest_fields(manifest_path, fields):
    """Replace the given key/value lines of a generated manifest."""
    with open(manifest_path) as f:
        lines = f.readlines()
    with open(manifest_path, "w") as f:
        for line in lines:
            key = line.split("\t", 1)[0]
            f.write(f"{key}\t{fields[key]}\n" if key in fields else line)


def backfill_genome(tolid, assembly_name, sample, ba_code, ear_path, workdir):
    """Process one genome; return True if its JSON was written to pending."""
    # BA.<PROJECT>_<MATERIAL>; project codes contain no underscore.
    project, material = ba_code.removeprefix("BA.").split("_", 1)

    if not ear_path or not os.path.isfile(ear_path):
        warn(f"{tolid} ({ba_code}) has no EAR report ({ear_path or 'not set in NGL-BI'}), skipping")
        return False

    manifest_path = os.path.join(workdir, f"{tolid}.manifest.txt")
    # Already submitted: the FASTA is never read, but generate_manifest requires
    # an existing gzipped file. An empty gzip stream satisfies its check.
    placeholder_fasta = os.path.join(workdir, f"{tolid}.fa.gz")
    gzip.open(placeholder_fasta, "wb").close()
    manifest_cmd = [
        "generate_manifest",
        "--project", project,
        "--material", material,
        "--fasta", placeholder_fasta,
        "--output", manifest_path,
    ]
    if subprocess.run(manifest_cmd).returncode != 0:
        warn(f"{tolid} ({ba_code}): manifest generation failed, skipping")
        return False
    # generate_manifest always writes "<tolid>.1"; use what ENA actually holds.
    set_manifest_fields(manifest_path, {"SAMPLE": sample, "ASSEMBLYNAME": assembly_name})

    ear_cmd = ["submit_ear", "--project", project, "--material", material, "--ear", ear_path]
    if subprocess.run(ear_cmd).returncode != 0:
        warn(f"{tolid} ({ba_code}): EAR upload failed, skipping")
        return False

    if not write_pending_bytesea_json(project, material, manifest_path, tolid):
        warn(f"{tolid} ({ba_code}): BYTESea JSON creation failed")
        return False
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Backfill BYTESea with the genomes already submitted to the ENA."
    )
    parser.add_argument(
        "--study",
        default=DEFAULT_STUDY,
        help=f"ENA study listing the assemblies (default: {DEFAULT_STUDY})",
    )
    parser.add_argument(
        "--tolid",
        action="append",
        help="Only process this ToLID (repeatable); default: every assembly of the study",
    )
    args = parser.parse_args()

    assemblies = fetch_ena_assemblies(args.study)
    analyses = fetch_ngl_analyses()
    tolids = args.tolid or sorted(assemblies)

    written, skipped = [], []
    with tempfile.TemporaryDirectory() as workdir:
        for tolid in tolids:
            print(f"=== {tolid} ===", file=sys.stderr)
            if tolid not in assemblies:
                warn(f"{tolid} is not an assembly of {args.study}, skipping")
                skipped.append(tolid)
                continue
            matches = analyses.get(tolid, [])
            if not matches:
                warn(f"{tolid} not found in NGL-BI project {NGL_PROJECT}, skipping")
                skipped.append(tolid)
                continue
            if len(matches) > 1:
                codes = ", ".join(code for code, _ in matches)
                warn(f"{tolid}: several {NGL_PROJECT} analyses ({codes}), skipping")
                skipped.append(tolid)
                continue
            ba_code, ear_path = matches[0]
            assembly_name, sample = assemblies[tolid]
            if backfill_genome(tolid, assembly_name, sample, ba_code, ear_path, workdir):
                written.append(tolid)
            else:
                skipped.append(tolid)

    print(
        f"\n{len(written)} JSON written to {BYTESEA_PENDING_DIRECTORY}, "
        f"{len(skipped)} skipped",
        file=sys.stderr,
    )
    if skipped:
        print(f"Skipped: {' '.join(skipped)}", file=sys.stderr)
    return 0 if not skipped else 1


if __name__ == "__main__":
    sys.exit(main())
