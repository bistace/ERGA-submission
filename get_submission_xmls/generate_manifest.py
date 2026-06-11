#!/usr/bin/env python3
import argparse
import csv
import glob
import requests
import sys


# Manifest fields that are constant for every ATLASea genome submission.
STATIC_FIELDS = {
    "ASSEMBLY_TYPE": "clone or isolate",
    "MINGAPLENGTH": "100",
    "MOLECULETYPE": "genomic DNA",
    "CHROMOSOME_LIST": "chr_list.txt.gz",
    "UNLOCALISED_LIST": "unloc_list.txt.gz",
}

# Downstream programs in the PROGRAM field, appended after the assembler.
DOWNSTREAM_PROGRAMS = "Purge_dups,Yahs"

# Assembler detected from the selectedAssembly path in the Assembly
# Decontamination treatment. Each entry maps a lowercase substring to its
# manifest label; the first match (in order) wins.
ASSEMBLERS = [("flye", "Flye"), ("hifiasm", "Hifiasm"), ("nextdenovo", "Nextdenovo")]
DEFAULT_ASSEMBLER = "Hifiasm"

# PLATFORM is built from the long-read technology plus the Hi-C kit(s) used,
# e.g. "PacBio,Arima", "ONT,OmniC" or "PacBio,Arima,OmniC".
LONG_READ_ONT = "ONT"
LONG_READ_PACBIO = "PacBio"

# NGL-BI readset typeCode that identifies an Oxford Nanopore readset.
ONT_TYPECODE = "rsnanopore"

# Hi-C kits, keyed by their hicLibqc "pairs" key, in manifest output order.
HIC_KITS = [("arima", "Arima"), ("omnic", "OmniC")]


def main():
    parser = argparse.ArgumentParser(description="Generate an ENA manifest file from NGL-BI data.")
    parser.add_argument("--project", required=True, help="Project code")
    parser.add_argument("--material", required=True, help="Material code")
    parser.add_argument(
        "--fasta", required=True, help="Assembly FASTA file (e.g. laThaTest1.fa.gz)"
    )
    parser.add_argument(
        "--output",
        default="manifest.txt",
        help="Path to the manifest file to write (default: manifest.txt)",
    )

    args = parser.parse_args()

    assembly_json = get_assembly_json(args.project, args.material)
    if assembly_json is None:
        sys.exit(1)

    study, assembly_name, coverage = extract_fields(assembly_json)
    platform = determine_platform(assembly_json)
    program = determine_program(assembly_json)
    description = get_description(args.project, args.material)

    manifest = build_manifest(
        study, assembly_name, coverage, program, platform, description, args.fasta
    )
    write_manifest(manifest, args.output)

    print(f"Manifest written to {args.output}", file=sys.stderr)
    print(f"  STUDY:        {study}", file=sys.stderr)
    print(f"  ASSEMBLYNAME: {assembly_name}", file=sys.stderr)
    print(f"  COVERAGE:     {coverage}", file=sys.stderr)
    print(f"  PROGRAM:      {program}", file=sys.stderr)
    print(f"  PLATFORM:     {platform}", file=sys.stderr)
    print(f"  DESCRIPTION:  {'set from EAR data' if description else 'blank'}", file=sys.stderr)
    print("  SAMPLE left blank for manual entry.", file=sys.stderr)


# --- NGL-BI API ---


def get_assembly_json(project_code: str, material_code: str) -> dict:
    """Fetch assembly data from the NGL-BI API."""
    headers = {"User-Agent": "bot"}
    url = f"http://ngl-bi.genoscope.cns.fr/api/analyses/BA.{project_code}_{material_code}"
    res = requests.get(url, headers=headers)
    if res.status_code == 404:
        print(
            f"ERROR: Project {project_code}/{material_code} has no assembly treatment!",
            file=sys.stderr,
        )
        return None
    res.raise_for_status()
    return res.json()


def get_readset_typecode(readset_code: str) -> str:
    """Fetch a single readset from NGL-BI and return its typeCode (or "" if unavailable)."""
    headers = {"User-Agent": "bot"}
    url = f"http://ngl-bi.genoscope.cns.fr/api/readsets/{readset_code}"
    res = requests.get(url, headers=headers)
    if res.status_code == 404:
        print(f"WARNING: readset {readset_code} not found in NGL-BI", file=sys.stderr)
        return ""
    res.raise_for_status()
    return res.json().get("typeCode", "")


def determine_long_read(assembly: dict) -> str:
    """Return the long-read technology: ONT if any nanopore readset is used, else PacBio."""
    readset_codes = assembly.get("readSetCodes", [])
    if not readset_codes:
        print(
            f"WARNING: no readsets found in NGL-BI, defaulting long-read to {LONG_READ_PACBIO}",
            file=sys.stderr,
        )
        return LONG_READ_PACBIO

    has_ont = any(get_readset_typecode(rs) == ONT_TYPECODE for rs in readset_codes)
    return LONG_READ_ONT if has_ont else LONG_READ_PACBIO


def determine_hic_kits(assembly: dict) -> list:
    """Return the Hi-C kit(s) (Arima and/or OmniC) found in the Hi-C LibQC treatment."""
    hic = assembly.get("treatments", {}).get("hicLibqc")
    if not hic:
        print(
            "WARNING: hicLibqc treatment not found in NGL-BI, no Hi-C kit added to PLATFORM",
            file=sys.stderr,
        )
        return []

    pairs = hic.get("pairs", {})
    kits = [label for key, label in HIC_KITS if pairs.get(key, {}).get("value")]
    if not kits:
        print(
            "WARNING: no Arima or OmniC entries in hicLibqc treatment, no Hi-C kit added to PLATFORM",
            file=sys.stderr,
        )
    return kits


def determine_platform(assembly: dict) -> str:
    """Build the PLATFORM value from the long-read technology and the Hi-C kit(s)."""
    components = [determine_long_read(assembly)] + determine_hic_kits(assembly)
    return ",".join(components)


def determine_assembler(assembly: dict) -> str:
    """Detect the assembler from the selectedAssembly path in the Assembly Decontamination treatment.

    The selectedAssembly value is a path to the chosen assembly FASTA, e.g.
    ".../best_formatted/Flye_fmt.fasta". The assembler name is matched
    case-insensitively against the known assemblers. Falls back to
    DEFAULT_ASSEMBLER with a warning if the field is missing or unrecognised.
    """
    decontamination = assembly.get("treatments", {}).get("assemblyDecontamination")
    if not decontamination:
        print(
            "WARNING: assemblyDecontamination treatment not found in NGL-BI, "
            f"defaulting assembler to {DEFAULT_ASSEMBLER}",
            file=sys.stderr,
        )
        return DEFAULT_ASSEMBLER

    selected = decontamination.get("pairs", {}).get("selectedAssembly", {}).get("value")
    if not selected:
        print(
            "WARNING: selectedAssembly not found in Assembly Decontamination treatment, "
            f"defaulting assembler to {DEFAULT_ASSEMBLER}",
            file=sys.stderr,
        )
        return DEFAULT_ASSEMBLER

    selected_lower = selected.lower()
    for needle, label in ASSEMBLERS:
        if needle in selected_lower:
            return label

    print(
        f"WARNING: could not detect a known assembler in selectedAssembly '{selected}', "
        f"defaulting assembler to {DEFAULT_ASSEMBLER}",
        file=sys.stderr,
    )
    return DEFAULT_ASSEMBLER


def determine_program(assembly: dict) -> str:
    """Build the PROGRAM value from the detected assembler and the downstream programs."""
    return f"{determine_assembler(assembly)},{DOWNSTREAM_PROGRAMS}"


def extract_fields(assembly: dict):
    """Extract study accession, assembly name and LR coverage from NGL-BI.

    Any field that cannot be found is left blank and a warning is emitted.
    """
    properties = assembly.get("properties", {})

    study = properties.get("primaryAssemblyProjectAccession", {}).get("value")
    if not study:
        print(
            "WARNING: primaryAssemblyProjectAccession not found in NGL-BI, leaving STUDY blank",
            file=sys.stderr,
        )
        study = ""

    tolid = properties.get("tolid", {}).get("value")
    if not tolid:
        print("WARNING: tolid not found in NGL-BI, leaving ASSEMBLYNAME blank", file=sys.stderr)
        assembly_name = ""
    else:
        assembly_name = f"{tolid}.1"

    coverage = extract_lr_coverage(assembly)

    return study, assembly_name, coverage


def extract_lr_coverage(assembly: dict) -> str:
    """Extract and round the LR reads coverage from the Genome Contigs Scaffolding treatment.

    Returns the rounded coverage as a string, or "" with a warning if not found.
    """
    treatments = assembly.get("treatments", {})
    scaffolding = treatments.get("genomeContigsScaffolding")
    if not scaffolding:
        print(
            "WARNING: genomeContigsScaffolding treatment not found in NGL-BI, leaving COVERAGE blank",
            file=sys.stderr,
        )
        return ""

    lr_reads = scaffolding.get("pairs", {}).get("lr_reads", {}).get("value")
    if not lr_reads:
        print(
            "WARNING: LR reads coverage not found in Genome Contigs Scaffolding treatment, leaving COVERAGE blank",
            file=sys.stderr,
        )
        return ""

    coverage = lr_reads[0].get("coverage")
    if coverage is None:
        print(
            "WARNING: coverage value missing for LR reads, leaving COVERAGE blank", file=sys.stderr
        )
        return ""

    return str(round(coverage))


# --- EAR data description ---


def get_description(project_code: str, material_code: str) -> str:
    """Read the curation notes from the EAR data CSV produced during curation.

    Returns the curation_notes text, or "" with a warning if the file cannot be
    reached or the field is missing.
    """
    pattern = (
        f"/env/cns/proj/projet_{project_code}/scratch/CORRECTED_SCAFFOLDING/"
        f"{project_code}_{material_code}_*/review/"
        f"{project_code}_{material_code}_review_ear_data.csv"
    )
    matches = glob.glob(pattern)
    if not matches:
        print(
            f"WARNING: could not reach EAR data file matching {pattern}, leaving DESCRIPTION blank",
            file=sys.stderr,
        )
        return ""
    if len(matches) > 1:
        print(
            f"WARNING: multiple EAR data files match {pattern}, using {matches[0]}", file=sys.stderr
        )

    path = matches[0]
    try:
        with open(path, newline="") as f:
            for row in csv.reader(f, delimiter=";"):
                if len(row) >= 2 and row[0] == "curation_notes":
                    return row[1].strip()
    except OSError as e:
        print(f"WARNING: could not read {path}: {e}, leaving DESCRIPTION blank", file=sys.stderr)
        return ""

    print(
        f"WARNING: curation_notes not found in {path}, leaving DESCRIPTION blank", file=sys.stderr
    )
    return ""


# --- Manifest writing ---


def build_manifest(
    study: str,
    assembly_name: str,
    coverage: str,
    program: str,
    platform: str,
    description: str,
    fasta: str,
) -> list:
    """Build the ordered list of (key, value) manifest lines."""
    return [
        ("STUDY", study),
        ("SAMPLE", ""),
        ("ASSEMBLYNAME", assembly_name),
        ("ASSEMBLY_TYPE", STATIC_FIELDS["ASSEMBLY_TYPE"]),
        ("COVERAGE", coverage),
        ("PROGRAM", program),
        ("PLATFORM", platform),
        ("MINGAPLENGTH", STATIC_FIELDS["MINGAPLENGTH"]),
        ("MOLECULETYPE", STATIC_FIELDS["MOLECULETYPE"]),
        ("DESCRIPTION", description),
        ("RUN_REF", ""),
        ("FASTA", fasta),
        ("CHROMOSOME_LIST", STATIC_FIELDS["CHROMOSOME_LIST"]),
        ("UNLOCALISED_LIST", STATIC_FIELDS["UNLOCALISED_LIST"]),
    ]


def write_manifest(manifest: list, path: str):
    """Write manifest lines as tab-separated key/value pairs."""
    with open(path, "w") as f:
        for key, value in manifest:
            f.write(f"{key}\t{value}\n")


if __name__ == "__main__":
    main()
