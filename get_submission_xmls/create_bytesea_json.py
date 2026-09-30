#!/usr/bin/env python3

import argparse
import json
import re
import sys
from datetime import datetime
from pathlib import Path
from urllib.parse import quote
from zoneinfo import ZoneInfo

import requests

NGL_BI_BASE_URL = "http://ngl-bi.genoscope.cns.fr"
NGL_SQ_BASE_URL = "http://ngl-sq.genoscope.cns.fr"
REQUEST_HEADERS = {"User-Agent": "bot", "Accept": "application/json"}
REQUEST_TIMEOUT = 60
SUBMISSION_TIMEZONE = ZoneInfo("Europe/Paris")
BUSCO_RELATIVE_DIRECTORY = Path("eukaryota") / "Busco_geno_eukaryota"
BUSCO_FILENAME_PATTERN = "short_summary.specific.*.Busco_geno*.json"
GFASTATS_FIELDS = {
    "# scaffolds": ("nb_scaffolds", int),
    "Total scaffold length": ("size", int),
    "GC content %": ("gc_content", float),
    "Scaffold N50": ("n50_scaffolds", int),
    "Scaffold N90": ("n90_scaffolds", int),
    "Scaffold L50": ("l50_scaffolds", int),
    "Scaffold L90": ("l90_scaffolds", int),
    "# contigs": ("nb_contigs", int),
    "Contig N50": ("n50_contigs", int),
    "Contig N90": ("n90_contigs", int),
    "Contig L50": ("l50_contigs", int),
    "Contig L90": ("l90_contigs", int),
}
# Metrics without which the submission JSON is not usable; every other gfastats
# metric is reported as missing and left out of the output.
GFASTATS_MANDATORY_FIELDS = ("size", "gc_content")
# NGL-BI readset typeCodes producing long reads. Every other known typeCode is
# Illumina, where libProcessTypeCode tells Hi-C apart from whole-genome reads.
LONG_READ_READSET_TYPES = ("rspacbio", "rsnanopore")
ILLUMINA_READSET_TYPE = "rsillumina"
HIC_LIB_PROCESS_TYPE = "DF"
# Output order of the sequencing entries.
SEQUENCING_TYPE_ORDER = ("short_read", "long_read", "hi_c")

ANALYSIS_INCLUDES = (
    "sampleCodes",
    "readSetCodes",
    "masterReadSetCodes",
    "properties.umbrellaProjectAccession",
    "properties.sequencingProjectAccession",
    "properties.primaryAssemblyProjectAccession",
    "properties.tolid",
    "properties.assemblyToDownloadVersion",
    "treatments.reviewing.pairs.completion",
    "treatments.reviewing.pairs.merquryScore",
    "treatments.reviewing.pairs.scoreBuscoEuk",
    "treatments.reviewing.pairs.scoreBuscoTaxon",
    "treatments.reviewing.pairs.taxonBusco",
    "treatments.reviewing.pairs.resultDirectory",
)

READSET_INCLUDES = (
    "code",
    "typeCode",
    "sampleCode",
    "runSequencingStartDate",
    "sampleOnContainer.properties.libProcessTypeCode",
)

BUSCO_SCORE_PATTERN = re.compile(
    r"^C:\d+(?:\.\d+)?%\["
    r"S:(?P<s>\d+(?:\.\d+)?)%,"
    r"D:(?P<d>\d+(?:\.\d+)?)%\],"
    r"F:(?P<f>\d+(?:\.\d+)?)%,"
    r"M:(?P<m>\d+(?:\.\d+)?)%,"
    r"n:(?P<n>\d+)"
    r"(?:,E:\d+(?:\.\d+)?%)?$"
)


class SubmissionDataError(Exception):
    pass


def warn(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr)


def mandatory_pair_value(container: dict, key: str, source: str, mandatory_field: str):
    pair = container.get(key)
    if not isinstance(pair, dict) or pair.get("value") is None:
        raise SubmissionDataError(
            f"{source} is missing {key}.value; {mandatory_field} could not be "
            "retrieved and is mandatory in the submission JSON"
        )
    return pair["value"]


def optional_pair_value(container: dict, key: str, source: str, skipped: str):
    pair = container.get(key)
    if not isinstance(pair, dict) or pair.get("value") is None:
        warn(f"{source} is missing {key}.value; skipping {skipped}")
        return None
    return pair["value"]


def optional_pair_float(container: dict, key: str, source: str, skipped: str):
    value = optional_pair_value(container, key, source, skipped)
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        warn(f"{source} has an invalid {key}.value ({value!r}); skipping {skipped}")
        return None


def optional_busco_score(container: dict, key: str, source: str, skipped: str):
    value = optional_pair_value(container, key, source, skipped)
    if value is None:
        return None

    match = BUSCO_SCORE_PATTERN.fullmatch(value) if isinstance(value, str) else None
    if match is None:
        warn(f"{source} has an invalid {key}.value ({value!r}); skipping {skipped}")
        return None

    return {
        "s": float(match.group("s")),
        "d": float(match.group("d")),
        "m": float(match.group("m")),
        "f": float(match.group("f")),
        "nbgenes": int(match.group("n")),
    }


def fetch_analysis(session: requests.Session, project: str, material: str) -> dict:
    analysis_code = f"BA.{project}_{material}"
    url = f"{NGL_BI_BASE_URL}/api/analyses/{quote(analysis_code, safe='')}"
    params = [("includes", field) for field in ANALYSIS_INCLUDES]
    response = session.get(url, params=params, timeout=REQUEST_TIMEOUT)
    if response.status_code == 404:
        raise SubmissionDataError(f"NGL-BI analysis {analysis_code} was not found")
    response.raise_for_status()
    return response.json()


def extract_analysis_values(analysis: dict) -> dict:
    properties = analysis.get("properties")
    if not isinstance(properties, dict):
        raise SubmissionDataError("NGL-BI analysis is missing properties")

    treatments = analysis.get("treatments")
    reviewing = treatments.get("reviewing") if isinstance(treatments, dict) else None
    pairs = reviewing.get("pairs") if isinstance(reviewing, dict) else None
    if not isinstance(pairs, dict):
        warn("NGL-BI analysis is missing the reviewing treatment")
        pairs = {}

    sample_codes = analysis.get("sampleCodes")
    if not isinstance(sample_codes, list) or len(sample_codes) != 1:
        count = len(sample_codes) if isinstance(sample_codes, list) else 0
        raise SubmissionDataError(
            f"NGL-BI analysis must link exactly one sample; found {count}"
        )

    result_directory_pair = pairs.get("resultDirectory")
    result_directory = (
        result_directory_pair.get("value")
        if isinstance(result_directory_pair, dict)
        else None
    )
    reviewing_source = "NGL-BI reviewing treatment"

    readset_codes = analysis.get("readSetCodes")
    if not isinstance(readset_codes, list):
        readset_codes = []
    if not readset_codes:
        warn("NGL-BI analysis has no readset; skipping sequencing")

    master_readset_codes = analysis.get("masterReadSetCodes")
    if not isinstance(master_readset_codes, list):
        master_readset_codes = []
    if readset_codes and not master_readset_codes:
        warn("NGL-BI analysis has no master readset; no sequencing entry will be main")

    tolid = mandatory_pair_value(
        properties, "tolid", "NGL-BI properties", "ear_report"
    )

    return {
        "bioproject_umbrella": optional_pair_value(
            properties,
            "umbrellaProjectAccession",
            "NGL-BI properties",
            "bioproject_umbrella",
        ),
        "bioproject_assembly": mandatory_pair_value(
            properties,
            "primaryAssemblyProjectAccession",
            "NGL-BI properties",
            "bioproject_assembly",
        ),
        "bioproject_reads": mandatory_pair_value(
            properties,
            "sequencingProjectAccession",
            "NGL-BI properties",
            "bioproject_reads",
        ),
        "assembly_name": mandatory_pair_value(
            properties,
            "assemblyToDownloadVersion",
            "NGL-BI properties",
            "assembly_name",
        ),
        "ear_report": f"EARs/{tolid}_EAR.pdf",
        "sample_code": sample_codes[0],
        "readset_codes": readset_codes,
        "master_readset_codes": master_readset_codes,
        "busco_euk": optional_busco_score(
            pairs, "scoreBuscoEuk", reviewing_source, "the busco_euk_* fields"
        ),
        "busco_lin": optional_busco_score(
            pairs, "scoreBuscoTaxon", reviewing_source, "the busco_lin_* fields"
        ),
        "busco_lin_dataset": optional_pair_value(
            pairs, "taxonBusco", reviewing_source, "busco_lin_dataset"
        ),
        "merqury_completion": optional_pair_float(
            pairs, "completion", reviewing_source, "merqury_completion"
        ),
        "merqury_score": optional_pair_float(
            pairs, "merquryScore", reviewing_source, "merqury_score"
        ),
        "result_directory": result_directory,
    }


def fetch_readsets(session: requests.Session, readset_codes: list) -> list:
    if not readset_codes:
        return []

    params = [("codes", code) for code in readset_codes]
    params += [("includes", field) for field in READSET_INCLUDES]
    response = session.get(
        f"{NGL_BI_BASE_URL}/api/readsets", params=params, timeout=REQUEST_TIMEOUT
    )
    response.raise_for_status()
    readsets = response.json()
    if not isinstance(readsets, list):
        raise SubmissionDataError("NGL-BI readset query did not return a list")

    missing = set(readset_codes) - {readset.get("code") for readset in readsets}
    if missing:
        warn(f"NGL-BI did not return readset(s) {', '.join(sorted(missing))}")
    return readsets


def fetch_specimen_codes(session: requests.Session, sample_codes: list) -> dict:
    """Map each NGL-SQ sample code to its specimen code ("Code unique individu")."""
    params = [("codes", code) for code in sample_codes]
    params += [("includes", "code"), ("includes", "properties.individualNumber")]
    response = session.get(
        f"{NGL_SQ_BASE_URL}/api/samples", params=params, timeout=REQUEST_TIMEOUT
    )
    response.raise_for_status()
    samples = response.json()
    if not isinstance(samples, list):
        raise SubmissionDataError("NGL-SQ sample query did not return a list")

    missing = set(sample_codes) - {sample.get("code") for sample in samples}
    if missing:
        warn(f"NGL-SQ did not return sample(s) {', '.join(sorted(missing))}")

    specimens = {}
    for sample in samples:
        properties = sample.get("properties")
        pair = (
            properties.get("individualNumber") if isinstance(properties, dict) else None
        )
        value = pair.get("value") if isinstance(pair, dict) else None
        if value is not None:
            specimens[sample.get("code")] = value
    return specimens


def resolve_busco_path(override: Path | None, result_directory: str | None) -> Path | None:
    skipped = "busco_version and busco_dataset_version"
    if override is not None:
        if not override.is_file():
            raise SubmissionDataError(f"BUSCO log was not found: {override}")
        return override

    if not result_directory:
        warn(
            "NGL-BI reviewing treatment is missing resultDirectory.value; "
            f"skipping {skipped}"
        )
        return None

    busco_directory = Path(result_directory) / BUSCO_RELATIVE_DIRECTORY
    matches = list(busco_directory.glob(BUSCO_FILENAME_PATTERN))
    if not matches:
        warn(
            f"no BUSCO log matching {BUSCO_FILENAME_PATTERN} in {busco_directory}; "
            f"skipping {skipped}"
        )
        return None
    if len(matches) > 1:
        warn(
            f"multiple BUSCO logs matching {BUSCO_FILENAME_PATTERN} in "
            f"{busco_directory}; skipping {skipped}"
        )
        return None
    return matches[0]


def parse_busco_log(path: Path) -> dict:
    try:
        with path.open(encoding="utf-8") as busco_file:
            log = json.load(busco_file)
    except (OSError, ValueError) as error:
        warn(
            f"BUSCO log {path} could not be read ({error}); "
            "skipping busco_version and busco_dataset_version"
        )
        return {}

    versions = log.get("versions") if isinstance(log, dict) else None
    version = versions.get("busco") if isinstance(versions, dict) else None
    if not isinstance(version, str) or not version:
        warn(f"BUSCO log is missing versions.busco ({path}); skipping busco_version")
        version = None

    parameters = log.get("parameters") if isinstance(log, dict) else None
    dataset_version = (
        parameters.get("datasets_version") if isinstance(parameters, dict) else None
    )
    if not isinstance(dataset_version, str) or not dataset_version:
        warn(
            f"BUSCO log is missing parameters.datasets_version ({path}); "
            "skipping busco_dataset_version"
        )
        dataset_version = None

    return {"version": version, "dataset_version": dataset_version}


def resolve_gfastats_path(
    override: Path | None,
    result_directory: str | None,
    lineage_dataset: str | None,
) -> Path:
    if override is not None:
        if not override.is_file():
            raise SubmissionDataError(f"gfastats file was not found: {override}")
        return override

    if not result_directory:
        raise SubmissionDataError(
            "NGL-BI reviewing treatment is missing resultDirectory.value"
        )
    if not lineage_dataset:
        raise SubmissionDataError("NGL-BI reviewing treatment is missing taxonBusco.value")

    path = (
        Path(result_directory)
        / lineage_dataset
        / "Results"
        / "gfastats_assembly.txt"
    )
    if not path.is_file():
        raise SubmissionDataError(f"gfastats file was not found: {path}")
    return path


def parse_gfastats(path: Path) -> dict:
    metrics = {}
    invalid_fields = set()
    with path.open(encoding="utf-8") as gfastats_file:
        for line in gfastats_file:
            label, separator, raw_value = line.strip().partition(":")
            if not separator or label not in GFASTATS_FIELDS:
                continue

            output_field, converter = GFASTATS_FIELDS[label]
            value = raw_value.strip()
            try:
                metrics[output_field] = converter(value)
            except ValueError as error:
                message = f"invalid {label} value in gfastats file {path}: {value!r}"
                if output_field in GFASTATS_MANDATORY_FIELDS:
                    raise SubmissionDataError(message) from error
                warn(f"{message}; skipping {output_field}")
                invalid_fields.add(output_field)

    missing_mandatory = [
        output_field
        for output_field in GFASTATS_MANDATORY_FIELDS
        if output_field not in metrics
    ]
    if missing_mandatory:
        raise SubmissionDataError(
            f"gfastats file is missing required metric(s): {', '.join(missing_mandatory)}"
        )

    for output_field, _ in GFASTATS_FIELDS.values():
        if output_field not in metrics and output_field not in invalid_fields:
            warn(f"gfastats file {path} is missing {output_field}; skipping it")
    return metrics


def parse_manifest(path: Path) -> dict:
    fields = {}
    with path.open(encoding="utf-8") as manifest_file:
        for line in manifest_file:
            parts = line.strip().split(None, 1)
            if len(parts) == 2:
                fields[parts[0]] = parts[1]
    return fields


def extract_ena_sample_code(manifest: dict, assembly_name: str) -> str | None:
    """Read the manifest SAMPLE field, cross-checking ASSEMBLYNAME against NGL-BI."""
    manifest_assembly_name = manifest.get("ASSEMBLYNAME")
    if manifest_assembly_name and manifest_assembly_name != assembly_name:
        warn(
            f"manifest ASSEMBLYNAME ({manifest_assembly_name}) does not match the "
            f"NGL-BI assembly name ({assembly_name}); NGL-BI wins"
        )

    ena_sample_code = manifest.get("SAMPLE")
    if not ena_sample_code:
        warn("manifest is missing SAMPLE; skipping ena_sample_code")
        return None
    return ena_sample_code


def build_submission_json(
    ngl: dict,
    main_specimen_code: str,
    sequencing: list,
    ena_sample_code: str | None,
    busco: dict,
    gfastats: dict,
) -> dict:
    euk = ngl["busco_euk"] or {}
    lineage = ngl["busco_lin"] or {}

    submission = {
        "main": True,
        "submission_date": datetime.now(SUBMISSION_TIMEZONE).isoformat(
            timespec="milliseconds"
        ),
        "assembly_level": "chromosome",
        "bioproject_umbrella": ngl["bioproject_umbrella"],
        "bioproject_assembly": ngl["bioproject_assembly"],
        "bioproject_reads": ngl["bioproject_reads"],
        "atlaseaid_prefix": "",
        "main_specimen_code": main_specimen_code,
        "assembly_accession": "",
        "busco_euk_genome_s": euk.get("s"),
        "busco_euk_genome_d": euk.get("d"),
        "busco_euk_genome_m": euk.get("m"),
        "busco_euk_genome_f": euk.get("f"),
        "busco_euk_dataset": "eukaryota" if euk else None,
        "busco_euk_nbgenes": euk.get("nbgenes"),
        "busco_lin_genome_s": lineage.get("s"),
        "busco_lin_genome_d": lineage.get("d"),
        "busco_lin_genome_m": lineage.get("m"),
        "busco_lin_genome_f": lineage.get("f"),
        "busco_lin_dataset": ngl["busco_lin_dataset"],
        "busco_lin_nbgenes": lineage.get("nbgenes"),
        "merqury_completion": ngl["merqury_completion"],
        "merqury_score": ngl["merqury_score"],
        "ena_sample_code": ena_sample_code,
        "assembly_name": ngl["assembly_name"],
        "ear_report": ngl["ear_report"],
        "sequencing": sequencing or None,
        "busco_version": busco.get("version"),
        "busco_dataset_version": busco.get("dataset_version"),
        "nb_scaffolds": gfastats.get("nb_scaffolds"),
        "size": gfastats["size"],
        "gc_content": gfastats["gc_content"],
        "n50_scaffolds": gfastats.get("n50_scaffolds"),
        "n90_scaffolds": gfastats.get("n90_scaffolds"),
        "l50_scaffolds": gfastats.get("l50_scaffolds"),
        "l90_scaffolds": gfastats.get("l90_scaffolds"),
        "nb_contigs": gfastats.get("nb_contigs"),
        "n50_contigs": gfastats.get("n50_contigs"),
        "n90_contigs": gfastats.get("n90_contigs"),
        "l50_contigs": gfastats.get("l50_contigs"),
        "l90_contigs": gfastats.get("l90_contigs"),
    }
    return {key: value for key, value in submission.items() if value is not None}


def readset_sequencing_type(readset: dict) -> str | None:
    type_code = readset.get("typeCode")
    if type_code in LONG_READ_READSET_TYPES:
        return "long_read"
    if type_code != ILLUMINA_READSET_TYPE:
        warn(
            f"readset {readset.get('code')} has an unknown typeCode "
            f"({type_code!r}); skipping it"
        )
        return None

    container = readset.get("sampleOnContainer")
    properties = container.get("properties") if isinstance(container, dict) else None
    pair = (
        properties.get("libProcessTypeCode") if isinstance(properties, dict) else None
    )
    process_type = pair.get("value") if isinstance(pair, dict) else None
    return "hi_c" if process_type == HIC_LIB_PROCESS_TYPE else "short_read"


def build_sequencing(
    readsets: list, master_readset_codes: list, specimens: dict
) -> list:
    """Build one sequencing entry per sequencing type and specimen."""
    masters = set(master_readset_codes)
    groups = {}
    for readset in readsets:
        sequencing_type = readset_sequencing_type(readset)
        if sequencing_type is None:
            continue

        sample_code = readset.get("sampleCode")
        specimen = specimens.get(sample_code)
        if specimen is None:
            warn(
                f"readset {readset.get('code')} sample {sample_code} has no "
                "specimen code; skipping it"
            )
            continue

        group = groups.setdefault(
            (sequencing_type, specimen), {"dates": [], "main": False}
        )
        run_date = readset.get("runSequencingStartDate")
        if run_date is None:
            warn(f"readset {readset.get('code')} is missing runSequencingStartDate")
        else:
            group["dates"].append(run_date)
        group["main"] = group["main"] or readset.get("code") in masters

    sequencing = []
    for key in sorted(groups, key=lambda k: (SEQUENCING_TYPE_ORDER.index(k[0]), k[1])):
        sequencing_type, specimen = key
        group = groups[key]
        entry = {"sequencing_type": sequencing_type}
        if group["dates"]:
            # NGL stores run start dates as a local midnight; drop the offset to
            # keep the rendered value a plain date and time.
            start = datetime.fromtimestamp(
                min(group["dates"]) / 1000, SUBMISSION_TIMEZONE
            )
            entry["creation_date"] = start.replace(tzinfo=None).isoformat(
                timespec="seconds"
            )
        entry["main"] = group["main"]
        entry["specimens"] = [specimen]
        sequencing.append(entry)
    return sequencing


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate a submission JSON file from NGL-BI and NGL-SQ data."
    )
    parser.add_argument("--project", required=True, help="NGL project code")
    parser.add_argument("--material", required=True, help="NGL material code")
    parser.add_argument(
        "--manifest",
        type=Path,
        help="Path to the ENA manifest file; only its SAMPLE field is used",
    )
    parser.add_argument(
        "--busco",
        type=Path,
        help=(
            "Path to a BUSCO JSON log; defaults to the reviewing result directory"
        ),
    )
    parser.add_argument(
        "--gfastats",
        type=Path,
        help=(
            "Path to gfastats output; defaults to the lineage dataset directory"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("submission.json"),
        help="Path to write (default: submission.json)",
    )
    args = parser.parse_args()

    try:
        with requests.Session() as session:
            session.headers.update(REQUEST_HEADERS)
            analysis = fetch_analysis(session, args.project, args.material)
            ngl = extract_analysis_values(analysis)
            readsets = fetch_readsets(session, ngl["readset_codes"])
            sample_codes = {ngl["sample_code"]} | {
                readset["sampleCode"]
                for readset in readsets
                if readset.get("sampleCode")
            }
            specimens = fetch_specimen_codes(session, sorted(sample_codes))
        if args.manifest is None:
            warn("no manifest provided; skipping ena_sample_code")
            ena_sample_code = None
        else:
            ena_sample_code = extract_ena_sample_code(
                parse_manifest(args.manifest), ngl["assembly_name"]
            )
        main_specimen_code = specimens.get(ngl["sample_code"])
        if main_specimen_code is None:
            raise SubmissionDataError(
                f"NGL-SQ sample {ngl['sample_code']} is missing "
                "properties.individualNumber.value"
            )
        sequencing = build_sequencing(
            readsets, ngl["master_readset_codes"], specimens
        )
        busco_path = resolve_busco_path(args.busco, ngl["result_directory"])
        busco = parse_busco_log(busco_path) if busco_path is not None else {}
        gfastats_path = resolve_gfastats_path(
            args.gfastats, ngl["result_directory"], ngl["busco_lin_dataset"]
        )
        gfastats = parse_gfastats(gfastats_path)
        submission = build_submission_json(
            ngl, main_specimen_code, sequencing, ena_sample_code, busco, gfastats
        )
        with args.output.open("w", encoding="utf-8") as output_file:
            json.dump(submission, output_file, indent=4)
            output_file.write("\n")
    except (OSError, requests.RequestException, SubmissionDataError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print(f"Submission JSON written to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
