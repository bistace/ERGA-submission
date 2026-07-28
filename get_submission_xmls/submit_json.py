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
SEQUENCING_PLATFORMS = ("PacBio", "ONT", "Arima", "OmniC")
LONG_READ_PLATFORMS = ("PacBio", "ONT")

ANALYSIS_INCLUDES = (
    "sampleCodes",
    "properties.umbrellaProjectAccession",
    "properties.sequencingProjectAccession",
    "properties.primaryAssemblyProjectAccession",
    "treatments.reviewing.pairs.completion",
    "treatments.reviewing.pairs.merquryScore",
    "treatments.reviewing.pairs.scoreBuscoEuk",
    "treatments.reviewing.pairs.scoreBuscoTaxon",
    "treatments.reviewing.pairs.taxonBusco",
    "treatments.reviewing.pairs.resultDirectory",
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


def mandatory_pair_value(container: dict, key: str, source: str):
    pair = container.get(key)
    if not isinstance(pair, dict) or pair.get("value") is None:
        raise SubmissionDataError(f"{source} is missing {key}.value")
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

    return {
        "bioproject_umbrella": optional_pair_value(
            properties,
            "umbrellaProjectAccession",
            "NGL-BI properties",
            "bioproject_umbrella",
        ),
        "bioproject_assembly": mandatory_pair_value(
            properties, "primaryAssemblyProjectAccession", "NGL-BI properties"
        ),
        "bioproject_reads": mandatory_pair_value(
            properties, "sequencingProjectAccession", "NGL-BI properties"
        ),
        "sample_code": sample_codes[0],
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


def fetch_main_specimen_code(session: requests.Session, sample_code: str) -> str:
    params = (
        ("codes", sample_code),
        ("includes", "code"),
        ("includes", "properties.individualNumber"),
    )
    response = session.get(
        f"{NGL_SQ_BASE_URL}/api/samples", params=params, timeout=REQUEST_TIMEOUT
    )
    response.raise_for_status()
    samples = response.json()

    if not isinstance(samples, list) or len(samples) != 1:
        count = len(samples) if isinstance(samples, list) else 0
        raise SubmissionDataError(
            f"NGL-SQ query for {sample_code} returned {count} samples"
        )

    returned_code = samples[0].get("code")
    if returned_code != sample_code:
        raise SubmissionDataError(
            f"NGL-SQ returned sample {returned_code!r} instead of {sample_code!r}"
        )

    properties = samples[0].get("properties")
    if not isinstance(properties, dict):
        raise SubmissionDataError(f"NGL-SQ sample {sample_code} is missing properties")
    return mandatory_pair_value(
        properties, "individualNumber", f"NGL-SQ sample {sample_code} properties"
    )


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


def extract_sequencing(raw_platforms: str | None) -> list | None:
    if not raw_platforms:
        warn("manifest is missing PLATFORM; skipping sequencing")
        return None

    platforms = [platform.strip() for platform in raw_platforms.split(",")]
    if not any(platform in LONG_READ_PLATFORMS for platform in platforms):
        warn(f"manifest PLATFORM ({raw_platforms!r}) contains neither PacBio nor ONT")

    sequencing = [
        platform for platform in platforms if platform in SEQUENCING_PLATFORMS
    ]
    if not sequencing:
        warn(
            f"manifest PLATFORM ({raw_platforms!r}) has no known platform; "
            "skipping sequencing"
        )
        return None
    return sequencing


def extract_manifest_values(manifest: dict) -> dict:
    assembly_name = manifest.get("ASSEMBLYNAME")
    if not assembly_name:
        raise SubmissionDataError("manifest is missing required field: ASSEMBLYNAME")

    ena_sample_code = manifest.get("SAMPLE")
    if not ena_sample_code:
        warn("manifest is missing SAMPLE; skipping ena_sample_code")
        ena_sample_code = None

    tolid = assembly_name.split(".", 1)[0]

    return {
        "ena_sample_code": ena_sample_code,
        "assembly_name": assembly_name,
        "ear_report": f"EARs/{tolid}_EAR.pdf",
        "sequencing": extract_sequencing(manifest.get("PLATFORM")),
    }


def build_submission_json(
    ngl: dict,
    main_specimen_code: str,
    manifest: dict,
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
        "ena_sample_code": manifest["ena_sample_code"],
        "assembly_name": manifest["assembly_name"],
        "ear_report": manifest["ear_report"],
        "sequencing": manifest["sequencing"],
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


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate a submission JSON file from NGL-BI and NGL-SQ data."
    )
    parser.add_argument("--project", required=True, help="NGL project code")
    parser.add_argument("--material", required=True, help="NGL material code")
    parser.add_argument(
        "--manifest", required=True, type=Path, help="Path to the ENA manifest file"
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
        manifest = extract_manifest_values(parse_manifest(args.manifest))
        with requests.Session() as session:
            session.headers.update(REQUEST_HEADERS)
            analysis = fetch_analysis(session, args.project, args.material)
            ngl = extract_analysis_values(analysis)
            main_specimen_code = fetch_main_specimen_code(session, ngl["sample_code"])
        busco_path = resolve_busco_path(args.busco, ngl["result_directory"])
        busco = parse_busco_log(busco_path) if busco_path is not None else {}
        gfastats_path = resolve_gfastats_path(
            args.gfastats, ngl["result_directory"], ngl["busco_lin_dataset"]
        )
        gfastats = parse_gfastats(gfastats_path)
        submission = build_submission_json(
            ngl, main_specimen_code, manifest, busco, gfastats
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
