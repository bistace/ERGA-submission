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


def required_pair_value(container: dict, key: str, source: str):
    pair = container.get(key)
    if not isinstance(pair, dict) or pair.get("value") is None:
        raise SubmissionDataError(f"{source} is missing {key}.value")
    return pair["value"]


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
        raise SubmissionDataError("NGL-BI analysis is missing the reviewing treatment")

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

    return {
        "bioproject_umbrella": required_pair_value(
            properties, "umbrellaProjectAccession", "NGL-BI properties"
        ),
        "bioproject_assembly": required_pair_value(
            properties, "primaryAssemblyProjectAccession", "NGL-BI properties"
        ),
        "bioproject_reads": required_pair_value(
            properties, "sequencingProjectAccession", "NGL-BI properties"
        ),
        "sample_code": sample_codes[0],
        "busco_euk": parse_busco_score(
            required_pair_value(pairs, "scoreBuscoEuk", "NGL-BI reviewing treatment")
        ),
        "busco_lin": parse_busco_score(
            required_pair_value(pairs, "scoreBuscoTaxon", "NGL-BI reviewing treatment")
        ),
        "busco_lin_dataset": required_pair_value(
            pairs, "taxonBusco", "NGL-BI reviewing treatment"
        ),
        "merqury_completion": float(
            required_pair_value(pairs, "completion", "NGL-BI reviewing treatment")
        ),
        "merqury_score": float(
            required_pair_value(pairs, "merquryScore", "NGL-BI reviewing treatment")
        ),
        "result_directory": result_directory,
    }


def parse_busco_score(score: str) -> dict:
    if not isinstance(score, str):
        raise SubmissionDataError(f"invalid NGL-BI BUSCO score: {score!r}")

    match = BUSCO_SCORE_PATTERN.fullmatch(score)
    if match is None:
        raise SubmissionDataError(f"invalid NGL-BI BUSCO score: {score!r}")

    return {
        "s": float(match.group("s")),
        "d": float(match.group("d")),
        "m": float(match.group("m")),
        "f": float(match.group("f")),
        "nbgenes": int(match.group("n")),
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
    return required_pair_value(
        properties, "individualNumber", f"NGL-SQ sample {sample_code} properties"
    )


def resolve_busco_path(override: Path | None, result_directory: str | None) -> Path:
    if override is not None:
        if not override.is_file():
            raise SubmissionDataError(f"BUSCO log was not found: {override}")
        return override

    if not result_directory:
        raise SubmissionDataError(
            "NGL-BI reviewing treatment is missing resultDirectory.value"
        )

    busco_directory = Path(result_directory) / BUSCO_RELATIVE_DIRECTORY
    matches = list(busco_directory.glob(BUSCO_FILENAME_PATTERN))
    if not matches:
        raise SubmissionDataError(
            f"no BUSCO log matching {BUSCO_FILENAME_PATTERN} in {busco_directory}"
        )
    if len(matches) > 1:
        raise SubmissionDataError(
            f"multiple BUSCO logs matching {BUSCO_FILENAME_PATTERN} in {busco_directory}"
        )
    return matches[0]


def parse_busco_log(path: Path) -> dict:
    with path.open(encoding="utf-8") as busco_file:
        log = json.load(busco_file)

    versions = log.get("versions") if isinstance(log, dict) else None
    version = versions.get("busco") if isinstance(versions, dict) else None
    if not isinstance(version, str) or not version:
        raise SubmissionDataError(f"BUSCO log is missing versions.busco: {path}")

    parameters = log.get("parameters")
    dataset_version = (
        parameters.get("datasets_version") if isinstance(parameters, dict) else None
    )
    if not isinstance(dataset_version, str) or not dataset_version:
        raise SubmissionDataError(
            f"BUSCO log is missing parameters.datasets_version: {path}"
        )

    return {"version": version, "dataset_version": dataset_version}


def resolve_gfastats_path(
    override: Path | None,
    result_directory: str | None,
    lineage_dataset: str,
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
    with path.open(encoding="utf-8") as gfastats_file:
        for line in gfastats_file:
            label, separator, raw_value = line.strip().partition(":")
            if not separator or label not in GFASTATS_FIELDS:
                continue

            output_field, converter = GFASTATS_FIELDS[label]
            try:
                metrics[output_field] = converter(raw_value.strip())
            except ValueError as error:
                raise SubmissionDataError(
                    f"invalid {label} value in gfastats file {path}: {raw_value.strip()!r}"
                ) from error

    missing_fields = [
        output_field
        for output_field, _ in GFASTATS_FIELDS.values()
        if output_field not in metrics
    ]
    if missing_fields:
        raise SubmissionDataError(
            f"gfastats file is missing required metric(s): {', '.join(missing_fields)}"
        )
    return metrics


def parse_manifest(path: Path) -> dict:
    fields = {}
    with path.open(encoding="utf-8") as manifest_file:
        for line in manifest_file:
            parts = line.strip().split(None, 1)
            if len(parts) == 2:
                fields[parts[0]] = parts[1]
    return fields


def extract_manifest_values(manifest: dict) -> dict:
    required_fields = ("SAMPLE", "ASSEMBLYNAME", "PLATFORM")
    missing_fields = [field for field in required_fields if not manifest.get(field)]
    if missing_fields:
        raise SubmissionDataError(
            f"manifest is missing required field(s): {', '.join(missing_fields)}"
        )

    assembly_name = manifest["ASSEMBLYNAME"]
    tolid = assembly_name.split(".", 1)[0]
    platforms = [platform.strip() for platform in manifest["PLATFORM"].split(",")]
    if not any(platform in {"PacBio", "ONT"} for platform in platforms):
        raise SubmissionDataError(
            "manifest PLATFORM must contain at least one of PacBio or ONT"
        )
    sequencing = [
        platform
        for platform in platforms
        if platform in {"PacBio", "ONT", "Arima", "OmniC"}
    ]

    return {
        "ena_sample_code": manifest["SAMPLE"],
        "assembly_name": assembly_name,
        "ear_report": f"EARs/{tolid}_EAR.pdf",
        "sequencing": sequencing,
    }


def build_submission_json(
    ngl: dict,
    main_specimen_code: str,
    manifest: dict,
    busco: dict,
    gfastats: dict,
) -> dict:
    euk = ngl["busco_euk"]
    lineage = ngl["busco_lin"]

    return {
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
        "busco_euk_genome_s": euk["s"],
        "busco_euk_genome_d": euk["d"],
        "busco_euk_genome_m": euk["m"],
        "busco_euk_genome_f": euk["f"],
        "busco_euk_dataset": "eukaryota",
        "busco_euk_nbgenes": euk["nbgenes"],
        "busco_lin_genome_s": lineage["s"],
        "busco_lin_genome_d": lineage["d"],
        "busco_lin_genome_m": lineage["m"],
        "busco_lin_genome_f": lineage["f"],
        "busco_lin_dataset": ngl["busco_lin_dataset"],
        "busco_lin_nbgenes": lineage["nbgenes"],
        "merqury_completion": ngl["merqury_completion"],
        "merqury_score": ngl["merqury_score"],
        "ena_sample_code": manifest["ena_sample_code"],
        "assembly_name": manifest["assembly_name"],
        "ear_report": manifest["ear_report"],
        "sequencing": manifest["sequencing"],
        "busco_version": busco["version"],
        "busco_dataset_version": busco["dataset_version"],
        "nb_scaffolds": gfastats["nb_scaffolds"],
        "size": gfastats["size"],
        "gc_content": gfastats["gc_content"],
        "n50_scaffolds": gfastats["n50_scaffolds"],
        "n90_scaffolds": gfastats["n90_scaffolds"],
        "l50_scaffolds": gfastats["l50_scaffolds"],
        "l90_scaffolds": gfastats["l90_scaffolds"],
        "nb_contigs": gfastats["nb_contigs"],
        "n50_contigs": gfastats["n50_contigs"],
        "n90_contigs": gfastats["n90_contigs"],
        "l50_contigs": gfastats["l50_contigs"],
        "l90_contigs": gfastats["l90_contigs"],
    }


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
        busco = parse_busco_log(busco_path)
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
