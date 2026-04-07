#!/usr/bin/env python3
import argparse
import configparser
import ngs_workflow.env
import ngl.analyses
import os
import requests
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET


os.environ["CONFFILE"] = "/env/atelier/ngs_ba/cns/conf/prod_ba.conf"
ngs_workflow.env.load_conf_file()


def read_credentials(filename=os.path.join(os.environ["HOME"], ".EBI/ebi.ini")):
    config = configparser.ConfigParser()
    config.read(filename)
    account = config.get("Credentials", "account")
    password = config.get("Credentials", "password")
    return account, password


def main():
    parser = argparse.ArgumentParser(
        description="Submit a genome to the ENA and validate its input files."
    )
    parser.add_argument("--project", required=True, help="Project code")
    parser.add_argument("--material", required=True, help="Material code")
    parser.add_argument("--manifest", required=True, help="Path to ENA manifest file")

    args = parser.parse_args()

    manifest = parse_manifest(args.manifest)
    study, sample, assembly_name = extract_manifest_fields(manifest)

    assembly_json = get_assembly_json(args.project, args.material)
    if assembly_json is None:
        print("Failed to get NGS-BA json from NGL.")
        sys.exit(1)
    sample_json = get_sample_json(args.project, args.material)
    if sample_json is None:
        print("Failed to get sample json from NG.")
        sys.exit(1)

    ebi_taxid = get_sample_taxid(sample)
    ngl_study, ngl_tolid, ngl_taxid = extract_ngl_fields(assembly_json, sample_json)
    validate(study, assembly_name, ebi_taxid, ngl_study, ngl_tolid, ngl_taxid)

    print(f"Manifest study: {study}")
    print(f"NGL-BI study:   {ngl_study}")
    print(f"Manifest assembly name: {assembly_name}")
    print(f"NGL-BI ToLID:           {ngl_tolid}")
    print(f"Manifest sample: {sample}")
    print(f"EBI taxid:       {ebi_taxid}")
    print(f"NGL taxid:       {ngl_taxid}")

    cred_path = os.path.join(os.environ["HOME"], ".EBI/ebi.ini")
    if not os.path.exists(cred_path):
        print(f"ERROR: credentials not found at path '{cred_path}'", file=sys.stderr)
        sys.exit(1)

    account, password = read_credentials(cred_path)
    webin_cli_jar = download_webin_cli()
    try:
        # submit_genome(webin_cli_jar, args.manifest, account, password)
        # update_ngl(args.project, args.material, assembly_name)
        pass
    finally:
        if os.path.exists(webin_cli_jar):
            os.remove(webin_cli_jar)


# --- Manifest parsing ---


def parse_manifest(path):
    """Parse an ENA manifest file (tab or space separated) into a dict."""
    fields = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split(None, 1)
            if len(parts) >= 2:
                fields[parts[0]] = parts[1]
    return fields


def extract_manifest_fields(manifest):
    """Extract and validate required fields from a parsed manifest."""
    study = manifest.get("STUDY")
    sample = manifest.get("SAMPLE")
    assembly_name = manifest.get("ASSEMBLYNAME")

    if not study:
        print("ERROR: STUDY field not found in manifest", file=sys.stderr)
        sys.exit(1)
    if not sample:
        print("ERROR: SAMPLE field not found in manifest", file=sys.stderr)
        sys.exit(1)
    if not assembly_name:
        print("ERROR: ASSEMBLYNAME field not found in manifest", file=sys.stderr)
        sys.exit(1)

    return study, sample, assembly_name


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


def get_sample_json(project_code: str, material_code: str) -> dict:
    """Fetch sample data from the NGL-BI API."""
    headers = {"User-Agent": "bot"}
    url = f"http://ngl-bi.genoscope.cns.fr/api/samples/{project_code}_{material_code}"
    res = requests.get(url, headers=headers)
    if res.status_code == 404:
        print(
            f"ERROR: Sample {project_code}_{material_code} not found in NGL-BI!",
            file=sys.stderr,
        )
        return None
    res.raise_for_status()
    return res.json()


def extract_ngl_fields(assembly, sample):
    """Extract study accession, tolid and taxid from an NGL-BI assembly response."""
    properties = assembly.get("properties", {})

    ngl_study = properties.get("primaryAssemblyProjectAccession", {}).get("value")
    ngl_tolid = properties.get("tolid", {}).get("value")
    ngl_taxid = sample.get("taxonCode")

    if not ngl_study:
        print("ERROR: primaryAssemblyProjectAccession not found in NGL-BI", file=sys.stderr)
        sys.exit(1)
    if not ngl_tolid:
        print("ERROR: tolid not found in NGL-BI", file=sys.stderr)
        sys.exit(1)
    if not ngl_taxid:
        print("ERROR: taxonCode not found in NGL-BI", file=sys.stderr)
        sys.exit(1)

    return ngl_study, ngl_tolid, ngl_taxid


def update_ngl(project_code: str, material_code: str, assembly_name: str):
    """Update assembly name and downloaded from NCBI checkbox"""
    code = f"BA.{project_code}_{material_code}"
    ngl.analyses.update_downloaded_from_ncbi(code, False)
    ngl.analyses.update_assembly_to_download_version(code, assembly_name)


# --- EBI API ---


def get_sample_taxid(sample_accession: str) -> str:
    """Fetch the taxid for a given sample accession from the ENA browser API."""
    url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{sample_accession}?download=true&includeLinks=false"
    response = requests.get(url)
    response.raise_for_status()

    root = ET.fromstring(response.content)
    taxon_id_elem = root.find(".//SAMPLE_NAME/TAXON_ID")
    if taxon_id_elem is None or not taxon_id_elem.text:
        print(
            f"ERROR: TAXON_ID not found in EBI XML for sample '{sample_accession}'",
            file=sys.stderr,
        )
        sys.exit(1)

    return taxon_id_elem.text


# --- Validation ---


def validate(study, assembly_name, taxid, ngl_study, ngl_tolid, ngl_taxid):
    """Validate that manifest fields match NGL-BI data."""
    if ngl_study != study:
        print(
            f"ERROR: STUDY mismatch: manifest has '{study}' but NGL-BI has '{ngl_study}'",
            file=sys.stderr,
        )
        sys.exit(1)

    if not assembly_name.startswith(ngl_tolid):
        print(
            f"ERROR: ASSEMBLYNAME '{assembly_name}' does not contain tolid '{ngl_tolid}'",
            file=sys.stderr,
        )
        sys.exit(1)

    if taxid != ngl_taxid:
        print(
            f"ERROR: TAXID mismatch: sample has '{taxid}' but NGL-BI has '{ngl_taxid}'",
            file=sys.stderr,
        )
        sys.exit(1)


# --- Webin CLI ---


def download_webin_cli():
    """Download the latest webin-cli JAR from GitHub to a unique temp path."""
    api_url = "https://api.github.com/repos/enasequence/webin-cli/releases/latest"
    res = requests.get(api_url)
    res.raise_for_status()
    release = res.json()

    jar_asset = None
    for asset in release.get("assets", []):
        if asset["name"].endswith(".jar"):
            jar_asset = asset
            break

    if jar_asset is None:
        print("ERROR: no JAR asset found in the latest webin-cli release", file=sys.stderr)
        sys.exit(1)

    fd, jar_path = tempfile.mkstemp(prefix="webin-cli-", suffix=".jar")
    os.close(fd)

    print(f"Downloading {jar_asset['name']} to {jar_path} ...", file=sys.stderr)
    download = requests.get(jar_asset["browser_download_url"], stream=True)
    download.raise_for_status()
    with open(jar_path, "wb") as f:
        for chunk in download.iter_content(chunk_size=8192):
            f.write(chunk)

    print(f"Downloaded webin-cli to {jar_path}", file=sys.stderr)
    return jar_path


# --- Submission ---


def submit_genome(webin_cli_jar, manifest_path, account, password):
    """Submit a genome to the ENA using webin-cli."""
    cmd = [
        "java",
        "-jar",
        webin_cli_jar,
        "-manifest",
        manifest_path,
        "-context",
        "genome",
        "-userName",
        account,
        "-password",
        password,
        "-submit",
        "-ascp",
    ]
    print(f" => Submitting: {' '.join(cmd)}", file=sys.stderr)

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    print("STDOUT:\n", out.decode("utf-8"), file=sys.stderr)
    print("STDERR:\n", err.decode("utf-8"), file=sys.stderr)

    if p.returncode != 0:
        print(f"ERROR: webin-cli exited with code {p.returncode}", file=sys.stderr)
        sys.exit(1)

    print("Submission was successful", file=sys.stderr)


if __name__ == "__main__":
    main()
