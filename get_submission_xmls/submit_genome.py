#!/usr/bin/env python3
import argparse
import configparser
import ngs_workflow.env
import ngl.analyses
import os
import subprocess
import sys


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
    study, assembly_name = extract_manifest_fields(manifest)

    assembly = get_assembly_json(args.project, args.material)
    if assembly is None:
        sys.exit(1)

    ngl_study, ngl_tolid = extract_ngl_fields(assembly)
    validate(study, assembly_name, ngl_study, ngl_tolid)

    print(f"Manifest study: {study}")
    print(f"NGL-BI study:   {ngl_study}")
    print(f"Manifest assembly name: {assembly_name}")
    print(f"NGL-BI ToLID:           {ngl_tolid}")

    cred_path = os.path.join(os.environ["HOME"], ".EBI/ebi.ini")
    if not os.path.exists(cred_path):
        print(f"ERROR: credentials not found at path '{cred_path}'", file=sys.stderr)
        sys.exit(1)

    account, password = read_credentials(cred_path)
    webin_cli_jar = download_webin_cli()

    # submit_genome(webin_cli_jar, args.manifest, account, password)
    update_ngl(args.project, args.material, assembly_name)


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
    assembly_name = manifest.get("ASSEMBLYNAME")

    if not study:
        print("ERROR: STUDY field not found in manifest", file=sys.stderr)
        sys.exit(1)
    if not assembly_name:
        print("ERROR: ASSEMBLYNAME field not found in manifest", file=sys.stderr)
        sys.exit(1)

    return study, assembly_name


# --- NGL-BI API ---


def get_assembly_json(project_code: str, material_code: str) -> dict:
    """Fetch assembly data from the NGL-BI API."""
    return ngl.get_from_bi(f"analyses/BA.{project_code}_{material_code}")


def extract_ngl_fields(assembly):
    """Extract study accession and tolid from an NGL-BI assembly response."""
    properties = assembly.get("properties", {})

    ngl_study = properties.get("primaryAssemblyProjectAccession", {}).get("value")
    ngl_tolid = properties.get("tolid", {}).get("value")

    if not ngl_study:
        print("ERROR: primaryAssemblyProjectAccession not found in NGL-BI", file=sys.stderr)
        sys.exit(1)
    if not ngl_tolid:
        print("ERROR: tolid not found in NGL-BI", file=sys.stderr)
        sys.exit(1)

    return ngl_study, ngl_tolid


def update_ngl(project_code: str, material_code: str, assembly_name: str):
    """Update assembly name and downloaded from NCBI checkbox"""
    code = f"BA.{project_code}_{material_code}"
    ngl.analyses.update_downloaded_from_ncbi(code, False)
    ngl.analyses.update_assembly_to_download_version(code, assembly_name)


# --- Validation ---


def validate(study, assembly_name, ngl_study, ngl_tolid):
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


# --- Webin CLI ---


def download_webin_cli(dest_dir="/tmp"):
    """Download the latest webin-cli JAR from GitHub into dest_dir."""
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

    jar_path = os.path.join(dest_dir, jar_asset["name"])
    if os.path.exists(jar_path):
        print(f"webin-cli already present at {jar_path}", file=sys.stderr)
        return jar_path

    print(f"Downloading {jar_asset['name']} ...", file=sys.stderr)
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
