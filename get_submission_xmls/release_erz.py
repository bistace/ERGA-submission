#!/usr/bin/env python3
import argparse
import configparser
import os
import sys
import tempfile
import subprocess
from pathlib import Path
import xml.etree.ElementTree as ET


def read_credentials(filename: Path | str = Path.home() / ".EBI" / "ebi.ini"):
    config = configparser.ConfigParser()
    if not os.path.exists(filename):
        print(f"ERROR: credentials not found at path '{filename}'", file=sys.stderr)
        sys.exit(1)
    config.read(filename)
    try:
        account = config.get("Credentials", "account")
        password = config.get("Credentials", "password")
    except Exception as e:
        print(f"ERROR: invalid credentials file '{filename}': {e}", file=sys.stderr)
        sys.exit(1)
    return account, password


def build_release_submission_xml(target_accession: str) -> str:
    sub = ET.Element("SUBMISSION")
    actions = ET.SubElement(sub, "ACTIONS")
    action = ET.SubElement(actions, "ACTION")
    ET.SubElement(action, "RELEASE", {"target": target_accession})
    return ET.tostring(sub, encoding="utf-8").decode("utf-8")


def post_submission_xml(xml_content: str, account: str, password: str, test: bool = False) -> tuple[bool, str, str]:
    # Write XML to a temp file because ENA drop-box expects multipart form with @file
    tmpdir = tempfile.mkdtemp(prefix="ena_release_")
    xml_path = os.path.join(tmpdir, "submission.xml")
    with open(xml_path, "w", encoding="utf-8") as f:
        f.write(xml_content)

    url = "https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/" if test else "https://www.ebi.ac.uk/ena/submit/drop-box/submit/"

    # Use curl consistently with other scripts in this repo
    curl_cmd = (
        f'curl -sS -u "{account}:{password}" '
        f'-F "SUBMISSION=@{xml_path}" '
        f"{url}"
    )

    proc = subprocess.Popen(
        curl_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    out, err = proc.communicate()
    stdout = out.decode("utf-8", errors="replace").strip()
    stderr = err.decode("utf-8", errors="replace").strip()

    if proc.returncode != 0:
        return False, stdout, stderr or f"curl exited with code {proc.returncode}"

    # Try to parse the receipt XML
    try:
        receipt = ET.fromstring(stdout)
        success_attr = receipt.get("success")
        success = (success_attr is not None and success_attr.lower() == "true")
        if not success:
            # Collect error messages if present
            errors = []
            for err_node in receipt.findall('.//ERROR'):
                msg = (err_node.text or "").strip()
                if msg:
                    errors.append(msg)
            return False, stdout, " | ".join(errors) if errors else "ENA receipt reported success=false"
        return True, stdout, ""
    except ET.ParseError:
        # Not XML; return raw outputs
        return False, stdout, "Failed to parse ENA receipt (not XML)"


def parse_accessions(arg_value: str) -> list[str]:
    # If it's a file that exists, read one per line
    if os.path.isfile(arg_value):
        accs: list[str] = []
        with open(arg_value, "r", encoding="utf-8") as f:
            for line in f:
                tok = line.strip()
                if tok:
                    accs.append(tok)
        return accs

    # If it's a comma-separated list
    if "," in arg_value:
        return [tok.strip() for tok in arg_value.split(",") if tok.strip()]

    # Otherwise, assume a single accession
    return [arg_value.strip()]


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Release one or more ERZ accessions at ENA by submitting a RELEASE action.\n"
            "Input can be: a single ERZ accession, a comma-separated list, or a file path with one ERZ per line."
        )
    )
    parser.add_argument(
        "input",
        help=(
            "ERZ accession, comma-separated list of ERZ accessions, or a file path containing ERZ accessions (one per line)"
        ),
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Use ENA dev endpoint (no actual release). Default: use production",
    )
    args = parser.parse_args()

    accessions = parse_accessions(args.input)
    if not accessions:
        print("No ERZ accessions provided.")
        sys.exit(1)

    # Simple validation/normalization
    normalized = []
    for a in accessions:
        a = a.strip()
        if not a:
            continue
        normalized.append(a)
    accessions = normalized

    account, password = read_credentials()

    overall_success = True
    for acc in accessions:
        xml = build_release_submission_xml(acc)
        ok, stdout, err_msg = post_submission_xml(xml, account, password, test=args.test)
        if ok:
            print(f"{acc}: success")
        else:
            overall_success = False
            if err_msg:
                print(f"{acc}: failure: {err_msg}")
            else:
                print(f"{acc}: failure")

    sys.exit(0 if overall_success else 2)


if __name__ == "__main__":
    main()

