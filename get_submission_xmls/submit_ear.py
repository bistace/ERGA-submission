#!/usr/bin/env python3
import argparse
import os
import sys

import boto3
import requests
from botocore.exceptions import BotoCoreError, ClientError
from dotenv import load_dotenv


# Target S3 bucket and prefix for EAR reports. Uploading with the "EARs/"
# prefix creates the folder implicitly (S3 has no real directories).
S3_BUCKET = "incoming"
S3_PREFIX = "EARs"


def main():
    parser = argparse.ArgumentParser(
        description="Upload an EAR report PDF to the bytesea S3 bucket."
    )
    parser.add_argument("--project", help="Project code")
    parser.add_argument("--material", help="Material code")
    parser.add_argument(
        "--ear",
        help="Path to an EAR report PDF to upload directly, bypassing the NGL-BI lookup",
    )

    args = parser.parse_args()

    if args.ear:
        ear_path = args.ear
    else:
        if not args.project or not args.material:
            print(
                "ERROR: provide --ear, or both --project and --material",
                file=sys.stderr,
            )
            sys.exit(1)
        assembly_json = get_assembly_json(args.project, args.material)
        if assembly_json is None:
            sys.exit(1)
        ear_path = extract_ear_report(assembly_json)

    if not os.path.isfile(ear_path):
        print(f"ERROR: EAR report file not found at '{ear_path}'", file=sys.stderr)
        sys.exit(1)

    load_dotenv()
    upload_ear(ear_path)


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


def extract_ear_report(assembly: dict) -> str:
    """Extract the EAR report path from the Reviewing treatment in NGL-BI."""
    reviewing = assembly.get("treatments", {}).get("reviewing")
    if not reviewing:
        print("ERROR: reviewing treatment not found in NGL-BI", file=sys.stderr)
        sys.exit(1)

    ear_path = reviewing.get("pairs", {}).get("earReport", {}).get("value")
    if not ear_path:
        print(
            "ERROR: earReport not found in the Reviewing treatment in NGL-BI",
            file=sys.stderr,
        )
        sys.exit(1)

    return ear_path


# --- S3 upload ---


def upload_ear(ear_path: str):
    """Upload the EAR report to s3://incoming/EARs/<basename>."""
    key = f"{S3_PREFIX}/{os.path.basename(ear_path)}"
    client = boto3.client("s3", endpoint_url=os.environ.get("AWS_ENDPOINT_URL"))

    print(f"Uploading {ear_path} to s3://{S3_BUCKET}/{key} ...", file=sys.stderr)
    try:
        client.upload_file(ear_path, S3_BUCKET, key)
    except (BotoCoreError, ClientError) as e:
        print(f"ERROR: failed to upload EAR report to S3: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Uploaded EAR report to s3://{S3_BUCKET}/{key}", file=sys.stderr)


if __name__ == "__main__":
    main()
