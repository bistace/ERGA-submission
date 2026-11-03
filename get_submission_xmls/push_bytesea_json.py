#!/usr/bin/env python3
"""Push BytESea assembly metrics JSON files to the Genoby assembly_metrics route.

Each file is sent on its own as {"payload": [<json>], "dry": <bool>}. A file that
the server accepts is moved from pending to the completed directory, unless
--dry-run is set.
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

import requests

ENDPOINT = "https://genome.atlasea.fr/genoby/genoscope/assembly_metrics/"
PENDING_DIRECTORY = Path("/env/cns/proj/projet_DKC/nda_bytesea/pending")
COMPLETED_DIRECTORY = Path("/env/cns/proj/projet_DKC/nda_bytesea/completed")
DEFAULT_TOKEN_FILE = Path.home() / ".bytesea_token" / "token.txt"


def push(session: requests.Session, path: Path, dry: bool) -> bool:
    try:
        document = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        print(f"{path}: cannot read JSON: {error}", file=sys.stderr)
        return False

    response = session.post(ENDPOINT, json={"payload": [document], "dry": dry}, timeout=120)
    try:
        body = response.json()
    except ValueError:
        body = response.text

    metrics = body.get("atlaseadb_assembly_metrics", {}) if isinstance(body, dict) else {}
    accepted = response.ok and metrics.get("response", {}).get("success") is True
    stream = sys.stdout if accepted else sys.stderr
    print(f"{path}: HTTP {response.status_code}", file=stream)
    print(json.dumps(body, indent=2) if isinstance(body, dict) else body, file=stream)
    return accepted


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "files",
        nargs="*",
        type=Path,
        help=f"JSON files to push (default: every *.json in {PENDING_DIRECTORY})",
    )
    parser.add_argument("--token-file", type=Path, default=DEFAULT_TOKEN_FILE)
    parser.add_argument("--completed-dir", type=Path, default=COMPLETED_DIRECTORY)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="ask the server to validate without writing; files are not moved",
    )
    args = parser.parse_args()

    files = args.files or sorted(PENDING_DIRECTORY.glob("*.json"))
    if not files:
        print("No JSON file to push.", file=sys.stderr)
        return 1

    try:
        token = args.token_file.read_text().strip()
    except OSError as error:
        print(f"Cannot read token: {error}", file=sys.stderr)
        return 1

    session = requests.Session()
    session.headers["Authorization"] = f"Bearer {token}"

    failures = 0
    for path in files:
        if not push(session, path, args.dry_run):
            failures += 1
            continue
        if not args.dry_run:
            args.completed_dir.mkdir(parents=True, exist_ok=True)
            destination = args.completed_dir / path.name
            shutil.move(path, destination)
            print(f"Moved {path} -> {destination}")

    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
