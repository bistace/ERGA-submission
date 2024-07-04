#!/usr/bin/env python
import argparse
import configparser
import jinja2
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

from collections import OrderedDict
from datetime import datetime
from xml.dom import minidom

# Author: Jessica Gomez-Garrido, CNAG.
# Contact email: jessica.gomez@cnag.eu
# Date:20230602

script_loc = os.path.dirname(os.path.abspath(sys.argv[0]))
env = jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.join(script_loc, "templates/")))


def get_attributes(root, parent, child, attr, **element):
    child = root.createElement(attr)
    if element:
        for key in element:
            if element[key] != "":
                child.setAttribute(key, element[key])
            else:
                text = root.createTextNode(key)
                child.appendChild(text)
    parent.appendChild(child)
    return child


def get_xml(project, center, species, tolid_pref, description, children):

    projects = root.createElement("PROJECT")
    projects.setAttribute("center_name", center)

    projects.setAttribute("alias", alias)
    xml.appendChild(projects)

    attr = OrderedDict()
    attr["NAME"] = tolid_pref
    attr["TITLE"] = species
    attr["DESCRIPTION"] = description

    for key in attr:
        attributes = root.createElement(key)
        text = root.createTextNode(attr[key])
        attributes.appendChild(text)
        projects.appendChild(attributes)

    attributes = get_attributes(root, projects, attributes, "UMBRELLA_PROJECT")
    organism = root.createElement("ORGANISM")
    attributes.appendChild(organism)
    for key in org_attr:
        org = root.createElement(key)
        text = root.createTextNode(org_attr[key])
        org.appendChild(text)
        organism.appendChild(org)

    if children:
        attributes = get_attributes(root, projects, attributes, "RELATED_PROJECTS")
        for key in children:
            seqp = get_attributes(root, attributes, attributes, "RELATED_PROJECT")
            accessions = get_attributes(
                root, seqp, seqp, "CHILD_PROJECT", **{"accession": key}
            )

    if args.project == "CBP" or project == "EASI" or project == "ERGA-BGE":
        keyword = {}
        keyword["TAG"] = "Keyword"
        keyword["VALUE"] = project
        attributes = ""
        study_attr = ""
        attributes = get_attributes(root, projects, attributes, "PROJECT_ATTRIBUTES")
        study_attr = get_attributes(
            root, attributes, study_attr, "PROJECT_ATTRIBUTE", **keyword
        )


def read_credentials(filename=os.path.join(os.environ["HOME"], ".EBI/ebi.ini")):
    config = configparser.ConfigParser()
    config.read(filename)
    account = config.get('Credentials', 'account')
    password = config.get('Credentials', 'password')
    return account, password


def generate_submission_xml(release=False):
    sub = ET.Element('SUBMISSION')
    actions = ET.SubElement(sub, 'ACTIONS')
    action = ET.SubElement(actions, 'ACTION')
    ET.SubElement(action, 'ADD')
    if release:
        action = ET.SubElement(actions, 'ACTION')
        ET.SubElement(action, "HOLD", {"HoldUntilDate": datetime.now().strftime("%Y-%m-%d")})

    xml_string = minidom.parseString(ET.tostring(sub)).toprettyxml(indent="\t", encoding="utf-8")
    with open("submission.xml", 'wb') as xml_file:
        xml_file.write(xml_string)


def submit_study(xml_path, test=True, release=False):
    account, password = read_credentials()
    generate_submission_xml(release)

    url = ""
    if test:
        url = '"https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/"'
    else:
        url = '"https://www.ebi.ac.uk/ena/submit/drop-box/submit"'

    curl_command = f"curl -u {account}:{password} " \
        f"""-F "SUBMISSION=@submission.xml" -F "PROJECT=@{xml_path}" """\
        f'{url}'
    print(f" => Submitting: {curl_command}", file=sys.stderr)

    p = subprocess.Popen(
        curl_command, shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    out, err = p.communicate()

    receipt = ET.fromstring(out.decode("utf-8"))
    success = receipt.get("success")
    if success == "false":
        print(f"Error running curl command, return code: {p.returncode}", file=sys.stderr)
        print("STDOUT:", out.decode("utf-8"), file=sys.stderr)
        print("STDERR:", err.decode("utf-8"), file=sys.stderr)
        sys.exit(1)

    # Print the project ID in case of success
    for child in receipt:
        if child.tag == "PROJECT":
            print(child.get("accession"))
            break

    if test:
        print("Test submission was successfull", file=sys.stderr)
    else:
        print("Submission was successfull", file=sys.stderr)

    print("STDOUT: \n", out.decode("utf-8"), file=sys.stderr)
    with open(xml_path.replace(".xml", ".receipt.xml"), "w") as fout:
        print("STDOUT: \n", out.decode("utf-8"), file=fout)

    print("STDERR: \n", err.decode("utf-8"), file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-p",
        "--project",
        default="ERGA-BGE",
        choices=["ERGA-BGE", "CBP", "ERGA-pilot", "EASI", "ATLASea", "other"],
        help="project",
    )
    parser.add_argument("-c", "--center", default="Genoscope", help="center name (Default: Genoscope)")
    parser.add_argument("-n", "--name", required=False, help="Species common name")
    parser.add_argument(
        "--sample-ambassador",
        required=False,
        help="Sample ambassador (only for ERGA-pilot projects)",
    )
    parser.add_argument("-t", "--tolid", required=True, help="tolid")
    parser.add_argument(
        "-s", "--species", required=True, help="species scientific name"
    )
    parser.add_argument("-x", "--taxon_id", required=True, help="species taxon_id")
    parser.add_argument(
        "-a",
        "--children_accessions",
        required=True,
        nargs="+",
        help="Children projects accessions separated by spaces (Example: '-a PRJEB1 PRJEB2 PRJEB3')",
    )
    parser.add_argument(
        "--submit", 
        dest="commit", action="store_true", required=False, 
        help="Do an actual submission if the test is successfull")
    parser.add_argument("--release", dest="release", action="store_true", required=False, help="Set the study to public")

    args = parser.parse_args()
    root = minidom.Document()
    xml = root.createElement("PROJECT_SET")
    root.appendChild(xml)

    tolid_pref = args.tolid
    if args.project == "ERGA-pilot":
        alias = args.name
        description_template = "pilot_umbrella_description.txt"
        if not args.sample_ambassador:
            exit("Required sample ambassador details for ERGA-pilot projects")
        else:
            sample_ambassador = args.sample_ambassador
    elif args.project == "CBP":
        description_template = "cbp_umbrella_description.txt"
        alias = (
            "cbp-"
            + tolid_pref
            + "-study-umbrella-"
            + datetime.now().strftime("%Y-%m-%d")
        )
    elif args.project == "ERGA-BGE":
        alias = (
            "erga-bge-"
            + tolid_pref
            + "-study-umbrella-"
            + datetime.now().strftime("%Y-%m-%d")
        )
        description_template = "bge_umbrella_description.txt"
    elif args.project == "ATLASea":
        alias = tolid_pref + "-study-umbrella-" + datetime.now().strftime("%Y-%m-%d")
        description_template = "atlasea_umbrella_description.txt"
    else:
        alias = tolid_pref + "-study-umbrella-" + datetime.now().strftime("%Y-%m-%d")
        description_template = "other_umbrella_description.txt"

    org_attr = {}
    org_attr["TAXON_ID"] = args.taxon_id
    org_attr["SCIENTIFIC_NAME"] = args.species

    description = env.get_template(description_template).render(
        species=args.species, alias=alias, sample_ambassador=args.sample_ambassador
    )

    get_xml(
        args.project,
        args.center,
        args.species,
        tolid_pref,
        description,
        args.children_accessions,
    )
    xml_str = root.toprettyxml(indent="\t")

    save_path_file = args.species.replace(" ", "_") + ".umbrella.xml"
    with open(save_path_file, "w") as f:
        f.write(xml_str)

    submit_study(save_path_file, not args.commit, args.release)
