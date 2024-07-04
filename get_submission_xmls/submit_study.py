#!/usr/bin/env python
import argparse
import configparser
import jinja2
import os
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
env = jinja2.Environment(
    loader=jinja2.FileSystemLoader(os.path.join(script_loc, "templates/"))
)


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


def get_studies(
    project,
    center,
    cname,
    tolid,
    species,
    sample_coordinator,
    study_type,
    locus_tag,
):
    study_title = ""
    if study_type == "assembly":
        study_title = "Genome assembly of " + species
    elif study_type == "sequencing":
        study_title = "Sequencing data of " + species

    description = ""
    study_name = tolid
    if study_type == "assembly":
        alias = cname.replace(" ", "_") + "_genome_assembly"
        study_title = env.get_template("assembly_title.txt").render(
            species=species, cname=cname, tolid=tolid
        )
        if project == "ERGA-pilot":
            description_template = "pilot_assembly_description.txt"
        elif project == "ERGA-BGE":
            alias = (
                "erga-bge-" + tolid + "_primary-" + datetime.now().strftime("%Y-%m-%d")
            )
            study_title = env.get_template("bge_assembly_title.txt").render(
                species=species, tolid=tolid
            )
            description_template = "bge_assembly_description.txt"
        elif project == "ATLASea":
            description_template = "atlasea_assembly_description.txt"
            alias = (
                "atlasea-" + tolid + "_primary-" + datetime.now().strftime("%Y-%m-%d")
            )
            study_title = env.get_template("bge_assembly_title.txt").render(
                species=species, tolid=tolid
            )
        description = env.get_template(description_template).render(
            species=species,
            cname=cname,
            sample_coordinator=sample_coordinator,
        )
    elif study_type == "sequencing":
        alias = cname.replace(" ", "_") + "_sequencing_data"
        study_register[tolid] = alias
        description = (
            "This project collects the sequencing data generated for "
            + species
            + " (common name "
            + cname
            + ")"
        )
        if project == "ERGA-pilot":
            description_template = "pilot_data_description.txt"
        elif project == "ERGA-BGE":
            alias = (
                "erga-bge-"
                + tolid
                + "-study-rawdata-"
                + datetime.now().strftime("%Y-%m-%d")
            )
            study_title = env.get_template("bge_data_title.txt").render(
                species=species, data=study_type
            )
            description_template = "bge_data_description.txt"
        elif project == "ATLASea":
            description_template = "atlasea_data_description.txt"
            study_title = env.get_template("bge_data_title.txt").render(species=species)
            alias = ("atlasea-" + tolid + "-study-rawdata-" + datetime.now().strftime("%Y-%m-%d"))
        else:
            description_template = "other_data_description.txt"
        study_register[tolid] = alias
        description = env.get_template(description_template).render(
            species=species,
            cname=cname,
            sample_coordinator=sample_coordinator,
            data=study_type,
        )

    get_study_xml(
        project,
        center,
        alias,
        study_name,
        study_title,
        description,
        study_type,
        locus_tag,
    )


def get_study_xml(
    project, center, alias, study_name, study_title, description, study_type, locus_tag
):

    projects = ""
    elements = {}
    elements["center_name"] = center
    elements["alias"] = alias
    projects = get_attributes(root["study"], study_xml, projects, "PROJECT", **elements)

    attr = OrderedDict()
    attr["NAME"] = study_name
    attr["TITLE"] = study_title
    attr["DESCRIPTION"] = description

    for key in attr:
        attributes = root["study"].createElement(key)
        text = root["study"].createTextNode(attr[key])
        attributes.appendChild(text)
        projects.appendChild(attributes)

    attributes = ""
    attributes = get_attributes(
        root["study"], projects, attributes, "SUBMISSION_PROJECT"
    )

    seqp = root["study"].createElement("SEQUENCING_PROJECT")
    if study_type == "assembly" and locus_tag != "-":
        loc = ""
        loc = get_attributes(
            root["study"], seqp, loc, "LOCUS_TAG_PREFIX", **{locus_tag: ""}
        )
    attributes.appendChild(seqp)

    if project == "CBP" or project == "EASI" or project == "ERGA-BGE":
        keyword = {}
        keyword["TAG"] = "Keyword"
        keyword["VALUE"] = project
        attributes = ""
        study_attr = ""
        attributes = get_attributes(
            root["study"], projects, attributes, "PROJECT_ATTRIBUTES"
        )
        study_attr = get_attributes(
            root["study"], attributes, study_attr, "PROJECT_ATTRIBUTE", **keyword
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
    parser.add_argument(
        "-c", "--center", default="Genoscope", help="center name (Default: Genoscope)"
    )
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
    parser.add_argument("-x", "--taxon-id", required=True, help="species taxon_id")
    parser.add_argument(
        "-l", "--locus-tag", required=False, default="-", 
        help="Locus tag to register for the study (Default: '-' [No locus tag])"
    )
    parser.add_argument(
        "--study-type",
        required=True,
        choices=["assembly", "sequencing"],
        help="Study type",
    )
    parser.add_argument("--submit", dest="commit", action="store_true", required=False, help="Do an actual submission if the test is successfull")
    parser.add_argument("--release", dest="release", action="store_true", required=False, help="Set the study to public")
    args = parser.parse_args()

    cred_path = os.path.join(os.environ["HOME"], ".EBI/ebi.ini")
    if not os.path.exists(cred_path):
        print(f"ERROR: credentials not found at path '{cred_path}'")

    cred_path = os.path.join(os.environ["HOME"], ".EBI/ebi.ini")
    if not os.path.exists(cred_path):
        print(f"ERROR: credentials not found at path '{cred_path}'")

    root = {}
    root["study"] = minidom.Document()
    study_xml = root["study"].createElement("PROJECT_SET")
    root["study"].appendChild(study_xml)
    study_register = {}

    center = args.center
    species = args.species
    tolid = args.tolid
    cname = args.name
    sample_coordinator = args.sample_ambassador
    locus_tag = args.locus_tag
    study_type = args.study_type

    get_studies(
        args.project,
        center,
        cname,
        tolid,
        species,
        sample_coordinator,
        study_type,
        locus_tag,
    )

    save_path = ""
    for i in root:
        xml_str = root[i].toprettyxml(indent="\t")
        save_path_file = species.replace(" ", "_") + "." + i + "." + study_type + ".xml"
        save_path = save_path_file
        with open(save_path_file, "w") as f:
            f.write(xml_str)

    submit_study(save_path, not args.commit, args.release)
