#!/usr/bin/env python
import os
from os import path
import json
import argparse
import sys
from xml.dom import minidom
from collections import OrderedDict
import jinja2
from datetime import datetime
import re
import pandas as pd

# Author: Jessica Gomez-Garrido, CNAG.
# Contact email: jessica.gomez@cnag.eu
# Date:20230602

script_loc = os.path.dirname(sys.argv[0])
env = jinja2.Environment(loader=jinja2.FileSystemLoader(script_loc + "/templates/"))


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
    use,
    locus_tag,
):
    study_title = species + " " + study_type
    description = ""
    study_name = tolid
    if study_type == "genome assembly":
        alias = cname.replace(" ", "_") + "_genome_assembly"
        study_title = env.get_template("assembly_title.txt").render(
            species=species, cname=cname, tolid=tolid
        )
        if project == "ERGA-pilot":
            description_template = "pilot_assembly_description.txt"
        elif project == "CBP":
            description_template = "cbp_assembly_description.txt"
        elif project == "ERGA-BGE":
            alias = (
                "erga-bge-" + tolid + "_primary-" + datetime.now().strftime("%Y-%m-%d")
            )
            study_title = env.get_template("bge_assembly_title.txt").render(
                species=species, tolid=tolid
            )
            description_template = "bge_assembly_description.txt"
        else:
            description_template = "other_assembly_description.txt"
        description = env.get_template(description_template).render(
            species=species,
            cname=cname,
            sample_coordinator=sample_coordinator,
            use=use.lower(),
        )
    elif study_type == "alternate assembly":
        alt_use = use
        if alternate_annot == "no":
            alt_use = "assembly"
        alias = cname.replace(" ", "_") + "_alternate_genome_assembly"
        study_name += ", alternate haplotype"
        study_title = env.get_template("alternate_assembly_title.txt").render(
            species=species, cname=cname, tolid=tolid
        )
        if project == "ERGA-pilot":
            description_template = "pilot_alternate_assembly_description.txt"
        elif project == "CBP":
            description_template = "cbp_alternate_assembly_description.txt"
        elif project == "ERGA-BGE":
            alias = (
                "erga-bge-"
                + tolid
                + "_alternate-"
                + datetime.now().strftime("%Y-%m-%d")
            )
            study_title = env.get_template("bge_alternate_assembly_title.txt").render(
                species=species, tolid=tolid
            )
            description_template = "bge_alternate_assembly_description.txt"
        else:
            description_template = "other_alternate_assembly_description.txt"
        description = env.get_template(description_template).render(
            species=species,
            cname=cname,
            sample_coordinator=sample_coordinator,
            use=alt_use.lower(),
        )
    elif study_type == "resequencing Data":
        alias = cname.replace(" ", "_") + "_resequencing_data"
        study_register[tolid_pref] = alias
        description = (
            "This project collects the "
            + study_type
            + " generated for "
            + species
            + " (common name "
            + cname
            + ")"
        )
    else:
        alias = cname.replace(" ", "_") + "_data"
        study_title = env.get_template("data_title.txt").render(
            species=species, cname=cname, data=study_type
        )
        study_name = tolid_pref
        if project == "ERGA-pilot":
            description_template = "pilot_data_description.txt"
        elif project == "CBP":
            description_template = "cbp_data_description.txt"
        elif project == "ERGA-BGE":
            alias = (
                "erga-bge-"
                + tolid_pref
                + "-study-rawdata-"
                + datetime.now().strftime("%Y-%m-%d")
            )
            study_title = env.get_template("bge_data_title.txt").render(
                species=species, data=study_type
            )
            description_template = "bge_data_description.txt"
        else:
            description_template = "other_data_description.txt"
        study_register[tolid_pref] = alias
        description = env.get_template(description_template).render(
            species=species,
            cname=cname,
            sample_coordinator=sample_coordinator,
            data=study_type,
            use=use.lower(),
        )

    if "all" in args.xml or "study" in args.xml:
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
    if study_type == "genome assembly" and locus_tag != "-":
        loc = ""
        loc = get_attributes(
            root["study"], seqp, loc, "LOCUS_TAG_PREFIX", **{locus_tag: ""}
        )

    if study_type == "alternate assembly" and alternate_annot == "yes":
        hlocus_tag = locus_tag + "H2"
        loc = ""
        loc = get_attributes(
            root["study"], seqp, loc, "LOCUS_TAG_PREFIX", **{hlocus_tag: ""}
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--files", required=True, help="TABLE file with appropriate headers"
    )
    parser.add_argument(
        "-p",
        "--project",
        default="ERGA-BGE",
        choices=["ERGA-BGE", "CBP", "ERGA-pilot", "EASI", "ATLASEA", "other"],
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
    parser.add_argument("-l", "--locus_tag", required=False, help="Locus tag to register for the study")
    args = parser.parse_args()

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

    get_studies(
        args.project,
        center,
        cname,
        tolid,
        species,
        sample_coordinator,
        type,
        locus_tag,
    )

    for i in root:
        xml_str = root[i].toprettyxml(indent="\t")
        save_path_file = args.out_prefix + "." + i + ".xml"
        with open(save_path_file, "w") as f:
            f.write(xml_str)
