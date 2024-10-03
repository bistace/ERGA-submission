#!/usr/bin/env python3

import requests
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom
import argparse
import configparser
import subprocess
import sys
import os
import re
from pathlib import Path

DEFAULT_CHECKLIST = "ERC000011"

def parse_sample_attributes(sample):
    attributes = {}
    for attr in sample.find('SAMPLE_ATTRIBUTES'):
        tag = attr.find('TAG').text
        value = attr.find('VALUE').text
        units = attr.find('UNITS')
        if units is not None:
            units = units.text
        attributes[tag] = (value, units)
    return attributes

def create_virtual_sample(common_attributes, sample_codes, common_taxon_id, common_scientific_name, alias, center, checklist, mandatory_fields):
    sample_attributes = {"alias": alias}
    if center:
        sample_attributes["center_name"] = center
    sample = ET.Element("SAMPLE", sample_attributes)

    # Adding the virtual sample description in the TITLE field
    title = ET.SubElement(sample, "TITLE")
    title.text = ("This sample is a virtual sample of assembled raw reads from multiple physical "
                  "samples of genome and is composed of physical samples {}.".format(', '.join(sample_codes)))

    sample_name = ET.SubElement(sample, "SAMPLE_NAME")
    taxon_id = ET.SubElement(sample_name, "TAXON_ID")
    taxon_id.text = common_taxon_id
    scientific_name = ET.SubElement(sample_name, "SCIENTIFIC_NAME")
    scientific_name.text = common_scientific_name

    sample_attributes = ET.SubElement(sample, "SAMPLE_ATTRIBUTES")

    # Adding checklist attribute
    if checklist:
        sample_attribute = ET.SubElement(sample_attributes, "SAMPLE_ATTRIBUTE")
        tag_elem = ET.SubElement(sample_attribute, "TAG")
        tag_elem.text = "ENA-CHECKLIST"
        value_elem = ET.SubElement(sample_attribute, "VALUE")
        value_elem.text = checklist

    # Adding common attributes and handling mandatory fields
    for tag, (value, units) in common_attributes.items():
        if tag == "ENA-CHECKLIST":
            continue  # Skip this iteration if the tag is "ENA-CHECKLIST"

        sample_attribute = ET.SubElement(sample_attributes, "SAMPLE_ATTRIBUTE")
        tag_elem = ET.SubElement(sample_attribute, "TAG")
        tag_elem.text = tag
        value_elem = ET.SubElement(sample_attribute, "VALUE")
        value_elem.text = value
        if units:
            units_elem = ET.SubElement(sample_attribute, "UNITS")
            units_elem.text = units

    # Add mandatory fields with default value if they are missing
    for mandatory_field in mandatory_fields:
        if mandatory_field not in common_attributes:
            sample_attribute = ET.SubElement(sample_attributes, "SAMPLE_ATTRIBUTE")
            tag_elem = ET.SubElement(sample_attribute, "TAG")
            tag_elem.text = mandatory_field
            value_elem = ET.SubElement(sample_attribute, "VALUE")
            value_elem.text = "missing: synthetic construct"

    return sample

def create_submission_xml():
    submission = ET.Element("SUBMISSION")
    actions = ET.SubElement(submission, "ACTIONS")
    action = ET.SubElement(actions, "ACTION")
    add = ET.SubElement(action, "ADD")
    return submission

def create_release_xml(sample_code):
    submission = ET.Element("SUBMISSION")
    actions = ET.SubElement(submission, "ACTIONS")
    modify_action = ET.SubElement(actions, "ACTION")
    ET.SubElement(modify_action, "MODIFY")
    release_action = ET.SubElement(actions, "ACTION")
    release = ET.SubElement(release_action, "RELEASE", target=sample_code)
    return submission

def download_sample_xml(output_dir, sample_code):
    print(f"Downloading XML for sample code {sample_code}")
    url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{sample_code}?download=true&includeLinks=false"
    response = requests.get(url)
    if response.status_code == 200:
        print(f"Successfully downloaded XML for sample code {sample_code}")
        # Write response to a file
        sample_xml_path = os.path.join(output_dir, f"{sample_code}.xml")
        with open(sample_xml_path, 'w') as f:
            f.write(response.text)
        return ET.ElementTree(ET.fromstring(response.content))
    else:
        raise Exception(f"Failed to download XML for sample code {sample_code}")

def prettify_xml(elem):
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

# Check if about the EBI account
def read_credentials():
    home = Path.home()
    credentials_file = home / '.EBI' / 'ebi.ini'

    if not credentials_file.exists():
        print(f"Error: Credentials file {credentials_file} not found.", file=sys.stderr)
        exit(1)

    config = configparser.ConfigParser()
    config.read(credentials_file)
    account = config.get('Credentials', 'account')
    password = config.get('Credentials', 'password')
    return account, password

def submit_to_ebi(account, password, submission, samples, submit, output_dir):
    url = "https://www.ebi.ac.uk/ena/submit/drop-box/submit/" if submit else "https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/"
    curl_command = f'curl -u "{account}:{password}" -F "SUBMISSION=@{submission}"'
    if samples:
        curl_command += f' -F "SAMPLE=@{samples}"'
    curl_command += f' {url}'
    print(f" => Executing command: {curl_command}", file=sys.stderr)
    try:
        result = subprocess.run(curl_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("=============== EBI's response =================\n", result.stdout.decode('utf-8'), "\n==============================================")

        # Write response to a file
        response_output_path = os.path.join(output_dir, "submission_response.xml")
        with open(response_output_path, 'w') as response_file:
            response_file.write(result.stdout.decode('utf-8'))

        return result.stdout.decode('utf-8')

    except subprocess.CalledProcessError as e:
        print(f"Error running curl command: {e}", file=sys.stderr)
        print("STDERR:", e.stderr.decode('utf-8'), file=sys.stderr)
        return None

def extract_sample_code_from_response(response):
    print("Extracting sample code from response")
    # Try to extract from success response
    match = re.search(r'<SAMPLE accession="(ERS\d+)"', response)
    if match:
        print(f"Found sample code: {match.group(1)}")
        return match.group(1)

    # Try to extract from error message
    match = re.search(r'accession: "(ERS\d+)"', response)
    if match:
        print(f"Found sample code: {match.group(1)}")
        return match.group(1)

    print("No sample code found in response")
    return None

def mandatory_fields_and_units(checklist_accession):
    mandatory_fields = []
    fields_with_units = {}
    recommended_fields = []
    # URL to download the sample checklist XML from ENA
    url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{checklist_accession}"

    try:
        # Download the XML file
        response = requests.get(url)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Unknown ENA checklist {checklist_accession} => Error downloading XML file: {e}", file=sys.stderr)
        exit(1)

    # Parse the XML content
    root = ET.fromstring(response.content)

    # Extract name of the checklist
    desc = root.find(".//DESCRIPTOR")
    name = desc.find("NAME").text

    # Extract mandatory fields
    for field in root.findall(".//FIELD"):
        field_name = field.find("NAME").text

        units_element = field.find("UNITS")
        if units_element is not None:
            unit_element = units_element.find("UNIT")
            unit_value = unit_element.text if unit_element is not None else None
            fields_with_units[field_name] = unit_value

        is_mandatory = field.find("MANDATORY").text == 'mandatory'
        if is_mandatory:
            mandatory_fields.append(field_name)
        is_recommended = field.find("MANDATORY").text == 'recommended'
        if is_recommended:
            recommended_fields.append(field_name)

    return name, mandatory_fields, fields_with_units, recommended_fields

def main(sample_codes, output_dir, submit, release, checklist, alias, center, force):
    print("Starting main function")

    # Set up the authentication for the HTTP request
    account, password = read_credentials()

    if not release and os.path.exists(output_dir):
        print(f"Error: Output directory '{output_dir}' already exists. Please specify a new directory.", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if release:
        print("Processing release")
        # Read the sample code from the previous submission response
        response_output_path = os.path.join(output_dir, "submission_response.xml")
        with open(response_output_path, 'r') as response_file:
            response = response_file.read()
            sample_code = extract_sample_code_from_response(response)

        if not sample_code:
            raise ValueError("No sample code found in submission response.")

        # Create release XML
        release_set = create_release_xml(sample_code)

        # Write the release XML to the output directory
        release_output_path = os.path.join(output_dir, "release.xml")
        with open(release_output_path, 'w') as f:
            f.write(prettify_xml(release_set))

        # Submit release XML to EBI
        virtual_sample_xml_path = os.path.join(output_dir, "virtual_sample.xml")
        submit_to_ebi(account, password, release_output_path, virtual_sample_xml_path, submit, output_dir)

    else:
        print("Processing submission")
        all_attributes = []
        sample_aliases = []
        common_taxon_id = None
        common_scientific_name = None
        checklist_set = set()
        warnings = []

        for sample_code in sample_codes:
            tree = download_sample_xml(output_dir, sample_code)
            root = tree.getroot()
            sample = root.find('SAMPLE')

            # Extract sample alias
            sample_aliases.append(sample.attrib['alias'])

            # Extract and check taxon_id and scientific_name
            sample_name = sample.find('SAMPLE_NAME')
            taxon_id = sample_name.find('TAXON_ID').text
            scientific_name = sample_name.find('SCIENTIFIC_NAME').text

            if common_taxon_id is None:
                common_taxon_id = taxon_id
            elif common_taxon_id != taxon_id:
                raise ValueError(f"Taxon ID mismatch for sample {sample_code}: {taxon_id} vs {common_taxon_id}")

            if common_scientific_name is None:
                common_scientific_name = scientific_name
            elif common_scientific_name != scientific_name:
                raise ValueError(f"Scientific name mismatch for sample {sample_code}: {scientific_name} vs {common_scientific_name}")

            # Parse sample attributes
            attributes = parse_sample_attributes(sample)
            all_attributes.append(attributes)

            # Collect checklists
            checklist_value = attributes.get("ENA-CHECKLIST", (None, None))[0]
            if checklist_value:
                checklist_set.add(checklist_value)

        # Check for checklist inconsistencies
        if checklist_set and len(checklist_set) > 1 and not checklist:
            raise ValueError("Inconsistent checklists found in input samples and no checklist provided.")

        # Find common attributes
        common_attributes = {}
        if all_attributes:
            first_attributes = all_attributes[0]
            for tag in first_attributes:
                if all(tag in attributes and first_attributes[tag] == attributes[tag] for attributes in all_attributes):
                    common_attributes[tag] = first_attributes[tag]

        # Use provided checklist or the common checklist from attributes
        if checklist:
            final_checklist = checklist
        elif len(checklist_set) == 1:
            final_checklist = checklist_set.pop()
        else:
            final_checklist = DEFAULT_CHECKLIST
            warnings.append("Using default checklist due to inconsistencies or missing checklist in input samples.")

        # Fetch mandatory fields for the checklist
        _, mandatory_fields, fields_with_units, recommended_fields = mandatory_fields_and_units(final_checklist)

        # Create a unique alias for the virtual sample if not provided
        virtual_sample_alias = alias if alias else "virtual_sample_" + "_".join(sample_codes)

        # Display warnings and switch to test mode if not forced
        if warnings and submit and not force:
            print("\n".join(warnings))
            print("Warnings detected, switching to test mode.")
            submit = False

        # Create virtual sample
        virtual_sample = create_virtual_sample(common_attributes, sample_codes, common_taxon_id, common_scientific_name, virtual_sample_alias, center, final_checklist, mandatory_fields)

        # Create the final XML structure for the samples
        sample_set = ET.Element("SAMPLE_SET")
        sample_set.append(virtual_sample)

        # Write the sample XML to the output directory
        samples_output_path = os.path.join(output_dir, "virtual_sample.xml")
        with open(samples_output_path, 'w') as f:
            f.write(prettify_xml(sample_set))

        # Create the submission XML
        submission_set = create_submission_xml()

        # Write the submission XML to the output directory
        submission_output_path = os.path.join(output_dir, "submission.xml")
        with open(submission_output_path, 'w') as f:
            f.write(prettify_xml(submission_set))

        # Submit to EBI
        response = submit_to_ebi(account, password, submission_output_path, samples_output_path, submit, output_dir)

        if response:
            sample_code = extract_sample_code_from_response(response)
            if sample_code:
                print(f"Sample submitted with code: {sample_code}")
            else:
                print("Sample already exists.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a virtual sample XML from multiple sample codes and submit to EBI.")
    parser.add_argument("sample_codes", nargs='*', help="List of sample codes.")
    parser.add_argument("--out", dest="output_dir", required=True, help="Output directory for XML files.")
    parser.add_argument("--submit", action='store_true', help="Submit to the EBI production server instead of the test server.")
    parser.add_argument("--release", action='store_true', help="Release a previously submitted virtual sample, do not provide again the sample codes but only the output directory.")
    parser.add_argument("--checklist", help="ENA checklist to use for the virtual sample. If the --checklist is not provided it will use the checklist found in the initial samples, if common. Otherwise ERC000011 (the default checklist) will be used. If --checklist is set, it will override any other checklists")
    parser.add_argument("--alias", help="Alias for the virtual sample.")
    parser.add_argument("--center", help="Center name for the virtual sample.")
    parser.add_argument("--force", action='store_true', help="Force submission despite warnings.")

    args = parser.parse_args()

    if args.release and not args.sample_codes:
        main([], args.output_dir, args.submit, args.release, args.checklist, args.alias, args.center, args.force)
    elif args.sample_codes:
        main(args.sample_codes, args.output_dir, args.submit, args.release, args.checklist, args.alias, args.center, args.force)
    else:
        parser.error("Either sample codes must be provided or --release must be specified.")
