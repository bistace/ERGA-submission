#!/bin/bash

module unload python
module load python/3.10
pip install --prefix . --ignore-installed -r requirements.txt

module unload python
module load python/3.11
pip install --prefix . --ignore-installed -r requirements.txt

mkdir -p bin
cd bin
ln -s ../get_submission_xmls/generate_manifest.py generate_manifest
ln -s ../get_submission_xmls/release_erz.py release_erz
ln -s ../get_submission_xmls/submit_genome.py submit_genome
ln -s ../get_submission_xmls/submit_study.py submit_study
ln -s ../get_submission_xmls/submit_umbrella.py submit_umbrella
ln -s ../get_submission_xmls/submit_studies
ln -s ../get_submission_xmls/templates
