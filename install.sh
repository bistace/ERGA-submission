#!/bin/bash

module unload python
module load python/3.10
pip install --prefix . --ignore-installed -r requirements.txt

module unload python
module load python/3.11
pip install --prefix . --ignore-installed -r requirements.txt

mkdir -p bin
cd bin
ln -s ../get_submission_xmls/submit_study.py
ln -s ../get_submission_xmls/submit_umbrella.py
ln -s ../get_submission_xmls/submit_studies
ln -s ../get_submission_xmls/templates
