#!/bin/bash

conda env create --yes -f bakir_analyis_environment.yml
conda activate bakir_analysis
pip install ./bakir

conda env create --yes -f bakir/bakir-env.yml
conda activate bakir
pip install ./bakir

conda deactivate

# Immunannot env
conda env create --yes -f HPRC-Immunanot-annotations/immunannot-env.yaml
mkdir Immuannot-data
wget https://zenodo.org/records/8372992/files/Data-2023Oct27.tar.gz -O Immuannot-data/refData-2023Oct27.tar.gz
tar xzvf Immuannot-data/refData-2023Oct27.tar.gz -C Immuannot-data/

# SKIRT env
conda env create --yes -f HPRC-Skirt-annotations/skirt-env.yml

# reconfigure nextflow.config
conda activate bakir_analysis
conda_env_path=$(which python | awk -F'/bakir_analysis' '{print $1}')
sed -i "s|/data/fordmk/miniconda3/envs|$conda_env_path|g" nextflow.config

# get assembly data
nextflow run HPRC-assemblies-annotations/download.nf -profile bakir_analysis

# generate BAKIR annotations
nextflow run HPRC-assemblies-annotations.nf -profile bakir
python trace_file_means.py bakir-trace.txt > bakir-resources-mean.txt

# generate Immunannot annotations
nextflow run HPRC-Immunannot-annotations.nf -profile immunannot

# generate SKIRT annotations
nextflow run HPRC-Skirt-annotations.nf -profile skirt

# Generate notebook with figures and tables
papermill --kernel python3 figures-data-generation.ipynb figures-data-generation.ipynb