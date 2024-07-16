#!/bin/bash

cd datasets/
wget -r -np -nd -R "index.html" https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182308/suppl/GSE182308_RAW.tar
tar -xvf *.tar
rm *.tar
#gzip *.gz
