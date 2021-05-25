#! /usr/bin/env bash 

set -e

for pdf_path in ../images/four-leaf-*grid.pdf
do
    path_prefix="${pdf_path/\.pdf/}"
    if [[ $path_prefix != *cropped ]]
    then
        cropped_path="${pdf_path/\.pdf/-cropped\.pdf}"
        pdfcrop "$pdf_path" "$cropped_path"
    fi
done
