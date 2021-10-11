#! /usr/bin/env bash 

set -e

# Crop all tree PDFs created by P4 + TreeGram
# The PDFs already exist from running P4 script, so no need to compile from tex
for tex_path in ../trees/*.tex
do
    path_prefix="${tex_path/\.tex/}"
    if [[ $path_prefix != *tikz ]]
    then
        pdf_path="${tex_path/\.tex/\.pdf}"
        if [ -e "$pdf_path" ]
        then
            cropped_path="${pdf_path/\.pdf/-cropped\.pdf}"
            pdfcrop "$pdf_path" "$cropped_path"
        fi
    fi
done

for tex_path in ../trees/4-leaf-*.tex
do
    path_prefix="${tex_path/\.tex/}"
    if [[ $path_prefix != *tikz ]]
    then
        pdf_path="${tex_path/\.tex/-cropped\.pdf}"
        blank_path="${pdf_path/\.pdf/-blank\.pdf}"
        convert "$pdf_path" -threshold -1 -alpha off "$blank_path"
    fi
done
