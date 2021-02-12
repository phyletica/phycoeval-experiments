#! /usr/bin/env bash 

set -e

for tex_path in ../trees/*.tex
do
    path_prefix="${tex_path/\.tex/}"
    if [[ $path_prefix != *tikz ]]
    then
        pdf_path="${tex_path/\.tex/\.pdf}"
        cropped_path="${pdf_path/\.pdf/-cropped\.pdf}"
        pdfcrop "$pdf_path" "$cropped_path"
    fi
done

for tex_path in ../trees/4-leaf-tree*.tex
do
    path_prefix="${tex_path/\.tex/}"
    if [[ $path_prefix != *tikz ]]
    then
        pdf_path="${tex_path/\.tex/-cropped\.pdf}"
        blank_path="${pdf_path/\.pdf/-blank\.pdf}"
        convert "$pdf_path" -threshold -1 -alpha off "$blank_path"
    fi
done
