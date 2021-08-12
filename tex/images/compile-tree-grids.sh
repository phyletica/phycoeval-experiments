#! /usr/bin/env bash 

set -e

for tex_path in ./four-leaf-*grid.tex
do
    path_prefix="${tex_path/\.tex/}"
    pdf_path="${path_prefix}.pdf"
    cropped_path="${path_prefix}-cropped.pdf"
    png_path="${path_prefix}-cropped.png"

    latexmk -C "$tex_path"
    latexmk -pdf "$tex_path"
    pdfcrop "$pdf_path" "$cropped_path"
    latexmk -C "$tex_path"

    convert -density 600 "$cropped_path" -flatten -strip -resize 2048 -transparent white "PNG8:${png_path}"
done
