#! /usr/bin/env bash 

set -e

# Plot trees with P4 + TreeGram
./plot-generalized-tree.py

# I couldn't figure out how to change the tree with TreeGram, so below is hack
# to use sed to create coloared trees
highlight_color="myorange"
trees_to_highlight=( \
    "../trees/4-leaf-tree-5-bare.tikz.tex" \
    "../trees/4-leaf-tree-5.tikz.tex" \
    "../trees/4-leaf-labeled-tree-bare-5-0.tikz.tex" \
    "../trees/4-leaf-labeled-tree-bare-5-2.tikz.tex" \
    "../trees/4-leaf-labeled-tree-bare-5-1.tikz.tex" \
)

for tree_path in ${trees_to_highlight[@]}
do
    new_tree_path="${tree_path/\.tikz\.tex/-highlight\.tikz\.tex}"
    template_path="${tree_path/\.tikz\.tex/\.tex}"
    new_template_path="${template_path/\.tex/-highlight\.tex}"
    sed -e "s/\[black\]/\[${highlight_color}\]/g" -e "s/\[black,line/\[${highlight_color},line/g" "$tree_path" > "$new_tree_path"
    sed -e "s|$(basename ${tree_path/\.tex/})|$(basename ${new_tree_path/\.tex/})|g" "$template_path" > "$new_template_path"
    (
        cd "$(dirname $new_template_path)"
        new_template_file="$(basename $new_template_path)"
        latexmk -C "$new_template_file"
        latexmk -pdf "$new_template_file"
        pdf_path="${new_template_file/\.tex/\.pdf}"
        cropped_pdf_path="${new_template_file/\.tex/-cropped\.pdf}"
        pdfcrop "$pdf_path" "$cropped_pdf_path"
        latexmk -C "$new_template_file"
    )
done
