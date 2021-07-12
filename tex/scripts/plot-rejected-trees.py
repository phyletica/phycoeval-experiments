#! /usr/bin/env p4

import os
import sys
import glob
import copy

from PIL import Image

from gram import treegram
from gram import gram 

var.punctuation = var.phylip_punctuation
# var.nexus_getAllCommandComments = True
# var.nexus_readBeastTreeCommandComments=True


def plot_tree_grid(rejected_trees,
        retained_trees,
        num_trees = 4,
        scale = 1.0,
        yscale = 1.0,
        x_space = 1.0,
        y_space = 5.5,
        label_y = 4.5,
        label_size = "Large",
        base_name = "rejected-tree-grid",
        dir_name = "gram"):
    reject_label_idx = 97
    retain_label_idx = reject_label_idx + num_trees
    current_x = 0.0
    max_width = max(
            rejected_trees[0].root_height * scale,
            retained_trees[0].root_height * scale)
    rejected_buffer = (max_width - (rejected_trees[0].root_height * scale)) / 2.0
    retained_buffer = (max_width - (retained_trees[0].root_height * scale)) / 2.0

    g_main = gram.Gram()
    g_main.baseName = base_name
    g_main.dirName = dir_name
    label_text = g_main.text(chr(reject_label_idx).upper(), 0.0, 0.0)
    reject_label_idx += 1
    label_text.textSize = "Large"
    g_main.gX = current_x

    col_label = gram.Gram()
    col_text = col_label.text("Rejected trees", 0.0, 0.0)
    col_text.rotate = 90
    col_text.textSize = "Large"
    col_text.anchor = "center"
    col_label.gX = -1.0
    col_label.gY = -(y_space / 2.2)
    g_main.grams.append(col_label)

    col_label = gram.Gram()
    col_text = col_label.text("Retained trees", 0.0, 0.0)
    col_text.rotate = 90
    col_text.textSize = "Large"
    col_text.anchor = "center"
    col_label.gX = -1.0
    col_label.gY = -y_space - (y_space / 2.2)
    g_main.grams.append(col_label)

    tg = treegram.TreeGram(rejected_trees[0],
            scale = scale,
            yScale = yscale)
    tg.gX = current_x + rejected_buffer
    tg.gY = -label_y
    g_main.grams.append(tg)

    tg = treegram.TreeGram(retained_trees[0],
            scale = scale,
            yScale = yscale)
    tg.setScaleBar(length=0.2, xOffset=0.0, yOffset=-0.6)
    tg.gX = current_x + retained_buffer
    tg.gY = -label_y - y_space
    g_main.grams.append(tg)

    label = gram.Gram()
    label_text = label.text(chr(retain_label_idx).upper(), 0.0, 0.0)
    retain_label_idx += 1
    label_text.textSize = "Large"
    label.gX = current_x
    label.gY = -y_space
    g_main.grams.append(label)

    current_x += (x_space + max_width)
    for i in range(1, num_trees):
        max_width = max(
                rejected_trees[i].root_height * scale,
                retained_trees[i].root_height * scale)
        rejected_buffer = (max_width - (rejected_trees[i].root_height * scale)) / 2.0
        retained_buffer = (max_width - (retained_trees[i].root_height * scale)) / 2.0

        tg = treegram.TreeGram(rejected_trees[i],
                scale = scale,
                yScale = yscale)
        # tg.setScaleBar(length=0.1, xOffset=0.0, yOffset=-0.6)
        tg.gX = current_x + rejected_buffer
        tg.gY = -label_y
        g_main.grams.append(tg)

        label = gram.Gram()
        label_text = label.text(chr(reject_label_idx).upper(), 0.0, 0.0)
        reject_label_idx += 1
        label_text.textSize = "Large"
        label.gX = current_x
        # label.gY = label_y
        g_main.grams.append(label)

        tg = treegram.TreeGram(retained_trees[i],
                scale = scale,
                yScale = yscale)
        # tg.setScaleBar(length=0.1, xOffset=0.0, yOffset=-0.6)
        tg.gX = current_x + retained_buffer
        tg.gY = -label_y - y_space
        g_main.grams.append(tg)

        label = gram.Gram()
        label_text = label.text(chr(retain_label_idx).upper(), 0.0, 0.0)
        retain_label_idx += 1
        label_text.textSize = "Large"
        label.gX = current_x
        label.gY = -y_space
        g_main.grams.append(label)

        max_height = max(rejected_trees[i].root_height, retained_trees[i].root_height)
        max_width = max_height * scale

        current_x += (x_space + max_width)
    g_main.epdf()

def get_root_height(tree):
    node = None
    for l in tree.iterLeavesNoRoot():
        node = l
        break
    height = 0.0
    while node.parent:
        height += node.br.len
        node = node.parent
    return height

def get_trees_by_length(path):
    var.trees = []
    read(path)
    for tree in var.trees:
        root_height = get_root_height(tree)
        tree_len = 0.0
        for node in tree.iterNodesNoRoot():
            tree_len += node.br.len / root_height
        tree.total_length = tree_len
        tree.root_height = root_height
    return sorted(var.trees, key = lambda x: x.total_length, reverse = False)

def main_cli():
    tip_labels = {
            "sp1" : "A",
            "sp2" : "B",
            "sp3" : "C",
            "sp4" : "D",
            "sp5" : "E",
            "sp6" : "F",
            "sp7" : "G",
            "sp8" : "H",
            "sp9" : "I",
            }

    rejected_bif_trees = get_trees_by_length("rejected-bif-trees.nex")
    retained_bif_trees = get_trees_by_length("retained-bif-trees.nex")

    for tree in rejected_bif_trees:
        for leaf in tree.iterLeavesNoRoot():
            leaf.name = tip_labels[leaf.name]
    for tree in retained_bif_trees:
        for leaf in tree.iterLeavesNoRoot():
            leaf.name = tip_labels[leaf.name]

    plot_tree_grid(rejected_bif_trees,
            retained_bif_trees,
            num_trees = 3,
            scale = 10.0,
            yscale = 0.45,
            x_space = 1.0,
            y_space = 4.2,
            label_y = 3.5,
            label_size = "Large",
            base_name = "rejected-bif-trees",
            dir_name = "../trees")


if __name__ ==  '__main__':
    main_cli()
