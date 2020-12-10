#! /usr/bin/env p4

import os
import sys
import glob
import copy

from PIL import Image
import ete3

from gram import treegram
from gram import gram 

var.punctuation = var.phylip_punctuation

class TreeInfo(object):
    IMAGES = {
            'gv': 'gekko-vittatus-1-green-shadow.png',
            'ov': 'gekko-vittatus-2-orange-shadow.png',
            'bv': 'gekko-vittatus-3-blue-shadow.png',
            'yv': 'gekko-vittatus-4-yellow-shadow.png',
            'tv': 'gekko-vittatus-5-teal-shadow.png',
            'av': 'gekko-vittatus-6-auburn-shadow.png',
            'gg': 'gekko-gecko-1-green-shadow.png',
            'og': 'gekko-gecko-2-orange-shadow.png',
            'bg': 'gekko-gecko-3-blue-shadow.png',
            'yg': 'gekko-gecko-4-yellow-shadow.png',
            'tg': 'gekko-gecko-5-teal-shadow.png',
            'ag': 'gekko-gecko-6-auburn-shadow.png',
            'gp': 'gecko-pixabay-cc0-1-green-shadow.png',
            'op': 'gecko-pixabay-cc0-2-orange-shadow.png',
            'bp': 'gecko-pixabay-cc0-3-blue-shadow.png',
            'yp': 'gecko-pixabay-cc0-4-yellow-shadow.png',
            'tp': 'gecko-pixabay-cc0-5-teal-shadow.png',
            'ap': 'gecko-pixabay-cc0-6-auburn-shadow.png',
            'p': 'gecko-pixabay-cc0.png',
            }
    NODE_LABELS = {
            'T1': '$T_1$',
            'T2': '$T_2$',
            'T3': '$T_3$',
            }
    _TREE = '(((({yv}:{lyv},{tv}:{ltv}):{lytv},{av}:{lav}):{lv},((({yp}:{lyp},{tp}:{ltp}):{lytp},{ap}:{lap}):{lp},(({yg}:{lyg},{tg}:{ltg}):{lytg},{ag}:{lag}):{lg}):{lgp}))root:{lroot};'

    LEAF_STYLE = gram.GramText("x")
    LEAF_STYLE.name = "myleaf"
    LEAF_STYLE.anchor = "west"
    # LEAF_STYLE.textHeight = 0.19
    # LEAF_STYLE.textDepth = 0.19
    LEAF_STYLE.textSize = "normalsize"
    LEAF_STYLE.yShift = -0.2

    def __init__(self, node_heights = None, height_multiplier = 1.0, tree = None):
        self.raw_node_heights = node_heights
        if not node_heights:
            self.raw_node_heights = {
                    'ytv': 0.05,
                    'ytp': 0.08,
                    'ytg': 0.05,
                    'ytav': 0.08,
                    'ytap': 0.08,
                    'ytag': 0.05,
                    'pg': 0.2,
                    'root': 0.2,
                    # 'ytv': 0.05,
                    # 'ytp': 0.11,
                    # 'ytg': 0.03,
                    # 'ytav': 0.06,
                    # 'ytap': 0.14,
                    # 'ytag': 0.12,
                    # 'pg': 0.18,
                    # 'root': 0.2,
                    }
        self.node_heights = {}
        for k, v in self.raw_node_heights.items():
            self.node_heights[k] = v * height_multiplier
        self.height_multiplier = height_multiplier
        if not tree:
            self.tree = self._TREE
        else:
            self.tree = tree

    def _get_stem_lengths(self):
        return {
                'lgp': self.node_heights['root'] - self.node_heights['pg'],
                'lv': self.node_heights['root'] - self.node_heights['ytav'],
                'lg': self.node_heights['pg'] - self.node_heights['ytag'],
                'lp': self.node_heights['pg'] - self.node_heights['ytap'],
                'lytv': self.node_heights['ytav'] - self.node_heights['ytv'],
                'lytg': self.node_heights['ytag'] - self.node_heights['ytg'],
                'lytp': self.node_heights['ytap'] - self.node_heights['ytp'],
                'lav': self.node_heights['ytav'],
                'lag': self.node_heights['ytag'],
                'lap': self.node_heights['ytap'],
                'lyv': self.node_heights['ytv'],
                'lyg': self.node_heights['ytg'],
                'lyp': self.node_heights['ytp'],
                'ltv': self.node_heights['ytv'],
                'ltg': self.node_heights['ytg'],
                'ltp': self.node_heights['ytp'],
                'lroot': self.node_heights['root'] * 0.08,
                }

    stem_lengths = property(_get_stem_lengths)

    def __str__(self):
        d = {}
        d.update(self.IMAGES)
        d.update(self.node_heights)
        d.update(self.stem_lengths)
        return self.tree.format(**d)

    def get_image_path(self, leaf_name):
        return os.path.join(os.path.pardir, 'images', 'phylopics', leaf_name)

    def plot_tree(self,
            base_name = 'tree',
            dir_name = 'gram',
            scale = 1.0,
            yscale = 1.0,
            image_width = 12,
            vmargin = 0,
            lmargin = 15,
            phylo_line_weight = 'ultra thick',
            event_line_weight = 'thick',
            event_line_color = 'black!50',
            event_line_style = 'dashed',
            title_str = None,
            show_title = False,
            title_font_size = 'huge',
            show_events = False,
            event_indices = None,
            star_indices = None,
            star_color = 'red',
            show_event_labels = False,
            remove_event_labels = False,
            event_labels_top = True,
            event_font_size = 'LARGE',
            show_node_labels = False,
            show_tip_images = True,
            event_label_plots = False,
            show_times = False,
            time_font_size = 'large',
            plot_width = 18,
            plot_x_shift = 0.0,
            plot_y_shift = -1.8,
            plot_suffix = "-vln.pdf",
            ):
        read(str(self))
        sys.stderr.write("{0}\n".format(str(self)))
        tree = var.trees[-1]
        nodes_to_shift = []
        for n in tree.iterNodes():
            if n.name in self.NODE_LABELS:
                if show_node_labels:
                    n.name = self.NODE_LABELS[n.name]
                    nodes_to_shift.append(n)
                else:
                    n.name = None
            if n.isLeaf:
                try:
                    img = Image.open(self.get_image_path(n.name))
                except IOError as e:
                    if n.name == 'root':
                        n.name = None
                        continue
                    raise e
                node_name = n.name
                n.name = r"\hspace{{{hspace}mm}}\includegraphics[width={width}mm]{{{path}}}".format(
                        hspace = lmargin,
                        width = image_width,
                        path = self.get_image_path(node_name))
                if not show_tip_images:
                    n.name = " "
        node_v = tree.node(2)
        node_p = tree.node(8)
        node_g = tree.node(13)
        node_root = tree.node(1)
        if event_label_plots:
            node_root.name = r"\includegraphics[width={plot_width}mm]{{../../results/fixed-gen-gen-node-probs-root{plot_suffix}}}".format(
                    plot_width = plot_width,
                    plot_suffix = plot_suffix)
            node_v.name = r"\includegraphics[width={plot_width}mm]{{../../results/fixed-gen-gen-node-probs-12-3{plot_suffix}}}".format(
                    plot_width = plot_width,
                    plot_suffix = plot_suffix)
            node_p.name = r"\includegraphics[width={plot_width}mm]{{../../results/fixed-gen-gen-node-probs-456{plot_suffix}}}".format(
                    plot_width = plot_width,
                    plot_suffix = plot_suffix)
            node_g.name = r"\includegraphics[width={plot_width}mm]{{../../results/fixed-gen-gen-node-probs-789{plot_suffix}}}".format(
                    plot_width = plot_width,
                    plot_suffix = plot_suffix)
        tg = treegram.TreeGram(tree, scale = scale, yScale = yscale) #, showNodeNums = True)
        tg.styleDict[self.LEAF_STYLE.name] = self.LEAF_STYLE
        for p4_node in tree.iterLeavesNoRoot():
            p4_node.label.myStyle = self.LEAF_STYLE.name
        for n in nodes_to_shift:
            # tg.tree.node(n.nodeNum).label.xShift = -0.1
            # tg.tree.node(n.nodeNum).label.yShift = 0.05
            tg.tree.node(n.nodeNum).label.myStyle = 'node right'
        tg.latexUsePackages.append('sfmath')
        tg.latexUsePackages.append('color')
        tg.latexOtherPreambleCommands.extend([
                "\definecolor{mygreen}{RGB}{50,162,81}",
                "\definecolor{myorange}{RGB}{255,127,15}",
                "\definecolor{myblue}{RGB}{60,183,204}",
                "\definecolor{myyellow}{RGB}{255,217,74}",
                "\definecolor{myteal}{RGB}{57,115,124}",
                "\definecolor{myauburn}{RGB}{184,90,13}",
                ])
        # tg.latexOtherPreambleCommands()
        tg.internalNodeLabelSize = 'small' # 'tiny' is default
        tg.tgDefaultLineThickness = phylo_line_weight
        tg.baseName = base_name
        tg.dirName = dir_name
        # tg.grid(0, 0, 8, 5)
        heights = sorted(set(self.node_heights.values()))
        show_indices = [i for i, h in enumerate(heights)]
        if event_indices:
            show_indices = [i for i in event_indices]
        if not star_indices:
            star_indices = []
        lines = []
        line_labels = []
        top = 8.41*yscale
        bottom = -0.3
        horizontal_center = scale*((self.node_heights['root'] + 0.1)/2.0)
        if show_title:
            if title_str:
                title = tg.text(title_str, horizontal_center, top)
            else:
                title = tg.text('$n_{{\\tau}} = {0}$'.format(len(heights)), horizontal_center, top)
            title.anchor = "south"
            title.textSize = title_font_size
        if show_times:
            root_x = scale * 0.1
            root_y = top / 2.0
            t = tg.text("\\sffamily {0}".format(self.raw_node_heights['root']), root_x, root_y)
            t.textSize = time_font_size
            t.anchor = "west"
            # t.xShift = 0.1
        height_xs = []
        for i, h in enumerate(heights):
            x = scale*((self.node_heights['root'] + 0.1) - h)
            height_xs.append(x)
            if not i in show_indices:
                continue
            l = tg.line(x, bottom, x, top)
            l.lineStyle = event_line_style
            l.lineThickness = event_line_weight
            l.colour = event_line_color
            if (i in star_indices) and star_color:
                l.colour = star_color
            if event_label_plots and (i > 1):
                l.colour = 'black!00'
            if not remove_event_labels:
                event_label = ""
                if star_indices:
                    if i in star_indices:
                        # event_label = '{{\\color{{{}}}$\\mathbf{{*}}$}}'.format(star_color)
                    # else:
                        event_label = ' ' 
                else:
                    event_label = '$\\tau_{{\\scriptscriptstyle {0}}}$'.format(i+1) 
                if event_label_plots:
                    if i == 0:
                        event_label = r"\includegraphics[width={plot_width}mm]{{../../results/fixed-gen-gen-height-probs-12-789{plot_suffix}}}".format(
                                plot_width = plot_width,
                                plot_suffix = plot_suffix)
                    elif i == 1:
                        event_label = r"\includegraphics[width={plot_width}mm]{{../../results/fixed-gen-gen-height-probs-123-456{plot_suffix}}}".format(
                                plot_width = plot_width,
                                plot_suffix = plot_suffix)
                    else:
                        event_label = ""
                if event_labels_top:
                    t = tg.text(event_label, x, top)
                    if show_times:
                        t2 = tg.text("\\sffamily {0}".format(h / self.height_multiplier), x, bottom)
                        t2.textSize = time_font_size
                        t2.anchor = "south east"
                else:
                    t = tg.text(event_label, x, bottom)
                    if show_times:
                        t2 = tg.text("\\sffamily {0}".format(h / self.height_multiplier), x, top)
                        t2.textSize = time_font_size
                        t2.anchor = "north east"
                t.textSize = event_font_size
                if event_labels_top:
                    t.anchor = "south"
                    t.yShift = -0.1
                else:
                    t.anchor = "north"
                    t.yShift = 0.1
                # if event_label_plots and (i == 0):
                #     t.xShift = 0.2
                # if event_label_plots and (i == 1):
                #     t.xShift = -0.3
            line_labels.append(t)
            lines.append(l)
        if not show_events:
            for l in lines:
                l.colour = 'black!00'
        if not show_event_labels:
            for l in line_labels:
                l.colour = 'black!00'
        if event_label_plots:
            node_root.label.anchor = "east"
            node_root.label.xShift = plot_x_shift
            node_root.label.yShift = plot_y_shift
            node_v.label.anchor = "east"
            node_v.label.xShift = plot_x_shift
            node_v.label.yShift = plot_y_shift
            node_p.label.anchor = "east"
            node_p.label.xShift = plot_x_shift
            node_p.label.yShift = plot_y_shift
            node_g.label.anchor = "east"
            node_g.label.xShift = plot_x_shift
            node_g.label.yShift = plot_y_shift
        tg.epdf()


class RevJumpMoveTreeInfo(object):
    NODE_LABELS = {
            'n1': '\\sffamily I',
            'n2': '\\sffamily H',
            'n3': '\\sffamily G',
            'n4': '\\sffamily F',
            'n5': '\\sffamily E',
            'n6': '\\sffamily D',
            'n7': '\\sffamily C',
            'n8': '\\sffamily B',
            'n9': '\\sffamily A',
            'n12': '$t_1$',
            'n34': '$t_2$',
            'n1234': '$t_3$',
            'n56': '$t_4$',
            'n789': '$t_5$',
            'root': '$t_6$',
            }
    _TREE = '(((n1:{l1},n2:{l2})n12:{l12},(n3:{l3},n4:{l4})n34:{l34})n1234:{l1234},((n5:{l5},n6:{l6})n56:{l56},((n7:{l7},n8:{l8})n78:{l78},n9:{l9})n789:{l789})n56789:{l56789})root:{lroot};'

    read("((A:0.1,B:0.1)AB:0.1,C:0.2)root:0.0;")
    _DUMMY_TREE = var.trees[-1]
    _DUMMY_TG = treegram.TreeGram(_DUMMY_TREE)
    _DUMMY_TG.render()
    LEAF_STYLE = copy.deepcopy(_DUMMY_TG.styleDict['leaf'])
    LEAF_STYLE.name = "rotated_leaf"
    LEAF_STYLE.anchor = "north"
    # LEAF_STYLE.textHeight = 0.19
    # LEAF_STYLE.textDepth = 0.19
    LEAF_STYLE.textSize = "normalsize"
    LEAF_STYLE.rotate = 90
    # LEAF_STYLE.yShift = -0.2

    INTERNAL_NODE_STYLE = copy.deepcopy(_DUMMY_TG.styleDict['node right'])
    INTERNAL_NODE_STYLE.name = "rotated_internal"
    INTERNAL_NODE_STYLE.anchor = "north"
    INTERNAL_NODE_STYLE.rotate = 90
    INTERNAL_NODE_STYLE.textSize = "normalsize"

    def __init__(self, node_heights = None, height_multiplier = 1.0, tree = None):
        self.raw_node_heights = node_heights
        if not node_heights:
            self.raw_node_heights = {
                    '12'    : 0.2,
                    '34'    : 0.4,
                    '1234'  : 0.7,
                    '56'    : 0.4,
                    '78'    : 0.2,
                    '789'   : 0.2,
                    '56789' : 1.0,
                    'root'  : 1.0,
                    }
        self.node_heights = {}
        for k, v in self.raw_node_heights.items():
            self.node_heights[k] = v * height_multiplier
        self.height_multiplier = height_multiplier
        if not tree:
            self.tree = self._TREE
        else:
            self.tree = tree

    def _get_stem_lengths(self):
        return {
                'l1234': self.node_heights['root'] - self.node_heights['1234'],
                'l56789': self.node_heights['root'] - self.node_heights['56789'],
                'l12': self.node_heights['1234'] - self.node_heights['12'],
                'l34': self.node_heights['1234'] - self.node_heights['34'],
                'l56': self.node_heights['56789'] - self.node_heights['56'],
                'l789': self.node_heights['56789'] - self.node_heights['789'],
                'l78': self.node_heights['789'] - self.node_heights['78'],
                'l1': self.node_heights['12'],
                'l2': self.node_heights['12'],
                'l3': self.node_heights['34'],
                'l4': self.node_heights['34'],
                'l5': self.node_heights['56'],
                'l6': self.node_heights['56'],
                'l7': self.node_heights['78'],
                'l8': self.node_heights['78'],
                'l9': self.node_heights['789'],
                'lroot': self.node_heights['root'] * 0.08,
                }

    def _get_event_labels(self):
        height_set = set(self.node_heights.values())
        height_set.add(0.0)
        heights = sorted(height_set)
        labels = []
        for i, h in enumerate(heights):
            l = '$\\tau_{{\\scriptscriptstyle {0}}}$'.format(i+1) 
            labels.append(l)
        return labels

    stem_lengths = property(_get_stem_lengths)

    def __str__(self):
        d = {}
        # d.update(self.node_heights)
        d.update(self.stem_lengths)
        return self.tree.format(**d)

    def plot_tree(self,
            base_name = 'tree',
            dir_name = 'gram',
            scale = 1.0,
            yscale = 1.0,
            vmargin = 0,
            lmargin = 15,
            phylo_line_weight = 'ultra thick',
            event_line_weight = 'thick',
            event_line_color = 'black!50',
            event_line_style = 'dashed',
            title_str = None,
            show_title = False,
            title_font_size = 'huge',
            show_events = False,
            splittable_indices = [],
            event_indices_to_labels = None,
            show_event_labels = False,
            remove_event_labels = False,
            event_labels_top = True,
            event_font_size = 'LARGE',
            show_node_labels = False,
            show_tip_labels = True,
            time_font_size = 'large',
            ):
        read(str(self))
        sys.stderr.write("{0}\n".format(str(self)))
        tree = var.trees[-1]
        nodes_to_shift = []
        for n in tree.iterNodes():
            if n.isLeaf:
                n.name = self.NODE_LABELS[n.name]
                if not show_tip_labels:
                    n.name = " "
                continue
            if n.name in self.NODE_LABELS:
                if show_node_labels:
                    n.name = self.NODE_LABELS[n.name]
                    nodes_to_shift.append(n)
                else:
                    n.name = None
            else:
                n.name = None
        node_root = tree.node(self.NODE_LABELS["root"])
        tg = treegram.TreeGram(tree, scale = scale, yScale = yscale) #, showNodeNums = True)
        tg.styleDict[self.LEAF_STYLE.name] = self.LEAF_STYLE
        tg.styleDict[self.INTERNAL_NODE_STYLE.name] = self.INTERNAL_NODE_STYLE
        for p4_node in tree.iterLeavesNoRoot():
            p4_node.label.myStyle = self.LEAF_STYLE.name
        for n in nodes_to_shift:
            tg.tree.node(n.nodeNum).label.myStyle = 'rotated_internal'
        tg.tree.node(node_root.nodeNum).label.xShift = -0.5
        tg.latexUsePackages.append('sfmath')
        tg.latexUsePackages.append('color')
        tg.latexOtherPreambleCommands.extend([
                "\definecolor{mygreen}{RGB}{50,162,81}",
                "\definecolor{myorange}{RGB}{255,127,15}",
                "\definecolor{myblue}{RGB}{60,183,204}",
                "\definecolor{myyellow}{RGB}{255,217,74}",
                "\definecolor{myteal}{RGB}{57,115,124}",
                "\definecolor{myauburn}{RGB}{184,90,13}",
                ])
        # tg.latexOtherPreambleCommands()
        tg.internalNodeLabelSize = 'small' # 'tiny' is default
        tg.tgDefaultLineThickness = phylo_line_weight
        tg.baseName = base_name
        tg.dirName = dir_name
        # tg.grid(0, 0, 8, 5)
        height_set = set(self.node_heights.values())
        height_set.add(0.0)
        heights = sorted(height_set)
        show_indices = [i for i, h in enumerate(heights)]
        event_labels = self._get_event_labels()
        if event_indices_to_labels:
            show_indices = sorted(event_indices_to_labels.keys())
            event_labels = [event_indices_to_labels[i] for i in show_indices]
        lines = []
        line_labels = []
        top = 8.41*yscale
        bottom = -0.3
        horizontal_center = scale*((self.node_heights['root'] + 0.1)/2.0)
        if show_title:
            if title_str:
                title = tg.text(title_str, horizontal_center, top)
            else:
                title = tg.text('$n_{{\\tau}} = {0}$'.format(len(heights)), horizontal_center, top)
            title.anchor = "south"
            title.textSize = title_font_size
        height_xs = []
        for i, h in enumerate(heights):
            x = scale*((self.node_heights['root'] + 0.0) - h)
            height_xs.append(x)
            if not i in show_indices:
                continue
            l = tg.line(x, bottom, x, top)
            l.lineStyle = event_line_style
            l.lineThickness = event_line_weight
            l.colour = event_line_color
            if not remove_event_labels:
                event_label = event_labels[i]
                if event_labels_top:
                    t = tg.text(event_label, x, top)
                else:
                    t = tg.text(event_label, x, bottom)
                t.textSize = event_font_size
                if event_labels_top:
                    t.anchor = "west"
                    # t.yShift = -0.1
                    t.rotate = 90
                else:
                    t.anchor = "east"
                    # t.yShift = 0.1
                    t.rotate = 90
            if i in splittable_indices:
                if event_labels_top:
                    split_t = tg.text('*', x, bottom)
                    # split_t.yShift = 0.1
                    split_t.anchor = "north"
                else:
                    split_t = tg.text('*', x, top)
                    split_t.yShift = -0.3
                    split_t.anchor = "south"
                split_t.textSize = "LARGE"
                line_labels.append(split_t)
            line_labels.append(t)
            lines.append(l)
        if not show_events:
            for l in lines:
                l.colour = 'black!00'
        if not show_event_labels:
            for l in line_labels:
                l.colour = 'black!00'
        tg.epdf()


def main_cli():
    out_dir = os.path.join(os.path.pardir, 'trees')
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    t = RevJumpMoveTreeInfo(
            node_heights = {
                    '12'    : 0.2,
                    '34'    : 0.4,
                    '1234'  : 0.7,
                    '56'    : 0.4,
                    '78'    : 0.2,
                    '789'   : 0.2,
                    '56789' : 1.0,
                    'root'  : 1.0,
                    },
            height_multiplier = 7.0)
    event_indices_to_labels = {
            0: r'$\tau_{\scriptscriptstyle 0}$',
            1: r'$\tau_{\scriptscriptstyle 1} \sim \textrm{\sffamily Beta}(\tau_{\scriptscriptstyle 0}, \tau_{\scriptscriptstyle 3}, \alpha, \beta)$',
            2: r'$\tau_{\scriptscriptstyle 2} \sim \textrm{\sffamily Beta}(\tau_{\scriptscriptstyle 0}, \tau_{\scriptscriptstyle 3}, \alpha, \beta)$',
            3: r'$\tau_{\scriptscriptstyle 3} \sim \textrm{\sffamily Beta}(\tau_{\scriptscriptstyle 0}, \tau_{\scriptscriptstyle 4}, \alpha, \beta)$',
            4: r'$\tau_{\scriptscriptstyle 4} \sim \textrm{\sffamily Gamma}(k, \theta)$',
            }
    t.plot_tree(
            base_name = 'rj-move-tree',
            dir_name = out_dir,
            scale = 1.0,
            yscale = 1.0,
            vmargin = 0,
            lmargin = 0,
            phylo_line_weight = 'ultra thick',
            event_line_weight = 'thick',
            event_line_color = 'black!50',
            event_line_style = 'dashed',
            title_str = None,
            show_title = False,
            title_font_size = 'huge',
            show_events = True,
            splittable_indices = [1, 2, 4],
            event_indices_to_labels = event_indices_to_labels,
            show_event_labels = True,
            remove_event_labels = False,
            event_labels_top = True,
            event_font_size = 'large',
            show_node_labels = True,
            show_tip_labels = True,
            time_font_size = 'large')

    t = TreeInfo(
            node_heights = {
                    'ytv': 0.05,
                    'ytp': 0.08,
                    'ytg': 0.05,
                    'ytav': 0.08,
                    'ytap': 0.08,
                    'ytag': 0.05,
                    'pg': 0.2,
                    'root': 0.2,
                    },
            height_multiplier = 40.0)
    t.plot_tree(
            base_name = 'gecko-generalized-tree',
            dir_name = out_dir,
            scale = 1.0,
            yscale = 0.8,
            image_width = 12,
            vmargin = 0,
            lmargin = 0,
            phylo_line_weight = 'ultra thick',
            event_line_weight = 'thick',
            event_line_color = 'black!50',
            event_line_style = 'dashed',
            show_title = False,
            title_str = '{\\sffamily True history}',
            show_events = True,
            show_event_labels = True,
            event_labels_top = True,
            show_node_labels = True,
            event_label_plots = True,
            plot_width = 16,
            plot_y_shift = -1.55,
            plot_suffix = "-vln.pdf")
    t.plot_tree(
            base_name = 'gecko-generalized-tree-all-sites',
            dir_name = out_dir,
            scale = 1.0,
            yscale = 0.63,
            image_width = 10,
            vmargin = 0,
            lmargin = 0,
            phylo_line_weight = 'ultra thick',
            event_line_weight = 'thick',
            event_line_color = 'black!50',
            event_line_style = 'dashed',
            show_title = False,
            title_str = '{\\sffamily True history}',
            show_events = True,
            event_indices = (0, 1),
            show_event_labels = True,
            event_labels_top = True,
            show_node_labels = True,
            event_label_plots = True,
            # show_tip_images = False,
            show_times = True,
            time_font_size = 'small',
            plot_width = 10,
            plot_x_shift = -0.02,
            plot_y_shift = 0.17,
            plot_suffix = "-all-sites-hist-cropped.pdf")
            # plot_suffix = "-all-sites-vln-cropped.pdf")

    t = TreeInfo(
            node_heights = {
                    'ytv': 0.05,
                    'ytp': 0.11,
                    'ytg': 0.03,
                    'ytav': 0.06,
                    'ytap': 0.14,
                    'ytag': 0.12,
                    'pg': 0.18,
                    'root': 0.2,
                    },
            height_multiplier = 40.0,
            tree = '(((({yv}:{lyv},{tv}:{ltv}):{lytv},{av}:{lav}):{lv},(({yp}:{lap},({tp}:{ltp},{ap}:{lyp}):{lytp}):{lp},(({yg}:{lyg},{tg}:{ltg}):{lytg},{ag}:{lag}):{lg}):{lgp}))root:{lroot};'
            )
    t.plot_tree(
            base_name = 'gecko-bifurcating-tree',
            dir_name = out_dir,
            scale = 1.0,
            yscale = 0.72,
            image_width = 12,
            vmargin = 0,
            lmargin = 15,
            phylo_line_weight = 'ultra thick',
            event_line_weight = 'thick',
            event_line_color = 'black!50',
            event_line_style = 'dashed',
            show_title = False,
            title_str = '{\\sffamily Current tree model}',
            show_events = True,
            show_event_labels = True,
            event_labels_top = False,
            show_node_labels = True,
            show_tip_images = False)

if __name__ ==  '__main__':
    main_cli()

