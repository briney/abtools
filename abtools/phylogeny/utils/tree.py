#!/usr/bin/python
# filename: tree.py


#
# Copyright (c) 2015 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


import os
import subprocess as sp

import ete2



def make_tree(alignment, timepoints, delimiter, is_aa, scale, branch_vert_margin,
    	fontsize, show_name, tree_orientation, show_scale=False):
    '''
    Builds a tree file (using FastTree) from a sequence alignment in FASTA format

    Input
    path to a FASTA-formatted sequence alignment

    Output
    path to a Newick-formatted tree file
    '''
    tree = alignment.replace('_aligned.aln', '_tree.nw')
    tree = fast_tree(alignment, tree, is_aa)
    make_figure(tree, timepoints, delimiter, scale, branch_vert_margin,
    	fontsize, show_name, tree_orientation, show_scale=show_scale)



def fast_tree(alignment, tree, is_aa, show_output=False):
    if is_aa:
        ft_cmd = 'fasttree {} > {}'.format(alignment, tree)
    else:
        ft_cmd = 'fasttree -nt {} > {}'.format(alignment, tree)
    ft = sp.Popen(ft_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = ft.communicate()
    if show_output:
    	print(ft_cmd)
    	print(stdout)
    	print(stderr)
    return tree


def make_figure(tree, timepoints, delimiter, scale, branch_vert_margin,
		fontsize, show_name, tree_orientation, show_scale=False):
    fig = tree.replace('_tree.nw', '_tree.pdf')
    orders = {tp.name: tp.order for tp in timepoints}
    colors = {tp.name: tp.color for tp in timepoints}
    # settins for name showing
    if show_name == 'none':
        show_name = []
    if show_name == 'all':
        show_name = ['mab', 'root', 'input']
    elif show_name == 'no-root':
        show_name = ['input', 'mab']
    elif type(show_name) in [str, unicode]:
        show_name = [show_name, ]
    # make the tree
    t = ete2.Tree(tree)
    t.set_outgroup(t&"root")
    # style the nodes based on timepoint
    for node in t.traverse():
        earliest = get_earliest_leaf(node.get_leaf_names(), orders, delimiter)
        color = colors[earliest]
        node_type = get_node_type(node.name)
        style = ete2.NodeStyle()
        style['size'] = 0
        style['vt_line_width'] = 1.0
        style['hz_line_width'] = 1.0
        style['vt_line_color'] = color
        style['hz_line_color'] = color
        style['vt_line_type'] = 0
        style['hz_line_type'] = 0
        if node_type in show_name:
            if node_type in ['mab', 'input']:
                name = ' ' + delimiter.join(node.name.split(delimiter)[1:])
            else:
                name = ' ' + node.name
            tf = ete2.TextFace(name)
            tf.fsize = fontsize
            node.add_face(tf, column=0)
            style['fgcolor'] = '#000000'
        node.set_style(style)
    # style the full tree
    # root = (t&"root")
    # nearest_to_root, distance = root.get_closest_leaf()
    # root_node = t.get_common_ancestor(root, nearest_to_root)
    t.dist = 0
    ts = ete2.TreeStyle()
    ts.orientation = tree_orientation
    ts.show_leaf_name = False
    if scale:
        ts.scale = int(scale)
    if branch_vert_margin:
        ts.branch_vertical_margin = float(branch_vert_margin)
    ts.show_scale = False
    # ladderize
    t.ladderize()
    # render the tree
    t.render(fig, tree_style=ts)


def get_node_type(node_name):
    if node_name == 'root':
        return 'root'
    if node_name.startswith('mab'):
        return 'mab'
    if node_name == 'NoName':
        return 'inner'
    return 'input'


def get_earliest_leaf(leaves, order, delimiter):
    counts = {}
    for leaf in leaves:
        tp = leaf.split(delimiter)[0]
        counts[tp] = counts[tp] + 1 if tp in counts else 1
    total = sum(counts.values())
    if 'root' in counts:
        return 'root'
    timepoints = sorted(counts.keys(), key=lambda x: order[x])
    for tp in timepoints:
        if 100. * counts[tp] / total >= 5:
            return tp
