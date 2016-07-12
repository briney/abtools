#!/usr/bin/env python
# filename: color.py


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


from __future__ import print_function, division

import seaborn as sns

import numpy as np

from matplotlib import pyplot as plt
from matplotlib import cm, colors


cmaps = {'heatmap': sns.diverging_palette(240, 10, as_cmap=True)}


def cmap_from_color(color, dark=False):
    '''
    Generates a matplotlib colormap from a single color.

    Colormap will be built, by default, from white to ``color``.

    Args:

        color: Can be one of several things:

            1. Hex code
            2. HTML color name
            3. RGB tuple

        dark (bool): If ``True``, colormap will be built from ``color`` to
            black. Default is ``False``, which builds a colormap from
            white to ``color``.

    Returns:

        colormap: A matplotlib colormap

    '''
    if dark:
        return sns.dark_palette(color, as_cmap=True)
    else:
        return sns.light_palette(color, as_cmap=True)


def hex_to_rgb(hex_string):
    rgb = colors.hex2color(hex_string)
    return tuple([int(255 * x) for x in rgb])


def rgb_to_hex(rgb_tuple):
    div = 1 if all([v <= 1.0 for v in rgb_tuple]) else 255
    return colors.rgb2hex([1.0 * x / div for x in rgb_tuple])


def hls(n_colors, hue=0.01, lightness=0.6, saturation=0.65):
    return sns.hls_palette(n_colors, h=hue, l=lightness, s=saturation)


def husl(n_colors, hue=0.01, saturation=0.9, lightness=0.65):
    return sns.husl_palette(n_colors, h=hue, s=saturation, l=lightness)


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=256):
    """
    Truncates a colormap, such that the new colormap consists of
    ``cmap[minval:maxval]``.

    If maxval is larger than minval, the truncated colormap will be reversed.

    Args:

    	cmap (colormap): Colormap to be truncated

    	minval (float): Lower bound. Should be a float betwee 0 and 1.

    	maxval (float): Upper bound. Should be a float between 0 and 1

    	n (int): Number of colormap steps. Default is ``256``.

    Returns:

    	colormap


    http://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib
    """
    cmap = get_cmap(cmap)
    name = "%s-trunc-%.2g-%.2g" % (cmap.name, minval, maxval)
    return colors.LinearSegmentedColormap.from_list(
        name, cmap(np.linspace(minval, maxval, n)))


def stack_colormap(lower, upper, n=256):
    """
    Stacks two colormaps (``lower`` and ``upper``) such that
    low half -> ``lower`` colors, high half -> ``upper`` colors

    Args:

    	lower (colormap): colormap for the lower half of the stacked colormap.

    	upper (colormap): colormap for the upper half of the stacked colormap.

    	n (int): Number of colormap steps. Default is ``256``.
    """
    A = get_cmap(lower)
    B = get_cmap(upper)
    name = "%s-%s" % (A.name, B.name)
    lin = np.linspace(0, 1, n)
    return array_cmap(np.vstack((A(lin), B(lin))), name, n=n)


def get_cmap(cmap=None, name=None, from_color=None, dark=False, n=256):
    # """
    # Generates a matplotlib colormap.

    # cmap can be one of several things:
    #     - a name ('Blues', 'BuGn_r') of a built-in colormap
    #     - a cmap
    #     - a filename, np.loadtxt() n x 3 or 4  ints 0..255 or floats 0..1
    #     - a numpy array
    # See http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps or in IPython, plt.cm.<tab>

    # An optional name for the colormap can be provided (name), as well as the number
    # of cmap steps (n).

    # Alternatively, to make a colormap using a single color, provide the color
    # via from_color. By default, the supplied color will be the dark end of
    # the cmap, with white as the lightest color. To reverse (use the input
    # color as the lighest value and black as the darkest), set dark=True.

    # Returns a matplotlib colormap object.
    # """
    if from_color is not None:
        return cmap_from_color(from_color, dark)
    elif cmap is None:
        err = 'You must provide either cmap or from_color'
        raise RuntimeError(err)
    if isinstance(cmap, colors.Colormap):
        return cmap
    if isinstance(cmap, basestring):
        if cmap in cm.cmap_d:
            return plt.get_cmap(cmap)  # "Blues" ...
        A = np.loadtxt(cmap, delimiter=None)  # None: white space
        name = name or cmap.split("/")[-1].split(".")[0]  # .../xx.csv -> xx
    else:
        A = cmap  # numpy array or array-like
    return array_cmap(A, name, n=n)
