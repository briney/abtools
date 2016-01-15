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
	if dark:
		return sns.dark_palette(color, as_cmap=True)
	else:
		return sns.light_palette(color, as_cmap=True)


def hls(n_colors, hue=0.01, lightness=0.6, saturation=0.65):
	return sns.hls_palette(n_colors, h=hue, l=lightness, s=saturation)


def husl(n_colors, hue=0.01, saturation=0.9, lightness=0.65):
	return sns.husl_palette(n_colors, h=hue, s=saturation, l=lightness)


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=256):
	"""
	Truncates the provided colormap, such that the new colormap consists of
	cmap[minval:maxval].

	minval and maxval should be floats 0..1

	if maxval is larger than minval, the subset will be reversed.

	http://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib
	"""
	cmap = get_cmap(cmap)
	name = "%s-trunc-%.2g-%.2g" % (cmap.name, minval, maxval)
	return colors.LinearSegmentedColormap.from_list(
		name, cmap(np.linspace(minval, maxval, n)))


def stack_colormap(A, B, n=256):
	"""
	Stacks two colormaps (A, B) such that
	low half -> A colors, high half -> B colors

	Optionally, provide the number of steps (n, default is 256)
	"""
	A = get_cmap(A)
	B = get_cmap(B)
	name = "%s-%s" % (A.name, B.name)
	lin = np.linspace(0, 1, n)
	return array_cmap(np.vstack((A(lin), B(lin))), name, n=n)


def get_cmap(cmap=None, name=None, from_color=None, dark=False, n=256):
	"""
	Generates a matplotlib colormap.

	cmap can be one of several things:
		- a name ('Blues', 'BuGn_r') of a built-in colormap
		- a cmap
		- a filename, np.loadtxt() n x 3 or 4  ints 0..255 or floats 0..1
		- a numpy array
	See http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps or in IPython, plt.cm.<tab>

	An optional name for the colormap can be provided (name), as well as the number
	of cmap steps (n).

	Alternatively, to make a colormap using a single color, provide the color
	via from_color. By default, the supplied color will be the dark end of
	the cmap, with white as the lightest color. To reverse (use the input
	color as the lighest value and black as the darkest), set dark=True.

	Returns a matplotlib colormap object.
	"""
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
