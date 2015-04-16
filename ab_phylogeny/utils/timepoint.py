#!/usr/bin/python
# filename: timepoint.py



###########################################################################
#
# Copyright (c) 2014 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT)
#
###########################################################################



class Timepoint(object):
	"""
	Stores and manipulates timepoint information.
	"""
	def __init__(self, tp, order, color):
		super(Timepoint, self).__init__()
		self.name = tp
		self.order = int(order)
		self.raw_color = color
		self.color = _get_color()


	def _get_color(self):
		if _is_rgb():
			return _convert_to_hex()
		return self.raw_color


	def _is_rbg(self):
		color = self.raw_color
		if color.split(',') == 3:
			return True
		return False


	def _convert_to_hex(self):
		r = float(self.raw_color.split(',')[0].replace('(', '').strip())
		g = float(self.raw_color.split(',')[1].strip())
		b = float(self.raw_color.split(',')[2].replace(')', '').strip())
		if 0 < round(sum([r, b, g]), 2) <= 3:
			rbg = [256 * r, 256 * b, 256 * g]
		else:
			rbg = [r, b, g]
		return '#%02x%02x%02x' % rgb
