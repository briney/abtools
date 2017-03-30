#!/usr/bin/python
# filename: timepoint.py


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


class Timepoint(object):
    """
    Stores and manipulates timepoint information.
    """
    def __init__(self, tp, order, color):
        super(Timepoint, self).__init__()
        self.name = tp
        self.order = int(order)
        self.raw_color = color
        self.color = self._get_color()


    def _get_color(self):
        if self._is_rbg():
            return self._convert_to_hex()
        return self.raw_color


    def _is_rbg(self):
        color = self.raw_color
        if type(color) == tuple:
            return True
        elif len(self.raw_color.split(',')) == 3:
            self.raw_color = self._convert_to_tuple()
            return True
        return False


    def _convert_to_tuple(self):
        c = self.raw_color.split(',')
        r = float(c[0].replace('(', '').strip())
        g = float(c[1].strip())
        b = float(c[2].replace(')', '').strip())
        return (r, g, b)


    def _convert_to_hex(self):
        r = float(self.raw_color[0])
        g = float(self.raw_color[1])
        b = float(self.raw_color[2])
        if 0 < round(sum([r, b, g]), 2) <= 3:
            r, b, g = (255 * r, 255 * b, 255 * g)
        return '#%02x%02x%02x' % (r, b, g)
