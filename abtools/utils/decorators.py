#!/usr/bin/env python
# filename: decorators.py


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


def lazy_property(func):
    '''
    Wraps a property to provide lazy evaluation. Eliminates boilerplate.
    Also provides for setting and deleting the property.

    Use as you would use the @property decorator::

        # OLD:
        class MyClass():
            def __init__():
                self._compute = None

            @property
            def compute(self):
                if self._compute is None:
                    # computationally intense stuff
                    # ...
                    # ...
                    self._compute = result
                return self._compute

            @compute.setter
            def compute(self, value):
                self._compute = value


        # NEW:
        class MyClass():

            def __init__():
                pass

            @lazy_property
            def compute(self):
                # computationally intense stuff
                # ...
                # ...
                return result

    .. note:

        Properties wrapped with ``lazy_property`` are only evaluated once.
        If the instance state changes, lazy properties will not be automatically
        re-evaulated and the update must be explicitly called for::

            c = MyClass(data)
            prop = c.lazy_property

            # If you update some data that affects c.lazy_property
            c.data = new_data

            # c.lazy_property won't change
            prop == c.lazy_property  # TRUE

            # If you want to update c.lazy_property, you can delete it, which will
            # force it to be recomputed (with the new data) the next time you use it
            del c.lazy_property
            new_prop = c.lazy_property
            new_prop == prop  # FALSE
    '''
    attr_name = '_lazy_' + func.__name__

    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, func(self))
        return getattr(self, attr_name)

    @_lazy_property.deleter
    def _lazy_property(self):
        if hasattr(self, attr_name):
            delattr(self, attr_name)

    @_lazy_property.setter
    def _lazy_property(self, value):
        setattr(self, attr_name, value)

    return _lazy_property


def coroutine(func):
    '''
    Initializes a coroutine -- essentially it just takes a
    generator function and calls generator.next() to get
    things going.
    '''
    def start(*args, **kwargs):
        cr = func(*args, **kwargs)
        cr.next()
        return cr
    return start
