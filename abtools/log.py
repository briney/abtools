#!/usr/bin/env python
# filename: log.py


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


from __future__ import print_function

import logging
import os

# from .pipeline import make_dir


def setup_logging(logfile, print_log_location=True, debug=False):
    '''
    Set up logging using the built-in ``logging`` package.

    A stream handler is added to all logs, so that logs at or above
    ``logging.INFO`` level are printed to screen as well as written
    to the log file.

    Arguments:

        logfile (str): Path to the log file. If the parent directory
            does not exist, it will be created. Required.

        print_log_location (bool): If ``True``, the log path will be
            written to the log upon initialization. Default is ``True``.

        debug (bool): If true, the log level will be set to ``logging.DEBUG``.
            If ``False``, the log level will be set to ``logging.INFO``.
            Default is ``False``.
    '''
    log_dir = os.path.dirname(logfile)
    make_dir(log_dir)
    fmt = '[%(levelname)s] %(name)s %(asctime)s %(message)s'
    if debug:
        logging.basicConfig(filename=logfile,
                            filemode='w',
                            format=fmt,
                            level=logging.DEBUG)
    else:
        logging.basicConfig(filename=logfile,
                            filemode='w',
                            format=fmt,
                            level=logging.INFO)
    logger = logging.getLogger('log')
    logger = add_stream_handler(logger)
    if print_log_location:
        logger.info('LOG LOCATION: {}'.format(logfile))


def get_logger(name=None):
    '''
    Get a logging handle.

    As with ``setup_logging``, a stream handler is added to the
    log handle.

    Arguments:

        name (str): Name of the log handle. Default is ``None``.
    '''
    logger = logging.getLogger(name)
    if len(logger.handlers) == 0:
        logger = add_stream_handler(logger)
    return logger


def add_stream_handler(logger):
    formatter = logging.Formatter("%(message)s")
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)
    return logger


def make_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)
