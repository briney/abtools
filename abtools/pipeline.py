#!/usr/bin/env python
# filename: pipeline.py


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


import glob
import os

from . import log


def initialize(log_file, project_dir=None, debug=False):
    '''
    Initializes an AbTools pipeline.

    Initialization includes printing the AbTools splash, setting up logging,
    creating the project directory, and logging both the project directory
    and the log location.

    Args:

        log_file (str): Path to the log file. Required.

        project_dir (str): Path to the project directory. If not provided,
            the project directory won't be created and the location won't be logged.

        debug (bool): If ``True``, the logging level will be set to ``logging.DEBUG``.
            Default is ``FALSE``, which logs at ``logging.INFO``.

    Returns:

        logger
    '''
    _print_splash()
    log.setup_logging(log_file, print_log_location=False, debug=debug)
    logger = log.get_logger('pipeline')
    if project_dir is not None:
        make_dir(os.path.normpath(project_dir))
        logger.info('PROJECT DIRECTORY: {}'.format(project_dir))
        logger.info('')
    logger.info('LOG LOCATION: {}'.format(log_file))
    print('')
    return logger


def make_dir(d):
    '''
    Makes a directory, if it doesn't already exist.

    Args:

        d (str): Path to a directory.
    '''
    if not os.path.exists(d):
        os.makedirs(d)


def list_files(d, extension=None):
    '''
    Lists files in a given directory.

    Args:

        d (str): Path to a directory.

        extension (str): If supplied, only files that contain the
            specificied extension will be returned. Default is ``False``,
            which returns all files in ``d``.

    Returns:

        list: A sorted list of file paths.
    '''
    if os.path.isdir(d):
        expanded_dir = os.path.expanduser(d)
        files = sorted(glob.glob(expanded_dir + '/*'))
    else:
        files = [d, ]
    if extension is not None:
        if type(extension) in [str, unicode]:
            extension = [extension, ]
        files = [f for f in files if any([f.split('.')[-1] in extension,
                                          f.split('.')[-1].upper() in extension,
                                          f.split('.')[-1].lower() in extension])]
    return files


def _print_splash():
    splash = '''
     _    _   _____           _       ____  _            _ _            
    / \  | |_|_   _|__   ___ | |___  |  _ \(_)_ __   ___| (_)_ __   ___ 
   / _ \ | '_ \| |/ _ \ / _ \| / __| | |_) | | '_ \ / _ \ | | '_ \ / _ \\
  / ___ \| |_) | | (_) | (_) | \__ \ |  __/| | |_) |  __/ | | | | |  __/
 /_/   \_\_.__/|_|\___/ \___/|_|___/ |_|   |_| .__/ \___|_|_|_| |_|\___|
                                             |_|                        '''
    # print('')
    print(splash)
    print('')
