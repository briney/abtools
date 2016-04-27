#!/usr/bin/env python
# filename: s3.py


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
import tarfile

from abtools import log


def compress_and_upload(data, compressed_file, s3_path,
    method='gz', delete=False, access_key=None, secret_key=None):
    '''
    Compresses data and uploads to S3.

    S3 upload uses ``s3cmd``, so you must either:

        1) Manually configure ``s3cmd`` prior to use (typically using ``s3cmd --configure``).

        2) Configure ``s3cmd`` using ``s3.configure()``.

        3) Pass your access key and secret key to ``compress_and_upload``, which will automatically configure s3cmd.

    .. note:

        ``s3cmd`` configuration only needs to be done once per computer,
        which means that relaunching a cloud instance or Docker image will
        require re-configuration of ``s3cmd``.

    Args:

        data: Can be one of three things:

            1) Path to a single file

            2) Path to a directory

            3) A list of one or more paths to files or directories

        compressed_file (str): Path to the compressed file. Required.

        s3_path (str): The S3 path, with the filename omitted. The S3 filename
          will be the basename of the ``compressed_file``. For example::

            compress_and_upload(data='/path/to/data',
                                compressed_file='/path/to/compressed.tar.gz',
                                s3_path='s3://my_bucket/path/to/')

          will result in an uploaded S3 path of ``s3://my_bucket/path/to/compressed.tar.gz``

        method (str): Compression method. Options are ``'gz'`` (gzip) or ``'bz2'`` (bzip2).
            Default is ``'gz'``.

        delete (bool): If ``True``, the ``compressed_file`` will be deleted after upload
            to S3. Default is ``False``.

        access_key (str): AWS access key.

        secret_key (str): AWS secret key.
    '''
    logger = log.get_logger('s3')
    if all([access_key, secret_key]):
        configure(access_key=access_key, secret_key=secret_key, logger=logger)
    compress(data, compressed_file, compress=method, logger=logger)
    put(compressed_file, s3_path, logger=logger)
    if delete:
        os.unlink(compressed_file)


def put(f, s3_path, logger=None):
    '''
    Uploads a single file to S3, using s3cmd.

    Args:

        f (str): Path to a single file.

        s3_path (str): The S3 path, with the filename omitted. The S3 filename
            will be the basename of the ``f``. For example::

                put(f='/path/to/myfile.tar.gz', s3_path='s3://my_bucket/path/to/')

            will result in an uploaded S3 path of ``s3://my_bucket/path/to/myfile.tar.gz``
    '''
    if not logger:
        logger = log.get_logger('s3')
    fname = os.path.basename(f)
    target = os.path.join(s3_path, fname)
    s3cmd_cline = 's3cmd put {} {}'.format(f, target)
    print_put_info(fname, target, logger)
    s3cmd = sp.Popen(s3cmd_cline,
                     stdout=sp.PIPE,
                     stderr=sp.PIPE,
                     shell=True)
    stdout, stderr = s3cmd.communicate()


def print_put_info(fname, target, logger):
    logger.info('')
    logger.info('')
    logger.info('')
    logger.info('-' * 25)
    logger.info('UPLOADING TO S3')
    logger.info('-' * 25)
    logger.info('')
    logger.info('File: {}'.format(fname))
    logger.info('Target S3 location: {}'.format(target))



def compress(d, output, compress='gz', logger=None):
    '''
    Creates a compressed/uncompressed tar file.

    Args:

        d: Can be one of three things:

            1. the path to a single file, as a string

            2. the path to a single directory, as a string

            3. an iterable of file or directory paths

        output (str): Output file path.

        compress: Compression method. Options are ``'gz'`` (gzip),
            ``'bz2'`` (bzip2) and ``'none'`` (uncompressed). Default is ``'gz'``.
    '''
    if not logger:
        logger = log.get_logger('s3')
    if type(d) not in [list, tuple]:
        d = [d, ]
    d = [os.path.expanduser(_d) for _d in d]
    print_compress_info(d, output, compress, logger)
    if compress.lower() == 'none':
        compress = ''
    elif compress.lower() not in ['gz', 'bz2']:
        logger.info('Compression option ("{}") is invalid.\nFalling back to uncompressed.'.format(compress))
        compress = ''
    output = os.path.expanduser(output)
    tar = tarfile.open(output, 'w:{}'.format(compress))
    for obj in d:
        tar.add(obj)
    tar.close()
    return output


def print_compress_info(d, output, compress, logger):
    if not logger:
        logger = log.get_logger('s3')
    dirs = [obj for obj in d if os.path.isdir(obj)]
    files = [obj for obj in d if os.path.isfile(obj)]
    logger.info('')
    logger.info('')
    logger.info('')
    logger.info('-' * 25)
    logger.info('COMPRESSING DATA')
    logger.info('-' * 25)
    logger.info('')
    logger.info('Ouptut file: {}'.format(output))
    logger.info('Compression: {}'.format(compress.lower()))
    if dirs:
        d = 'directories' if len(dirs) > 1 else 'directory'
        logger.info('Found {} {} to compress: {}'.format(len(dirs), d,
                                                         ', '.join(dirs)))
    if files:
        f = 'files' if len(files) > 1 else 'file'
        logger.info('Found {} {} to compress: {}'.format(len(files), f,
                                                         ', '.join(files)))


def configure(access_key=None, secret_key=None, logger=None):
    '''
    Configures s3cmd prior to first use.

    If no arguments are provided, you will be prompted to enter
    the access key and secret key interactively.

    Args:

        access_key (str): AWS access key

        secret_key (str): AWS secret key
    '''
    if not logger:
        logger = log.get_logger('s3')
    if not all([access_key, secret_key]):
        logger.info('')
        access_key = raw_input('AWS Access Key: ')
        secret_key = raw_input('AWS Secret Key: ')
    _write_config(access_key, secret_key)
    logger.info('')
    logger.info('Completed writing S3 config file.')
    logger.info('')


def _write_config(access_key, secret_key):
    cfg_string = '[default]\n'
    cfg_string += 'access_key = {}\n'.format(access_key)
    cfg_string += 'secret_key = {}\n'.format(secret_key)
    cfg_string += CONFIG_DEFAULTS
    cfg_file = os.path.expanduser('~/.s3cfg')
    open(cfg_file, 'w').write(cfg_string)


CONFIG_DEFAULTS = '''
access_token =
add_encoding_exts =
add_headers =
bucket_location = US
cache_file =
cloudfront_host = cloudfront.amazonaws.com
default_mime_type = binary/octet-stream
delay_updates = False
delete_after = False
delete_after_fetch = False
delete_removed = False
dry_run = False
enable_multipart = True
encoding = UTF-8
encrypt = False
expiry_date =
expiry_days =
expiry_prefix =
follow_symlinks = False
force = False
get_continue = False
gpg_command = /usr/bin/gpg
gpg_decrypt = %(gpg_command)s -d --verbose --no-use-agent --batch --yes --passphrase-fd %(passphrase_fd)s -o %(output_file)s %(input_file)s
gpg_encrypt = %(gpg_command)s -c --verbose --no-use-agent --batch --yes --passphrase-fd %(passphrase_fd)s -o %(output_file)s %(input_file)s
gpg_passphrase =
guess_mime_type = True
host_base = s3.amazonaws.com
host_bucket = %(bucket)s.s3.amazonaws.com
human_readable_sizes = False
ignore_failed_copy = False
invalidate_default_index_on_cf = False
invalidate_default_index_root_on_cf = True
invalidate_on_cf = False
list_md5 = False
log_target_prefix =
max_delete = -1
mime_type =
multipart_chunk_size_mb = 250
preserve_attrs = True
progress_meter = True
proxy_host =
proxy_port = 0
put_continue = False
recursive = False
recv_chunk = 4096
reduced_redundancy = False
restore_days = 1
send_chunk = 4096
server_side_encryption = False
simpledb_host = sdb.amazonaws.com
skip_existing = False
socket_timeout = 300
urlencoding_mode = normal
use_https = True
use_mime_magic = True
verbosity = WARNING
website_endpoint = http://%(bucket)s.s3-website-% (location)s.amazonaws.com/
website_error =
website_index = index.html
'''
