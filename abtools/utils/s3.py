#!/usr/bin/env python
# filename: scpcr_demultiplexing.py


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


import logging
import os
import subprocess as sp
import tarfile

logger = logging.getLogger('s3')


def put(f, s3_path):
	'''
	Uploads a file to S3, using s3cmd.

	Inputs: path to a single file, typically compressed using s3.compress() and
	the s3 path into which the file should be uploaded. The s3 path should not already
	include file name.
	'''
	fname = os.path.basename(f)
	target = os.path.join(s3_path, fname)
	s3cmd_cline = 's3cmd put {} {}'.format(f, target)
	print_put_info(fname, target)
	s3cmd = sp.Popen(s3cmd_cline,
					 stdout=sp.PIPE,
					 stderr=sp.PIPE,
					 shell=True)
	stdout, stderr = s3cmd.communicate()
	logger.info('')
	logger.info('')


def print_put_info(fname, target):
	logger.info('-' * 25)
	logger.info('UPLOADING TO S3')
	logger.info('-' * 25)
	logger.info('')
	logger.info('File: {}'.format(fname))
	logger.info('Target S3 location: {}'.format(target))



def compress(d, output, compress='gz'):
	'''
	Creates a compressed tar file.

	::d:: can be one of four things:
		1) the path to a single file, as a string
		2) the path to a single directory, as a string
		3) an iterable of files
		4) an iterable of directories

	::compress:: can be one of three things:
		1) 'gz' for gzip compression (default)
		2) 'bz2' for bzip2 compression
		3) 'none' for no compression
	'''
	if type(d) == str:
		d = [d, ]
	print_compress_info(d, output, compress)
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


def print_compress_info(d, output, compress):
	dirs = [obj for obj in d if os.path.isdir(obj)]
	files = [obj for obj in d if os.path.isfile(obj)]
	logger.info('-' * 25)
	logger.info('COMPRESSING DATA')
	logger.info('-' * 25)
	logger.info('')
	logger.info('Ouptut file: {}'.format(output))
	logger.info('Compression = {}'.format(compress.lower()))
	logger.info('Found {} directories to compress: {}'.format(len(dirs), ', '.join(dirs)))
	logger.info('Found {} files to compress: {}'.format(len(files), ', '.join(files)))
	logger.info('')


def configure():
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
