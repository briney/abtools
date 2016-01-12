#!/usr/bin/env python
# filename: mongodb.py


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
import subprocess as sp

from pymongo import MongoClient

from abtools.utils import log


def get_db(db, ip='localhost', port=27017, user=None, password=None):
	'''
	Returns a pymongo Database object.

	Inputs are the database name (::db::) and optional connection options.
	Default connection options are:
		ip: "localhost"
		port: 27017
		user: None
		password: None

	Note that ::user:: and ::password:: are only required when connecting to a
	MongoDB database that has authentication enabled.
	'''
	import urllib
	if user and password:
		pwd = urllib.quote_plus(password)
		uri = 'mongodb://{}:{}@{}:{}'.format(user, pwd, ip, port)
		conn = MongoClient(uri)
	else:
		conn = MongoClient(ip, port)
	return conn[db]


def get_collections(db, collection=None, prefix=None, suffix=None):
	'''
	Returns a sorted list of collection names found in ::db::

	Inputs are a pymongo Database object (::db::) and, optionally,
	a prefix and/or suffix.

	If prefix or suffix are provided, only collections
	containing the prefix/suffix will be returned.
	'''
	if collection is not None:
		return [collection, ]
	collections = db.collection_names(include_system_collections=False)
	if prefix:
		collections = [c for c in collections if c.startswith(prefix)]
	if suffix:
		collections = [c for c in collections if c.endswith(prefix)]
	return sorted(collections)


def rename_collection(db, collection, new_name):
	'''
	Renames a MongoDB collection.

	Inputs are a pymongo database connection object, the
	name of the collection to be renamed, and the new name.

	The new name can be provided either as a string or as a function
	that will be applied to the old collection name. If a function is
	provided and calling the function on the current collection name
	results in an empty string, the collection will not be renamed.
	'''
	if isinstance(new_name, function):
		_new = new_name(collection)
		if _new == '':
			return
	else:
		_new = new_name
	c = db[collection]
	c.rename(_new)


def update(field, value, db, collection, match=None):
	'''
	Updates records to set ::field:: equal to ::value:: for all
	records that meet ::match:: criteria.

	::db:: should be a pymongo Database object.
	::collection:: should be the name of a collection in db, as a string
	::match:: should be a dict containing the pymongo query criteria, like:
		{'seq_id': {'$in': ['a', 'b', 'c']}, 'size': {'$gte': 25}, 'prod': 'yes'}
	'''
	c = db[collection]
	match = match if match is not None else {}
	c.update(match, {'$set': {field: value}}, multi=True)
	# below is for MongoDB 3.0+
	# c.update_many(match, {'$set': {field: value}})


def unset(db, collection, field, match=None):
	'''
	Removes a field from all records in ::collection:: that meet
	::match:: criteria

	Inputs are a pymongo Database object, a collection name, and an optional
	::match::, a dict of the format:
		{'seq_id': {'$in': ['a', 'b', 'c']}, 'size': {'$gte': 25}, 'prod': 'yes'}
	'''
	c = db[collection]
	match = match if match is not None else {}
	c.update(match, {'$unset': {field: ''}}, multi=True)
	# below is for MongoDB 3.0+
	# c.update_many(match, {'$unset': {field: ''}})


def mongoimport(json, database,
				ip='localhost', port=27017,
				user=None, password=None,
				delim='_', delim1=None, delim2=None,
				delim_occurance=1, delim1_occurance=1, delim2_occurance=1):
	'''
	Performs mongoimport on one or more json files.

	Required inputs:
		::json:: can be one of several things:
			-- a single JSON file
			-- an iterable (list or tuple) of one or more JSON files
			-- a directory, containing one or more JSON files
		::database:: is the name of the MongoDB database into which
			the JSON files will be imported
	'''
	logger = log.get_logger('mongodb')
	_print_mongoimport_info(logger)
	if type(json) in (list, tuple):
		pass
	elif os.path.isdir(json):
		from abtools.utils.pipeline import list_files
		json = list_files(json)
	else:
		json = [json, ]
	jsons = sorted([os.path.expanduser(j) for j in json if j.endswith('.json')])
	collections = _get_import_collections(jsons, delim, delim_occurance,
										  delim1, delim1_occurance,
										  delim2, delim2_occurance)
	logger.info('Found {} files to import'.format(len(jsons)))
	logger.info('')
	for i, (json_file, collection) in enumerate(zip(jsons, collections)):
		logger.info('[ {} ] {} --> {}'.format(i + 1, os.path.basename(json_file), collection))
		# logger.info("Performing mongoimport on {}.".format(os.path.basename(json_file)))
		# logger.info("Importing the file into collection {}.".format(collection))
		if all([user, password]):
			host = '--host {} --port {} -username {} -password {}'.format(ip, port, user, password)
		else:
			host = '--host {} --port {}'.format(ip, port)
		mongo_cmd = "mongoimport {} --db {} --collection {} --file {}".format(
			host, database, collection, json_file)
		mongo = sp.Popen(mongo_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
		stdout, stderr = mongo.communicate()


def index(db, collection, fields, directions=None, desc=False):
	'''
	Builds a simple (single field) or complex (multiple fields) index
	on a single collection in a MongoDB database.

	::db:: should be a pymongo connection object to the desired MongoDB database
	::collection:: should be the name of the collection to be indexed
	::fields:: is the field(s) to be indexed, and can be one of two things:
		1) a single field, as a string
		2) one or more fields, as an iterable (list or tuple)

	By default, all fields will be indexed in ascending order. Setting ::desc:: to True
	will result in the index being built in DESCENDING order. For multi-directional indexes,
	a list of pymongo direction objects (pymongo.ASCENDING and pymongo.DESCENDING) can be provided,
	in the same order as the list of fields to be indexed.
	'''
	import pymongo
	if type(fields) == str:
		fields = [fields, ]
	if directions is None:
		_dir = pymongo.DESCENDING if desc else pymongo.ASCENDING
		directions = [_dir] * len(fields)
	field_tuples = zip(fields, directions)
	coll = db[collection]
	coll.create_index(field_tuples)


def remove_padding(db, collection, field='padding'):
	'''
	Removes a padding field from all records in ::collection::

	Inputs are a pymongo Database object, a collection name, and an optional
	padding ::field:: name (default is 'padding').
	'''
	unset(db, collection, field=field)
	# c = db[collection]
	# c.update({}, {'$unset': {field: ''}}, multi=True)


def _get_import_collections(jsons, delim, delim_occurance,
							delim1, delim1_occurance,
							delim2, delim2_occurance):
	jnames = [os.path.basename(j) for j in jsons]
	if not all([delim1, delim2]):
		collections = [delim.join(j.split(delim)[:delim_occurance]) for j in jnames]
	else:
		pre_colls = [delim1.join(j.split(delim1)[delim1_occurance:]) for j in jnames]
		collections = [delim2.join(j.split(delim2)[:delim2_occurance]) for j in pre_colls]
	return collections


def _print_mongoimport_info(logger):
	logger.info('')
	logger.info('')
	logger.info('')
	logger.info('-' * 25)
	logger.info('MONGOIMPORT')
	logger.info('-' * 25)
	logger.info('')


def _print_remove_padding():
	logger = log.get_logger('mongodb')
	logger.info('Removing MongoDB padding...')
