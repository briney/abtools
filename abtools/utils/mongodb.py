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
import subprocess as sp

from pymongo import MongoClient


logger = logging.getLogger('mongodb')


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
	if collection:
		return [collection, ]
	collections = db.collection_names(include_system_collections=False)
	if prefix:
		collections = [c for c in collections if c.startswith(prefix)]
	if suffix:
		collections = [c for c in collections if c.endswith(prefix)]
	return sorted(collections)



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
	if lookup_field and match:
		if type(match) != list:
			match = [match, ]
		query_match = {match_field: {'$in': match}}
	else:
		query_match = {}
	c.update_many(query_match, {field: value}, multi=True)
	conn.close()


def mongoimport(jfile, database, collection, ip='localhost', port=27017, user=None, password=None):
	u = " -u {}".format(user) if user else ""
	p = " -p {}".format(password) if password else ""
	# user_password = "{}{} --authenticationDatabase admin".format(username, password)
	mongo_cmd = "mongoimport --host {}:{}{}{} --db {} --collection {} --file {}".format(
		ip, port, u, p, database, collection, jfile)
	mongo = sp.Popen(mongo_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
	stdout, stderr = mongo.communicate()
	return stdout


def index(db, collection, fields, asc=True):
	import pymongo
	# if not fields:
	# 	logger.critical('ENSURE-INDEX EXCEPTION: either "field" or "fields" must be provided.')
	# 	raise RuntimeError('Either "field" or "fields" must be provided to index')
	_dir = pymongo.ASCENDING if asc else pymongo.DESCENDING
	_dirs = [_dir] * len(fields)
	_fields = zip(fields, _dirs)
	coll = db[collection]
	coll.create_index(_fields)


def remove_padding(db, collection, field='padding'):
	'''
	Removes a padding field from all records in ::collection::

	Inputs are a pymongo Database object, a collection name, and an optional
	padding ::field:: name (default is 'padding').
	'''
	c = db[collection]
	print_remove_padding()
	c.update({}, {'$unset': {field: ''}}, multi=True)
