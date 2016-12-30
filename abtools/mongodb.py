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
import platform
import os
import subprocess as sp

from pymongo import MongoClient

from abtools import log


def get_connection(ip='localhost', port=27017, user=None, password=None):
    '''
    Returns a pymongo MongoClient object.

    .. note:

        Both ``user`` and ``password`` are required when connecting to a MongoDB
        database that has authentication enabled.

    Arguments:

        ip (str): IP address of the MongoDB server. Default is ``localhost``.

        port (int): Port of the MongoDB server. Default is ``27017``.

        user (str): Username, if authentication is enabled on the MongoDB database.
            Default is ``None``, which results in requesting the connection
            without authentication.

        password (str): Password, if authentication is enabled on the MongoDB database.
            Default is ``None``, which results in requesting the connection
            without authentication.
    '''
    if platform.system().lower() == 'darwin':
        connect = False
    else:
        connect = True
    if user and password:
        import urllib
        pwd = urllib.quote_plus(password)
        uri = 'mongodb://{}:{}@{}:{}'.format(user, pwd, ip, port)
        return MongoClient(uri, connect=connect)
    return MongoClient(ip, port, connect=connect)


def get_db(db, ip='localhost', port=27017, user=None, password=None):
    '''
    Returns a pymongo Database object.

    .. note:

        Both ``user`` and ``password`` are required when connecting to a MongoDB
        database that has authentication enabled.

    Arguments:

        db (str): Name of the MongoDB database. Required.

        ip (str): IP address of the MongoDB server. Default is ``localhost``.

        port (int): Port of the MongoDB server. Default is ``27017``.

        user (str): Username, if authentication is enabled on the MongoDB database.
            Default is ``None``, which results in requesting the connection
            without authentication.

        password (str): Password, if authentication is enabled on the MongoDB database.
            Default is ``None``, which results in requesting the connection
            without authentication.
    '''
    if platform.system().lower() == 'darwin':
        connect = False
    else:
        connect = True
    if user and password:
        import urllib
        pwd = urllib.quote_plus(password)
        uri = 'mongodb://{}:{}@{}:{}'.format(user, pwd, ip, port)
        conn = MongoClient(uri, connect=connect)
    else:
        conn = MongoClient(ip, port, connect=connect)
    return conn[db]


def get_collections(db, collection=None, prefix=None, suffix=None):
    '''
    Returns a sorted list of collection names found in ``db``.

    Arguments:

        db (Database): A pymongo Database object. Can be obtained
            with ``get_db``.

        collection (str): Name of a collection. If the collection is
            present in the MongoDB database, a single-element list will
            be returned with the collecion name. If not, an empty list
            will be returned. This option is primarly included to allow
            for quick checking to see if a collection name is present.
            Default is None, which results in this option being ignored.

        prefix (str): If supplied, only collections that begin with
            ``prefix`` will be returned.

        suffix (str): If supplied, only collections that end with
            ``suffix`` will be returned.

    Returns:

        list: A sorted list of collection names.
    '''
    if collection is not None:
        return [collection, ]
    collections = db.collection_names(include_system_collections=False)
    if prefix is not None:
        collections = [c for c in collections if c.startswith(prefix)]
    if suffix is not None:
        collections = [c for c in collections if c.endswith(suffix)]
    return sorted(collections)


def rename_collection(db, collection, new_name):
    '''
    Renames a MongoDB collection.

    Arguments:

        db (Database): A pymongo Database object. Can be obtained
            with ``get_db``.

        collection (str): Name of the collection to be renamed.

        new_name (str, func): ``new_name`` can be one of two things::

            1. The new collection name, as a string.
            2. A function which, when passed the current collection name,
                returns the new collection name. If the function
                returns an empty string, the collection will not be
                renamed.
    '''
    if hasattr(new_name, '__call__'):
        _new = new_name(collection)
        if _new == '':
            return
    else:
        _new = new_name
    c = db[collection]
    c.rename(_new)


def update(field, value, db, collection, match=None):
    '''
    Updates MongoDB documents.

    Sets ``field`` equal to ``value`` for all documents that
    meet ``match`` criteria.

    Arguments:

        field (str): Field to update.

        value (str): Update value.

        db (Database): A pymongo Database object.

        collection (str): Collection name.

        match (dict): A dictionary containing the match criteria, for example::

            {'seq_id': {'$in': ['a', 'b', 'c']}, 'cdr3_len': {'$gte': 18}}
    '''
    c = db[collection]
    match = match if match is not None else {}
    # check MongoDB version to use appropriate update command
    if db.client.server_info()['version'].startswith('2'):
        c.update(match, {'$set': {field: value}}, multi=True)
    else:
        c.update_many(match, {'$set': {field: value}})


def unset(db, collection, field, match=None):
    '''
    Removes ``field`` from all records in ``collection`` that meet
    ``match`` criteria.

    Arguments:

        field (str): Field to be removed.

        db (Database): A pymongo Database object.

        collection (str): Collection name.

        match (dict): A dictionary containing the match criteria, for example::

            {'seq_id': {'$in': ['a', 'b', 'c']}, 'cdr3_len': {'$gte': 18}}
    '''
    c = db[collection]
    match = match if match is not None else {}
    # check MongoDB version to use appropriate update command
    if db.client.server_info()['version'].startswith('2'):
        c.update(match, {'$unset': {field: ''}}, multi=True)
    else:
        c.update_many(match, {'$unset': {field: ''}})


def mongoimport(json, database,
                ip='localhost', port=27017,
                user=None, password=None,
                delim='_', delim1=None, delim2=None,
                delim_occurance=1, delim1_occurance=1, delim2_occurance=1):
    '''
    Performs mongoimport on one or more json files.

    Args:

        json: Can be one of several things:

            - path to a single JSON file
            - an iterable (list or tuple) of one or more JSON file paths
            - path to a directory containing one or more JSON files

        database (str): Name of the database into which the JSON files
            will be imported

        ip (str): IP address of the MongoDB server. Default is ``localhost``.

        port (int): Port of the MongoDB database. Default is ``27017``.

        user (str): Username for the MongoDB database, if authentication is enabled.
            Default is ``None``, which results in attempting connection without
            authentication.

        password (str): Password for the MongoDB database, if authentication is enabled.
            Default is ``None``, which results in attempting connection without
            authentication.

        delim (str): Delimiter, when generating collection names using a single delimiter.
            Default is ``_``

        delim_occurance (int): Occurance at which to split filename when using a
            single delimiter. Default is ``1``

        delim1 (str): Left delimiter when splitting with two delimiters. Default is None.

        delim1_occurance (int): Occurance of ``delim1`` at which to split filename.
            Default is ``1``

        delim2 (str): Right delimiter when splitting with two delimiters. Default is None.

        delim2_occurance (int): Occurance of ``delim2`` at which to split filename.
            Default is ``1``
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


def index(db, collection, fields, directions=None, desc=False, background=False):
    '''
    Builds a simple (single field) or complex (multiple fields) index
    on a single collection in a MongoDB database.

    Args:

        db (Database): A pymongo Database object.

        collection (str): Collection name.

        fields: Can be one of two things:

            - the name of a single field, as a string
            - an iterable (list/tuple) of one or more field names

        desc (bool): If ``True``, all indexes will be created in descending order.
            Default is ``False``.

        directions (list): For complex indexes for which you'd like to have
            different indexing directions (ascending for some fields, descending
            for others), you can pass a list of pymongo direction objects (
            pymongo.ASCENDING and pymongo.DESCENDING), in the same order as the
            list of fields to be indexed. Must be the same length as the list
            of index fields. Default is ``None``.

        background (bool): If ``True``, the indexing operation will be processed
            in the background. When performing background indexes, the MongoDB
            database will not be locked.
    '''
    import pymongo
    if type(fields) == str:
        fields = [fields, ]
    if directions is None:
        _dir = pymongo.DESCENDING if desc else pymongo.ASCENDING
        directions = [_dir] * len(fields)
    field_tuples = zip(fields, directions)
    coll = db[collection]
    coll.create_index(field_tuples, background=background)


def remove_padding(db, collection, field='padding'):
    '''
    Removes a padding field.

    Args:

        db (Database): A pymongo Database object.

        collection (str): Collection name

        field (str): Name of the padding field. Default is ``padding``
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
