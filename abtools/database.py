#!/usr/local/bin/python
# filename: databases.py

#
# Copyright (c) 2017 Bryan Briney
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


import cPickle as pickle
import os
import sqlite3
import time


class KeyValueDB(object):
    """A SQLite database for storing key/value pairs.

    Args:

        name (str): Name for the database file. Required, except for in-memory databases.

        dir (str): Directory for the database. Required, except for in-memory databases.

        pickle_values (bool): If ``True``, values will be pickled (with ``cPickle``) before
            entry into the database. Necessary if you want to store Python objects as values.
            Default is ``False``.

        in_memory (bool): If ``True``, the KeyValueDB will be stored in memory rather than on disk.
            Default is ``False``.

    """
    def __init__(self, name=None, dir=None, pickle_values=True, in_memory=False):
        super(Database, self).__init__()
        self._name = name
        self._dir = dir
        self._path = None
        self.table_name = 'data'
        self.pickle_values = pickle_values
        self.structure = [('key', 'text'), ('value', 'text')]
        self.initialized = False
        self._connection = None
        self._cursor = None
        self._create_table_cmd = None
        self._insert_cmd = None
        if not self.initialized:
            self.create_table()

    def __getitem__(self, key):
        return self.find_one(key)

    def __setitem__(self, key, value):
        return self.insert_one(key, value)


    @property
    def name(self):
        if self._name is None:
            self._name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
        return self._name

    @name.setter
    def name(self, name):
        self._name = name


    @property
    def dir(self):
        if self._dir is None:
            if self.in_memory:
                self._dir = None
            else:
                self._dir = '/tmp'
        return self._dir

    @dir.setter
    def dir(self, dir):
        self._dir = dir


    @property
    def path(self):
        if self._path is None:
            if self.in_memory:
                self._path = ':memory:'
            else:
                self._path = os.path.join(self.dir, self.name)
        return self._path


    @property
    def connection(self):
        if self._connection is None:
            self._connection = sqlite3.connect(self.path)
        return self._connection

    @connection.setter
    def connection(self, connection):
        self._connection = connection


    @property
    def cursor(self):
        if self._cursor is None:
            self._cursor = self.connection.cursor()
        return self._cursor

    @cursor.setter
    def cursor(self, cursor):
        self._cursor = cursor


    @property
    def create_table_cmd(self):
        if self._create_table_cmd is None:
            field_string = ', '.join([' '.join(s) for s in self.structure])
            self._create_table_cmd = 'CREATE TABLE {} ({})'.format(self.table_name,
                field_string)
        return self._create_table_cmd


    @property
    def insert_cmd(self):
        if self._insert_cmd is None:
            self._insert_cmd = 'INSERT INTO {} VALUES ({})'.format(self.table_name,
                ','.join(['?'] * len(self.structure)))
        return self._insert_cmd


    def commit(self):
        '''
        Commits changes to the database.
        '''
        self.connection.commit()


    def close(self):
        '''
        Closes the database connection, but leaves the KeyValueDB object intact.
        The connection will automatically be re-established if necessary. Primarily designed
        to allow the KeyValueDB object to be passed as an argument via ``multiprocessing``,
        as open SQLite connection objects cannot be pickled.
        '''
        self.connection.close()
        del self._cursor
        del self._connection
        self._cursor = None
        self._connection = None


    def terminate(self):
        '''
        Terminates the database and, for databases stored on disk, removes the database file.
        All data will be lost and reconnection will not be possible.
        '''
        self.close()
        if not self.in_memory:
            os.unlink(self.path)


    def create_table(self):
        self.cursor.execute('DROP TABLE IF EXISTS {}'.format(self.table_name))
        self.cursor.execute(self.create_table_cmd)
        self.initialized = True


    def insert_one(self, key, value):
        '''
        Inserts a single key/value pair

        Args:

            key (str): The key

            value: The value. Can be an arbitrary object if the database was initialised
                with ``pickle_values``, otherwise should be a ``str``, ``int`` or ``float``.
        '''
        if self.pickle_values:
            value = pickle.dumps(value, protocol=0)
        with self.connection as conn:
            conn.execute(self.insert_cmd, (key, value))


    def insert_many(self, data):
        '''
        Inserts multiple key/value pairs.

        Args:

            data: a list/generator of iterables (lists or tuples) with each iterable containing
                a single key/value pair.
        '''
        if self.pickle_values:
            data = ((d[0], pickle.dumps(d[1], protocol=0)) for d in data)
        with self.connection as conn:
            conn.executemany(self.insert_cmd, data)


    def find_one(self, key):
        '''
        Finds one value with a key that matches ``key``.

        Args:

            key (str): a single key

        Returns:

            a single value (unpickled if necessary)
        '''
        self.cursor.execute(
            '''SELECT data.key, data.value
               FROM data
               WHERE data.key LIKE ?''', (key, ))
        if self.pickle_values:
            return pickle.loads(str(self.cursor.fetchone()[1]))
        return self.cursor.fetchone()[1]


    def find(self, keys):
        '''
        Finds all values with a key that match ``keys``. If ``keys`` is a list, values from
        any records with a key that matches any value in ``keys`` will be returned.

        Args:

            keys: a single key (as a ``str``) or iterable (``list`` or ``tuple``)
                containing one or more keys

        Returns:

            list: values, unpickled if necessary
        '''
        if type(keys) in [str, unicode]:
            keys = [keys, ]
        results = []
        for chunk in self.chunker(keys):
            result_chunk = self.cursor.execute(
                '''SELECT data.key, data.value
                   FROM data
                   WHERE data.key IN ({})'''.format(','.join('?' * len(chunk))), chunk)
            results.extend(result_chunk)
        if self.pickle_values:
            return [pickle.loads(str(r[1])) for r in results]
        return [r[1] for r in results]


    def find_all(self):
        '''
        Returns all values in the database.

        Returns:

            list: all values, unpickled if necessary
        '''
        results = self.cursor.execute(
            '''SELECT data.value
               FROM data''')
        if self.pickle_values:
            return [pickle.loads(str(r[0])) for r in results]
        return [r[0] for r in results]


    def delete(self, keys):
        '''
        Deletes all records with a key that match ``keys``. If ``keys`` is a list, any records
        with a key that matches any value in ``keys`` will be deleted.

        Args:

            keys: a single key (as a ``str``) or iterable (``list`` or ``tuple``)
                containing one or more keys
        '''
        if type(keys) in [str, unicode]:
            keys = [keys, ]
        with self.connection as conn:
            conn.executemany(
                '''DELETE
                   FROM data
                   WHERE data.key == ?''', keys)


    def index(self, field='key'):
        '''
        Indexes the database

        Args:

            field (str): The field on which to create the index. Default is 'key'.
        '''
        index_name = field + '_index'
        self.cursor.execute('CREATE INDEX {} ON {} ({})'.format(index_name,
            self.table_name,
            field))


    @staticmethod
    def chunker(l, n=900):
        '''
        Yields successive n-sized chunks from l.
        '''
        for i in xrange(0, len(l), n):
            yield l[i:i + n]
