#!/usr/bin/env python
# filename: ssh_tunnel.py


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


# Borrowed heavily from Paramiko:
# https://github.com/paramiko/paramiko/blob/master/demos/forward.py


import getpass
import os
import socket
import select
try:
    import SocketServer
except ImportError:
    import socketserver as SocketServer

import argparse
import sys

import paramiko

SSH_PORT = 22
DEFAULT_PORT = 27017


class ForwardServer(socketserver.ThreadingTCPServer):
    daemon_threads = True
    allow_reuse_address = True


class Handler(socketserver.BaseRequestHandler):
    def handle(self):
        try:
            chan = self.ssh_transport.open_channel(
                'direct-tcpip',
                (self.chain_host, self.chain_port),
                self.request.getpeername())
        except Exception as e:
            print('Incoming request to {}:{} failed: {}'.format(
                self.chain_host,
                self.chain_port,
                repr(e)))
            return
        if chan is None:
            print('Incoming request to {}:{} was rejected by the SSH server.'.format(
                self.chain_host,
                self.chain_port))
            return

        print('Connected!  Tunnel open {} -> {} -> {}'.format(
            self.request.getpeername(),
            chan.getpeername(),
            (self.chain_host, self.chain_port)))
        while True:
            r, w, x = select.select([self.request, chan], [], [])
            if self.request in r:
                data = self.request.recv(1024)
                if len(data) == 0:
                    break
                chan.send(data)
            if chan in r:
                data = chan.recv(1024)
                if len(data) == 0:
                    break
                self.request.send(data)
        peername = self.request.getpeername()
        chan.close()
        self.request.close()
        print('Tunnel closed from {}'.format(peername,))


def forward_tunnel(local_port, remote_host, remote_port, transport):
    # this is convoluted, but necessary to configure things for the Handler
    # object.  (SocketServer doesn't give Handlers any way to access the outer
    # server.)
    class SubHander (Handler):
        chain_host = remote_host
        chain_port = remote_port
        ssh_transport = transport
    ForwardServer(('', local_port), SubHander).serve_forever()


def get_host_port(spec, default_port):
    "parse 'hostname:22' into a host and port, with the port optional"
    args = (spec.split(':', 1) + [default_port])[:2]
    args[1] = int(args[1])
    return args[0], args[1]


HELP = """\
Set up a forward tunnel across an SSH server, using paramiko. A local port
(given with -p) is forwarded across an SSH session to an address:port from
the SSH server. This is similar to the openssh -L option.
"""


def parse_arguments():
    parser = argparse.ArgumentParser(usage='usage: ssh_tunnel [options] <ssh-server>[:<server-port>]',
                            version='ssh_tunnel', description=HELP)
    parser.add_argument('ssh_server', nargs=1,
                        help='SSH server, as <ssh-server>[:<ssh-port>]')
    parser.add_argument('-p', '--local-port', action='store', type=int, dest='port',
                        default=DEFAULT_PORT,
                        help='local port to forward (default: {})'.format(DEFAULT_PORT))
    parser.add_argument('-u', '--user', action='store', type=str, dest='user',
                        default=getpass.getuser(),
                        help='username for SSH authentication (default: {})'.format(getpass.getuser()))
    parser.add_argument('-K', '--key', action='store', type=str, dest='keyfile',
                        default=None,
                        help='private key file to use for SSH authentication')
    parser.add_argument('--no-key', action='store_false', dest='look_for_keys', default=True,
                        help="don't look for or use a private key file")
    parser.add_argument('-P', '--password', dest='readpass', default=False, action='store_true',
                        help='Use a password for SSH.')
    parser.add_argument('-r', '--remote', action='store', required=True, type=str, dest='remote',
                        default=None, metavar='host:port',
                        help='remote host and port to forward to')
    args = parser.parse_args()

    if args.remote is None:
        parser.error('Remote address required (-r).')

    server_host, server_port = get_host_port(args.ssh_server[0], SSH_PORT)
    remote_host, remote_port = get_host_port(args.remote, args.port)
    return args, (server_host, server_port), (remote_host, remote_port)


def main():
    args, server, remote = parse_arguments()
    if args.readpass:
        password = getpass.getpass('Enter SSH password: ')

    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.WarningPolicy())

    print('Connecting to ssh host {}:{} ...'.format(server[0], server[1]))
    try:
        client.connect(server[0], server[1], username=args.user, key_filename=args.keyfile,
                       look_for_keys=args.look_for_keys, password=password)
    except Exception as e:
        print('*** Failed to connect to {}:{}: {}'.format(server[0], server[1], e))
        sys.exit(1)

    print('Now forwarding port {} to {}:{} ...'.format(args.port, remote[0], remote[1]))
    try:
        forward_tunnel(args.port, remote[0], remote[1], client.get_transport())
    except KeyboardInterrupt:
        print('C-c: Port forwarding stopped.')
        sys.exit(0)


if __name__ == '__main__':
    main()
