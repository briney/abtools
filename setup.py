import os

try:
    from setuptools import setup
except ImportError:
    sys.exit('ERROR: setuptools is required.\n')

try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements
# try:
#     from pip.req import parse_requirements
# except ImportError:
#     sys.exit('ERROR: pip is required.\n')


if os.environ.get('READTHEDOCS', None):
    # Set empty install_requires to get install to work on readthedocs
    install_requires = []
else:
    req_file = 'requirements.txt'
    try:
        reqs = parse_requirements(req_file, session=False)
    except TypeError:
        reqs = parse_requirements(req_file)
    try:
        install_requires = [str(r.req) for r in reqs]
    except AttributeError:
        install_requires = [str(r.requirement) for r in reqs]


config = {
    'description': 'Utilities for analysis of antibody NGS data',
    'author': 'Bryan Briney',
    'url': 'https://www.github.com/briney/abtools',
    # 'download_url': 'www.github.com/briney/abtools/',
    'author_email': 'briney@scripps.edu',
    'version': '0.2.0',
    'install_requires': install_requires,
    'packages': ['abtools'],
    'scripts': ['bin/abcompare',
                'bin/abcorrect',
                'bin/abfinder',
                'bin/abstats',
                'bin/ssh_tunnel'],
    'name': 'abtools',
    'include_package_data': True
}

setup(**config)
