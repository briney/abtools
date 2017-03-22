import os

try:
    from setuptools import setup
except ImportError:
    sys.exit('ERROR: setuptools is required.\n')


try:
    from pip.req import parse_requirements
except ImportError:
    sys.exit('ERROR: pip is required.\n')


if os.environ.get('READTHEDOCS', None):
    # Set empty install_requires to get install to work on readthedocs
    install_requires = []
else:
    req_file = 'requirements.txt'
    try:
        reqs = parse_requirements(req_file, session=False)
    except TypeError:
        reqs = parse_requirements(req_file)
    install_requires = [str(r.req) for r in reqs]


config = {
    'description': 'Utilities for analysis of antibody NGS data',
    'author': 'Bryan Briney',
    'url': 'https://www.github.com/briney/abtools',
    # 'download_url': 'www.github.com/briney/abtools/',
    'author_email': 'briney@scripps.edu',
    'version': '0.1.5',
    'install_requires': install_requires,
    'packages': ['abtools'],
    'scripts': ['bin/abcompare',
                'bin/abcorrect',
                'bin/abfinder',
                'bin/abphylogeny',
                'bin/abstats',
                'bin/ssh_tunnel'],
    'name': 'abtools',
    'include_package_data': True
}

setup(**config)
