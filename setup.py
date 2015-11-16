from setuptools import setup

config = {
	'description': 'AbTools',
	'author': 'Bryan Briney',
	'url': 'www.github.com/briney/abtools/',
	# 'download_url': 'www.github.com/briney/abstar/',
	'author_email': 'briney@scripps.edu',
	'version': '0.1.0',
	'install_requires': ['nose',
						 'biopython',
						 'celery',
						 'ete2',
						 'matplotlib',
						 'numpy',
						 'nwalign',
						 'pandas',
						 'pymongo',
						 'seaborn'],
	'packages': ['abtools'],
	'scripts': ['bin/abtools'],
	'name': 'abstar',
	'include_package_data': True
}

setup(**config)
