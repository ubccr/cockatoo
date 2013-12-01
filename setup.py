try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'cockatoo',
    'author': 'Andrew E. Bruno',
    'url': '',
    'download_url': '',
    'author_email': 'aebruno2@buffalo.edu',
    'version': '0.0.1',
    'install_requires': ['nose'],
    'packages': ['cockatoo'],
    'scripts': [],
    'name': 'cockatoo'
}

setup(**config)
