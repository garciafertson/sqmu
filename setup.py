from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

exec(open('sqmu/version.py').read()) # loads __version__

setup(name='sqm exploratory utilities',
      version=__version__,
      author='fernando garcia-guevara',
      description='sqm is a set of simple tools to analyse the results from squeezemeta',
      long_description=readme,
      license='GPL3+',
      keywords="",
      packages=find_packages(exclude='docs'),
      install_requires=('keggrest >= 0.1.1 '),
      #setup_requires=['nose>=1.0'],
      #test_suite='nose.collector',
      url='https://github.com/garciafertson/sqmu',
      scripts=['bin/sqmu.py'],
      #data_files=[
      #    ('share', ['share/18S.hmm']),
      #],
)
