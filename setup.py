from setuptools import setup, find_packages
setup(
    name="wrfcmaq2inmap",
    version="0.1",
    packages=find_packages(),
    scripts = ['bin/wrfcmaq2inmap',],
    package_data = {'ancillary': ['ancillary/*'],
      'scripts': ['scripts/*']},
    python_requires='>3.5',
    setup_requires=['numpy>=1.12','netCDF4>=1.2.9','pandas'],
    install_requires=['numpy>=1.12','netCDF4>=1.2.9','pandas'],
    author_email='beidler.james@epa.gov'
)

