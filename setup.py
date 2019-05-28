from setuptools import setup

"""
Description of how to make a python package

https://python-packaging.readthedocs.io/en/latest/command-line-scripts.html

"""


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='PModel',
      version='0.0.1',
      long_description=readme(),
      description='Automate protein molecular dynamics analysis',
      url='https://github.com/zhenglz/PModel',
      author='zhenglz',
      author_email='zhenglz@outlook.com',
      license='GPL-3.0',
      install_requires=[
          'numpy',
          'pandas',
          'openmm',
          'mdtraj',
      ],
      include_package_data=True,
      zip_safe=False,
      python_requires='>=3.5',
      )
