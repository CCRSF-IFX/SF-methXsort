# import ez_setup
# import sys
# ez_setup.use_setuptools()
from setuptools import setup

# from mpld3
# from bwameth
def get_version(path):
    """Get the version info from the mpld3 package without importing it"""
    import ast

    with open(path) as init_file:
        module = ast.parse(init_file.read())

    version = (ast.literal_eval(node.value) for node in ast.walk(module)
               if isinstance(node, ast.Assign)
               and node.targets[0].id == "__version__")
    try:
        return next(version)
    except StopIteration:
        raise ValueError("version could not be located")

install_requires = ['toolshed>=0.4.5', 'xengsort>=2.0.9']

setup(name='methXsort',
      version=get_version("methXsort.py"),
      description="Sort bisulfite-/enzymatic-converted reads from xenograft experiment",
      py_modules=['methXsort'],
      author="Shaojun Xie",
      author_email="xies4@nih.gov",
      license="MIT",
      install_requires=install_requires,
      long_description=open('README.md').read(),
      classifiers=[
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'Programming Language :: Python :: 3'
      ],
      scripts=['methXsort.py']
)