import os
from setuptools import find_packages
from setuptools import setup
import re


VERSIONFILE=os.path.join('traitarm', '_version.py')
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else: 
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


setup(name='traitarm',
        version = verstr,
        description='Model-T - predictive genotype-phenotype models',
        url = 'http://github.com/aweimann/traitar-model',
        author='Aaron Weimann',
        author_email='weimann@hhu.de',
        license='GNU General Public License, version 3 (GPL-3.0)',
        packages= ['traitarm'],
        include_package_data = True,
        scripts = ['bin/build_edge_matrix_likelihood', 'bin/discretize_likelihood_recon', 'bin/gainLoss.VR01.266.dRep', 'bin/learn', 'bin/merge_gain_loss', 'bin/prune_ncbi_tree', 'bin/summary2gainLoss_input', 'bin/traitarm', 'bin/summary2gainLoss_input', 'bin/write_gainLoss_config'],
        zip_safe=False,
        install_requires = ["scikit-learn > 0.18.1", "traitar >= 1.1.2", "ete2", "dendropy", "seaborn"])
