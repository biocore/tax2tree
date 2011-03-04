#!/usr/bin/env python
from __future__ import division
from distutils.core import setup
import re

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The tax2tree project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"
 
long_description = """Taxonomy decoration tools for phylogenetic trees

http://tax2tree.sourceforge.net

"""
try:
    import cogent
except ImportError:
    print "PyCogent not installed but required. (Is it installed? Is it in the current users $PYTHONPATH or site-packages?) See http://pycogent.sourceforge.net."
    exit(1)
    
pycogent_version = tuple([int(v) \
        for v in re.split("[^\d]", cogent.__version__) if v.isdigit()])
        
if pycogent_version < (1,4,1):
    print "PyCogent >= 1.6.0 required, but %s is installed." % cogent.__version__
    exit(1)
    
setup(name='tax2tree',
      version=__version__,
      description='Taxonomy to tree decoration tools',
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url='http://tax2tree.sourceforge.net',
      packages=['t2t'],
      scripts=['scripts/t2t_workflow'],
      long_description=long_description,
)
