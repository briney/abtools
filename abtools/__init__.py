import warnings
warnings.filterwarnings('ignore')

import matplotlib as mpl
mpl.use('Agg', warn=False)

# import _compare as compare
from _correct import run as correct
from _finder import run as finder
import _phylogeny as phylogeny
# import _stats as stats
