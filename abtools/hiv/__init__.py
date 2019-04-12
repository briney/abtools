import os

CATNAP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CATNAP_data')

from .bnabs import get_bnab, get_bnabs, get_bnab_dict
