"""
diadem
Copyright 2019 William E. Fondrie

See the README for detailed documentation and examples.
"""
from pkg_resources import get_distribution, DistributionNotFound

__version__ = get_distribution(__name__).version
