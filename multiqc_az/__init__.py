from pkg_resources import get_distribution
from multiqc.utils import config

__version__ = get_distribution("multiqc_az").version
config.multiqc_az_version = __version__

