import logging
import re
import sys

import loguru
from loguru import logger

fmt = ("<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
       "<level>{level: <8}</level> | "
       "<level>{message}</level>")
def filter_color(data: dict) -> bool:
    data['message'] = re.sub(r'\033[^m]*m', '', data['message'])
    return True

# Hack loguru's FileDateFormatter
from loguru._file_sink import FileDateFormatter
def new_date_format(self, spec):
    if not spec:
        spec = "%Y-%m-%d"  # We don't care about hours and minutes and seconds when we rotate every day
    return self.datetime.__format__(spec)
FileDateFormatter.__format__ = new_date_format
logger.remove(0)  # Remove default stderr handler
# logger.add(sys.stdout, level='INFO', format=fmt, colorize=True)
logger.add(sys.stderr, level='DEBUG', format=fmt, colorize=True)
logger.add('logs/verdant_guardian.log', level='DEBUG', format=fmt, colorize=False, rotation='00:00', compression='tar.gz', filter=filter_color)

# disable other libs' logging
logging.getLogger('urllib3.connectionpool').setLevel(logging.INFO)
logging.getLogger('chardet.charsetprober').setLevel(logging.INFO)
logging.getLogger('chardet.universaldetector').setLevel(logging.INFO)
