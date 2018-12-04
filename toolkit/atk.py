import click
import logging
import sys

from recombination import recombination
from serotypes import serotypes
from utils import log_to_stdout

L = logging.getLogger(__name__)

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
def atk():
    pass

atk.add_command(serotypes)
atk.add_command(recombination)


if __name__ == "__main__":
    log_to_stdout(logging.INFO)
    try:
        atk()
    except Exception as e:
        L.error('\nError: {}'.format(e))
        sys.exit(1)