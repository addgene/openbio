import click
import logging
import sys

from serotypes import serotypes
from utils import log_to_stdout

L = logging.getLogger(__name__)

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
def atk():
    pass

atk.add_command(serotypes)


if __name__ == "__main__":
    log_to_stdout(logging.INFO)
    try:
        atk()
    except ValueError as e:
        L.error('\nError: ' + e.message)
        sys.exit(1)