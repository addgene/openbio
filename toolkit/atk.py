import click
import logging

from count_spacers import count_spacers
from recombination import recombination
from serotypes import serotypes
from utils import log_to_stdout

L = logging.getLogger(__name__)

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
def atk():
    pass

atk.add_command(serotypes)
atk.add_command(recombination)
atk.add_command(count_spacers)


if __name__ == "__main__":
    log_to_stdout(logging.INFO)
    try:
        atk()
    except Exception as e:
        L.exception(e)
        raise SystemExit('\nCommand terminated with an error: ' + str(e))
