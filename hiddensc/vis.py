import os
from typing import Optional, List

import IPython
import IPython.display as ipython_display
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def peek_df(df: pd.DataFrame, n_peek: int = 2) -> None:
    print(f'DF shape: {df.shape}')
    print(f'Columns: {list(df.columns)}')
    ipython_display.display(df.head(n_peek))

def visual_settings() -> None:
    """"Change matplotlib settings."""
    sns.set_style("white")
    sns.set_style('ticks')
    sns.set_context("paper", font_scale=2.25)
    sns.set_palette(sns.color_palette('bright'))

    params = {'savefig.dpi': 100,
              'lines.linewidth': 3,
              'axes.linewidth': 2.5,
              'savefig.dpi': 300,
              'xtick.major.width': 2.5,
              'ytick.major.width': 2.5,
              'xtick.minor.width': 1,
              'ytick.minor.width': 1,
              'font.weight': 'medium',
              'figure.figsize': (12, 8)
              }

    mpl.rcParams.update(params)
    IPython.display.set_matplotlib_formats('retina')

    pd.set_option("display.precision", 3)
    pd.set_option('display.float_format', lambda x: '%.3f' % x)

def save_figure(fname: str, exts: Optional[List[str]] = None) -> None:
    """Save figure in multiple formats."""
    if '.' in fname:
        raise ValueError('Expected no . in fname={fname}')
    exts = exts or ['svg', 'png']
    for ext in exts:
        fname = f'{fname}.{ext}'
        plt.savefig(fname, dpi=300, transparent=True, bbox_inches='tight')

def save_fig_results(fname: str) -> None:
    """An opinionated file structure for figures."""
    save_figure(fname, exts=['png'])
    dirname, basename = os.path.dirname(fname), os.path.basename(fname)
    for ext in ['pdf', 'svg']:
        save_figure(os.path.join(dirname, basename), exts=[ext])
