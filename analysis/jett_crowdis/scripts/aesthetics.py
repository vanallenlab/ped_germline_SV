"""Shared settings for controlling figure aesthetics. Borrowed from MX"""
import functools

import palettable
import seaborn as sns
from cycler import cycler
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import sys


custom_rcParams = {
    "figure.figsize": (4, 4),
    "font.family": "Arial",
    "font.size": 12,
    "font.weight": "regular",
    "axes.labelsize": 13,
    "axes.formatter.useoffset": False,
    "axes.formatter.limits": (-4, 4),
    "axes.titlesize": 14,
    "axes.grid": False,
    "legend.edgecolor": "none",
    "legend.fancybox": False,
    "legend.fontsize": 12,
    "legend.frameon": False,
    "legend.framealpha": 0,
    "legend.facecolor": "none",
    "legend.loc": "center left",
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "savefig.facecolor" : "white",
    "savefig.edgecolor": "white",
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.prop_cycle": cycler(
        color=palettable.cartocolors.qualitative.Bold_10.mpl_colors
    ),
}

illustrator_rcParams = {
    "figure.figsize": (8, 3),
    "font.family": "Arial",
    "font.size": 6,
    "font.weight": "regular",
    "axes.labelsize": 6,
    "axes.formatter.useoffset": False,
    "axes.formatter.limits": (-4, 4),
    "axes.titlesize": 6,
    "axes.grid": False,
    "legend.edgecolor": "none",
    "legend.fancybox": False,
    "legend.fontsize": 6,
    "legend.frameon": False,
    "legend.framealpha": 0,
    "legend.facecolor": "none",
    "legend.loc": "center left",
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "savefig.facecolor" : "white",
    "savefig.edgecolor": "white",
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "axes.prop_cycle": cycler(
        color=palettable.cartocolors.qualitative.Bold_10.mpl_colors
    )
}

def activate_paper_rcParams(activate='default'):
    
    if activate == 'default':
        rcParams.update(custom_rcParams)
        sns.set_context("paper", rc=custom_rcParams)
    elif activate == 'illustrator':
        rcParams.update(illustrator_rcParams)
        sns.set_context("paper", rc=illustrator_rcParams)


def strip_axis(ax, x="strip", y="strip", despine=True):
    """Given an axis, strips it to (most/all) axis elements and associated spines

    strip means remove the entire axis
    label means label the axis but remove the spine
    ignore means do not touch the axis

    Params:
    -------
    ax: axis object to strip
    x: operation to perform on the x axis. One of ['strip', 'label', 'ignore']
    y: operation to perform on the y axis. One of ['strip', 'label', 'ignore']
    """

    if despine:
        bottom = x in ["strip", "label"]
        left = y in ["strip", "label"]
        sns.despine(ax=ax, bottom=bottom, left=left)

    if x == "label":
        ax.tick_params(axis="x", length=0)
    elif x == "strip":
        ax.get_xaxis().set_visible(False)

    if y == "label":
        ax.tick_params(axis="y", length=0)
    elif y == "strip":
        ax.get_yaxis().set_visible(False)

    return ax

