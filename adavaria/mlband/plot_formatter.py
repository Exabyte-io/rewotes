# plot_formatting.py
import matplotlib.pyplot as plt
import matplotlib
# import seaborn as sns
import pandas as pd

# Defining all the parameters for the plots
font_size = 9
golden_ratio = (5 ** .5 - 1) / 2
mm_to_inch = 0.0393701
one_column = 90 * mm_to_inch
two_column = 190 * mm_to_inch

figsize = {
    'one_column': (one_column, one_column * golden_ratio),
    'two_column': (two_column, two_column * golden_ratio),
    'one_column_square': (one_column, one_column),
    'two_column_square': (two_column, two_column),
    'one_column_wide': (one_column, one_column * golden_ratio * 2),
    'two_column_wide': (two_column, two_column * golden_ratio * 2),
}
figsize['default'] = figsize['one_column']

# Assigning the parameters to the matplotlib and seaborn
params = {
    "legend.fontsize": font_size,
    "xtick.labelsize": font_size,
    "ytick.labelsize": font_size,
    'font.size': font_size,
    'axes.labelsize': font_size,
    'axes.titlesize': font_size + 1,
    'grid.color': 'silver',
    'grid.linestyle': '--',
    'grid.linewidth': 1,
    'grid.alpha': 0.45,
    'axes.grid': False,
    'xtick.bottom': True,
    'ytick.left': True,
}
matplotlib.rcParams.update(params)

plt.rcParams['xtick.bottom'] = True  # Show x-axis ticks at the bottom
plt.rcParams['ytick.left'] = True  # Show y-axis ticks on the left

matplotlib.rc('font', family='sans-serif')
matplotlib.rc('font', serif='Helvetica')
matplotlib.rc('text', usetex='false')
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.titlesize'] = font_size + 1

# # Update the parameters for tick labels format to be comma separated
# plt.rcParams['axes.formatter.use_locale'] = True
# plt.rcParams['axes.formatter.use_mathtext'] = True

plt.rcParams['figure.figsize'] = figsize['default']
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['savefig.bbox'] = 'tight'
