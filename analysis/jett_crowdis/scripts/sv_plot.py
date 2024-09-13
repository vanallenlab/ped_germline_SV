import matplotlib.pyplot as plt
from matplotlib import patches
import palettable

def plot_genomic_elements(ax, df, window=None, features = ('exon', 'CDS'),
                          feature_span = {'exon': 0.5, 'CDS': 1},
                          feature_colors = None,
                          label='top'):
    '''
    Plots genomic elements on an axis object. Defaults to exons and CDSes,
    but other elements can be specified.
    
    Params:
    -------
    ax: axis object to plot on. The y-axis will span -1 to 1.
    df: dataframe containing genomic element. Must include 'feature', 'start', and 'end'.
    window: the x axis limits. Elements outside this window will not be plotted.
    features: which features to plot. Defaults to 'exon' and 'CDS'.
    feature_span: relative height of the features. If extra features are included, but their height is not
        specified, it defaults to 1 (the full y span).
    feature_colors: color of the features. Defaults to color wheel.
    label: whether features should be labelled with a numeric value. Only labels the exon
        and start codons. 'top' for top, 'bottom' for bottom labeling. Start codon requires strand.
    '''
    
    for ft in features:
        if ft not in feature_span:
            feature_span[ft] = 1
            
    if feature_colors is None:
        feature_colors = {}

    # add genomic features
    exon_count = 0
    for index, row in df.iterrows():
        ft = row["feature"]
        start, end = row[['start', 'end']].values
        
        if ft == 'exon':
            exon_count += 1
        
        # don't show elements outside the window, clip those that overlap it
        if window is not None:
            if end < window[0] or start > window[1]:
                continue
            if end > window[1]:
                end = window[1]
            if start < window[0]:
                start = window[0]
        
        if ft not in features:
            continue
        
        # get feature aesthetics
        low, high = -1 * feature_span[ft], feature_span[ft]
        ft_color = feature_colors.get(ft)
        
        patch = patches.Rectangle(
                (start, low),
                end - start,
                high - low,
                color=ft_color,
                fill=True,
                zorder=3,
            )
        ax.add_patch(patch)
        
        if label and ft in ['exon', 'start_codon']:
            y = 1.15 if label == 'top' else -1.15
            if ft == 'exon':
                s = exon_count
            elif ft == 'start_codon':
                s = '>' if row['strand'] == '+' else '<'
                
            va = 'bottom' if label == 'top' else 'top'
            ax.text(x = (start + end) / 2, y = y, s = s, ha = 'center', va = va)
        
    # some reformatting
    ax.set_ylim([-1, 1])
    ax.axhline(y=0, color="black", linewidth=0.75)
    
    return ax
    
    
    
    