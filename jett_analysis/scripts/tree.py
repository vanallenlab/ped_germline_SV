import networkx as nx
import pandas as pd
import numpy as np


def classify_relationship(cat_i, cat_j, parent_dict):
    
    '''Given two categories and a dictionary explaining the relationships between
    category components, return 1 if category i is a parent of category j.'''

    cat_i_components = cat_i.split('.')
    cat_j_components = cat_j.split('.')
    
    # dicts are ordered in python 3.7
    subparents = []
    for cat, c_i, c_j in zip(parent_dict.keys(), cat_i_components, cat_j_components):
        if c_i == parent_dict[cat] or c_i == c_j:
            subparents.append(1)
        else:
            subparents.append(0)
    
    return all(subparents)


def transfer_edges(g, u, v):
    '''Given a kept node u and a dropped node v, transfer all of v's edges to u
    if they don't already exist.'''
    
    u_edges = [e[1] for e in g.edges(u)]
    v_edges = [e[1] for e in g.edges(v)]
    for n in v_edges:
        if n not in u_edges:
            g.add_edge(u, n)

    return g

def trim_tree_by_similarities(g, node_df, node_positions, preserved_cats, jaccard, threshold, keep_metric):
    '''Takes a tree created by hierarchy and removes nodes that are similar above some threshold
    via jaccard. Preserved categories are never dropped. Jaccard is the jaccard matrix. 
    Threshold is the threshold for jaccard similarity. Keep metric is used to remove nodes--the node with the higher
    keep_metric will be kept.'''
    
    # get the node levels
    node_pos_df = pd.DataFrame(node_positions, index = ['x', 'y']).T
    unique_y = {pos: i for i, pos in enumerate(sorted(set(node_pos_df['y']), reverse = True))}
    node_pos_df['y_level'] = node_pos_df['y'].replace(unique_y)

    # for each node level, remove nodes that have a jaccard similarity greater than some value
    # we keep the node with the highest metric or if it is in our preserved categories.
    nodes_to_remove = []
    for level in unique_y.values():
        level_nodes = node_pos_df[node_pos_df['y_level'] == level]

        node_similarities = jaccard.loc[level_nodes.index, level_nodes.index]

        for (i, j), sim in np.ndenumerate(node_similarities):
            cat_i = node_similarities.index[i]
            cat_j = node_similarities.index[j]

            # if more similar than a threshold, drop it
            if cat_i != cat_j and sim > threshold:

                # if both are preserved, we drop neither
                if cat_i in preserved_cats and cat_j in preserved_cats:
                    continue
                elif cat_i in preserved_cats and cat_j not in preserved_cats:
                    nodes_to_remove.append(cat_j)
                    g = transfer_edges(g, cat_i, cat_j)
                    continue
                elif cat_j in preserved_cats and cat_i not in preserved_cats:
                    nodes_to_remove.append(cat_i)
                    g = transfer_edges(g, cat_j, cat_i)
                    continue

                # if neither are in the preserved categories, drop the one with the metric
                sorted_by_p = node_df.loc[[cat_i, cat_j]].sort_values(by = keep_metric, ascending = False)
                kept, dropped = sorted_by_p.index
                nodes_to_remove.append(dropped)
                g = transfer_edges(g, kept, dropped)

    nodes_to_remove = set(nodes_to_remove)
    
    # remove these nodes
    g.remove_nodes_from(nodes_to_remove)
    
    # get the node positions of these
    node_pos_df = node_pos_df.loc[list(g.nodes)]
    
    # remake node positions
    node_positions = {node: tuple(row[['x', 'y']].values) for node, row in node_pos_df.iterrows()}
    
    return g, node_positions

def tighten_tree(g, node_df, node_positions, x_threshold, keep_metric):
    '''Tightens a tree by enforcing x distances that are very close to each other
    on different levels to be the same ("straightening" the tree). Also removes nodes
    that end up in the exact same position, keeping the one with the highest metric.'''
    
    # get the node levels
    node_pos_df = pd.DataFrame(node_positions, index = ['x', 'y']).T
    unique_y = {pos: i for i, pos in enumerate(sorted(set(node_pos_df['y']), reverse = True))}
    node_pos_df['y_level'] = node_pos_df['y'].replace(unique_y)
             
    # get the x-values of the nodes and their distances from each other
    x_vals = node_pos_df['x'].values
    difference_matrix = np.abs(x_vals[:, np.newaxis] - x_vals)
    diff_df = pd.DataFrame(difference_matrix, index = node_pos_df.index, columns = node_pos_df.index)

    # we "tighten" the tree from the top down. Given a category, force all "lower" nodes that are within some
    # x distance to be the average of the x positions.
    already_tightened = []
    for level in unique_y.values():
        level_nodes = node_pos_df[node_pos_df['y_level'] == level]
        deeper_nodes = node_pos_df[node_pos_df['y_level'] > level]

        for node in level_nodes.index:
            x_dist_deeper = diff_df.loc[node, deeper_nodes.index]

            # use some threshold
            x_dist_deeper_cats = list(x_dist_deeper[x_dist_deeper < x_threshold].index)

            # remove nodes that have already been tightened
            x_dist_deeper_cats = [c for c in x_dist_deeper_cats if c not in already_tightened]

            # average all of them
            xs = node_pos_df['x'].loc[[node] + x_dist_deeper_cats]
            tightened_xs = xs.mean()
            node_pos_df.loc[[node] + x_dist_deeper_cats, 'x'] = tightened_xs

            already_tightened += x_dist_deeper_cats
            
    # next, we remove nodes that are directly overlapping each other, keeping the one that
    # has the highest metric
    kept_nodes = []
    for level in unique_y.values():
        level_nodes = node_pos_df[node_pos_df['y_level'] == level]
        
        # get the overlapping x positions
        dups = level_nodes.duplicated(subset = 'x', keep = False).values
        
        # keep the nondups
        nondups = level_nodes[~dups]
        kept_nodes += list(nondups.index)
        
        # drop the overlapping nodes
        overlapping = level_nodes[dups]
        for x in set(overlapping['x']):
            x_nodes = overlapping[overlapping['x'] == x]
            metric_values = node_df.loc[x_nodes.index].sort_values(by = keep_metric, ascending = False)
            kept_n = metric_values.index[0]
            kept_nodes.append(kept_n)
            
            # when a node is dropped, the kept node must inherit all those nodes connections
            for dropped_node in metric_values.index[1:]:
                g = transfer_edges(g, kept_n, dropped_node)
            
    # drop nodes from graph
    nodes_to_remove = set(g.nodes) - set(kept_nodes)
    g.remove_nodes_from(nodes_to_remove)
    
    # get the node positions of these
    node_pos_df = node_pos_df.loc[list(g.nodes)]
    
    # remake node positions
    node_positions = {node: tuple(row[['x', 'y']].values) for node, row in node_pos_df.iterrows()}
    
    return g, node_positions