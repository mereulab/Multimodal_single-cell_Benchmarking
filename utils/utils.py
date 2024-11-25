

import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
rcParams['figure.figsize'] = (12,12)


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import matplotlib.patheffects as path_effects


sc.settings.verbosity =0

def adjust_labels(coordinates, label_size, scale_factor,categories, max_iterations=5, displacement=1, neighborhood_radius=1.0):
    """
    Adjust label positions to minimize overlaps by moving labels towards areas of lower density.
    Once a label's position is adjusted and it has no overlaps, its position is fixed.

    :param coordinates: numpy array of shape (n_labels, 2) with initial label coordinates.
    :param label_sizes: numpy array of shape (n_labels, 2) with width and height of each label.
    :param scale_factor: scaling factor to adjust label sizes to UMAP coordinate space.
    :param max_iterations: maximum number of iterations for adjustment.
    :param displacement: step size for moving labels.
    :param neighborhood_radius: radius to consider for neighborhood density.
    :return: adjusted coordinates.
    """
    fixed_labels = [False] * len(coordinates)  # Track whether a label's position is fixed

    for _ in range(max_iterations):
        for i in range(len(coordinates)):
            if fixed_labels[i]:
                continue  # Skip this label if its position is already fixed

            direction = np.zeros(2)
            has_overlap = False
            for j in range(len(coordinates)):
                if not fixed_labels[j] and i != j and check_overlap(coordinates[i], coordinates[j],label_size, scale_factor, pad = label_size +20):
                    # print(categories[i],' ', categories[j], 'OVERLAP' )
                    has_overlap = True
                    direction += coordinates[i] - coordinates[j]

            if has_overlap:
                # Normalize and move label
                if np.linalg.norm(direction) != 0:
                    direction /= np.linalg.norm(direction)
                    less_dense_direction = find_less_dense_direction(coordinates, i, neighborhood_radius)
                    direction += less_dense_direction
                    direction /= np.linalg.norm(direction)
                    coordinates[i] += direction * displacement
            else:
                fixed_labels[i] = True  # Mark this label as fixed

        # Check if all labels are fixed, and if so, break the loop
        if all(fixed_labels):
            break

    return coordinates

def check_overlap(coord1, coord2, size, scale_factor=1.0, pad = 25):
    """
    Check if two labels overlap, accounting for scale differences between UMAP coordinates and label sizes.

    :param coord1: coordinates of the first label.
    :param coord2: coordinates of the second label.
    :param size1: size of the first label.
    :param size2: size of the second label.
    :param scale_factor: scaling factor to adjust label sizes to UMAP coordinate space.
    :return: True if labels overlap, False otherwise.
    """
    distance = np.linalg.norm(coord1 - coord2) - pad/scale_factor
    adjusted_size = size / scale_factor
    return distance < adjusted_size

def find_less_dense_direction(coordinates, index, radius):
    """
    Find a direction towards a less dense area for a label.

    :param coordinates: numpy array of label coordinates.
    :param index: index of the label to move.
    :param radius: radius to consider for neighborhood density.
    :return: a unit vector pointing towards a less dense area.
    """
    center = coordinates[index]
    density_vector = np.zeros(2)
    for i, coord in enumerate(coordinates):
        if i != index and np.linalg.norm(coord - center) < radius:
            density_vector += center - coord

    if np.linalg.norm(density_vector) == 0:
        return np.random.rand(2) - 0.5  # Random direction if density is uniform
    else:
        return density_vector / np.linalg.norm(density_vector)
def label_size_to_font_size(label_size):
    # This is a placeholder function. You need to define how to convert label size to font size
    return label_size / 2 

def shorten_line(start, end, shorten_by):
    direction = end - start
    norm_direction = direction / np.linalg.norm(direction)
    return end - norm_direction * shorten_by

# VARIABLES




def umap_refined(adata,
                umap,
                var,
                 size=10,
                label_size = 50,
                width_in_inches = 18,
                height_in_inches = 18,
                max_iterations=20,
                displacement=0.5):
    x_min, y_min = np.min(adata.obsm[umap], axis=0)
    x_max, y_max = np.max(adata.obsm[umap], axis=0)

    umap_width = x_max - x_min
    umap_height = y_max - y_min


    fig = plt.figure(figsize=(width_in_inches, height_in_inches))
    dpi = fig.dpi  # Get the DPI of the figure

    # Convert figure size to pixels
    plot_width_pixels = width_in_inches * dpi
    plot_height_pixels = height_in_inches * dpi

    # You might need to adjust for margins and paddings
    # These are approximate values and might need tweaking
    left_margin = 0.1 * plot_width_pixels
    right_margin = 0.1 * plot_width_pixels
    top_margin = 0.1 * plot_height_pixels
    bottom_margin = 0.1 * plot_height_pixels

    effective_plot_width = plot_width_pixels - (left_margin + right_margin)
    effective_plot_height = plot_height_pixels - (top_margin + bottom_margin)


    scale_factor_x = effective_plot_width / umap_width
    scale_factor_y = effective_plot_height / umap_height

    # Use the smaller scale factor to maintain aspect ratio
    scale_factor = min(scale_factor_x, scale_factor_y)
    print('Figure size=',(width_in_inches, height_in_inches),'dpi=',dpi, 'Scale factor =', scale_factor)

    ax = plt.subplot()
    sc.pl.embedding(adata, basis = umap, color = [var], frameon = False, size = size, show=False, legend_fontoutline=2, ax=ax )


    label_colors = adata.uns[f'{var}_colors']
    # Get initial label positions (centroids of clusters)
    categories = adata.obs[var].cat.categories  # Replace with your category column
    initial_coords = np.array([np.median(adata.obsm[umap][adata.obs[var] == cat], axis=0) for cat in categories])
    # for i in range(len(categories)):
        # print(categories[i], initial_coords[i])
  # Replace with your estimated sizes

    adjusted_coords = adjust_labels(initial_coords.copy(), label_size, scale_factor = scale_factor,categories = categories, max_iterations=2, displacement=0.5)

    scaled_initial_coords = initial_coords 
    scaled_adjusted_coords = adjusted_coords 

    for i, cat in enumerate(categories):
        start_coord = scaled_initial_coords[i]
        end_coord = scaled_adjusted_coords[i]

        # Calculate new endpoint shortened by half the font size
        font_size = label_size_to_font_size(label_size / 2)
        new_end_coord = shorten_line(start_coord, end_coord, font_size/scale_factor)

        # Draw line
        ax.plot([start_coord[0], new_end_coord[0]], 
                [start_coord[1], new_end_coord[1]], 
                color='black', linestyle='-', linewidth=1)

        # Plot label
        text = ax.text(end_coord[0], end_coord[1], cat, 
                       ha='center', va='center', 
                       fontsize=font_size, color=label_colors[i])

        # Add black outline to text
        text.set_path_effects([path_effects.Stroke(linewidth=4, foreground='white'),
                               path_effects.Stroke(linewidth=2, foreground='black'),
                               path_effects.Normal()])
        




import seaborn as sns
def comparison_heatmap(adata, key1, key2, label_1=None, label_2=None, cmap = 'Reds', annot = True, figsize=(7,7)):
    if label_1==None:
        label_1=key1
    if label_2==None:
        label_2=key2
    expected_df = adata.obs[[key1,key2]].groupby(by=[key2,key1]).size().reset_index(name = 'count')
    counts = np.array(expected_df['count'].tolist())
    df = pd.DataFrame(counts.reshape(((len(adata.obs[key2].cat.categories),len(adata.obs[key1].cat.categories)))), index = expected_df[key2].unique(), columns = expected_df[key1].unique())
    if annot ==True:
        annot_ = df.astype(int)
        sc.settings.set_figure_params(figsize=figsize, color_map='inferno')
    else:
        annot_=None
        sc.settings.set_figure_params(figsize=figsize, color_map='inferno')
    s = sns.heatmap(df/np.sum(df,axis = 0), cbar_kws={'label': '% cell shared between annotations',"shrink": .5}, cmap=cmap, vmax=1, vmin=0, annot = annot_,  fmt='.7g',center=0.5,square=True, linewidths=.5)
    s.set_ylabel(label_2, fontsize=12)
    s.set_xlabel(label_1, fontsize = 12)
    # plt.show()
    return df

def comparison_heatmap_percent(adata, key1, key2, label_1=None, label_2=None, cmap = 'Reds', annot = True, figsize=(7,7)):
    if label_1==None:
        label_1=key1
    if label_2==None:
        label_2=key2
    expected_df = adata.obs[[key1,key2]].groupby(by=[key2,key1]).size().reset_index(name = 'count')
    counts = np.array(expected_df['count'].tolist())
    df = pd.DataFrame(counts.reshape(((len(adata.obs[key2].cat.categories),len(adata.obs[key1].cat.categories)))), index = expected_df[key2].unique(), columns = expected_df[key1].unique())
    if annot ==True:
        annot_ = df/np.sum(df,axis = 0)
        sc.settings.set_figure_params(figsize=figsize, color_map='inferno')
    else:
        annot_=None
        sc.settings.set_figure_params(figsize=figsize, color_map='inferno')
    s = sns.heatmap(df/np.sum(df,axis = 0), cbar_kws={'label': '% cell shared between annotations',"shrink": .5}, cmap=cmap, vmax=1, vmin=0, annot = annot_,  fmt='.1f',center=0.5,square=True, linewidths=.5)
    s.set_ylabel(label_2, fontsize=12)
    s.set_xlabel(label_1, fontsize = 12)
    # plt.show()
    return df