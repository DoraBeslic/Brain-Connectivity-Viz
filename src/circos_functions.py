import numpy as np
import pandas as pd
from pycirclize import Circos
import matplotlib.pyplot as plt


def calc_weighted_deg_bound(conn_matrix):
    '''
    Calculates greatest weighted connectivity degree of a connectivity matrix. 
    - Note: Does not preserve sign of greatest weighted connectivity degree.

    Parameters:
    - conn_matrix (np.array): Effective/directed or functional connectivity matrix.
    
    Returns:
    - bound (float): Greatest weighted connectivity degree.
    '''

    # Assume that a symmetrical matrix is a functional connectivity matrix, otherwise it is directed
    if np.array_equal(conn_matrix, conn_matrix.T):
        deg = np.abs(np.sum(conn_matrix, axis=0))
        bound = deg.max()
    else:
        in_deg = np.abs(np.sum(conn_matrix, axis=0))
        out_deg = np.abs(np.sum(conn_matrix, axis=1))
        bound = max(in_deg.max(), out_deg.max())    
    return bound


def calc_min_max_conn_strength(conn_matrix):
    '''
    Calculates the minimum and maximum connectivity strengths of a connectivity matrix.
    - Note: Does not preserve signs of connectivity strengths.
    
    Parameters:
    - conn_matrix (np.array): Effective/directed or functional connectivity matrix.
    
    Returns:
    - conn_strength_min, conn_strength_max (tuple of two float values): Minimum and maximum connectivity strengths.
    '''
    conn_matrix_copy = conn_matrix.copy()
    conn_matrix_copy[conn_matrix_copy == 0] = np.nan
    conn_strength_min = np.nanmin(np.abs(conn_matrix_copy))
    conn_strength_max = np.nanmax(np.abs(conn_matrix_copy))
    return conn_strength_min, conn_strength_max


def plot_degree_tracks(conn_matrix, df_roi_list, deg_bound_lower, deg_bound_upper):
    '''
    Creates a Circos plot labelled with brain regions of interest (ROIs) and adds heatmap track(s) for weighted connectivity degree of each region.
    - For functional connectivty matrices: one track for weighted degree.
    - For directed/effective connectivity matrices: first track for weighted in-degree, second track for weighted out-degree.
        
    Parameters:
    - conn_matrix (np.array): Effective/directed or functional connectivity matrix.
    - df_roi_list (pd.DataFrame): ROI information including atlas, network, label, hemisphere, and X, Y, Z coordinates.
    - deg_bound_lower (float): Lowest weighted degree value.
    - deg_bound_upper (float): Highest weighted degree value. 
    
    Returns:
    - circos (Circos): Circos plot with labelled ROIs and heatmap track(s) of weighted degree.
    '''

    # Determine if connectivity matrix is directed
    # Assume that a symmetrical matrix is a functional connectivity matrix, otherwise it is directed
    directed = False if np.array_equal(conn_matrix, conn_matrix.T) else True
    
    # Define color mapping for different regions
    name2color = {
        'tian_subcortex-right': 'limegreen', 'tian_subcortex-left': 'limegreen', 'brainstem-right': 'deepskyblue',
        'brainstem-left': 'deepskyblue', 'brainstem-medial': 'deepskyblue', 'cerebellum-medial': 'mediumpurple',
        'Vis-right': '#781286', 'SomMot-right': '#4682B4', 'DorsAttn-right': '#00760E', 'SalVentAttn-right': '#C43AFA',
        'Limbic-right': '#DCF8A4', 'Cont-right': '#E69422', 'Default-right': '#CD3E4E', 'Vis-left': '#781286',
        'SomMot-left': '#4682B4', 'DorsAttn-left': '#00760E', 'SalVentAttn-left': '#C43AFA', 'Limbic-left': '#DCF8A4',
        'Cont-left': '#E69422', 'Default-left': '#CD3E4E'
    }

    # Initialize sectors dictionary
    sectors = {}

    # Create sectors for right hemisphere schaefer ROIs
    for network in np.flip(df_roi_list[df_roi_list['atlas'] == 'schaefer']['network'].unique()):
        network_r_labels = df_roi_list[(df_roi_list['network'] == network) & (df_roi_list['hemisphere'] == 'right')]['label'].tolist()
        sectors[f'{network}-right'] = len(network_r_labels)

    # Create sectors for right hemisphere non-schaefer ROIs
    for atlas in df_roi_list['atlas'].unique():
        if atlas != 'schaefer':
            atlas_r_labels = df_roi_list[(df_roi_list['atlas'] == atlas) & (df_roi_list['hemisphere'] == 'right')]['label'].tolist()
            sectors[f'{atlas}-right'] = len(atlas_r_labels)

    # Create sector for medial ROIs
    for atlas in df_roi_list['atlas'].unique():
        if atlas != 'schaefer':
            atlas_m_labels = df_roi_list[(df_roi_list['atlas'] == atlas) & (df_roi_list['hemisphere'] == 'medial')]['label'].tolist()
            sectors[f'{atlas}-medial'] = len(atlas_m_labels)

    # Create sectors for left hemisphere non-schaefer ROIs
    for atlas in np.flip(df_roi_list['atlas'].unique()):
        if atlas != 'schaefer':
            atlas_l_labels = df_roi_list[(df_roi_list['atlas'] == atlas) & (df_roi_list['hemisphere'] == 'left')]['label'].tolist()
            sectors[f'{atlas}-left'] = len(atlas_l_labels)

    # Create sectors for left hemisphere schaefer ROIs
    for network in df_roi_list[df_roi_list['atlas'] == 'schaefer']['network'].unique():
        network_l_labels = df_roi_list[(df_roi_list['network'] == network) & (df_roi_list['hemisphere'] == 'left')]['label'].tolist()
        sectors[f'{network}-left'] = len(network_l_labels)

    # Remove sectors with zero labels
    sectors = {key: value for key, value in sectors.items() if value != 0}
    
    # Create Circos plot
    circos = Circos(sectors, space=3)

    # Create tracks and add heatmaps for each sector
    for sector in circos.sectors:
        track1 = sector.add_track((95, 100))
        track1.axis(fc=name2color[sector.name])
        pos_list = list(range(0, int(track1.size)))
        sub, hemisphere = sector.name.split('-')

        if sub in df_roi_list['atlas'].unique():
            labels = df_roi_list[(df_roi_list['atlas'] == sub) & (df_roi_list['hemisphere'] == hemisphere)]['label'].tolist()
        else:
            labels = df_roi_list[(df_roi_list['network'] == sub) & (df_roi_list['hemisphere'] == hemisphere)]['label'].tolist()
        
        if hemisphere == 'left':
            labels = np.flip(labels)
            
        track1.xticks(
            pos_list,
            labels,
            label_size=6,
            label_orientation='vertical',
            text_kws={'fontweight': 'demibold'}
        )
        
        idx1 = df_roi_list[df_roi_list['label'].isin(labels)].index
        idx2 = df_roi_list.index
        vmin, vmax = deg_bound_lower, deg_bound_upper

        if directed:
            in_deg = np.sum(conn_matrix[np.ix_(idx2, idx1)], axis=0)
            out_deg = np.sum(conn_matrix[np.ix_(idx1, idx2)], axis=1)
            track2, track3 = sector.add_track((90, 95)), sector.add_track((85, 90))
            track2.heatmap(in_deg, vmin=vmin, vmax=vmax, rect_kws={'ec': 'black', 'lw': 0.5})
            track3.heatmap(out_deg, vmin=vmin, vmax=vmax, rect_kws={'ec': 'black', 'lw': 0.5})
        else:
            deg = np.sum(conn_matrix[np.ix_(idx2, idx1)], axis=0)
            track2 = sector.add_track((86, 93))
            track2.heatmap(deg, vmin=vmin, vmax=vmax, rect_kws={'ec': 'black', 'lw': 0.5})
    return circos  


def plot_links(circos, conn_matrix, df_roi_list, conn_strength_min, conn_strength_max):
    '''
    Plots connectivity links between ROIs on a Circos plot.
    - Connectivity strength represented by width and transparency of links. 
    
    Parameters:
    - circos (Circos): Circos plot with labelled ROIs. 
    - conn_matrix (np.array): Effective/directed or functional connectivity matrix.
    - df_roi_list (pd.DataFrame): ROI information including atlas, network, label, hemisphere, and X, Y, Z coordinates.
    - conn_strength_min (float): Minimum connectivity strength.
    - conn_strength_max (float): Maximum connectivity strength. 
    - color (str or matplotlib.colors): color of links.
    
    Returns:
    - circos (Circos): Circos plot with links between ROIs. 
    '''
    
    # Determine if connectivity matrix is directed
    # Assume that a symmetrical matrix is a functional connectivity matrix, otherwise it is directed
    directed = False if np.array_equal(conn_matrix, conn_matrix.T) else True
    
    # For every pair of ROIs with a non-zero connectivity value, plot a link between them
    for i in range(conn_matrix.shape[0]):
        for j in range(conn_matrix.shape[1]):
            if conn_matrix[i, j] != 0:
                
                # Get sector and position in sector of first ROI
                if i in df_roi_list[df_roi_list['atlas'] == 'schaefer'].index.tolist():
                    sector1_net = df_roi_list.iloc[i]['network']
                    sector1_hem = df_roi_list.iloc[i]['hemisphere']
                    sector1 = f'{sector1_net}-{sector1_hem}'
                    df_roi_list_sector1 = df_roi_list[(df_roi_list['network'] == sector1_net) & (df_roi_list['hemisphere'] == sector1_hem)].reset_index()
                    pos1 = df_roi_list_sector1[df_roi_list_sector1['index'] == i].index[0]
                else:
                    sector1_atlas = df_roi_list.iloc[i]['atlas']
                    sector1_hem = df_roi_list.iloc[i]['hemisphere']
                    sector1 = f'{sector1_atlas}-{sector1_hem}'
                    df_roi_list_sector1 = df_roi_list[(df_roi_list['atlas'] == sector1_atlas) & (df_roi_list['hemisphere'] == sector1_hem)].reset_index()
                    pos1 = df_roi_list_sector1[df_roi_list_sector1['index'] == i].index[0]
                
                # Get sector and position in sector of second ROI
                if j in df_roi_list[df_roi_list['atlas'] == 'schaefer'].index.tolist():
                    sector2_net = df_roi_list.iloc[j]['network']
                    sector2_hem = df_roi_list.iloc[j]['hemisphere']
                    sector2 = f'{sector2_net}-{sector2_hem}'
                    df_roi_list_sector2 = df_roi_list[(df_roi_list['network'] == sector2_net) & (df_roi_list['hemisphere'] == sector2_hem)].reset_index()
                    pos2 = df_roi_list_sector2[df_roi_list_sector2['index'] == j].index[0]
                else:
                    sector2_atlas = df_roi_list.iloc[j]['atlas']
                    sector2_hem = df_roi_list.iloc[j]['hemisphere']
                    sector2 = f'{sector2_atlas}-{sector2_hem}'
                    df_roi_list_sector2 = df_roi_list[(df_roi_list['atlas'] == sector2_atlas) & (df_roi_list['hemisphere'] == sector2_hem)].reset_index()
                    pos2 = df_roi_list_sector2[df_roi_list_sector2['index'] == j].index[0]
                
                # Calculate normalized and scaled connectivity strength 
                norm_strength = (np.abs(conn_matrix[i, j]) - conn_strength_min) / (conn_strength_max - conn_strength_min)
                scaled_strength = norm_strength * 0.8 + 0.1

                # Set color of link based on connectivity sign
                color = 'red' if conn_matrix[i, j] > 0 else 'blue'
                
                # Plot link - set transparency and line width of link to scaled connectivity strength
                circos.link_line(
                    (sector1, pos1), 
                    (sector2, pos2), 
                    direction=1 if directed else 0,
                    color=color,
                    alpha=scaled_strength,
                    lw=scaled_strength
                )    
    return circos


def create_legend_handles(conn_strength_min, conn_strength_max, color):
    '''
    Creates legend handles for connectivity strengths represented by links on a Circos plot.
    
    Parameters:
    - conn_strength_min (float): Minimum connectivity strength.
    - conn_strength_max (float): Maximum connectivity strength. 
    - color (str or matplotlib.colors): Color of handles (should be the same color as links).
    
    Returns:
    - handles (list of matplotlib.lines.Line2D): List of legend handles.
    '''
    legend_values = np.linspace(0.1, 0.9, 4)
    handles = []
    for alpha_value in legend_values:
        weight = (alpha_value - 0.1) / 0.8 * (conn_strength_max - conn_strength_min) + conn_strength_min
        handle = plt.Line2D([0], [0], color=color, alpha=alpha_value, linewidth=1, label=f'{weight:.2f}')
        handles.append(handle)
    return handles