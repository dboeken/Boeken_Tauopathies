from lib2to3.pgen2.pgen import DFAState
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import functools
import matplotlib

from smma.src.utilities import closest_node, find_pairs
from smma.src.visualise import plot_colocalisation

from loguru import logger

logger.info('Import OK')

input_path = 'results/spot_detection/colocalisation/colocalisation_summary.csv'
input_path2 = 'results/spot_detection/colocalisation/colocalisation_spots.csv'
output_folder = 'results/spot_detection/colocalisation/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


font = {'family': 'normal',
        'weight': 'normal',
        'size': 14}
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['figure.dpi'] = 600


############
input_path2 = 'results/spot_detection/colocalisation/colocalisation_spots.csv'

spots = pd.read_csv(f'{input_path2}')
spots.drop([col for col in spots.columns.tolist() 
            if 'Unnamed : ' in col], axis=1, inplace=True)

spots['position'] = spots['fov'].str[:4]

#spots[['capture','sample', 'detect']] = spots['sample'].str.split('_', expand=True)

spots[['disease', 'patient']] = spots['sample'].str.split('-', expand=True)


spots2 = spots[spots['channel']==641]
spots3 = spots2#[spots2[disease']=='AD']


for (detect, channel), df in spots3.groupby(['disease', 'channel']):
    fig,ax = plt. subplots()
    sns.histplot(
    data=df, #[df['smoothed_Length'] > 50],
    x='mean_intensity', 
    hue= 'coloc?',
    # common_norm=False
    # palette=palette,
    )
    plt.xlabel('Intensity')
    plt.title(f' {detect} {channel} AT8-T181 colocalisation ') 
    plt.xlim (400,5000)



sns.set_theme (style="ticks", font_scale=1.6)
fig.axes = plt.subplots(1, 4, figsize=(20, 4))
fig.tight_layout()

#print (np. shape (axes))
axesReorder = [axes[1], axes[2], axes[3], axes[0]]
for x, (day, df) in enumerate(spots3.groupby(['detect', 'day', 'channel'])):
    f = sns.histplot (
        data=df,
        x=' mean_intensity',
        hue= 'coloc?',
        ax=axesReorder[x],
    )
axesReorder[x].set_title(f'{day}')
axesReorder[x].set(xlabel='Intensity', ylabel='Count')
axesReorder[x].set_xlim(400, 8000)
axesReorder[x].set_ylim(00, 1800)
fig.tight_layout()
plt.savefig(f'{output_folder}coloc_brightness.svg')

#########


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
input_path2 = 'results/spot_detection/colocalisation/colocalisation_spots.csv'
spots = pd.read_csv(input_path2)
spots.drop([col for col in spots.columns.tolist() if 'Unnamed : ' in col], axis=1, inplace=True)

spots['position'] = spots['fov'].str[:4]
spots[['disease', 'patient']] = spots['sample'].str.split('-', expand=True)

# Filter for channel 641
spots2 = spots[spots['channel'] == 641]
spots3 = spots2

# Specify the disease order
disease_order = ['PSP', 'Pick', 'CBD', 'AD', 'HC']

# Create subplots: one row, one column per disease
#fig, axes = plt.subplots(1, len(disease_order), figsize=(5 * len(disease_order), 5), sharey=True)
fig, axes = plt.subplots(1, len(disease_order), figsize=(20, 3), sharey=True)


# Define the custom palette for the 'coloc?' column
custom_palette = {0: '#b4bfd0ff', 1: '#6c6e7eff'}

# Loop over the ordered diseases and plot each on its corresponding subplot
for ax, disease in zip(axes, disease_order):
    # Filter the dataframe for the current disease
    df_disease = spots3[spots3['disease'] == disease]
    
    # If there is no data for a disease, print a message on that subplot
    if df_disease.empty:
        ax.text(0.5, 0.5, 'No data available', horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes)
        ax.set_title(f"{disease} - 641")
        ax.set_xlabel('Intensity')
        ax.set_xlim(400, 5000)
        continue

    # Plot the histogram of mean_intensity with a hue for 'coloc?' using the custom palette
    sns.histplot(data=df_disease, x='mean_intensity', hue='coloc?', palette=custom_palette, ax=ax)
    ax.set_xlabel('Intensity')
    ax.set_title(f"{disease} - 641")
    ax.set_xlim(400, 4000)

plt.tight_layout()
plt.savefig(f'{output_folder}brightness_T181.svg')
plt.show()






#####

# Read in summary FOV data
coloc_spots = pd.read_csv(f'{input_path}')
coloc_spots.drop([col for col in coloc_spots.columns.tolist()
                  if 'Unnamed: ' in col], axis=1, inplace=True)

coloc_spots['position'] = coloc_spots['fov'].str[:4]
#coloc_spots['spots_count'] = coloc_spots['spots_count'].fillna(0)
coloc_spots['sampleID'
      ] = coloc_spots['sample']
coloc_spots[['sample', 'patient']
      ] = coloc_spots['sample'].str.split('-', expand=True)


# Read in summary FOV data
coloc_spots2 = pd.read_csv(f'{input_path2}')
coloc_spots2.drop([col for col in coloc_spots2.columns.tolist()
                  if 'Unnamed: ' in col], axis=1, inplace=True)

coloc_spots2['position'] = coloc_spots2['fov'].str[:4]
#coloc_spots['spots_count'] = coloc_spots['spots_count'].fillna(0)
coloc_spots2['sampleID'
      ] = coloc_spots2['sample']
coloc_spots2[['sample', 'patient']
      ] = coloc_spots2['sample'].str.split('-', expand=True)




#calculate mean
meandf = coloc_spots.groupby(
    ['sample', 'position', 'channel', 'detect', 'capture', 'layout', 'sampleID']).median().reset_index()

# sample_dict = {'AD1': 'AD', 'AD2': 'AD', 'AD3': 'AD',
#                'AD4': 'AD', 'AD5': 'AD', 'AD6': 'CRL', 'AD7': 'CRL', 'AD8': 'CRL', 'AD9': 'CRL', 'AD10': 'CRL',}
antibody_dict = {488: 'T181', 641: 'AT8', 'colocalised': 'colocalised'}

# meandf['disease_state'] = meandf['sample'].map(sample_dict)


#only colocalised images
# filtered_meandf = meandf[meandf['detect'] == 'AT8 r T181 b'].copy()
filtered_meandf = meandf

#remove one set of channels for coloc spots
filtered_coloc = filtered_meandf[filtered_meandf['channel'] == 488].copy()
filtered_coloc = filtered_coloc[[ 'channel', 'sample',
                                 'detect', 'coloc_spots', 'sampleID']].rename(columns={'coloc_spots': 'total_spots'})
filtered_coloc['channel'] = 'colocalised'


for_plotting = pd.concat([filtered_meandf[[
                         'channel', 'sample', 'detect', 'total_spots', 'sampleID']], filtered_coloc])
for_plotting.to_csv(f'{output_folder}coloc.csv')
#for_plotting['antibody_state'] = for_plotting['channel'].map(antibody_dict)

palette = {
    'T181': '#17A398',
    488: '#17A398',
    'AT8': '#b0185E',
    641: '#b0185E',
    'colocalised': '#e3b504',
}

palette2 = {
    'T181': 'black',
    'AT8': 'black',
    'colocalised': 'black',
}

sns.set_theme(style="ticks", font_scale=1.6)
fig, axes = plt.subplots(1, 2, figsize=(20, 4))
fig.tight_layout()

for x, (capture, df) in enumerate(for_plotting.groupby('detect')):
    f = sns.barplot(
        data=df,
        x='sampleID',
        y='total_spots',
        hue='channel',
        #palette=palette,
        ax=axes[x],
        capsize=.15,
        errwidth=0.7,
        saturation=0.8,
        alpha=0.9
    )
    axes[x].tick_params(axis='x', rotation=90)
    sns.stripplot(
        data=df,
        x='sampleID',
        y='total_spots',
        hue='channel',
        ax=axes[x],
        #palette=palette2,
        dodge=True,
        s=5,
        alpha=0.8

    )

    axes[x].set_title(capture)
    axes[x].set(xlabel='Disease State', ylabel='# of spots per FOV')
    #plt.legend(f, ['T181', 'AT8', 'colocalised'])

    #handles, labels = axes[x].get_legend_handles_labels()
    axes[0].legend('')
    #axes[2].legend('')
    #axes[3].legend('')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    axes[1].legend(by_label.values(), by_label.keys(),
                   bbox_to_anchor=(5.0, 6.0))

    plt.tight_layout()

    sns.move_legend(
        f, "lower center",
        bbox_to_anchor=(-0.3, -0.45), ncol=3, title=None, frameon=False,
    )
    plt.savefig(f'{output_folder}coloc.svg')


#change error width and capsize, fade out, alpha, change size scatterplots use s =


#plotting


for_plotting_proportion = filtered_meandf[[
    'channel', 'sample', 'detect', 'proportion_coloc', 'chance_proportion_coloc', 'sampleID']].copy()
for_plotting_proportion.to_csv(f'{output_folder}coloc_prportion.csv')
# for_plotting_proportion['antibody_state'] = for_plotting_proportion['channel'].map(
#     antibody_dict)
#sns.set(font_scale=1)
sns.set_theme(style="ticks", font_scale=1.4)

fig, axes = plt.subplots(1, 2, figsize=(14, 14))
fig.tight_layout()
for x, (capture, df) in enumerate(for_plotting_proportion.groupby('detect')):
    f = sns.barplot(
        data=df,
        x='sampleID',
        y='proportion_coloc',
        hue='channel',
        palette=palette,
        ax=axes[x],
        capsize=.15,
        errwidth=0.7,
        saturation=0.8,
        alpha=0.9
    )
    axes[x].tick_params(axis='x', rotation=90)
    # sns.stripplot(
    #     data=df,
    #     x='disease_state',
    #     y='proportion_coloc',
    #     hue='antibody_state',
    #     ax=axes[x],
    #     palette=palette2,
    #     dodge=True,
    #     s=5,
    #     alpha=0.8

    # )
    axes[x].set_title(capture)
    axes[x].set(xlabel='Disease State', ylabel='Fraction of spots coloaclised')
    axes[x].set_ylim(0, 100)

    axes[x].axhline(df['chance_proportion_coloc'].mean(),
                    linestyle='--', color='lightgrey')
    logger.info(f"Chance av.: {df['chance_proportion_coloc'].mean()}")

    # axes[x].legend('')
    # handles, labels = plt.gca().get_legend_handles_labels()
    # by_label = dict(zip(labels, handles))
    axes[x].legend(by_label.values(), by_label.keys(),
                   bbox_to_anchor=(1, -0.8))

    plt.tight_layout()

    # sns.move_legend(
    #     f, "lower center",
    #     bbox_to_anchor=(-0.3, -0.45), ncol=3, title=None, frameon=False,
    # )
plt.savefig(f'{output_folder}proportion_coloc.svg')

###### mean



sns.set_theme(style="ticks", font_scale=1.4)

fig, axes = plt.subplots(1, 2, figsize=(14, 14))
fig.tight_layout()
for x, (capture, df) in enumerate(for_plotting_proportion.groupby('detect')):
    f = sns.barplot(
        data=df,
        x='sample',
        y='proportion_coloc',
        hue='channel',
        palette=palette,
        ax=axes[x],
        capsize=.15,
        errwidth=0.7,
        saturation=0.8,
        alpha=0.9
    )
    axes[x].tick_params(axis='x', rotation=90)
    # sns.stripplot(
    #     data=df,
    #     x='disease_state',
    #     y='proportion_coloc',
    #     hue='antibody_state',
    #     ax=axes[x],
    #     palette=palette2,
    #     dodge=True,
    #     s=5,
    #     alpha=0.8

    # )
    axes[x].set_title(capture)
    axes[x].set(xlabel='Disease State', ylabel='Fraction of spots coloaclised')
    axes[x].set_ylim(0, 100)

    axes[x].axhline(df['chance_proportion_coloc'].mean(),
                    linestyle='--', color='lightgrey')
    logger.info(f"Chance av.: {df['chance_proportion_coloc'].mean()}")

    # axes[x].legend('')
    # handles, labels = plt.gca().get_legend_handles_labels()
    # by_label = dict(zip(labels, handles))
    axes[x].legend(by_label.values(), by_label.keys(),
                   bbox_to_anchor=(1, -0.8))

    plt.tight_layout()
    plt.savefig(f'{output_folder}prop_coloc.svg')





for_plotting_intensity = filtered_meandf[['sampleID', 'channel', 'sample','capture', 'mean_intensity-coloc', 'mean_intensity-noncoloc']].copy()
for_plotting_intensity['antibody_state'] = for_plotting_intensity['channel'].map(antibody_dict)

melted = pd.melt(for_plotting_intensity, id_vars=['sampleID', 'channel', 'antibody_state', 'sample', 'capture'], value_vars=[
                 'mean_intensity-coloc', 'mean_intensity-noncoloc'], value_name='coloc_intensity', var_name='coloc_status')

melted['coloc_status'] = melted['coloc_status'].str.replace(
    'mean_intensity-', '')
melted['key'] = melted['coloc_status'] + \
    '_' + melted['antibody_state'].astype(str)

palette3 = {
    'coloc_T181': '#e3b504',
    'coloc_AT8': '#e3b504',
    'noncoloc_T181': '#17A398',
    'noncoloc_AT8': '#b0185E',
    'colocalised': '#e3b504',
}



sns.set_theme(style="ticks", font_scale=1.4)
fig, axes = plt.subplots(1, 2, figsize=(14, 4))
fig.tight_layout()
for x, (capture, df) in enumerate(melted.groupby(['capture'])):
    f = sns.barplot(
        data=df,
        x='sample',
        y='coloc_intensity',
        hue='key',
        palette=palette3,
        ax=axes[x],
        capsize=.15,
        errwidth=0.7,
        #width = 15,
        saturation=0.8,
        alpha=0.9,
        hue_order=['noncoloc_T181', 'coloc_T181',
                   'noncoloc_AT8', 'coloc_AT8']
    )

    sns.stripplot(
        data=df,
        x='sample',
        y='coloc_intensity',
        hue='key',
        ax=axes[x],
        palette=palette3,
        dodge=True,
        s=5,
        alpha=0.8,
        hue_order=['noncoloc_T181', 'coloc_T181', 'noncoloc_AT8', 'coloc_AT8']

    )
    axes[x].set_title(capture)
    axes[x].set(xlabel='Disease State', ylabel='Mean intensity')
axes[0].legend('')
axes[1].legend()
plt.tight_layout()
sns.move_legend(
    f, "lower center",
    bbox_to_anchor=(-0.4, -0.5), ncol=3, title=None, frameon=False,
)
axes[0].legend('')
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
axes[1].legend(by_label.values(), by_label.keys(), bbox_to_anchor=(2.0, 1.0))

plt.tight_layout()

sns.move_legend(
    f, "lower center",
    bbox_to_anchor=(-0.3, -0.4), ncol=3, title=None, frameon=False,
)


###plot ratio for intensity coloc vs non-coloc

for_plotting_intensity = filtered_meandf[['sampleID', 'channel', 'sample', 'capture', 'mean_intensity-coloc',  'mean_intensity-noncoloc']].copy()
for_plotting_intensity['antibody_state'] = for_plotting_intensity['channel'].map(
    antibody_dict)
for_plotting_intensity['intensity_ratio'] = for_plotting_intensity['mean_intensity-coloc'] / \
    for_plotting_intensity['mean_intensity-noncoloc']
for_plotting_intensity['log2_intensity_ratio'] = np.log2(
    for_plotting_intensity['intensity_ratio'])


palette3 = {
    'T181': '#17A398',
    'AT8': '#b0185E',
    'colocalised': '#e3b504',
}


sns.set_theme(style="ticks", font_scale=1.4)
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
fig.tight_layout()
for x, (capture, df) in enumerate(for_plotting_intensity.groupby(['capture'])):
    f = sns.barplot(
        data=df,
        x='sample',
        y='intensity_ratio',
        hue='antibody_state',
        palette=palette3,
        ax=axes[x],
        capsize=.15,
        errwidth=0.7,
        #width = 15,
        saturation=0.8,
        alpha=0.9,
        order=['AD', 'CBD', 'HC', 'PSP', 'Pick']
    )

    sns.stripplot(
        data=df,
        x='sample',
        y='intensity_ratio',
        hue='antibody_state',
        ax=axes[x],
        palette=palette2,
        dodge=True,
        s=5,
        alpha=0.8,
         order=['AD', 'CBD', 'HC', 'PSP', 'Pick']

    )
    axes[x].set_title(capture)
    axes[x].set(xlabel='Disease State', ylabel='log2 fold change intensity')
axes[0].legend()
axes[1].legend()
plt.tight_layout()
sns.move_legend(
    f, "lower center",
    bbox_to_anchor=(-0.3, -0.4), ncol=3, title=None, frameon=False,
)

axes[0].legend('')
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
axes[1].legend(by_label.values(), by_label.keys(), bbox_to_anchor=(2.0, 1.0))

plt.tight_layout()

sns.move_legend(
    f, "lower center",
    bbox_to_anchor=(-0.3, -0.45), ncol=3, title=None, frameon=False,
)
plt.savefig(f'{output_folder}log2_intensity.svg')

#1 sample t test or 2 way anova


#read in df
#
# do calcualtions
### average FOVs: group by position and colour
### average triplicates - one replicate will be in a different csv
##### do for each detection Ab, wavelength, disease vs CRL
##### do for number, number of coloc
#
# make plots






# Visualise distribution of sizes
    
for (detect, channel), df in coloc_spots2.groupby(['sample', 'channel']):
    fig, ax = plt.subplots()
    sns.ecdfplot(
        data=df,
        x='mean_intensity',
        hue='coloc?',
        # common_norm=False
        #palette=palette,
    )
    ax.set_xlabel('Mean intensity')
    plt.title(f'{detect} {channel}')
    #ax.set_xlim(left=0, right=1000)
    #plt.xlabel('Particle length (nm)')
    #plt.title(f'{detect} {protein}')
    plt.savefig(f'{output_folder}{detect}_brightness_dist.png')
    plt.savefig(f'{output_folder}{detect}_brightness_dist.svg')
    
    
    
for (detect, channel), df in coloc_spots2.groupby(['sample', 'channel']):
    fig, ax = plt.subplots()
    sns.histplot(
        data=df,
        x='mean_intensity',
        hue='coloc?',
        # common_norm=False
        #palette=palette,
        #kde=True,
        
        #element='step',
        #stat="density", common_norm=False,
    )
    ax.set_xlabel('Mean intensity')
    plt.title(f'{detect} {channel}')
    

for (detect, channel), df in coloc_spots2.groupby(['sample', 'channel']):
    fig, ax = plt.subplots()
    sns.kdeplot(
        data=df,
        x='mean_intensity',
        hue='coloc?',
        fill = True,
        # common_norm=False
        #palette=palette,
        #de=True
    )
    plt.title(f'{detect} {channel}')
    ax.set_xlim(left=0, right=7000)
    ax.set_xlabel('Mean intensity')
    plt.savefig(f'{output_folder}{detect}_kde_dist.svg')


#############
    
coloc_spots2_fil = coloc_spots2[coloc_spots2['coloc?']==1].copy().reset_index()

coloc_spots2_fil = coloc_spots2_fil.drop_duplicates(subset=["well_info", "pair_id"])

channel_641 = coloc_spots2_fil[coloc_spots2_fil["channel"] == 641]
channel_488 = coloc_spots2_fil[coloc_spots2_fil["channel"] == 488]

merged_df = pd.merge(
    channel_641,
    channel_488,
    on=["fov", "pair_id"],
    suffixes=("_641", "_488")
)

# Calculate the ratio of mean_intensity and its natural logarithm
merged_df["ratio"] = merged_df["mean_intensity_488"] / merged_df["mean_intensity_641"]
merged_df["ln_ratio"] = merged_df["ratio"].apply(lambda x: pd.np.log(x) if x > 0 else None)

result_df = merged_df[
    [
        "index_641",
        "fov",
        "sample_641",
        "capture_641",
        "detect_641",
        "well_info_641",
        "ratio",
        "ln_ratio",
    ]
].rename(
    columns={
        "index_641": "index",
        "sample_641": "sample",
        "capture_641": "capture",
        "detect_641": "detect",
        "well_info_641": "well_info",
    }
)

result_df['sampleID'
      ] = result_df['sample']
result_df[['sample', 'patient']
      ] = result_df['sample'].str.split('-', expand=True)



sns.set_theme(style="ticks", font_scale=1.4)
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

fig, ax = plt.subplots()
fig.tight_layout()
for (detect), df in result_df.groupby(['sampleID']):
    fig, ax = plt.subplots()
    sns.histplot(
        data=df,
        x='ln_ratio',
        #hue='coloc?',
        #fill = True,
        # common_norm=False
        #palette=palette,
        #de=True
    )
    plt.title(f'{detect}')
    
    
import matplotlib.pyplot as plt
import seaborn as sns
import math

# Determine the grid dimensions
n_rows = 6
n_cols = 5

# Create a figure with 5x6 panels
fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 15))  # Adjust figsize as needed
axes = axes.flatten()  # Flatten the 2D array of axes into a 1D array for easier iteration

# Plot each sample in a separate panel
for i, (detect, df) in enumerate(result_df.groupby('sampleID')):
    sns.histplot(
        data=df,
        x='ln_ratio',
        ax=axes[i]  # Assign each plot to a specific axis
    )
    axes[i].set_title(f'{detect}')  # Set the title for each subplot
    axes[i].set_xlim(-5, 5)

# Hide any unused axes
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

# Adjust layout to prevent overlap
plt.tight_layout()
plt.show()




import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import seaborn as sns

# Define the sum of Gaussians function
def gaussian_mixture(x, *params):
    n = len(params) // 3  # Number of Gaussians
    y = np.zeros_like(x)
    for i in range(n):
        A = params[i * 3]
        mu = params[i * 3 + 1]
        sigma = params[i * 3 + 2]
        y += A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))
    return y

# Fit histograms for each detect
n_rows, n_cols = 6, 5
fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 15))
axes = axes.flatten()

for i, (detect, df) in enumerate(result_df.groupby('sample')):
    # Get histogram data
    hist, bin_edges = np.histogram(df['ln_ratio'].dropna(), bins=30, range=(-4, 4), density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Initial guess for parameters (2 Gaussians as an example)
    initial_guess = [0.5, -1, 0.5, 0.5, 1, 0.5]
    
    try:
        # Fit the data
        popt, _ = curve_fit(gaussian_mixture, bin_centers, hist, p0=initial_guess)
        
        # Plot histogram
        sns.histplot(data=df, x='ln_ratio', bins=30, ax=axes[i], stat="density")
        axes[i].set_xlim(-5, 5)
        axes[i].set_title(f'{detect}')
        
        # Overlay the fit
        x_fit = np.linspace(-4, 4, 500)
        y_fit = gaussian_mixture(x_fit, *popt)
        axes[i].plot(x_fit, y_fit, 'r-')  # No label added
    except RuntimeError:
        axes[i].text(0.5, 0.5, "Fit failed", horizontalalignment='center', transform=axes[i].transAxes)

# Hide unused axes
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

plt.tight_layout()
plt.show()
#####

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.mixture import GaussianMixture

# Fit histograms for each detect
n_rows, n_cols = 6, 5
fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 15))
axes = axes.flatten()

for i, (detect, df) in enumerate(result_df.groupby('sampleID')):
    # Extract the ln_ratio values (dropping NaNs)
    data = df['ln_ratio'].dropna().values.reshape(-1, 1)
    
    # Fit a Gaussian Mixture Model (e.g., 2 components)
    n_components = 3  # You can increase this number for more complex distributions
    gmm = GaussianMixture(n_components=n_components, random_state=42)
    gmm.fit(data)
    
    # Create histogram
    sns.histplot(data=data.flatten(), bins=30, stat="density", ax=axes[i])
    axes[i].set_xlim(-5, 5)
    axes[i].set_title(f'{detect}')
    
    # Generate x values for plotting the fit
    x_fit = np.linspace(-4, 4, 500).reshape(-1, 1)
    y_fit = np.exp(gmm.score_samples(x_fit))  # GMM outputs log-density, so we exponentiate
    
    # Overlay the GMM fit
    axes[i].plot(x_fit, y_fit, color='red')
    
# Hide unused axes
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

plt.tight_layout()
plt.show()




####

import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import seaborn as sns

# Define the number of rows and columns for the grid
n_rows, n_cols = 5, 6

# Create subplots
fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 15))
axes = axes.flatten()

for i, (detect, df) in enumerate(result_df.groupby('sampleID')):
    data = df['ln_ratio'].dropna().values.reshape(-1, 1)

    # Fit Gaussian Mixture Model (start with 2 components)
    gmm = GaussianMixture(n_components=3, random_state=0)
    gmm.fit(data)

    # Extract parameters
    means = gmm.means_.flatten()
    covariances = np.sqrt(gmm.covariances_).flatten()
    weights = gmm.weights_.flatten()

    # Plot histogram
    sns.histplot(data=data.flatten(), bins=30, stat="density", ax=axes[i])
    axes[i].set_xlim(-4, 4)
    axes[i].set_title(f'{detect}')

    # Plot each Gaussian component
    x_fit = np.linspace(-4, 4, 500)
    y_fit = np.zeros_like(x_fit)
    for mean, cov, weight in zip(means, covariances, weights):
        y = weight * (1 / (np.sqrt(2 * np.pi) * cov)) * np.exp(-0.5 * ((x_fit - mean) / cov) ** 2)
        y_fit += y
        axes[i].plot(x_fit, y, linestyle="--", label="Subpopulation")

    # Plot the sum of Gaussians
    axes[i].plot(x_fit, y_fit, color="red")

# Hide unused axes
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

plt.tight_layout()
plt.show()








from sklearn.mixture import GaussianMixture
import pandas as pd

# Initialize a list to store results
results = []

# Loop through each sample in the DataFrame
for detect, df in result_df.groupby('sampleID'):
    # Extract the 'ln_ratio' data, dropping NaN values
    data = df['ln_ratio'].dropna().values.reshape(-1, 1)
    
    # Fit GMMs for 2, 3, and 4 components
    for n_components in range(2, 5):
        # Fit the Gaussian Mixture Model
        gmm = GaussianMixture(n_components=n_components, random_state=0)
        gmm.fit(data)
        
        # Store results for each component
        for i in range(n_components):
            results.append({
                'sampleID': detect,
                'n_components': n_components,
                'BIC': gmm.bic(data),
                'AIC': gmm.aic(data),
                'mean': gmm.means_[i][0],
                'variance': gmm.covariances_[i][0][0],
                'weight': gmm.weights_[i]
            })

# Convert the results list into a DataFrame
results_df = pd.DataFrame(results)

# Sort results by sampleID and number of components for clarity
results_df = results_df.sort_values(by=['sampleID', 'n_components'])

# Save or display the results
results_df.to_csv("gmm_validation_results2.csv", index=False)
print("Results saved to 'gmm_validation_results.csv'.")



#####

import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import seaborn as sns
import pandas as pd

# Load the optimal components data (replace with your actual file path if needed)
optimal_components_df = pd.read_csv('/Users/dorotheaboeken/Documents/GitHub/240530_Tauopathies_AT8_T181/results/spot_detection/colocalisation/Optimal_Number_of_Components_per_Sample__Updated_.csv')

# Define the number of rows and columns for the grid
n_rows, n_cols = 6,5

# Create subplots
fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 15))
axes = axes.flatten()

# Iterate over each sample and its optimal number of components
for i, (detect, df) in enumerate(result_df.groupby('sampleID')):
    # Get the optimal number of components for this sample
    optimal_components = optimal_components_df.loc[
        optimal_components_df['sampleID'] == detect, 'n_components'
    ].values[0]
    
    data = df['ln_ratio'].dropna().values.reshape(-1, 1)

    # Fit Gaussian Mixture Model with the optimal number of components
    gmm = GaussianMixture(n_components=optimal_components, random_state=0)
    gmm.fit(data)

    # Extract parameters
    means = gmm.means_.flatten()
    covariances = np.sqrt(gmm.covariances_).flatten()
    weights = gmm.weights_.flatten()

    # Plot histogram
    sns.histplot(data=data.flatten(), bins=30, stat="density", ax=axes[i])
    axes[i].set_xlim(-5, 5)
    axes[i].set_title(f'{detect}')

    # Plot each Gaussian component
    x_fit = np.linspace(-4, 4, 500)
    y_fit = np.zeros_like(x_fit)
    for mean, cov, weight in zip(means, covariances, weights):
        y = weight * (1 / (np.sqrt(2 * np.pi) * cov)) * np.exp(-0.5 * ((x_fit - mean) / cov) ** 2)
        y_fit += y
        axes[i].plot(x_fit, y, linestyle="--")  # Subpopulation

    # Plot the sum of Gaussians
    axes[i].plot(x_fit, y_fit, color="red")

# Hide unused axes
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

plt.tight_layout()
plt.show()





#### averaged


# Initialize a list to store results
results = []
print("Unique values in well_info for channel_641:", results_df["sample"].unique())
# Loop through each sample in the DataFrame
for detect, df in result_df.groupby('sample'):
    # Extract the 'ln_ratio' data, dropping NaN values
    data = df['ln_ratio'].dropna().values.reshape(-1, 1)
    
    # Fit GMMs for 2, 3, and 4 components
    for n_components in range(1, 10):
        # Fit the Gaussian Mixture Model
        gmm = GaussianMixture(n_components=n_components, random_state=0)
        gmm.fit(data)
        
        # Store results for each component
        for i in range(n_components):
            results.append({
                'sample': detect,
                'n_components': n_components,
                'BIC': gmm.bic(data),
                'AIC': gmm.aic(data),
                'mean': gmm.means_[i][0],
                'variance': gmm.covariances_[i][0][0],
                'weight': gmm.weights_[i]
            })

# Convert the results list into a DataFrame
results_df = pd.DataFrame(results)

# Sort results by sampleID and number of components for clarity
results_df = results_df.sort_values(by=['sample', 'n_components'])

# Save or display the results
results_df.to_csv("gmm_validation_resultsaveraged.csv", index=False)
print("Results saved to 'gmm_validation_results.csv'.")