import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from loguru import logger
logger.info('Import OK')

input_folder = 'results/super-res/measurements/'
output_folder = 'results/super-res/summary/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


font = {'family': 'normal',
        'weight': 'normal',
        'size': 40}
#matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['figure.dpi'] = 400



# Read in summary measurements
properties = pd.read_csv(f'{input_folder}properties_compiled.csv')
properties = properties[properties['smoothed_length'] > 30].copy()

for df in [properties]:
    df.drop([col for col in df.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)



# Visualise particle count per FOV
measure_count = properties.groupby(['well_info', 'sample', 'slide_position', 'detect', 'disease_state']).count()['area'].reset_index().rename(columns={'area': 'count'})

palette = {
    'AD': '#DE4968',
    'HC': '#FE9F6D',
    'CBD': '#8C2981',
    'Pick': '#3B0F70',
    'BSA': 'lightgrey',
    'PSP': '#000004',
}


# ================= Mean parameters =================

# Measure mean length 
measure_mean_length = properties.groupby(['well_info', 'sample', 'slide_position', 'detect', 'disease_state']).mean()['smoothed_length'].reset_index()
measure_mean_length.to_csv(f'{output_folder}measure_mean_length.csv')

# Measure mean area 
measure_mean_area = properties[properties['area'] > 0].groupby(['well_info', 'sample', 'slide_position', 'detect', 'disease_state']).mean()['area'].reset_index()
measure_mean_area.to_csv(f'{output_folder}area.csv')

# Measure mean eccentricity 
measure_mean_eccentricity = properties.groupby(['well_info', 'sample', 'slide_position', 'detect', 'disease_state']).mean().reset_index()


# ================= Cumulative plots =================
diseases = properties["disease_state"].dropna().unique()
plt.figure(figsize=(8, 6))

for disease in diseases:
    disease_data = properties[properties["disease_state"] == disease]
    patients = disease_data["patient"].dropna().unique()
    
    x_min = disease_data["smoothed_length"].min()
    x_max = disease_data["smoothed_length"].max()
    x_vals = np.linspace(x_min, x_max, 200)
    
    # Collect the ECDF (fraction ≤ x) for each patient 
    y_matrix = []
    for pat in patients:
        pat_data = disease_data[disease_data["patient"] == pat]["smoothed_length"].dropna()
        sorted_vals = np.sort(pat_data)
        n = len(sorted_vals)

        def fraction_leq(x):
            idx = np.searchsorted(sorted_vals, x, side='right')
            return idx / n
        
        y_ecdf = [fraction_leq(xi) for xi in x_vals]
        y_matrix.append(y_ecdf)
    
    if len(y_matrix) == 0:
        continue
    
    y_array = np.array(y_matrix)
    
    y_mean = y_array.mean(axis=0)
    y_std = y_array.std(axis=0)
    
    # Plot the mean ECDF
    plt.plot(x_vals, y_mean, label=f"{disease}")
    plt.fill_between(x_vals, y_mean - y_std, y_mean + y_std, alpha=0.2)

plt.xlabel("smoothed_length")
plt.ylabel("ECDF")
plt.title("Cumulative length")
plt.legend()
plt.show()



# ================= Localisation density =================
properties['density'] = properties['#locs']/ properties['area']
properties_long_aggregates = properties[properties['smoothed_length'] > 300].groupby(['well_info', 'sample', 'slide_position', 'detect', 'disease_state']).mean()['density'].reset_index()

for (detect), df in properties_long_aggregates.groupby(['detect']):

    fig, ax = plt.subplots()
    sns.boxplot(
        data=df.groupby(['sample', 'disease_state']).mean().reset_index(),
        x='disease_state',
        y='density',
        #marker="$\circ$", ec="face",
        palette =palette,
        #color=palette[protein],
        order=['PSP', 'Pick', 'CBD', 'AD', 'HC'],
        # order=[]
    )
    plt.title('Loc density of aggregates > 200 nm')
    plt.ylabel('#locs/px')
    plt.xlabel('Disease')
    ax.set_ylim(bottom=0, top=13)
  

properties_short_aggregates = properties[properties['smoothed_length'] < 100].groupby(['well_info', 'sample', 'slide_position', 'detect', 'disease_state']).mean()['density'].reset_index()

for (detect), df in properties_short_aggregates.groupby(['detect']):

    fig, ax = plt.subplots()
    sns.boxplot(
        data=df.groupby(['sample', 'disease_state']).mean().reset_index(),
        x='disease_state',
        y='density',
        #marker="$\circ$", ec="face",
        palette =palette,
        #color=palette[protein],
        order=['PSP', 'Pick', 'CBD', 'AD', 'HC'],
        # order=[]
    )
    plt.title('Loc density of aggregates < 200 nm')
    plt.ylabel('#locs/px')
    plt.xlabel('Disease')
    ax.set_ylim(bottom=0, top=13)
    
#ratio loc density
properties_merged = pd.merge(properties_long_aggregates, properties_short_aggregates, on = ['well_info', 'sample', 'slide_position', 'detect', 'disease_state'])#properties_long[properties_long['density']].copy()
properties_merged['density_ratio']= properties_merged['density_x']/properties_merged['density_y']


# ===================================================================================
# ================= Proportion calculation long/fibrillar aggregates =================
# ===================================================================================

# =================Setting parameters=================

thresholds = {
    'branch-distance_max': 200,
    'branch-distance_min': 50,
    'eccentricity_max': 0.9,
    'eccentricity_min': 0.6,
    '#locs_max': 30,
    '#locs_min': 10,
    '#locs_density_max': 0.7,
    '#locs_density_min': 0.3
}


# =================Processing functions=================
def thresholding(data, parameter):
    parameter_cat = parameter + '_cat'
    data[parameter_cat] = ['high' if val > thresholds[parameter]
                           else ('low' if val < 100 else 'medium')for val, detect in data[[parameter, 'detect']].values]


def proportion_calc(data, parameter_cat):

    proportion = (data.groupby(['capture', 'sample', 'slide_position', 'detect', 'disease_state', parameter_cat]).count(
    )['label'] / data.groupby(['capture', 'sample', 'slide_position', 'detect', 'disease_state']).count()['label']).reset_index()
    proportion['label'] = proportion['label'] * 100
    proportion = pd.pivot(
        proportion,
        index=['capture', 'sample', 'slide_position',
               'detect', 'disease_state'],
        columns=parameter_cat,
        values='label'
    ).fillna(0).reset_index()

    proportion_plotting = proportion.groupby(
        ['capture', 'sample', 'detect', 'disease_state']).mean().reset_index()

    return proportion_plotting

# =================Proportion calculation length=================
for_plotting = properties
for_plotting['branch-distance'] = for_plotting['smoothed_length']
parameters = ['branch-distance']
for_plotting['label'] = for_plotting.reset_index().index
#for_plotting.index.names = ['label']

for parameter in parameters:
    #for_plotting = thresholding(for_plotting, parameter)
    parameter_cat = parameter + '_cat'
    for_plotting[parameter_cat] = ['high' if val > thresholds[parameter + '_max']
                                       else ('low' if val < thresholds[parameter + '_min'] else 'medium')for val in for_plotting[parameter].values]

proportion_branch_distance = proportion_calc(for_plotting, 'branch-distance_cat')

# =================Proportion calculation eccentricity=================
for_plotting_ecc = properties
parameters = ['eccentricity']
for_plotting_ecc['label'] = for_plotting_ecc.reset_index().index
#for_plotting.index.names = ['label']

for parameter in parameters:
    #for_plotting = thresholding(for_plotting, parameter)
    parameter_cat = parameter + '_cat'
    for_plotting_ecc[parameter_cat] = ['high' if val > thresholds[parameter + '_max']
                                       else ('low' if val < thresholds[parameter + '_min'] else 'medium')for val in for_plotting_ecc[parameter].values]

proportion_ecc = proportion_calc(for_plotting_ecc, 'eccentricity_cat')


proportion_ecc.to_csv(f'{output_folder}proportion_ecc.csv')
proportion_branch_distance.to_csv(f'{output_folder}proportion_long.csv')

