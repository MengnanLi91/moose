# Model validation

import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_data_and_relative_error(experimental_df, numerical_df, labels, dpi):
    # Ensure that the dataframes have the same index and columns
    if not experimental_df.index.equals(numerical_df.index) or not experimental_df.columns.equals(numerical_df.columns):
        raise ValueError("Experimental and numerical dataframes must have the same index and columns")

    # Calculate relative error
    relative_error_df = (experimental_df - numerical_df).abs() / experimental_df.abs()

    # Set global font size
    plt.rcParams.update({'font.size': 14})

    for label in labels:
        # Plot experimental vs numerical data
        plt.figure(figsize=(12, 6))

        plt.plot(experimental_df[f'void_fraction_{label}'], experimental_df[f'vg_{label}'], '-o', label=label, markersize=12, linewidth=3, markerfacecolor='none')
        plt.plot(numerical_df[f'void_fraction_{label}'], numerical_df[f'vg_{label}'], '--s', label='MOOSE', markersize=12, linewidth=3, markerfacecolor='none')

        plt.xlabel(r'Void Fraction[-]')
        plt.ylabel(r'Gas velocity($\mathrm{m/s}$)')
        plt.ylim([0,1])
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join('results/', f'gas_velocity_{label}.png'), dpi=dpi)


    for label in labels:
        # Plot relative error
        plt.figure(figsize=(12, 6))

        plt.plot(experimental_df[f'void_fraction_{label}'], relative_error_df[f'vg_{label}'],'-o')

        plt.xlabel('Void Fraction')
        plt.ylabel(f'Relative Error for {label}')
        plt.grid(True)
        plt.savefig(os.path.join('results/', f'relative_error_{label}.png'), dpi=dpi)

# Bhagwat2016
experimental_data = {'void_fraction_Bhagwat2016': [0.455, 0.475, 0.55, 0.6], 'vg_Bhagwat2016': [0.21, 0.58, 0.42, 0.65]}
numerical_data = {'void_fraction_Bhagwat2016': [0.455, 0.475, 0.55, 0.6], 'vg_Bhagwat2016': [0.2095, 0.5804, 0.4176, 0.6480]}
experimental_df = pd.DataFrame(experimental_data)
numerical_df = pd.DataFrame(numerical_data)

# Hibiki2001
Hibiki2001_data = pd.DataFrame({'void_fraction_Hibiki2001': [0.049, 0.092, 0.192], 'vg_Hibiki2001': [0.561, 0.642, 0.672]})
numerical_data2 =pd.DataFrame( {'void_fraction_Hibiki2001': [0.049, 0.092, 0.192], 'vg_Hibiki2001': [0.605, 0.646, 0.635]})
experimental_df=pd.concat([experimental_df,Hibiki2001_data])
numerical_df=pd.concat([numerical_df,numerical_data2])

plot_data_and_relative_error(experimental_df, numerical_df, {'Bhagwat2016','Hibiki2001'}, 400)
