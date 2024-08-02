# Model validation

import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_data_and_relative_error(experimental_df, numerical_df, label, dpi):
    # Ensure that the dataframes have the same index and columns
    if not experimental_df.index.equals(numerical_df.index) or not experimental_df.columns.equals(numerical_df.columns):
        raise ValueError("Experimental and numerical dataframes must have the same index and columns")

    # Calculate relative error
    relative_error_df = (experimental_df - numerical_df).abs() / experimental_df.abs()

    # Set global font size
    plt.rcParams.update({'font.size': 14})

    # Plot experimental vs numerical data
    plt.figure(figsize=(12, 6))

    plt.plot(experimental_df['void_fraction'], experimental_df['vg'], '-o', label=label)
    plt.plot(numerical_df['void_fraction'], numerical_df['vg'], '--x', label='MOOSE')

    plt.xlabel('Void Fraction')
    plt.ylabel('Gas velocity(m/s)')
    plt.ylim([0,1])
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join('results/', f'gas_velocity_{label}.png'), dpi=dpi)

    # Plot relative error
    plt.figure(figsize=(12, 6))


    plt.plot(experimental_df['void_fraction'], relative_error_df['vg'],'-o')

    plt.xlabel('Void Fraction')
    plt.ylabel('Relative Error')
    plt.ylim([0,0.01])
    plt.grid(True)
    plt.savefig(os.path.join('results/', f'relative_error_{label}.png'), dpi=dpi)

# Bhagwat2016
experimental_data = {'void_fraction': [0.455, 0.475, 0.55, 0.6], 'vg': [0.21, 0.58, 0.42, 0.65]}
numerical_data = {'void_fraction': [0.455, 0.475, 0.55, 0.6], 'vg': [0.2095, 0.5804, 0.4176, 0.6480]}
experimental_df = pd.DataFrame(experimental_data)
numerical_df = pd.DataFrame(numerical_data)

plot_data_and_relative_error(experimental_df, numerical_df, 'Bhagwat2016', 400)

