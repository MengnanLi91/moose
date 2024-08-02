import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Function to read CSV file and select specific columns
def read_and_select_columns(file_path, columns):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file_path)

    # Select the specified columns
    selected_columns_df = df[columns]

    return selected_columns_df

# Function to compute the GCI index
def compute_gci(fine, medium, coarse, fine_size, medium_size, coarse_size):
    # Calculate the apparent order of accuracy (p)
    r21 = medium_size / fine_size
    r32 = coarse_size/ medium_size

    epsilon21 = medium - fine
    epsilon32 = coarse - medium

    p = np.log(abs(epsilon21 / epsilon32)) / np.log(r21)

    exact_solution21 = (r21**p*fine - medium)/(r21**p-1)
    exact_solution32 = (r32**p*medium - coarse)/(r32**p-1)

    epsilon21_a = abs((fine-medium)/fine)
    epsilon21_ext = abs((exact_solution21-fine)/exact_solution21)

    epsilon32_a = abs((medium-coarse)/medium)
    epsilon32_ext = abs((exact_solution32-medium)/exact_solution32)

    # Calculate the GCI for medium and fine grids
    F_safety = 1.25
    gci_fine_medium = F_safety *epsilon21_a/(r21**p-1)
    gci_medium_coarse = F_safety *epsilon32_a/(r32**p-1)

    return p, gci_fine_medium, gci_medium_coarse, epsilon21, epsilon32

# Define the paths to your CSV files
csv_files = {
    'fine': 'refine_2.csv',
    'medium': 'refine_1.csv',
    'coarse': 'refine_0.csv'
}

# Define the columns you want to select from each CSV file
columns_to_select = ['grid_size', 'vg_x', 'vslip_value']

# Read and select columns from each CSV file
fine_df = read_and_select_columns(csv_files['fine'], columns_to_select)
medium_df = read_and_select_columns(csv_files['medium'], columns_to_select)
coarse_df = read_and_select_columns(csv_files['coarse'], columns_to_select)

# Extract the required data from the selected DataFrames
fine_grid_size = fine_df['grid_size'].iloc[0]
medium_grid_size = medium_df['grid_size'].iloc[0]
coarse_grid_size = coarse_df['grid_size'].iloc[0]

for var in ['vg_x', 'vslip_value']:
  fine_solution = fine_df[var].iloc[0]
  medium_solution = medium_df[var].iloc[0]
  coarse_solution = coarse_df[var].iloc[0]


  # Compute GCI
  p, gci_fine_medium, gci_medium_coarse, epsilon21, epsilon32 = compute_gci(fine_solution, medium_solution, coarse_solution,
                                                      fine_grid_size, medium_grid_size, coarse_grid_size)

  print(f"Apparent order of accuracy (p): {p:.2f}")
  print(f"GCI for fine to medium grids: {gci_fine_medium:.2%}")
  print(f"GCI for medium to coarse grids: {gci_medium_coarse:.2%}")

  # Prepare data for plotting
  grid_sizes = [coarse_grid_size, medium_grid_size]
  errors = [epsilon32, epsilon21]  # Error for fine grid is assumed to be 0

  # Plot grid size vs. error
  plt.figure(figsize=(10, 6))
  plt.plot(grid_sizes, errors, marker='o')
  plt.xlabel('Grid Size')
  plt.ylabel('Error')
  plt.title(var)
  plt.grid(True)
  plt.gca().invert_xaxis()  # Invert x-axis to show decreasing grid size from left to right
  plt.show()
