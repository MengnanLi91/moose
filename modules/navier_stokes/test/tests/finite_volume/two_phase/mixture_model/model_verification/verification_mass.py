# The mixture model verification
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import os.path

def read_data(file_name):

  my_path = os.path.abspath(os.path.dirname(__file__))
  path = os.path.join(my_path, file_name)

  MASS_data = pd.read_csv(path)

  return MASS_data

def analytical_solution_void(data):
  analytical_mass_data = pd.DataFrame(data['time'])
  analytical_mass_data['void_fraction'] = 0.05 * (np.sin(2 * np.pi * analytical_mass_data['time'] - np.pi / 2) + 1)

  return analytical_mass_data

def plot_void(data, label, marker='o', linestyle='none'):
  plt.plot(data['time'], data['void_fraction'], label=label,linestyle=linestyle, marker=marker)

def plot_norm(x, y, order, label, marker='x', linestyle='none'):

  xx = np.linspace(np.min(x), np.max(x), 30)
  if order == '1':
    a = y[1]/x[1]
    plt.loglog(xx, a*xx, label='first order', linestyle=':')
  if order == '2':
    a = y[1]/(x[1]*x[1])
    plt.loglog(xx, a*xx*xx, label='second order', linestyle=':')
  plt.loglog(x,y, label=label,linestyle=linestyle, marker=marker)

def L2_norm(data, analytical):
  difference = data - analytical
  l2_norm = np.linalg.norm(difference)
  l_infinite = np.max(np.abs(difference))
  return l2_norm, l_infinite



def main():
  file1 = "./dt=0.0125.csv"
  file2 = "./dt=0.025.csv"
  file3 = "./dt=0.1.csv"

  data1 = read_data(file1)
  data2 = read_data(file2)
  data3 = read_data(file3)
  ana_data = analytical_solution_void(data2)

  plot_void(data2, 'MOOSE')
  plot_void(ana_data, 'Analytical', linestyle='-', marker='none')

 # plt.title('Mass Transient')
  plt.xlabel('time(s)')
  plt.ylabel('Void fracton[-]')
  plt.legend(loc='best')
  plt.savefig('results/MassTransient.png', dpi=200)
  plt.show()

  dt = [0.05, 0.025, 0.0125]
  L2_1, Lini_1 = L2_norm(data1['void_fraction'], analytical_solution_void(data1)['void_fraction'])
  L2_2, Lini_2= L2_norm(data2['void_fraction'], ana_data['void_fraction'])
  L2_3, Lini_3 = L2_norm(data3['void_fraction'], analytical_solution_void(data3)['void_fraction'])

  L2_norm_void = []
  L2_norm_void.append(L2_1)
  L2_norm_void.append(L2_2)
  L2_norm_void.append(L2_3)

  Linfini_norm_void = []
  Linfini_norm_void.append(Lini_1)
  Linfini_norm_void.append(Lini_2)
  Linfini_norm_void.append(Lini_3)

  plot_norm(dt, L2_norm_void, order='2', label='$l_2$', marker='x', linestyle='none')
 # plot_norm(dt, Linfini_norm_void, order='2', label='$l_{\infty}$', marker='x', linestyle='none')

  #plt.title('Mass Transient Convergence')
  plt.xlabel('Timestep size(s)')
  plt.ylabel('$l_2$ error of void [-]')
  plt.ylim((pow(10,-15),pow(10,-13)))
  plt.legend(loc='best')
  plt.savefig('results/MassTransientConvergence.png', dpi=200)
  plt.show()

if __name__ == "__main__":
    main()
