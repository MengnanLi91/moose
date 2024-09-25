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

def analytical_solution_vel(data):
  analytical_data = pd.DataFrame(data['time'])
  analytical_data['outlet_vel'] = 1.992-0.9921*np.cos(np.pi*analytical_data['time'])

  return analytical_data

def plot_vel(data, label, marker='o', linestyle='none'):
  plt.plot(data['time'], data['outlet_vel'], label=label,linestyle=linestyle, marker=marker)

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
  file1 = "./moment_dt=0.00625.csv"
  file2 = "./moment_dt=0.0125.csv"
  file3 = "./moment_dt=0.05.csv"

  data1 = read_data(file1)
  data2 = read_data(file2)
  data3 = read_data(file3)
  ana_data = analytical_solution_vel(data1)

  plot_vel(data2, 'MOOSE')
  plot_vel(ana_data, 'Analytical', linestyle='-', marker='none')

  #plt.title('Momentum Transient')
  plt.xlabel('time(s)')
  plt.ylabel('velocity($m/s$)')
  plt.legend(loc='best')
  plt.savefig('results/MomentumTransient.png', dpi=200)
  plt.show()

  dt = [0.05, 0.025, 0.0125]
  L2_1, Lini_1 = L2_norm(data1['outlet_vel'], analytical_solution_vel(data1)['outlet_vel'])
  L2_2, Lini_2= L2_norm(data2['outlet_vel'], analytical_solution_vel(data2)['outlet_vel'])
  L2_3, Lini_3 = L2_norm(data3['outlet_vel'], analytical_solution_vel(data3)['outlet_vel'])

  L2_norm_vel = []
  L2_norm_vel.append(L2_1)
  L2_norm_vel.append(L2_2)
  L2_norm_vel.append(L2_3)

  Linfini_norm_vel = []
  Linfini_norm_vel.append(Lini_1)
  Linfini_norm_vel.append(Lini_2)
  Linfini_norm_vel.append(Lini_3)

  plot_norm(dt, L2_norm_vel, order='1', label='$l_2$', marker='x', linestyle='none')
 # plot_norm(dt, Linfini_norm_vel, order='1', label='$l_{\infty}$', marker='x', linestyle='none')

  #plt.title('Momentum Transient Convergence')
  plt.xlabel('Timestep size(s)')
  plt.ylabel('$l_2$ error of velocity [$m/s$]')
  plt.legend(loc='best')
  plt.ylim((pow(10,-2),pow(10,-0.5)))
  plt.savefig('results/MomentumTransientConvergence.png', dpi=200)
  plt.show()

if __name__ == "__main__":
    main()
