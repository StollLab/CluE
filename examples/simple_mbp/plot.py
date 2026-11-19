from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import os
import re
import sys

def main():
  for d in sys.argv:
    if os.path.isdir(d):
      plot(d)

def plot(clue_path):
  file_path = f'{clue_path}/signal.csv'
  assert os.path.isfile(file_path)

  data = load_clue_data(file_path)

  fig,ax = plt.subplots(1,1,figsize=(6,6))

  ax.plot(2*data['tau'],np.real(data['signal_2']),color='tab:green')

  ax.set_xlabel(r'$2\tau$ (µs)')
  ax.set_ylabel('relative echo amplitude')
  plt.savefig(f'{clue_path}/fig.png')

def load_clue_data(file_path):
  # A signal X is saved as the columns re_X and im_X; out['X'] is the
  # complex signal. Other columns, such as tau and tau2, are kept as they are.
  df = pd.read_csv(file_path)
  out = {}
  for k,v in df.items():
    if k.startswith('re_'):
      name = k[len('re_'):]
      out[name] = np.array(v) + 1j*np.array(df['im_' + name])
    elif not k.startswith('im_'):
      out[k] = np.array(v)

  return out
main()
