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
  
  ax.plot(2*data['tau'],data['signal_2'],color='tab:green')

  ax.set_xlabel(r'$2\tau$ (µs)')
  ax.set_ylabel('relative echo amplitude')
  plt.savefig(f'{clue_path}/fig.png')

def load_clue_data(file_path):

  df = pd.read_csv(file_path)
  out = {}
  for k,v in df.items():
    if type(v[0])==type('str'):
      signal = np.array(
        (v.str.replace('i','j')).apply(lambda z: complex(z))
      )
      out[k] = np.real(signal)
    else:
      out[k] = np.real(v)


  return out
main()  
