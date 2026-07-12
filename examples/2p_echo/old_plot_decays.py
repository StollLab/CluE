import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import re
import sys
import os

color_palette = [
  [222/255,66/255,91/255,1],
  [5/255,135/255,0, 1],
  [84/255,93/255,187/255,1],
  [255/255,166/255,0/255,1],
  [221/255,119/255,136/255,1],
  [102/255,153/255,204/255,1],
  [144/255,238/255,144/255,1],
  [253/255,230/255,47/255,1],
]
#-------------------------------------------------------------------------------
def main():
  if len(sys.argv) < 2:
    print("Please supply a target directory:")
    return

  folder = sys.argv[1]
  if not os.path.isdir(folder):
    print("Please supply a valid directory:")
    return

  time_file = f"{folder}/tau_axis.csv"
  if not os.path.isfile(time_file):
    print(f"\"{folder}\" does not contain tau_axis.csv")

  sig_file = f"{folder}/signal.csv"
  if not os.path.isfile(time_file):
    print(f"\"{folder}\" does not contain signal.csv")

  time = read_tau_axis(time_file)
  signals = read_signal_file(sig_file)

  fig, (ax) = plt.subplots(1,1,figsize=(8, 6) );
  
  legend = []
  for ii,signal in enumerate(signals):
    
    color = color_palette[ii % 8]
    ax.plot(time,np.real(signal),
        linewidth=1,color=color);
    
    cluster_size = ii + 1
    legend.append( f"{cluster_size}-CCE") 

  ax.set_ylim(0,1)
  ax.set_xlim(0,time[-1])
  ax.legend(legend)
  ax.set_xlabel(r"$2\tau$ (μs)")
  ax.set_ylabel("Re(signal)")
  plt.savefig("fig.png")

#-------------------------------------------------------------------------------
def read_tau_axis(csv_file):
  data = pd.read_csv(csv_file);
  t = data['tau_axis'];
  return np.array(t);
#-------------------------------------------------------------------------------
def read_signal_file(csv_file):
  convert_to_python = {'i':'j', '\+-':'-'};

  data = pd.read_csv(csv_file);
  order = 1
  key = "signal"

  signals = []
  
  while True:
    new_key = f"signal_{order}"
    if new_key in data:
      v = (data[new_key].replace(convert_to_python,regex=True)).apply(
          lambda z: complex(z));
      signals.append(np.array(v))
      key = new_key
      order += 1
    else:
      break;

  return signals
#-------------------------------------------------------------------------------
def stretched_exponential(t,V0,TM,xi):
  return V0*np.exp( -(t/TM)**xi )
#-------------------------------------------------------------------------------

#===============================================================================
if __name__ == "__main__":
  main()
#===============================================================================
