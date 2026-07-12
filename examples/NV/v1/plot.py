from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

SRC = 'CluE-out'
def main():
  time_axis, signal = load_clue_data(SRC)

  fig,ax = plt.subplots(1,1,figsize=(6,6))

#  ax.plot(time_axis,np.real(signal))
#  ax.plot(time_axis,np.imag(signal))
  ax.plot(time_axis,np.abs(signal))

#  ax.legend(['Re(v)','Im(v)','|v|'])
  ax.set_xlabel('time (μs)')
  ax.set_ylabel('signal (μs)')
  ax.set_ylim([0,1.2])

  plt.savefig('fig.png')
  plt.close()

def load_clue_data(src: str) -> (np.ndarray,np.ndarray):  
  df = pd.read_csv(f'{src}/tau_axis.csv')

  return np.array(df['tau1_axis']), load_signal(src)

def load_signal(src: str) -> np.ndarray:  
  df = pd.read_csv(f'{src}/signal.csv')

  key0 = 'signal_'
  key = ''
  for ii in range(1,100):
   new_key = f'signal_{ii}'
   if new_key in df.keys():
     key = new_key
   else:
     break  

  signal = np.array(
      (df[key].str.replace('i','j')).apply(lambda z: complex(z))
  );   

  return signal

if __name__ == '__main__':
  main()
