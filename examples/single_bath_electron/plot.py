from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

PI = np.pi
MUB = 9.2740100783e-24
MU0 = 1.25663706212e-6
ELECTRON_G = 2.00231930436256
PLANCK = 6.62607015e-34
HBAR = (PLANCK/2.0/PI)

SRC = 'CluE-out'
def main():
  time_axis, signal = load_clue_data(SRC)

  fig,ax = plt.subplots(1,1,figsize=(6,6))

  ax.plot(time_axis,np.real(signal))
#  ax.plot(time_axis,np.imag(signal))
#  ax.plot(time_axis,np.abs(signal))
  
  ref_sig = get_expected_signal(time_axis)
  ax.plot(time_axis,ref_sig,linestyle='--')

  rmsd = np.sqrt(np.mean(np.abs(ref_sig - signal)**2))
  ax.set_title(f'RMSD = {rmsd}')
  ax.set_ylim([0,1.2])

#  ax.legend(['Re(v)','Im(v)','|v|'])
  ax.legend(['CluE','analytic'])
  ax.set_xlabel('time (μs)')
  ax.set_ylabel('signal')

  plt.savefig('fig.png')
  plt.close()

def load_clue_data(src: str) -> (np.ndarray,np.ndarray):  
  df = pd.read_csv(f'{src}/time_axis.csv')

  return np.array(df['time_axis']), load_signal(src)

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

def get_expected_signal(t: np.ndarray):
  k = -(MU0/4/PI)*(ELECTRON_G*MUB)**2/HBAR

  theta = 0
  ori_fac = (3*np.cos(theta)**2 -1)
  r = 30e-10
  azz = k*ori_fac/r**3
  omega = -azz/2
  print(omega)
  print(0.5*omega/PI)
  return 0.75 + 0.25*np.cos(omega*t*1e-6)

if __name__ == '__main__':
  main()
