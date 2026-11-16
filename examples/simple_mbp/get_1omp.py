import chilife

def main():
  
  protein = chilife.fetch('1omp', save=True)
  
  site = 66
  
  spin_label = chilife.SpinLabel('R1C', site=site, use_H=True, protein=protein)
  
  spin_labeled_protein = chilife.mutate(protein,spin_label)
  
  chilife.save = chilife.save(f'1omp_{site}R1.pdb', spin_labeled_protein,
                              conect = True)

main()  


