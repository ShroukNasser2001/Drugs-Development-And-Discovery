import pandas as pd
df = pd.read_csv('bioactivity_data_preprocessed.csv')
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
def lipinski(smiles, verbose=False):

    moldata= []
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem) 
        moldata.append(mol)
       
    baseData= np.arange(1,1)
    i=0  
    for mol in moldata:        
       
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
           
        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])   
    
        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1      
    
    columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]   
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)
    
    return descriptors
df_lipinski = lipinski(df.canonical_smiles)
print("Calculate the Molecular weight(MW),Octanol-water partition coefficient(Log P),Hydrogen bond donors and acceptors:")
print("----------------------------------------------------------------------------------------------------------")
print(df_lipinski)
print("***********************************************************")
print("***********************************************************")
print("View the table with molecule_chembl_id, canonical_smiles, bioactivity_class, standard_value:")
print("----------------------------------------------------------------------------------------")
print(df)
print("***********************************************************")
print("***********************************************************")
print("View the table in which all of the above: ")
print("---------------------------------------")
df_combined = pd.concat([df,df_lipinski], axis=1)
print(df_combined)
print("***********************************************************")
print("***********************************************************")
import numpy as np

def pIC50(input):
    pIC50 = []

    for i in input['standard_value_norm']:
        molar = i*(10**-9) # Converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)
        
    return x
df_combined.standard_value.describe()
def norm_value(input):
    norm = []

    for i in input['standard_value']:
        if i > 100000000:
          i = 100000000
        norm.append(i)

    input['standard_value_norm'] = norm
    x = input.drop('standard_value', 1)
        
    return x
df_norm = norm_value(df_combined)
print("View all data with standard_value_norm: ")
print("-------------------------------------")
print(df_norm)
print("***********************************************************")
print("***********************************************************")
df_norm.standard_value_norm.describe()
df_final = pIC50(df_norm)
print("Add the PIC50 and delete the standard_value_norm:")
print("-------------------------------------------------")
print(df_final)
print("**********************************************************")
print("**********************************************************")
df_2class = df_final[df_final.bioactivity_class != 'intermediate']
print("view all data but remove intermediate in bioactivity_class:")
print("-----------------------------------------------------------")
print(df_2class)
print("**********************************************************")
print("**********************************************************")
import seaborn as sns
sns.set(style='ticks')
import matplotlib.pyplot as plt
print("1-compare between frequency and bio activity class (active and inactive):")
print("2-compare between(Log P) and (MW) in bio activity class")
print("3-compare between(PIC50) and bio activity class :")
print("4-Discription of (PIC50)")
print("5-(Compare between (MW) and (bioactivity)) and (Discription of (MW))")
print("6-(Compare between (LogP) and (bioactivity)) and (Discription of (LogP))")
print("7-(Compare between (NumHDonors) and (bioactivity)) and (Discription of (NumHDonors))")
print("8-(Compare between (NumHAcceptors) and (bioactivity)) and (Discription of (NumHAcceptors))")
print("******************************************************************************************")
print("*******************************************************************************************")
plt.figure(figsize=(5.5, 5.5))

sns.countplot(x='bioactivity_class', data=df_2class, edgecolor='black')

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('Frequency', fontsize=14, fontweight='bold')

plt.savefig('plot_bioactivity_class.pdf')
plt.figure(figsize=(5.5, 5.5))

sns.scatterplot(x='MW', y='LogP', data=df_2class, hue='bioactivity_class', size='pIC50', edgecolor='black', alpha=0.7)
plt.xlabel('MW', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
plt.savefig('plot_MW_vs_LogP.pdf')
plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'pIC50', data = df_2class)
plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('pIC50 value', fontsize=14, fontweight='bold')
plt.savefig('plot_ic50.pdf')
def mannwhitney(descriptor, verbose=False):
  from numpy.random import seed
  from numpy.random import randn
  from scipy.stats import mannwhitneyu
  seed(1)
  selection = [descriptor, 'bioactivity_class']
  df = df_2class[selection]
  active = df[df.bioactivity_class == 'active']
  active = active[descriptor]

  selection = [descriptor, 'bioactivity_class']
  df = df_2class[selection]
  inactive = df[df.bioactivity_class == 'inactive']
  inactive = inactive[descriptor]

  stat, p = mannwhitneyu(active, inactive)

  alpha = 0.05
  if p > alpha:
    interpretation = 'Same distribution (fail to reject H0)'
  else:
    interpretation = 'Different distribution (reject H0)'
  
  results = pd.DataFrame({'Descriptor':descriptor,
                          'Statistics':stat,
                          'p':p,
                          'alpha':alpha,
                          'Interpretation':interpretation}, index=[0])
  filename = 'mannwhitneyu_' + descriptor + '.csv'
  results.to_csv(filename)
  return results
print("Discription of (PIC50)")
print(mannwhitney('pIC50'))
print("****************************************************************")
plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'MW', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('MW', fontsize=14, fontweight='bold')

plt.savefig('plot_MW.pdf')
print("Discription of (MW)")
print(mannwhitney('MW'))
print("****************************************************************")
plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'LogP', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')

plt.savefig('plot_LogP.pdf')
print("Discription of (LOgP)")
print(mannwhitney('LogP'))
print("****************************************************************")
plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'NumHDonors', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('NumHDonors', fontsize=14, fontweight='bold')

plt.savefig('plot_NumHDonors.pdf')
print("Discription of (NumHDonors)")
print(mannwhitney('NumHDonors'))
print("****************************************************************")
plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'NumHAcceptors', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('NumHAcceptors', fontsize=14, fontweight='bold')

plt.savefig('plot_NumHAcceptors.pdf')
print("Discription of (NumHAcceptors)")
print(mannwhitney('NumHAcceptors'))
print("****************************************************************")
print("All drugs are effective and acceptable in (WHO)")
print("****************************************************************")
print("****************************************************************")


# In[ ]:





# In[ ]:




