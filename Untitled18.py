#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#This code to indtall the data
import pandas as pd
from chembl_webresource_client.new_client import new_client
target = new_client.target
Search = input("Enter your disease: ")
target_query = target.search(Search)
targets = pd.DataFrame.from_dict(target_query)
print(targets)
id_num = int(input("Enter your ID: "))
selected_target = targets.target_chembl_id[id_num]
print(selected_target)
activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
df = pd.DataFrame.from_dict(res)
result=int(input("Enter the number of result: "))
print(df.head(result))
df.standard_type.unique()
df2 = df[df.standard_value.notna()]
print(df2)
bioactivity_class = []
for i in df2.standard_value:
  if float(i) >= 10000:
    bioactivity_class.append("inactive")
  elif float(i) <= 1000:
    bioactivity_class.append("active")
  else:
    bioactivity_class.append("intermediate")
mol_cid = []
for i in df2.molecule_chembl_id:
  mol_cid.append(i)
canonical_smiles = []
for i in df2.canonical_smiles:
  canonical_smiles.append(i)
standard_value = []
for i in df2.standard_value:
  standard_value.append(i)
data_tuples = list(zip(mol_cid, canonical_smiles, bioactivity_class, standard_value))
df3 = pd.DataFrame( data_tuples,  columns=['molecule_chembl_id', 'canonical_smiles', 'bioactivity_class', 'standard_value'])
print(df3)


# In[ ]:





# In[ ]:





# In[ ]:




