import pickle
td_path = r'D:\openclaw\intelligence-augmentation\models\Geneformer\geneformer\token_dictionary_gc104M.pkl'
name_path = r'D:\openclaw\intelligence-augmentation\models\Geneformer\geneformer\gene_name_id_dict_gc104M.pkl'
with open(td_path, 'rb') as f: td = pickle.load(f)
with open(name_path, 'rb') as f: name_dict = pickle.load(f)
print('Token dict keys sample:', list(td.keys())[:5])
print('Name dict sample:', list(name_dict.items())[:5])
targets = ['CADM2','NRXN1','NLGN1','APP','FOXO3']
for t in targets:
    ens = name_dict.get(t, 'MISSING')
    tok = td.get(ens, 'MISSING') if ens != 'MISSING' else 'N/A'
    print(f'{t} -> {ens} -> token {tok}')
