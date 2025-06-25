#%%
from ftplib import FTP
import os
import pandas as pd
import gzip
import shutil
from pyteomics import fasta, parser
import digest_find_siblings as sib

#uniprot_server = "ftp.uniprot.org"
#uniprot_dir = "pub/databases/uniprot"

#uniprot_server = "ftp.ebi.ac.uk"
#uniprot_dir = "pub/databases/uniprot"

uniprot_server = "ftp.expasy.org"
uniprot_dir = "databases/uniprot"

#%%

def download_file_ftp(host, remote_path, local_path, username='anonymous', password=''):
    with FTP(host) as ftp:
        ftp.login(user=username, passwd=password)
        with open(local_path, 'wb') as f:
            ftp.retrbinary(f'RETR {remote_path}', f.write)

#%% download proteomes README
os.makedirs('./work', exist_ok=True)
readme_file = './work/proteomes_README'

# download the current proteomes README from UniProt FTP server
download_file_ftp(uniprot_server, f'{uniprot_dir}/current_release/knowledgebase/reference_proteomes/README', readme_file)

# %% read in the information from the README file
in_rel_info = False
in_data = False
release_str = None
df_proteomes = None

with open(readme_file, 'r') as f:
    for line in f.readlines():
        if in_data:
            data = line.strip().split('\t')
            if len(data) == len(df_proteomes.columns):
                df_proteomes.loc[len(df_proteomes)] = data
            else:
                if len(data) <= 1:
                    in_data = False
            
        elif in_rel_info:
            if line.startswith('Release '):
                release = line
            elif line.startswith('Proteome_ID'):
                headers = line.strip().split('\t')
                df_proteomes = pd.DataFrame(columns=headers)
                in_data = True
        elif line.startswith('========================================================================'):
            in_rel_info = True

df_proteomes['#(1)'] = df_proteomes['#(1)'].astype(int)
df_proteomes['#(2)'] = df_proteomes['#(2)'].astype(int)

print(f"Found {len(df_proteomes)} proteomes in release {release.strip()}")

# %% cutting for testing
df_proteomes = df_proteomes[df_proteomes['SUPERREGNUM'] == 'eukaryota'][:5]

# %%
# go through the proteomes, download them and calculate the sibling information

proteome_count = 0
for index, row in df_proteomes.iterrows():
    proteome_id = row['Proteome_ID']
    superregnum = row['SUPERREGNUM'].capitalize()
    tax_id = row['Tax_ID']
    print(f"Downloading proteome {proteome_id}")

    fasta_file = f"./work/{proteome_id}_{tax_id}.fasta"
    siblings_file = f"./work/{proteome_id}_{tax_id}_siblings.txt"

    if row['#(1)'] > 0:
        ftp_str = f"{uniprot_dir}/current_release/knowledgebase/reference_proteomes/{superregnum}/{proteome_id}/{proteome_id}_{tax_id}.fasta.gz"
        proteome_file = f"./work/{proteome_id}_{tax_id}.fasta.gz"
        print(f"    {row['#(1)']} sequences")
        download_file_ftp(uniprot_server, ftp_str, proteome_file)
        with gzip.open(proteome_file, 'rb') as f_in, open(fasta_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(proteome_file)


    if row['#(2)'] > 0:
        ftp_str = f"{uniprot_dir}/current_release/knowledgebase/reference_proteomes/{superregnum}/{proteome_id}/{proteome_id}_{tax_id}_additional.fasta.gz"
        additional_proteome_file = f"./work/{proteome_id}_{tax_id}_additional.fasta.gz"
        print(f"    {row['#(2)']} additionals")
        download_file_ftp(uniprot_server, ftp_str, additional_proteome_file)
        with gzip.open(additional_proteome_file, 'rb') as f_in, open(fasta_file, 'ab') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(additional_proteome_file)
    
    peps_by_len = sib.digest_fasta_keep_with_leucines(fasta_file)

    count_siblings = 0
    count_siblings_per_length = {}
    for length, peps in peps_by_len.items():
        count_siblings_per_length[length] = 0
        for group_pep in list(peps.keys()):
            if len(peps[group_pep]) == 1:
                del peps[group_pep]
            else:
                count_siblings += len(peps[group_pep])
                count_siblings_per_length[length] += len(peps[group_pep])
    
    with open(siblings_file, 'w') as f:
        f.write(f"Total number of sibling peptides: {count_siblings}\n")
        print(f"Total number of sibling peptides: {count_siblings}")

        f.write("============================================================\n")
        for length in sorted(count_siblings_per_length.keys()):
            if (count_siblings_per_length[length] > 0):
                f.write(f"{length}\t{count_siblings_per_length[length]}\n")

        f.write("============================================================\n")
        for length in sorted(peps_by_len.keys()):
            for group_pep in peps_by_len[length]:
                f.write(f"{peps_by_len[length][group_pep]}\n")
    
    os.remove(fasta_file)

    proteome_count += 1
    print(f"Processed proteome {proteome_count} of {len(df_proteomes)}")


# %%
