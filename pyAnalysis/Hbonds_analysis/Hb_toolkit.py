import re
from pathlib import Path
import pandas as pd

def get_match(text_for_match):
    '''
    text_for_match: raw H-bond files from chimeraX

    returns: a list containing matches
    '''
    pat = r"/(A|B|C) *(\w+) *(\d+) *(\S+) */(A|B|C) *(\w+) *(\d+) *(\S+) *\w+ \w+ *(\d.\d+) *"

    match = re.findall(pat, text_for_match)

    match_list = []

    for tuples in match:
        match_list.append(list(tuples))

    return match_list

def to_df(matched_list):
    '''
    matched_list: A list containing matches using regular expression. Returned value from get_match()

    returns: dataframe
    '''

    df = pd.DataFrame(matched_list, columns = ['donor_chain', 'donor_residue_ID', 'donor_seq', 'donor_atom', 'acceptor_chain','acceptor_residue_ID', 'acceptor_seq','acceptor_atom','distance'])

    df_dtype_corrected = df.copy()
    # Convert residue_seq in df_dtype_corrected to int64 data type
    df_dtype_corrected["donor_seq"]    = df_dtype_corrected["donor_seq"].astype("int64")
    df_dtype_corrected["acceptor_seq"] = df_dtype_corrected["acceptor_seq"].astype("int64")
    df_dtype_corrected["distance"]     = df_dtype_corrected["distance"].astype("float64")

    return df_dtype_corrected

def valid_resi(df):
    '''
    df: dataframe
    
    returns: A dataframe contains only 20 common amino acids and the chain A
    '''

    amino_acids = ["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
    condition = (df["donor_residue_ID"].isin(amino_acids)) & (df["acceptor_residue_ID"].isin(amino_acids))
    
    df_valid = df.loc[condition]

    # select only chain A

    df_valid = df_valid.loc[(df_valid["donor_chain"] == "A") & (df_valid["acceptor_chain"] == "A")]

    return df_valid

def input_to_df(input_text):
    '''
    This function calls above funtions.
    
    input_text: raw H-bond files from chimeraX

    returns: cleaned dataframe
    '''

    match_text = get_match(input_text)
    df         = to_df(match_text)
    df_valid   = valid_resi(df)

    return df_valid

# Define a function for domain merge

def merge_domains(df, df_domain):
    '''
    df: dataframe

    df_domain: A dataframe containing domain information of each residue

    returns: A merged dataframe
    '''

    # Duplicate domain columns for merging
    df_domain["donor_domain"]    = df_domain["domain"]
    df_domain["acceptor_domain"] = df_domain["domain"]
    df_domain["donor_seq"]       = df_domain["residue_seq"]
    df_domain["acceptor_seq"]    = df_domain["residue_seq"]

    # Merge domain info based on donor_seq and acceptor_seq
    df_merge_donor_domain = pd.merge(df, df_domain[["donor_seq","donor_domain"]],
                                     on="donor_seq",
                                     how="left")
    
    df_merge_acceptor_domain = pd.merge(df_merge_donor_domain, df_domain[["acceptor_seq","acceptor_domain"]],
                                     on="acceptor_seq",
                                     how="left")  

    return df_merge_acceptor_domain

if __name__ == "__main__":
    main()