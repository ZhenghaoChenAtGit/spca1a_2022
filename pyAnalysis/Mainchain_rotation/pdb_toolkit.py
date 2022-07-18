import os
from pdbecif.mmcif_io import CifFileReader
import pandas as pd
import numpy as np
import re

def extract_coor(cifs, entry_ID):
    '''
    cifs: cif file read by CifFileReader(), For example:
          cifs = CifFileReader().read(f"6LLE.cif")
    
    entry_ID: the ID of cif file. For example: "6LLE"

    returns: a dataframe containing coordination 
    '''
    coordinates = cifs[entry_ID]["_atom_site"]
    df_coor = pd.DataFrame.from_dict(coordinates, orient="index").transpose() 
    # some may lack info in columns, use transpose incase  
    df_renamed = df_coor.rename(columns = {"auth_asym_id"  : "residue_chain",
                                           "label_atom_id" : "atom_ID",
                                           "label_comp_id" : "residue_ID",
                                           "auth_seq_id"   : "residue_seq",
                                           "Cartn_x"       : "coor_x" ,
                                           "Cartn_y"       : "coor_y" ,
                                           "Cartn_z"       : "coor_z" })
    
    column_for_coor   = ["residue_chain",
                         "residue_ID",
                         "residue_seq",
                         "atom_ID",
                         "coor_x" ,
                         "coor_y" ,
                         "coor_z" ]
    # select columns
    df_coor_filtered = df_renamed[column_for_coor].copy()
    # Convert residue_seq in df_coor_filtered to int64 data type
    df_coor_filtered["residue_seq"] = df_coor_filtered["residue_seq"].astype("int64")
    df_coor_filtered["coor_x"] = df_coor_filtered["coor_x"].astype("float64")
    df_coor_filtered["coor_y"] = df_coor_filtered["coor_y"].astype("float64")
    df_coor_filtered["coor_z"] = df_coor_filtered["coor_z"].astype("float64")

    return df_coor_filtered

def valid_resi(df):
    '''
    df: a dataframe

    returns: the dataframe containing atoms belonging to 20 amino acids
    '''

    amino_acids = ["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
    condition = df["residue_ID"].isin(amino_acids)
    df_valid = df.loc[condition]

    return df_valid


def alpha_C(df):
    '''
    df: a dataframe

    returns: a dataframe containing only alpha carbon atoms
    '''

    df_alpha = df.loc[df["atom_ID"] == "CA"]
    df_renemad = df_alpha.rename(columns = {"coor_x" : "alpha_x",
                                            "coor_y" : "alpha_y",
                                            "coor_z" : "alpha_z",})
    return df_renemad


if __name__ == "__main__":
    main()