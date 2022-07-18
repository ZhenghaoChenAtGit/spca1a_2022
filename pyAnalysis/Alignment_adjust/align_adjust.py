##################################################
# This script does not align sequences.          #
# This script merely adjust the position of "-"  #
# and add missing residues for displaying        #
# multiple sequence alignment.                   #
#                                                #
# Missing residues:                              #
# SPCA1a : 68-70   "---"    to  "DEP"            #
#          205-210 "------" to  "ATNGDL"         #
# SERCA2b: 502-507 "P---SS" to  "PSRT-S"         #
##################################################


with open("spca1_s2b.dnas", "r") as seq_set1:
    content_set1 = seq_set1.readlines()

# >6LLE.A (A refers to chain A)
set1_s2b = content_set1[1]
# >Mb_Ca_2.cif.A
set1_spca1 = content_set1[3]


with open("s2b_s1a.dnas", "r") as seq_set2:
    content_set2 = seq_set2.readlines()

# >6LLE.A
set2_s2b = content_set2[1]
# >3AR2.A
set2_s1a = content_set2[3]

with open("s2b_pmr1.dnas", "r") as seq_set3:
    content_set3 = seq_set3.readlines()

# >6LLE.A
set3_s2b = content_set3[1]
# >PMR1_Alphafold.pdb.A
set3_pmr1 = content_set3[3]

set1_s2b   = list(set1_s2b)
set1_spca1 = list(set1_spca1)
set2_s2b   = list(set2_s2b)
set2_s1a   = list(set2_s1a)
set3_s2b   = list(set3_s2b)
set3_pmr1  = list(set3_pmr1)

# Add "-" to make lengths the same

for i in range(len(set3_s2b)):
    if set3_s2b[i] != set2_s2b[i]:
        # When found differnt s2b seq, insert a "-" to other sets.
        if   set3_s2b[i] == "-":
            set2_s2b.insert(i, "-")
            set2_s1a.insert(i, "-")
        elif set2_s2b[i] == "-":
            set3_s2b.insert(i, "-")
            set3_pmr1.insert(i, "-")
        i += 1
    else:
        pass

for i in range(len(set3_s2b)):
    if set3_s2b[i] != set1_s2b[i]:
        if   set3_s2b[i] == "-":
            set1_s2b.insert(i, "-")
            set1_spca1.insert(i, "-")
        elif set1_s2b[i] == "-":
            set3_s2b.insert(i, "-")
            set3_pmr1.insert(i, "-")
        i += 1
    else:
        pass

for i in range(len(set3_s2b)):
    if set3_s2b[i] != set2_s2b[i]:
        if   set3_s2b[i] == "-":
            set2_s2b.insert(i, "-")
            set2_s1a.insert(i, "-")
        elif set2_s2b[i] == "-":
            set3_s2b.insert(i, "-")
            set3_pmr1.insert(i, "-")
        i += 1
    else:
        pass

# Add missing residues 68-70 "DEP", 205-210 "ATNGDL" in SPCA1a
set1_spca1[47:50]   = ["D", "E", "P"]
set1_spca1[200:206] = ["A", "T", "N", "G", "D", "L"]

# Add missing residues 68-70 "DEP", 205-210 "ATNGDL" in SPCA1a
set1_s2b[516:522]   = ["P", "S", "R", "T", "-", "S"]

# convert list to string
set1_s2b_op   = "".join(set1_s2b)
set1_spca1_op = "".join(set1_spca1)
set2_s2b_op   = "".join(set2_s2b)
set2_s1a_op   = "".join(set2_s1a)
set3_s2b_op   = "".join(set3_s2b)
set3_pmr1_op  = "".join(set3_pmr1)


# Add the following to the indicated sequence, add "-" to make their lengths the same.
# N-ter
spca1_n = "------------MKVARFQKIPNGENETMIPVLT"
pmr1_n  = "MSDNPFNASLLDEDSNREREILDATAEALSKPSP"

# C-ter
spca1_c = "IQKHVSSTSSSFLEV-----------------------------------"
pmr1_c  = "EDSTYFSNV-----------------------------------------"
s1a_c   = "G-------------------------------------------------"
s2b_c   = "PGKECVQPATKSCSFSACTDGISWPFVLLIMPLVIWVYSTDTNFSDMFWS"

adjusted_spca1 = spca1_n + set1_spca1_op + spca1_c
adjusted_pmr1  = pmr1_n  + set3_pmr1_op  + pmr1_c
adjusted_s1a   = "-"*34  + set2_s1a_op   + s1a_c
adjusted_s2b   = "-"*34  + set1_s2b_op   + s2b_c

# Summary
summary = f'''>SPCA1a
{adjusted_spca1}
>PMR1
{adjusted_pmr1}
>SERCA1a
{adjusted_s1a}
>SERCA2b
{adjusted_s2b}'''

# Write out the summary
with open("align_output.fasta", 'w') as out:
    out.write(summary)