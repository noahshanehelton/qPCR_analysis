# +
#import numpy as np
#import pandas as pd 
#from scipy import stats

def import_and_tidy_data(file_path):
    '''
    Import csv file with the following headings:
    Gene
    Condition
    Replicate
    Ct1
    Ct2
    Ct3

    This will import the csv file and average technical Ct replicates 
    '''
    df = pd.read_csv(file_path)
    df['Ct_value'] = df[['Ct1', 'Ct2', 'Ct3']].mean(axis=1)
    return df
    
def primer_effiency_calc(slope):
    '''
    Equation for primer effieciency 
    '''
    exponent = -1/slope
    efficiency = (10**(exponent)-1 ) * 100
    return round(efficiency,2)



def primer_efficiency(file_path, gene):
    """
    Input csv file with the following headings 
    Gene
    Dilution
    Ct1
    Ct2
    Ct3

    This function will import csv file and return the primer efficiency 
    for a gene of interest by converting dilutions to log_dilutions, performing
    linear regression to calculate the slope, and doing the primer effeciency 
    calculation on the data
    """
    df = import_and_tidy_data(file_path)
    df = df[df['Gene'] == gene]
    df['log_dilution'] = np.log10(df['Dilution'])
    slope, intercept, r_value, p_value, std_err = stats.linregress(df['log_dilution'], df['Ct_value'])
    return primer_effiency_calc(slope)
#print(primer_efficiency('data/dilution_testdata.csv', 'GOI'))

def delta_Ct(GOI, control_gene):
    dCT = GOI - control
    return dCT

def delta_delta_Ct(Gene, experimental_condition, control_condition):
    ddCt = delta_Ct(Gene, control) - delta_Ct(Gene, control)
    return True

def pfaffl(GOI, control_gene, experimental_condition, control_condition, E_GOI, E_control):
    ratio = (E_GOI**(delta_Ct(GOI, experimental_condition, control_condition)))
    
    return df

def polysome_analyis(df):
    return percent_df
# -


