
import numpy as np
import pandas as pd 
from scipy import stats

def import_and_tidy_data(file_path):
    '''
    Import csv file with the following headings:
    Gene
    Condition or Dilution or Condition and Fraction 
    Replicate
    Ct1
    Ct2
    Ct3

    This will import the csv file and average technical Ct replicates 
    '''
    df = pd.read_csv(file_path)

    expected_headers1 = ['Gene', 'Condition', 'Replicate', 'Ct1', 'Ct2', 'Ct3'] #for GER tests
    expected_headers2 = ['Gene', 'Dilution', 'Replicate', 'Ct1', 'Ct2', 'Ct3']  #for dilution series and    
    expected_headers3 = ['Gene', 'Fraction', 'Condition', 'Ct1', 'Ct2', 'Ct3'] #for polysome profiling qPCR

    df['Ct_value'] = df[['Ct1', 'Ct2', 'Ct3']].mean(axis=1)    
    df['SEM'] = df[['Ct1', 'Ct2', 'Ct3']].sem(axis=1)

    return df
    
def primer_effiency_calc(slope):
    '''
    Equation for primer effieciency 
    '''
    exponent = -1/slope
    efficiency = (10**(exponent)-1 ) * 100
    return round(efficiency,2)



def primer_efficiency(df, gene):
    """
    Input dataframe from import and tidy data 

    This function will import csv file and return the primer efficiency 
    for a gene of interest by converting dilutions to log_dilutions, performing
    linear regression to calculate the slope, and doing the primer effeciency 
    calculation on the data
    """
    #Convert dilution series to log10
    df['log_dilution'] = np.log10(df['Dilution'])
    
    #Use Scipy to calculate linear regression of dilution series
    slope, intercept, r_value, p_value, std_err = stats.linregress(df['log_dilution'], df['Ct_value'])

    #Convert data to a pandas dataframe
    primer_efficiency_df = pd.DataFrame({
    'Gene': gene,
    'Slope': slope,
    'Intercept': intercept,
    'Error': std_err,
    'R': r_value, 
    'p_value': p_value,
    'Primer Efficiency': [primer_effiency_calc(slope)]
        })
    
    return primer_efficiency_df

def delta_Ct(GOI, control_gene):
    '''
    Goal: 

    Input: 

    Output: 
    '''
    dCT = GOI - control_gene

    return dCT

def delta_delta_Ct(Gene, experimental_condition, control_condition):
    '''
    Goal: 

    Input: 

    Output: 
    '''
    ddCt = delta_Ct(Gene, control) - delta_Ct(Gene, control)

    return True

def pfaffl_calc(deltaCT1, deltaCT2, efficiency1, efficiency2):
    '''
    Goal: 

    Input: 

    Output: 
    '''
    efficiency1 = (efficiency1/100) + 1
    efficiency2 = (efficiency2/100) + 1
    ratio = (efficiency1**deltaCT1/efficiency2**deltaCT2)
    return ratio

def pfaffl(GOI, control_gene, experimental_condition, control_condition, E_GOI, E_control, df):
    '''
    Goal: 

    Input: 

    Output: 


    
    '''
    #Get control averages 
    average_GOI = df[(df['Condition'] == control_condition) & (df['Gene'] == control_gene)]['Ct_value'].mean()
    average_control = df[(df['Condition'] == control_condition) & (df['Gene'] == GOI)]['Ct_value'].mean()

    # Create a new columns 
    df['DeltaCt'] = None  # Initialize DeltaCt column with None values

    df.loc[df['Gene'] == control_gene, 'DeltaCt'] = df['Ct_value'] - average_control
    df.loc[df['Gene'] == GOI, 'DeltaCt'] = df['Ct_value'] - average_GOI
    df = df[['Gene', 'Condition', 'DeltaCt','Replicate']]
    df_goi = df[df['Gene'] == 'GOI']  # DataFrame containing rows where 'Gene' is 'GOI'
    df_control = df[df['Gene'] == 'Control']  # DataFrame containing rows where 'Gene' is 'Control'
    df_goi.reset_index(drop=True, inplace=True)
    df_control.reset_index(drop=True, inplace=True)

    GER1 = []
    GER2 = []
    GER3 = []
    for i in range(len(df_goi)):
        if (df_goi.loc[i, 'Replicate'] == df_control.loc[i, 'Replicate']) and (df_goi.loc[i, 'Condition'] == df_control.loc[i, 'Condition']):

            if df_goi.loc[i, 'Replicate'] == 1: 
                GER1.append(pfaffl_calc(df_goi['DeltaCt'][i], df_control['DeltaCt'][i], E_GOI, E_control)) 
            if df_goi.loc[i, 'Replicate'] == 2: 
                GER2.append(pfaffl_calc(df_goi['DeltaCt'][i], df_control['DeltaCt'][i], E_GOI, E_control)) 
            if df_goi.loc[i, 'Replicate'] == 3: 
                GER3.append(pfaffl_calc(df_goi['DeltaCt'][i], df_control['DeltaCt'][i], E_GOI, E_control)) 
    
    data = [GER1, GER2, GER3]
    sem_value = stats.sem(data, axis=0)
    mean_value = np.mean(data, axis=0)



    pfaffl_df = pd.DataFrame({
         'Gene': [GOI, GOI],
         'Condition': [control_condition, experimental_condition],
         'Gene Expression Ratio 1':GER1 , 
         'Gene Expression Ratio 2':GER2 , 
         'Gene Expression Ratio 3':GER3 , 
         'Average GER':mean_value , 
         'SEM GER': sem_value
     }) 

    return pfaffl_df


#polysome_test = import_and_tidy_data('data/polysome_profile_testdata.csv')
#print(polysome_test)

def calculate_polysome_diff(reference_cts, cts):
    return reference_cts - cts
    
def calculate_2_delta_ct_polysome(delta_ct):
    return 2 ** delta_ct

def calculate_percentages(total, column):
    return (column * 100) / total if total != 0 else 0

def polysome_profiling_analysis(df, gene, condition, reps):
    
    subset_df = df[(df['Condition'] == condition) & (df['Gene'] == gene)].copy()

    # Calculate the baseline CT values for Fraction 1 to be used as reference
    baseline_cts = subset_df.loc[subset_df['Fraction'] == 1, ['Ct1', 'Ct2', 'Ct3']].mean()

    # Calculate delta Ct, 2^-deltaCT, and percentages for each Ct replicate
    for r in range(1, reps+1):
        replicate = f'Ct{r}'
        delta_ct_col = f'deltaCtR{r}'
        two_delta_ct_col = f'2DeltaCtR{r}'
        percent_col = f'Percent in fraction R{r}'
        
        # Calculate delta Ct
        delta_ct = calculate_polysome_diff(baseline_cts[replicate], subset_df[replicate])
        subset_df.loc[:, delta_ct_col] = delta_ct
        
        # Calculate 2^-deltaCT
        subset_df.loc[:, two_delta_ct_col] = calculate_2_delta_ct_polysome(delta_ct)
        
        # Calculate the percentage for each replicate
        sum_2_delta_ct = subset_df[two_delta_ct_col].sum()
        subset_df.loc[:, percent_col] = calculate_percentages(sum_2_delta_ct, subset_df[two_delta_ct_col])
        print(subset_df[two_delta_ct_col])
    
    
    
    subset_df = subset_df[subset_df.columns[subset_df.columns.str.contains('percent', case=False) | (subset_df.columns == 'fraction')]]
    row_averages = subset_df.mean(axis=1)
    row_sem = subset_df.sem(axis=1)

    subset_df['Average Percent in Fraction'] = row_averages
    subset_df['SEM Percent in Fraction'] = row_sem

    # Return the modified dataframe
    return subset_df.reset_index(drop=True)


#polysome_profiling_analysis(polysome_test, 'GOI','Untreated', 3) 
