# +
def import_and_tidy_data(file_path):
    return df


def primer_efficiency_calc(slope):
    efficiency = (10^(-1/slope)-1 ) * 100
    return efficiency 

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
