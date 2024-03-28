# +
import pandas as pd
import matplotlib.pyplot as plt

def plot_efficiency_graph(df, df_primer_efficiency, gene):
    # Filter the data for the specific gene
    gene_data = df[df['Gene'] == gene]
    
    # Filter the primer efficiency data for the specific gene
    gene_efficiency_data = df_primer_efficiency[df_primer_efficiency['Gene'] == gene]
    
    # Extracting relevant information
    log_dilution = gene_data['log_dilution']
    ct_value = gene_data['Ct_value']
    ct1 = gene_data['Ct1']
    ct2 = gene_data['Ct2']
    ct3 = gene_data['Ct3']
    slope = gene_efficiency_data['Slope'].values[0]  # Assuming there's only one row for each gene
    r_value = gene_efficiency_data['R'].values[0]  # Assuming there's only one row for each gene
    primer_efficiency = gene_efficiency_data['Primer Efficiency'].values[0]  # Assuming there's only one row for each gene
    sem = gene_data['SEM']  # Assuming SEM column contains standard error values
    
    # Plot the data
    plt.plot(log_dilution, ct_value, marker='o', color='blue', label=f'{gene} Ct_value')
    plt.scatter(log_dilution, ct1, color='red', label='Ct1', alpha=0.3, s=20, zorder=5)
    plt.scatter(log_dilution, ct2, color='green', label='Ct2', alpha=0.3, s=20, zorder=5)
    plt.scatter(log_dilution, ct3, color='orange', label='Ct3', alpha=0.3, s=20, zorder=5)
    
    # Plot the main dot (average Ct_value) with error bars (standard error)
    plt.errorbar(log_dilution, ct_value, yerr=sem, fmt='o', color='blue', ecolor='black', capsize=5)
    
    # Add labels
    plt.xlabel('log_dilution')
    plt.ylabel('Ct_value')

    # Add custom legend for slope, R_value, and primer efficiency
    legend_text = f'Slope: {slope:.2f}\nR_value: {r_value:.2f}\nPrimer Efficiency: {primer_efficiency:.2f}'
    plt.text(0.65, 0.95, legend_text, transform=plt.gca().transAxes, fontsize=10, verticalalignment='top', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=0.5'))

    # Show the plot
    plt.show()

# Example usage:
# plot_efficiency_graph(df_gene_data, df_primer_efficiency, 'GeneA')

