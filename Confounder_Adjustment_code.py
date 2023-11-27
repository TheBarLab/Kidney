import pandas as pd
from scipy import stats
import numpy as np

# Function to calculate partial correlations
def partial_correlation(r_xy, r_xz, r_yz):
    return (r_xy - r_xz * r_yz) / np.sqrt((1 - r_xz**2) * (1 - r_yz**2))

# Loading the data
file_path = "/path/to/your/file.csv"  # Replace with the actual path to your file
data = pd.read_csv(file_path)

# Variables of interest
variables_of_interest = ['Interstitium_Fibrosis', 'eGFR']
# Epigenetic clocks - adjust according to analysis
epigenetic_clocks = ['DNAmPhenoAge', 'DNAmB2M', 'DNAmGDF15', 'DNAmLeptin', 'DNAmPAI1', 'DNAmGrimAge']
# Confounders
confounders = ['Age', 'Gender', 'BMI', 'glucose', 'BP_systole', 'BP_diastole', 'glom_sclerosis_old', 'Race']

# DataFrame to hold the results with separate adjustments for each confounder
detailed_results = pd.DataFrame(columns=['Variable', 'Clock', 'R', 'p-value'] + [f'Adjusted R ({conf})' for conf in confounders] + [f'Adjusted p-value ({conf})' for conf in confounders])

# Calculating correlations with separate adjustments for each confounder
for variable in variables_of_interest:
    for clock in epigenetic_clocks:
        # Removing missing values for the unadjusted correlation
        clean_data_unadjusted = data[[variable, clock]].dropna()
        x = clean_data_unadjusted[variable]
        y = clean_data_unadjusted[clock]
        r, p_value = stats.pearsonr(x, y)

        # Storing results for this combination
        result_row = {'Variable': variable, 'Clock': clock, 'R': r, 'p-value': p_value}

        # Adjusting for each confounder separately
        for confounder in confounders:
            # Removing missing values for the adjusted correlation
            clean_data_adjusted = data[[variable, clock, confounder]].dropna()
            x_adj = clean_data_adjusted[variable]
            y_adj = clean_data_adjusted[clock]
            z = clean_data_adjusted[confounder]
            r_xz, _ = stats.pearsonr(x_adj, z)
            r_yz, _ = stats.pearsonr(y_adj, z)
            adjusted_r = partial_correlation(r, r_xz, r_yz)
            adjusted_p_value = p_value  # Placeholder value for adjusted p-value
            result_row[f'Adjusted R ({confounder})'] = adjusted_r
            result_row[f'Adjusted p-value ({confounder})'] = adjusted_p_value

        # Adding to detailed results
        detailed_results = detailed_results.append(result_row, ignore_index=True)


