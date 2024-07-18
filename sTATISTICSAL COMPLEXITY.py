import numpy as np
import pandas as pd
from collections import defaultdict
import os

def calculate_statistical_complexity(time_series, memory_length=5, tolerance=0.05):
    def get_probability_distributions(time_series, memory_length):
        history_counts = defaultdict(int)
        future_counts = defaultdict(lambda: defaultdict(int))
        
        for i in range(len(time_series) - memory_length):
            history = tuple(time_series[i:i+memory_length])
            future = time_series[i+memory_length]
            history_counts[history] += 1
            future_counts[history][future] += 1
        
        history_probabilities = {}
        for history, count in history_counts.items():
            future_probabilities = {future: count / history_counts[history] for future, count in future_counts[history].items()}
            history_probabilities[history] = future_probabilities
        
        return history_probabilities

    def refine_states(prob_distributions, tolerance):
        refined_states = {}
        state_index = 0
        
        for history, prob_dist in prob_distributions.items():
            found = False
            for state, state_dist in refined_states.items():
                if np.allclose(list(prob_dist.values()), list(state_dist.values()), atol=tolerance):
                    refined_states[state].update(prob_dist)
                    found = True
                    break
            
            if not found:
                refined_states[state_index] = prob_dist
                state_index += 1
        
        return refined_states

    def shannon_entropy(prob_dist):
        return -sum(p * np.log2(p) for p in prob_dist if p > 0)

    print(f"Original time series: {time_series[:10]} ...")  # Debug: show first 10 values
    prob_distributions = get_probability_distributions(time_series, memory_length)
    print(f"Probability Distributions: {prob_distributions}")  # Debug statement
    refined_states = refine_states(prob_distributions, tolerance)
    print(f"Refined States: {refined_states}")  # Debug statement
    
    if not refined_states:
        print("No refined states found.")
        return 0
    
    state_probabilities = [sum(probs.values()) for probs in refined_states.values()]
    
    # Normalize state probabilities
    total_probability = sum(state_probabilities)
    if total_probability > 0:
        state_probabilities = [p / total_probability for p in state_probabilities]
    else:
        print("Total probability is zero.")
        state_probabilities = [0 for p in state_probabilities]
    
    print("State Probabilities:", state_probabilities)  # Debug statement
    
    complexity = shannon_entropy(state_probabilities)
    
    return complexity

def read_time_series(file_path):
    try:
        data = pd.read_csv(file_path, header=None)
        # Ensure only numeric data is returned
        numeric_data = data.apply(pd.to_numeric, errors='coerce').dropna().values.flatten()
        return numeric_data
    except FileNotFoundError:
        print(f"Error: The file at {file_path} was not found.")
        return None

def binarize_time_series(time_series):
    median = np.median(time_series)
    binarized_series = (time_series > median).astype(int)
    return binarized_series

def calculate_complexity_from_file(file_path):
    time_series = read_time_series(file_path)
    if time_series is None or len(time_series) == 0:
        print(f"No valid data found in {file_path}")
        return None
    binarized_series = binarize_time_series(time_series)
    print(f"Binarized time series: {binarized_series[:10]} ...")  # Debug: show first 10 values of binarized series
    complexity = calculate_statistical_complexity(binarized_series)
    return complexity

def process_folder(folder_path):
    complexities = {}
    for filename in os.listdir(folder_path):
        if filename.endswith(".csv"):
            file_path = os.path.join(folder_path, filename)
            complexity = calculate_complexity_from_file(file_path)
            if complexity is not None:
                complexities[filename] = complexity
    return complexities

# Example usage
folder_path = r"C:\Users\prsyu\OneDrive\Bidlung\University\M.S. Leiden University\M.S. Neuroscience (Research)\Thesis\CopBET\CopBET\data"  # Replace with your actual folder path
complexities = process_folder(folder_path)

for filename, complexity in complexities.items():
    print(f'{filename}: Statistical Complexity: {complexity}')



#----#
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from collections import defaultdict

# Define the statistical complexity values for each file
data = {
    "dat-basel_LAM_sub-01_task-Amphetamine_ses-1_atlas-yeo17.csv": 3.0,
    "dat-basel_LAM_sub-01_task-LSD_ses-3_atlas-yeo17.csv": 3.459431618637298,
    "dat-basel_LAM_sub-01_task-MDMA_ses-4_atlas-yeo17.csv": 3.321928094887362,
    "dat-basel_LAM_sub-01_task-Placebo_ses-2_atlas-yeo17.csv": 3.321928094887362,
    "dat-basel_LAM_sub-02_task-Amphetamine_ses-3_atlas-yeo17.csv": 3.0,
    "dat-basel_LAM_sub-02_task-LSD_ses-4_atlas-yeo17.csv": 3.169925001442312,
    "dat-basel_LAM_sub-02_task-MDMA_ses-2_atlas-yeo17.csv": 3.0,
    "dat-basel_LAM_sub-02_task-Placebo_ses-1_atlas-yeo17.csv": 3.584962500721156,
    "dat-basel_LAM_sub-03_task-Amphetamine_ses-4_atlas-yeo17.csv": 3.321928094887362,
    "dat-basel_LAM_sub-03_task-LSD_ses-2_atlas-yeo17.csv": 3.459431618637298,
    "dat-basel_LAM_sub-03_task-MDMA_ses-1_atlas-yeo17.csv": 3.169925001442312,
    "dat-basel_LAM_sub-03_task-Placebo_ses-3_atlas-yeo17.csv": 3.0,
    "dat-basel_LAM_sub-04_task-Amphetamine_ses-3_atlas-yeo17.csv": 3.321928094887362,
    "dat-basel_LAM_sub-04_task-LSD_ses-4_atlas-yeo17.csv": 3.321928094887362,
    "dat-basel_LAM_sub-04_task-MDMA_ses-1_atlas-yeo17.csv": 3.459431618637298,
    "dat-basel_LAM_sub-04_task-Placebo_ses-2_atlas-yeo17.csv": 3.584962500721156,
}

# Initialize dictionaries to hold the sums and counts for each drug
drug_sums = defaultdict(float)
drug_counts = defaultdict(int)
drug_values = defaultdict(list)

# Aggregate the sums and counts for each drug
for filename, complexity in data.items():
    drug = filename.split("_task-")[1].split("_ses")[0]
    drug_sums[drug] += complexity
    drug_counts[drug] += 1
    drug_values[drug].append(complexity)

# Calculate the averages for each drug
drug_averages = {drug: drug_sums[drug] / drug_counts[drug] for drug in drug_sums}

# Print the results
print("Average Statistical Complexity for each drug:")
for drug, average in drug_averages.items():
    print(f"{drug}: {average:.4f}")

# Create a boxplot of the values
plt.figure(figsize=(10, 6))
plt.boxplot([drug_values[drug] for drug in drugs], labels=drugs)
plt.xlabel('Drug')
plt.ylabel('Statistical Complexity')
plt.title('Statistical Complexity for Each Drug')
plt.show()

# Perform statistical tests
# Using ANOVA to compare means across multiple groups
f_val, p_val = stats.f_oneway(*[drug_values[drug] for drug in drugs])
print(f"\nANOVA results: F-value = {f_val:.4f}, p-value = {p_val:.4f}")

# If ANOVA is significant, perform post-hoc t-tests
if p_val < 0.05:
    print("\nPost-hoc t-tests (Bonferroni corrected):")
    for i in range(len(drugs)):
        for j in range(i + 1, len(drugs)):
            drug1, drug2 = drugs[i], drugs[j]
            t_val, p_val = stats.ttest_ind(drug_values[drug1], drug_values[drug2])
            # Bonferroni correction
            p_val_corrected = p_val * (len(drugs) * (len(drugs) - 1) / 2)
            print(f"Comparison between {drug1} and {drug2}: t-value = {t_val:.4f}, p-value = {p_val_corrected:.4f}")

# Comparing average entropies directly
comparisons = []
for i in range(len(drugs)):
    for j in range(i + 1, len(drugs)):
        drug1, drug2 = drugs[i], drugs[j]
        diff = abs(drug_averages[drug1] - drug_averages[drug2])
        comparisons.append((drug1, drug2, diff))

# Print the comparisons
print("\nComparisons of Statistical Complexity between drugs:")
for drug1, drug2, diff in comparisons:
    print(f"Difference between {drug1} and {drug2}: {diff:.4f}")
