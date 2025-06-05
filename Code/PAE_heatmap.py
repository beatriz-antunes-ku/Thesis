#from: https://github.com/Ash100/Biopython/blob/main/AF3_Results_Visualization.ipynb
#@title # **Predicted Alignment Error (PAE) Heatmap 🔥**
import json
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# ---Prepare Data---
# Function to load JSON data from a file
def extract_data(json_file_path):
    try:
        # Open the JSON file and load the data
        with open(json_file_path, 'r') as file:
            data = json.load(file)
        return data

    except Exception as e:
        print("Error reading the JSON file:", e)

# Function to calculate chain lengths
def calculate_chain_lengths(token_res_ids):
    chain_lengths = []
    current_chain_length = 1  # Initialize with 1 to account for the first residue
    previous_residue = token_res_ids[0]

    for residue in token_res_ids[1:]:
        if residue == previous_residue + 1:
            current_chain_length += 1
        else:
            chain_lengths.append(current_chain_length)
            current_chain_length = 1
        previous_residue = residue

    # Add length of the last chain
    chain_lengths.append(current_chain_length)

    return chain_lengths

# Put the data into arrays
data = extract_data("/Users/rnq676/Documents/thesis/KRAS - PIK3CG/fold_1he8/fold_1he8_full_data_4.json")
pae_data = data['pae']
token_res_ids = data['token_res_ids']

print('Number of residues: {}'.format(len(token_res_ids)))

# Calculate entity lengths
chain_lengths = calculate_chain_lengths(token_res_ids)
print('Entity lengths:', chain_lengths)

# Calculate total number of residues
total_residues = sum(chain_lengths)

# Create heatmap
plt.figure(figsize=(8, 6))

# Plot the heatmap for the entire data
pae_heatmap = sns.heatmap(pae_data, cmap='viridis', cbar_kws={'label': 'Expected Position Error (Ångströms)'}, square=True)

# Customize plot
plt.title('Predicted Alignment Error (PAE) Heatmap - 1he8 model 4')
plt.xlabel('Scored residue')
plt.ylabel('Aligned residue')

# Set tick positions to cover the full range of data
tick_positions = np.arange(total_residues)
tick_labels = [str(i+1) for i in range(total_residues)]  # Adjust for 1-based indexing

# Adjust tick frequency
num_ticks = 10  # Change this value to adjust the number of ticks
tick_step = total_residues // num_ticks
if tick_step == 0:  # Ensure at least one tick per interval
    tick_step = 1
plt.xticks(tick_positions[::tick_step], tick_labels[::tick_step], rotation=45)
plt.yticks(tick_positions[::tick_step], tick_labels[::tick_step])

# Show color legend
cbar = pae_heatmap.collections[0].colorbar
cbar.set_label('Expected Position Error (Ångströms)')

# Show the plot
plt.show()

pae_heatmap_image = pae_heatmap.get_figure()
