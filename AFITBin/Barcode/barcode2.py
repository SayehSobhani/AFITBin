import pandas as pd
import matplotlib.pyplot as plt

# Replace 'file_path.csv' with the path to your CSV file containing k-mers
file_path = 'your_file_path.csv'

# Load the CSV file into a Pandas DataFrame
data = pd.read_csv(file_path)

# Assuming the CSV file has a column named 'Sequence' containing the sequences
sequences = data['Sequence']

# Assigning unique colors to different k-mers for visualization
unique_kmers = set(sequences.str.split(',').explode().unique())
color_map = {kmer: f'C{i}' for i, kmer in enumerate(unique_kmers)}

# Generating the barcode representation
plt.figure(figsize=(10, 6))

for i, sequence in enumerate(sequences):
    y = [i] * len(sequence.split(','))
    x = [j for j in range(len(sequence.split(',')))]
    colors = [color_map[kmer] for kmer in sequence.split(',')]
    plt.scatter(x, y, color=colors, marker='|')

# Customize plot labels and appearance
plt.xlabel('Position in Sequence')
plt.ylabel('Sequence')
plt.title('Barcode Representation of Sequences based on k-mers')
plt.yticks(range(len(sequences)), range(1, len(sequences) + 1))
plt.ylim(-1, len(sequences))
plt.grid(axis='y')
plt.tight_layout()

# Display or save the plot
plt.show()
# Uncomment the following line to save the plot to a file
# plt.savefig('barcode_plot.png')

