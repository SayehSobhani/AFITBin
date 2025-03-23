import pandas as pd
import matplotlib.pyplot as plt

# Load your CSV file containing k-mer frequencies (replace 'filename.csv' with your file name)
data = pd.read_csv('/home/sayeh/fasta_BMC/fol/results.csv', index_col=0)  # Assuming the first column contains sequence IDs

# Plotting the barcode
plt.figure(figsize=(10, 6))  # Adjust the figure size as needed
plt.imshow(data.T, cmap='viridis', aspect='auto')
plt.xlabel('Sequences')  # Label for x-axis (if needed)
plt.ylabel('K-mer Frequencies')  # Label for y-axis (if needed)
plt.title('K-mer Barcode')  # Title of the plot (if needed)
plt.colorbar(label='Frequency')  # Add a color bar indicating frequency values

# Show or save the plot
plt.tight_layout()  # Adjust layout for better appearance
plt.show()  # Show the plot

# If you want to save the plot as an image file (e.g., PNG), uncomment the line below
# plt.savefig('barcode_plot.png', dpi=300)  # Specify the file name and DPI (adjust as needed)

