

### AFITBin

Description: Briefly describe the project and its primary goals.

#### Files Included:

1. **AFITBIN.py**: Contains code related to command-line argument parsing using argparse.
2. **binning.py**: Includes functions for clustering sequences and guessing cluster numbers using KMeans.
3. **composition.py**: Defines functions for calculating novel and nucleotide frequencies from DNA sequences.
4. **validation.py**: Contains functions for plotting t-SNE visualizations and distance distributions.

### Dependencies:

- Python 3.x
- Necessary Python libraries:
  - BioPython (`Bio`)
  - NumPy (`numpy`)
  - Pandas (`pandas`)
  - SciKit-Learn (`sklearn`)
  - Matplotlib (`matplotlib`)
  - Seaborn (`seaborn`) (for some visualization functions)
  - Plotly (`plotly`) (for interactive visualizations)

### Usage Instructions:

1. **Installation:**
    - Make sure Python 3.x is installed on your system.
    - Install necessary libraries using `pip`:
        ```bash
        pip install biopython numpy pandas scikit-learn matplotlib seaborn plotly
        ```

2. **Running the code:**
    - Each file contains distinct functions. Import necessary functions into your Python script or Jupyter Notebook.
    - Ensure proper data input is provided as per each function's requirements.
    - Refer to the inline comments/documentation in each function for guidance on parameters and usage.
    - Example usage might look like this:
        ```python
        # Import necessary functions
        from AFITBIN import *
        from binning import guess_clusters_num
        from composition import novel_frequency, nucleotide_frequency, combine_frequency
        from validation import plot_tsne, plot_dist

        # Your data - replace this with your actual data
        sequences = [...]  # Your list of DNA sequences or the path to your sequences file

        # Example usage of functions
        clusters_num = guess_clusters_num(sequences)
        novel_freq = novel_frequency(sequences)
        nucleo_freq = nucleotide_frequency(sequences)
        combined_freq = combine_frequency(sequences)
        ```



