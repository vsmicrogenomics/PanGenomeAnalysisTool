# PanGenomeAnalysisTool
PanGenomeAnalysisTool: A Python script for pan-genome analysis, generating plots, and statistical insights. Analyze gene presence and absence in multiple genomes effortlessly.

This Python script is designed for analyzing pan-genomes, specifically for estimating the parameters \(k\) and \(\gamma\) for the Heap Law equation: \(n = \kappa N^\gamma\), where \(n\) is the number of pangenome genes and \(N\) is the number of genomes. It takes the `gene_presence_absence.Rtab` output of the Roary tool as input and provides insights into whether the pan-genome is open or closed based on the estimated \(\gamma\) value.

## Features:
- Estimate the parameters \(k\) and \(\gamma\) for the Heap Law equation.
- Determine if the pan-genome is open or closed based on the \(\gamma\) value.
- Generate pan-genome and modified core-genome plots.
- Save statistics and plots as high-resolution images for publication.

## Usage

1. **Installation**:
   - Ensure you have Python 3.x installed on your system.
   - Install required Python packages using pip:
     ```bash
     pip install numpy matplotlib scipy
     ```

2. **Running the Script**:
   - Clone this repository to your local machine.
   - Navigate to the directory containing `pan_genome_analysis.py`.

3. **Command-line Usage**:
   - Run the script using the following command-line arguments:
     ```bash
     python pan_genome_analysis.py -f input_file -i iterations -o output_dir
     ```
     - `-f`, `--input_file`: Path to the `gene_presence_absence.Rtab` output file of the Roary tool.
     - `-i`, `--iterations`: Number of iterations for analysis (default: 10).
     - `-o`, `--output_dir`: Directory path for output files.

4. **Output**:
   - The script will generate the following outputs:
     - Pan-genome and modified core-genome plots in high-resolution image formats (e.g., PNG, PDF) saved in the specified output directory.
     - A text file named `pan_genome_statistics.txt` containing the values of \(k\), \(\gamma\), and the pan-genome status (open or closed).

5. **Interpreting Results**:
   - The script estimates \(k\) and \(\gamma\) for the Heap Law equation and determines whether the pan-genome is open or closed based on the \(\gamma\) value.

5. **Test run**:
  To demonstrate the usage of the `pan_genome_analysis.py` script, we provide a test directory with input and output subdirectories with test input and its output files. You can reproduce the analysis using the following command:

   ```bash
   python pan_genome_analysis.py -f input/gene_presence_absence.Rtab -o output -i 10

# Acknowledgements

This script is designed to process the output files of the Roary tool and analyze pan-genomes. We would like to acknowledge the developers of the Roary tool for their contribution to the field of comparative genomics. The Roary tool is a valuable resource for pan-genome analysis, and its documentation is available at [Roary GitHub Repository](https://github.com/sanger-pathogens/Roary).

If you are using the `pan_genome_analysis.py` script for your research, please consider citing it as follows:

```plaintext
Sharma, V. (2024). pan_genome_analysis.py [Python script]. Retrieved from https://github.com/vsmicrogenomics/PanGenomeAnalysisTool


**Acknowledgments**:
    - This script utilizes libraries such as NumPy, Matplotlib, and SciPy for data analysis and visualization.

Feel free to use and contribute to this pan-genome analysis tool for your research or bioinformatics projects.

