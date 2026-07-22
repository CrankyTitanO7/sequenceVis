### Overview
**SequenceVis** is a tool designed to bridge the gap between position-specific frequency matrices (from CSVs) and sequence-based visualization tools. It generates a foreground dataset reflecting specific motif probabilities and a constrained random background dataset.
https://colab.research.google.com/github/CrankyTitanO7/sequenceVis/blob/main/colab/sequenceVis.ipynb

### Usage Instructions
1. **Prepare your CSV**: Ensure your CSV has a `titles` column for amino acids and columns labeled `-5` to `+4` containing frequency or probability data.
2. **Upload**: Run the **Upload Your Own CSV** cell to import your data. The script will automatically rename it to `input.csv` for processing.
3. **Run Logic**: Execute the **Run the Conversion** cells. This will:
    * Generate 1,000 foreground sequences based on your matrix.
    * Generate 5,000 background sequences (random 10-mers, position 0 fixed to S/T).
4. **Download**: Use the final cell to download `fg.csv` and `bg.csv` for use in logo generators like WebLogo or IceLogo.

### Summary of Output Files
*   **`fg.csv`**: The primary output for motif visualization.
*   **`bg.csv`**: Used as a reference background to highlight enrichment.
*   **`fg_full.csv`**: Useful for debugging, showing the specific amino acid chosen at every position for every generated sequence.
