# Conotoxin Structure-Based Clustering Analysis

## Overview
This project analyzes the structural diversity of conotoxins, small peptides produced by marine cone snails, using AI-guided 3D protein structure predictions and clustering techniques. The pipeline leverages AlphaFold-generated PDB structures, PCA, and clustering to classify conotoxins based on structural features rather than sequence similarity.  

The analysis helps uncover evolutionary patterns, functional similarities, and structure-function relationships among conotoxins.

---

## Features
- Reads 20 representative conotoxin PDB files.  
- Extracts structural features:  
  - Atom count  
  - Radius of gyration  
- Performs PCA for dimensionality reduction.  
- Uses K-Means clustering to classify conotoxins.  
- Computes silhouette scores for clustering validation.  
- Generates multiple plots:
  - PCA scatter plot  
  - Bubble plot of structural features  
  - Silhouette plot  
  - Feature heatmap  
  - PCA explained variance  
  - Species-level cluster distribution  
  - Species PCA centroids  
  - Feature boxplots per species  
  - RMSD heatmap for pairwise structural similarity  
- Saves all results to CSV (`results/results_species.csv`) for downstream analysis.

---

## Requirements
- Python 3.10+  
- Packages:
  ```bash
  pip install numpy pandas matplotlib seaborn scikit-learn biopython

OR

conda create -n conotoxin python=3.10
conda activate conotoxin
pip install numpy pandas matplotlib seaborn scikit-learn biopython

## File Structure
Conotoxin-main/
├── main.py                    # Main script to run analysis
├── visualization/
│   └── pca_plot.py            # Plotting functions
├── pdb_files/                 # 20 conotoxin PDB files
├── results/                   # CSV output: results_species.csv
└── plots/                     # All generated plots

## Usage
1. Place all conotoxin PDB files in the pdb_files/ directory.
2. Run the main script:
    
   python main.py

3. The script outputs:
    - results/results_species.csv – table of features, PCA coordinates, cluster IDs, silhouette scores, and species.
    - plots/ – folder containing all charts and heatmaps.

## CSV Output (results_species.csv) contains the following columns:
    - Name – PDB ID or conotoxin name
    - Species – Cone snail species
    - Atom_Count – Number of atoms
    - Radius_of_Gyration – Structural compactness
    - PC1, PC2 – Principal component coordinates
    - Cluster_ID – Assigned cluster from K-Means
    - Silhouette – Silhouette score per sample

Notes:
    - The PCA scatter plot and bubble plot include clear labels and legends.
    - RMSD heatmap shows pairwise structural similarity.
    - Boxplots and species-level centroids summarize clustering across cone snail species.
    - Legends use real molecule symbols to improve clarity and reduce overlapping text in plots.
    - Silhouette score provides a measure of clustering quality; higher is better.
    - Plots are saved in plots/, CSV results in results/.