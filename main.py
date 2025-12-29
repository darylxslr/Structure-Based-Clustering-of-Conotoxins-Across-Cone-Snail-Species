import os
import numpy as np
import pandas as pd
from features.extractor import extract_features
from preprocessing.scaler import scale_features
from models.pca_model import run_pca
from models.clustering_model import kmeans_clustering
from evaluation.silhouette import compute_silhouette
from visualization.pca_plot import generate_all_plots, generate_rmsd_matrix, plot_rmsd_heatmap, compute_rmsd

# -------------------------
# Configuration
# -------------------------
PDB_DIR = "pdb_files"
OUTPUT_DIR = "plots"
RESULTS_DIR = "results"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# Mapping of PDB files to cone snail species
species_mapping = {
    "1AKG":"Conus magus", "1F3K":"Conus geographus", "1FU3":"Conus textile", "1FYG":"Conus striatus",
    "1K64":"Conus catus", "1MII":"Conus purpurascens", "1MTQ":"Conus tulipa", "1P1P":"Conus marmoreus",
    "1QMW":"Conus obscurus", "2H8S":"Conus ventricosus", "2JQC":"Conus imperialis", "2YYF":"Conus victoriae",
    "8K3M":"Conus ermineus", "8K3N":"Conus radiatus", "EVIA":"Conus ermineus", "GI":"Conus geographus",
    "GIIIA":"Conus geographus", "ImI":"Conus imperialis", "MVIIA":"Conus magus", "PVIIA":"Conus purpurascens"
}

# -------------------------
# 1. Load PDB files & extract features
# -------------------------
files = sorted(os.listdir(PDB_DIR))
names = []
features = []

for file in files:
    if file.endswith(".pdb"):
        names.append(file.replace(".pdb",""))
        feat = extract_features(os.path.join(PDB_DIR, file))
        features.append(feat)

features = np.array(features)

if len(features) != 20:
    print(f"Warning: Expected 20 PDB files, found {len(features)}")
else:
    print("All 20 PDB files successfully read.")

# -------------------------
# 2. Scale features
# -------------------------
scaled_features = scale_features(features)

# -------------------------
# 3. PCA
# -------------------------
pca_data, pca_model = run_pca(scaled_features)

# -------------------------
# 4. K-Means clustering
# -------------------------
labels, _ = kmeans_clustering(scaled_features, n_clusters=4)

# -------------------------
# 5. Silhouette score
# -------------------------
silhouette_avg = compute_silhouette(scaled_features, labels)
print(f"Silhouette Score: {silhouette_avg:.3f}")

# -------------------------
# 6. Generate all plots
# -------------------------
generate_all_plots(
    features=features,
    scaled_features=scaled_features,
    pca_data=pca_data,
    pca_model=pca_model,
    labels=labels,
    names=names,
    silhouette_avg=silhouette_avg,
    species_mapping=species_mapping,
    output_dir=OUTPUT_DIR
)

# -------------------------
# 7. Compute RMSD & heatmap
# -------------------------
rmsd_matrix = generate_rmsd_matrix(PDB_DIR, names)
plot_rmsd_heatmap(rmsd_matrix, names, species_mapping, output_dir=OUTPUT_DIR)
print("RMSD heatmap generated and saved to 'plots/'.")

# -------------------------
# 8. Save results_species.csv (with RMSD)
# -------------------------
results_species = pd.DataFrame({
    "Conotoxin": names,
    "Species": [species_mapping.get(n,"Unknown") for n in names],
    "Cluster": labels,
    "PC1": pca_data[:,0],
    "PC2": pca_data[:,1],
    "Atom Count": features[:,0],
    "Radius of Gyration": features[:,1]
})

# Add RMSD values to all others
for i, name in enumerate(names):
    results_species[[f"RMSD_to_{n}" for n in names]] = rmsd_matrix

results_species.to_csv(os.path.join(RESULTS_DIR,"results_species.csv"), index=False)
print(f"results_species.csv saved to '{RESULTS_DIR}/'.")
