import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import silhouette_samples
from Bio.PDB import PDBParser

sns.set(style="whitegrid", font_scale=1.1)

# -------------------------
# Main plotting function
# -------------------------
def generate_all_plots(features, scaled_features, pca_data, pca_model, labels, names, silhouette_avg, species_mapping, output_dir="plots", results_dir="results"):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    cmap = "tab10"

    # -------------------------------
    # 1. PCA Scatter Plot
    # -------------------------------
    plt.figure(figsize=(10,8))
    scatter = plt.scatter(
        pca_data[:,0], pca_data[:,1], 
        c=labels, cmap=cmap, s=150, edgecolor="black", alpha=0.9
    )

    for i, name in enumerate(names):
        plt.text(
            pca_data[i,0]+0.02, pca_data[i,1]+0.02,
            name, fontsize=10, weight='bold',
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='gray', boxstyle='round,pad=0.2')
        )

    plt.xlabel("Principal Component 1", fontsize=12, weight='bold')
    plt.ylabel("Principal Component 2", fontsize=12, weight='bold')
    plt.title("Structure-Based Clustering of Conotoxins (PCA)", fontsize=14, weight='bold')
    plt.colorbar(scatter, label="Cluster ID")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"01_pca_scatter.png"), dpi=300)
    plt.close()

    # -------------------------------
    # 2. Bubble Plot
    # -------------------------------
    plt.figure(figsize=(10,8))
    sizes = (features[:,0] - np.min(features[:,0])) / (np.max(features[:,0])-np.min(features[:,0])+1e-5)*300 + 50
    scatter = plt.scatter(
        features[:,0], features[:,1],
        s=sizes, c=labels, cmap=cmap, edgecolor="black", alpha=0.8
    )
    for i, name in enumerate(names):
        plt.text(
            features[i,0]+1, features[i,1]+0.02,
            name, fontsize=9, weight='bold',
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='gray', boxstyle='round,pad=0.2')
        )
    plt.xlabel("Atom Count", fontsize=12, weight='bold')
    plt.ylabel("Radius of Gyration", fontsize=12, weight='bold')
    plt.title("Structural Features of Conotoxins (Bubble Plot)", fontsize=14, weight='bold')
    plt.colorbar(scatter, label="Cluster ID")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"02_feature_bubble.png"), dpi=300)
    plt.close()

    # -------------------------------
    # 3. Silhouette Plot
    # -------------------------------
    if len(np.unique(labels)) > 1:
        sample_silhouette = silhouette_samples(scaled_features, labels)
        plt.figure(figsize=(8,6))
        y_lower = 10
        for i in np.unique(labels):
            ith_vals = sample_silhouette[labels==i]
            ith_vals.sort()
            size_cluster_i = ith_vals.shape[0]
            y_upper = y_lower + size_cluster_i
            plt.barh(range(y_lower,y_upper), ith_vals, edgecolor="black", height=1.0, alpha=0.8)
            y_lower = y_upper + 10
        plt.axvline(x=sample_silhouette.mean(), linestyle="--", color="red",
                    label=f"Average = {sample_silhouette.mean():.2f}")
        plt.xlabel("Silhouette Coefficient", fontsize=12, weight='bold')
        plt.ylabel("Conotoxins", fontsize=12, weight='bold')
        plt.title("Silhouette Plot (Clustering Validation)", fontsize=14, weight='bold')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir,"03_silhouette.png"), dpi=300)
        plt.close()

    # -------------------------------
    # 4. Feature Heatmap
    # -------------------------------
    df_features = pd.DataFrame(features, index=names, columns=["Atom Count","Radius of Gyration"])
    plt.figure(figsize=(10,6))
    sns.heatmap(df_features, annot=True, fmt=".1f", cmap="viridis", linewidths=0.5, annot_kws={"fontsize":10})
    plt.title("Structural Feature Heatmap of Conotoxins", fontsize=14, weight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"04_feature_heatmap.png"), dpi=300)
    plt.close()

    # -------------------------------
    # 5. PCA Explained Variance
    # -------------------------------
    plt.figure(figsize=(7,5))
    cum_var = np.cumsum(pca_model.explained_variance_ratio_)
    plt.plot(range(1,len(cum_var)+1), cum_var, marker="o", linestyle="-", color="blue")
    plt.xticks(range(1,len(cum_var)+1))
    plt.xlabel("Number of Principal Components", fontsize=12, weight='bold')
    plt.ylabel("Cumulative Explained Variance", fontsize=12, weight='bold')
    plt.title("PCA Explained Variance", fontsize=14, weight='bold')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"05_pca_variance.png"), dpi=300)
    plt.close()

    # -------------------------------
    # 6. Species-Level Cluster Distribution
    # -------------------------------
    plot_cluster_distribution_by_species(names, labels, species_mapping, output_dir)

    # -------------------------------
    # 7. Species PCA Centroids
    # -------------------------------
    plot_species_pca_centroids(pca_data, names, labels, species_mapping, output_dir)

    # -------------------------------
    # 8. Feature Boxplots per Species
    # -------------------------------
    plot_species_feature_boxplots(features, names, species_mapping, output_dir)

    # -------------------------------
    # 9. Save results_species.csv
    # -------------------------------
    df_results = pd.DataFrame({
        "Conotoxin": names,
        "Cluster": labels,
        "Species": [species_mapping.get(n,"Unknown") for n in names],
        "PC1": pca_data[:,0],
        "PC2": pca_data[:,1],
        "Atom Count": features[:,0],
        "Radius of Gyration": features[:,1],
        "Silhouette": silhouette_avg
    })
    df_results.to_csv(os.path.join(results_dir,"results_species.csv"), index=False)
    print(f"results_species.csv saved to '{results_dir}/'")

# -------------------------
# Supporting Plot Functions
# -------------------------
def plot_cluster_distribution_by_species(names, labels, species_mapping, output_dir="plots"):
    species_list = [species_mapping.get(name,"Unknown") for name in names]
    df = pd.DataFrame({"Name":names, "Species":species_list, "Cluster":labels})
    cluster_counts = df.groupby(["Species","Cluster"]).size().unstack(fill_value=0)
    ax = cluster_counts.plot(kind="bar", stacked=True, colormap="tab10", figsize=(12,6))
    plt.ylabel("Number of Conotoxins", fontsize=12, weight='bold')
    plt.xlabel("Cone Snail Species", fontsize=12, weight='bold')
    plt.title("Structure-Based Clustering Across Species", fontsize=14, weight='bold')
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"06_cluster_by_species.png"), dpi=300)
    plt.close()

def plot_species_pca_centroids(pca_data, names, labels, species_mapping, output_dir="plots"):
    df = pd.DataFrame({
        "PC1": pca_data[:,0], 
        "PC2": pca_data[:,1],
        "Name": names, 
        "Species": [species_mapping.get(n,"Unknown") for n in names],
        "Cluster": labels
    })
    species_centroids = df.groupby("Species")[["PC1","PC2"]].mean()
    plt.figure(figsize=(10,8))
    
    # scatter all points colored by cluster
    plt.scatter(df["PC1"], df["PC2"], c=labels, cmap="tab10", s=80, alpha=0.6, edgecolor="black")
    
    # Add a dummy scatter for legend
    plt.scatter([], [], c='gray', alpha=0.6, s=80, edgecolor='black', label="Individual Conotoxins")
    
    # plot centroids with legend only
    for species, row in species_centroids.iterrows():
        plt.scatter(row.PC1, row.PC2, s=200, marker="X", label=species, edgecolor='black')
    
    plt.xlabel("PC1", fontsize=12, weight='bold')
    plt.ylabel("PC2", fontsize=12, weight='bold')
    plt.title("Species-Level PCA Centroids", fontsize=14, weight='bold')
    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"07_species_pca_centroids.png"), dpi=300)
    plt.close()

def plot_species_feature_boxplots(features, names, species_mapping, output_dir="plots"):
    df = pd.DataFrame({"Atom Count":features[:,0],
                       "Radius of Gyration":features[:,1],
                       "Species":[species_mapping.get(n,"Unknown") for n in names]})

    plt.figure(figsize=(12,5))
    sns.boxplot(x="Species", y="Atom Count", data=df, palette="Set2")
    plt.title("Atom Count Distribution Across Species", fontsize=14, weight='bold')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"08_species_atomcount.png"), dpi=300)
    plt.close()

    plt.figure(figsize=(12,5))
    sns.boxplot(x="Species", y="Radius of Gyration", data=df, palette="Set3")
    plt.title("Radius of Gyration Distribution Across Species", fontsize=14, weight='bold')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"09_species_radius.png"), dpi=300)
    plt.close()

# -------------------------
# RMSD Computation & Heatmap
# -------------------------
def compute_rmsd(pdb1_path, pdb2_path):
    parser = PDBParser(QUIET=True)
    struct1 = parser.get_structure("p1", pdb1_path)
    struct2 = parser.get_structure("p2", pdb2_path)
    coords1 = np.array([atom.coord for atom in struct1.get_atoms()])
    coords2 = np.array([atom.coord for atom in struct2.get_atoms()])
    coords1 -= coords1.mean(axis=0)
    coords2 -= coords2.mean(axis=0)
    n = min(len(coords1), len(coords2))
    rmsd = np.sqrt(np.mean(np.sum((coords1[:n]-coords2[:n])**2, axis=1)))
    return rmsd

def generate_rmsd_matrix(pdb_dir, names):
    n = len(names)
    rmsd_matrix = np.zeros((n,n))
    pdb_paths = [os.path.join(pdb_dir, f"{name}.pdb") for name in names]
    for i in range(n):
        for j in range(i, n):
            rmsd = compute_rmsd(pdb_paths[i], pdb_paths[j])
            rmsd_matrix[i,j] = rmsd_matrix[j,i] = rmsd
    return rmsd_matrix

def plot_rmsd_heatmap(rmsd_matrix, names, species_mapping, output_dir="plots"):
    os.makedirs(output_dir, exist_ok=True)
    names_sorted = sorted(names, key=lambda x: species_mapping.get(x,"Unknown"))
    idx = [names.index(n) for n in names_sorted]
    plt.figure(figsize=(12,10))
    sns.heatmap(rmsd_matrix[np.ix_(idx,idx)], xticklabels=names_sorted, yticklabels=names_sorted,
                cmap="magma", annot=True, fmt=".2f", annot_kws={"fontsize":10})
    plt.title("Pairwise Structural Similarity (RMSD)", fontsize=14, weight='bold')
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"10_rmsd_heatmap.png"), dpi=300)
    plt.close()
