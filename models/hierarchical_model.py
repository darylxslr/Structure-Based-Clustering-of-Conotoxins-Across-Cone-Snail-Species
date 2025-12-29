from scipy.cluster.hierarchy import linkage

def hierarchical_clustering(features):
    """
    Perform hierarchical clustering using Ward linkage
    """
    linkage_matrix = linkage(features, method="ward")
    return linkage_matrix
