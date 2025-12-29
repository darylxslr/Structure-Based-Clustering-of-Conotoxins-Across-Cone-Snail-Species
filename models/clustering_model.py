from sklearn.cluster import KMeans
import numpy as np

def kmeans_clustering(features, n_clusters=4):
    model = KMeans(
        n_clusters=n_clusters,
        random_state=42,
        n_init=20
    )
    labels = model.fit_predict(features)

    # SAFETY CHECK: avoid single-cluster collapse
    if len(np.unique(labels)) == 1:
        model = KMeans(n_clusters=2, random_state=42, n_init=20)
        labels = model.fit_predict(features)

    return labels, model
