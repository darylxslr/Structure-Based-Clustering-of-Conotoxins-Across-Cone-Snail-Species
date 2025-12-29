from sklearn.metrics import silhouette_score

def compute_silhouette(features, labels):
    if len(set(labels)) < 2:
        return -1
    return silhouette_score(features, labels)
