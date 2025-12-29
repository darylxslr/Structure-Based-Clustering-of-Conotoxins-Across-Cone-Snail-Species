from sklearn.decomposition import PCA

def run_pca(features, n_components=2):
    pca = PCA(n_components=n_components)
    reduced = pca.fit_transform(features)
    return reduced, pca
