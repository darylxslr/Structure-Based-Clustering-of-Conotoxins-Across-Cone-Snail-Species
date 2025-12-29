from sklearn.preprocessing import StandardScaler

def scale_features(features):
    scaler = StandardScaler()
    return scaler.fit_transform(features)
