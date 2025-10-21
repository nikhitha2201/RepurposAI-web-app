def predict_targets(compound_features):
    # Placeholder ML function
    import numpy as np
    import tensorflow as tf

    # Example gene targets
    TARGET_GENES = ["TP53", "EGFR", "BRCA1", "MTOR", "VEGFA"]

    if not isinstance(compound_features, np.ndarray):
        compound_features = np.random.rand(1024)

    # small mock neural network 
    model = tf.keras.Sequential([
        tf.keras.layers.Dense(64, activation="relu", input_shape=(compound_features.shape[0],)),
        tf.keras.layers.Dense(32, activation="relu"),
        tf.keras.layers.Dense(len(TARGET_GENES), activation="softmax")
    ])

    # Select targets with probability above 0.2
    threshold = 0.2
    selected_targets = [TARGET_GENES[i] for i, p in enumerate(preds) if p > threshold]

    return selected_targets
