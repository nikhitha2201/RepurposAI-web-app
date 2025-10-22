import tensorflow as tf
import numpy as np
from typing import List, Dict

def build_model(input_dim: int, n_classes: int) -> tf.keras.Model:
    model = tf.keras.Sequential([
        tf.keras.layers.Input(shape=(input_dim,)),
        tf.keras.layers.Dense(1024, activation="relu"),
        tf.keras.layers.Dropout(0.3),
        tf.keras.layers.Dense(512, activation="relu"),
        tf.keras.layers.Dropout(0.2),
        tf.keras.layers.Dense(n_classes, activation="softmax")
    ])
    
    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3),
        loss="sparse_categorical_crossentropy",
        metrics=[
            "accuracy",
            tf.keras.metrics.SparseTopKCategoricalAccuracy(k=3, name="top_3_acc")
        ]
    )
    
    return model

def predict_targets(model: tf.keras.Model, compound_features: np.ndarray, label_map: Dict[int, str]) -> List[Dict[str, float]]:
    probs = model.predict(compound_features)  # shape: (n_samples, n_classes)
    all_predictions = []
    for prob_vec in probs:
        pred_dict = {label_map[i]: float(prob) for i, prob in enumerate(prob_vec)}
        all_predictions.append(pred_dict)
    return all_predictions
