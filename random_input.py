import pandas as pd
import numpy as np
import biology as bio

cols = list(range(-5, 5))  # -5 through 3

rng = np.random.default_rng(seed=None)

# Generate a 20x9 matrix all at once
raw = rng.uniform(1, 100, size=(len(bio.AMINO_ALPH), len(cols)))
normalized = (raw / raw.sum(axis=0)) * 20

# Create DataFrame
df = pd.DataFrame(normalized, columns=cols)
df.insert(0, 'titles', bio.AMINO_ALPH)

df.to_csv('input.csv', index=False)