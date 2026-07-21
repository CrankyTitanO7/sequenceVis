import pandas as pd
import numpy as np

from biology import AMINO_ALPH, AMINOS

# constants and formulae
PATH = "input.csv"
SIGFIGS = 2
SEED = 42
TRUERAND = True  # enable truly random seed

# temporary input constants
FG_SIZE = 1000
BG_SIZE = 5000
CHAIN_LENGTH = 10


def norm(x, value):
    return (x / AMINOS) * value


def open_csv():
    df = pd.read_csv(PATH)
    return df


def trim(df):
    # TODO: write trim logic, probably first column and first row
    return df


def normalize(subset, val):
    # Vectorized math is significantly faster than .map() or .apply()
    return (subset / AMINOS) * val


def csv_to_quantity():
    # BUG FIX 1: Added () to actually call open_csv()
    df = open_csv()

    # BUG FIX 2: Separate numeric columns to prevent TypeErrors on string/text columns
    numeric_cols = df.select_dtypes(include=[np.number]).columns

    df_normalized = df.copy()
    df_normalized[numeric_cols] = normalize(df[numeric_cols], FG_SIZE)

    # BUG FIX 3: Apply the SIGFIGS rounding specified in your constants/comments
    return df_normalized.round(SIGFIGS)


def bg(seq_len=CHAIN_LENGTH):  # BUG FIX 4: Use defined CHAIN_LENGTH constant
    rng = np.random.default_rng(seed=None if TRUERAND else SEED)

    char_matrix = rng.choice(AMINO_ALPH, size=(BG_SIZE, seq_len))
    sequences = ["".join(row) for row in char_matrix]

    return {"bg": sequences}


# tests
fg_dat = csv_to_quantity()
print("////////////// foreground generator tests //////////////")
print(fg_dat["A"][:3] if "P" in fg_dat.columns else fg_dat.head(3))

bg_dat = bg()
print("////////////// background generator tests //////////////")
print(bg_dat["bg"][:3])
print("//////////////       tests  concluded     //////////////")