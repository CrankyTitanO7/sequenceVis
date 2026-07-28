import pandas as pd
import numpy as np

from biology import AMINO_ALPH, AMINOS

# constants and formulae
PATH = "real_input.csv"
BGEX_PATH = "bginput.csv"
SIGFIGS = 2
SEED = 42
TRUERAND = True  # enable truly random seed

rng = np.random.default_rng(seed=None if TRUERAND else SEED)

# temporary input constants
FG_SIZE = 1000
BG_SIZE = 5000


def norm(x, value):
    return (x / AMINOS) * value


def fill_blanks(df):
    """
    Replaces missing values (NaN/None) and empty/whitespace-only string cells with 0.
    """
    return df.replace(r'^\s*$', 0, regex=True).fillna(0)


def open_csv(path):
    df = pd.read_csv(path)
    return fill_blanks(df)


def normalize(subset, val):
    return (subset / AMINOS) * val


def csv_to_quantity(path, length):
    df = open_csv(path)
    # Auto-detect chain positions based on the presence of numeric data
    numeric_cols = df.select_dtypes(include=[np.number]).columns

    df_normalized = df.copy()
    df_normalized[numeric_cols] = normalize(df[numeric_cols], length)

    return df_normalized.round(SIGFIGS)


def gen(path, length):
    quant = csv_to_quantity(path, length)
    numeric_cols = quant.select_dtypes(include=[np.number]).columns

    # Determine amino acid labels from a 'titles' column, the first column if non-numeric, or default alphabet
    if "titles" in quant.columns:
        amino_acids = quant["titles"].values
    elif quant.columns[0] not in numeric_cols:
        amino_acids = quant[quant.columns[0]].values
    else:
        amino_acids = AMINO_ALPH

    gen_dict = {}

    # Iterate over automatically detected numeric columns directly (avoids header matching)
    for col in numeric_cols:
        counts = quant[col].round().astype(int)

        col_pool = []
        for aa, count in zip(amino_acids, counts):
            col_pool.extend([aa] * count)

        # Handle rounding adjustments
        if len(col_pool) < length:
            col_pool.extend(
                rng.choice(amino_acids, size=length - len(col_pool))
            )
        elif len(col_pool) > length:
            col_pool = col_pool[:length]

        col_arr = np.array(col_pool)
        rng.shuffle(col_arr)

        gen_dict[col] = col_arr

    return pd.DataFrame(gen_dict)


def generate_compiled_sequence(path, length):
    df = gen(path, length)

    # Since fg() outputs a DataFrame exclusively containing the valid position columns,
    # we can safely join characters across all columns for each row without filtering headers.
    compiled_seqs = df.astype(str).agg("".join, axis=1)

    # Insert compiled sequences into position 0 (the very first column)
    df.insert(0, "sequence", compiled_seqs)

    return df

def fg():
    return generate_compiled_sequence(PATH, FG_SIZE)
def bg():
    return generate_compiled_sequence(BGEX_PATH, BG_SIZE)


# tests
fg_dat = fg()
print("////////////// foreground generator tests //////////////")
print(fg_dat.head())

bg_dat = bg()
print("\n////////////// background generator tests //////////////")
print(bg_dat.head())
print("//////////////       tests  concluded     //////////////\n")


# export to csv
bg_dat.to_csv("bg_full.csv", index=False)
fg_dat.to_csv("fg_full.csv", index=False)
bg_dat.iloc[:, [0]].to_csv("bg.csv", index=False)
fg_dat.iloc[:, [0]].to_csv("fg.csv", index=False)