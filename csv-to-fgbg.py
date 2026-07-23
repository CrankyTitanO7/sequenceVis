import pandas as pd
import numpy as np

from biology import AMINO_ALPH, AMINOS

# constants and formulae
PATH = "real_input.csv"
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


def open_csv():
    df = pd.read_csv(PATH)
    return fill_blanks(df)


def normalize(subset, val):
    return (subset / AMINOS) * val


def csv_to_quantity():
    df = open_csv()
    # Auto-detect chain positions based on the presence of numeric data
    numeric_cols = df.select_dtypes(include=[np.number]).columns

    df_normalized = df.copy()
    df_normalized[numeric_cols] = normalize(df[numeric_cols], FG_SIZE)

    return df_normalized.round(SIGFIGS)


def fg():
    quant = csv_to_quantity()
    numeric_cols = quant.select_dtypes(include=[np.number]).columns

    # Determine amino acid labels from a 'titles' column, the first column if non-numeric, or default alphabet
    if "titles" in quant.columns:
        amino_acids = quant["titles"].values
    elif quant.columns[0] not in numeric_cols:
        amino_acids = quant[quant.columns[0]].values
    else:
        amino_acids = AMINO_ALPH

    fg_dict = {}

    # Iterate over automatically detected numeric columns directly (avoids header matching)
    for col in numeric_cols:
        counts = quant[col].round().astype(int)

        col_pool = []
        for aa, count in zip(amino_acids, counts):
            col_pool.extend([aa] * count)

        # Handle rounding adjustments
        if len(col_pool) < FG_SIZE:
            col_pool.extend(
                rng.choice(amino_acids, size=FG_SIZE - len(col_pool))
            )
        elif len(col_pool) > FG_SIZE:
            col_pool = col_pool[:FG_SIZE]

        col_arr = np.array(col_pool)
        rng.shuffle(col_arr)

        fg_dict[col] = col_arr

    return pd.DataFrame(fg_dict)


def fg_generate_compiled_sequence():
    df = fg()

    # Since fg() outputs a DataFrame exclusively containing the valid position columns,
    # we can safely join characters across all columns for each row without filtering headers.
    compiled_seqs = df.astype(str).agg("".join, axis=1)

    # Insert compiled sequences into position 0 (the very first column)
    df.insert(0, "sequence", compiled_seqs)

    return df


def bg(seq_len=None):
    if seq_len is None:
        # Dynamically auto-detect chain length from the input CSV width
        df = open_csv()
        seq_len = len(df.select_dtypes(include=[np.number]).columns)

    char_matrix = rng.choice(AMINO_ALPH, size=(BG_SIZE, seq_len))

    # Override position 0 to only pick S or T
    # Position 0 is assumed to be exactly in the middle of the sequence
    pos_zero_idx = seq_len // 2
    char_matrix[:, pos_zero_idx] = rng.choice(["S", "T"], size=BG_SIZE)

    sequences = ["".join(row) for row in char_matrix]

    return pd.DataFrame({"bg": sequences})


# tests
fg_dat = fg_generate_compiled_sequence()
print("////////////// foreground generator tests //////////////")
print(fg_dat.head())

bg_dat = bg()
print("\n////////////// background generator tests //////////////")
print(bg_dat["bg"].head(10).tolist(), end="... and etc\n")
print("//////////////       tests  concluded     //////////////\n")


# export to csv
bg_dat.to_csv("bg.csv", index=False)
fg_dat.to_csv("fg_full.csv", index=False)
fg_dat.iloc[:, [0]].to_csv("fg.csv", index=False)