# csv to fg/bg
# converts a csv matrix of csv data into foreground/background data
# csv table is position (from -5 to +4) by amino acid 
# each amino acid at each position is converted into a probability out of 20, then 
# multiplied into a larger set of size FG_SIZE
# 
# the way this works is that suppose probability 4/20 of R at position -5 (20%). then this is converted into
# value out of 1000 (200) and 200 of the  new 1000 items in the new column -5 become R 
# (in other words, "1000 choose 200" items in column -5 become R)
# 
# bg is a random generation of amino acids of length -5 to +4, inclusive


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
CHAIN_LENGTH = 10


def norm(x, value):
    return (x / AMINOS) * value

def fill_blanks(df):
    """
    Replaces missing values (NaN/None) and empty/whitespace-only string cells with 0.
    """
    # 1. Replace empty strings or space-only strings with 0
    # 2. Fill standard pandas NaN/None values with 0
    return df.replace(r'^\s*$', 0, regex=True).fillna(0)

def open_csv():
    df = pd.read_csv(PATH)
    return fill_blanks(df)


def normalize(subset, val):
    return (subset / AMINOS) * val


def csv_to_quantity():
    df = open_csv()
    numeric_cols = df.select_dtypes(include=[np.number]).columns

    df_normalized = df.copy()
    df_normalized[numeric_cols] = normalize(df[numeric_cols], FG_SIZE)

    return df_normalized.round(SIGFIGS)


def fg():
    quant = csv_to_quantity()

    # Determine amino acid labels from a 'titles' column or default alphabet
    if "titles" in quant.columns:
        amino_acids = quant["titles"].values
    else:
        amino_acids = AMINO_ALPH

    # Identify position columns (e.g., -5 to 4 or string representations '-5' to '4')
    # Strip both leading '-' and '+' signs before checking .isdigit()
    pos_cols = [col for col in quant.columns if str(col).lstrip("-+").isdigit()]

    # Ensure they are sorted numerically (-5 to +4)
    pos_cols = sorted(pos_cols, key=lambda x: int(x))

    fg_dict = {}

    for col in pos_cols:
        # Convert quantities to integer row counts
        counts = quant[col].round().astype(int)

        # Build an exact pool of amino acids matching target quantities
        col_pool = []
        for aa, count in zip(amino_acids, counts):
            col_pool.extend([aa] * count)

        # Handle any minor floating point rounding mismatches to ensure exactly FG_SIZE
        if len(col_pool) < FG_SIZE:
            col_pool.extend(
                rng.choice(amino_acids, size=FG_SIZE - len(col_pool))
            )
        elif len(col_pool) > FG_SIZE:
            col_pool = col_pool[:FG_SIZE]

        # Convert to numpy array and shuffle randomly
        col_arr = np.array(col_pool)
        rng.shuffle(col_arr)

        fg_dict[col] = col_arr

    return pd.DataFrame(fg_dict)

# foreground, but the big long chain is put all together at the first column
def fg_generate_compiled_sequence():
    df = fg()

    # Identify and sort position columns to ensure correct N-to-C terminal order (-5 to +4)
    pos_cols = sorted(
        [col for col in df.columns if str(col).lstrip("-").isdigit()],
        key=lambda x: int(x),
    )

    # Join characters across columns for each row into a single sequence string
    compiled_seqs = df[pos_cols].astype(str).agg("".join, axis=1)

    # Insert compiled sequences into position 0 (the very first column)
    df.insert(0, "sequence", compiled_seqs)

    return df


def bg(seq_len=CHAIN_LENGTH):  
    # 1. Generate random matrix for all positions
    char_matrix = rng.choice(AMINO_ALPH, size=(BG_SIZE, seq_len))
    
    # 2. Override position 0 (index 5) to only pick S or T
    pos_zero_idx = 5  
    char_matrix[:, pos_zero_idx] = rng.choice(["S", "T"], size=BG_SIZE)
    
    # 3. Join into strings
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
fg_dat.iloc[:, [0]].to_csv('fg.csv', index=False)