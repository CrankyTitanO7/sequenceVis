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
    elif quant.columns[0].startswith("Unnamed"):
        amino_acids = quant[quant.columns[0]].values
    else:
        amino_acids = AMINO_ALPH

    # Identify position columns (-5 to +4)
    pos_cols = [col for col in quant.columns if str(col).lstrip("-+").isdigit()]
    pos_cols = sorted(pos_cols, key=lambda x: int(x))

    fg_dict = {}

    for col in pos_cols:
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

    # Strip both '-' AND '+' so all 10 position columns are selected
    pos_cols = sorted(
        [col for col in df.columns if str(col).lstrip("-+").isdigit()],
        key=lambda x: int(x),
    )

    # Join characters across all 10 columns for each row
    compiled_seqs = df[pos_cols].astype(str).agg("".join, axis=1)

    # Insert compiled sequences into position 0 (the very first column)
    df.insert(0, "sequence", compiled_seqs)

    return df


def bg(seq_len=CHAIN_LENGTH):
    char_matrix = rng.choice(AMINO_ALPH, size=(BG_SIZE, seq_len))

    # Override position 0 (index 5) to only pick S or T
    pos_zero_idx = 5
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