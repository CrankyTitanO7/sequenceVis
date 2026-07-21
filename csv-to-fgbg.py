import pandas as pd
import numpy as np
import biology as bio

from biology import AMINO_ALPH, AMINOS

# constants and formulae
PATH = "data.csv"
SIGFIGS = 2
SEED = 42
TRUERAND = True # enable truly random seed

# temporary input constants
FG_SIZE = 1000
BG_SIZE = 5000
CHAIN_LENGTH = 10
#FG_CHAIN_LENGTH = 10
#BG_CHAIN_LENGTH = 10

def norm (x, value): 
    return ((x/AMINOS) * value)

# returns a pandas dataframe
def open_csv():
    df = pd.read_csv(PATH)
    return df

# trim header and etc
def trim(): 
    # TODO: write trim logic, probably first column and first row. return trimmed data, and first column in order as a list
    return

# returns a normalized dataset from input (dataset, new sum)
def normalize(subset, val): 
    # return subset.apply(lambda x: norm(x, val))
    return subset.map(lambda x: norm(x, val))

# normalized dataset + rounded
def csv_to_quantity ():
    df = open_csv
    normal = normalize(df, FG_SIZE)
    return normal


def bg () :
    # Initialize the standard NumPy random generator
    rng = np.random.default_rng(seed=SEED) 
    if TRUERAND : rng = np.random.default_rng(seed=None) 

    num_rows = BG_SIZE
    data = {
        "bg" : rng.choice(AMINO_ALPH, size=num_rows)
    }
    return data
