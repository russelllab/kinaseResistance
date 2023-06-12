import fetchData
from tqdm import tqdm
import pickle
import pandas as pd
from cls import Kinase, Mutation
import argparse
import gzip, sys

import pandas as pd

data = {
    'Name': ['John', 'Emily', 'Michael'],
    'Age': ['25', '32', '19'],
    'City': ['New York', 'San Francisco', 'Chicago']
}

df = pd.DataFrame(data)

# Set the display width for each column
# pd.set_option('display.max_colwidth', max(len(value) for value in df.values.flatten()))

# Print the formatted DataFrame
print(df.to_string(index=False))

import shutil

terminal_width, terminal_height = shutil.get_terminal_size()
print("Terminal width:", terminal_width)
print("Terminal height:", terminal_height)
