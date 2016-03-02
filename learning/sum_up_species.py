import argparse
import pandas as pd

def run(args):
    # parse input
    df = pd.read_csv(args.input_file, sep="\t")
    # find indizes of same species
    #same_species = find_same_species(df)
    # call method
    if args.method == "UN":
        df = calculate_union(df)
    else:
        df = calculate_majority(df)
    # write to output file
    df.to_csv(args.output, sep="\t")

def find_same_species(df):
    same_species = {}
    for idx, sp1 in df[:-1].iterrows():
        for id, sp2 in df[idx+1:].iterrows():
            if sp1[0] == sp2[0]:
                try:
                    same_species[sp1[0]].append(id)
                except KeyError:
                    same_species[sp1[0]] = [idx,id]
    return same_species

def calculate_union(df):
    df = df.groupby(df.columns[0])
    df = df.sum()
    return df

def calculate_majority(df):
    df = df.groupby(df.columns[0])
    df = df.max()
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='input file name')
    parser.add_argument('-o', '--output', default= 'corrected_df.txt',
                       help='choose putput file name. Default: corrected_df.txt')
    parser.add_argument('-m', '--method', default= 'UN', choices=['UN', 'MV'],
                       help='choose method to unite results. Union or Majority Vote')
    args = parser.parse_args()
    run(args)
