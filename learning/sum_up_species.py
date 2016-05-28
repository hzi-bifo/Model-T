import argparse
import pandas as pd

def run(args):
    # parse input
    df = pd.read_csv(args.input_file, sep="\t")
    df ["Refseq IDs"] = df["Unnamed: 0"]
    df = df.set_index('Unnamed: 0')
    df.index = df.index.values.astype('string')

    # create dictionary mapping refseq to species IDs
    ref2sp = create_ref2sp_dict(args.mapping_file)
    print ref2sp
    # map refseq ids to sp ids
    df = map_ref2sp(ref2sp, df)
    print df
    # call method
    if args.method == "UN":
        df = calculate_union(df)
    else:
        df = calculate_majority(df)
    # replace ids again
    #df = df.set_index("Refseq IDs")
    # replace RefSeq ID column
    del df["Refseq IDs"]
    # write to output file
    df.loc[list(set(ref2sp.values())), :].to_csv(args.output, sep="\t")

def create_ref2sp_dict(map):
    ref2sp = {}
    with open(map, "r") as mapfile:
        for line in mapfile:
            tmp = line.split("\t")
            sp = tmp[0]
            for el in tmp[1].replace("\n", "").split(","):
                ref2sp[el.strip()] = sp.strip()
    return ref2sp

def map_ref2sp(ref2sp, df):
    df.rename(ref2sp, inplace=True)
    return df

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
    df = df.groupby(df.index)
    df = df.sum()
    return df

def calculate_majority(df):
    df = df.groupby(df.index)
    df = df.max()
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='input file name')
    parser.add_argument('mapping_file', help='mapping file name that maps refseq IDs to species IDs')
    parser.add_argument('-o', '--output', default= 'corrected_df.txt',
                       help='choose putput file name. Default: corrected_df.txt')
    parser.add_argument('-m', '--method', default= 'UN', choices=['UN', 'MV'],
                       help='choose method to unite results. Union or Majority Vote')
    args = parser.parse_args()
    run(args)
