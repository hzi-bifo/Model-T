import joblib
import pandas as ps
def extract_bias(pickled_f, out_f):
    l = joblib.load(pickled_f)
    bias_cparams = []
    for m in l:
        bias_cparams.append((m.C, m.intercept_[0]))
    ps.DataFrame(bias_cparams).to_csv(out_f, sep = "\t", index = None, header = None)

import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser("extract the bias term from pickeld SVM models")
    parser.add_argument("pickled_f", help = "pickled model file")
    parser.add_argument("out_f", help = "out file for the bias term")
    args = parser.parse_args()
    extract_bias(args.pickled_f, args.out_f) 
