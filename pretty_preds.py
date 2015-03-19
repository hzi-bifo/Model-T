PFAM_PTS_F = "../stol_2_NCBI20140115_candidatus/pfam_pts_names_nl.txt"
import pandas as ps
def flatten_df(df1, df2, out):
    with open(out, 'w') as f:
        for i in set(df1.index).union(set(df2.index)):
            for j in set(df1.columns).union(set(df2.columns)):
                if i in df1.index and j in df1.columns and  not ps.isnull(df1.loc[i,j]):
                    f.write("\t".join([i, j, str(df1.loc[i,j]), "phypats plus max likelihood reconstruction"]) + "\n")
                if i in df2.index and j in df2.columns and not ps.isnull(df2.loc[i,j]):
                    f.write("\t".join([i, j, str(df2.loc[i,j]), "phypats"]) + "\n")

def flatten_preds(phypat_f, ml_plus_phypat_f, sample_names_f, out_f):
    #read in draft genomes phenotype predictions 
    ml_plus_phypat = ps.read_csv(ml_plus_phypat_f, sep = "\t", index_col = 0)
    ml_plus_phypat.columns = ml_plus_phypat.columns.astype('int') + 1
    phypats = ps.read_csv(phypat_f, sep = "\t", index_col = 0)
    phypats.columns = phypats.columns.astype('int') + 207
    bin_names = []
    #read in names of bins
    with open(sample_names_f) as f:
        for l in f:
            bin_names.append(l.split("/")[-1].split(".")[0].strip())
    #read in names of phenotypes
    pt_names = ps.read_csv(PFAM_PTS_F, sep = "\t", header = None, index_col = 0)
    ml_plus_phypat.columns = pt_names.loc[ml_plus_phypat.columns, 1]
    phypats.columns = pt_names.loc[phypats.columns, 1 ]
    ml_plus_phypat.index = bin_names
    phypats.index = bin_names
    flatten_df(ml_plus_phypat, phypats, out_f)



    

if __name__ == '__main__':
    import getopt
    import sys
    if len(sys.argv) == 1:
        print """USAGE: python %s
-p <phypat_preds> a file with the predictions generated with phypats only models 
-m <ml+phypat_preds> a file with the predictions from the phypat plus ml models 
-s <sample name file> containing the name of the samples   
-o <out file> for the pretty predictions
""" % (sys.argv[0])
        sys.exit(1)
    optlist, args = getopt.getopt(sys.argv[1:], "p:m:s:o:")
    phypat_f = ml_plus_phypat_f = sample_names_f = out_f = None
    try:
        for o, a in optlist:
            if o == "-p":
                phypat_f = a 
            if o == "-m":
                ml_plus_phypat_f = a
            if o == "-s":
                sample_names_f = a
            if o == "-o":
                out_f = a
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    flatten_preds(phypat_f, ml_plus_phypat_f, sample_names_f, out_f)
    

                 
