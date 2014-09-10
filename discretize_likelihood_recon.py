"""reads in a max likelihood reconstruction matrix and discretizes the entries based on a likelihood threshold"""
import os.path
import pandas
import sys
def threshold_matrix(dir1, t, outdir, loss_dir2=None):
    for f1 in os.listdir(dir1):
        m1 = pandas.read_csv(os.path.join(dir1, f1), sep="\t", index_col=0, header=None) 
        if not loss_dir2 is None:
            if os.path.isfile(os.path.join(loss_dir2, f1)):
                m2 = pandas.read_csv(os.path.join(loss_dir2, f1), sep="\t", index_col=0, header=None) 
                m1[m1>=t] = 1
                m1[m1<t] = 0
                m2[m2>=t] = 1
                m2[m2<t] = 0
                m = m1 + m2
                print m.shape
                print m[m>1]
                if any(m > 1):
                    sys.stderr.write("matrix m1 and m2 incompatible with threshold %s\n"%t)
                    raise Exception
                m.to_csv(os.path.join(outdir, f1), sep="\t", header=None, index=True)
            else: 
                raise Exception
        else:
            m1[m1>=t] = 1
            m1[m1<t] = 0
            m1.to_csv(os.path.join(outdir, f1), sep="\t", header=None, index=True)


if __name__=="__main__":
    #test script
    #threshold_matrix("/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/likelihood/gain_prob", t=0.5, outdir="/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/likelihood/gain_prob0.5/")
    threshold_matrix("/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/likelihood/gain_prob",loss_dir2="/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/likelihood/loss_prob", t=0.5, outdir="/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/likelihood/gain_loss_prob0.5/")
