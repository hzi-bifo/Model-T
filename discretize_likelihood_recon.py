"""reads in a max likelihood reconstruction matrix and discretizes the entries based on a likelihood threshold"""
import os.path
import pandas
import sys
import numpy as np

def single_matrix(m, f,  t, outdir):
    """discretize a single matrix""" 
    print "processing file %s"%(f)
    m[m>=t] = 1
    m[m<t] = 0
    s = m.shape
    m = m[~((m.iloc[:,m.shape[1] - 1] != 0) & (m.iloc[:,m.shape[1] - 1] < t))]
    print "there are %s cases with probable gain and loss event for the same edge" % (m > 1).sum().sum()
    print "there are %s probable phenotype events" % (m.iloc[:, m.shape[1] - 1] >= 1).sum()
    print "there were %s samples removed due to unclear status" % (s[0] - m.shape[0])
    m.to_csv(os.path.join(outdir, f), sep="\t", header=None, index=True)
    return m

def threshold_matrix(dir1, t, outdir, loss_dir2=None):
    """discretize one or combine two matrices into one discretized matrix"""
    if loss_dir2 is None:
        for f1 in os.listdir(dir1):
            m1 = pandas.read_csv(os.path.join(dir1, f1), sep="\t", index_col=0, header=None) 
            return single_matrix(m1, f1, t, outdir)

    if not loss_dir2 is None:
        m_set = set(os.listdir(dir1)) | set(os.listdir(loss_dir2))
        for m_f in m_set:
            if os.path.isfile(os.path.join(dir1, m_f)):
                m1 = pandas.read_csv(os.path.join(dir1, m_f), sep="\t", index_col=0, header=None) 
            else: 
                #no gain event file found, go into single matrix case
                m2 = pandas.read_csv(os.path.join(loss_dir2, m_f), sep="\t", index_col=0, header=None) 
                return single_matrix(m2, m_f, t, outdir)
                continue 
            if os.path.isfile(os.path.join(loss_dir2, m_f)):
                m2 = pandas.read_csv(os.path.join(loss_dir2, m_f), sep="\t", index_col=0, header=None) 
            else: 
                #no loss event file found, go into single matrix case
                return single_matrix(m1, m_f, t, outdir)
                continue
            print "processing file %s/%s"%(loss_dir2, m_f)
            #drop condition for m1 where the pt is non-zero but does not exceed the threshold
            drop_m1 = ~((m1.iloc[:,m1.shape[1] - 1] != 0) & (m1.iloc[:,m1.shape[1] - 1] < t)) 
            #set all elements greater or equal the threshold to 1
            m1[m1>=t] = 1
            #set all elements smaller than the threshold to 0 
            m1[m1<t] = 0
            #drop condition for m2
            drop_m2 = ~((m2.iloc[:,m2.shape[1] - 1] != 0) & (m2.iloc[:,m2.shape[1] - 1] < t)) 
            m2[m2>=t] = 1
            m2[m2<t] = 0
            m = m1 + m2
            s = m.shape
            #write to file all those samples that have been dropped due to unclear phenotype status
            #TODO not yet tested
            f = open(os.path.join(outdir, "%s_%s" %(m_f, "dropped_samples.txt")), 'w')
            for i in m[~(drop_m1 & drop_m2 | (m.iloc[:,m.shape[1] - 1] >= 1))].index:
                f.write("%s\n"%i)
            f.close()
            #combine the two conditions
            m = m[drop_m1 & drop_m2 | (m.iloc[:,m.shape[1] - 1] >= 1)]
            print "there are %s cases with probable gain and loss event for the same edge" % (m > 1).sum().sum()
            print "there are %s probable phenotype events" % (m.iloc[:, m.shape[1] - 1] >= 1).sum()
            print "there were %s samples removed due to unclear status" % (s[0] - m.shape[0])
            #for debugging purposes only
            #if any(m > 1):
                #get the elements that are ambigous
                #for i in m.index:
                #    if any(m.loc[i] > 1):
                #        print m.loc[i][m.loc[i] > 1]
                #        print m1.loc[i][m1.loc[i] > t]
                #        print m2.loc[i][m2.loc[i] > t]
                #sys.stderr.write("matrix m1 and m2 incompatible with threshold %s\n"%t)
                #raise Exception
            m.to_csv(os.path.join(outdir, m_f), sep="\t", header=None, index=True)
            return m


if __name__=="__main__":
    #test script
    import getopt
    import sys
    if len(sys.argv) == 1:
        print """USAGE: python %s
-d <input dir> with edge matrices corresponding to the likelihood of gain or loss events
-t <threshold> threshold to discretize the likelihood matrices i.e. transform them into actual reconstruction events i.e. m [ m > threshold ] = 1 ; m [ m < threshold] = 0 
-o <out dir> for the discretized matrices 
-f <optional second directory> with gain or loss likelihood matrices that shall be combined with those given in the input dir specified by the -d option
        """ % (sys.argv[0])
        sys.exit(1)
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "d:t:o:f:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    dir1 = None
    t = None
    out = None
    dir2 = None
    for o, a in optlist:
        if o == "-d":
            dir1 = a
        if o == "-t":
            t = float(a)
        if o == "-o":
            out = a 
            if os.path.exists(out):
                sys.stderr.write("output directory %s already exists; delete and rerun\n"%a)
                sys.exit(1)
            else:
                os.mkdir(out)
        if o == "-f":
            dir2 = a
    threshold_matrix(dir1, t, out, dir2)
    #threshold_matrix("/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/likelihood/gain_prob", t=0.5, outdir="/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/likelihood/gain_prob0.5/")
    #threshold_matrix("/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/likelihood/gain_prob",loss_dir2="/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/likelihood/loss_prob", t=0.5, outdir="/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/likelihood/gain_loss_prob0.5/")
