
"""reads in a max likelihood reconstruction matrix and discretizes the entries based on a likelihood threshold"""
import os.path
import pandas
import sys
import numpy as np

def single_matrix(m, f, outdir):
    """discretize a single matrix""" 
    print "processing file %s"%(f)
    m.to_csv(os.path.join(outdir, f), sep="\t", header=None, index=True)
    return m

def threshold_matrix(dir1, outdir, loss_dir2=None, is_internal = False):
    """discretize one or combine two matrices into one discretized matrix"""
    if loss_dir2 is None:
        for f1 in os.listdir(dir1):
            m1 = pandas.read_csv(os.path.join(dir1, f1), sep="\t", index_col=0, header=None) 
            m =  single_matrix(m1, f1, outdir)
            if is_internal:
                return m

    if not loss_dir2 is None:
        m_set = set(os.listdir(dir1)) | set(os.listdir(loss_dir2))
        for m_f in m_set:
            if os.path.isfile(os.path.join(dir1, m_f)):
                m1 = pandas.read_csv(os.path.join(dir1, m_f), sep="\t", index_col=0, header=None) 
            else: 
                #no gain event file found, go into single matrix case
                m2 = pandas.read_csv(os.path.join(loss_dir2, m_f), sep="\t", index_col=0, header=None) 
                m =  single_matrix(m2, m_f, outdir)
                if is_internal:
                    return m
                continue 
            if os.path.isfile(os.path.join(loss_dir2, m_f)):
                m2 = pandas.read_csv(os.path.join(loss_dir2, m_f), sep="\t", index_col=0, header=None) 
            else: 
                #no loss event file found, go into single matrix case
                m = single_matrix(m1, m_f, outdir)
                if is_internal:
                    return m
                continue
            print "processing file %s/%s"%(loss_dir2, m_f)
            #combine gain and loss events
            m = m1 + (1 - m1) * m2 
            if is_internal:
                return m
            else:
                m.to_csv(os.path.join(outdir, m_f), sep="\t", header=None, index=True)


if __name__=="__main__":
    #test script
    import getopt
    import sys
    if len(sys.argv) == 1:
        print """USAGE: python %s
-d <input dir> with edge matrices corresponding to the likelihood of gain or loss events
-f <second directory> with gain or loss likelihood matrices that shall be combined with those given in the input dir specified by the -d option
-o <out dir> for the joined gain loss matrices 
        """ % (sys.argv[0])
        sys.exit(1)
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "d:o:f:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    dir1 = None
    t = None
    out = None
    dir2 = None
    discretize_pt_only = False 
    for o, a in optlist:
        if o == "-d":
            dir1 = a
        if o == "-o":
            out = a 
            if os.path.exists(out):
                sys.stderr.write("output directory %s already exists; delete and rerun\n"%a)
                sys.exit(1)
            else:
                os.mkdir(out)
        if o == "-f":
            dir2 = a
    threshold_matrix(dir1,  out, dir2)
