"""reads in a max likelihood reconstruction matrix and discretizes the entries based on a likelihood threshold"""
import os.path
import pandas as pd

def single_matrix(m, pt,  t,  outdir, discretize_pt_only = False):
    """discretize a single matrix""" 
    if discretize_pt_only:
        m.loc[m.iloc[:, m.shape[1] - 1] >= t, m.columns[m.shape[1] - 1]] = 1
    else: 
        m[m >= t] = 1
    s = m.shape
    #write to disk complete matrix
    m.to_csv(os.path.join(outdir, "pt%s.dat"%pt), sep="\t")
    #drop edges / samples that do not exceed the phenotype threshold 
    m = m[~((m.iloc[:,m.shape[1] - 1] != 0) & (m.iloc[:,m.shape[1] - 1] < t))]
    if not discretize_pt_only:
        m[m < t] = 0
        #print "there are %s cases with probable gain and loss event for the same edge" % (m > 1).sum().sum()
    #print "there are %s probable phenotype events" % (m.iloc[:, m.shape[1] - 1] >= 1).sum()
    #print "there were %s samples removed due to unclear status" % (s[0] - m.shape[0])
    m.to_csv(os.path.join(outdir, "pt%s_pruned.dat"%pt), sep="\t")
    return m

def threshold_matrix(dir1, t, outdir, pts, loss_dir2 = None, discretize_pt_only = False,  is_internal = False, are_continuous_features = False):
    """discretize one or combine two matrices into one discretized matrix"""
    if not os.path.exists(outdir):
        #if script is run in parallel the directory might have been created already
        try:
            os.mkdir(outdir)
        except OSError:
            pass
    if loss_dir2 is None:
        for pt in pts:
            m1 = pd.read_csv(os.path.join(dir1, "pt%s.dat" % pt), sep="\t", index_col=0) 
            m =  single_matrix(m1, pt, t,   outdir, discretize_pt_only = discretize_pt_only)
            if is_internal:
                return m

    if not loss_dir2 is None:
        #m_set = set(os.listdir(dir1)) | set(os.listdir(loss_dir2))
        for pt in pts:
            if os.path.isfile(os.path.join(dir1, "pt%s.dat" % pt)):
                m1 = pd.read_csv(os.path.join(dir1, "pt%s.dat" % pt), sep="\t", index_col=0) 
            else: 
                #no gain event file found, go into single matrix case
                m2 = pd.read_csv(os.path.join(loss_dir2, "pt%s.dat" % pt), sep="\t", index_col=0) 
                m =  single_matrix(m2, pt, t, outdir, discretize_pt_only = discretize_pt_only)
                if is_internal:
                    return m
                continue 
            if os.path.isfile(os.path.join(loss_dir2, "pt%s.dat" % pt)):
                m2 = pd.read_csv(os.path.join(loss_dir2, "pt%s.dat" % pt), sep="\t", index_col=0) 
            else: 
                #no loss event file found, go into single matrix case
                m = single_matrix(m1, pt, t, outdir, discretize_pt_only = discretize_pt_only)
                if is_internal:
                    return m
                continue
            #print "processing file %s/pt%s.dat"%(loss_dir2, pt)
            #drop condition for m1 where the pt is non-zero but does not exceed the threshold
            #join gains and losses
            if are_continuous_features:
                print "hier"
                m = m1
                #join phenotype probabilities by counter probabilities p1 + (1- p1) * p2
                m.loc[:, pt] = m1.loc[:, pt] + (1 - m1.loc[:, pt]) * m2 .loc[:, pt]
            #if features are probabilities join by counter probabilities p1 + (1- p1) * p2
            else:
                m = m1 + (1 - m1) * m2 
            #if features are continuous sum continuous features
            drop_m = ~((m.iloc[:,m.shape[1] - 1] != 0) & (m.iloc[:,m.shape[1] - 1] < t)) 
            #drop condition for m2
            #set all elements greater or equal the threshold to 1
            if discretize_pt_only:
                m.loc[m.iloc[:, m.shape[1] - 1] >= t, m.columns[m.shape[1] - 1]] = 1
            else:
                m[m>=t] = 1
                m[m<t] = 0
            s = m.shape
            #write to disk complete matrix
            m.to_csv(os.path.join(outdir, "pt%s.dat" % pt), sep="\t")
            #write to file all those samples that have been dropped due to unclear phenotype status
            f = open(os.path.join(outdir, "%s_%s" %(pt, "dropped_samples.txt")), 'w')
            for i in m[~(drop_m | (m.iloc[:,m.shape[1] - 1] >= 1))].index:
                f.write("%s\n"%i)
            f.close()
            #combine the two conditions
            m = m[drop_m | (m.iloc[:,m.shape[1] - 1] >= 1)]
            if discretize_pt_only:
                pass
                #print "there are %s cases with probable gain and loss event for the same edge" % (m > 1).sum().sum()
            #print "there are %s probable phenotype events" % (m.iloc[:, m.shape[1] - 1] >= 1).sum()
            #print "there were %s samples removed due to unclear status" % (s[0] - m.shape[0])
            #TODO only write to file if the call is external i.e. in the bulk discretization run
            m.to_csv(os.path.join(outdir, "pt%s_pruned.dat" % pt), sep="\t")
            if is_internal:
                return m
