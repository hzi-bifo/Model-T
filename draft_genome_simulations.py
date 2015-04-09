"""script to simulate genomes with incomplete gene sets from fully sequenced genomes"""
import predict
import domtblout2gene
import random 
import math
import pandas as ps
import sys
import tarfile
import nested_cv
from scipy.stats import nanmean


UID2FN = "/net/metagenomics/data/genomes/public/NCBI20140115_FAAs_uid2fn.txt"
class SimulateDraftGenome:

    def __init__(self, sample_points, no_samples, refseq_annot_fs, refseq_f, pfam_acc_f, model_dir, train_pfam_f, sp2str_f):
        self.sample_points = sample_points
        self.no_samples = no_samples
        self.model_dir = model_dir
        self.train_pfam_f = train_pfam_f
        fs = [i.strip() for i in open(refseq_annot_fs, 'r').readlines()]
        annot, self.gene2hmm, self.rs2f = domtblout2gene.gene2hmm(fs, pfam_acc_f)
        #read the species to strain mapping and infer a species2gene dictionary  
        self.sp2gene = self.sp2str(sp2str_f)
        #print self.sp2gene.keys()
        #print annot.sum()
        #print self.rs2f
        #self.rs2gene = ps.read_csv(gene_length_f, sep = "\t") 
        #self.rs2gene.index = ["RefSeqProId+%s" % i for i in rs2gene.index]

    def sp2str(self, sp2str_f):
        """read in species to strain mapping and infer a species2gene dictonary analogously to the gene2hmm dictionary"""
        sp2str = ps.read_csv(sp2str_f, sep = "\t", index_col = 3, header = None) 
        sp2gene = {}
        for sp in sp2str.index:
            if type(sp2str.loc[sp,]) == ps.Series:
                sp2gene[sp] = self.gene2hmm["RefSeqProId+%s" % sp2str.loc[sp, 2]]
            else:
                for rs in sp2str.loc[sp, 2]:
                    if sp not in sp2gene:
                        sp2gene[sp] = self.gene2hmm["RefSeqProId+%s" % rs]
                    else: 
                        sp2gene[sp] = dict(sp2gene[sp].items() + self.gene2hmm["RefSeqProId+%s" % rs].items())
        return sp2gene
            

    def simulate_draft_genomes(self):
        """create for each fully sequenced genome, x number of samples for y number of sample points"""
        #read in pfam accs
        pfams = ps.read_csv(pfam_acc_f, index_col = 0, header = None).index
        #sample Pfam matrices at different samples points
        full_m = ps.DataFrame(ps.np.zeros(shape = (len(self.sp2gene) * len(self.sample_points) * self.no_samples, len(pfams))))
        full_m.columns = pfams
        #set row names
        rns = ["%s_%s_%s" % (self.sp2gene.keys()[j], i, k) for j in range(len(self.sp2gene)) for i in range(self.no_samples) for k in self.sample_points]
        full_m.index = rns
        for j in range(len(self.sp2gene)):
            print "sample", self.sp2gene.keys()[j]
            for i in range(self.no_samples):
                #reverse the list of sample points
                sample_points_t = self.sample_points[::-1]
                #create a matrix to store the simulated draft genomes 
                m = ps.DataFrame(ps.np.zeros(shape = (len(self.sample_points), len(pfams) )))
                #rows for the different samples points
                m.index = self.sample_points
                #columns for the different pfam families
                m.columns = pfams
                #start sampling
                cur_sample_point = sample_points_t.pop()
                #create a vector to store the genome annotations
                cur_sample_series = m.iloc[0, ].copy()
                #draw a samples of genes with the same size 
                genes = random.sample(self.sp2gene[self.sp2gene.keys()[j]], len(self.sp2gene[self.sp2gene.keys()[j]])) 
                print "total number of genes", len(genes)
                #draw random sample
                
                #retrieve vector of pfams 
                for k in range(len(genes)):
                    #check if the current sampling point is complete
                    if k + 1 >= cur_sample_point * len(genes):
                        #print i, sample_points_t, cur_sample_point
                        #if so assign that vector to the matrix
                        print "sampling point at", cur_sample_point, k
                        m.loc[cur_sample_point, ] = cur_sample_series
                        #stop when all sampling points have been processed
                        if sample_points_t == []:
                            break
                        cur_sample_point = sample_points_t.pop()
                    for pf in self.sp2gene[self.sp2gene.keys()[j]][genes[k]]:
                        cur_sample_series.loc[pf] += 1
                #print j * len(self.sample_points) * self.no_samples + len(self.sample_points) * i
                #print (j * len(self.sample_points) * self.no_samples + len(self.sample_points) * i + len(self.sample_points))
                #print m.shape
                #print full_m.shape
                full_m.iloc[(j * len(self.sample_points) * self.no_samples + len(self.sample_points) * i) : (j * len(self.sample_points) * self.no_samples + len(self.sample_points) * i + len(self.sample_points)), ] = m.values
                #if self.gene2hmm.keys()[j] not in pf_sampled:
                #    pf_sampled[self.gene2hmm.keys()[j]] = [m]
                #else: 
                #    pf_sampled[self.gene2hmm.keys()[j]].append(m)
        #ps.set_option("display.max_columns",101)
        return full_m

   
    def evaluate_draft_genome_preds(self, preds, pt_gold_f, pt2id_f, outdir):
        """read in gold standard pt predictions"""
        #macro accuracy evaluation per sampling point averaged over all phenotypes; this might be skewed due to very few test samples for some phenotypes
        #evaluation per sampling point per phenotype -> influence on specific phenotypes
        #micro accuracy pooled evaluation of all predictions
        #matrix: 
        #Rf1s1 10%
        #Rf1s1 20%
        #Rf1s2 10%
        #Rf1s2 20%
        #Rf2s1 10%
        pt_gold_df = ps.read_csv(pt_gold_f, sep = "\t", index_col = 0)
        #print pt_gold_df.columns
        #read in pt2id file
        pt2id = ps.read_csv(pt2id_f, index_col = 0, sep = "\t", header = None) 
        pt_gold_df.columns = pt2id.index
        #restrict to pt that were selected 
        pt_gold_sel = pt_gold_df.loc[:, preds.columns]
        #print pt_gold_sel
        #restrict to pt that have a model available and restrict to the species in the species2gene dict
        pt_gold_sel_model = pt_gold_sel.loc[self.sp2gene.keys(), preds.columns]
        #get boolean matrix that indices the non-missing values
        pt_gold_missing_m = ~ps.isnull(pt_gold_sel.loc[self.sp2gene.keys(), preds.columns])
        #get phenotypes that are completely absent in the gold standard
        pt_gold_absent = (~pt_gold_missing_m).apply(ps.np.all)
        #print "pt_gold_absent", pt_gold_absent
        #if not pt_gold_missing_m.all().all():
        #    print "no labels found for the input phenotypes in the gold standard"
        #    sys.exit(1)
        #unfold array and index missing phenotypes
        perf_sampling_points_pt = ps.DataFrame(ps.np.zeros(shape = (len(sampling_points), 4 * pt_gold_sel.shape[1])))
        perf_sampling_points_pt.index = sampling_points
        #perf_sampling_points.columns = ["recall positive class", "recall negative class", "balanced accuracy", "precision"]
        for i in range(len(sampling_points)):
            print "processing sampling point", sampling_points[i]
            #get confusion matrix at each sampling point i.e. 10%, 20%, ..., 100%
            #matrix to store the confusion matrices for each phenotype and sub sample
            perf_sample_pt  = ps.DataFrame(ps.np.zeros(shape = (no_samples, 4 * pt_gold_sel.shape[1])))
            for j in range(no_samples):
                #get indices  for each species at the same sampling point for the same sub sample
                loc_pred = [i + len(sampling_points) * j + len(sampling_points) * no_samples * k for k in range(len(self.sp2gene))] 
                labels = pt_gold_df 
                #print "sampling_point", sampling_points[i], "sample", j
                #print "loc_pred", loc_pred
                #get subset of predictions that correspond to the current sub sample and sampling point
                preds_ij = preds.iloc[loc_pred, ]
                preds_ij.index = pt_gold_sel_model.index
                #print "predictions", preds_ij
                #print "pt_gold_missing row", pt_gold_missing_m
                #retrieve confusion matrix for every phenotype
                for pt_i in range(preds_ij.shape[1]):
                        perf_sample_pt.iloc[j, pt_i * 4 : (pt_i * 4 + 4)] += nested_cv.nested_cv.confusion_m(pt_gold_sel_model.iloc[:, pt_i][pt_gold_missing_m.iloc[:, pt_i].values], preds_ij.iloc[:, pt_i][pt_gold_missing_m.iloc[:, pt_i].values])
                        #print "phenotype", pt_i, "confusion m", nested_cv.nested_cv.confusion_m(pt_gold_sel_model.iloc[:, pt_i][pt_gold_missing_m.iloc[:, pt_i].values], preds_ij.iloc[:, i][pt_gold_missing_m.iloc[:, pt_i].values])
            #sum up samples per phenotype
            perf_sampling_points_pt.iloc[i, ] = perf_sample_pt.sum(axis = 0).values
            #print "perf_sample_pt sample point", i, perf_sample_pt.sum(axis = 1)
        #print "perf_sampling_points_pt", perf_sampling_points_pt
        #micro accuracy
        conf_per_sp = ps.DataFrame(ps.np.zeros(shape = (len(sampling_points), 4)))
        conf_per_sp.index = sampling_points
        conf_per_sp.columns = ["TN", "FP", "FN", "TP"]
        #store performance measures per sampling point per phenotype
        perf_per_pt = ps.np.zeros(shape = (len(sampling_points), pt_gold_sel_model.shape[1], 4))
        for sp in range(len(sampling_points)):
            for pt_i in range(pt_gold_sel_model.shape[1]):
                conf_pt = perf_sampling_points_pt.iloc[sp, pt_i * 4 : (pt_i * 4 + 4)]
                #print "conf_pt", conf_pt
                conf_per_sp.iloc[sp, :] += conf_pt.values
                rec_pos = nested_cv.nested_cv.recall_pos_conf(conf_pt)
                rec_neg = nested_cv.nested_cv.recall_neg_conf(conf_pt)
                prec = nested_cv.nested_cv.precision_conf(conf_pt)
                bacc = nested_cv.nested_cv.bacc(rec_pos, rec_neg)
                #print rec_pos, rec_neg, bacc, prec
                perf_per_pt[sp, pt_i, :] = ps.np.array([rec_pos, rec_neg, bacc, prec])
        #print "conf_per_sp", conf_per_sp
        #print "perf_per_pt", perf_per_pt
        perf_colnames = ["recall positive class", "recall negative class", "balanced accuracy", "precision"]
        #get macro measures 
        macro_per_sp = ps.DataFrame(nanmean(perf_per_pt[:, (~pt_gold_absent).values, : ], axis = 1))
        macro_per_sp.index = sampling_points
        macro_per_sp.columns = perf_colnames 
        #get micro measures 
        micro_per_sp = conf_per_sp.apply(lambda x: ps.Series([nested_cv.nested_cv.recall_pos_conf(x), nested_cv.nested_cv.recall_neg_conf(x), nested_cv.nested_cv.bacc(nested_cv.nested_cv.recall_pos_conf(x), nested_cv.nested_cv.recall_neg_conf(x)), nested_cv.nested_cv.precision_conf(x)]), axis = 1)
        micro_per_sp.columns = perf_colnames 
        #print "micro", micro_per_sp
        #write results to disk
        micro_per_sp.to_csv("%s/micro.tsv"%outdir, sep = "\t")
        macro_per_sp.to_csv("%s/macro.tsv"%outdir, sep = "\t")
        #write individual phenotype performances to disk per sampling point
        for i in range(len(sampling_points)):
            df = ps.DataFrame(perf_per_pt[i, :, :])
            df.index = pt_gold_sel_model.columns
            df.columns = perf_colnames
            df[~pt_gold_absent].to_csv("%s/perf_sp%s.tsv"%(outdir, sampling_points[i]), sep = "\t")


         
                

    def predict_draft_genomes(self,data, pt1, pt2, k, mask_pts = None):
        """for all the phenotypes and all simulated genomes get phenotype predictions"""
        self.pt_has_model = dict((i, False) for i in range(pt1, pt2 + 1)) 
        #restrict the data matrix to the pfams from the training data set
        data_train_pf = data.loc[:, ps.read_csv(self.train_pfam_f, index_col = 0, header = None).index]
        #instaniate a matrix for the individual predictions
        preds = ps.DataFrame(ps.np.zeros(shape = (len(self.sp2gene) * len(self.sample_points) * self.no_samples ,pt2 + 1 - pt1)))
        preds.columns = range(pt1, pt2 + 1)
        for pt in range(pt1, pt2 + 1):
            maj_pred = predict.majority_predict(pt, self.model_dir, data_train_pf, k)
            maj_preds_aggr = maj_pred.apply(lambda x: 1 if (x > 0).sum() > k/2  else -1, axis = 1)
            #if there is a model available fill in the phenotype prediction into the matrix
            if not len(maj_preds_aggr) == 0:
                self.pt_has_model[pt] = True
                preds.loc[:, pt] = maj_preds_aggr 
        if not mask_pts is None:
            for pt in mask_pts:
                self.pt_has_model[pt] = False
        return preds.loc[:, [self.pt_has_model[i] for i in range(pt1, pt2 + 1)]]

    def get_genes_targz(self, gff_targz):
        rs2gene = {}
        #read refseq id (uid) to directory mapping
        uid2fn = ps.read_csv(UID2FN, sep = "\t", index_col = 0, header = None) 
        uid2fn.index = ["RefSeqProId+%s" % i for i in uid2fn.index]
        print self.sp2gene
        with tarfile.open(gff_targz, mode = 'r:gz') as tf:
            for rs in self.gene2hmm:
                for ce in self.rs2f[rs]:
                    print "extracted", uid2fn
                    tar = tf.getmember("./%s/%s" % (uid2fn.loc[rs, ].iloc[0], ce))
                    rs2gene[refseq] = ps.read_csv(tf.extractfile(tar), sep = "\t", index_col = 0, commentchar = '#').index
        return rs2gene

 
if __name__ == '__main__':
    import getopt
    if len(sys.argv) == 1:
        print """USAGE: python %s
-s <number of samples> number of samples to be drawn for each genome
-p <[sampling point 1, sampling point 2, ...]> sequence of ordered relative number of genes to be retained e.g. -p 0.1,0.2,0.5,1.0 
-g Pfam (gene/protein) family accession files
-f <gene counts> which contains a bunch of gff file for every species
-m <model dir> directory with the liblinear models 
-a <file> with refseq chromosomal annotation files 
-r <file> with refseq ids 
-b <phenotype_range> e.g. 5000:5010
-c <gold standard matrix> with the phenotypes
-d <pfams> that occur in the training set 
-i <mapping file> that maps phenotypes to the internally used ids
-e <mapping file> that maps RefSeqProIds to species
-o <outdir> for macro, micro accuracy etc. 
""" % (sys.argv[0])
        sys.exit(1)
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "s:p:g:m:r:a:b:c:d:i:e:o:k:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    no_samples = sampling_points = gff_targz = pfam_acc_f = model_dir = refseq_annot_fs = refseq_f = tr_pfams_f = pt_gold_f = pt2id_f = sp2str_f = outdir = mask_pts = None
    for o, a in optlist:
        if o == "-b":
            pt1, pt2 = [int(i) for i in a.split("-")]
        if o == "-s":
            no_samples = int(a)
        if o == "-p":
            sampling_points = [float(i) for i in a.split(",")] 
        if o == "-g":
            pfam_acc_f = a
        #if o == "-f":
        #    gene_length_f = a
        if o == "-m":
            model_dir = a
        if o == "-a":
            refseq_annot_fs = a
        if o == "-r":
            refseq_f = a
        if o == "-c":
            pt_gold_f = a
        if o == "-d":
            tr_pfams_f = a
        if o == "-i":
            pt2id_f = a
        if o == "-e":
            sp2str_f = a
        if o == "-k":
            mask_pts = [int(i) for i in a.split(",")]
        if o == "-o":
            outdir = a
    import os
    if os.path.exists(outdir):
        sys.stderr.write("outdir %s already exists"%outdir)
        sys.exit(1)
    else:
        os.mkdir(outdir)

    simulator = SimulateDraftGenome(sampling_points, no_samples, refseq_annot_fs, refseq_f, pfam_acc_f,  model_dir, tr_pfams_f, sp2str_f)
    data = simulator.simulate_draft_genomes()
    preds = simulator.predict_draft_genomes(data, pt1, pt2, 5, mask_pts)
    simulator.evaluate_draft_genome_preds(preds, pt_gold_f, pt2id_f, outdir)
