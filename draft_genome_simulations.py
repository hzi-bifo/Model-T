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

    def __init__(self, sample_points, no_samples, refseq_annot_fs, refseq_f, pfam_acc_f, model_tar, sp2str_f, min_samples):
        self.sample_points = sample_points
        self.no_samples = no_samples
        self.model_tar = model_tar
        #minimum number of samples per phenotype + and - class 
        self.min_samples = min_samples
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
            

    def simulate_draft_genomes(self, pfam_acc_f):
        """create for each fully sequenced genome, x number of samples for y number of sample points"""
        #read in pfam accs
        pfams = ps.read_csv(pfam_acc_f, header=None, index_col = 0, sep = "\t").index
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
                #print "total number of genes", len(genes)
                #draw random sample
                
                #retrieve vector of pfams 
                for k in range(len(genes)):
                    #check if the current sampling point is complete
                    if k + 1 >= cur_sample_point * len(genes):
                        #print i, sample_points_t, cur_sample_point
                        #if so assign that vector to the matrix
                        #print "sampling point at", cur_sample_point, k
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

    def predict_draft_genomes(self,data, pt1, pt2, k = 5, mask_pts = None, adjust_bias = False):
        """for all the phenotypes and all simulated genomes get phenotype predictions"""
        self.pt_has_model = dict((i, False) for i in range(pt1, pt2 + 1)) 
        #restrict the data matrix to the pfams from the training data set
        pfam_pts_mapping = ps.read_csv(self.model_tar.extractfile("pfam_pts_names_nl_desc.txt"), header=None, index_col = 0, sep = "\t")
        #full_m.columns = pfam_pts_mapping.loc[full_m.columns, 0] - 1
        data_train_pf = data.loc[:, pfam_pts_mapping.iloc[:-93, 0]]
        data_train_pf.columns = pfam_pts_mapping.index[:-93].values - 1
        #instaniate a matrix for the individual predictions
        preds = ps.DataFrame(ps.np.zeros(shape = (len(self.sp2gene) * len(self.sample_points) * self.no_samples ,pt2 + 1 - pt1)))
        preds.columns = range(pt1, pt2 + 1)
        #print preds
        for pt_i in range(pt2 + 1 - pt1):
            #predict each sampling point independently to allow adjusting the bias term
            for i in range(len(self.sample_points)):
                data_train_pf_sp = data_train_pf.iloc[[i + j * len(self.sample_points)  for j in range(self.no_samples * len(self.sp2gene))],:]
                if adjust_bias:
                    maj_pred = predict.majority_predict(range(pt1, pt2 + 1)[pt_i], self.model_tar, data_train_pf_sp, k, bias_weight = self.sample_points[i])
                else:
                    maj_pred = predict.majority_predict(range(pt1, pt2 + 1)[pt_i], self.model_tar, data_train_pf_sp, k, bias_weight = 1)
                maj_preds_aggr = maj_pred.apply(lambda x: 1 if (x > 0).sum() > k/2  else -1, axis = 1)
                #if there is a model available fill in the phenotype prediction into the matrix
                if not len(maj_preds_aggr) == 0:
                    self.pt_has_model[range(pt1, pt2 + 1)[pt_i]] = True
                    #set the corresponding entry in the prediction matrix, make sure that this is the value not a vector, will default to na (NASTY BUG)
                    preds.iloc[[i + j * len(self.sample_points)  for j in range(self.no_samples * len(self.sp2gene))], pt_i] = maj_preds_aggr.values
        if not mask_pts is None:
            for pt in mask_pts:
                self.pt_has_model[pt] = False
        return preds.loc[:, [self.pt_has_model[i] for i in range(pt1, pt2 + 1)]]
   
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
        #get phenotypes that have less than x number of pos samples and y number of negative samples
        pt_gold_too_few_samples = pt_gold_sel_model.apply(lambda x: ps.Series(((x[~ps.isnull(x) & (x < 0)].sum() <= -self.min_samples) & (x[~ps.isnull(x) & (x > 0)].sum() >= self.min_samples))))
        #print "pt_gold_absent", pt_gold_absent
        #if not pt_gold_missing_m.all().all():
        #    print "no labels found for the input phenotypes in the gold standard"
        #    sys.exit(1)
        #unfold array and index missing phenotypes
        perf_sampling_points_pt = ps.DataFrame(ps.np.zeros(shape = (len(self.sample_points), 4 * pt_gold_sel.shape[1])))
        perf_sampling_points_pt.index = self.sample_points
        #perf_sampling_points.columns = ["recall positive class", "recall negative class", "balanced accuracy", "precision"]
        for i in range(len(self.sample_points)):
            #print "processing sampling point", sampling_points[i]
            #get confusion matrix at each sampling point i.e. 10%, 20%, ..., 100%
            #matrix to store the confusion matrices for each phenotype and sub sample
            perf_sample_pt  = ps.DataFrame(ps.np.zeros(shape = (self.no_samples, 4 * pt_gold_sel.shape[1])))
            for j in range(self.no_samples):
                #get indices  for each species at the same sampling point for the same sub sample
                loc_pred = [i + len(self.sample_points) * j + len(self.sample_points) * self.no_samples * k for k in range(len(self.sp2gene))] 
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
        conf_per_sp = ps.DataFrame(ps.np.zeros(shape = (len(self.sample_points), 4)))
        conf_per_sp.index = self.sample_points
        conf_per_sp.columns = ["TN", "FP", "FN", "TP"]
        #store performance measures per sampling point per phenotype
        perf_per_pt = ps.np.zeros(shape = (len(self.sample_points), pt_gold_sel_model.shape[1], 4))
        for sp in range(len(self.sample_points)):
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
        macro_per_sp = ps.DataFrame(nanmean(perf_per_pt[:, ((~pt_gold_absent) & pt_gold_too_few_samples).values[0], : ], axis = 1))
        macro_per_sp.index = self.sample_points
        macro_per_sp.columns = perf_colnames 
        #get micro measures 
        micro_per_sp = conf_per_sp.apply(lambda x: ps.Series([nested_cv.nested_cv.recall_pos_conf(x), nested_cv.nested_cv.recall_neg_conf(x), nested_cv.nested_cv.bacc(nested_cv.nested_cv.recall_pos_conf(x), nested_cv.nested_cv.recall_neg_conf(x)), nested_cv.nested_cv.precision_conf(x)]), axis = 1)
        micro_per_sp.columns = perf_colnames 
        #print "micro", micro_per_sp
        #write results to disk
        micro_per_sp.to_csv("%s/micro.tsv"%outdir, sep = "\t")
        macro_per_sp.to_csv("%s/macro.tsv"%outdir, sep = "\t")
        #write individual phenotype performances to disk per sampling point
        for i in range(len(self.sample_points)):
            df = ps.DataFrame(perf_per_pt[i, :, :])
            df.index = pt_gold_sel_model.columns
            df.columns = perf_colnames
            df[~pt_gold_absent].to_csv("%s/perf_sp%s.tsv"%(outdir, self.sample_points[i]), sep = "\t")


         
                


    def get_genes_targz(self, gff_targz):
        rs2gene = {}
        #read refseq id (uid) to directory mapping
        uid2fn = ps.read_csv(UID2FN, sep = "\t", index_col = 0, header = None) 
        uid2fn.index = ["RefSeqProId+%s" % i for i in uid2fn.index]
        #print self.sp2gene
        with tarfile.open(gff_targz, mode = 'r:gz') as tf:
            for rs in self.gene2hmm:
                for ce in self.rs2f[rs]:
                    print "extracted", uid2fn
                    tar = tf.getmember("./%s/%s" % (uid2fn.loc[rs, ].iloc[0], ce))
                    rs2gene[refseq] = ps.read_csv(tf.extractfile(tar), sep = "\t", index_col = 0, commentchar = '#').index
        return rs2gene

 
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser("combine the misclassified samples of different phenotypes into data matrices")
    parser.add_argument("sample_points", help='<[sampling point 1, sampling point 2, ...]> sequence of ordered relative number of genes to be retained e.g. -p 0.1,0.2,0.5,1.0')
    parser.add_argument("pfam_acc_f", help='Pfam (gene/protein) family accession files')
    parser.add_argument("no_samples", help='number of samples to be drawn for each genome', type = int)
    parser.add_argument("model_tar", help='archive with the liblinear models')
    parser.add_argument("refseq_annot_fs", help='directory with refseq chromosomal annotation files ')
    parser.add_argument("refseq_id_f", help='file with refseq ids')
    parser.add_argument("pt_range", help='range of phenotypes e.g. 8476-8568')
    parser.add_argument("pt_gold_f", help='gold standard matrix with the phenotypes')
    parser.add_argument("pt2id_f", help='out_dir for macro, micro accuracy etc.')
    parser.add_argument("outdir", help='out_dir for macro, micro accuracy etc.')
    parser.add_argument("sp2str_f", help='mapping file that maps RefSeqProIds to species')
    parser.add_argument("min_samples", help='number of samples required for the phenotype + and - class to be considered in the macro accuracy computation', type = int)
    parser.add_argument("mask_pts", help='mask those phenotypes e.g. 8476,8521,8533')
    parser.add_argument("--adjust_bias", help='if set adjust the bias term according to the percentage of genes retained', action = 'store_true')
    args = parser.parse_args()
    pt1, pt2 = [int(i) for i in args.pt_range.split("-")]
    args.sample_points = [float(i) for i in args.sample_points.split(",")] 
    import os
    if os.path.exists(args.outdir):
        sys.stderr.write("outdir %s already exists"%args.outdir)
        sys.exit(1)
    else:
        os.mkdir(args.outdir)
    
    simulator = SimulateDraftGenome(args.sample_points, args.no_samples, args.refseq_annot_fs, args.refseq_id_f, args.pfam_acc_f,  tarfile.open(args.model_tar, mode = "r:gz"),  args.sp2str_f, args.min_samples)
    data = simulator.simulate_draft_genomes(args.pfam_acc_f)
    preds = simulator.predict_draft_genomes(data, pt1, pt2, 5, args.mask_pts, args.adjust_bias)
    simulator.evaluate_draft_genome_preds(preds, args.pt_gold_f, args.pt2id_f, args.outdir)
