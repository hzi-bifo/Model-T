"""script to simulate genomes with incomplete gene sets from fully sequenced genomes"""
import predict
import domtblout2gene
import random 
import math
import pandas as ps
import sys
import tarfile
import nested_cv


UID2FN = "/net/metagenomics/data/genomes/public/NCBI20140115_FAAs_uid2fn.txt"
class SimulateDraftGenome:

    def __init__(self, sample_points, no_samples, refseq_annot_fs, refseq_f, pfam_acc_f, model_dir, train_pfam_f):
        self.sample_points = sample_points
        self.no_samples = no_samples
        self.model_dir = model_dir
        self.train_pfam_f = train_pfam_f
        fs = [i.strip() for i in open(refseq_annot_fs, 'r').readlines()]
        annot, self.gene2hmm, self.rs2f = domtblout2gene.gene2hmm(fs, pfam_acc_f)
        print annot.sum()
        print self.rs2f
        #self.rs2gene = ps.read_csv(gene_length_f, sep = "\t") 
        #self.rs2gene.index = ["RefSeqProId+%s" % i for i in rs2gene.index]

    def simulate_draft_genomes(self):
        """create for each fully sequenced genome, x number of samples for y number of sample points"""
        #read in pfam accs
        pfams = ps.read_csv(pfam_acc_f, index_col = 0, header = None).index
        #sample Pfam matrices at different samples points
        full_m = ps.DataFrame(ps.np.zeros(shape = (len(self.gene2hmm) * len(self.sample_points) * self.no_samples, len(pfams))))
        full_m.columns = pfams
        #pf_sampled = {}
        for j in range(len(self.gene2hmm)):
            for i in range(self.no_samples):
                sample_points_t = self.sample_points[::-1]
                m = ps.DataFrame(ps.np.zeros(shape = (len(self.sample_points), len(pfams) )))
                m.index = self.sample_points
                m.columns = pfams
                cur_sample_point = sample_points_t.pop()
                cur_sample_series = m.iloc[0,].copy()
                genes = random.sample(self.gene2hmm[self.gene2hmm.keys()[j]], len(self.gene2hmm[self.gene2hmm.keys()[j]])) 
                #draw random sample
                print "draw sample", self.gene2hmm.keys()[j]
                #retrieve vector of pfams 
                for k in range(len(genes)):
                    if k > cur_sample_point * len(genes):
                        print i, sample_points_t, cur_sample_point
                        m.loc[cur_sample_point, ] = cur_sample_series
                        if sample_points_t == []:
                            break
                        cur_sample_point = sample_points_t.pop()
                    for pf in self.gene2hmm[self.gene2hmm.keys()[j]][genes[k]]:
                        cur_sample_series.loc[pf] += 1
                print j * len(self.sample_points) * self.no_samples + len(self.sample_points) * i
                print (j * len(self.sample_points) * self.no_samples + len(self.sample_points) * i + len(self.sample_points))
                print m.shape
                print full_m.shape
                full_m.iloc[(j * len(self.sample_points) * self.no_samples + len(self.sample_points) * i) : (j * len(self.sample_points) * self.no_samples + len(self.sample_points) * i + len(self.sample_points)), ] = m.values
                #if self.gene2hmm.keys()[j] not in pf_sampled:
                #    pf_sampled[self.gene2hmm.keys()[j]] = [m]
                #else: 
                #    pf_sampled[self.gene2hmm.keys()[j]].append(m)
        return full_m

   
    def evaluate_draft_genome_preds(self, preds, pt_gold_f, pt2id_f):
        """read in gold standard pt predictions"""
        print preds
        #evaluation per sampling point averaged over all phenotypes; this might be skewed due to very few test samples for some phenotypes
        #evalution over all predictions per phenotype
        #evaluation per sapling point per phenotype -> influence on specific phenotypes
        #matrix: 
        #Rf1s1 10%
        #Rf1s1 20%
        #Rf1s2 10%
        #Rf1s2 20%
        #Rf2s1 10%
        pt_gold_df = ps.read_csv(pt_gold_f, sep = "\t", index_col = 0)
        print pt_gold_df.columns
        #read in pt2id file
        pt2id = ps.read_csv(pt2id_f, index_col = 0, sep = "\t", header = None) 
        pt_gold_df.columns = pt2id.index
        #restrict to pt that were selected 
        print pt_gold_df
        pt_gold_sel = pt_gold_df.loc[:, preds.columns]
        #restrict to pt that have a model avialble
        pt_gold_sel_model = pt_gold_sel[[self.pt_has_model[i] for i in range(pt2, pt1)]]
        pt_gold_array = pt_gold_sel_model.values.flatten(order = 'C')
        pt_gold_missing = ps.notnull(pt_gold_array)
        if pt_gold_missing.all():
            print "no labels found for the input phenotypes in the gold standard"
            sys.exit(1)
        #unfold array and index missing phenotypes
        for i in range(len(sampling_points)):
            #get accuracy at each sampling point
            accs = ps.np.zeros(shape = no_samples)
            for j in range(no_samples):
                loc_pred = [i + i * j + i * j * k for k in range(len(self.gene2hmm))] 
                labels = pt_gold_df 
                print preds.iloc[loc_pred, ]
                print loc_pred
                rec_pos = nested_cv.nested_cv.recall_pos(pt_gold_array[pt_gold_missing], preds.iloc[loc_pred, ].values.flatten()[pt_gold_missing])
                rec_neg = nested_cv.nested_cv.recall_neg(pt_gold_array[pt_gold_missing], preds.iloc[loc_pred, ].values.flatten()[pt_gold_missing])
                print rec_pos
                print rec_neg

                

    def predict_draft_genomes(self,data, pt1, pt2, k):
        """"""
        self.pt_has_model = dict((i, False) for i in range(pt2, pt1)) 
        #restrict the data matrix to the pfams from the training data set
        data_train_pf = data.loc[:, ps.read_csv(self.train_pfam_f, index_col = 0, header = None).index]
        #instaniate a matrix for the individual predictions
        preds = ps.DataFrame(ps.np.zeros(shape = (len(self.gene2hmm) * len(self.sample_points) * self.no_samples ,pt2 + 1 - pt1)))
        preds.columns = range(pt1, pt2 + 1)
        for pt in range(pt1, pt2 + 1):
            maj_pred = predict.majority_predict(pt, self.model_dir, data_train_pf, k)
            maj_preds_aggr = maj_pred.apply(lambda x: (x > 0).sum() > k/2, axis = 1)
            #print maj_preds_aggr
            #if there is a model available fill in the phenotype prediction into the matrix
            if not maj_preds_aggr is None:
                self.pt_has_model[pt] = True
                preds.loc[:, pt] = maj_preds_aggr 
            #preds.loc[:, pt] = 
        return preds.loc[:, [self.pt_has_model[i] for i in range(pt1, pt2 + 1)]]

    def get_genes_targz(self, gff_targz):
        rs2gene = {}
        #read refseq id (uid) to directory mapping
        uid2fn = ps.read_csv(UID2FN, sep = "\t", index_col = 0, header = None) 
        uid2fn.index = ["RefSeqProId+%s" % i for i in uid2fn.index]
        print self.gene2hmm
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
""" % (sys.argv[0])
        sys.exit(1)
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "s:p:g:m:r:a:b:c:d:i:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    no_samples = sampling_points = gff_targz = pfam_acc_f = model_dir = refseq_annot_fs = refseq_f = tr_pfams_f = pt_gold_f = pt2id_f = None
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

    simulator = SimulateDraftGenome(sampling_points, no_samples, refseq_annot_fs, refseq_f, pfam_acc_f,  model_dir, tr_pfams_f)
    data = simulator.simulate_draft_genomes()
    preds = simulator.predict_draft_genomes(data, pt1, pt2, 5)
    simulator.evaluate_draft_genome_preds(preds, pt_gold_f, pt2id_f)