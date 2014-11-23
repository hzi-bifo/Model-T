import pandas as ps
import numpy as np
from ete2 import TreeNode, Tree
#import dendropy as dp
import sys

class node_anc_rec:

    def __init__(self, tree_f, majority_feat_f, gain_m_f, loss_m_f, pfam_pts_f, pfam_ids_f, pt_id, feat_list):
        #read in gain and loss event matrices 
        self.gain_m = ps.read_csv(gain_m_f, sep = "\t", index_col = 0, header = None)
        self.loss_m = ps.read_csv(loss_m_f, sep = "\t", index_col = 0, header = None)
        #reindex gain_m and loss_m by getting rid of the internal nodes
        self.gain_m.index = ["_".join((i.split("_")[0], i.split("_")[len(i.split("_")) - 1])) for i in self.gain_m.index.values]
        #########print ["_".join((i.split("_")[0], i.split("_")[len(i.split("_")) - 1])) for i in self.gain_m.index.values]
        #sys.exit(1)
        self.loss_m.index = ["_".join((i.split("_")[0], i.split("_")[len(i.split("_")) - 1])) for i in self.loss_m.index.values]
        #self.loss_m.drop(8477, inplace = True, axis = 1)
        #self.gain_m.drop(8477, inplace = True, axis = 1)
        #self.loss_m.rename(columns = {8478:8477}, inplace = True)
        #self.gain_m.rename(columns = {8478:8477}, inplace = True)
        #read in majority features
        if not majority_feat_f is None:
            self.feats = ps.read_csv(majority_feat_f, index_col = 1, sep = "\t", header = True)
        else:
            self.feats = ps.DataFrame(feat_list)
            self.feats.index = self.feats.iloc[:, 0]
        pf2id = ps.read_csv(pfam_ids_f, sep = "\t", header = None, index_col = 1)
        self.feat_ids = pf2id.loc[self.feats.index.values, 0] 
        self.feat_ids = self.feat_ids
        #prune full tree to pt tree
        self.pt_tree = self.prune(tree_f, self.gain_m.index.values)
        #print self.pt_tree.write(format = 8)
        #print self.pt_tree
        #node state matrix
        node_names = [i.name for i in self.pt_tree.traverse(strategy = 'preorder')]
        self.node_state_m = ps.DataFrame(np.zeros(shape = (len(node_names), self.feats.shape[0] + 1)))
        #node edge state matrix
        self.node_edge_state_m = ps.DataFrame(np.zeros(shape = (len(node_names), self.feats.shape[0] + 1)))
        self.node_state_m.index = node_names
        self.node_edge_state_m.index = node_names
        self.node_state_m.columns = self.feat_ids.append(ps.Series([8477]))
        self.node_edge_state_m.columns = self.feat_ids.append(ps.Series([8477]))
        #set leave node states
        pfam_pts_m = ps.read_csv(pfam_pts_f, sep = "\t",  index_col = 0, header = None)
        #convert int index into string index
        pfam_pts_m.index = [str(i) for i in pfam_pts_m.index.values]
        #set node states of the leave nodes
        self.node_state_m.loc[self.pt_tree.get_leaf_names(), self.feat_ids] = pfam_pts_m.loc[self.pt_tree.get_leaf_names(), self.feat_ids]
        self.node_state_m.loc[self.pt_tree.get_leaf_names(), 8477] = ps.Series(pfam_pts_m.loc[self.pt_tree.get_leaf_names(), pt_id + 1], dtype = 'int')
        #print self.node_state_m.loc[:, 8477] 


    def reconstruct(self):
        #traverse tree to determine the internal node states based on the edge events
        self.get_internal_states(self.pt_tree, 8477)
        for i in self.feat_ids:
            self.get_internal_states(self.pt_tree, i)
        return self.pt_tree, self.node_edge_state_m, self.node_state_m
    
    def get_internal_states(self, node, feat_id):
        if node.is_leaf():
            #get edge events
            n_g_ev = self.gain_m.loc[node.name + "_" + node.up.name, feat_id]
            n_l_ev = self.loss_m.loc[node.name + "_" + node.up.name, feat_id]
            #loss event
            if n_l_ev != 0 and n_g_ev == 0: 
                self.node_edge_state_m.loc[node.name, feat_id] = -n_l_ev
            #gain event
            if n_g_ev != 0 and n_l_ev == 0: 
                self.node_edge_state_m.loc[node.name, feat_id] = n_g_ev
            #conflicting event
            if n_g_ev != 0 and n_l_ev != 0:
                self.node_edge_state_m.loc[node.name, feat_id] = np.NaN 
            #no event
            if n_g_ev == 0 and n_l_ev == 0:
                self.node_edge_state_m.loc[node.name, feat_id] = 0 
            #print self.node_edge_state_m.loc[node.name, feat_id] 
            #print node.name
            #print n_g_ev
            #print n_l_ev
            return 
        l, r  = node.get_children()
        self.get_internal_states(l, feat_id)
        self.get_internal_states(r, feat_id)
        if node.is_root():
            return 
        #get states of left and right nodes
        l_state = self.node_state_m.loc[l.name, feat_id]
        r_state = self.node_state_m.loc[r.name, feat_id]
        #check events on branches from l to current node and r to current node and assign state of current node accordingly
        l_ev = self.node_edge_state_m.loc[l.name, feat_id]
        r_ev = self.node_edge_state_m.loc[r.name, feat_id]
        n_g_ev = self.gain_m.loc[node.name + "_" + node.up.name, feat_id]
        n_l_ev = self.loss_m.loc[node.name + "_" + node.up.name, feat_id]
        #set edge event of node 
        #loss event
        #print l_state, r_state
        #print l_ev, r_ev
        #print l.name, r.name
        if n_l_ev != 0 and n_g_ev == 0: 
            self.node_edge_state_m.loc[node.name, feat_id] = -n_l_ev
        #gain event
        if n_g_ev != 0 and n_l_ev == 0: 
            self.node_edge_state_m.loc[node.name, feat_id] = n_g_ev
        #conflicting event
        if n_g_ev != 0 and n_l_ev != 0:
            self.node_edge_state_m.loc[node.name, feat_id] = np.NaN 
        #no event
        if n_g_ev == 0 and n_l_ev == 0:
            self.node_edge_state_m.loc[node.name, feat_id] = 0 
        #set node event
        #node states of children is ambiguous
        if l_state == -1 or r_state == -1:
            self.node_state_m.loc[node.name, feat_id] = -1
        #check if states of children correspond
        if l_state - l_ev == r_state - r_ev:
            self.node_state_m.loc[node.name, feat_id] = r_state - r_ev
        #node states of children do not correspond
        else:
            self.node_state_m.loc[node.name, feat_id] = np.NaN

    
    
    def prune(self, tree_f, edges): 
        pt_leaves = []
        #read in tree
        ncbi_ids = tree_f
        t = TreeNode(tree_f,1)
        leaves_full_tree = set(t.get_leaf_names())
        #print t.write(format = 1)
        #print leaves_full_tree
        #get leaves nodes in full tree that are also part of the phenotype tree
        for e in edges:
            if e.split("_")[0] in leaves_full_tree:
                pt_leaves.append(e.split("_")[0])
        t.prune(pt_leaves)
        return t
    

if __name__ == '__main__':
    import getopt
    if len(sys.argv) == 1:
        print """USAGE: python %s
-g matrix with gain events
-l matrix with loss events
-t tree with all taxa
-f features to be mapped on the tree
-a feature ids 
-p phyletic pattern matrix
-t query phenotype id 
> <out file> with gain and loss events
        """ % (sys.argv[0])
        sys.exit(1)
    try:
        optlist, args = getopt.getopt(
                sys.argv[1:], "g:l:t:f:p:a:q:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    g = l = f = t = p = pfam_ids_f  = q = None

    for o, a in optlist:
        if o == "-g":
            g = a
        if o == "-l":
            l = a
        if o == "-t":
            t = a
        if o == "-f":
            f = a
        if o == "-p":
            p = a
        if o == "-a":
            pfam_ids_f = a
        if o == "-q":
            q = int(a)
    node_recon = node_anc_rec(t, f, g, l, p, pfam_ids_f, q)
    node_recon.reconstruct()
    
    
