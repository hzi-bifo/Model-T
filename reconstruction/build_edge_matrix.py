import sys
import os.path
import numpy as np
import pandas as pd
import dendropy as dp

class build_edge_matrix:

    def __init__(self,tree_f, format, phyn_f, phym_f, event_f):
        t = dp.Tree()
        self.tree = t
        t.read_from_path(tree_f, format, suppress_internal_node_taxa=False)
        self.phy = pd.read_csv(pf_pt_f, index_col = 0, sep = "\t")
        self.char2ev = self.read_events(event_f)
        self.missing_characters = ("?", "-", "-1.0")
        self.event_f=event_f

    def get_edges(self, pt):
        edges =[]
        #recurse through tree to get edges
        tl, tr = self.tree.seed_node.child_nodes()
        root_update = self.get_edges_h( tl, tr,self.tree.seed_node, edges, pt)
        if type(root_update) == type([]):
            #the root is obsolete but needs to be kept
            for e in edges:
                if root_update[0] in e:
                    e += root_update[1:(len(root_update))]
        return edges

    def get_edges_h(self, l, r,node, edges, pt):
        if l.is_leaf() and r.is_leaf():
            return self.update(l, r, node, edges, pt)
        elif l.is_leaf():
            rl, rr = r.child_nodes()
            return self.update(l, self.get_edges_h(rl, rr, r, edges, pt), node, edges, pt)
        elif r.is_leaf():
            ll, lr = l.child_nodes()
            return self.update(self.get_edges_h(ll, lr, l, edges, pt), r, node, edges, pt)
        else:
            ll, lr = l.child_nodes()
            rl, rr = r.child_nodes()
            return self.update(self.get_edges_h(ll, lr, l, edges, pt), self.get_edges_h(rl, rr,r , edges, pt), node, edges, pt)


    def update(self, l, r, node, edges, pt):
        #both children stem from missing values (NAs)
        l1, r1 = node.child_nodes()
        if l is None and r is None:
            return None
        #node l is ancestor of NA node return the right node
        elif l is None:
            if type(r) == type([]):
                r.append(node.taxon.label)
                return r
            elif type(r) == type(True):
                return [r1.taxon.label, node.taxon.label]
            elif type(r) == type(node):
                if str(self.phy.loc[r.taxon.label, pt]) not in self.missing_characters:
                    return [r.taxon.label, node.taxon.label]
                else:
                    return None
        #node r is ancestor of NA node return the left node
        elif r is None:
            if type(l) == type([]):
                l.append(node.taxon.label)
                return l
            elif type(l) == type(True):
                return [l1.taxon.label, node.taxon.label]
            elif type(l) == type(node):
                if str(self.phy.loc[l.taxon.label, pt]) not in self.missing_characters:
                    return [l.taxon.label, node.taxon.label]
                else:
                    return None
        #one node is list the other leave
        elif type(l) == type(node) and type(r) == type([]):
            if str(self.phy.loc[l.taxon.label, pt]) not in self.missing_characters:
                r.append(node.taxon.label)
                edges.append(r)
                edges.append([l.taxon.label, node.taxon.label])
                return True
            else:
                r.append(node.taxon.label)
                return r
        elif type(r) == type(node) and type(l) == type([]):
            if str(self.phy.loc[r.taxon.label, pt]) not in self.missing_characters:
                l.append(node.taxon.label)
                edges.append(l)
                edges.append([r.taxon.label, node.taxon.label])
                return True
            else:
                l.append(node.taxon.label)
                return l
        #both nodes derive from non NA nodes
        elif type(l) == type([]) and type(r) == type([]):
            l.append(node.taxon.label)
            r.append(node.taxon.label)
            edges.append(l)
            edges.append(r)
            return True
        #one of the nodes is sane and isn't an ancestor of a NA node
        elif type(l) == type(True) or type(r) == type(True):
            #check if left node is either type node or list of nodes
            if type(l) != type(True):
                #check if left node is leaf node
                if type(l) == type(l1) and l.is_leaf():
                    if str(self.phy.loc[l.taxon.label, pt]) not in self.missing_characters:
                        edges.append([l.taxon.label, node.taxon.label])
                        edges.append([r1.taxon.label, node.taxon.label])
                        return True
                    else: return [ r1.taxon.label, node.taxon.label]
                #check if left node is type list
                elif type(l) == type([]):
                    l.append(node.taxon.label)
                    #l.append(node.taxon.label)
                    edges.append(l)
                    edges.append([r1.taxon.label, node.taxon.label])
                    return True
            #check if right node is either type node or list of nodes
            elif type(r) != type(True):
                #check if right node is leave node
                if type(r) == type(r1) and r.is_leaf():
                    if str(self.phy.loc[r.taxon.label, pt]) not in self.missing_characters:
                        edges.append([r.taxon.label, node.taxon.label])
                        edges.append([l1.taxon.label, node.taxon.label])
                        return True
                    else: return [l1.taxon.label, node.taxon.label]
                #check if right node is type list
                elif type(r) == type([]):
                    #r.append(r1.taxon.label)
                    r.append(node.taxon.label)
                    edges.append(r)
                    edges.append([l1.taxon.label, node.taxon.label])
                    return True
            else:
                edges.append([l1.taxon.label,node.taxon.label])
                edges.append([r1.taxon.label, node.taxon.label])
                return True
        #both nodes are leave nodes
        elif l.is_leaf() and r.is_leaf():
            if str(self.phy.loc[l.taxon.label, pt]) not in self.missing_characters and str(self.phy.loc[r.taxon.label, pt]) not in self.missing_characters:
                edges.append([l.taxon.label,node.taxon.label])
                edges.append([r.taxon.label,node.taxon.label])
                return True
            elif str(self.phy.loc[l.taxon.label, pt]) not in self.missing_characters:
                return [l.taxon.label, node.taxon.label]
            elif str(self.phy.loc[r.taxon.label, pt]) not in self.missing_characters:
                return [r.taxon.label, node.taxon.label]
            else:
                return None


    def read_events(self,event_f):
        """read parsimony events for all characters"""
        char2ev = {}
        f = open(event_f, 'r')
        #skip header
        f.readline()
        for l in f:
            char, n1,n2,etype, ecount = l.strip().split("\t")
            if not int(char) in char2ev:
                char2ev[int(char)] = [(n1,n2,etype,ecount)]
            else: char2ev[int(char)].append((n1,n2,etype,ecount))
        f.close()
        return char2ev

    def get_pt_dict(self,edges):
        """for a given pt character get a mapping from nodes to edges """
        node2edge = {}
        for e in edges:
            for node in e:
                if not node in node2edge:
                    node2edge[node] = set()
                    node2edge[node].add(tuple(e))
                else:
                    node2edge[node].add(tuple(e))

        return node2edge

    def map_events(self,node2edge, gt_start, gt_end, char2ev, pt, edges):
        """map events from all  genes in gt_range to pt edges and add pt states too"""
        edge2char2val = {}
        for e in edges:
            edge2char2val[tuple(e)] = {}
        #map/match genotype to phenotype edges
        for gt in range(gt_start, gt_end):
            #TODO account for what happen if there are no events
            if not gt in char2ev:
                continue
            for ev in char2ev[gt]:
                #map gt event to pt event and discard if not possible
                if ev[0] in node2edge:
                    edge1 = node2edge[ev[0]]
                else: continue
                if ev[1] in node2edge:
                    edge2 = node2edge[ev[1]]
                else: continue
                #print edge1, edge2
                isn = edge1.intersection(edge2)
                if len(isn) == 0:
                    print "event",ev
                    print "edge1",edge1
                    print "edge2",edge2
                isn = isn.pop()
                if isn in edge2char2val and gt in edge2char2val[isn]:
                    if ev[2] == 'loss':
                        edge2char2val[isn][gt] -= int(ev[3])
                    else: edge2char2val[isn][gt] += int(ev[3])
                elif isn in edge2char2val and gt not in edge2char2val[isn]:
                    if ev[2] == 'loss':
                        edge2char2val[isn][gt] = -int(ev[3])
                    else: edge2char2val[isn][gt] = int(ev[3])
        #account for phenotype edges
        #print "edge2char2val", edge2char2val
        for ev in char2ev[pt]:
            #print ev
            edge1 = node2edge[ev[0]]
            edge2 = node2edge[ev[1]]
            isn = edge1.intersection(edge2).pop()
            if ev[2] == 'loss':
                edge2char2val[isn][pt]=-int(ev[3])
            else: edge2char2val[isn][pt]=int(ev[3])
        return edge2char2val

    def get_edge_m(self,edge2char2val, edges, feats, pt, out_f, is_internal = False):
        """generate for each edge a vector  of all characters"""
        #out_fo = open(out_f, 'w')
        if not pt in feats:
            out_m = pd.DataFrame(np.zeros(shape = (len(edges), len(feats) + 1)))
            out_m.columns = feats + [pt]
        else:
            out_m = pd.DataFrame(np.zeros(shape = (len(edges), len(feats))))
            out_m.columns = feats 
        out_m.index = ["_".join(e) for e in edges]
        #only consider the phenotype if it's not part of the genotypes ergo no missing values are involved
        for e in edges:
            #s="%s\t"%str("_".join(e))
            for gt in feats:
                if gt in edge2char2val[tuple(e)]:
                    out_m.loc["_".join(e), gt] = edge2char2val[tuple(e)][gt]
                    #s+="%s\t"%edge2char2val[tuple(e)][gt]
                #else:
                #    s+="0\t"
            #only consider the phenotype if it's not part of the genotypes ergo no missing values are involved
            if pt not in feats and pt in edge2char2val[tuple(e)]:
                #print "this shouldn't happen"
                #sys.exit(0)
                out_m.loc["_".join(e), pt] = edge2char2val[tuple(e)][pt]
                #s+="%s\n"%edge2char2val[tuple(e)][pt]
                #else: s+="0\n"
                #out_fo.write(s)
        out_m.to_csv(out_f, sep = "\t")
        if is_internal: 
            return out_m 


    def get_all_edge_m(self, feats, pts, out_dir, is_internal = False):
        """for all phenotypes generate a edge based matrix"""
        for pt in pts:
            #check if phenotype has any events associated with it
            if pt not in self.char2ev:
                print "phenotype", pt, "has no events associated with it. Skipping.."
                continue
            edges = self.get_edges(pt)
            pt_dict = self.get_pt_dict(edges)
            #plus one because we want to include column gt_end
            edge2char2val = self.map_events(pt_dict, feats, self.char2ev, pt, edges)
            m = self.get_edge_m(edge2char2val, edges, feats, pt, "%s/pt%s.dat"%(out_dir,pt), is_internal = is_internal)
            if is_internal:
                return m

if __name__ == '__main__':
    parser = argparse.ArgumentParser("reconstruct likelihood matrix from gainLoss output")
    parser.add_argument("tree", help='phylogenetic tree')
    parser.add_argument("phypat_pt_f", help='phylogenetic tree')
    parser.add_argument("--format", default = "newick", help='phylogenetic tree format')
    parser.add_argument("event_f", help='gainLoss tabular output file')
    parse.add_argument("feature_f", help = "list of features used")
    parse.add_argument("phenotype_f", help = "list of phenotypes")
    parse.add_argument("outdir", help = "output directory")
    a = parser.parse_args()
    #check if the directory already exists
    if os.path.exists(out):
        sys.stderr.write("output directory %s already exists; delete and rerun\n"%a)
        sys.exit(1)
    else:
        os.mkdir(a.outdir)
    bem = build_edge_matrix(a.tree, a.format, a.phypat_pt_f, a.event_f)
    #read in features 
    feats = pd.read_csv(feat_f, sep = "\t", index_col = 0).index.tolist()
    #read in phenotypes 
    pts = pd.read_csv(pt_f, sep = "\t", index_col = 0).index.astype('string').tolist()
    bem.get_all_edge_m(feats, pts, args.outdir)
