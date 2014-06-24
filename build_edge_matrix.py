import numpy
import sys
import dendropy as dp
class build_edge_matrix:

    def __init__(self,tree_f, format, phyletic_f, event_f):
        t = dp.Tree()
        self.tree = t
        t.read_from_path(tree_f, format, suppress_internal_node_taxa=False)
        self.phy = dp.StandardCharacterMatrix.get_from_path(phyletic_f, "fasta")
        self.char2ev = self.read_events(event_f)
        self.missing_characters = ("?", "-")
        self.event_f=event_f

    def get_edges(self, pt):
        edges =[]
        #recurse through tree to get edges
        tl, tr = self.tree.seed_node.child_nodes()
        self.get_edges_h( tl, tr,self.tree.seed_node, edges, pt)
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


    def update(self, l, r,node, edges, pt):
        #both l and r are leaves
        print node, l, r
        if l is None and r is None:
            return None
        elif l is None:
            return r.append(node.taxon.label)
        elif r is None:
            return l.append(node.taxon.label)
        elif l == True or r == True and not l.is_leaf():
            l1, r1 = node.child_nodes()
            if type(l) != type(True):
                l.append(node.taxon.label)
                edges.append(l)
                edges.append((r1.taxon.label,node.taxon.label))
            if type(r) != type(True):
                r.append(node.taxon.label)
                edges.append((l1.taxon.label,node.taxon.label))
                edges.append(r)
            else:
                edges.append((l1.taxon.label,node.taxon.label))
                edges.append((r1.taxon.label, node.taxon.label))
            return True

        if l.is_leaf() and r.is_leaf():
            if str(self.phy[l.taxon.label][pt]) not in self.missing_characters and str(self.phy[r.taxon.label][pt]) not in self.missing_characters:
                edges.append((l.taxon.label,node.taxon.label))
                edges.append((r.taxon.label,node.taxon.label))
                return True
            elif str(self.phy[l.taxon.label][pt]) not in self.missing_characters:
                return [l.taxon.label, node.taxon.label]
            elif str(self.phy[r.taxon.label][pt]) not in self.missing_characters:
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
        return char2ev

    def get_pt_dict(self,edges, char2ev):
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
        """map events from all  genes in gt_range to pt edges and add pt states to"""
        edge2char2val = {}
        for e in edges:
            edge2char2val[tuple(e)] = {}
        #map/match genotype to phenotype edges
        for gt in range(gt_start, gt_end):
            #TODO account for what happen if there are no events
            for ev in char2ev[gt]:
                #map gt event to pt event and discard if not possible
                print ev[0],ev[1]
                if ev[0] in node2edge:
                    edge1 = node2edge[ev[0]]
                else: continue
                if ev[1] in node2edge:
                    edge2 = node2edge[ev[1]]
                else: continue
                print edge1, edge2
                print edge1.intersection(edge2)
                edge2char2val[edge1.intersection(edge2).pop()][gt]=ev[2]
        #account for phenotype edges
        print "edge2char2val", edge2char2val
        for ev in char2ev[pt]:
            edge1 = node2edge[ev[0]]
            edge2 = node2edge[ev[1]]
            edge2char2val[edge1.intersection(edge2).pop()][gt]=ev[2]
        return edge2char2val

    def get_edge_m(self,edge2char2val, edges, gt_start, gt_end,pt, out_f):
        """generate for each edge a vector  of all characters"""
        out_fo = open(out_f, 'w')
        for e in edges:
            s="%s\t"%str(e)
            for gt in range(gt_start, gt_end):
                if gt in edge2char2val[tuple(e)]:
                    s+="%s\t"%edge2char2val[tuple(e)][gt]
                else:
                    s+="0\t"
            if pt in edge2char2val[tuple(e)]:
                s+="%s\n"%edge2char2val[tuple(e)][pt]
            else: s+="0\n"
            out_fo.write(s)


if __name__ == '__main__':
    import getopt
    if len(sys.argv) == 1:
        print """USAGE: python %s
-t <tree> in a dendropy compliant format
-f <format> above format either one of nexus, newick, nexml, phylip
-p <phyletic pattern> phyletic_patterns in fasta format
-e <events> file with parsimony events
-g <range of genotypes> e.g. 1-8400
-pt <range of phenotypes> to consider e.g 8550-8560
> <out file> with gain and loss events
        """ % (sys.argv[0])
    t = "/net/metagenomics/projects/phenotypes_20130523/code/build_edge_m/sampletree.newick"
    p = "/net/metagenomics/projects/phenotypes_20130523/code/build_edge_m/charm.fasta"
    e = "/net/metagenomics/projects/phenotypes_20130523/code/build_edge_m/events.txt"
    f = "newick"
    bem = build_edge_matrix(t,f,p,e)
    #test edge generation
    e = bem.get_edges(0)
    print e
    for i in e:
        s=""
        for j in i:
            s+="\t%s"%j
        print s

    #test edge matrix generation
    char2ev =  bem.char2ev
    pt_dict = bem.get_pt_dict(e, char2ev)
    print pt_dict
    edge2char2val = bem.map_events(pt_dict, 1,2, char2ev,0, e)
    #print edge2char2val
    bem.get_edge_m(edge2char2val, e, 0,1,0, "out_matrix.test.txt")
    sys.exit(1)

    try:
        optlist, args = getopt.getopt(sys.argv[1:], "t:f:p:e:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    t=None
    f=None
    p=None
    e = None
    for o, a in optlist:
        if o == "-t":
            t = a
        if o == "-f":
            f = a
        if o == "-p":
            p = a
        if o == "-e":
            e = a
    bem = build_edge_matrix(t,f,p,e)

