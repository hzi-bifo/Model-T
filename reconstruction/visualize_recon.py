import pandas as pd
import ete2
from ete2 import Tree, faces, AttrFace, TreeStyle
import pylab

kelly_colors_hex = [
    0xFFB300, # Vivid Yellow
    0x803E75, # Strong Purple
    0xFF6800, # Vivid Orange
    0xA6BDD7, # Very Light Blue
    0xC10020, # Vivid Red
    0xCEA262, # Grayish Yellow
    0x817066, # Medium Gray

    # The following don't work well for people with defective color vision
    0x007D34, # Vivid Green
    0xF6768E, # Strong Purplish Pink
    0x00538A, # Strong Blue
    0xFF7A5C, # Strong Yellowish Pink
    0x53377A, # Strong Violet
    0xFF8E00, # Vivid Orange Yellow
    0xB32851, # Strong Purplish Red
    0xF4C800, # Vivid Greenish Yellow
    0x7F180D, # Strong Reddish Brown
    0x93AA00, # Vivid Yellowish Green
    0x593315, # Deep Yellowish Brown
    0xF13A13, # Vivid Reddish Orange
    0x232C16, # Dark Olive Green
    ]

def my_layout(node):
    if node.is_leaf():
         # If terminal node, draws its name
         name_face = AttrFace("name")
    else:
         # If internal node, draws label with smaller font size
         name_face = AttrFace("name", fsize=10)
    # Adds the name face to the image at the preferred position
    faces.add_face_to_node(name_face, node, column=0, position="branch-right")
    

def get_style():
    ts = TreeStyle()
    # Do not add leaf names automatically
    ts.show_leaf_name = False
    ts.show_scale = True 
    ts.force_topology = True
    # Use my custom layout
    ts.layout_fn = my_layout
    return ts

def plot_tree(pt_tree, target_node, out):
    #pt_tree, feats, pf2color = get_tree(phenotype = phenotype, feat_list = "top_cor", is_ml_plus_phypat = True, target_node = target_node)
    pt_tree.dist = 0
    target = pt_tree.search_nodes(name = target_node)[0]
    target.render("%%inline",  tree_style = get_style())
    target.render(out + '_tree.svg',  tree_style = get_style())
    target.render(out + '_tree.png', tree_style = get_style())
    return target, feats, pf2color

def plot_legend(feats, out, pf2color,  pf_desc = False, pf_acc = True, include_class = False):
    fig = pylab.figure()
    figlegend = pylab.figure(figsize = (9, 6))
    ax = fig.add_subplot(111)
    x = [0,1]
    lines = [ax.plot(x, pd.np.ones(len(x)), 'o', color = "#%06x" % (pf2color[feats.index[i]]))[0] for i in range(feats.shape[0])]
    #labels= ["%s" %(feats.loc[:,"Pfam_acc"].iloc[i]) for i in range(feats.shape[0])]
    labels= [i for i in feats.index]
    #if include_class:
    #    labels= ["%s %s" %(labels[i], feats.loc[:, "class"].iloc[i]) for i in range(len(labels))]
    #if pf_desc:
    #    labels = ["%s %s" % (labels[i], pf2short_desc.loc[feats.loc[:,"Pfam_acc"].iloc[i], 1]) for i in range(len(labels))]
    #if pf_acc:
    #    labels = ["%s %s" % (labels[i], pf2acc.loc[feats.loc[:,"Pfam_acc"].iloc[i], 1]) for i in range(len(labels))]

    figlegend.legend(lines, labels, markerscale = 2.5, numpoints = 1, frameon = False)
    #fig.show()
    fig.tight_layout()
    figlegend.savefig(out + "_legend.svg")
    figlegend.savefig(out + "_legend.png")
    return figlegend

def get_tree(phenotype, tree, gain_recon, loss_recon, node_recon, pfam_mapping, feat_list, sample_mapping, threshold = 0.5, target_node = None):
    #read target feats
    feats = pd.read_csv(feat_list, index_col = 0, sep = "\t")
    pt_tree = ete2.Tree(tree, format = 1)
    sample_mapping = pd.read_csv(sample_mapping, index_col = 0, sep = "\t")
    #read node and edge reconstruction matrices
    node_recon = pd.read_csv(node_recon, sep = "\t", index_col = 0)
    gain_recon = pd.read_csv(gain_recon, sep = "\t", index_col = 0)
    gain_recon.index = ["_".join((i.split("_")[0], i.split("_")[len(i.split("_")) - 1])) for i in gain_recon.index.values]
    loss_recon = pd.read_csv(loss_recon, sep = "\t", index_col = 0)
    loss_recon.index = ["_".join((i.split("_")[0], i.split("_")[len(i.split("_")) - 1])) for i in loss_recon.index.values]
    #prune to target node
    if target_node is not None:
        pt_tree = pt_tree.search_nodes(name = target_node)[0]
    node2name = dict((i.name, i.name) for i in pt_tree.traverse(strategy = 'preorder'))
    pfams_with_event = set()
    pfam2color = {}
    #set the style of the branches and nodes according to the posterior probability
    top10_feats = feats.iloc[:10,]
    for n in pt_tree.traverse():
        ns = node_recon.loc[n.name, phenotype] 
        style = ete2.NodeStyle()
        style["shape"] = 'square'
        style['size'] = 10
        if pd.isnull(ns):
            style['fgcolor'] = 'grey'
        elif ns < threshold:
            style['fgcolor'] = 'darkred'
        else:
            style['fgcolor'] = 'green'
        if not n.name == "N1":
            branch_id = n.name + "_" + n.up.name
            #print gain_recon.loc[branch_id, phenotype], loss_recon.loc[branch_id, phenotype]
            if gain_recon.loc[branch_id, phenotype] > threshold:
                style["hz_line_type"] = 1
                style["hz_line_color"] = 'green' 
                style["hz_line_width"] = 3
            elif loss_recon.loc[branch_id, phenotype] > threshold:
                style["hz_line_type"] = 1
                style["hz_line_color"] = 'red'
                style["hz_line_width"] = 3
            else:
                style["hz_line_type"] = 0
                style["hz_line_color"] = 'black'
            n.set_style(style)
            #set species name instead of tax id
            if n.name in sample_mapping.index:
                node2name[n.name] = sample_mapping.loc[n.name,][0]
            #add majority feature gains and losses
            events = []
            for i in range(top10_feats.shape[0]): 
                cf = faces.CircleFace(radius = 8, style = "circle", color = kelly_colors_hex[i])
                pfam2color[top10_feats.index[i]] = kelly_colors_hex[i]
                #gain events
                #print gain_recon.columns
                if gain_recon.loc[branch_id, top10_feats.index[i]] > threshold:
                    pfam2color[top10_feats.index[i]] = kelly_colors_hex[i]
                    tf = faces.TextFace("+")
                    events.append(tf)
                    pfams_with_event.add(node_recon.index[i])
                    events.append(cf)
                #loss events
                elif loss_recon.loc[branch_id, top10_feats.index[i]] > threshold:
                    pfam2color[top10_feats.index[i]] = kelly_colors_hex[i]
                    tf = faces.TextFace("-")
                    events.append(tf)
                    pfams_with_event.add(node_recon.index[i])
                    events.append(cf)
            for i in range(len(events)):
                n.add_face(events[i], column = i, position = "branch-top")
    for n in pt_tree.traverse():
        if n.name in node2name:
            n.name = node2name[n.name]
    #print top10_feats.loc[:,"Pfam_acc"].values
    #print list(pfams_with_event)
    #filtered_pfams = filter(lambda i: i in list(pfams_with_event), top10_feats.loc[:,"Pfam_acc"].values)
    #print filtered_pfams
    #filtered_ids = pt_gt2id.loc[filtered_pfams, 0] - 1
    #print filtered_ids
    #top10_feats_with_event = top10_feats.loc[filtered_ids,]
    return pt_tree, top10_feats, pfam2color

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("""build node to feat matrix from gainLoss input""")
    parser.add_argument("node_recon", help = "gainLoss node ancestral character state reconstruction") 
    parser.add_argument("gain_recon", help = "gainLoss node ancestral character state reconstruction") 
    parser.add_argument("loss_recon", help = "gainLoss node ancestral character state reconstruction") 
    parser.add_argument("tree", help = "tree with internal nodes labeled") 
    parser.add_argument("pfam_mapping",  help = "feature mapping/list")
    parser.add_argument("feat_list",  help = "list of features")
    parser.add_argument("--target_node", default = "N1", help = "list of features")
    parser.add_argument("phenotype",  help = "target phenotype")
    parser.add_argument("threshold", type = float,  help = "treshold to call genotype/phenotype events")
    parser.add_argument("sample_mapping",  help = "mapping between sample ids and names")
    parser.add_argument("out",  help = "output file")
    a = parser.parse_args()
    pt_tree, feats, pf2color = get_tree(node_recon = a.node_recon, gain_recon = a.gain_recon, loss_recon = a.loss_recon, pfam_mapping = a.pfam_mapping, tree = a.tree, feat_list = a.feat_list, phenotype = a.phenotype, target_node = a.target_node, threshold = a.threshold, sample_mapping = a.sample_mapping)
    plot_tree(pt_tree, a.target_node, a.out)
    plot_legend(feats, a.out, pf2color) 
