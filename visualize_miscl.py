import pandas as ps
import numpy as np
import ete2
from ete2 import TreeStyle, faces, AttrFace, CircleFace
import os
import colorsys

def map_bacc2color(val, minval, maxval):
    # convert val in range minval..maxval to the range 0..120 degrees which
    # correspond to the colors red..green in the HSV colorspace
    h = (float(val-minval) / (maxval-minval)) * 120
    # convert hsv color (h,1,1) to its rgb equivalent
    # note: the hsv_to_rgb() function expects h to be in the range 0..1 not 0..360
    return  '#%02x%02x%02x' %tuple(map(lambda x: int(x*255), colorsys.hsv_to_rgb(h/360, 1., 1.)))


def miscl_layout(node):
    if node.is_leaf():
        name_face = AttrFace("miscl_full")
        faces.add_face_to_node(name_face, node, column=0, position="branch-right")
    else:
        name_face = AttrFace("miscl_full", fsize = 10)
        faces.add_face_to_node(name_face, node, column=0, position="branch-right")
    if not ps.isnull(miscl_phypat_ml.loc[node.name, '0']): 
        c_mp = CircleFace( radius = 10, style="sphere", color = map_bacc2color(miscl_phypat_ml.loc[node.name, '0'], bacc_min, bacc_max))
        c_p = CircleFace( radius = 10, style="sphere", color = map_bacc2color(miscl_phypat.loc[node.name, '0'], bacc_min, bacc_max))
    else:
        c_mp = CircleFace( radius = 10, style="sphere", color = 'grey')
        c_p = CircleFace( radius = 10, style="sphere", color = 'grey')
    c_p.opacity = 0.7
    c_mp.opacity = 0.7
 
    faces.add_face_to_node(c_mp, node, column=1, position="float")
    faces.add_face_to_node(c_p, node, column=2, position="float")

if __name__ == "__main__":
    import argparse
    FULL = False
    parser = argparse.ArgumentParser("visualize ancestral classification performance")
    parser.add_argument("phypat", help='the input directory with the extant misclassified statistics for phypat ')
    parser.add_argument("phypat_ml", help='the input directory with the extant misclassified statistics for phypat_ml')
    parser.add_argument("-f", "--more_info", help = "add statistics about false positives / false negatives to the nodes", default = False, type = bool)
    parser.add_argument("-r", "--rank", help = "specify rank at which to cut the tree", default = "genus")
    parser.add_argument("out_f", help='output image path')
    args = parser.parse_args()
    miscl_phypat = ps.read_csv(args.phypat, index_col = 0, sep = "\t")
    print miscl_phypat
    miscl_phypat_ml = ps.read_csv(args.phypat_ml, index_col = 0, sep = "\t")
    t = ete2.Tree("gideon_tree.nwck")
    t.scientific_name = "Bacteria"
    t.rank = "Kingdom"
    bacc_max = max(miscl_phypat['0'].max(), miscl_phypat_ml['0'].max())
    bacc_min = min(miscl_phypat['0'].min(), miscl_phypat_ml['0'].min())
    #prune tree to genus level
    #get genus nodes
    genera = []
    for i in t.traverse():
        if i.rank == args.rank:
            genera.append(i)
    t.prune(genera)
    #join scientific name and miscl statistics
    for n in t.traverse():
        if args.more_info:
            n.miscl_full = "%s | %i %i %i %i %s | %i %i %i %i %s" % tuple([n.scientific_name] + [round(miscl_phypat_ml.loc[n.name,].iloc[i], 2) for i in range(5)] + [round(miscl_phypat.loc[n.name,].iloc[i], 2) for i in range(5)])
        else:
            print miscl_phypat_ml.columns
            n.miscl_full = "%s | %s | %s" % (n.scientific_name , round(miscl_phypat_ml.loc[n.name,'0'], 2),  round(miscl_phypat.loc[n.name,'0'], 2))

    ts = TreeStyle()
    #ts.rotation = 90
    ts.layout_fn = miscl_layout
    ts.show_leaf_name = False
    t.render(args.out_f, tree_style = ts)    



