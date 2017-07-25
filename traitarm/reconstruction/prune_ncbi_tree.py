from ete2 import TreeNode, Tree
import sys
#read ncbi tree from saved newick file
def prune_ncbi_tree(full_tree_fn, ncbi_ids_fn, format_out, outfile, newick_extended = False,  do_name_internal = False, failed_to_map = None, resolve_polytomy = False):
    t= TreeNode(full_tree_fn,1)
    #read in ncbi ids
    if not ncbi_ids_fn is None:
        tax_id_f = open(ncbi_ids_fn, "r")
        tax_ids =[line.strip() for line in tax_id_f]
        tax_ids_exists = dict((tax_id, False) for tax_id in tax_ids)
        tax_ids_existing = []
        #check if all tax ids are in the ncbi tree
        for node in t.traverse('preorder'): 
            if tax_ids_exists.has_key(node.name):
                tax_ids_exists[node.name] = True 
                tax_ids_existing.append(node)
        if not failed_to_map is None:
            with open(failed_to_map, 'w') as ftm:
                for tax_id in tax_ids_exists:
                    if tax_ids_exists[tax_id]:
                        ftm.write("%s\n"% tax_id)
        t.prune(tax_ids_existing)
    if resolve_polytomy:
        t.resolve_polytomy(recursive = True)
    if do_name_internal:
        n = 2
        t.rank = "bacteria"
        for node in t.traverse(strategy = 'preorder'):
            if node.is_root():
                node.name = "N1"  
            #if node.name in tax_ids:
            #    node.children = []
            elif not node.is_leaf():
               node.name = "N%s" %n
               n += 1
    if newick_extended:
        t.write(features = ["name", "scientific_name", "rank"], outfile=outfile, format=format_out)
    else:
        t.write(outfile=outfile,format=format_out, format_root_node = True)



def convert_ncbi_tree(infile,outfile):
    """convert tree from format 1 to format 8"""
    t = TreeNode(infile,1)
    t.write(features = ["name", "scientific_name", "rank"], outfile=outfile, format=1)
