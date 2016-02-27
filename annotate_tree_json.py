# vim: fdm=indent
'''
author:     Fabio Zanini
date:       27/02/16
content:    Combine a bare newick tree with an annotation table into a JSON
            annotated tree.
'''
def annotate_and_convert(tree_file, tree_format,
                         annotation_file,
                         output_file):
    '''Annotate tree and convert it to JSON
    
    Arguments:
       tree_file (string): file name of the input tree
       tree_format (string or None): format of the input tree (e.g. 'newick')
       annotation_file (string): file name of the input annotation table
       output_file (string): output file name
    '''

    def combine_tree_annotations(tree_file, tree_format, annotation_file):
        '''Combine a bare tree with annotations'''
        import pandas as pd
        from Bio import Phylo
    
        if tree_format is None:
            tree_format = tree_file.split('.')[-1]
    
        tree = Phylo.read(tree_file, tree_format)
        anno = pd.read_csv(annotation_file, sep='\t', index_col='name')
    
        metadata = list(anno.columns)
        for node in tree.find_elements():
            if node.name in anno.index:
                for key, value in anno.loc[node.name].items():
                    setattr(node, key, value)
        
        return {'tree': tree, 'metadata_nodes': metadata}
    
    
    def write_tree_to_json(tree, output_file,
                           metadata_tree=[],
                           metadata_nodes=[],
                           children_attrname="children"):
        '''Write Biopython tree to JSON'''
    
        tree_dict = {}
        metadata_tree = [m for m in metadata_tree if m not in ['tree']]
        for field in metadata_tree:
            if hasattr(tree, field):
                tree_dict[field] = getattr(tree, field)
    
    
        not_metadata = ['clades', children_attrname]
        metadata_nodes = [m for m in metadata_nodes if m not in not_metadata]
        def convert_to_dict(node):
            d = {}
            for field in metadata_nodes:
                if hasattr(node, field):
                    d[field] = getattr(node, field)
            d[children_attrname] = [convert_to_dict(c) for c in node.clades]
            return d
    
        tree_dict['tree'] = convert_to_dict(tree.root)
    
        import json
        with open(output_file, 'w') as handle:
            json.dump(tree_dict, handle, indent=1)


    data = combine_tree_annotations(tree_file, tree_format, annotation_file)
    write_tree_to_json(data['tree'], output_file,
                       metadata_tree=[],
                       metadata_nodes=data['metadata_nodes'])


# Script
if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='Annotate tree into JSON')
    parser.add_argument('treefile',
                        help='Input filename with the tree')
    parser.add_argument('annotationfile',
                        help='Input filename with the annotations')
    parser.add_argument('outputfile',
                        help='Output filename (JSON)')
    parser.add_argument('--treeformat', default=None,
                        help='Format of the input tree file')

    args = parser.parse_args()

    annotate_and_convert(args.treefile, args.treeformat,
                         args.annotationfile, args.outputfile)
