
#-------------------------------------------------------------------------------
# Name:        mesh_enrichment
# Purpose:     save mesh terms and performs gene set enrichment analysis
#
# Author:      Francesca
#
# Created:     30/04/2014

#-------------------------------------------------------------------------------
#!/usr/bin/env python

"""
Usage: python mesh_enrichment.py -p parameter_file.txt
-h help
-p parameter file				*[No default value]
"""


import urllib
import xml.dom.minidom
import json
import scipy.stats as sc
import sys
import pickle
import os
import utils

def __counts2tfidf__(g2a):
    """ Computes Term Frequency - Inverse Document Frequency, which reflects the specific importance of a MeSH term for a gene
    :param g2a: dictionary gene terms.
    """
    global tcounts, df, idf
    # term document frequency
    terms = set([])
    for ts in g2a.values():
        terms.update(set(ts.keys()))
    # (document frequency); for a specific term, how many genes are annotated to it
    df = dict([(t, sum([True for g, ts in g2a.items() if t in ts])) for t in terms])
    # gene term counts; total number of counts for each gene
    tcounts = dict([(g, sum([c for _, c in ts.items()])) for g, ts in g2a.items()])
    # inverse document frequency
    n = float(len(g2a))
    idf = dict([(t, log(n / df[t])) for t in terms])
    return dict([ (g, dict([(t, float(c)/tcounts[g] * idf[t]) for t,c in ts.items()]) ) for g, ts in g2a.items()])

def normalize(tfidf):
    """ Normalize vectors (val:count) in the dictionary such that vectors have unit length """
    n = {}
    for g, td in tfidf.items():
        norm = sqrt(sum(map(lambda x: x*x, td.values())))
        n[g] = dict([(t, v/norm) for t, v in td.items()])
    return n

def update_mesh():
    """ Updates the json file containing gene-MeSH annotations: # for each gene, saves the associated MeSH, and for each MeSH, the number of reference papers
    Current MeSH terms and frequencies are downloaded for the whole genome from http://gene2mesh.ncibi.org/ """
    path = os.getcwd()
    filein = open(os.path.join(path, 'genes_human_types.txt'))
    genes = [line.rstrip().split()[0] for line in filein.readlines()[1:] if line.rstrip().split()[1] == 'protein-coding']

    gene2mesh = {}
    g2mpart = {}
    for k, gs in enumerate(genes_cod):
        url = "http://gene2mesh.ncibi.org/fetch?genesymbol=%s&taxid=9606&limit=10000" %gs
        data = urllib.urlopen(url).read()
        dom = xml.dom.minidom.parseString(data)
        res_xml = dom.getElementsByTagName("Result")

        for i, el in enumerate(res_xml):
            h = el.getElementsByTagName("MeSH") # mesh
            j = h[0].childNodes[0].getElementsByTagName("Name")
            name = str(j[0].childNodes[0].data)
            p = el.getElementsByTagName("DocumentSet")
            np = len(p[0].childNodes) # number of papers
            dictg = gene2mesh.get(gs, {}) # save
            dictg[name] = np
            gene2mesh[gs] = dictg
            g2mpart[gs] = dictg

    # save dictionary
    path = os.getcwd()
    path_cbm = '/restricted/projectnb/montilab-p/CBMrepositoryData/annot'
    res_folder = os.path.join(path_cbm, "mesh_annotation")
    if not(os.path.exists(res_folder)):
        os.makedirs(res_folder)
    outfile = '%s/%s' %(res_folder, "gene2mesh_human")
    json.dump(gene2mesh, open(outfile, 'w'))

    norm = normalize(gene2mesh)
    norm_final = preprocessing(norm)
    # save normalized dict
    outfile_n = '%s/%s' %(res_folder, "gene2mesh_human_normalized")
    json.dump(norm_final, open(outfile_n, 'w'))

    return norm_final


def normalize(annot):
    tfidf = counts2tfidf(annot)
    norm = normalize(tfidf)
    outfile = '%s/%s' %(res_folder, "gene2mesh_human")
    json.dump(data, open(outfile, 'w'))
    return norm


def reverse_annot(anndict):
    """ Returns a dictionary with MeSH terms as keys and genes as values """
    # inverse dictionary (terms: genes: values)
    anninv = {}
    for kg, tvals in anndict.items():
        for tv, n in tvals.items():
            dict_t = anninv.get(tv, {})
            dict_t.update({kg: n})
            anninv[tv] = dict_t

    # save dict
    path = os.getcwd()
    res_folder = os.path.join(path, "mesh_annotation")
    outfile = '%s/%s' %(res_folder, "mesh2genes_human_normalized")
    json.dump(anninv, open(outfile, 'w'))
    return anninv


def preprocessing():
    """ Removes very general annotation terms, i.e. with frequency above the 99.5th percentile of the distribution of frequencies """
    # remove too general terms
    all_length = [(tt, len(n_inv[tt])) for tt in n_inv.items()]
    perc = sc.scoreatpercentile([v[1] for v in all_length], 99.5)
    keys_general = [v[0] for v in all_length if v[1] > perc]

    n_inv_ok = dict([(k, v) for k, v in n_inv.items() if k not in keys_general])
    n_ok = {}
    for g, td in norm.items():
        new_n = []
        for tname, tvalue in td.items():
            if tname not in keys_general:
                new_n.append((tname, tvalue))
        n_ok[g] = dict(new_n)


def main():
    """ Computes MeSH-based enrichment. Input folder with gene lists are obtained from the parameter file. For each list, the output is a ranked set of MeSH terms, enrichment pvalues, FDR and sum of TF-IDF values (json format)
  """

    list_args = sys.argv[1:]
    print list_args
    path = os.getcwd()

    if '-h' in list_args:
        print __doc__
        sys.exit(0)
    elif len(list_args) <2 or '-p' not in list_args:
        print __doc__
        sys.exit(0)
    else:
        params, parameter_file, updated = utils.update_parameters(list_args)
        folder = params['data_folder']
        filelist_str = params['filelist']
        folder_res = params['output_folder']
        print params

        if params['update'] == 'TRUE':
            norm = update_mesh()
            n_inv = reverse_dict(norm)
        else:
            res_folder = os.path.join(path, "mesh_annotation")
            json_data_norm = open('%s/%s' %(res_folder, "gene2mesh_human_normalized")).read()
            norm = json.loads(json_data_norm)
            reference = norm.keys()

            json_inv_norm = open('%s/%s' %(res_folder, "mesh2genes_human_normalized")).read()
            n_inv = json.loads(json_inv_norm)

        if filelist_str == 'all':
           filelist = os.listdir(folder)
        else:
             filelist = filelist_str.split(',')
        data = {}
        for filename in filelist:
            genes = [l.rstrip() for l in open(os.path.join(folder, filename)).readlines()]
            print len(genes)
            # genes in the list
            genesok = [gg for gg in genes if gg in norm.keys()]
            # save all terms
            terms = []
            for g in genesok:
                terms.extend(norm[g].keys())
            terms = list(set(terms))

            res = []
            for t in terms:
                # select terms (genesets) with at least 3 genes
                if len(n_inv[t]) > 3:
                    vals_t = [(v, k) for k,v in n_inv[t].items()]

                    # genes in the genome associated to this term
                    list_genes = n_inv[t].keys()
                    # our genes associated to this term
                    mappedGenes = list(set(list_genes).intersection(genesok))
                    pv = 1-sc.hypergeom.cdf(len(mappedGenes), len(reference), len(list_genes), len(genesok))
                    vals_t.sort(reverse = True)

                    sumrank = sum([float(i+1)/len(vals_t) for i,el in enumerate(vals_t) if el[1] in genesok])
                    res.append((pv, sumrank, len(mappedGenes), len(list_genes), t, ' '.join(mappedGenes)))

            res.sort()
            pv_corr = utils.FDR([pp[0] for pp in res])
            res_corr = [(p, res[i][0], 1.0/res[i][1], 1.0/res[i][2], 1.0/res[i][3], res[i][4], res[i][5]) for i, p in enumerate(pv_corr)]
            res_corr.sort()
            res_data = [['MeSH term', 'FDR', 'p-value', 'ranksum','count', 'reference', 'gene list']]

            for r in res_corr:
                res_data.append([r[5], r[0], r[1], 1.0/r[2],  1.0/r[3],  1.0/r[4], r[6]])
            data[filename[:-4]] = res_data


        res_folder = os.path.join(path, folder_res)
        if not(os.path.exists(res_folder)):
            os.makedirs(res_folder)
        outfile = '%s/%s' %(res_folder, params['output_file'])
        json.dump(data, open(outfile, 'w'))
        print "results saved in folder %s" %res_folder

if __name__ == '__main__':

    main()







