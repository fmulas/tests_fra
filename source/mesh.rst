mesh_enrichment.py
==============================
Loads MeSH terms from NCBI papers and performs MeSH-based gene set enrichment.


.. autofunction:: mesh_enrichment.__counts2tfidf__
.. autofunction:: mesh_enrichment.update_mesh
.. autofunction:: mesh_enrichment.reverse_annot
.. autofunction:: mesh_enrichment.preprocessing
.. autofunction:: mesh_enrichment.main


**List of parameters**

data_folder: where the lists of genes are, can contain an entire path

filelist: names of the input files separated by comma with no spaces, or all if entire folder

output_folder: name of the results folder

output_file: name of the resulting file

update: TRUE if the MeSH annotation needs update


**Examples**

*Python call*

# a) Using the parameters in the parameter file

python mesh_enrichment.py -p mesh_parameters_template.txt

# b) Updating some parameters in the parameter file

python mesh_enrichment.py -p mesh_parameters_template.txt --data_folder myfolder

*Calling this script from R*

The R script mesh_xls.R provides an interface with R to both MeSH enrichment and xls writing.
The function mesh_enrich calls the Python script mesh_enrichment. An example of use follows:

source('mesh_xls.R')

# a) Using the parameters in the parameter file

mesh_enrich(fname, par_file='mesh_parameters_template.txt')

# b) Updating some parameters in the parameter file

# folder containing the gene lists

fname <- 'genelists'

# list of files to be considered in the folder

flist <- 'genelist1.txt, genelist2.txt'

folder_out <- 'res_mesh_gl'

mesh_enrich(fname, listf= flist, output='mesh_enriched', output_f=folder_out, update='FALSE')

**Notes**

All needed python packages are available on scc, after loading epd/epd-2.7.3 (Traditional Enthought Python Distribution)