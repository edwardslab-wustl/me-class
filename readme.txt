********
ME-Class
********
ME-class is a tool to determine genes that show differential methylation changes
that are associated with differential expression changes.  Methylation data is
typically from whole-genome bisulfite sequencing (WGBS) and genome-wide
expression (typically RNA-Seq) data. Methods are provided to encode methylation
data in a variety of ways. However we highly recommend using the default TSS
method.  See Schlosberg et al. for more information.


******************
AUTHOR INFORMATION
******************
Chris Schlosberg chris.schlosberg@gmail.com 
John Edwards  jredwards@wustl.edu


*********************
COPYRIGHT AND LICENSE
*********************
ME-Class is copyrighted by Chris Schlosberg and John Edwards.

ME-Class is distributed under the GPL-3 licence. See license.txt distributed for
this file for the complete license information.


*************
INSTALLATION
*************
No special installation instructions.  Simply run each command as a script. i.e.
python <script-name>

ME-Class has been tested on python 2.7 and 3.5. 


************
DEPENDENCIES 
************
numpy (1.11.1)
scipy (0.17.1)
matplotlib (1.5.1)
sklearn (0.17.1)
mlpy (3.5.0)

Note: To simply this process, we highly recommend installing anaconda (4.1.1) or a
similar distribution which contains all the required packages for ME-class.


****
DEMO
****
(1) Run interpolation on included example files.
	
    python methylation_interpolation.py TSS -f TRUE -g ../example_dataset/refGene.hg19.21_12_15.txt ../example_dataset/E071_E079.bedgraph ../example_dataset/E071_E079.expr E071_E079
	python methylation_interpolation.py TSS -f TRUE -g ../example_dataset/refGene.hg19.21_12_15.txt ../example_dataset/E094_E095.bedgraph ../example_dataset/E094_E095.expr E094_E095
	python methylation_interpolation.py TSS -f TRUE -g ../example_dataset/refGene.hg19.21_12_15.txt ../example_dataset/E096_E097.bedgraph ../example_dataset/E096_E097.expr E096_E097
	python methylation_interpolation.py TSS -f TRUE -g ../example_dataset/refGene.hg19.21_12_15.txt ../example_dataset/E095_E096.bedgraph ../example_dataset/E095_E096.expr E095_E096

(2) Run classification module to perform leave-one-sample-out model for classification.
    python expression_classification.py TSS evaluate_samples.txt



***********
QUICK START
***********

To get started with ME-Class we recommend first running the demo dataset above.


(I) Input Files
     (a) Reference gene coordinate file.  Reference file for genes to be analyzed
         downloaded from the UCSC Genome Browser.

     (b) Methylation data for individual samples. Each sample comparison should
         be in a single file in standard tab-delimited bedgraph format: 
         <chr> <pos1> <pos2> <differential methylation level>
 
     (c) Expression data for each sample.  Each sample comparison should be in a
         single file in tab-delimited format:
         <gene_id>  <sample1 expression>  <sample2 expression>


(II) Run interpolations to prepare data for classification.
     - Use the script "methylation_interpolation.py".  
     - We have found the "TSS" method works the best. See Schlosberg et al. for
       more details about performance of different models. By default samples
       are scaled and normalized.

(III) Perform classification.
    - Resulting prediction files will consist of prediction of expression change 
        along with gene information.

************
USEFUL NOTES
************

(1) The interpolator prepares data for the resulting classification.  There are
several versions of the interpolator that match the different models used in the
interpolator.  Type "python methylation_interpolation.py -h" to see available
options.  There are specific options available for each model.  To see type
"python methylation_interpolation.py <model> -h".

(2) All model-specific options must be placed after the model name and before
the positional arguments (methylation and expression data).  E.g. python
methylation_interpolation.py <model> <model_specific_options> <bedgraph> <expr>
<tag> 






