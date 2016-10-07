### Custom container Class definitions

class Window:
    def __init__(self,low,high,bins,bp_window_size=None):
        self.low = low
        self.high = high
        self.len = self.high-self.low
        self.bins = bins
#        print "Window Bins: ",self.bins
        self.x = list()
        self.x_raw_meth = list()
        self.y_meth = list()
        self.y_cpg = list()
        self.y_raw_meth = list()
        if bp_window_size == None:        
            self.set_bp_window_size()        
        else: 
            self.bp_window_size = bp_window_size   
#        print "BP Window Size: ",self.bp_window_size
            
    def set_counter(self,counter):
        self.counter = counter
    
    def set_bp_window_size(self):
        self.bp_window_size = (self.len)/self.bins

    def scale_sample_x_coords(self):
        self.scale_sample_x = [((x-self.low)/float(self.len))+self.counter for x in self.sample_x]
        
    def scale_raw_x_coords(self):
        self.scale_x_raw_meth = [((x-self.low)/float(self.len))+self.counter for x in self.x_raw_meth]

    def sample_window(self):
        if self.bp_window_size != 0:
            self.downsample_window()
        else:
            self.sample_x = list()
            self.sample_y_meth = list()
            self.sample_y_cpg = list()

    def downsample_window(self):
        self.sample_x = downsample(self.x,self.bp_window_size,self.bins)
#        print "Downsample: Number of bins in sample_x: ",len(self.sample_x)
        self.sample_y_meth = downsample(self.y_meth,self.bp_window_size,self.bins)
        self.sample_y_cpg = downsample(self.y_cpg,self.bp_window_size,self.bins)

class Gene:
    def __init__(self,gene_id):
        ### Either the REFSEQ or ENSG ID        
        self.gene_id = gene_id
        ### List of all the transcripts
        self.transcripts = list()
        
    def set_location_information(self,chrom,strand,pos_low,pos_high,exons,common_name):        
        self.chrom = chrom
        self.strand = strand
        self.pos_low = pos_low
        self.pos_high = pos_high
        self.exons = exons
        self.common_name = common_name

    def set_half_tss_window(self,half_tss_window):
        self.half_tss_window = half_tss_window

    def set_expr(self,expr):
        self.expr = expr    
    
    def set_sample(self,sample):
        self.sample = sample
        
    def append_transcript(self,transcript):
        self.transcripts.append(transcript)

    def set_introns(self):
        self.introns = self.calculate_introns()        
        
    def calculate_introns(self):
        introns = list()    
        for i,exon_tup in enumerate(self.exons):
            try:
                i_tup = (exon_tup[1]+1,self.exons[i+1][0]-1) 
                introns.append(i_tup)
            except IndexError:
                continue
        return introns        
        
    def merge_gene_transcripts(self):
        common_names = list()        
        for i,t in enumerate(self.transcripts):
            common_names.append(t.common_name)            
            if i == 0:
                merged_exons = merge_intervals(t.exons)
            else:
                current_merged_exons = merged_exons+t.exons
                merged_exons = merge_intervals(current_merged_exons)                              
        pos_low = merged_exons[0][0]
        pos_high = merged_exons[-1][1]
        chrom = self.transcripts[0].chrom
        strand = self.transcripts[0].strand
        common_name = get_most_common_element(common_names)
        self.set_location_information(chrom,strand,pos_low,pos_high,merged_exons,common_name)
        self.set_introns()
        return
                        
class Transcript:
    def __init__(self,transcript_id,chrom,strand,pos_low,pos_high,exons,common_name,cds_tx_start_cmpl,cds_tx_end_cmpl):
        self.transcript_id = transcript_id
        self.chrom = chrom
        self.strand = strand
        self.pos_low = pos_low
        self.pos_high = pos_high
        self.exons = exons
        self.common_name = common_name
        self.cds_tx_start_cmpl = cds_tx_start_cmpl
        self.cds_tx_end_cmpl = cds_tx_end_cmpl
        
class FilterFlags:
    def __init__(self,annot_filetype="",only_simple_chrs=True,only_mrna=True,only_complete=False,only_uniq_tx_start=False,differential=False,collapse_on_common_name=True,exclude_lt_four_exons=False):
        ### Gene Annotation Format 
        self.annot_filetype = annot_filetype
        ### Transcript Level Filters        
        ## Only look at transcripts defined as complete. cds start and end well defined.
        self.only_complete = only_complete
        ### Gene Level Filters
        ## Only look at genes on chromosomes defined in valid_chr    
        self.only_simple_chrs = only_simple_chrs
        ## Only look at genes defined as protein coding (starting with NM_)
        self.only_mrna = only_mrna
        ## Only look at genes with unique TxStart sites
        self.only_uniq_tx_start = only_uniq_tx_start
        ## Identify that this interpolation is for differential
        self.differential = differential
        ## Collapse Annotations so that only common gene name represented
        self.collapse_on_common_name = collapse_on_common_name
        ## Exclude genes with < 4 exons to have straightforward comparison with ROI
        self.exclude_lt_four_exons=exclude_lt_four_exons


class FilterCounts:
    def __init__(self):
        ## Total number of entries        
        self.num_transcripts = 0
        ## Number of genes/transcripts defined as complete. cds start and end well defined.
        self.num_tx_complete = 0
        ## Number of genes/transcripts defined as incomplete. cds start and end ill defined.
        self.num_tx_incomplete = 0
        ## Number of unique gene ids
        self.num_uniq_gene_ids = 0 
        ## Number of genes ids with no transcripts after transcript filters
        self.num_gene_ids_no_tx_post_tx_filters = 0      
        ## Number of unique gene ids
        self.num_gene_ids_post_tx_filters = 0
        ## Number of genes on chromosomes defined in valid_chr    
        self.num_genes_simple_chrs = 0
        ## number of genes on chromosomes not defined in valid_chr    
        self.num_genes_other_chrs = 0
        ## Number of genes defined as protein coding (starting with NM_)
        self.num_genes_mrna = 0
         ## Number of genes defined as not protein coding (not starting with NM_)
        self.num_genes_notPC = 0
        ## Number of genes with unique TxStart sites
        self.num_genes_uniq_tx_start = 0
        ## Number of genes with unique TxStart sites
        self.num_genes_nonuniq_tx_start = 0
        ## Number of genes loaded into database for processing
        self.num_final_loaded_genes = 0
        
    def write_output(self,filter_flags):
        sys.stderr.write("Transcript Processing:\n")
        sys.stderr.write("\tIdentified %d transcripts\n" % (self.num_transcripts))
        if filter_flags.only_complete:
            sys.stderr.write("\tLoaded %d transcripts with well-defined CDS start and end\n" % (self.num_tx_complete))
            sys.stderr.write("\tLoaded %d transcripts with ill-defined CDS start and end\n" % (self.num_tx_incomplete))
        sys.stderr.write("Gene Processing:\n")        
        sys.stderr.write("\tIdentified %d unique gene annotations\n" % (self.num_uniq_gene_ids))
        if filter_flags.only_complete:
            sys.stderr.write("\tIdentified %d genes with no valid transcripts, post-transcript filter(s)\n" % (self.num_gene_ids_no_tx_post_tx_filters))
            sys.stderr.write("\tIdentified %d genes post-transcript filter(s)\n" % (self.num_gene_ids_post_tx_filters))
        if filter_flags.only_simple_chrs:
            sys.stderr.write("\tLoaded genes %d on simple chromosomes\n" % (self.num_genes_simple_chrs))        
            sys.stderr.write("\tLoaded genes %d on other chromosomes\n" % (self.num_genes_other_chrs))
        if filter_flags.only_mrna:
            sys.stderr.write("\tLoaded %d protein coding genes\n" % (self.num_genes_mrna))
            sys.stderr.write("\tLoaded %d non-protein coding genes\n" % (self.num_genes_notPC))
        if filter_flags.only_uniq_tx_start:
            sys.stderr.write("\tLoaded %d genes with unique Tx start site\n" % (self.num_genes_uniq_tx_start))
            sys.stderr.write("\tLoaded %d genes with non-unique Tx start site\n" % (self.num_genes_nonuniq_tx_start))
        sys.stderr.write("Loaded %d genes for cross-referencing\n" % (self.num_final_loaded_genes))
    
def main():
    ### For help call: python <script>.py -h
    
    ### Define valid chromosomes
    global valid_chr
    valid_chr = set(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'])
    ### Set up parser arguments
    global parser
    parser,subparsers = setup_parser_arguments()
    ### Grab argument variables
    args = parser.parse_args()
    ### Load gene annotations
    gene_annot,filter_flags = load_gene_annot(args.annot_file,args)
    ### Load sample files 
    sample_files = load_sample_info(args.bedgraph,args.expr,args.tag)    
    ### Load gene expression
    gene_list,differential = load_gene_expression_rpkm(sample_files,gene_annot)
    ### Release memory from holding ENSG database
    del gene_annot      
    ### Print Label List
    label_list_fh = print_label_list_header(gene_list,args.tag) 
    if args.interp_cpg_density:    
        ### Load up all CpG positions
        cpg_density = dict()
        cpg_density = load_cpg_perc_ref(args.cpg_density_file)   
    else:
        cpg_density = None
    if args.interp_type=="DMR":
        process_dmr(sample_files,gene_list,label_list_fh,differential,args)
    else:
        ### Calculate and print interpolation for each gene
        interpolate_methylation_cpgdensity_windows(sample_files,gene_list,cpg_density,label_list_fh,differential,args)    

### Function Definitions below here

### DMR Functions

def process_dmr(sample_files,gene_list,label_list_fh,differential,args):
    interp_type_tag, interp_type_header = setup_interp_type(args)
    meth_data_file = ".".join([args.tag,interp_type_tag,"meth","dat"])
    if os.path.isfile(meth_data_file):
        os.remove(meth_data_file)
    meth_fh = open(meth_data_file,'a')    
    meth_fh.write(interp_type_header)        
    gene_count = 0
    ### For each gene in each sample, interpolate methylation over whole gene
    for sample,val in sample_files.items():
        ### Load provided methylation file for given sample
        sample_meth = load_methylation_dmr(val[1])
        sys.stderr.write("Interpolating genes for %s\n" % (sample))
        ### Iterate through each gene to get methylation values over whole gene
        for gene in gene_list:
#            print gene.gene_id,gene.common_name,gene.strand,gene.pos_low,gene.pos_high
            if gene.sample != sample:
                continue
            gene_count +=1
            if gene_count%10==0:
                sys.stderr.write("\tNow on gene %d out of %d\n" % (gene_count,len(gene_list)))
            if gene.strand == "+":
                tss = gene.pos_low
            else:
                tss = gene.pos_high
            if not gene.chrom in sample_meth.keys():
                sys.stderr.write("\tWarning: %s has no DMRs assayed in %s for %s \n" % (gene.gene_id,val[1],gene.chrom))  
                continue  
            dmr_meth_vec = identify_closest_dmr(tss,sample_meth[gene.chrom],gene.strand)
            ### If moved past filters, print the gene            
            print_label_gene(gene,label_list_fh)   
            wl = " ".join([str(s) for s in dmr_meth_vec])+"\n"
            meth_fh.write(wl)
    return

### Interpolation Functions

def interpolate_methylation_cpgdensity_windows(sample_files,gene_list,cpg_density,label_list_fh,differential,args):        
    #Input: List of samples and genes
    #Output: Interpolated values for each gene at for defined gene representation.    
    ### Descriptions of the parameters encoded into the file names    
    interp_type_tag, interp_type_header = setup_interp_type(args)
    
    if args.meth_cpg_density_score:
        meth_cpg_density_data_file = ".".join([args.tag,interp_type_tag,"mcs","dat"])
        if os.path.isfile(meth_cpg_density_data_file):
            os.remove(meth_cpg_density_data_file)
        plot_dir_tag = "_".join([args.tag,interp_type_tag,"meth","cpg","score"])
    else:
        if args.interp_cpg_density:
            cpg_density_data_file = ".".join([args.tag,interp_type_tag,"cpg","dat"])
            if os.path.isfile(cpg_density_data_file):
                os.remove(cpg_density_data_file)
        meth_data_file = ".".join([args.tag,interp_type_tag,"meth","dat"])
        if os.path.isfile(meth_data_file):
            os.remove(meth_data_file)
        if args.interp_cpg_density:
            plot_dir_tag = "_".join([args.tag,interp_type_tag,"meth","cpg"])
        else:
            plot_dir_tag = "_".join([args.tag,interp_type_tag,"meth"])
    if args.plot:
        output_dir = plot_dir_tag+"_interp_curves"
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
    
    ### Output file handling. 
    ### Output files will exist with only a non-informative header 
    ### and whatever bins exist for that given gene.
    ### for TSS and WSG: this will always be the same number for each gene
    ### for SGF: this will vary for each gene, and will need to use the label file 
    ### and the actual data file to create a model for each gene. 
    ### Doing it this way makes it easy to concatenate files too. 

    ### Throw the parameters in the comment header at the top of the file    
    if args.meth_cpg_density_score:
        meth_cpg_fh = open(meth_cpg_density_data_file,'a')
        meth_cpg_fh.write(interp_type_header)
    else:
        meth_fh = open(meth_data_file,'a')    
        meth_fh.write(interp_type_header)
        if args.interp_cpg_density:
            cpg_density_fh = open(cpg_density_data_file,'a')
            cpg_density_fh.write(interp_type_header)

    gene_count = 0
    ### For each gene in each sample, interpolate methylation over whole gene
    for sample,val in sample_files.items():
        ### Load provided methylation file for given sample
        sample_meth = load_methylation(val[1])
        sys.stderr.write("Interpolating genes for %s\n" % (sample))
        
        ### Need to change this loop so that we're iterating once through the bedgraph file
        ### Loading the data, and then interpolating         
                
        ### Iterate through each gene to get methylation values over whole gene
        for gene in gene_list:
#            print gene.gene_id,gene.common_name,gene.strand,gene.pos_low,gene.pos_high
            if gene.sample != sample:
                continue
            gene_count +=1
            if gene_count%10==0:
                sys.stderr.write("\tNow on gene %d out of %d\n" % (gene_count,len(gene_list)))
            ### Define gene boundaries for interpolation of full gene
            ### then subset using windows 
            ### TSS option is the only time where this idea is broken, so base it off of the tss window for it
            lower, upper, interp_window_start, interp_window_stop = define_gene_boundaries(args,gene)

#            time1 = time.time()
#            print "Meth index before: %d" % (sample_meth[gene.chrom][0])
            x_raw_meth,y_raw_meth,start_index = get_position_value_from_range(lower,upper,sample_meth[gene.chrom])
            
            ### Minimum methylation assay filter over all CpGs      
            if len(x_raw_meth) < args.min_meth:
                sys.stderr.write("\tWarning: %s has < %d CpGs assayed in %s\n" % (gene.gene_id,args.min_meth,val[1]))  
                continue            
            
            ### need to identify genes where there are no anchor points.
            if identify_missing_anchors(x_raw_meth,interp_window_start,interp_window_stop,lower,upper):
                sys.stderr.write("\tWarning: %s has no measured CpGs in anchor window(s)\n" %(gene.gene_id))
                continue
                
            ### Need to inform where to start on the chromosome to save time interpolating          
            sample_meth[gene.chrom][0] = start_index            
#            print "Meth index after: %d" % (sample_meth[gene.chrom][0])
#            time2 = time.time()  
#            print '%s took %0.3f ms' % ("methylation lookup", (time2-time1)*1000.0)


            x = range(int(lower),int(upper)+1,1) 
            if args.flankNorm:
#                y_raw_meth = flank_zscore_normalize_values(x_raw_meth,y_raw_meth,x,interp_window_start,interp_window_stop,lower,upper,os.path.join(output_dir,sample+"_"+gene.gene_id))
                y_raw_meth = flank_zscore_normalize_values(x_raw_meth,y_raw_meth,x,interp_window_start,interp_window_stop,lower,upper)
            interp_y_meth = interpolate_defined_window(x_raw_meth,y_raw_meth,x,args.sigma,interp_window_start,interp_window_stop)
            if not args.flankNorm:            
                interp_y_meth = methylation_sanity_filter(interp_y_meth,differential)            
            ### CpG Density Interpolation Section
            if args.interp_cpg_density:            
                ### For each gene, get CpG Density values
#                time1 = time.time()
#                print "CpG density index before: %d" % (cpg_density[gene.chrom][0])
                raw_pos_cpg,raw_val_cpg,start_index = get_position_value_from_range(lower,upper,cpg_density[gene.chrom])   
                ### Need to inform where to start on the chromosome to save time interpolating          
                cpg_density[gene.chrom][0] = start_index
#                print "CpG density index after: %d" % (cpg_density[gene.chrom][0])                
#                time2 = time.time()
#                print '%s took %0.3f ms' % ("CpG density lookup", (time2-time1)*1000.0)
                ### Interpolate the CpG Density            
#                if args.flankNorm:
#                    raw_val_cpg = flank_zscore_normalize_values(raw_pos_cpg,raw_val_cpg,x,interp_window_start,interp_window_stop,lower,upper,sample+"_"+gene.gene_id)
                interp_y_cpg = interpolate_defined_window(raw_pos_cpg,raw_val_cpg,x,args.sigma,interp_window_start,interp_window_stop)                                
                interp_y_cpg = methylation_sanity_filter(interp_y_cpg,False)
            
#            time1 = time.time()
            ### Need to work in physical position space and normalize afterwards
            
            windows = setup_windows(args,gene)            
            
            ### This filter is checking to make sure that there are no genes that
            ### are shorter than the number of bins in any given window.            
            if args.interp_type == "WSG":
                short_flag = False
                for w in windows:
                    if w.bp_window_size == 0:
                        short_flag = True
                if short_flag:
                    sys.stderr.write("\tWarning: %s has < %d bp (# WSG bins)\n" % (gene.gene_id,args.wsg_bins))  
                    continue
            ### This filter is to check to make sure that there are no genes
            ### that are shorter than the downstream window boundary 
            if args.interp_type == "TSS":
                ### only 1 window for TSS                
                w = windows[0]
                if gene.strand == "+":
                    if w.high > gene.pos_high:
                        sys.stderr.write("\tWarning: %s has downstream window boundary > gene length\n" % (gene.gene_id))  
                        continue                           
                else: 
                    if w.low < gene.pos_low:
                        sys.stderr.write("\tWarning: %s has downstream window boundary > gene length\n" % (gene.gene_id))  
                        continue       
                
            ### Scale and paste windows together
            ### Currently, there are only contiguous windows, so just iterate through window list. 
            ### Average every bin bp
            ### Having output plots relative to both input base pairs and normalized so that there are even spacings for each window.
            ### Each window will be represented as 1 for the scaled plots.
            ### Will output both a scaled plot and a non-scaled (actual) plot
            if args.interp_cpg_density:
                x_scale_windows,x_windows_raw_meth,y_windows_raw_meth,\
                x_scale_windows_raw_meth,x_windows,y_windows_interp_meth,\
                y_windows_interp_cpg = sample_scale_windows(args,\
                windows,x,x_raw_meth,y_raw_meth,interp_y_meth,interp_y_cpg=interp_y_cpg)
            else:
                x_scale_windows,x_windows_raw_meth,y_windows_raw_meth,\
                x_scale_windows_raw_meth,x_windows,y_windows_interp_meth \
                = sample_scale_windows(args,windows,x,\
                x_raw_meth,y_raw_meth,interp_y_meth)
                
                
            ### Minimum methylation assay filter over all windowed regions      
            if len(x_windows_raw_meth) < args.min_meth:
                sys.stderr.write("\tWarning: %s has < %d CpGs assayed in %s\n" % (gene.gene_id,args.min_meth,val[1]))  
                continue
            ### Minimum methylation change filter (differential samples only)
            if not args.flankNorm and differential and max([math.fabs(z) for z in y_windows_raw_meth])<=args.min_meth_change:
                sys.stderr.write("\tWarning: %s has < %f methylation change in %s\n" % (gene.gene_id,args.min_meth_change,val[1]))
                continue
            ### If moved past filters, print the gene            
            print_label_gene(gene,label_list_fh)            
            
#            time2 = time.time()
#            print '%s took %0.3f ms' % ("actual interpolation/windowing", (time2-time1)*1000.0)
            
            #### Plotting section
            if args.plot:
                if args.interp_cpg_density:
                    plot_wrapper(differential,output_dir,gene,args,x_windows,x_scale_windows,x_scale_windows_raw_meth,y_windows_interp_meth,x_windows_raw_meth,y_windows_raw_meth,y_windows_interp_cpg=y_windows_interp_cpg)
                else:
                    plot_wrapper(differential,output_dir,gene,args,x_windows,x_scale_windows,x_scale_windows_raw_meth,y_windows_interp_meth,x_windows_raw_meth,y_windows_raw_meth)
            
            ### Need to flip each level of information according to strand to print         
            y_windows_interp_meth = flip_window(y_windows_interp_meth,gene.strand)
            if args.interp_cpg_density:
                y_windows_interp_cpg = flip_window(y_windows_interp_cpg,gene.strand)                        
            
            ### Number of windows for interpolation:
            ### For TSS interpolation: 1
            ### For WSG interpolation: 1(whole scaled gene)+2(proximal windows) = 3
            ### For SGF interpolation: 2+(# exons)+(# introns)
                       
            
            if args.interp_type == "ROI":
                if args.interp_cpg_density:
                    output_roi_data(cpg_density_fh,y_windows_interp_cpg,gene)
                output_roi_data(meth_fh,y_windows_interp_meth,gene)
            else:  
                if args.meth_cpg_density_score:                          
                    y_windows_interp_meth_cpg_score = normalize_methylation_cpg_density_wrapper(y_windows_interp_meth,y_windows_interp_cpg)
                    ### Write interpolated Methylation/CpG Density score to data file
                    wl = " ".join([str(s) for s in y_windows_interp_meth_cpg_score])+"\n"
                    meth_cpg_fh.write(wl)
                else:                
                    ### Now output curves to appropriate file:
                    if args.interp_cpg_density:            
                        ### Write interpolated CpG Density to data file
                        wl = " ".join([str(s) for s in y_windows_interp_cpg])+"\n"
        #                print "# bins in interp cgp: ",len(y_windows_interp_cpg)                
                        cpg_density_fh.write(wl)
        #            print "# bins in interp meth: ",len(y_windows_interp_meth)                
                    ### Write methylation to data file
        #            print len(y_windows_interp_meth),gene.pos_high-gene.pos_low 
                    wl = " ".join([str(s) for s in y_windows_interp_meth])+"\n"
                    meth_fh.write(wl)

        ### Free memory so only holding one sample at a time
        del sample_meth

    return    

def identify_missing_anchors(x_raw_meth,interp_window_start,interp_window_stop,lower,upper):
    first_pos = x_raw_meth[0]
    last_pos = x_raw_meth[-1]
    lower_anchor_cpg_exists = False
    upper_anchor_cpg_exists = False    
    if lower <= first_pos <= interp_window_start:
        lower_anchor_cpg_exists = True
    if interp_window_stop <= last_pos <= upper:
        upper_anchor_cpg_exists = True
    if lower_anchor_cpg_exists and upper_anchor_cpg_exists:
        return False
    else:
        return True
    
def output_roi_data(fh,window_vals,gene):
    upstream_vals = window_vals[:5]
    downstream_vals = window_vals[-5:]
    num_exons = len(gene.exons)
    num_introns = len(gene.introns)    
    num_internal_bins = num_exons+num_introns
    if num_internal_bins == 1:
        first_exon = window_vals[5]
        first_intron = "NA"
        internal_exon = "NA"
        internal_intron = "NA"   
        last_intron = "NA"
        last_exon = window_vals[5]
    elif num_internal_bins == 3:
        first_exon = window_vals[5]
        first_intron = window_vals[6]
        internal_exon = "NA"
        internal_intron = "NA"   
        last_intron = window_vals[6]
        last_exon = window_vals[7]
    elif num_internal_bins == 5:
        internal_vals = [x for j,x in enumerate(window_vals) if j>=5 and j<len(window_vals)-5]        
        first_exon = internal_vals[0]
        first_intron = internal_vals[1]
        internal_exon = internal_vals[2]
        internal_intron = "NA"   
        last_intron = internal_vals[3]
        last_exon = internal_vals[4]       
    else:
        internal_vals = [x for j,x in enumerate(window_vals) if j>=5 and j<len(window_vals)-5]        
        first_exon = internal_vals[0]
        first_intron = internal_vals[1]
        internal_internal_vals = [x for j,x in enumerate(internal_vals) if j > 1 and j < len(internal_vals)-2]
        internal_exon = np.mean([x for j,x in enumerate(internal_internal_vals) if j%2==0])
        internal_intron = np.mean([x for j,x in enumerate(internal_internal_vals) if j%2!=0])
        last_intron = internal_vals[-2]
        last_exon = internal_vals[-1] 
    wll = upstream_vals+[first_exon,first_intron,internal_exon,internal_intron,last_intron,last_exon]+downstream_vals  
    wl = " ".join([str(s) for s in wll])+"\n"
    fh.write(wl)
    return
        

def sample_scale_windows(args,windows,x,x_raw_meth,y_raw_meth,interp_y_meth,interp_y_cpg=None):
    x_scale_windows = list()
    x_raw_meth_set = set(x_raw_meth)
    x_windows_raw_meth = list()
    y_windows_raw_meth = list()            
    x_scale_windows_raw_meth = list()            
    x_windows = list()
    y_windows_interp_meth = list()
    if args.interp_cpg_density:
        y_windows_interp_cpg = list()
    ### This subsets methylation/cpg density to the windows of interest
    if args.interp_cpg_density:        
        windows = subset_windows_of_interest(args,windows,x,x_raw_meth_set,x_raw_meth,y_raw_meth,interp_y_meth,interp_y_cpg=interp_y_cpg)  
    else:
        windows = subset_windows_of_interest(args,windows,x,x_raw_meth_set,x_raw_meth,y_raw_meth,interp_y_meth)     
    ### This section samples and scales x-coordinates for each window                                    
    for k,w in enumerate(windows):               
        w.sample_window()               
        w.set_counter(k)                
        w.scale_sample_x_coords()
        w.scale_raw_x_coords()
        ### Concatenate everything together for output
        x_windows += w.sample_x
        x_scale_windows += w.scale_sample_x
        y_windows_interp_meth +=w.sample_y_meth
        if args.interp_cpg_density:
            y_windows_interp_cpg += w.sample_y_cpg
        x_windows_raw_meth += w.x_raw_meth
        x_scale_windows_raw_meth += w.scale_x_raw_meth                
        y_windows_raw_meth += w.y_raw_meth
    
    if args.interp_cpg_density:
        return (x_scale_windows,x_windows_raw_meth,y_windows_raw_meth,x_scale_windows_raw_meth,x_windows,y_windows_interp_meth,y_windows_interp_cpg)
    else:
        return (x_scale_windows,x_windows_raw_meth,y_windows_raw_meth,x_scale_windows_raw_meth,x_windows,y_windows_interp_meth)

def subset_windows_of_interest(args,windows,x,x_raw_meth_set,x_raw_meth,y_raw_meth,interp_y_meth,interp_y_cpg=None):
    i = 0
    ### We'll left bound the windows.
    for j,pos in enumerate(x):
        if pos < windows[0].low:
            continue
        if pos >= windows[-1].high:
            break     
        if pos == windows[i].high:
            i+=1               
        if pos in x_raw_meth_set:
            windows[i].x_raw_meth.append(pos)                  
            windows[i].y_raw_meth.append(y_raw_meth[x_raw_meth.index(pos)])      
        windows[i].x.append(pos)
        windows[i].y_meth.append(interp_y_meth[j])
        if args.interp_cpg_density:
            windows[i].y_cpg.append(interp_y_cpg[j])
    
    if args.average:
        for j,w in enumerate(windows):
            if len(w.y_raw_meth)==0:
                #We're doing this in order from left to right, so there's has to be a value on one side                
                ### Average and actual methylation not the same thing in terms of averaging.                
                surrounding_y_raw_meth = list()              
                try:
                    prev_y_raw_meth = windows[j-1].y_raw_meth
                    ## If raw methylation doesn't exist, use the estimate from previous                    
                    if len(prev_y_raw_meth)==0:
                        prev_y_raw_meth = windows[j-1].y_meth
                    surrounding_y_raw_meth+=prev_y_raw_meth
                except IndexError:
                    continue
                try:
                    next_y_raw_meth = windows[j+1].y_raw_meth
                    ## There's not going to be an estimate if multiple windows are missing.
                    ## So just use raw. 
                    surrounding_y_raw_meth+=next_y_raw_meth
                except IndexError:
                    continue
                avg_window_meth = np.mean(surrounding_y_raw_meth)                
            else:  
                avg_window_meth = np.mean(w.y_raw_meth)
            w.y_meth = [avg_window_meth for z in w.y_meth]

#    print "x low: ",x[0],"x high: ",x[-1]
#    for w in windows:
#        print "Window Low: ",w.low,"Window High: ",w.high
#        print "Len of window x: ", len(w.x)
#        print "Average of raw meth: ", np.mean(w.y_raw_meth)
    return windows

def setup_windows(args,gene):
    windows = list()
    ### Need to fix the problem that TSS and WG need actual coordinates
    ### and WSG and SGF need bins. 
    if args.interp_type == "TSS":
        ### Single window for TSS
        if gene.strand == "+":
            windows.append(Window(gene.pos_low-gene.half_tss_window+args.tss_offset,gene.pos_low+gene.half_tss_window+args.tss_offset,args.tss_bins))
        else:
            windows.append(Window(gene.pos_high-gene.half_tss_window-args.tss_offset,gene.pos_high+gene.half_tss_window-args.tss_offset,args.tss_bins))
    elif args.interp_type == "WG":
        ### 3 windows for whole scaled gene:
        ### -low proximal window
        ### -whole gene window
        ### -high proximal window
        if args.proximal_window % args.wg_bp == 0:                
            proximal_bins = args.proximal_window/args.wg_bp
        else:
            proximal_bins = (args.proximal_window/args.wg_bp) + 1
        if (gene.pos_high-gene.pos_low) % args.wg_bp == 0:
            wg_bins = (gene.pos_high-gene.pos_low)/args.wg_bp
        else:
            wg_bins = ((gene.pos_high-gene.pos_low)/args.wg_bp) + 1 
        if args.proximal_window > 0:
            windows.append(Window(gene.pos_low-args.proximal_window,gene.pos_low,proximal_bins))
        windows.append(Window(gene.pos_low,gene.pos_high,wg_bins,args.wg_bp))
        if args.proximal_window > 0:                
            windows.append(Window(gene.pos_high,gene.pos_high+args.proximal_window,proximal_bins))             
    elif args.interp_type == "WSG":
        ### 3 windows for whole scaled gene:
        ### -low proximal window
        ### -whole gene window
        ### -high proximal window  
        if args.proximal_window > 0:
            windows.append(Window(gene.pos_low-args.proximal_window,gene.pos_low,args.proximal_bins))
        windows.append(Window(gene.pos_low,gene.pos_high,args.wsg_bins))
        if args.proximal_window > 0:                
            windows.append(Window(gene.pos_high,gene.pos_high+args.proximal_window,args.proximal_bins))                
    elif args.interp_type == "SGF": # interp_type == "SGF"
        ### -Low proximal window
        if args.proximal_window > 0:
            windows.append(Window(gene.pos_low-args.proximal_window,gene.pos_low-1,args.proximal_bins))
        ### -treat every exon and intron as a window                
        for i,e in enumerate(gene.exons):
            ### Append exon window                    
            windows.append(Window(e[0],e[1],args.exon_bins))
            ### Append intron window
            try:
                windows.append(Window(gene.introns[i][0],gene.introns[i][1],args.intron_bins))
            except IndexError:
                continue
        ### -High proximal window
        if args.proximal_window > 0:
            windows.append(Window(gene.pos_high+1,gene.pos_high+args.proximal_window,args.proximal_bins)) 
    else: ### interp_type == "ROI"
        ### -5 Low proximal windows
        for j in reversed(range(1,6)):
            windows.append(Window(gene.pos_low-(j*args.proximal_window),gene.pos_low-((j-1)*args.proximal_window)-1,args.proximal_bins))
#            windows.append(Window(gene.pos_low-args.proximal_window,gene.pos_low-1,args.proximal_bins))
        ### -treat every exon and intron as a window                
        for i,e in enumerate(gene.exons):
            ### Append exon window                    
            windows.append(Window(e[0],e[1],args.exon_bins))
            ### Append intron window
            try:
                windows.append(Window(gene.introns[i][0],gene.introns[i][1],args.intron_bins))
            except IndexError:
                continue
        ### -High proximal windows
        for j in range(5):
            windows.append(Window(gene.pos_high+(j*args.proximal_window)+1,gene.pos_high+((j+1)*args.proximal_window),args.proximal_bins)) 
#            windows.append(Window(gene.pos_high+1,gene.pos_high+args.proximal_window,args.proximal_bins)) 
    return windows            
        


def define_gene_boundaries(args,gene):
    if args.interp_type == "TSS":
        half_tss_window = args.tss_window/2
        lower = gene.pos_low-(half_tss_window+args.anchor_window)
        upper = gene.pos_high+(half_tss_window+args.anchor_window)
        interp_window_start = gene.pos_low-half_tss_window
        interp_window_stop = gene.pos_high+half_tss_window  
        gene.set_half_tss_window(half_tss_window)                          
    else:
        lower = gene.pos_low-(args.proximal_window+args.anchor_window)
        upper = gene.pos_high+(args.proximal_window+args.anchor_window)
        interp_window_start = gene.pos_low-args.proximal_window
        interp_window_stop = gene.pos_high+args.proximal_window
    return (lower,upper,interp_window_start,interp_window_stop)

def setup_interp_type(args):
    if args.interp_type=="TSS":
        ### TSS parameters:
        #   window size: 2*proximal window
        #   tss_bins: tss_bins
        interp_type_tag = "tss"
        interp_type_header = "#tssWin:%d,tssBins:%d,tssOffset:%d\n" % (args.tss_window,args.tss_bins,args.tss_offset)
    elif args.interp_type =="WG":
        ### WSG parameters:
        #   proximal_window
        #   wg_bp
        interp_type_tag = "wg"
        interp_type_header = "#proxWin:%d,wgBp:%d\n" % (args.proximal_window,args.wg_bp)    
    elif args.interp_type =="WSG":
        ### WSG parameters:
        #   proximal_window
        #   proximal_bins   
        #   wsg_bins
        interp_type_tag = "wsg"
        interp_type_header = "#proxWin:%d,proxBins:%d,wsgBins:%d\n" % (args.proximal_window,args.proximal_bins,args.wsg_bins)
    elif args.interp_type == "SGF": ## interp_type == "SGF"
        ### SGF parameters
        #   proximal_window
        #   proximal_bins
        #   exon_bins
        #   intron_bins        
        interp_type_tag = "sgf"
        interp_type_header = "#proxWin:%d,proxBins:%d,exonBins:%d,intronBins:%d\n" % (args.proximal_window,args.proximal_bins,args.exon_bins,args.intron_bins)
    elif args.interp_type == "DMR": ## interp_type == "DMR"
        ## DMR parameters
        #   DSS filename
        interp_type_tag = "dmr"
        interp_type_header = "#diffMeth distDmrTss lenDmr\n"
    else: ### interp_type == "ROI"
        ### Adjust args parameters to create ROI classifier        
        args.proximal_window = 400
        args.proximal_bins = 1
        args.exon_bins = 1
        args.intron_bins = 1
        args.average = True
        interp_type_tag = "roi"
        interp_type_header = "#Up5 Up4 Up3 Up2 Up1 FirstEx FirstIn IntnEx IntnIn LastIn LastEx Dw1 Dw2 Dw3 Dw4 Dw5\n"
    return (interp_type_tag,interp_type_header) 

def plot_wrapper(differential,output_dir,gene,args,x_windows,x_scale_windows,x_scale_windows_raw_meth,y_windows_interp_meth,x_windows_raw_meth,y_windows_raw_meth,y_windows_interp_cpg=None):
    if differential:
        if args.flankNorm:
            ylim = [-4,4]
        else:
            ylim = [-1.05,1]
    else:
        ylim = [-0.05,1]
    if args.interp_type=="TSS":
        ### There is no scaled  rsion of the plot since there are only actual coordinates
        plotname_actual = os.path.join(output_dir,"_".join([gene.sample,gene.gene_id,gene.common_name])+".actual.png")           
        title = " ".join([str(s) for s in [gene.sample,gene.gene_id,gene.common_name,gene.expr]])\
        +"\n"+gene.chrom+":"+str(gene.pos_low)+"-"+str(gene.pos_high)+" ("+gene.strand+")"
        if args.interp_cpg_density:
            if args.meth_cpg_density_score:
                tss_plot(x_windows,y_windows_interp_meth,x_windows_raw_meth,y_windows_raw_meth,title,plotname_actual,gene,pd=True,ps=True,d=y_windows_interp_cpg,ylim=ylim)
            else:
                tss_plot(x_windows,y_windows_interp_meth,x_windows_raw_meth,y_windows_raw_meth,title,plotname_actual,gene,pd=True,ps=False,d=y_windows_interp_cpg,ylim=ylim)
        else:
            tss_plot(x_windows,y_windows_interp_meth,x_windows_raw_meth,y_windows_raw_meth,title,plotname_actual,gene,pd=False,ps=False,ylim=ylim)                
    else:                
        ### For the whole scaled gene, and whole scaled gene features, plot 2 versions:
        ### actual version: (this is what the interpolation would look like on the genome browser)
        ### scaled version: (this is the scaled version with unit length for each feature)
        plotname_actual = os.path.join(output_dir,"_".join([gene.sample,gene.gene_id,gene.common_name])+".actual.png")  
        plotname_scale = os.path.join(output_dir,"_".join([gene.sample,gene.gene_id,gene.common_name])+".scale.png")                   
        title = " ".join([str(s) for s in [gene.sample,gene.gene_id,gene.common_name,gene.expr]])\
        +"\n"+gene.chrom+":"+str(gene.pos_low)+"-"+str(gene.pos_high)+" ("+gene.strand+")"
        if args.interp_cpg_density:
            if args.meth_cpg_density_score:
                ##Plot actual                        
                gene_plot(x_windows,y_windows_interp_meth,x_windows_raw_meth,y_windows_raw_meth,title,plotname_actual,gene,pd=True,ps=True,d=y_windows_interp_cpg,ylim=ylim)       
                if args.interp_type == "WSG" or args.interp_type == "SGF" or args.interp_type=="ROI":                        
                    ##Plot scaled
                    gene_plot(x_scale_windows,y_windows_interp_meth,x_scale_windows_raw_meth,y_windows_raw_meth,title,plotname_scale,gene,pd=True,ps=True,d=y_windows_interp_cpg,interp_type=args.interp_type,scaled=True,ylim=ylim)
                
            else:
                ##Plot actual                        
                gene_plot(x_windows,y_windows_interp_meth,x_windows_raw_meth,y_windows_raw_meth,title,plotname_actual,gene,pd=True,ps=False,d=y_windows_interp_cpg,ylim=ylim)       
                if args.interp_type == "WSG" or args.interp_type == "SGF" or args.interp_type=="ROI":                        
                    ##Plot scaled
                    gene_plot(x_scale_windows,y_windows_interp_meth,x_scale_windows_raw_meth,y_windows_raw_meth,title,plotname_scale,gene,pd=True,ps=False,d=y_windows_interp_cpg,interp_type=args.interp_type,scaled=True,ylim=ylim)
        else:
            ##Plot actual                        
            gene_plot(x_windows,y_windows_interp_meth,x_windows_raw_meth,y_windows_raw_meth,title,plotname_actual,gene,ylim=ylim)       
            if args.interp_type == "WSG" or args.interp_type == "SGF" or args.interp_type == "ROI":                            
                ##Plot scaled
                gene_plot(x_scale_windows,y_windows_interp_meth,x_scale_windows_raw_meth,y_windows_raw_meth,title,plotname_scale,gene,scaled=True,interp_type=args.interp_type,ylim=ylim)

def plot_simple(x,y):
    plt.figure(figsize=(8,5))
    plt.plot(x,y,'k.')
    plt.savefig("test.png")
    plt.close()
    return
    
def tss_plot_high_vis(x,y,rx,ry,title,plotname,gene,pd=False,ps=False,d=[],ylim=[-0.05,1]):
    ### Need to somehow delinate features based on the data
    ### We know that the boundaries are based on unit length, so draw at 0,1,2,...  
    ### for now these will be vertical lines, but I want to have some sort of horizontal lines at the bottom      
    plt.figure(figsize=(18,6))    
    plt.plot(x,y,'b-',lw=5,color="blue",label="Interpolated % Methylation")
    if pd:    
        plt.plot(x,d,'-',lw=5,color="green",label="Interpolated CpG Density")
    plt.plot(rx,ry,'k.',ms=10,alpha=1,label="Raw Methylation")    

    if gene.strand == "-":
        tss = gene.pos_high
#        tes = gene.pos_low
    else:
        tss = gene.pos_low
#        tes = gene.pos_high
#    for i,exon in enumerate(gene.exons):    
#        plt.plot([exon[0],exon[1]],[ylim[0]+0.02,ylim[0]+0.02],'-',color='purple',lw=5)
#    for i,intron in enumerate(gene.introns):
#        plt.plot([intron[0],intron[1]],[ylim[0]+0.02,ylim[0]+0.02],'y-',lw=5)
#    plt.plot([tss,tss],[ylim[0]+0.04,ylim[0]],'g-',lw=2,label="TSS")
#    plt.plot([tes,tes],[ylim[0]+0.04,ylim[0]],'r-',lw=2,alpha=0.5,label="TES") 
    plt.plot([x[0],x[-1]],[0,0],'-',lw=1,color='gray',alpha=0.5)
    plt.plot([tss,tss],[ylim[0]+0.05,ylim[1]],'-',lw=1,color='gray',alpha=0.5)
    if gene.strand == "-":
        plt.xlim(x[-1],x[0])
    else:
        plt.xlim(x[0],x[-1])
    plt.ylim(ylim)
    plt.ylabel("Differential Fractional Methylation")
    plt.axis('off')
    plt.tight_layout()
#    plt.xlabel("Genomic Position (bp)")
#    plt.legend(loc="best",fontsize='xx-small')
#    plt.title(title)
    plt.savefig(plotname)
    plt.close()
    return    


def tss_plot(x,y,rx,ry,title,plotname,gene,pd=False,ps=False,d=[],ylim=[-0.05,1]):
    ### Need to somehow delinate features based on the data
    ### We know that the boundaries are based on unit length, so draw at 0,1,2,...  
    ### for now these will be vertical lines, but I want to have some sort of horizontal lines at the bottom      
    plt.figure(figsize=(8,5))    
    plt.plot(x,y,'b-',label="Interpolated % Methylation")
    if pd:    
        plt.plot(x,d,'-',color="green",label="Interpolated CpG Density")
    if ps:
        s = normalize_methylation_cpg_density_wrapper(y,d)
        plt.plot(x,s,'-',color="firebrick",label="Interpolated Methylation * CpG Density")
    plt.plot(rx,ry,'k.',ms=3,alpha=0.5,label="Raw Methylation")    

    if gene.strand == "-":
        tss = gene.pos_high
        tes = gene.pos_low
    else:
        tss = gene.pos_low
        tes = gene.pos_high
    for i,exon in enumerate(gene.exons):    
        plt.plot([exon[0],exon[1]],[ylim[0]+0.02,ylim[0]+0.02],'-',color='purple',lw=2)
    for i,intron in enumerate(gene.introns):
        plt.plot([intron[0],intron[1]],[ylim[0]+0.02,ylim[0]+0.02],'y-',lw=2)
    plt.plot([tss,tss],[ylim[0]+0.04,ylim[0]],'g-',lw=2,label="TSS")
    plt.plot([tes,tes],[ylim[0]+0.04,ylim[0]],'r-',lw=2,alpha=0.5,label="TES") 
    plt.plot([x[0],x[-1]],[0,0],'-',lw=1,color='gray',alpha=0.5)
    plt.plot([tss,tss],[ylim[0]+0.05,ylim[1]],'-',lw=1,color='gray',alpha=0.5)
    if gene.strand == "-":
        plt.xlim(x[-1],x[0])
    else:
        plt.xlim(x[0],x[-1])
    plt.ylim(ylim)
    plt.ylabel("% Methylation / CpG Density (200bp)")

    plt.xlabel("Genomic Position (bp)")
    plt.legend(loc="best",fontsize='xx-small')
    plt.title(title)
    plt.savefig(plotname)
    plt.close()
    return      



def gene_plot(x,y,rx,ry,title,plotname,gene,pd=False,ps=False,d=[],scaled=False,interp_type="WSG",ylim=[-0.05,1]):
    ### Need to somehow delinate features based on the data
    ### We know that the boundaries are based on unit length, so draw at 0,1,2,...  
    ### for now these will be vertical lines, but I want to have some sort of horizontal lines at the bottom      
    plt.figure(figsize=(8,5))    
    plt.plot(x,y,'b-',label="Interpolated % Methylation")
    if pd:    
        plt.plot(x,d,'-',color="green",label="Interpolated CpG Density")
    if ps:
        s = normalize_methylation_cpg_density_wrapper(y,d)
        plt.plot(x,s,'-',color="firebrick",label="Interpolated Methylation * CpG Density")        
    plt.plot(rx,ry,'k.',ms=3,alpha=0.5,label="Raw Methylation")    
    if scaled:
        gene_len = gene.pos_high-gene.pos_low        
        if interp_type=="WSG":
            pos_low = int(math.floor(x[0]))
            pos_high = int(math.ceil(x[-1]))
            if gene.strand == "-":
                tss = pos_high-1
                tes = pos_low+1
            else:
                tss = pos_low+1
                tes = pos_high-1
            for i,exon in enumerate(gene.exons): 
                exon_low_perc = float(exon[0]-gene.pos_low)/gene_len
                exon_high_perc = float(exon[1]-gene.pos_low)/gene_len
                plt.plot([exon_low_perc+1,exon_high_perc+1],[ylim[0]+0.02,ylim[0]+0.02],'-',color='purple',lw=2,alpha=0.5)
            for i,intron in enumerate(gene.introns):
                intron_low_perc = float(intron[0]-gene.pos_low)/gene_len
                intron_high_perc = float(intron[1]-gene.pos_low)/gene_len
                plt.plot([intron_low_perc+1,intron_high_perc+1],[ylim[0]+0.02,ylim[0]+0.02],'y-',lw=2,alpha=0.5)
            plt.plot([tss,tss],[ylim[0]+0.04,ylim[0]],'g-',lw=2,alpha=0.5,label="TSS")
            plt.plot([tes,tes],[ylim[0]+0.04,ylim[0]],'r-',lw=2,alpha=0.5,label="TES")                
        else:       
            pos_low = int(math.floor(x[0]))
            pos_high = int(math.ceil(x[-1]))
            if gene.strand == "-":
                tss = pos_high-1
                tes = pos_low+1
            else:
                tss = pos_low+1
                tes = pos_high-1
            for i in range(pos_low+1,pos_high-1):
                if i%2==1:            
                    plt.plot([i,i+1],[ylim[0]+0.02,ylim[0]+0.02],'-',color='purple',lw=2,alpha=0.5)
                else:
                    plt.plot([i,i+1],[ylim[0]+0.02,ylim[0]+0.02],'y-',lw=2,alpha=0.5)
            plt.plot([tss,tss],[ylim[0]+0.04,ylim[0]],'g-',lw=2,alpha=0.5,label="TSS")
            plt.plot([tes,tes],[ylim[0]+0.04,ylim[0]],'r-',lw=2,alpha=0.5,label="TES") 
    else:
        if gene.strand == "-":
            tss = gene.pos_high
            tes = gene.pos_low
        else:
            tss = gene.pos_low
            tes = gene.pos_high

        for i,exon in enumerate(gene.exons):    
            plt.plot([exon[0],exon[1]],[ylim[0]+0.02,ylim[0]+0.02],'-',color='purple',lw=2,alpha=0.5)
        for i,intron in enumerate(gene.introns):
            plt.plot([intron[0],intron[1]],[ylim[0]+0.02,ylim[0]+0.02],'y-',lw=2,alpha=0.5)
        plt.plot([tss,tss],[ylim[0]+0.04,ylim[0]],'g-',lw=2,alpha=0.5,label="TSS")
        plt.plot([tes,tes],[ylim[0]+0.04,ylim[0]],'r-',lw=2,alpha=0.5,label="TES")
    plt.plot([x[0],x[-1]],[0,0],'-',lw=1,color='gray',alpha=0.5)
    plt.ylim(ylim)
    if gene.strand == "-":
        plt.xlim(x[-1],x[0])
    else:
        plt.xlim(x[0],x[-1])
    plt.ylabel("% Methylation / CpG Density (200bp)")
    if scaled:
        plt.xlabel("Unit length bins")
    else: 
        plt.xlabel("Genomic Position (bp)")
    plt.legend(loc="best",fontsize='xx-small')
    plt.title(title)
    plt.savefig(plotname)
    plt.close()
    return        

def downsample(x,r,b):
    x = np.array(x)
    if (x.size % r) !=0:
        pad_size = math.ceil(float(x.size)/r)*r-x.size
        x_pad = np.append(x,np.zeros(pad_size)*np.NaN)
        ds = scipy.nanmean(x_pad.reshape(-1,r),axis=1)
    else:
        ds = x.reshape(-1,r).mean(axis=1)
    #Average the remaining values if indivisible by bin size.
    if len(ds)>b:
        avg_lbs = np.mean(ds[b-1:])
        ds = ds[:b]
        ds[-1] = avg_lbs
    return list(ds)
    
def upsample(x,b):
    ### this function works but I wouldn't     
    
    ## Since downsampling is a simple average,
    ## upsampling will be a simple linear interpolation
    ## we also know that every xi in x is exactly 1 bp apart
    ## this gives us a length to work from.
    x_n = np.linspace(0,len(x)-1,num=b)
    x_adj = np.arange(0,len(x))
    f = scipy.interpolate.interp1d(x_adj,x)
    return list(f(x_n))

def flip_window(y,strand):
    ### Check strand for bin order direction
    if strand == "+":
        flip_y = y
    ### Adjust bin order for negative strand
    else:
        flip_y = y[::-1]
    return flip_y

def flank_zscore_normalize_values(x_raw_meth,y_raw_meth,x,interp_window_start,interp_window_stop,lower,upper):
    flank_y = list()    
#    gene_y = list()
    for i,x in enumerate(x_raw_meth):
        if lower<=x<=interp_window_start or interp_window_stop<=x<=upper:
            flank_y.append(y_raw_meth[i])
        else:
#            gene_y.append(y_raw_meth[i])
            continue
#    plot_simple_histogram(flank_y,tag+".flank_y.png")
#    plot_simple_histogram(gene_y,tag+".gene_y.png")
    flank_y = np.array(flank_y)
    flank_mean = flank_y.mean()
    flank_std = flank_y.std()
    zscore_norm_y = (y_raw_meth-flank_mean)/flank_std
#    zscore_norm_flank_y = (flank_y-flank_mean)/flank_std
#    zscore_norm_gene_y = (gene_y-flank_mean)/flank_std
#    plot_simple_histogram(zscore_norm_flank_y,tag+".zscore_flank_y.png")
#    plot_simple_histogram(zscore_norm_gene_y,tag+".zscore_gene_y.png")
    return zscore_norm_y

def plot_simple_histogram(x,outname):    
    plt.figure()    
    sns.set_style("white")    
    sns.distplot(x,kde=True)
    plt.title(outname+" # CpGs: %d"%(len(x)))
    plt.savefig(outname)
    plt.close()
    return

def interpolate_defined_window(raw_pos,raw_val,x,sigma,start,stop):        
    if len(raw_pos) >= 2:
        ### Linear Interpolation                
#        pre_smooth_interp_y = scipy.interpolate.interp1d(raw_pos,raw_val,bounds_error=False,fill_value=0)(x)
#        interp_y = pre_smooth_interp_y
        ### PCHIP Interpolation
        pre_smooth_interp_y = scipy.interpolate.pchip_interpolate(raw_pos,raw_val,x)
        interp_y = scipy.ndimage.filters.gaussian_filter1d(pre_smooth_interp_y,sigma)          
    #Single CpG in region of interest
    elif len(raw_pos) == 1:
        interp_y = [float(raw_val[0]) for val in x]
    #No CpGs in region of interest
    #No good way to deal with this and cannot put NA into learning algorithm
    #Imputation also doesn't make any sense, but a non-informative line does
    else:
        interp_y = [float(0) for val in x]    
    return interp_y

def methylation_sanity_filter(y,diff):
    l = list()
    for el in y:
        if diff:
            if el > 1:
                l.append(1.0)
            elif el < -1:
                l.append(-1.0)
            else:
                l.append(el)
        else:
            if el > 1:
                l.append(1.0)
            elif el < 0:
                l.append(0.0)
            else:
                l.append(el)
    return l

def get_position_value_from_range(min_bound,max_bound,tup):
    sample_meth_list = tup[1]
    ### previous start index    
    psi = tup[0]
    #Input: bounds to search, list of methylation values 
    #Output: List of methylation values for given bounds 
    search_bound = 1000
    #Initialize return structure
    pos = list()
    val = list()
    #Do single search through chromosome for gene
    #Speedup for searching through chromosome
    len_list = len(sample_meth_list)
    cc = psi
    nsi = -1
    while cc != len_list:      
        cpg = sample_meth_list[cc]        
        #Haven't gotten to gene area yet
        if cpg[0] < min_bound-search_bound:
            cc+=1            
            continue
        #Already past gene area
        elif cpg[0] > max_bound+search_bound:
            break
        #In region of gene
        else:
            ## This only trips the first time you come in here.
            if nsi < 0:
                nsi = cc
            if check_between(cpg[0],min_bound,max_bound):
                pos.append(cpg[0])
                val.append(cpg[1])
            else:
                cc+=1                
                continue
        cc+=1       
    return (pos,val,nsi)

def get_max_min_from_list(roi_list):
    ret_list = list()
    for roi in roi_list:
        if roi[0] == -1 or roi[1] == -1:
            continue
        else:
            ret_list.append(roi[0])
            ret_list.append(roi[1])
    return (max(ret_list),min(ret_list))

def define_region(mid,left,right):
    return (mid+left,mid+right)

def check_overlap_range(x,y):
    return range(max(x[0], y[0]), min(x[-1], y[-1])+1)
    
def check_overlap(x,y):
    if max(x[1],y[1]) - min(x[0],y[0]) <= (x[1] - x[0]) + (y[1] - y[0]):
        return True
    else:
        return False

def merge_intervals(intervals):   
    ### Softens the bound for overlap. 
    ### Not only can the intervals overlap, but they can be close by
    ### this amount.    
    overlap_bound = 1     
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0]) 
    merged = list()
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            ### Here is where to add the 
            if higher[0]-overlap_bound <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)
            else:
                merged.append(higher)
    return merged        
        
def get_most_common_element(mylist):
    return max(set(mylist), key=mylist.count)

def check_between(check_val,lower,upper):
    if lower <= check_val <= upper:
        return True
    else:
        return False 
        
def normalize_methylation_cpg_density_wrapper(mlist,clist):
    score = list()    
    for i,m in enumerate(mlist):
        score.append(normalize_methlylation_cpg_density(m,clist[i]))
    return score
      
def normalize_methlylation_cpg_density(m,c):
    scaling_factor = float(100)    
    return m*c*scaling_factor

def identify_closest_dmr(tss,chr_meth,strand):
    dist_list = np.array([min(math.fabs(tss-t[0]),math.fabs(tss-t[1])) for t in chr_meth])
    mdi = np.argmin(dist_list)
    ## If the tss is within the DMR, the distance will be 0
    if chr_meth[mdi][0]<=tss<=chr_meth[mdi][1]:
        chr_meth_dist = 0
    else:
        chr_meth_dist_low = tss-chr_meth[mdi][0]
        chr_meth_dist_high = tss-chr_meth[mdi][1]
        ## Have to switch the directionality so that positive refers to downstream,
        ## negative refers to upstream. Also, strand dependent. 
        if strand == "+":
            if math.fabs(chr_meth_dist_low) < math.fabs(chr_meth_dist_high):
                chr_meth_dist = -(chr_meth_dist_low)
            else:
                chr_meth_dist = -(chr_meth_dist_high)
        else:
            if math.fabs(chr_meth_dist_low) < math.fabs(chr_meth_dist_high):
                chr_meth_dist = chr_meth_dist_low
            else:
                chr_meth_dist = chr_meth_dist_high
    dat_tup = (chr_meth[mdi][2],chr_meth_dist,chr_meth[mdi][3])
    return dat_tup
        
### Input Functions

def load_methylation(meth_file):
    #Input: bedms file
    #Output: Hash Table
    #   Key: chr    Value: list of CpG positions
    fh = open(meth_file,'r')
    chr_hash = dict()
    sys.stderr.write("Loading methylation values for %s\n" % (meth_file))
    for line in fh:
        ll = line.strip().split()
        chr = ll[0]
        if not chr in valid_chr:
            continue
        cpg_pos = int(ll[2])-1
        val = float(ll[3])
        tup = (cpg_pos,val)
        try:
            chr_hash[chr][1].append(tup)
        except KeyError:
            sys.stderr.write("%s," % (chr))
            ### hash table will have an index of the beginning index of the last looked up 
            ### gene, so that there is an index to start from rather the beginning of the list
            chr_hash[chr] = [0,list()]
            chr_hash[chr][1].append(tup)
    sys.stderr.write("\n")
    return chr_hash
    
def load_methylation_dmr(meth_file):
    #Input: DSS output file
    #Output: Hash Table
    #   Key: chr    Value: (start,end,diff_meth)
    fh = open(meth_file,'r')
    chr_hash = dict()
    sys.stderr.write("Loading DMR methylation values for %s\n" % (meth_file))
    i = 0    
    for line in fh:
        ll = line.strip().split()
        if i == 0:
            ci = ll.index("chr")
            si = ll.index("start")
            ei = ll.index("end")
            li = ll.index("length")
            dmi = ll.index("diff.Methy")
            i+=1
            continue
        else:
            chr = ll[ci]
            try:
                start = int(ll[si])
            except ValueError:
                start = int(float(ll[si]))
            try:
                end = int(ll[ei])
            except ValueError:
                end = int(float(ll[ei]))
            length = int(ll[li])
            diff_meth = float(ll[dmi])
            if not chr in valid_chr:
                continue
            tup = (start,end,diff_meth,length)
            try:
                chr_hash[chr].append(tup)
            except KeyError:
                sys.stderr.write("%s," % (chr))
                chr_hash[chr] = list()
                chr_hash[chr].append(tup)
    sys.stderr.write("\n")
    return chr_hash    

def load_gene_expression_rpkm(sample_files,gene_annot):
    # Input: sample_files hash table, gene_annot hash table (key: gene_id, value: Gene)
    # Output: List of Genes
    #   Sorted by 1) sample 2) chr 3) pos_low
    label_list = list()
    total_num_genes = 0
    differential = False
    for sample,files in sample_files.items():
        sys.stderr.write("Processing Sample: %s\n" % (sample))
        ### Perform an initial search to exclude genes where the expression is defined twice. 
        ### This would effectively inject pure noise into the classification/regression
        ### because it's the same methylation pattern with different expression.
        ### this is different than the multiple samples with the same gene because 
        ### there are potentially different methylation patterns in different samples.
        uniq_sample_gene_ids = [x.strip().split()[0] for x in open(files[0],'r')]
        dup_sample_gene_ids = set([k for k,v in Counter(uniq_sample_gene_ids).items() if v > 1]) 
        i = 0
        for line in open(files[0],'r'):
            if line.startswith("#"):
                continue
            ll = line.strip().split()
            gene_id = ll[0]
            expr1 = float(ll[1])
            try:
                expr2 = float(ll[2])
                if expr1 > expr2:
                    expr = -(expr1/expr2)
                elif expr2 > expr1:
                    expr = expr2/expr1
                else: #expr1 == expr2
                    expr = float(0)
                differential = True
            except IndexError:
                expr = expr1
            
            if gene_id in dup_sample_gene_ids:
                sys.stderr.write("\tWarning: %s found >1 in sample expression file: %s\n" % (gene_id,files[0]))
                continue
            try:
                gene = copy.deepcopy(gene_annot[gene_id])
                gene.set_expr(expr)
                gene.set_sample(sample)
            except KeyError:
                sys.stderr.write("\tWarning: %s not found in reference database\n" % (gene_id))
                continue
            label_list.append(gene)
            i+=1
        sys.stderr.write("Obtained %d genes for sample %s\n\n" % (i,sample))
        total_num_genes += i
    sys.stderr.write("Obtained %d genes for %d samples\n" % (total_num_genes,len(sample_files.keys())))
    label_list.sort(key = lambda j: (j.sample,j.chrom,j.pos_low))
    return (label_list,differential)         


def load_sample_info(bedgraph,expr,tag):
    ### This is a hack, and can be extended later to process multiple samples per run.    
    sample_files = dict()
    if os.path.exists(bedgraph) and os.path.exists(expr):
        tup = (expr,bedgraph)
        sample_files[tag] = tup
    else:
        sys.stderr.write("Specified files do not exist\n")
        callHelp(parser)
    sys.stderr.write("\tLoaded sample\n")
    return sample_files

def load_conf_file(conf_file):
    # Input: 3-column configuration file:
    #   <sample> <expr5_path> <bedms_path>
    # Output: Hash Table
    #   Key: <sample> Value: Tuple of (expr5_file,bedms_file)
    sample_files = dict()
    sys.stderr.write("Loading configuration file from %s\n" % (conf_file))
    i = 0
    for line in open(conf_file):
        ll = line.strip().split()
        if len(ll)!=3:
            sys.stderr.write("Incorrect format for %s\n" % (conf_file))
            callHelp(parser)
        samp_name = ll[0]
        expr_name = ll[1]
        bedms_name = ll[2]
        i+=1
        if os.path.exists(expr_name) and os.path.exists(bedms_name):
            sample_files[samp_name] = (expr_name,bedms_name)
        else:
            sys.stderr.write("Cannot find files for sample %s\n" % (samp_name))
            callHelp(parser)
    sys.stderr.write("\tLoaded %d samples\n" % (i))
    return sample_files 

def load_gene_annot(annot_file,args):
    sys.stderr.write("Processing Annotation from %s\n" % (annot_file))
    ## Get only gene ids with single ID per gene
    ## This is because we have a single expression value that needs to map 
    ## to the ID
    ## Define filtering flags and counts
    if args.interp_type == "ROI":
        filter_flags = FilterFlags(exclude_lt_four_exons=True)
    else:
        filter_flags = FilterFlags()
    filter_counts = FilterCounts()    
    ### Process as REFSEQ or ENSG
    first_id = open(annot_file,'r').readline().strip().split()[1]    
    if first_id.startswith("NM_") or first_id.startswith("NR_"):
        filter_flags.annot_filetype="REFSEQ"
        uniq_id_index = 1
    elif first_id.startswith("ENS"):
        filter_flags.annot_filetype="ENSG"
        filter_flags.only_mrna = False
        uniq_id_index = 12
    else:
        sys.stderr.write("Unidentified filetype: %s Please specify either refGene or ensGene\n")
        callHelp(parser)
    ### Process for obtaining unique, gene-level annotations:
    num_tx_ids, uniq_gene_ids = get_unique_gene_ids(annot_file,uniq_id_index)
    filter_counts.num_uniq_gene_ids = len(uniq_gene_ids)
    filter_counts.num_transcripts = num_tx_ids    
    ### Using hash of genes because each transcript is going to be referencing a gene id
    genes = {x:Gene(x) for x in uniq_gene_ids}
    ### Going through file and associating all transcripts with an individual gene id. 
    fh = open(annot_file,'r')
    for line in fh:
        ll = line.strip().split()
        if filter_flags.annot_filetype == "REFSEQ":
            g,transcript = process_refseq_line_gene_pred_ext(ll)
            ### REFSEQ has few transcripts per gene, and the first transcript in the file should be used.
            if len(genes[g].transcripts) == 0:
                genes[g].append_transcript(transcript)
        else: ## ENSG
            g,transcript = process_ensembl_line_gene_pred_ext(ll)
            ### ENSG has many transcripts that should all be merged.         
            genes[g].append_transcript(transcript)   
    ### Applying transcript filtering 
    # Complete CDS start/end sites
    if filter_flags.only_complete:
        sys.stderr.write("Filtering transcripts without complete CDS start & stop\n") 
        filter_non_complete_transcripts(genes,filter_counts)
    ### Merging all transcripts (post transcription filtering)
    sys.stderr.write("Merging transcripts into genes based on transcript filter(s)\n")
    merge_gene_transcripts_wrapper(genes,filter_counts)    
    ### Applying gene filtering    
    # Simple Chromosome name filtering
    if filter_flags.only_simple_chrs:
        sys.stderr.write("Filtering genes without simple chromosomes\n") 
        filter_non_simple_chromosomes(genes,filter_counts)
    # Protein coding filtering if REFSEQ
    if filter_flags.annot_filetype == "REFSEQ" and filter_flags.only_mrna:
        sys.stderr.write("Filtering non-protein coding genes\n") 
        filter_non_protein_coding_genes(genes,filter_counts)
    # Unique transcription start only
    if filter_flags.only_uniq_tx_start:
        sys.stderr.write("Filtering genes with non-unique Tx start sites\n")
        genes = filter_duplicate_tx_starts(genes,filter_counts)
    # Report 1 common gene per multiple REFSEQ IDs only
    if filter_flags.collapse_on_common_name:
        sys.stderr.write("Collapsing genes with multiple REFSEQ IDs per common gene name\n")
        genes = collapse_genes_by_common_name(genes,filter_counts)
    # Exclude genes with less than 4 exons
    if filter_flags.exclude_lt_four_exons:
        sys.stderr.write("Excluding genes with < 4 exons\n")
        genes = exclude_lt_exon_num(genes,4)
    filter_counts.num_final_loaded_genes = len(genes.keys())
    filter_counts.write_output(filter_flags)
    return (genes,filter_flags)

def exclude_lt_exon_num(genes,exon_num):
    genes_out = dict()
    for gid,gene in genes.items():
        if len(gene.exons) >= exon_num:
            genes_out[gid] = gene
    return genes_out

def sort_genes(genes):
    label_list = list()
    for gid,gene in genes.items():
        label_list.append(gene)
    label_list.sort(key = lambda j: (chrom_to_num(j.chrom),j.pos_low))
    return label_list

def chrom_to_num(chrom_name):
    lc = chrom_name.lstrip("chr") 
    if lc == "X":
        rc = 23
    elif lc == "Y":
        rc = 24
    else:
        rc = int(lc)
    return rc

def find_smallest_refseq_id(rsl):
    nl = list()
    for r in rsl:
        nl.append(int(r.lstrip("NM_")))
    return rsl[nl.index(min(nl))]   

def collapse_genes_by_common_name(in_gene_annot,filter_counts):
    cnd = dict()
    for gid,g in in_gene_annot.items():
        try:
            cnd[g.common_name].append(gid)
        except KeyError:
            cnd[g.common_name] = list()
            cnd[g.common_name].append(gid)
    out_gene_annot = dict()
    for gcn,gid_list in cnd.items():
        if len(gid_list) > 1:
            small_gid = find_smallest_refseq_id(gid_list)
            out_gene_annot[small_gid] = in_gene_annot[small_gid]
        else:
            out_gene_annot[gid_list[0]] = in_gene_annot[gid_list[0]]
    return out_gene_annot

def filter_non_protein_coding_genes(genes,filter_counts):
    pop_list = list()
    for gid in genes.keys():
        if not gid.startswith("NM_"):
            pop_list.append(gid)
            filter_counts.num_genes_notPC+=1
        else:
            filter_counts.num_genes_mrna+=1
    for gid in pop_list:
        del genes[gid]
    return

def filter_non_simple_chromosomes(genes,filter_counts):
    pop_list = list()    
    for gid,g in genes.items():
        if not g.chrom in valid_chr:
            pop_list.append(gid)
            filter_counts.num_genes_other_chrs+=1
        else:
            filter_counts.num_genes_simple_chrs+=1
    for gid in pop_list:
        del genes[gid]
    return

def filter_non_complete_transcripts(genes,filter_counts):
    for gid,g in genes.items():
        pop_list = list()
        for i,t in enumerate(g.transcripts):
            if t.cds_tx_start_cmpl != "cmpl" or t.cds_tx_end_cmpl != "cmpl":
                pop_list.append(i)
        ### Have to remove in reverse order
        for idx in sorted(pop_list, reverse=True):
            del g.transcripts[idx]
        filter_counts.num_tx_complete+=len(g.transcripts)
        filter_counts.num_tx_incomplete+=len(pop_list)   
    return

def merge_gene_transcripts_wrapper(genes,filter_counts):
    no_transcript_list = list()    
    for gid,g in genes.items():
        if len(g.transcripts) == 0:
            no_transcript_list.append(gid)
            filter_counts.num_gene_ids_no_tx_post_tx_filters+=1
        elif len(g.transcripts) == 1:
            filter_counts.num_gene_ids_post_tx_filters+=1
            t = g.transcripts[0]            
            g.set_location_information(t.chrom,t.strand,t.pos_low,t.pos_high,merge_intervals(t.exons),t.common_name)
            g.set_introns()
        else:
            g.merge_gene_transcripts()
            filter_counts.num_gene_ids_post_tx_filters+=1
    for gid in no_transcript_list:
        del genes[gid]
    return
    
def get_unique_gene_ids(filename,idx):
    gene_set = set()
    fh = open(filename,'r')
    i = 0    
    for line in fh:
        gene_set.add(line.strip().split()[idx])
        i+=1
    return (i,gene_set)

def filter_duplicate_tx_starts(in_gene_annot,filter_counts):
    chr_start_tups = list()
    dummy_list = list()
    for gid,g in in_gene_annot.items():
        if g.strand == "-":
            start = g.pos_high            
        else:
            start = g.pos_low
        chr_start_tups.append((g.chrom,start))
        dummy_list.append((g.gene_id,g.chrom,g.strand,g.pos_low,g.pos_high))
    dup_chr_start = set([k for k,v in Counter(chr_start_tups).items() if v > 1])
#    for k,v in Counter(chr_start_tups).items():
#        print k,v
#    for el in dummy_list:
#        print el
    out_gene_annot = dict()
    for gid,g in in_gene_annot.items():
        if g.strand == "-":
            start = g.pos_high
        else:
            start = g.pos_low
        if not (g.chrom,start) in dup_chr_start:
            out_gene_annot[gid] = g
            filter_counts.num_genes_uniq_tx_start+=1
        else:
            filter_counts.num_genes_nonuniq_tx_start+=1
    return out_gene_annot

def process_ensembl_line_gene_pred_ext(ll):
    tx_id = ll[1]
    chrom = ll[2]
    cds_start_stat = ll[13]
    cds_end_stat = ll[14]
    strand = ll[3]
    pos_low = int(ll[4])
    pos_high = int(ll[5])
    exon_starts = [int(x) for x in ll[9].split(',') if x != '']
    exon_stops = [int(x) for x in ll[10].split(',') if x != '']
    ## Exons will currently be a list of tuples of strandless starts and stops    
    exons = calculate_exons(exon_starts,exon_stops)
    gene_id = ll[12]
    common_name = "NA"
    return (gene_id,Transcript(tx_id,chrom,strand,pos_low,pos_high,exons,common_name,cds_start_stat,cds_end_stat)) 

def process_refseq_line_gene_pred_ext(ll):
    tx_id = ll[1]
    chrom = ll[2]
    cds_start_stat = ll[13]
    cds_end_stat = ll[14]
    strand = ll[3]
    pos_low = int(ll[4])
    pos_high = int(ll[5])
    exon_starts = [int(x) for x in ll[9].split(',') if x != '']
    exon_stops = [int(x) for x in ll[10].split(',') if x != '']
    ## Exons will currently be a list of tuples of strandless starts and stops    
    exons = calculate_exons(exon_starts,exon_stops)
    common_name = ll[12]
    return (tx_id,Transcript(tx_id,chrom,strand,pos_low,pos_high,exons,common_name,cds_start_stat,cds_end_stat))

def calculate_exons(exon_starts,exon_stops):
    exons = list()
    for i,start in enumerate(exon_starts):
        e_tup = (start,exon_stops[i])
        exons.append(e_tup)        
    return exons

def calculate_introns(exons):
    introns = list()    
    for i,exon_tup in enumerate(exons):
        try:
            i_tup = (exon_tup[1]+1,exons[i+1][0]-1) 
            introns.append(i_tup)
        except IndexError:
            continue
    return introns

def calculate_exons_introns(exon_starts,exon_stops):
    exons = list()
    introns = list()    
    for i,start in enumerate(exon_starts):
        e_tup = (start,exon_stops[i])
        exons.append(e_tup)        
        try:        
            i_tup = (exon_stops[i]+1,exon_starts[i+1]-1)
            introns.append(i_tup)
        except IndexError:
            continue
    return (exons,introns)     

def load_chr_lengths(genome_file,valid_chr):
    # Input: 3-column file
    #   <chr> <chr_length> <2-bit_file>
    # Output: Hash Table
    #   Key: <chr>  Value: <chr_length>
    sys.stderr.write("Loading Chromosome Annotation from %s\n" % (genome_file))
    chr_len = dict()
    i = 0
    for line in open(genome_file,'r'):
        ll = line.strip().split()
        chr = ll[0]
        if not chr in valid_chr:
            continue
        len = int(ll[1])
        chr_len[chr] = len
        i+=1
    sys.stderr.write("\tLoaded %d Chromosomes\n" % (i))
    return chr_len

def load_cpg_perc_ref(refFile):
    #Hash Table of Key: Chromosome, Value: Numpy Array of Positions
    ret_hash = dict()
    sys.stderr.write("Loading Reference CpG Information for: %s \n" % (refFile))
    for line in open(refFile,'r'):
        ll = line.strip().split()
        chrom = ll[0]
        pos = int(ll[1])
        perc = float(ll[2])
        tup = (pos,perc)
        try:
            ret_hash[chrom][1].append(tup)
        except KeyError:
            sys.stderr.write("%s," % (chrom))
            ### Loading index of the last start of the last gene that cpg density was on 
            ### to speed up the lookup of the cpg density values
            ret_hash[chrom] = [0,[tup]]
    sys.stderr.write("\n")
    return ret_hash

### Output Functions

def print_label_gene(g,outfile):
    exon_starts = ",".join([str(x[0]) for x in g.exons])
    exon_stops = ",".join([str(x[1]) for x in g.exons])
    wl = "\t".join([str(x) for x in [g.sample,g.gene_id,g.common_name,g.expr,g.chrom,g.pos_low,g.pos_high,g.strand,len(g.exons),exon_starts,exon_stops]])+"\n"
    outfile.write(wl)
    return

def print_label_list_header(gene_list,tag):
    outname = tag+".label"
    sys.stderr.write("Writing %s\n" % (outname))
    outfile = open(outname,'w')
    #Print Header
    header_list = ["SAMPLE","GENE_ID","GENE_NAME","EXPR","CHR","POS_LOW","POS_HIGH","STRAND","NUM_EXONS","EXON_STARTS","EXON_STOPS"] 
    outfile.write("#"+"\t".join(header_list)+"\n")
    return outfile

### Parser Functions

def setup_parser_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,description='''
For multiple types of gene level interpolation and scaling of WGBS methylation data.''')
    
    ### Required positional arguments
    ### Setup subparser    
    subparsers = parser.add_subparsers(dest='interp_type')
    parser.add_argument('bedgraph',help="Sample bedgraph or DMR (DSS) output file")
    parser.add_argument('expr',help="2- (single) or 3- (differential) column expression file with: <gene_id> <expression1> (<expression2>)") 
    parser.add_argument('tag',help="Shared prefix for all output files created")       
    
    ### TSS Subparser    
    sp_tss = subparsers.add_parser('TSS',help="TSS-centered (TSS)")
    sp_tss.add_argument('-c',"--interpCpGDensity",dest="interp_cpg_density",default=False,help="Interpolate CpG density as an orthogonal source of information")   
    sp_tss.add_argument('-g',"--genePredExt",dest="annot_file",default="/mnt/work2/genomeData/hg19/refGene.hg19.21_12_15.txt",help="REFSEQ or ENSG gene predictions (extended) format")
    sp_tss.add_argument('-d',"--cpgDensity",dest="cpg_density_file",default="/mnt/work2/genomeData/hg19/cpg_density_cpg_sites_all_hg19.window_200bp.txt",help="Interpolate CpG density as an additional curve for the given genes. Location of all CpGs in genome with density of CpGs within window. 3-Column format: <chr> <start> <density>")    
    sp_tss.add_argument('-p',"--plot",dest="plot",default=False,type=bool,help="Plot all interpolated curves in directory <tag>_interp_curves")
    sp_tss.add_argument('-a',"--anchorWindow",dest="anchor_window",default=100000,type=int,help="Window size for both sides of the interpolation for anchoring and flank normalization.")    
    sp_tss.add_argument('-f',"--flankNorm",dest="flankNorm",default=False,type=bool,help="Use z-score flank (anchor windows up/downstream) normalized values instead of raw methylation and CpG density.")
    sp_tss.add_argument('-v',"--average",dest="average",default=False,type=bool,help="Use average value of raw methylation within windows instead of interpolated values")
    sp_tss.add_argument('-y',"--minMeth",dest="min_meth",default=40,type=int,help="Minimum number of measured CpGs in whole gene region to interpolate")    
    sp_tss.add_argument('-q',"--minMethChange",dest="min_meth_change",default=0.2,type=float,help="Minimum methylation change in interpolant. (differential only)")   
    sp_tss.add_argument('-r',"--methCpGDensityScore",dest="meth_cpg_density_score",default=False,help="Report a combined methylation/CpG density score")
    sp_tss.add_argument('-t',"--tssWindow",dest="tss_window",default=10000,type=int,help="Window size for TSS-centered interpolation.") 
    sp_tss.add_argument('-o',"--tssOffset",dest="tss_offset",default=0,type=int,help="Offset for TSS-centered interpolation. -ve = upstream, +ve = downstream")
    sp_tss.add_argument('-n',"--tssBins",dest="tss_bins",default=500,type=int,help="TSS-centered default based on 20bp windows for 10kb window.")  
    sp_tss.add_argument('-s',"--sigma",dest="sigma",default=50,type=int,help="Sigma (in bp) for Gaussian smoothing filter")
    ### WG Subparser
    sp_wg = subparsers.add_parser('WG',help="Whole Gene Unscaled (WG)")
    sp_wg.add_argument('-c',"--interpCpGDensity",dest="interp_cpg_density",default=False,help="Interpolate CpG density as an orthogonal source of information")   
    sp_wg.add_argument('-g',"--genePredExt",dest="annot_file",default="/mnt/work2/genomeData/hg19/refGene.hg19.21_12_15.txt",help="REFSEQ or ENSG gene predictions (extended) format")
    sp_wg.add_argument('-d',"--cpgDensity",dest="cpg_density_file",default="/mnt/work2/genomeData/hg19/cpg_density_cpg_sites_all_hg19.window_200bp.txt",help="Interpolate CpG density as an additional curve for the given genes. Location of all CpGs in genome with density of CpGs within window. 3-Column format: <chr> <start> <density>")    
    sp_wg.add_argument('-p',"--plot",dest="plot",default=False,type=bool,help="Plot all interpolated curves in directory <tag>_interp_curves")
    sp_wg.add_argument('-a',"--anchorWindow",dest="anchor_window",default=100000,type=int,help="Window size for both sides of the interpolation for anchoring and flank normalization.")    
    sp_wg.add_argument('-f',"--flankNorm",dest="flankNorm",default=False,type=bool,help="Use z-score flank (anchor windows up/downstream) normalized values instead of raw methylation and CpG density.")
    sp_wg.add_argument('-v',"--average",dest="average",default=False,type=bool,help="Use average value of raw methylation within windows instead of interpolated values")
    sp_wg.add_argument('-y',"--minMeth",dest="min_meth",default=40,type=int,help="Minimum number of measured CpGs in whole gene region to interpolate")    
    sp_wg.add_argument('-q',"--minMethChange",dest="min_meth_change",default=0.2,type=float,help="Minimum methylation change in interpolant. (differential only)")      
    sp_wg.add_argument('-r',"--methCpGDensityScore",dest="meth_cpg_density_score",default=False,help="Report a combined methylation/CpG density score")
    sp_wg.add_argument('-x',"--proximalWindow",dest="proximal_window",default=5000,type=int,help="Window size for windows adjacent to the TSS and TES")
    sp_wg.add_argument('-b',"--wgbp",dest="wg_bp",default=40,type=int,help="Number of base pairs to downsample the interpolated curve") 
    sp_wg.add_argument('-s',"--sigma",dest="sigma",default=50,type=int,help="Sigma (in bp) for Gaussian smoothing filter")    
    ### WSG Subparser
    sp_wsg = subparsers.add_parser('WSG',help="Whole Gene Scaled (WSG)")
    sp_wsg.add_argument('-c',"--interpCpGDensity",dest="interp_cpg_density",default=False,help="Interpolate CpG density as an orthogonal source of information")   
    sp_wsg.add_argument('-g',"--genePredExt",dest="annot_file",default="/mnt/work2/genomeData/hg19/refGene.hg19.21_12_15.txt",help="REFSEQ or ENSG gene predictions (extended) format")
    sp_wsg.add_argument('-d',"--cpgDensity",dest="cpg_density_file",default="/mnt/work2/genomeData/hg19/cpg_density_cpg_sites_all_hg19.window_200bp.txt",help="Interpolate CpG density as an additional curve for the given genes. Location of all CpGs in genome with density of CpGs within window. 3-Column format: <chr> <start> <density>")    
    sp_wsg.add_argument('-p',"--plot",dest="plot",default=False,type=bool,help="Plot all interpolated curves in directory <tag>_interp_curves")
    sp_wsg.add_argument('-a',"--anchorWindow",dest="anchor_window",default=100000,type=int,help="Window size for both sides of the interpolation for anchoring and flank normalization.")    
    sp_wsg.add_argument('-f',"--flankNorm",dest="flankNorm",default=False,type=bool,help="Use z-score flank (anchor windows up/downstream) normalized values instead of raw methylation and CpG density.")
    sp_wsg.add_argument('-v',"--average",dest="average",default=False,type=bool,help="Use average value of raw methylation within windows instead of interpolated values")
    sp_wsg.add_argument('-y',"--minMeth",dest="min_meth",default=40,type=int,help="Minimum number of measured CpGs in whole gene region to interpolate")    
    sp_wsg.add_argument('-q',"--minMethChange",dest="min_meth_change",default=0.2,type=float,help="Minimum methylation change in interpolant. (differential only)")   
    sp_wsg.add_argument('-r',"--methCpGDensityScore",dest="meth_cpg_density_score",default=False,help="Report a combined methylation/CpG density score")
    sp_wsg.add_argument('-x',"--proximalWindow",dest="proximal_window",default=5000,type=int,help="Window size for windows adjacent to the TSS and TES")
    sp_wsg.add_argument('-m',"--proximalBins",dest="proximal_bins",default=125,type=int,help="Number of bins to divide the proximal window into.") 
    sp_wsg.add_argument('-w',"--wsgBins",dest="wsg_bins",default=500,type=int,help="Number of bins to divide the whole gene window into. Whole gene default (500) based on median size of protein coding gene in ENSG.") 
    sp_wsg.add_argument('-s',"--sigma",dest="sigma",default=50,type=int,help="Sigma (in bp) for Gaussian smoothing filter")        
    ### SGF Subparser
    sp_sgf = subparsers.add_parser('SGF',help="Scaled Gene Features (SGF)")
    sp_sgf.add_argument('-c',"--interpCpGDensity",dest="interp_cpg_density",default=False,help="Interpolate CpG density as an orthogonal source of information")   
    sp_sgf.add_argument('-g',"--genePredExt",dest="annot_file",default="/mnt/work2/genomeData/hg19/refGene.hg19.21_12_15.txt",help="REFSEQ or ENSG gene predictions (extended) format")
    sp_sgf.add_argument('-d',"--cpgDensity",dest="cpg_density_file",default="/mnt/work2/genomeData/hg19/cpg_density_cpg_sites_all_hg19.window_200bp.txt",help="Interpolate CpG density as an additional curve for the given genes. Location of all CpGs in genome with density of CpGs within window. 3-Column format: <chr> <start> <density>")    
    sp_sgf.add_argument('-p',"--plot",dest="plot",default=False,type=bool,help="Plot all interpolated curves in directory <tag>_interp_curves")
    sp_sgf.add_argument('-a',"--anchorWindow",dest="anchor_window",default=100000,type=int,help="Window size for both sides of the interpolation for anchoring and flank normalization.")    
    sp_sgf.add_argument('-f',"--flankNorm",dest="flankNorm",default=False,type=bool,help="Use z-score flank (anchor windows up/downstream) normalized values instead of raw methylation and CpG density.")
    sp_sgf.add_argument('-v',"--average",dest="average",default=False,type=bool,help="Use average value of raw methylation within windows instead of interpolated values")
    sp_sgf.add_argument('-y',"--minMeth",dest="min_meth",default=40,type=int,help="Minimum number of measured CpGs in whole gene region to interpolate")    
    sp_sgf.add_argument('-q',"--minMethChange",dest="min_meth_change",default=0.2,type=float,help="Minimum methylation change in interpolant. (differential only)")  
    sp_sgf.add_argument('-r',"--methCpGDensityScore",dest="meth_cpg_density_score",default=False,help="Report a combined methylation/CpG density score")
    sp_sgf.add_argument('-x',"--proximalWindow",dest="proximal_window",default=5000,type=int,help="Window size for windows adjacent to the TSS and TES")
    sp_sgf.add_argument('-m',"--proximalBins",dest="proximal_bins",default=125,type=int,help="Number of bins to divide the proximal window into") 
    sp_sgf.add_argument('-e',"--exonBins",dest="exon_bins",default=10,type=int,help="Number of bins to divide the scaled exons into")
    sp_sgf.add_argument('-i',"--intronBins",dest="intron_bins",default=30,type=int,help="Number of bins to divide the scaled introns into")
    sp_sgf.add_argument('-s',"--sigma",dest="sigma",default=50,type=int,help="Sigma (in bp) for Gaussian smoothing filter")    
    ### ROI Subparser
    sp_roi = subparsers.add_parser('ROI',help="Region of Interest (ROI) (Lou et al 2014) Implementation")    
    sp_roi.add_argument('-c',"--interpCpGDensity",dest="interp_cpg_density",default=False,help="Interpolate CpG density as an orthogonal source of information")   
    sp_roi.add_argument('-g',"--genePredExt",dest="annot_file",default="/mnt/work2/genomeData/hg19/refGene.hg19.21_12_15.txt",help="REFSEQ or ENSG gene predictions (extended) format")
    sp_roi.add_argument('-d',"--cpgDensity",dest="cpg_density_file",default="/mnt/work2/genomeData/hg19/cpg_density_cpg_sites_all_hg19.window_200bp.txt",help="Interpolate CpG density as an additional curve for the given genes. Location of all CpGs in genome with density of CpGs within window. 3-Column format: <chr> <start> <density>")    
    sp_roi.add_argument('-p',"--plot",dest="plot",default=False,type=bool,help="Plot all interpolated curves in directory <tag>_interp_curves")
    sp_roi.add_argument('-a',"--anchorWindow",dest="anchor_window",default=100000,type=int,help="Window size for both sides of the interpolation for anchoring and flank normalization.")    
    sp_roi.add_argument('-f',"--flankNorm",dest="flankNorm",default=False,type=bool,help="Use z-score flank (anchor windows up/downstream) normalized values instead of raw methylation and CpG density.")
    sp_roi.add_argument('-v',"--average",dest="average",default=False,type=bool,help="Use average value of raw methylation within windows instead of interpolated values")
    sp_roi.add_argument('-y',"--minMeth",dest="min_meth",default=40,type=int,help="Minimum number of measured CpGs in whole gene region to interpolate")    
    sp_roi.add_argument('-q',"--minMethChange",dest="min_meth_change",default=0.2,type=float,help="Minimum methylation change in interpolant. (differential only)")   
    sp_roi.add_argument('-r',"--methCpGDensityScore",dest="meth_cpg_density_score",default=False,help="Report a combined methylation/CpG density score")
    sp_roi.add_argument('-s',"--sigma",dest="sigma",default=50,type=int,help=argparse.SUPPRESS)  
    ### DMR Subparser
    sp_dmr = subparsers.add_parser('DMR',help="Differentially Methylated Region (DMR)")
    sp_dmr.add_argument('-c',"--interpCpGDensity",dest="interp_cpg_density",default=False,help="Interpolate CpG density as an orthogonal source of information")   
    sp_dmr.add_argument('-g',"--genePredExt",dest="annot_file",default="/mnt/work2/genomeData/hg19/refGene.hg19.21_12_15.txt",help="REFSEQ or ENSG gene predictions (extended) format")
    sp_dmr.add_argument('-d',"--cpgDensity",dest="cpg_density_file",default="/mnt/work2/genomeData/hg19/cpg_density_cpg_sites_all_hg19.window_200bp.txt",help="Interpolate CpG density as an additional curve for the given genes. Location of all CpGs in genome with density of CpGs within window. 3-Column format: <chr> <start> <density>")    
    sp_dmr.add_argument('-p',"--plot",dest="plot",default=False,type=bool,help="Plot all interpolated curves in directory <tag>_interp_curves")
    sp_dmr.add_argument('-a',"--anchorWindow",dest="anchor_window",default=100000,type=int,help="Window size for both sides of the interpolation for anchoring and flank normalization.")    
    sp_dmr.add_argument('-f',"--flankNorm",dest="flankNorm",default=False,type=bool,help="Use z-score flank (anchor windows up/downstream) normalized values instead of raw methylation and CpG density.")
    sp_dmr.add_argument('-v',"--average",dest="average",default=False,type=bool,help="Use average value of raw methylation within windows instead of interpolated values")
    sp_dmr.add_argument('-y',"--minMeth",dest="min_meth",default=40,type=int,help="Minimum number of measured CpGs in whole gene region to interpolate")    
    sp_dmr.add_argument('-q',"--minMethChange",dest="min_meth_change",default=0.2,type=float,help="Minimum methylation change in interpolant. (differential only)")
    sp_dmr.add_argument('-s',"--sigma",dest="sigma",default=50,type=int,help=argparse.SUPPRESS)
    return parser,subparsers


def callHelp(parser):
    parser.print_help()
    sys.exit()

if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')  
    import sys, os, argparse, copy, math, time
    from collections import Counter
    import scipy    
    import scipy.stats
    import scipy.interpolate
    import scipy.ndimage.filters
    import numpy as np  
    from matplotlib import pyplot as plt
    main()
