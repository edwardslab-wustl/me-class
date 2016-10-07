# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 16:42:29 2016

@author: cschlosberg
"""
   
class Sample:
    def __init__(self,sample_name,Y,P,M):
        self.sample_name = sample_name
        self.Y = np.asarray(Y)
        self.P = np.asarray(P)  
        self.M = np.asarray(M)        
        
    def set_roc(self,tpr,fpr,roc_auc,roc_curve):
        self.tpr = tpr
        self.fpr = fpr
        self.roc_auc = roc_auc
        self.roc_curve = roc_curve
        
    def set_pr(self,recall,precision,pr_auc,pr_curve):
        self.recall = recall
        self.precision = precision
        self.pr_auc = pr_auc
        self.pr_curve = pr_curve    
        
    def set_print_vals(self,print_rej_rate,print_acc):
        self.print_rej_rate = print_rej_rate
        self.print_acc = print_acc        
        
    def set_acc_thresh(self,acc,thresh):
        self.acc = acc
        self.thresh = thresh
        
        
def main():
    ### Set up parser arguments
    global parser
    global global_lw
    global axisSpacing
    parser = setup_parser_arguments()
    ### Grab argument variables
    args = parser.parse_args()   
    global_lw = args.lineWidth
    axisSpacing = args.axisSpacing
    method_results = list()
    for method in args.methods:
        method_results.append(process_samples(method))
    method_names = list()
    for method in args.methods:
        method_names.append(os.path.splitext(os.path.basename(method))[0])
    meth_tag = "_".join(method_names)
    tag = args.tag+"."+meth_tag
    method_results = calculate_performance_stats(method_results,method_names,tag,args.pred_prob)       
    plot_ROC(method_results,method_names,tag,args.legendLoc)

         
def plot_ROC(method_results,method_names,tag,legendLoc):
    """ Receiver Operator Characteristic Curve """
    rocs_df = list() 
#    print("ROC stats")
    for j,samples in enumerate(method_results):
        avg_auc_list = list()
        for sample in samples:
#            print(method_names[j],sample.roc_auc)
            avg_auc_list.append(sample.roc_auc)
        avg_auc = np.mean(avg_auc_list)
        std_auc = np.std(avg_auc_list)
        for sample in samples:
            for i,y in enumerate(sample.roc_curve):
                rocs_df.append([i,"%s (%.2f+/-%.2f)"%(method_names[j],avg_auc,std_auc),sample.sample_name,np.nan_to_num(y)])
    rocs_df = pd.DataFrame(rocs_df)
#    rocs_df.to_csv(tag+".roc.txt")

    fig = plt.figure()    
    sns.set(style='ticks',font_scale=4)
    sns.despine()
    plt.figure(figsize=(10,10))
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    sns.set_palette("husl", len(method_names))
    ax = sns.tsplot(rocs_df,time=0,condition=1,unit=2,value=3,err_style=None,linewidth=global_lw)
    sns.despine(fig)
    n = len(ax.xaxis.get_ticklabels())
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
#    plt.plot([x0,x1],[y0,y1],'k--')
    ax.set_xticklabels(np.linspace(0,1,n)) 
                
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    if legendLoc == "right":
        lgd=plt.legend(bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0., fontsize='xx-small', handlelength=1, handletextpad=0.3)
        #plt.savefig(tag+".roc.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        lgd=plt.legend(loc="best",fontsize="xx-small", handlelength=1, handletextpad=0.3)
        #plt.legend(loc="best",fontsize="xx-small")
        #plt.tight_layout()
        #plt.savefig(tag+".roc.png")
    plt.savefig(tag+".roc.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
    #plt.legend(loc="best",fontsize='xx-small')
    #plt.tight_layout()    
    #plt.savefig(tag+".roc.png")
    #JOHN
    #lgd=plt.legend(bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0., fontsize='xx-small', handlelength=1, handletextpad=0.3)
    #plt.savefig(tag+".roc.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close('all')

def calculate_performance_stats(method_results,method_names,tag,print_pred_prob):       
    """ Loop through .pred files and calculate results """    
    sfh = open(tag+".stats.txt",'w')
    spacing_curves=101    
    wl = "#"+"\t".join(['Sample','Method','Prob_Pred','Tot','TP','FP','TN','FN','ACC','PC','PPV','NPV'])+"\n"
    sfh.write(wl)
    for i,samples in enumerate(method_results):
        sys.stderr.write("Processing method: %s\n" % method_names[i])
        for sample in samples:  
            sys.stderr.write("\tProcessing sample: %s\n" % (sample.sample_name))
            fpr, tpr, thresholds = sklearn.metrics.roc_curve(sample.M,sample.P[:,1],pos_label=1)
            roc_auc = sklearn.metrics.auc(fpr,tpr)
            precision, recall, thresholds = sklearn.metrics.precision_recall_curve(sample.M,sample.P[:,1],pos_label=1) 
            pr_auc = sklearn.metrics.auc(recall,precision)
            print_acc, print_rejrate, acc, fixed_thresholds = calculate_accuracy_reject_rate(sample.M,sample.P,method_names[i],sample.sample_name,sfh,print_pred_prob)
            f = scipy.interpolate.interp1d(fpr,tpr)
            roc_curve = f(np.linspace(0,1,num=spacing_curves))
            f = scipy.interpolate.interp1d(recall,precision)
            pr_curve = f(np.linspace(0,1,num=spacing_curves))
            sample.set_roc(fpr,tpr,roc_auc,roc_curve)
            sample.set_pr(recall,precision,pr_auc,pr_curve)
            sample.set_print_vals(print_rejrate,print_acc)
            sample.set_acc_thresh(acc,fixed_thresholds)
    return method_results
            
def calculate_accuracy_reject_rate(y_true,y_score,method_name,sample_name,fh,print_pred_prob):    
    #['Sample','Method','Prob_Pred','Tot','TP','FP','TN','FN','ACC','PC','PPV','NPV']
    num_thresholds = 101
    pos_thresh = np.linspace(0.5,1,num_thresholds)[:(num_thresholds-1)]
    acc = list()
    fixed_thresholds = list()    
    for t in pos_thresh:
        t_true = list()
        t_score = list()
        for i,p in enumerate(y_score):
            max_p = max(p)
            if max_p >= t:
                t_true.append(y_true[i])
                t_score.append(y_score[i])
        t_pred = convert_score_pred(t_score)  
        t_rej_rate = len(t_true)/float(len(y_true))
        t_num_genes = int(t_rej_rate*len(y_true))
        t_acc = np.nan_to_num(sklearn.metrics.accuracy_score(t_true,t_pred))
        t_num_genes = len(t_true)
        if t == print_pred_prob:
            pred_prob_acc = t_acc
            pred_prob_rej_rate = t_rej_rate      
        acc.append(t_acc)
        fixed_thresholds.append(t)      
        cm = sklearn.metrics.confusion_matrix(t_true,t_pred)
        if cm.shape==(2,2):
            t_num_tp = cm[0][0]
            t_num_tn = cm[1][1]
            t_num_fp = cm[1][0]
            t_num_fn = cm[0][1]
            t_ppv = np.nan_to_num(t_num_tp/float(t_num_tp+t_num_fp))
            t_npv = np.nan_to_num(t_num_tn/float(t_num_tn+t_num_fn))
        elif cm.shape==(1,1):
            if all(t_true)==1:
                t_num_tp = t_num_genes
                t_num_tn = 0
                t_num_fp = 0
                t_num_fn = 0
                t_ppv = float(1)
                t_npv = float(0)
            else:
                t_num_tp = 0
                t_num_tn = t_num_genes
                t_num_fp = 0
                t_num_fn = 0
                t_ppv = float(0)
                t_npv = float(1)
        else:
            t_num_tp = 0
            t_num_tn = 0
            t_num_fp = 0
            t_num_fn = 0
            t_ppv = float(0)
            t_npv = float(0)
        wl = "\t".join([str(x) for x in [sample_name,method_name,t,t_num_genes,t_num_tp,t_num_fp,t_num_tn,t_num_fn,t_acc,t_rej_rate,t_ppv,t_npv]])+"\n"
        fh.write(wl)
    return pred_prob_acc,pred_prob_rej_rate,acc,fixed_thresholds


def convert_score_pred(score):
    ret_list = list()    
    for s in score:
        if np.argmax(s) == 0:
            ret_list.append(-1)
        else:
            ret_list.append(1)
    return ret_list

def calculate_accuracy(A):
    corr = A.count(1)
    return (float(corr)/len(A))*100
    
    
def process_samples(method_list):
    samples = list()
    for x in open(method_list,'r'):
        pred_file = x.strip()     
        if pred_file.endswith(".pred"):
            sample_name = pred_file.split(".")[0]   
            sys.stderr.write("\tPreparing evaluation for %s\n" %(sample_name))      
            Y,L,E,M,P = load_pred(pred_file)
            sample = Sample(sample_name,Y,P,M)
            samples.append(sample)
    return samples
    
def load_pred(filename):
    fh = open(filename,'r')
    Y = list()
    L = list()
    E = list()
    M = list()
    P = list()
    for line in fh:
        if line.startswith("#"):
            ll = line.lstrip("#").rstrip().split()
            pd = ll.index("PROB_DW")
            pu = ll.index("PROB_UP")
            ei = ll.index("EXPR")
            phi = ll.index("POS_HIGH")
            pli = ll.index("POS_LOW")            
            nei = ll.index("NUM_EXONS")
            continue
        ll = line.strip().split()
        expr = float(ll[ei])
        length = int(ll[phi])-int(ll[pli])
        num_exons = int(ll[nei])
        c = get_expression_binary_class(expr)
        pred_down = float(ll[pd])
        pred_up = float(ll[pu])
        Y.append(expr)
        M.append(c)
        L.append(length)
        E.append(num_exons)
        P.append([pred_down,pred_up])
    return (Y,L,E,M,P)  

def get_expression_binary_class(expr):
    if expr > 0:
        c = 1
    else:
        c = -1
    return c

def get_function_name(func):
    return str(func).split('(')[0]

def setup_parser_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,description='''
For plotting ROC for classification of gene expression based on interpolated methylation values over multiple samples and methods''')
    ### Required positional arguments
    parser.add_argument('tag',help="Output tag")
    parser.add_argument('methods',nargs='+',help="1 or more lists of <sample>.<method>.pred files")
    ### Optional arguments
    parser.add_argument('-p',"--pred_prob",dest="pred_prob",default=0.9,help="Prediction probability to use for printing box plots of accuracy and rejection rate")    
    parser.add_argument('--lineWidth',help="linewidth for plots, default=4",type=int, default=4)
    parser.add_argument('--legendLoc',help="location of legend for plots, default=best",choices=["best","right"],default="best")
    parser.add_argument('--axisSpacing',help="axis offset so you can see 0 pts, 0.02 works well. default=0.0",type=float, default=0)

    return parser

if __name__ == "__main__":
    import sys, os, argparse, linecache, operator, math
    import numpy as np

    import pandas as pd
    import scipy.stats
    import scipy.interpolate
    import sklearn.metrics
        
    from matplotlib import pyplot as plt
    import seaborn as sns
    import matplotlib as mpl
    main()    
    