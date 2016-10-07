# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 16:42:29 2016

@author: cschlosberg
"""

class Sample:
    def __init__(self,sample_name,X,Y,L,E,C,I,D,header):
        self.sample_name = sample_name
        self.individual_sample_names = set(sample_name.split("_"))
        ## Methylation Profiles        
        self.X = X 
        ## Expression Deciles           
        self.Y = Y
        ## log10 gene length
        self.L = L
        ## Number of exons
        self.E = E
        ## Mean CpG density over TSS region
        self.C = C    
        ## Gene ID
        self.I = I 
        ## label file info
        self.D = D
        ## label file header
        self.header = header
        
    def combine_features(self):
        F = list()        
        for i,x in enumerate(self.X):
            f = x
            F.append(f)
        self.F = F
        
    def reduce_memory(self):
        del self.X
        del self.L        
        del self.C
        del self.E
        
        
def main():
    ### Set up parser arguments
    global parser
    parser = setup_parser_arguments()
    ### Grab argument variables
    args = parser.parse_args()  
    eval_flag = os.path.isfile(args.eval_data)
    model_flag = os.path.isfile(args.model_data)
    ### Load data
    if eval_flag:        
        sys.stderr.write("Loading Evaluation Data\n")
        eval_samples = load_samples(args.eval_data,args.interp_type)
        if len(eval_samples) == 1:
            sys.exit("Error: Specify >=2 samples for evaluation\n")
    if model_flag:
        sys.stderr.write("Loading Model Data\n")
        model_samples = load_samples(args.model_data,args.interp_type)  

    ### Import statements
    sys.stderr.write("Loading packages\n")
    method = setup_methods(args)

    ### Evaluation 
    if eval_flag and not model_flag:
        loso_evaluation(eval_samples,method,args)    
    else:
        train_test_evaluation(model_samples,eval_samples,method,args)
    return
    
    
### Create prediction file that has label file too.     

    
def setup_methods(args):
    cc,k,nj,it = args.classifier,args.k,args.num_jobs,args.interp_type
    if it == "SGF" or it == "WG":
        sys.stderr.write("DTW or related curve similarity metric must be used to compare curves with different numbers of features. Continuing with DTW calculation.\n")
        return "DTW-NN"
    if cc == "RF":
        return sklearn.ensemble.RandomForestClassifier(n_estimators=1001,n_jobs=nj)
    elif cc == "LR":
        return sklearn.linear_model.LogisticRegressionCV(max_iter=1001,n_jobs=nj)
    elif cc == "GBCT":
        return sklearn.ensemble.GradientBoostingClassifier(n_estimators=1001,)
    elif cc == "NN":
        return sklearn.neighbors.KNeighborsClassifier(n_neighbors=k,n_jobs=nj)
    elif cc == "NB":
        return sklearn.naive_bayes.GaussianNB()
    else: ## DTW-NN
        return "DTW-NN"    

def train_test_evaluation(train_samples,test_samples,method,args):
    X_train = list()
    for ts in train_samples:
        X_train+=ts.F            
    X_train = np.asarray(X_train)       
    Y_train = list()
    for ts in train_samples:
        Y_train+=ts.Y
    Y_train = np.asarray(Y_train)
    I_train = list()
    for ts in train_samples:
        I_train+=ts.I
    I_train = np.asarray(I_train)
    if args.equal_class:
        X_train,Y_train,I_train = create_random_equivalent_training_classes(X_train,Y_train,I_train)
    for i,test_sample in enumerate(test_samples):
        sys.stderr.write("Preparing evaluation for %s\n" % (test_sample.sample_name))
        X_test = np.asarray(test_sample.F)     
        Y_test = np.asarray(test_sample.Y)
        I_test = np.asarray(test_sample.I)
        D_test = np.asarray(test_sample.D)
        sys.stderr.write("Sample: %s, Training Labels:\n" % (test_sample.sample_name))
        print_label_counter(Y_train)        
        sys.stderr.write("Sample: %s, Testing Labels:\n" % (test_sample.sample_name))
        print_label_counter(Y_test)     
        if method == "DTW-NN":
            expression_prediction_dtw_nn(X_train,Y_train,I_train,X_test,Y_test,I_test,D_test,test_sample.sample_name,test_sample.header,args)
        else:
            expression_prediction_train_test(X_train,Y_train,I_train,X_test,Y_test,I_test,D_test,test_sample.sample_name,test_sample.header,method,args)    
    return
    

def loso_evaluation(samples,method,args):  
    for i,sample in enumerate(samples):
        sys.stderr.write("Preparing evaluation for %s\n" % (sample.sample_name))
        if args.obs_sample:
            train_samples = [s for j,s in enumerate(samples) if i!=j]
        else:
            train_samples = list()
            for s in samples:
                include = True
                for idsn in s.individual_sample_names:
                    if idsn in sample.individual_sample_names:
                        include = False
                if include:
                    train_samples.append(s)
        X_train = list()
        for ts in train_samples:
            X_train+=ts.F            
        X_train = np.asarray(X_train)
        X_test = np.asarray(sample.F)        
        Y_train = list()
        for ts in train_samples:
            Y_train+=ts.Y
        Y_train = np.asarray(Y_train)
        Y_test = np.asarray(sample.Y)
        I_train = list()
        for ts in train_samples:
            I_train+=ts.I
        I_train = np.asarray(I_train)
        I_test = np.asarray(sample.I)
        D_test = np.asarray(sample.D)
        sys.stderr.write("Sample: %s, Training Labels:\n" % (sample.sample_name))
        print_label_counter(Y_train)        
        sys.stderr.write("Sample: %s, Testing Labels:\n" % (sample.sample_name))
        print_label_counter(Y_test)
        if args.equal_class:
            X_train,Y_train,I_train = create_random_equivalent_training_classes(X_train,Y_train,I_train)
        sys.stderr.write("Sample: %s, Training Labels:\n" % (sample.sample_name))
        print_label_counter(Y_train)                
        if method == "DTW-NN":
            expression_prediction_dtw_nn(X_train,Y_train,I_train,X_test,Y_test,I_test,D_test,sample.sample_name,sample.header,args)
        else:
            expression_prediction_loso(X_train,Y_train,I_train,X_test,Y_test,I_test,D_test,sample.sample_name,sample.header,method,args) 
    return

def create_random_equivalent_training_classes(X,Y,I):
    counter = collections.Counter(Y)
    min_count = min(counter.values())
#    print min_count
    pos_lab_idx = np.where(Y==1)[0]
#    print pos_lab_idx
    neg_lab_idx = np.where(Y==-1)[0]
#    print neg_lab_idx
    pos_lab_randsub_idx = np.random.choice(pos_lab_idx,min_count,replace=False)
    neg_lab_randsub_idx = np.random.choice(neg_lab_idx,min_count,replace=False)
    rand_sub_idx = np.sort(np.concatenate((pos_lab_randsub_idx,neg_lab_randsub_idx),axis=0))
    return (X[rand_sub_idx],Y[rand_sub_idx],I[rand_sub_idx])


def expression_prediction_loso(X_train,Y_train,I_train,X_test,Y_test,I_test,D_test,sample_name,header,method,args):               
    func_name = get_function_name(method)    
#    outname = ".".join([sample_name,func_name,"png"])
    sys.stderr.write("Training %s for %s\n" % (func_name,sample_name))
    ### Need to perform Leave-one-fold-out CV to compare.
    n_folds = args.strat_k
    Y_pred_prob = [list() for x in range(len(Y_test))]
    Y_pred = [0]*len(Y_test)
    ### Randomly split indexes
    kf = sklearn.cross_validation.KFold(len(Y_test),n_folds=n_folds,shuffle=True)
    for train_idx,test_idx in kf:    
        I_kf_test = I_test[test_idx]
        I_kf_test_set = set(I_kf_test.tolist())
        X_train_loo = list()
        Y_train_loo = list()
        ## Need to exclude any examples of exactly the gene        
        for j,i_train in enumerate(I_train):
            if args.obs_gene:
                X_train_loo.append(X_train[j])
                Y_train_loo.append(Y_train[j])
            else:
                if not i_train in I_kf_test_set:
                    X_train_loo.append(X_train[j])
                    Y_train_loo.append(Y_train[j])
        X_train_loo = np.array(X_train_loo)
        Y_train_loo = np.array(Y_train_loo)
        method.fit(X_train_loo,Y_train_loo)
        X_test_loo = X_test[test_idx]
        Y_pred_prob_el= method.predict_proba(X_test_loo)
        Y_pred_el = method.predict(X_test_loo)
        for j,y_pred in enumerate(Y_pred_el):
            Y_pred_prob[test_idx[j]] = Y_pred_prob_el[j]
            Y_pred[test_idx[j]] = y_pred            
    ofh = open(".".join([sample_name,func_name,"pred"]),'w')  
    ofh.write(header)
    for j,p in enumerate(Y_pred_prob):
        l = [str(x) for x in p]
        ofh.write("\t".join(l+list(D_test[j]))+"\n")
    ofh.close()
    return


def expression_prediction_dtw_nn(X_train,Y_train,I_train,X_test,Y_test,I_test,D_test,sample_name,header,args):  
#    k_start = 1
#    k_stop = 151
#    k_step = 2
    ik = args.k 
    k_array = [ik]
#    k_array = range(k_start,k_stop+1,k_step)
    sys.stderr.write("Testing DTW for %s\n" % (sample_name))

    Y_pred_prob,Y_pred = predict_nearest_neighbor_dynamic_time_warping(X_train,Y_train,I_train,X_test,I_test,k_array,args)    
    for i,k in enumerate(k_array):
        func_name = "dtw_nn%d" % k
        ofh = open(".".join([sample_name,func_name,"pred"]),'w')  
        ofh.write(header)
        for j,p in enumerate(Y_pred_prob):
            l = [str(x) for x in p[0]]
            ofh.write("\t".join(l+list(D_test[j]))+"\n")
        ofh.close()
    return
   
def predict_nearest_neighbor_dynamic_time_warping(X_train,Y_train,I_train,X_test,I_test,k_array,args):   
    sys.stderr.write("Calculating DTW distances\n")
    Y_pred = list()
    Y_pred_prob = list()   
    pool = multiprocessing.Pool()      
    for q,x_test in enumerate(X_test):
#        c_test = C_test[q]       
#        sys.stderr.write("Now testing: %s\n" % I_test[q]) 
        if q%50 == 0:    
            sys.stderr.write("\tNow on gene %d/%d\n" % (q,len(X_test)))   
        ## Calculating methylation DTW
        compare_iter = list()
        for x_train in X_train:        
            compare_iter.append((x_test,x_train))
        XD = pool.map(wrap_mlpy_dtw,compare_iter)
        XD = np.asarray(XD)
#        XD = max_normalize(np.asarray(XD))
        ## Calculating CpG density DTW
#        compare_iter = list()
#        for c_train in C_train:
#            compare_iter.append((c_test,c_train))
#        CD = pool.map(wrap_mlpy_dtw,compare_iter)
#        CD = np.asarray(CD)
#        CD = max_normalize(np.asarray(CD))
#        Z = XD+CD
        
#        all_min_idx = np.argsort(Z)[0:k_array[-1]] 
#        all_min_idx = np.argsort(XD)[0:k_array[-1]] 
        all_min_idx = np.argsort(XD)
#        print all_min_idx
#        print Z[all_min_idx]
#        print XD[all_min_idx]
#        print CD[all_min_idx]
#        print I_train[all_min_idx]
        Y_pred_k_array = list()
        Y_pred_prob_k_array = list() 
        I_train_min_idx = I_train[all_min_idx]
        for k in k_array:
            ## So testing example is not overfit on it's own methylation signature,
            ## The testing example cannot see an example of itself.
            if args.obs_gene:
                min_idx = all_min_idx[:k]
            else:
                gc = 0
                min_idx = list()
                for j,i_train in enumerate(I_train_min_idx):
                    if gc == k:
                        break
                    elif i_train != I_test[q]:
                        min_idx.append(all_min_idx[j])
                        gc+=1
                    else:
                        continue
                min_idx = np.array(min_idx)
            min_y = Y_train[min_idx]
#            print min_y
            num_pos = (min_y == 1).sum()
            num_neg = (min_y == -1).sum()
            Y_pred_prob_k_array.append((float(num_neg)/k,float(num_pos)/k))
            if num_pos >= num_neg:
                Y_pred_k_array.append(1)
            else:
                Y_pred_k_array.append(-1)
        Y_pred.append(Y_pred_k_array)
        Y_pred_prob.append(Y_pred_prob_k_array)
#    print Y_pred[0]
#    print Y_pred_prob[0]
    return Y_pred_prob,Y_pred
    
    
    
    
    
### Expression Prediction for training & testing
def expression_prediction_train_test(X_train,Y_train,I_train,X_test,Y_test,I_test,D_test,sample_name,header,method,args):               
    func_name = get_function_name(method)    
#    outname = ".".join([sample_name,func_name,"png"])
    sys.stderr.write("Training %s for %s\n" % (func_name,sample_name))
    ### Need to perform Leave-one-fold-out CV to compare.
    n_folds = args.strat_k
    Y_pred_prob = [list() for x in range(len(Y_test))]
    Y_pred = [0]*len(Y_test)
    ### Randomly split indexes
    kf = sklearn.cross_validation.KFold(len(Y_test),n_folds=n_folds,shuffle=True)
    for train_idx,test_idx in kf:    
        I_kf_test = I_test[test_idx]
        I_kf_test_set = set(I_kf_test.tolist())
        X_train_loo = list()
        Y_train_loo = list()
        ## Need to exclude any examples of exactly the gene        
        for j,i_train in enumerate(I_train):
            if args.obs_gene:
                X_train_loo.append(X_train[j])
                Y_train_loo.append(Y_train[j])
            else:
                if not i_train in I_kf_test_set:
                    X_train_loo.append(X_train[j])
                    Y_train_loo.append(Y_train[j])
        X_train_loo = np.array(X_train_loo)
        Y_train_loo = np.array(Y_train_loo)
        method.fit(X_train_loo,Y_train_loo)
        X_test_loo = X_test[test_idx]
        Y_pred_prob_el= method.predict_proba(X_test_loo)
        Y_pred_el = method.predict(X_test_loo)
        for j,y_pred in enumerate(Y_pred_el):
            Y_pred_prob[test_idx[j]] = Y_pred_prob_el[j]
            Y_pred[test_idx[j]] = y_pred            
    ofh = open(".".join([sample_name,func_name,"pred"]),'w')  
    ofh.write(header)
    for j,p in enumerate(Y_pred_prob):
        l = [str(x) for x in p]
        ofh.write("\t".join(l+list(D_test[j]))+"\n")
    ofh.close()
    return

def max_normalize(arr):
    return arr/max(arr)
   
def wrap_mlpy_dtw(t):
#    return mlpy.dtw_std(t[0],t[1],dist_only=True)/float(len(t[0])+len(t[1]))   
    return mlpy.dtw_std(t[0],t[1],dist_only=True)   
    
def load_methylation(index,filename):
    line = linecache.getline(filename,index)
    return [float(x) for x in line.strip().split()]

def load_samples(dat_file_list,interp_type):
    samples = list()
    for line in open(dat_file_list,'r'):
        ll = line.strip().split()
        if len(ll) != 2:
            sys.exit("Error: Please specify a file with 2 columns: <label> <dat>")
        elif not ll[0].endswith("label"):
            sys.exit("Error: Please specify a file with 2 columns: <label> <dat>")
        elif not ll[1].endswith("dat"):
            sys.exit("Error: Please specify a file with 2 columns: <label> <dat>")
        else: 
            pass
        label_path = ll[0]
        dat_path = ll[1]
        sample_name = os.path.basename(dat_path).split(".")[0]            
        sys.stderr.write("Processing sample %s\n" % (sample_name))            
        Y,L,E,I,D,header = load_multilabel_expr(label_path)
        if interp_type == "ROI":
            X = load_vals_roi(dat_path)
        elif interp_type == "DMR":
            X = load_vals_dmr(dat_path)
        else:
            X = load_vals(dat_path)
        C = list()
        sample = Sample(sample_name,X,Y,L,E,C,I,D,header)
        sample.combine_features()
        samples.append(sample)      
        print_label_counter(sample.Y)
        sample.reduce_memory()
    return samples
    
def print_label_counter(Y):
    counter = collections.Counter(Y)
    for k in sorted(counter.keys()):
        sys.stderr.write("Class: %s, #: %d\n" % (str(k),counter[k]))
    sys.stderr.write("\n")
    return

def load_vals(filename):
    sys.stderr.write("Loading %s\n" % filename) 
    fh = open(filename,'r')
    C = list()
    for line in fh:
        if line.startswith("#"):
            continue
        ll = [float(x) for x in line.strip().split()]
        C.append(ll)
    return C
        
def load_vals_dmr(filename):    
    fh = open(filename,'r')
    X = list()
    for line in fh:
        if line.startswith("#"):
            continue
        Xl = list()
        for x in line.strip().split():
            try:
                Xl.append(float(x))                
            except ValueError:
                Xl.append('NaN')
        X.append(Xl)
    imputer = sklearn.preprocessing.Imputer()
    X = np.array(X)
    return imputer.fit_transform(X).tolist()

def load_vals_roi(filename):    
    fh = open(filename,'r')
    X = list()
    for line in fh:
        if line.startswith("#"):
            continue
        Xl = list()
        for x in line.strip().split():
            try:
                Xl.append(float(x))                
            except ValueError:
                Xl.append('NaN')
        X.append(Xl)
    imputer = sklearn.preprocessing.Imputer()
    X = np.array(X)
    return imputer.fit_transform(X).tolist()     
        
def load_multilabel_expr(filename):
    fh = open(filename,'r')
    Ex = list()
    L = list()
    E = list()
    Y = list()
    I = list()
    D = list()
    for line in fh:
        if line.startswith("#"):
            ll = line.lstrip("#").rstrip().split()
            ii = ll.index("GENE_ID")            
            ei = ll.index("EXPR")
            xi = ll.index("NUM_EXONS")
            pli = ll.index("POS_LOW")            
            phi = ll.index("POS_HIGH")
            header = "#PROB_DW\tPROB_UP\t"+line.lstrip("#")
            continue
        ll = line.strip().split()
        D.append(ll)
        id_name = ll[ii]
        expr = float(ll[ei])
        E.append(int(ll[xi]))
        L.append(math.log(int(ll[phi])-int(ll[pli]),10))
        Ex.append(expr)
        I.append(id_name)
        if expr > 0:
            Y.append(1)
        else:
            Y.append(-1)
    return Y,L,E,I,D,header

def get_function_name(func):
    return str(func).split('(')[0]

def setup_parser_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,description='''
For classification of gene expression based on interpolated methylation values over multiple samples. Must specify at least either model or evaluation data. ''')
    ### Required positional arguments
    parser.add_argument('interp_type',choices={"TSS","WG","WSG","UGF","ROI","DMR"},help="Interpolation type of input .dat files. Must use DTW-NN for WG and UGF.")
    parser.add_argument('eval_data',default="",help="2 column with path to <sample>.label and <sample>.<interp_type>.meth.dat files.")
    parser.add_argument('-m','--model_data',default="",help="2 column with path to <sample>.label and <sample>.<interp_type>.meth.dat files.")
    parser.add_argument('-c','--classifier',choices={"RF","LR","GBCT","NN","DTW-NN","NB"},default="RF",help="Choice of classifier: Random Forest (RF), Logistic Regression (LR), Gradient Boosted Classification Trees (GBCT), Nearest Neighbor - L2 Distance (NN), Nearest Neighbors - Dynamic Time Warping (DTW-NN), Naive Bayes (NB)")
    parser.add_argument('-j','--num_jobs',default=1,type=int,help="Number of jobs for multithreaded methods. Use -1 for all cores.")    
    parser.add_argument('-g','--obs_gene',type=bool,default=False,help="Include examples of testing genes in the training set")
    parser.add_argument('-s','--obs_sample',type=bool,default=False,help="Include individual sample from training set if observed in testing set")
    parser.add_argument('-k','--k',default=21,type=int,help="k for k-Nearest Neighbor classifiers")
    parser.add_argument('-t','--strat_k',default=10,type=int,help="k for Stratified K-fold cross validation with default evalution framework")
    parser.add_argument('-q','--equal_class',default=True,type=bool,help="Randomly subsample the dominant class to create equal classes")
    
    return parser

if __name__ == "__main__":
    import sys, os, argparse, linecache, math, collections
    import sklearn.preprocessing
    import numpy as np
    import sklearn.cross_validation  
    import sklearn.ensemble
    import sklearn.linear_model    
    import sklearn.ensemble    
    import sklearn.neighbors    
    import sklearn.naive_bayes
    import multiprocessing
    #import mlpy    
    main()    
    
