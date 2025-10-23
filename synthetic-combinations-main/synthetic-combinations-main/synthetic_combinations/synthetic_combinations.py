from typing import List
import numpy as np
from celer import LassoCV
from tqdm import tqdm
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold
from scipy.linalg import hadamard

class Synthetic_Combinations():
    """
    Impute counterfactuals for combinatorial interventions.
    """
    
    def __init__(self,max_rank = None,spectral_t = None):
        """
        Parameters
        ----------
        max_rank : int 
        Perform truncated SVD on training data with this value as its rank 
        spectral_t : float 
        Perform truncated SVD on training data with (100*thresh)% of spectral energy retained. 
        If omitted, then the default value is chosen via Donoho & Gavish '14 paper. 
        """
        self.max_rank = max_rank
        self.spectral_t = spectral_t
    
    def horizontal_fit(self,observation_matrix,donor_unit_indices,horizontal_estimator = LassoCV(cv = 5)): 
        """
        Does the horizontal regression via the lasso with regularization parameter chosen by cross-validation
        
        Parameters
        ----------
        observation_matrix: numpy array
        Observed matrix with missing entries denoted by NA
        
        donor_unit_indices: list
        list of row indices that correspond to donor units
        """
        fourier_characteristic_matrix = hadamard(observation_matrix.shape[1]) # get matrix of fourier coefficients. 2^number of items
        for donor_unit in tqdm(donor_unit_indices): #perform horizontal regression for each donor unit
            donor_unit_outcomes = observation_matrix[donor_unit,:]
            non_nan_indices = np.argwhere(~np.isnan(donor_unit_outcomes)) #get non-missing entries
            non_nan_indices = [non_nan_indices[index][0] for index in range(len(non_nan_indices))] #convert to list
                
            donor_unit_fourier_characteristic_matrix = fourier_characteristic_matrix[non_nan_indices,:] #X matrix used for donor unit
            donor_unit_observed_outcomes = donor_unit_outcomes[non_nan_indices] #observed outcomes
            lasso_reg = horizontal_estimator.fit(donor_unit_fourier_characteristic_matrix,donor_unit_observed_outcomes) #fit lasso
               
            observation_matrix[donor_unit,:] = lasso_reg.predict(fourier_characteristic_matrix) #predict outcomes
        return observation_matrix
       
     
    def vertical_fit(self,horizontal_imputed_observation_matrix,donor_unit_indices,use_cv = True):
        """
        Does the vertical regression via PCR
        
        Parameters
        ----------
        horizontal_imputed_observation_matrix: numpy array
        Observed matrix with missing entries for donor units imputed. 
        
        donor_unit_indices: list
        list of row indices that correspond to donor units
        """
        N = horizontal_imputed_observation_matrix.shape[0] #Number of units
        indices = set([i for i in range(N)])
        non_donor_unit_indices = list(indices.difference(set(donor_unit_indices))) #get indices for non-donor units
        for n in tqdm(non_donor_unit_indices): #impute outcomes for all non-donor units
            n_outcomes = horizontal_imputed_observation_matrix[n,:]
            non_nan_indices = np.argwhere(~np.isnan(n_outcomes))
            non_nan_indices = [non_nan_indices[index][0] for index in range(len(non_nan_indices))]
        
            n_observed_outcomes = n_outcomes[non_nan_indices] #observed outcomes for non-donor unit 
            donor_unit_outcomes = horizontal_imputed_observation_matrix[donor_unit_indices,:] 
            donor_unit_n_outcomes = donor_unit_outcomes[:,non_nan_indices] #imputed outcomes for donor set corresponding to observations of non-donor unit
            if use_cv:
                w, u_rank, s_rank, v_rank = self._pcr_cv(donor_unit_n_outcomes.T,n_observed_outcomes) #get coefficients from PCR CV
            else:
                w, u_rank, s_rank, v_rank = self._pcr(donor_unit_n_outcomes.T,n_observed_outcomes) 
            n_preds = np.matmul(horizontal_imputed_observation_matrix[donor_unit_indices,:].T,w)
            horizontal_imputed_observation_matrix[n,:] = n_preds
        return horizontal_imputed_observation_matrix
    
    
        
      
    def _pcr(self, X, y):
        """
        principal component regression (PCR) 
        """
        (u, s, v) = np.linalg.svd(X, full_matrices=False) 
        if self.max_rank is not None: 
            rank = self.max_rank 
        elif self.spectral_t is not None: 
            rank = self._spectral_rank(s)
        else: 
            (m, n) = X.shape 
            rank = self._universal_rank(s, ratio=m/n) #choose num components according to Donoho & Gavish '14
        s_rank = s[:rank]
        u_rank = u[:, :rank]
        v_rank = v[:rank, :] 
        beta = ((v_rank.T/s_rank) @ u_rank.T) @ y #compute regression vector
        return (beta, u_rank, s_rank, v_rank)
    
    def _pcr_cv(self,X,y,num_folds = 5):
        """
        do principal component regression (PCR) choosing the number of components via CV 
        """
        (u, s, v) = np.linalg.svd(X, full_matrices=False)
        rank = len(s) #rank of X
        kf = KFold(n_splits=num_folds) #for splitting data
        cv_error = np.empty((rank,num_folds))
        cv_error[:] = np.nan
        fold_num = 0
        for i, (train_index, test_index) in enumerate(kf.split(X)): # K-fold CV split
            X_train = X[train_index,:] 
            y_train = y[train_index]
            X_test = X[test_index,:]
            y_test = y[test_index] 
            for r in range(1,rank+1): #iterate through all ranks
                (u_train, s_train, v_train) = np.linalg.svd(X_train, full_matrices=False)
                s_train_r = s_train[:r]
                u_train_r = u_train[:,:r]
                v_train_r = v_train[:r,:]
                beta_train_r = ((v_train_r.T/s_train_r) @ u_train_r.T) @ y_train
                test_preds = np.matmul(X_test,beta_train_r)
                cv_error[r-1,fold_num] = r2_score(y_test,test_preds)
            fold_num += 1
        average_error = np.nanmean(cv_error,axis = 0)
        opt_rank = np.argmax(average_error)
        
        s_rank = s[:opt_rank]
        u_rank = u[:, :opt_rank]
        v_rank = v[:opt_rank, :] 
        beta = ((v_rank.T/s_rank) @ u_rank.T) @ y
        return (beta, u_rank, s_rank, v_rank)
        
    
    def _spectral_rank(self, s):
        if self.spectral_t==1.0: 
            rank = len(s)
        else: 
            total_energy = (s**2).cumsum() / (s**2).sum()
            rank = list((total_energy>self.spectral_t)).index(True) + 1
        return rank

    def _universal_rank(self, s, ratio): 
        """
        retain all singular values above optimal threshold as per Donoho & Gavish '14:
        https://arxiv.org/pdf/1305.5870.pdf
        """ 
        omega = 0.56*ratio**3 - 0.95*ratio**2 + 1.43 + 1.82*ratio
        t = omega * np.median(s) 
        rank = max(len(s[s>t]), 1)
        return rank 