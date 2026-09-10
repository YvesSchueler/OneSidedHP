# -*- coding: utf-8 -*-
"""
Created on Fri May  5 13:41:37 2023

@author: elias
"""


import numpy as np
from scipy.optimize import minimize


"""
 Code for "On adjusting the one-sided Hodrick-Prescott filter" (2020) 
 by Elias Wolf (FU-Berlin), Frieder Mokinski (Deutsche Bundesbank), and
 Yves Schüler (Deutsche Bundesbank)

 Version date 2023/05/04
 If you encounter any bug, please mail Yves Schüler at yves.schueler (at) bundesbank.de

 This program estimates the cyclical component of the adjusted one-sided HP
 filter (HP-1s*) by estimating the two-sided HP filter (HP-2s) on an expanding sample and keeping the last observation.

 For the exact references that are cited in this program, please see
 the paper.

       Input:    
                   y     - a Tx1 vector (y_1,...,y_T)', where T is the
                   number of observations. 
                       
                   lm_2  - a scalar. This is the value of the smoothing
                   parameter usually employed for the two-sided HP filter
                   (HP-2s). This input is optional. The default is lm_2=1,600.

                   opt    - 0/1. If opt=1, the adjustment parameters are
                   obtained through the optimization routine (computationally more intensive). This input is optional. The default is
                   opt=0 (computationally less intensive).

                   sample - 0/1. If sample = 1, the adjustment takes sample size into account. May be used for large smoothing
                   parameters in combination with small sample sizes, for instance, if
                   lm_2 = 400,000 and T < 100. 

      Output:  
                   ycycle_adj  - a Tx1 vector with the adjusted extracted cyclical component
                  
                   lm_1        - a scalar. This is the adjusted value of the smoothing parameter used as an input to the one-sided HP filter

                   kappa       - a scalar. This is the mulitplicative scaling factor, with wich the cyclical component of the
                   one-sided HP filter is rescaled

       Example:
                   ycycle_adj, lm_1, kappa = adj1s_hpfilter(y,1600)
                   yields a Tx1 vector of the extracted cyclical component
                   using a smoothing parameter of 650 and scaling
                   parameter of size 1.1513 for the one-sided HP filter
                  
"""




def adj1s_hpfilter(y, lm_2, opt=True, sample=False):
    
    if isinstance(y, (np.ndarray)) == False:
        
        y = np.array(y)
    
    if len(y) == 0:

        raise Exception('Time series is missing or empty')

    if len(y.shape) > 1:

        raise Exception('y contains more than one time series: the program is designed for one series only')

    if y.shape[0] < 3:

        raise Exception('y contains less than 3 observations: please check your data input')
        
    if np.isnan(y).any() == True:

        raise Exception('y contains missing values: please remove')

    # Fitted polynomials for corresponding lambda and kappa values
    poly_lm = lambda lm_2: -5.071743 +  0.4080780*lm_2  -0.749537*lm_2**(1/2) + 14.50001*lm_2**(1/3) -50.36123*lm_2**(1/4) + 41.44369*lm_2**(1/5)

    poly_kappa = lambda lm_2: 0.984928 +  0.547209*(1/lm_2)  - 0.703636*(1/np.sqrt(lm_2)) + 3.291657*(1/(lm_2)**(1/3)) - 3.165415*(1/(lm_2)**(1/4)) + 1.759527 *(1/(lm_2)**(1/5))
            
    if opt == False:

        lm_1 = poly_lm(lm_2)
        kappa = poly_kappa(lm_2)

    
    if opt == True:
    
    # finds adjustment parameters lm_1 and kappa by harmonizing PTF of HP-2s
    # and HP-1s
    
        if sample == True:
            T_eff = y.size
        else:
            T_eff = 1000  #size of vector of filter weight for observation t=T
    
        def ptf_cost(params, lm_2, T):
    
            # Parameters to be optimized
            lm_1 = params[0]
            ratio = params[1]
    
            I = np.eye(T)
            Q_t = np.diff(np.eye(T), 2, axis=0)
            A_inv = np.linalg.inv(I + lm_1*np.dot(Q_t.T, Q_t))
    
            wts_one = A_inv[-1,:]
            
            # FFT of the weights to obtain PTF of the one-side HP-filter
            grid = 2*np.pi*np.asarray(range(int(np.floor(T/2))+1))/T
    
            k = np.asarray(range(-(T-1), 1))
            k = np.transpose(k)
            hpgain = []
    
            for i in range(1, len(grid)+1):
    
                hpgain.append(1 - np.ones(shape=(1, len(k))).dot(wts_one.T*np.exp(1j*grid[i-1]*k)))
    
            hp1_ptf = np.concatenate(ratio**2*np.abs(hpgain)**2)
    
            # PTF of the two-sided HP_filter 
            hp2_ptf = ((4*lm_2*(1-np.cos(grid))**2/(1+4*lm_2*(1-np.cos(grid))**2))**2)
    
            # Calculate loss function defined as the squared distance    
            cost = sum((hp1_ptf-hp2_ptf)**2)
    
            return(cost)
    
        # Initial values and optimization of the cost function
        init_vals = np.asarray([poly_lm(lm_2), poly_kappa(lm_2)])
        res = minimize(ptf_cost, init_vals, args=(lm_2, T_eff), method='BFGS',
                       options={'disp': True})
    
        # Collect optimal parameters
        lm_1 = res.x[0]
        kappa = res.x[1]
        
    
    def hp_one(lm_1, kappa, y):

        y = np.asarray(y)
        T_series =  len(y)
        cycle_os = []

        for i in range(3,T_series+1):
            
            I = np.eye(i)
            Q_t = np.diff(I, 2, axis=0)
            A_inv = np.linalg.inv(I + lm_1*np.dot(Q_t.T, Q_t))
            
            psi_t = y[i-1] - A_inv[-1,:].dot(y[:i])
            
            cycle_os.append(psi_t)

        return(kappa*np.asarray(cycle_os))
    
    # Apply the one-sided HP-Filter to the series
    ycycle_adj = hp_one(lm_1, kappa, y)

    return([ycycle_adj, lm_1, kappa])


