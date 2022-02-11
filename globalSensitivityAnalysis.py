import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt

# These methods are applicable to measurement data, where there is no control over all or some of the variables.
# There are three methods:
## 1. Multiple Linear regression
## 2. HDMR

def gsa_mlr(x,y):
    ''' Method:	Global sensitivity analysis with Multiple Linear Regression

        Author:	Jari Vepsalainen
        Updated:	02/02/2022

        References: 	[1] Groen, E. A., & Heijungs, R. (2017). Ignoring correlation in uncertainty and sensitivity analysis in life cycle assessment: what is the risk?. Environmental Impact Assessment Review, 62, 98-109.
                        [2] Xu, C., & Gertner, G. Z. (2008). Uncertainty and sensitivity analysis for models with correlated parameters. Reliability Engineering & System Safety, 93(10), 1563-1573.
                        [3] Li, L., Lu, Z., & Zhou, C. (2011). Importance analysis for models with correlated input variables by the state dependent parameters method. Computers & Mathematics with Applications, 62(12), 4547-4556.

        input: X (input(s)) and Y (output) as pandas dataframes.
        output: df, R2, bHat, RMSE
    '''
    X       = x.values
    Y       = y.values
    names   = x.columns.values
    
    N       = X.shape[0]  # Number of rows (observations)
    k       = X.shape[1]  # Number of columns (input factors)
    Yi      = np.array(Y) # of the output
    vHat    = np.zeros(k)
    vHatU   = np.zeros(k)
    vHatC   = np.zeros(k) 
    
    for i in range(k): 
        # Partial variance of each input parameter [1],[2] 
        R           = np.ones((N, 2))
        R[:,1]      = X[:,i]
        thetaHat    = (np.linalg.inv(R.T@R)) @ (R.T @ Yi)
        yHat        = thetaHat[0]+thetaHat[1]*X[:,i]
        vHat[i]     = (1/(N-1))*(np.sum((yHat - np.mean(Yi))**2))
        
        # Uncorrelated and correlated variance contributions [1],[2]
        C           = X
        C           = np.delete(C,i,1)
        Reta        = np.concatenate((np.ones((N,1)), C), axis=1) 
        etaHat      = (np.linalg.inv(Reta.T@Reta)) @ (Reta.T @ X[:,i])

        zHat        = X[:,i] - (etaHat[0]+(C@etaHat[1:]))  

        Ar          = np.concatenate((np.ones((N,1)), zHat[np.newaxis,:].T), axis=1)
        r           = (np.linalg.inv(Ar.T@Ar)) @ (Ar.T @ Yi)
        yHatU       = Ar@r

        vHatU[i]    = (1/(N-1))*(np.sum((yHatU - np.mean(Yi))**2))
        vHatC[i]    = vHat[i]-vHatU[i]
        
    # 2nd Order Correlation variance [3]
    cT      = np.array(list(itertools.permutations(np.arange(k), 2)))
    vHatlq  = np.zeros(cT.shape[0])
    for i in range(cT.shape[0]):
        D           = X
        D           = np.delete(D,cT[i],1) # Remove columns l and q
        Reta2       = np.concatenate((np.ones((N,1)), D), axis=1)
        etaHat2     = (np.linalg.inv(Reta2.T@Reta2)) @ (Reta2.T @ X[:,cT[i][0]])
        zHat2       = X[:,cT[i][0]] - (etaHat2[0]+(D@etaHat2[1:]))
            
        Ar2         = np.concatenate((np.ones((N,1)), zHat2[np.newaxis,:].T), axis=1) # should here be zhat2 or not?
        r2          = (np.linalg.inv(Ar2.T@Ar2)) @ (Ar2.T @ X[:,cT[i][1]])
        yHatlq      = Ar2@r2
        vHatlq[i]   = (1/(N-1))*(np.sum((yHatlq - np.mean(X[:,cT[i][1]]))**2))

    # Importance matrix [3]
    S       = np.diagflat(vHatU/np.var(Yi))
    for i in range(cT.shape[0]):
        S[cT[i][0],cT[i][1]]    =  vHatlq[i]/np.var(X[:,cT[i,1]])
    S       = (S+S.T) - np.eye(k)*np.diag(S) 
    np.fill_diagonal(S,0) # Set diagonal as zeros
    for i in range(k):
        S[:,i]=S[:,i]/np.sum(S[:,i])

    # Total Variance of the system and the R-squared value [2]
    R1      = np.concatenate((np.ones((N,1)),X), axis=1)
    bHat    = (np.linalg.inv(R1.T@R1)) @ (R1.T @ Yi)
    yHatL   = np.zeros(N)
    for i in range(N):
        yHatL[i]    = bHat[0]+sum(np.reshape(bHat[1:], bHat[1:].shape[0])*X[i,:])
    vHatL     = (1/(N-1))*(np.sum((yHatL - np.mean(Yi))**2)) 
    R2      = vHatL/np.var(Yi)
    RMSE    = np.sqrt(np.mean(np.power(Yi-yHatL,2)))

    # Sensitivity Indices [2][3]
    St   = vHat/np.var(Yi)
    Su   = vHatU/np.var(Yi)
    Sc   = vHatC/np.var(Yi)
    df   = pd.DataFrame(np.array([Su, Sc]).T, columns = ['Main effect', 'Corr. effect'], index = names)
    df.index.name="Factor"
    return df, R2, bHat, RMSE


from SALib.analyze import hdmr

def gsa_hdmr(x,y):
    '''
    Calculating the sensitivity indices based on RS-HDMR approach which considers correlated inputs.

    Based on the paper: Li, G., Rabitz, H., Yelvington, P. E., Oluwole, O. O., Bacon, F., Kolb, C. E., & Schoendorf, J. (2010). 
    Global sensitivity analysis for systems with independent and/or correlated inputs. The journal of physical chemistry A, 114(19), 6022-6032.
    
    input: X (input(s)) and Y (output) as pandas dataframes.
    output: df of sensitivity indices
    '''

    X       = x.values
    Y       = np.reshape(y.values, y.values.shape[0]) # of the output
    names   = x.columns.values
    k       = X.shape[1]  # Number of columns (input factors)

    bounds = np.zeros((k, 2))
    for i in range (k):
        bounds[i,:] = [min(X[:,i]), max(X[:,i])]
    
    problem = {
    'num_vars': k,
    'names': names,
    'bounds': bounds
    }

    Si = hdmr.analyze(problem, X, Y, print_to_console=False)

    return Si.to_df()

def plot_hdmr(df, X):
    # Sensitivity indicies
    k = X.shape[1]
    S = (df.S.values)[0:k]
    St = (df.ST.values)[0:k]
    S_err = (df.S_conf.values)[0:k]
    St_err = (df.ST_conf.values)[0:k]
    names=(df.index.values)[0:k]

    # First-order and total
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    X = np.arange(len(names))
    ax.bar(X + 0.00, St, color = 'b', width = 0.25, label='Total SI', yerr=St_err, capsize=5)
    ax.bar(X + 0.25, S, color = 'g', width = 0.25, label='First-order SI', yerr=S_err, capsize=5)
    plt.xticks(X + 0.25 / 2, names)
    plt.ylabel('Sensitivity Index [-]')
    plt.xlabel('Factor')
    ax.set_axisbelow(True)
    ax.grid(color='gray', linestyle='dashed')
    ax.legend()
    
    return S
