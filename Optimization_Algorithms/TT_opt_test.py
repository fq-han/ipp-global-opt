"""The demo of using ttopt. Simple example with cache usage.

We'll find the minimum for the simple analytic function of the form f(X), where
input X is the [samples, dimension] numpy array using the cache.

Run it from the root of the project as "python demo/cache.py".

As a result of the script work we expect the output in console like this:
"
...
Simple-5d  | evals=9.80e+03+3.01e+03 | t_cur=6.47e-02 | e_x=3.55e-02 e_y=1.03e-04
----------------------------------------------------------------------
Simple-5d  | evals=9.80e+03+3.01e+03 | t_all=7.61e-02 | e_x=3.55e-02 e_y=1.03e-04 
"

"""
import numpy as np
from scipy.optimize import rosen


from ttopt import TTOpt
from ttopt import ttopt_init

np.random.seed(42)
 
def f_example1(X):
    d = X.shape[1] 
    yy = 1 + np.sum((X[:, 0:]+0.1)**2/4000, axis=1) - np.prod(np.cos((X[:, 0:]+0.1)/np.sqrt(np.arange(1, d + 1))), axis = 1)
    return yy
# Levy
def f_example2(X):
    d = X.shape[1]  # Number of dimensions
    
    wj = 1 + (X[:, :-1] - 0.1) / 4
    sum_term = np.sum((wj - 1) ** 2 * (1 + 10 * np.sin(np.pi * (wj + 1)) ** 2), axis=1)
    
    w1 = 1 + (X[:, 0] - 0.1) / 4
    wd = 1 + (X[:, d-1] - 0.1) / 4 
    return (sum_term + np.sin(np.pi * w1) ** 2 +
          (wd - 1) ** 2 * (1 + np.sin(2 * np.pi * wd) ** 2)) 

# Rasrigin
def f_example3(X):
    d = X.shape[1]  # Number of dimensions
    yz_sum = np.sum((X + 0.1) ** 2 - 10 * np.cos(2 * np.pi * (X + 0.1)), axis=1)
    return (yz_sum + 10 * d) 

# Ackley
def f_example5(X):
    d = X.shape[1]  # Number of dimensions
    yz1 = np.sum((X + 0.0) ** 2, axis=1)
    yz2 = np.sum(np.cos(2 * np.pi * (X + 0.0)), axis=1)
    yz = -20 * np.exp(-0.2 * np.sqrt(yz1 / d)) - np.exp(yz2 / d) + 20 + np.exp(1)
    return yz 

# Corrugated spring
def f_example6(X):
    yx = 0
    xsum = np.sum(X**2, axis=1)
    yx = yx + 0.1*xsum - np.cos(5*np.sqrt(xsum))
    yx = yx + 1
    return yx

# Brown
def f_example8(X):
    X = np.array(X)  # Ensure X is a NumPy array
    d = X.shape[1]  # Number of dimensions
    
    # Handle potential invalid values during power operation
    with np.errstate(all='ignore'):
        exponent1 = 2 * (X[:, 1:] ** 2 + 1)
        exponent2 = 2 * (X[:, :-1] ** 2 + 1)
        
        term1 = np.power(X[:, :-1] ** 2, exponent1)
        term2 = np.power(X[:, 1:] ** 2, exponent2)
        
        yz = np.sum(term1 + term2, axis=1)
    
        return yz 
# Exponential
def f_example9(X):
    d = X.shape[1]  # Number of dimensions
    yz = np.sum(X ** 2, axis = 1)
    return (-np.exp(-0.5 * yz * d) +1)

# Trid Function
def f_example10(X):
    d = X.shape[1]  # Number of dimensions
    yx = (X[:, 0] - 1 + d) ** 2
    for jd in range(1, d):
        yx += (X[:, jd] - 1 + jd * (d + 1 - jd)) ** 2 - (X[:, jd - 1] + (jd - 1) * (d + 1 - jd + 1)) * (X[:, jd] + jd * (d + 1 - jd))
    yx += d * (d + 4) * (d - 1) / 6
    yx *= 2 * d
    return yx

# Discuss Function
def f_example11(X):
    d = X.shape[1]  # Number of dimensions
    N = 4
    yx = X[:, 0] ** 2 * N
    for jd in range(1, d):
        yx += X[:, jd] ** 2
    yx = yx 
    return yx

#Cosine mixture
def f_example12(X):
    d = X.shape[1]  # Number of dimensions
    yx = 0
    for jd in range(d):
        yx += -0.1 * np.cos(5 * np.pi * X[:, jd] / 20) - (X[:, jd] / 20) ** 2
    yx += 0.1 * d
    yx /= 2
    return yx

# Alpine 1 Function
def f_example13(X):
    d = X.shape[1]  # Number of dimensions
    yx = 0
    for jd in range(d):
        yx += np.abs(X[:, jd] * np.sin(X[:, jd]) + 0.1 * X[:, jd])
    return yx

# Example 14
def f_example14(X):
    d = X.shape[1]  # Number of dimensions
    yx = 0
    for jd in range(d):
        yx += 0.429 * X[:, jd] - 1.126 * X[:, jd] ** 2 - 0.143 * X[:, jd] ** 3 + 0.563 * X[:, jd] ** 4
    yx += 0.8490 * d
    return yx

def f_example15(xs):
    d = xs.shape[1]  # Number of dimensions
    yx = 0.5 * d
    for jd in range(d - 1):
        sx = xs[:, jd]**2 + xs[:, jd + 1]**2
        yx += (np.sin(np.sqrt(sx))**2 - 0.5) / (1 + 0.001 * sx)**2
    return yx
   
def f_example31(xs):
    # Define the covariance matrix
    # Sigma = np.array([[2.250, 0.300, 1.500, 2.250],
    #                   [0.300, 4.000, 3.500, 2.400],
    #                   [1.500, 3.500, 6.250, 6.000],
    #                   [2.250, 2.400, 6.000, 9.000]])
    # Initialize a d x d matrix with zeros
    d = 10
    Sigma = np.zeros((10, 10))

# Fill each column of A
    for j in range(1, d + 1):  # Python is 0-indexed, so range starts at 1 for equivalence
         Sigma[:, j - 1] = np.exp(-((np.arange(1, d + 1) - j) ** 2) / 2)
    # Number of assets
    n = Sigma.shape[0]

    # Target risk contribution (equal for all assets)
    target_RC = 1 / n

    # Initialize an array to store objective values for each row of xs
    y_values = []

    # Loop over each sample in xs (each row)
    for row in range(xs.shape[0]):
        # Extract current weight vector
        w = xs[row, :]
        
        # Compute portfolio risk as sqrt(w' * Sigma * w)
        portfolio_risk = 0
        for i in range(n):
            for j in range(n):
                portfolio_risk += w[i] * Sigma[i, j] * w[j]
        portfolio_risk = np.sqrt(portfolio_risk)
        
        # Compute risk contributions
        RC = np.zeros(n)
        for i in range(n):
            RC[i] = w[i] * sum(Sigma[i, j] * w[j] for j in range(n)) / portfolio_risk
        
        # Compute the objective function
        y = 0
        for i in range(n):
            y += (RC[i] - target_RC * portfolio_risk) ** 2
        
        # Add penalty term for the sum of weights
        y += 2*abs(1 - np.sum(w))
        
        # Store the objective value for this sample
        y_values.append(y)
    
    return np.array(y_values)

d = 10                     # Number of function dimensions:
rank = 10                    # Maximum TT-rank while cross-like iterations

# We initialize the TTOpt class instance with the correct parameters:
tto = TTOpt(
    f=f_example31,                    # Function for minimization. X is [samples, dim]
    d=d,                    # Number of function dimensions
    a=0.,                  # Grid lower bound (number or list of len d)
    b=1.,                  # Grid upper bound (number or list of len d)
    n=2**8,                 # Number of grid points (number or list of len d)   2**12*3 if 1e-2
    evals=1E+6,            # Number of function evaluations max 1E+6
    name='Simple',          # Function name for log (this is optional)
    x_opt_real=np.zeros(d), # Real value of x-minima (x; this is for test)
    y_opt_real=0.,          # Real value of y-minima (y=f(x); this is for test)
    with_cache=True,        # We save all requests into cache
    with_log=True)


# And now we launching the minimizer:
tto.optimize(rank)


# We can extract the results of the computation:
x = tto.x_opt          # The found value of the minimum of the function (x)
y = tto.y_opt          # The found value of the minimum of the function (y=f(x))
k_c = tto.k_cache      # Total number of cache usage (should be 0 in this demo)
k_e = tto.k_evals      # Total number of requests to func (is always = evals)
k_t = tto.k_total      # Total number of requests (k_cache + k_evals)
t_f = tto.t_evals_mean # Average time spent to real function call for 1 point
                       # ... (see "ttopt.py" and docs for more details)


# We log the final state:
print('-' * 70 + '\n' + tto.info() +'\n\n')

# We log the final state:
print('y_opt : ', y)
print('i_opt : ', x)
