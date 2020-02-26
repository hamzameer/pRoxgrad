# # overlapping group lasso structure test
# #
# # #
# ng=50  # number of groups
# g_size=100  # group size
# overlap=10  # number of variables overlapped between two consecutive groups
# n=1000      # sample size
# p=ng*(g_size-overlap)+overlap   # total number of variables
# gamma=ng/10   # regularization parameter for group penalty L2-norm
# lambda=ng/10  # regularization parameter for L1-norm
# mu = 1e-3
# #### Settings #########
# #
# toy_group = gentoy_group(n, ng, g_size, overlap)
# X = toy_group$X
# Y = toy_group$Y
# Tm = toy_group$Tm
# Tw = toy_group$Tw
# w = toy_group$w
# #
# # #
# pre = pre_group(Tm, Tw)
# C = pre$C
# g_idx = pre$g_idx
# CNorm = pre$Cnorm
# option = NULL
# C1 <- as.matrix(C)
# 
# solve = SPG(prob = 'group', Y = Y, X = X, gamma = gamma, lambda = lambda, C= C, CNorm = Cnorm, g_idx = g_idx)
# 
# 
# #
# 
# # option.maxiter=10000;
# # max iteration
# 
# # option.mu=1e-3;
# # Smoothing Parameter
# # Carefully choose mu, a larger mu may not converge while a smaller one
# # leads to slow convergence
# # mu = (eps/2D) | D = |G|/2  | From the paper
# 
# # option.verbose=true;
# # option.display_iter=10;
# # option.tol=1e-8;          # tolerance
# 
# # generate toy data
# # X: design matrix
# # Y: output
# # T: 0/1 sparse matrix: # of groups by # of features,
# # indicate each group contains which variables
# # Tw: # of groups by 1 column vector,  weight for each group
# # w: true regression coefficients
# #
# #
# 
# 
# #
# # prob='group'
# 
# # SPG with a pre-compuated Lipschitz constant for medium-scale problems
# 
# ### Function 4 | [grad_beta,grad_obj,grad_density,grad_iter,grad_time, L] = SPG_linesearch(prob, Y, X, gamma, lambda, C,option, g_idx);
# 
# # SPG with line search on Lipschitz constant for large-scale problems
