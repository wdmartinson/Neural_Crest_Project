import numpy as np
# import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances
from sklearn.cluster import KMeans

def compute_inertia(a, X):
    # W = [np.mean(pairwise_distances(X[a == c, :])) for c in np.unique(a)]
    W = [0.5*np.mean(pairwise_distances(X[a == c, :])**2) for c in np.unique(a)]
    # return np.mean(W)
    return np.sum(W)

def bounding_box(X):
    num_cols = X.shape[1]
    limits = [None]*num_cols*2 # min and max
    for k in range(num_cols):
        limits[2*k] = min(X,key=lambda a:a[k])[k]
        limits[2*k+1] = max(X,key=lambda a:a[k])[k]
    # ymin, ymax = min(X,key=lambda a:a[1])[1], max(X,key=lambda a:a[1])[1]
    # return xmin, xmax, ymin, ymax
    return limits

def compute_gap(clustering, data, k_max=5, n_references=50):
    if len(data.shape) == 1:
        data = data.reshape(-1, 1)
    reference = np.random.rand(*data.shape,n_references) # Generate n_references uniform random distributions
    # xmin, xmax, ymin, ymax = bounding_box(data)
    # reference[:,0] = xmin + (xmax-xmin)*reference[:,0]
    # reference[:,1] = ymin + (ymax-ymin)*reference[:,1]
    limits = bounding_box(data)
    num_cols = data.shape[1]
    reference_inertia = [None]*k_max
    std_err = [None]*k_max
    for k in range(1, k_max+1):
        local_inertia = [None]*n_references
        for ii in range(n_references):
            clustering.n_clusters = k
            for j in range(num_cols):
                reference[:,j,ii] = limits[2*j] + (limits[2*j+1]-limits[2*j])*reference[:,j,ii] # Make sure the reference distribution is in the same space as your data
            assignments = clustering.fit_predict(reference[:,:,ii])
            local_inertia[ii] = compute_inertia(assignments, reference[:,:,ii])
        reference_inertia[k-1] = np.mean(np.log(local_inertia))
        std_err[k-1] = np.sqrt(1+1.0/n_references)*np.std(np.log(local_inertia), ddof = 0)

    ondata_inertia = [None]*k_max
    for k in range(1, k_max+1):
        clustering.n_clusters = k
        assignments = clustering.fit_predict(data)
        ondata_inertia[k-1] = compute_inertia(assignments, data)

    gap = reference_inertia-np.log(ondata_inertia)
    return gap, reference_inertia, np.log(ondata_inertia), std_err


# Load in the positions of the cells at the final time point:
points = np.loadtxt('final_cell_positions_puncta_model.csv', delimiter = ',')
# plt.plot(points[:,0], points[:,1], 'o')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()

if len(points) <= 10:
    k_max = max(len(points)-1 , 1) # Upper limit of clusters is pre-set to 10, for performance reasons
else:
    k_max = 10 # Upper limit of clusters is pre-set to 10, for performance reasons
gap, reference_inertia, ondata_inertia, std_err = compute_gap(KMeans(), points, k_max=k_max, n_references=k_max)
std_err = np.asarray(std_err)

gap_statistic = np.argmax(gap)+1
file = open('gap_statistic_puncta_model.txt', 'w')
file.write(str(gap_statistic))
file.close()

# Get list of all gap values that satisfy gap(k)-gap(k+1)+std_err(k+1) >= 0
gap_std_1_err = np.zeros(len(gap))
for i in range(1,k_max):
    if gap[i-1]-gap[i]+std_err[i] >= 0:
        gap_std_1_err[i-1] = gap[i-1]

other_gap_statistic = next((i+1 for i, x in enumerate(gap_std_1_err) if x), None)
if other_gap_statistic == None:
    other_gap_statistic = gap_statistic
file = open('gap_std_1_error_statistic_puncta_model.txt', 'w')
file.write(str(other_gap_statistic))
file.close()

# plt.plot(points[:,0], points[:,1], 'o')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()
#
# plt.plot(range(1, k_max+1), reference_inertia,
#          '-o', label='reference')
# plt.plot(range(1, k_max+1), ondata_inertia,
#          '-o', label='data')
# plt.xlabel('k')
# plt.ylabel('log(inertia)')
# plt.show()
#
# plt.plot(range(1, k_max+1), gap, '-o')
# plt.ylabel('gap')
# plt.xlabel('k')
# plt.show()
#
# plt.errorbar(range(1, k_max+1), gap, fmt='-o', yerr = std_err)
# plt.ylabel('gap')
# plt.xlabel('k')
# plt.show()
