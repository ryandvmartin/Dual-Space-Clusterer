# (c) Ryan Martin 2017 under MIT license

import numpy as np
import numba
from scipy.spatial import cKDTree
from math import sin, cos, pi
from popdiff import population_difference_mod as popdiff_mod


class AgglomCluster:
    """
    Clustering by the Dual Space Search (DS) agglomeration

    Parameters
    ----------
    mvdata: np.ndarray
        the data containing variables in nvar x ndata
    locations: np.ndarray
        the data containing locations in nvar x ndata
    niter: int
        number of clustering iterations to run
    nnears: int
        number of nearest neighbors to consider in the spatial search
    stage1merge: int
        the number of neighbors merged in stage 1
    seed: int
        the seed to set the random state
    ansio: 5-tuple
        containing (ang1, ang2, ang3, r1, r2) following GSLIB rotation ordering

    .. codeauthor:: Ryan Martin - 2017
    """

    def __init__(self, mvdata, locations, niter=100, nnears=10, stage1merge=5, seed=69069,
                 aniso=(0, 0, 0, 1, 1)):
        """
        Init the pars of the class
        """
        # some checks
        assert mvdata.ndim == 2, "mvdata must be 2 dimensional with nd x dimension"
        assert mvdata.shape[0] > mvdata.shape[1], "mvdata.shape[0] must be > mvdata.shape[1]"
        assert locations.ndim == 2, "locations must be 2 dimensional with nd x dimension"
        assert locations.shape[0] > locations.shape[1], "mvdata.shape[0] must be > mvdata.shape[1]"
        assert isinstance(aniso, (tuple, list)), 'Pass either list or tuple for `aniso`'
        assert len(aniso) == 5, 'Must pass a 5-long tuple of (ang1, ang2, ang3, r1, r2) for `aniso`'
        # assign stuff to the class
        self.mvdata = mvdata
        self.locations = locations
        self.numpairings = None
        self.pars = dict(niter=niter, nnears=nnears, seed=seed,
                         aniso=aniso, stage1merge=stage1merge)

    def fit(self, target_nclus, minprop=0.01, maxprop=0.65, verbose=True, ipdiff=0):
        """
        Fit the clustering ensemble to the multivariate and spatial locations using the dual-stage
        search

        Parameters
        ----------
        target_nclus: int
            the target number of clusters to merge to, could be larger than the target number of
            clusters for the domain so that this parameter can be inferred from the pairings matrix
        minprop: float
            the minimum proportion allowed in a single cluster
        maxprop: float
            the maximum proprtion allowed in a single cluster
        verbose: bool
            whether or not to write a progress bar to the jupyter notebook
        ipdiff: int
            0: randomized(default), 1:centers, 2: energy distance, 3: m-energy distance,
            4: wards distance

        Returns
        -------
        clusterings: np.ndarray

        """
        seed = self.pars["seed"]
        niter = self.pars["niter"]
        nnears = self.pars["nnears"]
        stage1merge = self.pars["stage1merge"]
        aniso = self.pars["aniso"]
        mvdata = self.mvdata
        locations = self.locations
        self.nclus = target_nclus

        # dimensioning checks
        if hasattr(mvdata, 'values'):
            mvdata = mvdata.values
        ndata, nvar = mvdata.shape
        if hasattr(locations, 'values'):
            locations = locations.values
        ndata1, dim = locations.shape
        assert isinstance(target_nclus, int), "target_nclus must be an integer"

        # generate the random seeds
        np.random.RandomState(seed)
        iterable = np.random.randint(20000, high=50000, size=niter)
        self.clusterings = np.zeros((ndata, niter))
        irun = 0
        for seed in iterable:
            args = (mvdata, locations, nnears, stage1merge, target_nclus,
                    aniso, int(seed), minprop, maxprop)
            self.clusterings[:, irun] = agglomclus_single(*args)
            irun += 1

    def predict(self, target_nclus=None, method='ward'):
        """
        return the class labels from the clustering of the data using ensemble methods to get the
        single realization that best describes the clustering of the dataset
        method (str): `single`, `complete`, `average`, `ward`, `weighted`,
        """
        if target_nclus is None:
            target_nclus = self.nclus
        self.numpairings = _getpairings(self.clusterings)
        clustering = get_hierarchy(self.numpairings, target_nclus, method)
        return clustering


def agglomclus_single(mvdata, xyzlocs, nnears, num_take, target_nclus, aniso=(0, 0, 0, 1, 1),
                      rseed=69069, minprop=0.001, maxprop=0.999, ipdiff=0):
    """
    Find a single agglomerated clustering
    """
    from scipy.spatial import cKDTree
    np.random.RandomState(rseed)
    ndata, nvar = mvdata.shape
    ndata1, dim = xyzlocs.shape
    mvdata_t = mvdata.T
    assigned = np.zeros(ndata, dtype=bool)
    anisolocs = get_rot_coords(aniso[0:3], [1.0, *aniso[3:5]], xyzlocs)
    aniso_tree = cKDTree(anisolocs)
    assigned[:] = False
    idx = np.random.permutation(ndata)
    clusterlabels = np.zeros(ndata)
    iclus = 0
    # ------ primary merging stage, based on the spatial neighbors -----------------------------
    for ix in idx:
        # query the kdtree around the location ix
        _, didx = aniso_tree.query(anisolocs[ix, :], k=nnears + 1)
        didx = didx[1:]
        mvdist = np.zeros(nnears)
        # for each neighbor, compute the sqr eucldist and sort based on closest
        for i, jx in enumerate(didx):
            mvdist[i] = sqr_euclidean(mvdata[ix, :], mvdata[jx, :])
        didx = didx[np.argsort(mvdist)]
        # merge closest nearby points
        nassigned = 0
        for i, jx in enumerate(didx):
            if assigned[ix] and not assigned[jx]:
                clusterlabels[jx] = clusterlabels[ix]
                assigned[jx] = True
            elif assigned[jx] and not assigned[ix]:
                clusterlabels[ix] = clusterlabels[jx]
                assigned[ix] = True
            elif not assigned[ix] and not assigned[jx]:
                clusterlabels[ix] = clusterlabels[jx] = iclus
                assigned[ix] = assigned[jx] = True
                iclus += 1
            nassigned += 1
            if nassigned >= num_take:
                break
    # ------ secondary merging based on the MV differences of the populations ------------------
    nclus_left = 1e21
    while nclus_left > target_nclus + 1:
        clusters = np.unique(clusterlabels)
        nclus_left = len(clusters)
        iclus = np.max(clusters) + 1
        if ipdiff == 0:
            pdiff = np.random.randint(1, 4)
        else:
            pdiff = ipdiff
        # choose a random cluster center, store the data indexes to labels as clusix
        ix = np.random.randint(0, nclus_left)
        clusix = clusterlabels == clusters[ix]
        x = mvdata_t[:, clusix]
        cx = np.mean(x, axis=1)
        jxs = np.array([i for i in range(nclus_left) if i != ix])
        # pre-sort the populations based on the mean to limit the amount of pop-diffing
        centerdistances = []
        for jx in jxs:
            y = mvdata_t[:, clusterlabels == clusters[jx]]
            cy = np.mean(y, axis=1)
            centerdistances.append(sqr_euclidean(cx, cy))
        jxs = jxs[np.argsort(centerdistances)]
        mv_div_score = []
        for jx in jxs[:min(5, nclus_left)]:
            y = mvdata_t[:, clusterlabels == clusters[jx]]
            if pdiff == 1:
                cx = np.mean(x, axis=1)
                cy = np.mean(y, axis=1)
                diffs = sqr_euclidean(cx, cy)
            elif pdiff == 2:
                diffs = popdiff_mod.energy_statistic(x, y)
            elif pdiff == 3:
                diffs = popdiff_mod.energy_statistic_mdist(x, y, regconst=1e-4)
            elif pdiff == 4:
                diffs = popdiff_mod.wards_distance(x, y)
            elif pdiff == 5:
                diffs = popdiff_mod.kd_distance(x, y, band=0.65)
            else:
                print("ERROR: invalid popdiff")
            mv_div_score.append(diffs)
        # get the jx corresponding to the smalled MV distance
        sortidx = np.argsort(mv_div_score)
        jx = jxs[sortidx[0]]
        clusjx = clusterlabels == clusters[jx]
        # now we have to merge clusix and clusjx
        clusterlabels[clusix] = clusterlabels[clusjx] = iclus
    return clusterlabels


@numba.jit
def sqr_euclidean(pt1, pt2):
    """ n-dimensional squared distance between pt1 and pt2 """
    dis = 0.0
    for i in range(pt1.shape[0]):
        dis += (pt1[i] - pt2[i]) ** 2
    return dis


def get_rot_coords(a, r, coords):
    """
    Translation of the GSLIB setrot function, a are the rotation angles ang1, ang2, ang3 in degrees
    """
    arads = np.deg2rad(a)
    rotmat = np.zeros((3, 3))
    sina = sin(arads[0])
    sinb = sin(arads[1])
    sint = sin(arads[2])
    cosa = cos(arads[0])
    cosb = cos(arads[1])
    cost = cos(arads[2])
    # Construct the rotation matrix:
    r1 = r[0] / r[1]
    r2 = r[0] / r[2]
    rotmat[0, :] = [cosb * cosa, cosb * sina, -sinb]
    rotmat[1, :] = [r1 * (-cost * sina + sint * sinb * cosa),
                    r1 * (cost * cosa + sint * sinb * sina),
                    r1 * (sint * cosb)]
    rotmat[2, :] = [r2 * (sint * sina + cost * sinb * cosa),
                    r2 * (-sint * cosa + cost * sinb * sina),
                    r2 * (cost * cosb)]
    dim = coords.shape[1]
    return np.dot(rotmat[:dim, :dim], coords.T).T


@numba.jit
def total_wcss(clustering, mvdata):
    """
    Calculate the within cluster sum of squares

    Parameters:
        clustering (np.ndarray): 1D array of integer codes corresponding to a cluster assignment
        mvdata (np.ndarray): N-dimensional array (ndata x dim) array of values to compute the sum
            of squares

    Returns:
        sum_of_squares (float): the within-cluster-sum-of-squares
    """
    categories = np.unique(clustering)
    nclus = categories.shape[0]
    ndata, nvar = mvdata.shape
    wcss = np.zeros(nclus)
    for i in range(nclus):
        didx = clustering == categories[i]
        dslice = mvdata[didx, :]
        cx = np.mean(dslice, axis=0)
        for j in range(dslice.shape[0]):
            wcss[i] += sqr_euclidean(cx, dslice[j, :])
    return wcss.sum()


@numba.jit
def total_local_entropy(xyzlocs, clustering, knears=25):
    """
    Function to get the local entropy considering the number of nearest neighbors
    and the categories defined in the input array

    Parameters:
        xyzlocs (np.ndarray): locations of the sample points, MV space
        clustering (np.ndarray): cluster definitions (or categories) for each location
        knears (int): number of nearest neighbors in the search
    Returns:

    """
    categories = np.unique(clustering)
    ncats = categories.shape[0]
    assert xyzlocs.shape[0] > xyzlocs.shape[1], " xyzlocs must be ndata x dim dimensioned "
    ndata, dim = xyzlocs.shape
    search_tree = cKDTree(xyzlocs)
    pks = np.zeros(ncats)
    hsum = 0.0
    for iloc in range(ndata):
        _, idx = search_tree.query(xyzlocs[iloc, :], k=knears)
        for icat, cat in enumerate(categories):
            count = 0
            for lcat in clustering[idx]:
                if cat == lcat:
                    count += 1
            pks[icat] = np.maximum(count / knears, 1.0e-10)
        hsum += np.sum(pks * np.log(pks))
    return hsum


def calc_cluster_stats(mvdata, locations, clusdefs, nnears):
    # make sure things are properly dimensioned
    assert mvdata.shape[0] > mvdata.shape[1], "transpose mvdata"
    assert locations.shape[0] > locations.shape[1], "transpose locations!"
    wcss = total_wcss(clusdefs, mvdata)
    spentropy = total_local_entropy(locations, clusdefs, nnears)
    return wcss, spentropy


def label_cluster_stats(mvdata, locations, clusdefs, nnears, anisos, ax, top=True, fontsize=7,
                        coords=None, **kwargs):
    from functions import label_subplot
    if hasattr(mvdata, 'values'):
        mvdata = mvdata.values
    if hasattr(locations, 'values'):
        locations = locations.values
    if hasattr(clusdefs, 'values'):
        clusdefs = clusdefs.values
    wcss, spentropy = calc_cluster_stats(mvdata, locations, clusdefs, nnears)
    annot = 'WCSS: %.1f\nSE: %.1f' % (wcss, spentropy)
    if top:
        if coords is None:
            coords = (0.01, 0.99)
        va = 'top'
    else:
        if coords is None:
            coords = (0.01, 0.01)
        va = 'bottom'
    label_subplot(ax, labels=annot, fontsize=fontsize, annot_xy=coords, va=va, **kwargs)


@numba.jit(nopython=True)
def _getpairings(clusterings):
    """ counts the number of times i was paired with j over all realizations """
    ndata, nreal = clusterings.shape
    numpairings = np.zeros((ndata, ndata))
    for ireal in range(nreal):
        cluscodes = clusterings[:, ireal]
        for i in range(ndata):
            for j in range(ndata):
                if cluscodes[i] == cluscodes[j]:
                    numpairings[i, j] += 1
    numpairings = numpairings / float(nreal)
    return numpairings


def get_hierarchy(pairingsmatrix, target_nclus, method='ward'):
    """
    Return the hierarchical clustering of the nd x nd pairings matrix truncated to the
    target number of clusters
    """
    from scipy import cluster
    linkage = cluster.hierarchy.linkage(pairingsmatrix, method=method)
    clusdefs = cluster.hierarchy.fcluster(linkage, target_nclus, criterion='maxclust')
    return clusdefs
