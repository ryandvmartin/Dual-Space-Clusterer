# (c) Ryan Martin 2017 under MIT license

import numpy as np
import pygeostat as gs
import matplotlib.pyplot as plt


class AutocorrClus:
    """
    Parameters:
    mvdata: pd.DataFrame
        the multivariate dataframe
    locations: pd.DataFrame
        the Cartesian locations of the data
    target_nclus: int
        the number of clusters
    clustermethod: str
        can be one of 'kmeans', 'gmm', or 'hier'
    autocor: str
        the autocorrelation metric to consider, options are 'morans', 'getis', 'localmean', or
        'localvar'
    nnears: int
        the number of nearest neighbors in a local search to consider
    aniso_search: 5-tuple
        (ang1, ang2, ang3, r1, r2) anisotropic properties following the GSLIB rotation conventions

    .. codeauthor:: Ryan Martin - 17-10-2017
    """
    def __init__(self, mvdata, locations, target_nclus=None, clustermethod='kmeans',
                 autocor='getis', nnears=45, aniso_search=(0, 0, 0, 1, 1)):
        # parse the input, make sure arrays are nvar x ndata dimensioned
        if hasattr(mvdata, 'values'):
            mvdata = mvdata.values.T
        elif isinstance(mvdata, np.ndarray):
            if mvdata.shape[0] > mvdata.shape[1]:
                mvdata = mvdata.T
        if hasattr(locations, 'values'):
            locations = locations.values.T
        elif isinstance(locations, np.ndarray):
            if locations.shape[0] > locations.shape[1]:
                locations = locations.T

        # set all the parmeters to this object
        self.mvdata = mvdata
        self.locations = locations
        self.target_nclus = target_nclus
        self.clustermethod = clustermethod
        self.aniso_search = aniso_search
        self.nnears = nnears
        autocordict = dict(morans=1, getis=2, localvar=3, localmean=4)  # values for fortran
        self.autocor = autocordict[autocor]

    def fit_predict(self, target_nclus=None):
        """
        Call all the functions to get a cluster using the autocorrelation metrics and the data
        passed to the constructor of this class
        """
        from pandas import DataFrame
        import pygeostat as gs
        from autocorr import autocorr

        # initialize the spatial locations for autocorrelation
        autocorr.initautocorr(self.locations, self.aniso_search)
        # get the autocorrelation dataset
        tdata = autocorr.localautocorrelation(self.autocor, 1, self.nnears, self.mvdata).T

        tcols = ['T%i' % i for i in range(tdata.shape[1])]
        tdf = DataFrame(tdata, columns=tcols)
        tdata, _ = gs.nscore(tdf, tcols)
        tdata = tdata.values
        if target_nclus is not None:
            self.target_nclus = target_nclus
        assert self.target_nclus is not None, " set a number of clusters ! "
        clusdefs = cluster(self.target_nclus, tdata, algorithm=self.clustermethod)
        return clusdefs


class MDSPlot:
    """
    Wrapper for MDS plotting for a given dataset, passed on init
    Parameters:
        data (dataframe or np.ndarray): the data to embed with MDS
        mtype (str): the embedding type, valid are 'mds' and 'tsne'
        embedding (int): the number of dimensions to embed the coordinates in
        randstate (int): the random state for embedding
    """
    def __init__(self, data, mtype='MDS', embedding=2, randstate=69069, verbose=False,
                 autofit=True):
        self.data = data
        self.mtype = mtype
        self.embedding = embedding
        self.randstate = randstate
        self.verbose = verbose
        self.coords = None

    def embed(self, data=None):
        """
        Quick wrapper around sklearn manifold for MDS and TSNE
        Parameters:
            data (dataframe or np.ndarray):
        """
        from sklearn import manifold
        if data is None:
            data = self.data
        np.random.RandomState(self.randstate)
        if self.mtype.lower() == 'mds':
            model = manifold.MDS(n_components=self.embedding, random_state=self.randstate,
                                 verbose=self.verbose)
        if self.mtype.lower() == 'tsne':
            model = manifold.t_sne.TSNE(n_components=self.embedding, random_state=self.randstate,
                                        verbose=self.verbose)
        nrow, ncol = np.shape(data)
        if ncol > nrow:
            data = data.T
        model = model.fit(data)
        self.coords = model.embedding_  # the N-dimensional embedding of the data

    def plot_mds(self, colors='k', cmap='viridis', catdata=False, pltstyle='pt8',
                 ax=None, cax=None, figsize=(2.5, 3), s=15, lw=0.1, cbar=None, grid=True,
                 legend_loc='lower right', title=None, vlim=None, legstr='Cluster',
                 xlabel=None, ylabel=None):
        """
        Plotting for the MDS coordinates and the probabilities or categories
        """
        gs.set_style(pltstyle)
        coords = self.coords
        # setup the figure
        fig, ax, cax = gs.setup_plot(ax, cbar=cbar, cax=cax, figsize=figsize, grid=grid)
        if vlim is None and colors is not None and not isinstance(colors, str):
            vlim = (np.min(colors), np.max(colors))
        else:
            vlim = (None, None)
        # deal with non-array input
        if hasattr(colors, 'values'):
            colors = colors.values
        if len(np.unique(colors)) <= 12:
            catdata = True
        # plot categories
        if catdata:
            clusids = np.unique(colors)
            ncat = len(clusids)
            if cmap in gs.cat_palettes:
                cmap = gs.get_palette(cmap, ncat)
            else:
                cmap = gs.catcmapfromcontinuous(cmap, ncat)
            dump = colors.copy()
            for i in range(ncat):
                dump[colors == clusids[i]] = np.arange(ncat)[i]
            colors = dump
            ax.scatter(coords[:, 0], coords[:, 1], c=colors, cmap=cmap, s=s, lw=lw,
                       label=legstr)
            if isinstance(legend_loc, str):
                ax.legend(loc=legend_loc, scatterpoints=1, handletextpad=0.1)
            elif isinstance(legend_loc, tuple):
                ax.legend(loc='upper left', bbox_to_anchor=legend_loc, scatterpoints=1,
                          handletextpad=0.1)

        # plot continous data with a colorbar
        else:
            plot = ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=s, lw=lw, cmap=cmap,
                              vmin=vlim[0], vmax=vlim[1])
            if cbar:
                vlims, ticklocs, ticklabels = gs.get_contcbarargs(colors, 2, vlim, nticks=8)
                cbar = fig.colorbar(plot, cax=cax, ticks=ticklocs)
                cbar.ax.set_yticklabels(ticklabels, ha='left')

        applytickpad(ax, 1)
        addgridtoplt(ax)
        if ylabel is None:
            ax.set_ylabel('$MDS_2$')
        else:
            ax.set_ylabel(ylabel)
        if xlabel is None:
            ax.set_xlabel('$MDS_1$')
        else:
            ax.set_xlabel(xlabel)
        applylabelpad(ax, 0.5)
        if title:
            ax.set_title(title)
        return ax


def cluster(nclus, dfs, n_init=1, variables=None, algorithm='kmeans', ret_gmm=False):
    """
    Wrapper around Gaussian Mixture and Kmeans

    Parameters:
        nclus (int): number of clusters
        dfs (pd.DataFrame): the dataframe
        n_init (int): the number of times to initialize the clusterer
        variables (list): the list of variables to take from the dataframe
        algorithm (str): clustering type, `kmeans` or `gmm`
        ret_gmm (bool): either `true` for returning the mixture or false to return just the
            clustering

    Returns:
        clusdef +- gmm depending on ret_gmm
    """
    from sklearn.cluster import KMeans, hierarchical
    from sklearn.mixture import GaussianMixture
    # get the data for clustering:
    if variables is not None:
        nd_pts = dfs[variables].values
    else:
        nd_pts = dfs
    nclus = int(nclus)
    if algorithm.lower() == 'kmeans':
        clusdef = KMeans(n_clusters=nclus, n_init=n_init).fit(nd_pts).labels_
    elif algorithm.lower() == 'gmm':
        gmm = GaussianMixture(n_components=nclus, n_init=n_init)
        gmm = gmm.fit(nd_pts)
        clusdef = np.array(gmm.predict(nd_pts))
    elif algorithm.lower() == 'hier':
        hier = hierarchical.AgglomerativeClustering(n_clusters=nclus, linkage='ward')
        clusdef = hier.fit_predict(nd_pts)
    if ret_gmm:
        return clusdef, gmm
    else:
        return clusdef


def fixlabels(ax, tickpad=1.55, labelpad=0.45):
    applylabelpad(ax, labelpad)
    applytickpad(ax, tickpad)


def label_subplot(ax, fignum=None, labels=True, case='lower', annot_xy=(0.0, 1.03), ha='left',
                  va='bottom', fontsize=10, **kwargs):
    """
    Use the module variables and the default properties listed herein to label the axis given
    the ax and the figure number
    NEW:
    Optionally omit the figure number, and set kwarg `labels` to type str, will label the string at
    the given location!
    """
    """ set of subplot labels for label_subplot() """
    figurelabels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)',
                    '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)', '(p)',
                    '(q)', '(r)', '(s)', '(t)']
    figurenumbers = ['(i)', '(ii)', '(iii)', '(iv)', '(v)', '(vi)', '(vii)',
                     '(viii)', '(ix)', '(x)', '(xi)', '(xii)', '(xiii)', '(xiv)']
    if isinstance(labels, str):
        lbl = labels
    else:
        if labels:
            lbls = figurelabels
        else:
            lbls = figurenumbers
        if case == 'lower':
            lbl = lbls[fignum].lower()
        else:
            lbl = lbls[fignum].upper()
    ax.annotate(lbl, xy=annot_xy, ha=ha, va=va, xycoords='axes fraction', fontsize=fontsize,
                **kwargs)


def reclass_cluster_single(clus1, clus2):
    """
    Use the a pairwise similarity score to determine the best reclassification of cluster 2
    relative to clustering 1

    Note: directly calls the fortran codes which calculates the similarity metrics using
        permutations

    Parameters
    ----------
    clus1: np.ndarray
        the first array of clusterings
    clus2: np.ndarray
        the second array of clusterings

    Returns
    -------
    reclus: np.ndarray
        the relclassified clus2
    diff: float
        the matching score between clus1 and the final reclassified clusters ... not the jaccard!

    """
    from clus_compare import clus_compare_mod
    reclus, diff = clus_compare_mod.reclass_cluster_single(clus1, clus2)
    return reclus, diff


def simple_density_calc(x, y, ncomp, niter=50, seed=69069):
    from sklearn.mixture import GaussianMixture
    if hasattr(x, 'values'):
        x = x.values
    if hasattr(y, 'values'):
        y = y.values
    sample = np.c_[x, y].T
    gmm = GaussianMixture(n_components=ncomp, n_init=niter)
    gmm = gmm.fit(sample)
    return gmm.score_samples(sample)


def standardize(data, variables=None):
    """
    Subtract the mean and standardize by the stdev, returning a dataframe with the new variables,
    and the means and standard deviations used to normalize.
    """
    if variables is not None:
        data = data[variables]
    return_df = False
    if hasattr(data, 'values'):
        arr = data.values
        return_df = True
    ncol = np.size(arr, 1)
    means = []
    stds = []
    for i in range(ncol):
        mean = np.mean(arr[:, i])
        std = np.std(arr[:, i])
        arr[:, i] -= mean
        arr[:, i] /= std
        means.append(mean)
        stds.append(std)
    if return_df:
        import pandas as pd
        columns = ['std_%s' % var for var in data.columns]
        arr = pd.DataFrame(arr, columns=columns)
    return arr, means, stds


def addgridtoplt(ax=None):
    """ Add the standard pygeostat grid to the current plot"""
    if ax is None:
        fig = plt.gcf()
        ax = fig.axes[0]
    ax.grid(which='major', color='k', linestyle=':', lw=0.2)


def applylabelpad(ax, labelpad):
    ax.yaxis.labelpad = labelpad
    ax.xaxis.labelpad = labelpad


def applytickpad(ax, tickpad):
    ax.tick_params(axis='both', which='major', pad=tickpad)


