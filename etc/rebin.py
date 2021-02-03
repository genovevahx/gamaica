import numpy as np


def get_borders(bins):
    '''Quick and dirty function to convert bins to borders

    In the ETC, we assume that bins are specified by their central
    value. For flux and transmisson invariant rebinning however we
    need the borders instead. The borders are now simply assumed to be
    in the middle of the central bins, which is a good assumption as
    long as the bin width doesn't change too much between neighbored
    bins. The leftmost and rightmost borders are chosen so that the
    values are really central.

    Example::

        >>> import astropy.units as u
        >>> bins = np.arange(360., 365., 1.) * u.nm
        >>> print(bins)
        [360. 361. 362. 363. 364.] nm
        >>> print(get_borders(bins))
        [359.5 360.5 361.5 362.5 363.5 364.5] nm

    Parameters
    ----------
    bins : numpy.ndarray
        Central bin values

    Returns
    -------
    numpy.ndarray
        Border bin values, with length ``len(bins)+1``

    '''
    borders = np.resize(bins, len(bins)+1)
    borders[1:-1] = 0.5*(bins[1:] + bins[:-1])
    borders[0] = 1.5 * bins[0] - 0.5 * bins[1]
    borders[-1] = 1.5 * bins[-1] - 0.5 * bins[-2]
    return borders


def rebin_1d_flux_box(flux, org_bins, new_bins):
    '''Flux-invariant re-binning.

    This function takes n flux values and their central bins.  It
    rebins the flux into a sequence of m new flux values, defined by
    their central bins.

    The flux on the part of the bins that are covered by both original
    and new borders is preserved. If the new borders cover a smaller
    interval, the uncovered part of the flux will be lost, hence
    reducing the total flux. Bins in the new borders that are not
    covered by the original bins are set to 0 flux.

    This function should work for both, increasing and decreasing the
    number or size of bins.

    For this function, it is assumed, that the flux within one bin is
    uniformly distributed over the range of that bin. That means, re-
    binning to bins of significantly smaller size than in the original
    sequence will lead to plateaus on the new flux sequence.

    Preserved flux means, that if a bin with the borders a and b, a <
    b has a flux of f, any smaller bin (c, d) within the borders of a
    and b (a <= c < d <= b) will also the flux that falls into the
    smaller bin f*(d-c)/(b-a). Also, if the new borders completely
    cover the old borders, sum(flux) = sum(new_flux), safe for
    numerical errors.

    This function can become numerically instable if elements in
    org_bins are similar but not identical to elements in new_bins.

    Example ::

        >>> import astropy.units as u
        >>> flux = np.array([1., 5., 1., 1., 1.])
        >>> bin0 = np.arange(360., 365., 1.) * u.nm
        >>> bin1 = np.arange(360., 365., 2.) * u.nm
        >>> print(rebin_1d_flux_box(flux, bin0, bin1))
        [3.5 4.  1.5]

    Parameters
    ----------
    flux : np.ndarray
        Flux values as a 1-dim array with at least two entries.

    org_bins : np.ndarray
        Original bins (central values). The length must be equal to
        the length of the ``flux`` array.

    new_bins : np.ndarray
        New bins (central values). The length must be at least 2.

    Returns
    -------
    np.ndarray
        Rebinned flux values. The length is equal to the length
        of the ``new_bins`` parameter.

    '''
    org_borders = get_borders(org_bins)
    new_borders = get_borders(new_bins)
    
    # int_borders means interlinked borders. The interlinked borders are used
    # to calculate how the flux in the old bins is distributed to the new
    # bins. The interlinked borders contain both the new and old borders.
    # union1d returns a sorted array of unique values.

    int_borders = np.union1d(org_borders.astype(new_borders.dtype),
                             new_borders)
    
    # Make sure that only the intersecting area of both original and new
    # borders is used.
    int_borders = int_borders[(int_borders >= org_borders[0])
                              & (int_borders <= org_borders[-1])]
    int_borders = int_borders[(int_borders >= new_borders[0])
                              & (int_borders <= new_borders[-1])]

    # The centers of the interlinked bins are used to find the indices within
    # the original and new borders.
    int_centers = 0.5 * (int_borders[1:] + int_borders[:-1])
    org_ids = np.digitize(int_centers, org_borders, right=True) - 1
    new_ids = np.digitize(int_centers, new_borders, right=True) - 1

    # The interlinked bucket width is calculated relative to the original
    # bucket width. The list of interlinked borders defines buckets are always
    # entirely within one original bucket. The flux intensity in the new
    # buckets is then the old intensity times the new bucket size divided by
    # the old bucket size. So its a simple multiplication with a ratio of
    # bucket sizes.
    org_width = org_borders[1:] - org_borders[:-1]
    int_width = int_borders[1:] - int_borders[:-1]
    int_flux = flux[org_ids] * int_width / org_width[org_ids]

    # bincount sums all weight values according to the list of indices in x.
    # See examples in the documentation of numpy. minlength needs to be there
    # in order to cover bins of the new flux that are larger than the right
    # border of org_borders.
    return np.bincount(new_ids, weights=int_flux,
                       minlength=len(new_borders) - 1)


def rebin_1d_trans_box(trans, org_bins, new_bins):

    '''Transmission-invariant re-binning

    Example ::

        >>> import astropy.units as u
        >>> trans = np.array([0.8, 0.1, 0.8, 0.8, 0.8])
        >>> bin0 = np.arange(360., 365., 1.) * u.nm
        >>> bin1 = np.arange(360., 365., 2.) * u.nm
        >>> print(rebin_1d_trans_box(trans, bin0, bin1))
        [0.625 0.625 0.8  ]

    Parameters
    ----------
    trans : np.ndarray
        Transmission values as a 1-dim array with at least two entries.

    org_bins : np.ndarray
        Original bins (central values). The length must be equal to
        the length of the ``trans`` array.

    new_bins : np.ndarray
        New bins (central values). The length must be at least 2.

    Returns
    -------
    np.ndarray
        Rebinned transmission values. The length is equal to the length
        of the ``new_bins`` parameter.

    '''
    org_borders = get_borders(org_bins)
    new_borders = get_borders(new_bins)

    # org_binwidth = org_borders[1:] - org_borders[:-1]
    new_binwidth = new_borders[1:] - new_borders[:-1]

    # int_borders means interlinked borders. The interlinked borders
    # are used to calculate how the flux in the old bins is
    # distributed to the new bins. The interlinked borders contain
    # both the new and old borders.  union1d returns a sorted array of
    # unique values.
    int_borders = np.union1d(org_borders, new_borders)

    # Rather than only using the intersecting area, use only the
    # output area.  This means, one has to deal with index mismatch
    # later, which also takes care of edge effects.

    int_borders = int_borders[(int_borders >= new_borders[0])
                              & (int_borders <= new_borders[-1])]

    # The centers of the interlinked bins are used to find the indices
    # within the original and new borders.
    int_centers = 0.5 * (int_borders[1:] + int_borders[:-1])
    int_binwidth = int_borders[1:] - int_borders[:-1]
    org_ids = np.digitize(int_centers, org_borders, right=True) - 1
    new_ids = np.digitize(int_centers, new_borders, right=True) - 1

    # All indeices outside the new borders are set to the first and
    # last bin respectively.
    org_ids[org_ids == -1] = 0
    org_ids[org_ids == len(trans)] = len(trans) - 1

    # bincount sums all weight values according to the list of indices
    # in x.  See examples in the documentation of numpy. minlength
    # needs to be there in order to cover bins of the new flux that
    # are larger than the right border of org_borders.
    return np.bincount(new_ids,
                       weights=trans[org_ids] * int_binwidth,
                       minlength=len(new_borders) - 1) / new_binwidth
