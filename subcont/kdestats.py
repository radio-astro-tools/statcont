import numpy as np
from scipy.stats import gaussian_kde


def kde_max(kde, tol=1e-3):
    """Determine the maximum value of a KDE
    INPUT: scipy.stats.kde object"""
    # 2010-07-23 11:14 IJC: Created
    cenguess = np.median(kde.dataset,1)
    init_steps = np.std(kde.dataset,1)
    max_vals = np.zeros(kde.d)
    for dim in range(kde.d):
        this_step = init_steps[dim].copy()
        tempkde = gaussian_kde(kde.dataset[dim,:])
        this_cen = cenguess[dim].copy()
        dguess = 1.
        while dguess>tol:
            pos_guess = this_cen + this_step
            neg_guess = this_cen - this_step
            if tempkde(pos_guess)>tempkde(this_cen):
                dguess = abs(this_cen-pos_guess)
                this_cen = pos_guess
            elif tempkde(neg_guess)>tempkde(this_cen):
                dguess = abs(this_cen-neg_guess)
                this_cen = neg_guess
            else: 
                this_step *= 0.5
                dguess *=0.5
        max_vals[dim] = this_cen
    return max_vals

def findkdeval(kde, val, guess=None,tol=1e-6, maxiter=1000,verbose=False):
    """Find the "x"-value such that kde(x)=val +/- tol

    Without a specified guess for a monomodal distribution, tends to
    find the lower of the two possible values.

    Uses scipy.stats.gaussian_kde objects"""
    # 2010-07-23 11:49 IJC: Created
    if maxiter is None:
        maxiter = 1000
    
    if kde.d>1:
        retvals = np.zeros(kde.d,float)
        if hasattr(val,'__iter__'):
            pass
        else:
            val = np.tile(val,kde.d)

        for dim in range(kde.d):
            print dim
            tempkde = gaussian_kde(kde.dataset[dim,:])
            retvals[dim] = findkdeval(tempkde, val[dim],tol=tol, guess=guess,maxiter=maxiter)
    else:
        max_location = kde_max(kde)
        maxval = kde(max_location)
        if val>maxval:
            print "Value entered (%01.3g) is greater than KDE's maximum value (%01.3g)" % (val,maxval)
            return max_location
        
        iter = 0
        initial_step = np.std(kde.dataset)
        step = initial_step
        medval = np.median(kde.dataset)
        if guess is None:
            guess = medval - step
        else:
            pass

        initial_guess = guess # Save this, just in case
        kdeval = kde(guess)
        dval = abs(kdeval-val)
        while dval>tol and iter<=maxiter:
#            pos = min(guess + step, max_location)
            pos = guess + step
            neg = guess - step
            kdepos = kde(pos)
            kdeneg = kde(neg)
            if verbose and (iter/10)==(iter/10.):
                print "iter, val, pos, neg, kdepos,kdeneg,guess, step"
            if verbose:
                print "%i %01.6g %01.6g %01.6g %01.6g %01.6g %01.6g" %(iter, val, pos, neg,kdepos,kdeneg, step)
            dval_pos = abs(kdepos-val)
            dval_neg = abs(kdeneg-val)
            if dval_pos<dval: # we're closer than we were!
                dval = dval_pos
                guess = pos
            elif dval_neg<dval: # we're closer than we were!
                dval = dval_neg
                guess = neg
            elif dval_pos==dval_neg: # We must be really far away
                if initial_guess>medval:
                    guess = medval + initial_step
                else:
                    guess = medval - initial_step
            else: # neither guess is closer, but we're still not within tol
                step *= 0.5

            iter += 1
        retvals = guess
        if iter>maxiter:
            print "Exceeded maximum number of iterations (%i) without convergance" % maxiter

    return retvals
                
            
 

def conflevel(kde, frac, ftol=1e-6, tol=1e-6, usespline=False, verbose=False,maxiter=None):
    """Determine the lower and upper confidence levels required to
    enclose a given fraction 'frac' of a KDE object's dataset(s) to
    within a tolerance ftol."""
    # 2010-07-23 14:00 IJC: Created
    from scipy import interpolate

    if kde.d>1:
        ret = []
        for ii in range(kde.d):
            if verbose: print ii,'/',kde.d
            ret.append(conflevel(gaussian_kde(kde.dataset[ii]),frac,ftol=ftol,tol=tol,usespline=usespline,verbose=verbose,maxiter=maxiter))
        return ret
    else:
        median_value = np.median(kde.dataset)
        step = np.std(kde.dataset)
        low_limit = median_value - 2*step
        enclosed_fraction = np.Inf
        thismovewasup=True
        if usespline:
            spx = np.linspace(-5*step,5*step,1e3)+median_value
            sp = interpolate.UnivariateSpline(spx,kde(spx),k=3.0,s=0.0)
            sp.d = kde.d
            sp.dataset = kde.dataset
            kde0 = kde
            kde = sp

        while abs(enclosed_fraction-frac)>ftol:
            kdelow = kde(low_limit)
            guess = 2*median_value-low_limit
            if verbose:
                print "kdelow,guess",kdelow,guess
            high_limit = findkdeval(kde, kdelow, tol=tol, guess=guess,verbose=verbose,maxiter=maxiter)
            if usespline:
                enclosed_fraction = kde0.integrate_box_1d(low_limit, high_limit)
            else:
                enclosed_fraction = kde.integrate_box_1d(low_limit, high_limit)

            lastmovewasup = thismovewasup
            if enclosed_fraction > frac: # need to close in
                low_limit += step
                thismovewasup = True
            else: # need to open up limits
                low_limit += -step
                thismovewasup = False

            movingInSameDirection = lastmovewasup==thismovewasup
            if movingInSameDirection: # keep going
                pass
            else: # We passed it; turn around and home in
                step *= 0.5

    return low_limit, high_limit


def confmap(map, frac, **kw):
    """Return the confidence level of a 2D histogram or array that
    encloses the specified fraction of the total sum.

    :INPUTS:
      map : 1D or 2D numpy array
        Probability map (from hist2d or kde)

      frac : float, 0 <= frac <= 1
        desired fraction of enclosed energy of map

    :OPTIONS:
      ordinate : None or 1D array
        If 1D map, interpolates onto the desired value.  This could
        cause problems when you aren't just setting upper/lower
        limits....
    """
    # 2010-07-26 12:54 IJC: Created
    # 2011-11-05 14:29 IJMC: Fixed so it actually does what it's supposed to!
    from scipy.optimize import bisect

    def diffsum(level, map, ndesired):
        return ((1.0*map[map >= level].sum()/map.sum() - ndesired))

    if hasattr(frac,'__iter__'):
        return [confmap(map,thisfrac, **kw) for thisfrac in frac]

    #nx, ny = map.shape
    #ntot = map.size
    #n = int(ntot*frac)

    #guess = map.max()
    #dx = 10.*float((guess-map.min())/ntot)
    #thisn = map[map<=guess].sum()

    ret = bisect(diffsum, map.min(), map.max(), args=(map, frac))
    if kw.has_key('ordinate') and kw['ordinate'] is not None:
        sortind = np.argsort(map)
        ret = np.interp(ret, map[sortind], kw['ordinate'][sortind])

    return ret


def kdehist2(x, y, npts, xrange=None, yrange=None):
    """Generate a 2D histogram map from data, using Gaussian KDEs

    :INPUTS:
      x : seq
        X data

      y : seq
        Y data

      npts : int or 2-seq
        number of points across final histogram, or [nx, ny]

    :OPTIONAL_INPUTS:
      xrange : 2-seq
        [x_min, x_max] values for final histogram

      yrange : 2-seq
        [y_min, y_max] values for final histogram
        
    :RETURNS:
      [kdehist, xbins, ybins]

    :EXAMPLE:
      ::
 
        import kdestats as kde
        import numpy as np
        import pylab as py

        covmat = [[1., 1.5], [1.5, 4.]]
        xy = np.random.multivariate_normal([0, 0], covmat, [1e4])
        kdehist = kde.kdehist2(xy[:,0], xy[:,1], [30, 30])
        clevels = kde.confmap(kdehist[0], [.6827,.9545,.9973])

        py.figure()  # Plot 1-, 2-, and 3-sigma contours
        c = py.contour(kdehist[1], kdehist[2], kdehist[0], clevels)
    """
    # 2012-02-11 19:45 IJMC: Created

    if hasattr(npts, '__iter__'):
        if len(npts)==1:
            npts = [npts[0], npts[1]]
    else:
        npts = [npts, npts]

    # Generate KDE:
    thiskde = gaussian_kde([x, y])

    # Generate coordinates for KDE evaluation:
    if xrange is None:
        thisx0 = np.median(thiskde.dataset[0])
        thisdx0 = np.std(thiskde.dataset[0])
        thisx = np.linspace(-5*thisdx0,5*thisdx0,npts[0])+thisx0
    else:
        thisx = np.linspace(xrange[0], xrange[1], npts[0])

    if yrange is None:
        thisy0 = np.median(thiskde.dataset[1])
        thisdy0 = np.std(thiskde.dataset[1])
        thisy = np.linspace(-5*thisdy0,5*thisdy0,npts[1])+thisy0
    else:
        thisy = np.linspace(yrange[0], yrange[1], npts[1])

    thisxx,thisyy = np.meshgrid(thisx,thisy)
    thishist = thiskde([thisxx.ravel(),thisyy.ravel()]).reshape(npts[0],npts[0])

    return thishist, thisx, thisy
