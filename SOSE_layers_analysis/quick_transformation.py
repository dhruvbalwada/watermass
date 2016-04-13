import numpy as np

def surface_transformation_rate(rac, rho1, rho1levs,
                                rho2=None, rho2levs=None,
                                fields=[]):
    """Calculate surface transformation rate given
        rac: area weights
        rho: a surface density field
        rholevs: density levels (to define bins)
        args: various components of the density flux
    """

    drho1 = np.diff(rho1levs)
    Nbins1 = len(rho1levs)-1
    
    # mask anything outside the range
    rho1_idx, mask1 = get_rho_idx(rho1, rho1levs)
    assert rho1_idx.min()>=0
    assert rho1_idx.max()<Nbins1
        
    if rho2 is None:
        rho_idx = rho1_idx
        mask = mask1
        drho = drho1
        Nbins = Nbins1
    else:
        rho2_idx, mask2 = get_rho_idx(rho2, rho2levs)
        drho2 = np.diff(rho2levs)
        Nbins2 = len(rho2levs)-1
        assert rho2_idx.min()>=0
        assert rho2_idx.max()<Nbins2
        
        # turn 2d coords into 1d indexer
        # rho1_idx is the inner axis
        rho_idx = Nbins1*rho2_idx + rho1_idx
        mask = mask1 | mask2
        drho = np.outer(drho2, drho1).ravel()
        Nbins = Nbins1 * Nbins2
    
    res = []
    for a in fields:
        dens_flux = np.ma.masked_array(a * rac, mask)
        q = np.bincount(rho_idx.compressed(), weights=dens_flux.compressed(), minlength=Nbins) / drho
        if rho2 is not None:
            q = q.reshape([Nbins2, Nbins1])
        res.append(q)
    return np.array(res)
    
def get_rho_idx(rho, rholevs):
    rho_m = np.ma.masked_greater_equal(
            np.ma.masked_less(rho, rholevs.min()), rholevs.max())
    mask = rho_m.mask
    # don't want to compress because we want to preserve the size of the array
    #rho_idx = np.digitize(rho_m.compressed(), rholevs)-1
    Nrho = len(rholevs)
    rho_idx = np.ma.masked_outside(np.digitize(rho, rholevs)-1, 0, Nrho-2)
    return rho_idx, rho_idx.mask

def idx_2d_to_1d(rho1_idx, Nbins1, rho2_idx, Nbins2):
    """Convert a pair of indices to a single index"""
    
    
def average_field_in_rho(rac, rho, rholevs, a, avg=True):
    """Calculate surface transformation rate given
        rac: area weights
        rho: a surface density field
        rholevs: density levels (to define bins)
        args: various components of the density flux
    """

    drho = np.diff(rholevs)
    Nbins = len(rholevs)-1
    
    # mask anything outside the range
    rho_m = np.ma.masked_greater_equal(
                np.ma.masked_less(rho, rholevs.min()), rholevs.max())
    mask = rho_m.mask
    
    rho_idx = np.digitize(rho_m.compressed(), rholevs)-1
   
    # rho_idx==i means rholevs[i-1] <= rholevs < rholevs[i]
    
        
    rac_masked = np.ma.masked_array(rac, mask)
    area = np.bincount(rho_idx, weights=rac_masked.compressed(), minlength=Nbins)
    area_recip = np.where(area>0., area**-1, 0.)

    # area-weighted average
    res = np.bincount(rho_idx, weights=(a*rac_masked).compressed(), minlength=Nbins)
    if avg:
        res *= area_recip
    return res, area