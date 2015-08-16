import numpy as np
from matplotlib import pyplot as plt
from mitgcmdata import MITgcmmodel
from mitgcmdata import layers

def process_layers(ddir, iter0, diters, N, deltaTclock=900, layers_name='1RHO',
                  layers_trend_Hc_index=None):

    iters = range(iter0,iter0+N*diters,diters)
    DT = diters * deltaTclock
    m = MITgcmmodel.ModelInstance(ddir, default_iter=list(iters))


    layers_transport = m.rdmds('DiagLAYERS-transport', iters[1:])
    layers_diapycnal = m.rdmds('DiagLAYERS-diapycnal', iters[1:])
    if N > 2:
        layers_transport = layers_transport.mean(axis=0)
        layers_diapycnal = layers_diapycnal.mean(axis=0)
    
    
    uh, vh = layers_transport[0], layers_transport[3]
    print 'uh vh shape =', uh.shape, vh.shape
    pi_s = layers_transport[5]
    layers_trend = m.rdmds('DiagLAYERS-trend', iters)
    
    if layers_trend_Hc_index is not None:
        layers_trend = layers_trend[:, layers_trend_Hc_index].squeeze()
    print 'layers_trend.shape', layers_trend.shape
    
    # process data
    la = layers.LayersAnalyzer(m, layers_name=layers_name)
    dvol_dt = la.volume_trend(layers_trend, DT)
    if dvol_dt.ndim == 4:
        # average over time
        dvol_dt = dvol_dt.mean(axis=0)

    # everything but the advective tendencies
    diapycnal_vel = la.diapycnal_velocity(layers_diapycnal[np.r_[:3,5:8]])

    advective_vel = la.advective_flux_divergence(uh, vh)
    
    # "direct" form of advective vel
    advective_vel_direct = la.diapycnal_velocity(layers_diapycnal[np.r_[3,4,8,9]]).sum(axis=0)
    
    # "TOTTEND"
    diapycnal_vel_tottend = -la.diapycnal_velocity(layers_diapycnal[np.r_[10,11]]/(24*60*60)).sum(axis=0)

    print dvol_dt.shape, diapycnal_vel.sum(axis=0).shape, advective_vel.shape
    diapycnal_vel_numerical = dvol_dt - diapycnal_vel.sum(axis=0) - advective_vel
    diapycnal_vel_numerical_direct = dvol_dt - diapycnal_vel.sum(axis=0) - advective_vel_direct
    diapycnal_vel_numerical_tottend = dvol_dt - diapycnal_vel_tottend
    
    psi = la.calculate_moc(vh) / 1e6
    print psi.shape

    # something is wrong with this
    pi = pi_s.mean(axis=-1)
    
    return {'dvol_dt': dvol_dt,
            'diapycnal_vel': diapycnal_vel,
            'advective_vel': advective_vel,
            'advective_vel_direct': advective_vel_direct,
            'diapycnal_vel_numerical': diapycnal_vel_numerical,
            'diapycnal_vel_numerical_direct': diapycnal_vel_numerical_direct,
            'diapycnal_vel_tottend': diapycnal_vel_tottend,
            'diapycnal_vel_numerical_tottend': diapycnal_vel_numerical_tottend,            
            'psi': psi, 'pi': pi, 'la': la, 'm': m}

def plot_diapycnal_velocities(d, r=0.7e6, ylim = [28.2, 22], direct=False):
    
    if direct:
        diapyc_suff = '_direct'
    else:
        diapyc_suff = ''

    m = d['m']
    Yg = np.hstack([m.yg[0,:,0], 2*m.yg[0,-1,0]-m.yg[0,-2,0]])
    Xg = np.hstack([m.xg[0,0,:], 2*m.xg[0,0,-1]-m.xg[0,0,-2]])
    Yc = m.yc[0,:,0]
    Xc = m.xc[0,0,:]

    lb = d['la'].layers_bounds_w
    lbc = d['la'].layers_top

    plt.figure(figsize=(16,12))

    plt.subplot(321)
    plt.pcolormesh(Yg, lb, d['diapycnal_vel'].sum(axis=0).sum(axis=-1), cmap='RdBu_r')
    plt.clim([-r, r])
    plt.colorbar()
    print 'Yg.shape', Yg.shape
    print "d['psi'].shape", d['psi'].shape
    plt.contour(Yg[:-1], lbc, d['psi'], 10, colors='k')
    #plt.contour(Yg[:-1], lbc, pi, [0.05, 0.5, 0.95], colors='y')
    plt.ylim(ylim)
    plt.title('Explicit Diapycnal Velocity')

    plt.subplot(322)
    plt.pcolormesh(Yg, lb, d['dvol_dt'].sum(axis=-1), cmap='RdBu_r')
    plt.clim([-r, r])
    plt.colorbar()
    plt.contour(Yc, lbc, d['psi'], 10, colors='k')
    plt.ylim(ylim)
    plt.title('Volume Trend')

    plt.subplot(323)
    plt.pcolormesh(Yg, lb, d['advective_vel' + diapyc_suff].sum(axis=-1), cmap='RdBu_r')
    plt.clim([-r, r])
    plt.colorbar()
    plt.contour(Yc, lbc, d['psi'], 10, colors='k')
    plt.ylim(ylim)
    plt.title('Advective Diapycnal Velocity')

    plt.subplot(324)
    plt.pcolormesh(Yg, lb, d['diapycnal_vel_numerical' + diapyc_suff].sum(axis=-1), cmap='RdBu_r')
    plt.clim([-r, r])
    plt.colorbar()
    plt.contour(Yc, lbc, d['psi'], 10, colors='k')
    plt.ylim(ylim)
    plt.title('Numerical Diapycnal Velocity')
    
    plt.subplot(325)
    plt.pcolormesh(Yg, lb, d['diapycnal_vel_tottend'].sum(axis=-1), cmap='RdBu_r')
    plt.clim([-r, r])
    plt.colorbar()
    plt.contour(Yc, lbc, d['psi'], 10, colors='k')
    plt.ylim(ylim)
    plt.title('TOTTEND Diapycnal Velocity')

    plt.subplot(326)
    plt.pcolormesh(Yg, lb, d['diapycnal_vel_numerical_tottend'].sum(axis=-1), cmap='RdBu_r')
    plt.clim([-r, r])
    plt.colorbar()
    plt.contour(Yc, lbc, d['psi'], 10, colors='k')
    plt.ylim(ylim)
    plt.title('Numerical Diapycnal Velocity')
