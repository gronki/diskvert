# coding: utf-8



# def rxdtype(info):
#
#     names = [ 'n','z','zscale','rho','T' ]
#     formats = ['i4'] + 4 * ['f4']
#
#     def apcl(l, f = 'f4'):
#         names.append(l)
#         formats.append(f)
#
#     if info.has_corona:
#         apcl('Trad')
#
#     apcl('Frad')
#     apcl('kabs')
#     apcl('ksct')
#
#     if info.has_magnetic:
#         apcl('Pmag')
#         apcl('beta')
#
#     if info.has_conduction:
#         apcl('Fcond')
#         apcl('kcnd')
#
#     return dict(names = names, formats = formats)

def rxread(name = 'disk', frames = None):

    from pyminiconf import pyminiconf
    from numpy import loadtxt,unique
    from col2python import read_dtype

    with open('{n}.txt'.format(n = name),'r') as f:
        info = pyminiconf(f)

    with open('{n}.col'.format(n = name),'r') as f:
        dtype = read_dtype(f)

    data = {}

    framenrs = range(info.niter)
    if frames: framenrs = filter(lambda i: 0 <= i and i <= info.niter, frames)

    for i in unique(framenrs):
        data[i] = loadtxt('{n}.{i:03}.dat'.format(n = name, i = i),
            dtype = dtype,
        )

    return data, info


def axprep(plt,info):
    from matplotlib import rc
    rc('font', family = 'serif', size = 10)

    if info.has_magnetic:
        fig, axes0 = plt.subplots(2, 2, figsize = (8,8))
        axes = axes0.reshape(1,4)[0,:]
    else:
        fig, axes = plt.subplots(1, 3, figsize = (12,4))

    for ax in axes:
        pass
        # ax.set_xscale('log')
        # ax.set_xlabel('tau')
        # ax.set_xlim(1e4,1e-4)
        # ax.axvline(2.0/3.0, color = '#C7E3E6', linewidth = 0.8)

    ax = axes[0]
    ax.set_title('rho')
    ax.set_yscale('log')

    ax = axes[1]
    ax.set_title('T')
    ax.set_yscale('log')

    ax = axes[2]
    ax.set_title('Frad')

    if info.has_magnetic:
        ax = axes[3]
        ax.set_title('beta')
        ax.set_yscale('log')
        ax.set_ylim(1e3,1e-3)

    return fig, axes

def axplot(axes, data, info, **kwargs):
    xkey = 'h'
    # desnity
    axes[0].plot(data[xkey], data['rho'], **kwargs)
    # temperature
    if info.has_corona:
        axes[1].plot(data[xkey], data['trad'], linewidth = 0.8, **kwargs)
    axes[1].plot(data[xkey], data['temp'], **kwargs)
    # fluxes
    axes[2].plot(data[xkey], data['frad'], **kwargs)
    # magnetic field
    if info.has_magnetic:
        axes[3].plot(data[xkey], data['beta'], **kwargs)

def parseargv():
    from sys import argv
    name = argv[1]
    frames = [ int(arg) for arg in argv[2:] ] if len(argv) > 2 else None
    return name,frames

def main_plotmulti():
    import matplotlib.pyplot as plt
    name, frames = parseargv()
    data,info = rxread(name, frames)
    for i,d in data.items():
        fig, axes = axprep(plt, info)
        axplot(axes, d, info)
        plt.savefig('{0}.{1:03}.png'.format(name,i), dpi = 144)
        plt.close(fig)

def main_plotcumul():
    import matplotlib.pyplot as plt
    from matplotlib.cm import get_cmap
    name, frames = parseargv()
    cmap = get_cmap('inferno_r')
    data,info = rxread(name, frames)
    fig, axes = axprep(plt, info)
    for i,d in data.items():
        c = cmap(1.0 * (i + 1) / len(data))
        axplot(axes, d, info, color = c)
    plt.savefig('{0}.multi.png'.format(name), dpi = 144)
    plt.close(fig)

if __name__ == '__main__':
    main_plotmulti()
