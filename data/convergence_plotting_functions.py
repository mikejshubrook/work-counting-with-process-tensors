import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def load_dynamics(a, tp, sta, dt, p):
    filename  = f'dynamics-files/dyns-a{a}_G10.0_w25.0_e25.0_tp{tp}_sta{sta}_dt{dt}_p{p}.csv'
    try:
        df = pd.read_csv(filename)
        return df
    except:
        return None

### CONVERGENCE PLOTS FOR DYNAMICS

def plot_dynamics_convergence(alpha, tp, dt_list, p_list, sta_list=[0,1]):
    # dictionary to store dataframes with (dt,p) as key
    dfs = {}

    colors = ['blue', 'red']
    linestyles = ['-', '--', '-.', ':']

    for dt in dt_list:
        for p in p_list:
            for sta in sta_list:
                dfs[(dt, p, sta)] = pd.read_csv(f'dynamics-files/dyns-a{alpha}_G10.0_w25.0_e25.0_tp{tp}_sta{sta}_dt{dt}_p{p}.csv')
         
    # Create figure with three subplots
    fig, ax = plt.subplots(3, 2, figsize=(15, 8))  # match WCF size
    sta_markers = ['o', 's']

    # choose integers for slicing (n=0, m=-1 for the whole range)
    n = 0
    m = -1
    fnt=12 # fontsize
    ms = 0 # markersize
    # NO STA
    # Plot magnetization
    for i, dt in enumerate(dt_list):
        for j, p in enumerate(p_list):
            for sta in sta_list:
                ax[0,sta].plot(dfs[(dt,p, sta)]['tlist'][n:m], dfs[(dt,p, sta)]['mag'][n:m],
               label=f'dt={dt}, p={p}',
               color=colors[i],
               linestyle=linestyles[j],
               marker=sta_markers[sta],
               markersize=ms)

    ax[0,0].set_ylabel('Magnetization', fontsize=fnt)
    ax[0,1].set_ylabel('Magnetization', fontsize=fnt)

        # Plot coherence
    for i, dt in enumerate(dt_list):
        for j, p in enumerate(p_list):
            for sta in sta_list:
                ax[1,sta].plot(dfs[(dt,p, sta)]['tlist'][n:m], dfs[(dt,p, sta)]['coh'][n:m],
               label=f'dt={dt}, p={p}',
               color=colors[i],
               linestyle=linestyles[j],
               marker=sta_markers[sta],
               markersize=ms)

    ax[1,0].set_ylabel('Coherence', fontsize=fnt)
    ax[1,1].set_ylabel('Coherence', fontsize=fnt)

    # Plot overlap_target
    for i, dt in enumerate(dt_list):
        for j, p in enumerate(p_list):
            for sta in sta_list:
                ax[2, sta].plot(dfs[(dt,p, sta)]['tlist'][n:m], dfs[(dt,p, sta)]['overlap_target'][n:m],
               label=f'dt={dt}, p={p}',
               color=colors[i],
               linestyle=linestyles[j],
               marker=sta_markers[sta],
               markersize=ms)
    ax[2,0].set_ylabel('Overlap Target', fontsize=fnt)
    ax[2,1].set_ylabel('Overlap Target', fontsize=fnt)

    for axes in ax.flat:
        axes.legend(ncol=2, fontsize=fnt)
        axes.set_xlabel('Time', fontsize=fnt)

    ax[0,0].set_title('RAW')
    ax[0,1].set_title('STA')
    
    plt.suptitle(f'Convergence tests for dynamics: alpha={alpha}, tp={tp}', fontsize=16)
    plt.tight_layout()
    plt.show()
    return

def load_moments(a, tp, sta, dt, p, s):
    filename  = f'moments-files/moms-a{a}_G10.0_w25.0_e25.0_tp{tp}_sta{sta}_dt{dt}_p{p}_eq{s}.csv'
    df = pd.read_csv(filename)
    return df

### CONVERGENCE PLOTS FOR MOMENTS

def plot_convergence_moments(alpha, tp, dt_list, p_list, s, sta_list=[0,1]):
    # dictionary to store dataframes with (dt,p) as key
    dfs = {}

    colors = ['blue', 'red']
    linestyles = ['-', '--', '-.', ':']

    for dt in dt_list:
        for p in p_list:
            for sta in sta_list:
                dfs[(dt, p, sta)] = load_moments(alpha, tp, sta, dt, p, s)
            
    # Create figure with three subplots
    # fig, ax = plt.subplots(2, 2, figsize=(10, 6))

    # choose integers for slicing (n=0, m=-1 for the whole range)
    n = 0
    m = -1
    fnt=12 # fontsize
    ms = 0 # markersize

    # Get the length of moments data...
    n = len(dfs[(dt, p, sta)]['mean_real'][1::2])
    # ...and create a list of even integers: 2, 4, ..., 2n
    orders = [2 * i for i in range(1, n + 1)]

    # Create figure with three subplots
    fig, ax = plt.subplots(3, 2, figsize=(15, 8))

    # Mean
    for i, dt in enumerate(dt_list):
        for j, p in enumerate(p_list):
            for sta in sta_list:
                ax[0, sta].plot(orders, dfs[(dt,p, sta)]['mean_real'][1::2],
                        '*',
                        label=f'dt={dt}, p={p}',
                        color=colors[i],
                        linestyle=linestyles[j])
    ax[0,0].set_ylabel('Mean', fontsize=fnt)
    ax[0,1].set_ylabel('Mean', fontsize=fnt)          
        # variance
    for i, dt in enumerate(dt_list):
        for j, p in enumerate(p_list):
            for sta in sta_list:
                ax[1, sta].plot(orders, dfs[(dt,p, sta)]['variance_real'][1::2],
                        '*',
                        label=f'dt={dt}, p={p}',
                        color=colors[i],
                        linestyle=linestyles[j])
    ax[1,0].set_ylabel('Variance', fontsize=fnt)
    ax[1,1].set_ylabel('Variance', fontsize=fnt)            
        # Mean
    for i, dt in enumerate(dt_list):
        for j, p in enumerate(p_list):
            for sta in sta_list:
                ax[2, sta].plot(orders, dfs[(dt,p, sta)]['skewness_real'][1::2],
                        '*',
                        label=f'dt={dt}, p={p}',
                        color=colors[i],
                        linestyle=linestyles[j])
    ax[2,0].set_ylabel('Skewness', fontsize=fnt)
    ax[2,1].set_ylabel('Skewness', fontsize=fnt)            
    
    ax[0,0].set_title('RAW')
    ax[0,1].set_title('STA')
    
    for axes in ax.flat:
        axes.set_xlabel('Finite Difference Order', fontsize=fnt)        
        axes.legend(ncol=len(dt_list), fontsize=fnt)

    plt.suptitle(f'Convergence tests for moments: alpha={alpha}, tp={tp}', fontsize=16)
    plt.tight_layout()
    plt.show()
    return

def load_wcf(a, tp, sta, dt, p, s, X0, Xf):
    filename  = f'wcf-files-combined/WCF_a{a}_G10.0_w25.0_e25.0_t{tp}_sta{sta}_dt{dt}_p{p}_eq{s}_X{X0}_Xf{Xf}.txt'
    wcf = pd.read_csv(filename, sep='\t', header=None)
    return wcf

### CONVERGENCE PLOTS FOR MOMENTS

def plot_convergence_wcf(alpha, tp, dt_list, p_list, s, X0, Xf, sta_list=[0,1]):
    # dictionary to store dataframes with (dt,p) as key
    wcfs = {}

    colors = ['blue', 'red']
    linestyles = ['-', '--', '-.', ':']

    for dt in dt_list:
        for p in p_list:
            for sta in sta_list:
                wcfs[(dt, p, sta)] = load_wcf(alpha, tp, sta, dt, p, s, X0, Xf)
            
    # Create figure with three subplots
    # fig, ax = plt.subplots(2, 2, figsize=(10, 6))

    # choose integers for slicing (n=0, m=-1 for the whole range)
    n = 0
    m = -1
    fnt=12
    sta_markers = ['o', 's']
    chis = {}
    for dt in dt_list:
        chis[dt] = np.arange(0, len(wcfs[(dt, p, sta)][0]))*dt

    # Create figure with three subplots
    fig, ax = plt.subplots(2, 2, figsize=(15, 8))

    # Mean
    for i, dt in enumerate(dt_list):
        for j, p in enumerate(p_list):
            for sta in sta_list:
                ax[0, sta].plot(chis[dt], wcfs[(dt,p, sta)][0],
                        marker=sta_markers[sta],
                        label=f'dt={dt}, p={p}',
                        color=colors[i],
                        linestyle=linestyles[j])
    ax[0,0].set_ylabel('Real WCF', fontsize=fnt)
    ax[0,1].set_ylabel('Real WCF', fontsize=fnt)
        # variance
    for i, dt in enumerate(dt_list):
        for j, p in enumerate(p_list):
            for sta in sta_list:
                ax[1, sta].plot(chis[dt], wcfs[(dt,p, sta)][1],
                        marker=sta_markers[sta],
                        label=f'dt={dt}, p={p}',
                        color=colors[i],
                        linestyle=linestyles[j])
    ax[1,0].set_ylabel('Imag  WCF', fontsize=fnt)
    ax[1,1].set_ylabel('Imag  WCF', fontsize=fnt)
    ax[0,0].set_title('RAW')
    ax[0,1].set_title('STA')
    
    for axes in ax.flat:
        axes.set_xlabel('Counting parameter', fontsize=fnt)        
        axes.legend(ncol=len(dt_list), fontsize=fnt)

    plt.suptitle(f'Convergence tests for WCF: alpha={alpha}, tp={tp}', fontsize=16)
    plt.tight_layout()
    plt.show()
    return