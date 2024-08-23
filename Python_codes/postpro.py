import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.stats import gaussian_kde
from matplotlib.animation import FuncAnimation
import run_one_config as roc


def plot_pos_net(pos, opinion, state, CMap, NPop, linkThresh, show_net_bool=True, marker_size=40, bool_from_arr=False, time=0, link_color=(153/255, 255/255, 204/255), fig=None, ax=None, show_Msngrs=False):
    if bool_from_arr:
        pos = pos[:, :, time].reshape(2, -1)
        opinion = opinion[:, time].reshape(-1)
        state = state[:, time].reshape(-1)

    if fig is None:
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.set_axis_off()

    # make color based on normalized opinions and colormap CMap
    norm = plt.Normalize(0, 1)
    colors = CMap(norm(opinion))

    if show_net_bool:
        # Adjc = np.zeros((NPop, NPop))
        Adjc = roc.get_adjacency_matrix(pos, linkThresh)
        for jj in range(NPop - 1):
            for ii in range(jj + 1, NPop):
                if Adjc[jj, ii] == 1:
                    ax.plot([pos[0, jj], pos[0, ii]], [pos[1, jj], pos[1, ii]], color=link_color, linewidth=1, zorder=0)

    sc = ax.scatter(pos[0, state == 1], pos[1, state == 1], marker_size, colors[state == 1], edgecolor='none')
    if show_Msngrs:
        sc2 = ax.scatter(pos[0, state == 0], pos[1, state == 0], 1.5 * marker_size, colors[state == 0], edgecolor='r', linewidth=1.5)
    else:
        sc2 = ax.scatter(pos[0, state == 0], pos[1, state == 0], 1.5 * marker_size, colors[state == 0], edgecolor='none')

    return sc, sc2, fig, ax



def show_animation(result, params=-1):
    posArr = result['posArr']
    opinionArr = result['zpArr']
    stateArr = result['stArr']

    if params == -1:
        NPop = 100
        linkThresh = 0.15
    else:
        NPop = params['NPop']
        linkThresh = params['linkThresh']

    from palettable.scientific.diverging import Berlin_20
    CMap = Berlin_20.mpl_colormap   

    fig = plt.figure(figsize=(7, 6))

    pos_ax1 = [0.1, 0.1, 0.7, 0.7 * 7 / 6]
    pos_ax2 = [0.9, 0.1, 0.07, 0.7 * 7 / 6]

    ax1 = fig.add_axes(pos_ax1)
    ax2 = fig.add_axes(pos_ax2)

    ax1.set_xlim([-1, 1])
    ax1.set_ylim([-1, 1])

    nSkipPlot = 4
    nTrace = 5
    nSkipTrace = 10


    length_of_sim = opinionArr.shape[1] - 1

    for time in range(1, length_of_sim, nSkipPlot):
        ax1.clear()
        ax2.clear()

        # ax1.contour(xx, yy, zz, levels=[zEnv], linestyles='-.', linewidths=1.5, alpha=0.5)
        # ax1.contour(xx, yy, zz, linestyles='--', alpha=0.2)

        # scs_trace, scs_trace2 = [], []
        for i_trace in range(nTrace):
            time_trace = max(1, time - i_trace * nSkipTrace)
            marker_size = (1 - (i_trace - 1) / nTrace) * 40 / 4
            
            plot_pos_net(pos=posArr, opinion=opinionArr, state=stateArr, time=time_trace, CMap=CMap, NPop=NPop, 
                                         linkThresh=linkThresh, show_net_bool=False, marker_size=marker_size, bool_from_arr=True, fig=fig, ax=ax1, show_Msngrs=False)
            # sc, sc2, _, _ = pp.plot_pos_net(pos=posArr, opinion=opinionArr, state=stateArr, time=time_trace, CMap=CMap, NPop=NPop, linkThresh=linkThresh, show_net_bool=True, marker_size=marker_size, bool_from_arr=True)
            # scs_trace.append(sc)
            # scs_trace2.append(sc2)
        
        sc, sc2, _, _ = plot_pos_net(pos=posArr, opinion=opinionArr, state=stateArr, time=time, CMap=CMap, NPop=NPop, linkThresh=linkThresh, 
                                     show_net_bool=True, marker_size=40, bool_from_arr=True, fig=fig, ax=ax1, show_Msngrs=True)

        ax1.set_xlim([-1, 1]); ax1.set_ylim([-1, 1])
        ax1.set_title('Time = ' + str(time)) 

        opinions = opinionArr[:, time]

        # Kernel density estimation for the opinions
        kde = gaussian_kde(opinions)
        densities = kde(opinions)
        density_scale = densities / densities.max()  # Normalize density

        norm = plt.Normalize(0, 1)
        colors = CMap(norm(opinions))

        # Deterministic jittering based on position index
        def deterministic_jitter(index, scale=0.05):
            return (index % 2 - 0.5)/2 * scale + np.random.rand(1)/3 # Jitter in [-0.5*scale, 0.5*scale]

        # Sort opinions for consistent plotting
        sorted_indices = np.argsort(opinions)
        # print(sorted_indices)
        sorted_opinions = opinions[sorted_indices]
        sorted_colors = colors[sorted_indices]
        sorted_density_scale = density_scale[sorted_indices]

        


        # Scatter plot with deterministic jitter
        for i, (opinion, color) in enumerate(zip(sorted_opinions, sorted_colors)):
            plt.scatter(deterministic_jitter(i,sorted_density_scale[i]), opinion, color=color, s=2, alpha=0.6)


        # # Scatter plot with density-based offset
        for opinion, density, color in zip(opinions, density_scale, colors):
            ax2.scatter((density) * np.random.rand(1), opinion, color=color, s=2, alpha=0.6)


        ax2.set_ylim([0, 1.5])
        # set xticks off
        ax2.set_xticks([])

        plt.draw()
        # plt.show()
        plt.pause(0.00001)

    plt.show()


def __main__():
    print("Running one experiment from run_one_config.py", flush=True)
    result, params = roc.__main__()
    print("finished running one experiment")
    

    # # load from pickle file
    # file_name = 'data/single_run_2024_08_15__NPop_100_Arena_1.0__tf_1.0k__landFunc__cone__BasicMarkov__initENumMsngr_sensRang_0.4__i_p2e_14__i_p2m_3.pkl'
    # file_name = 'data/single_run_2024_08_15__NPop_100_Arena_1.0__tf_1.0k__landFunc__cone__BasicMarkov__initENumMsngr_sensRang_0.4__i_p2e_35__i_p2m_44.pkl'
    # result = pickle.load(open(file_name, 'rb'))
    # params = pickle.load(open(file_name[:-4] + '_Params.pkl', 'rb'))

    print("showing animation ...")
    show_animation(result, params=params)

    return 1


if __name__ == '__main__':
    __main__()




