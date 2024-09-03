"""
Title: Simulation Runner Functions and Script

Description: This script contains the functions needed to run simulations of the opinion dynamics model with Messengers using DMP. 
The "main" script runs an independent simulation with predetermined parameters. Then it shows an animation of the simulation run, using post-processing codes. 
The animation shows the evolution of the opinion dynamics and the movement of the agents in the environment. 
Depending on the parameters, you can observe how messengers break echo chambers in collective opinion dynamics with homophily.
You can change the parameters in the "main" function to run different simulations.

Author: Mohsen Raoufi

Contact: mohsenraoufi@icloud.com

Affiliation: Research Cluster of Excellence, "Science of Intelligence"

Date Created: August, 2024

Version: 1.0

Usage: Run this script with Python 3.8 or later. Ensure you create a "data" directory in the working directory before execution.

License: Distributed under the MIT License. See LICENSE.txt for more information.

Citation: Please cite our work if you use this code:
"Messengers: Breaking Echo Chambers in Collective Opinion Dynamics with Homophily"
"""

import numpy as np
from datetime import date
import pickle
import postpro as pp
import os

# added chmod
def get_adjacency_matrix(pos, linkThresh):
    """ Update the adjacency matrix based on the positions of agents and the link threshold.
    args:
        pos: Positions of agents
        linkThresh: Link threshold for the spatial network
    returns:
        The updated adjacency matrix """
    
    NPop = pos.shape[1]
    Adjc = np.zeros((NPop, NPop))
    for jj in range(NPop-1):
        for ii in range(jj+1, NPop):
            if np.linalg.norm(pos[:, jj] - pos[:, ii]) < linkThresh:
                Adjc[ii, jj] = 1
    Adjc = Adjc + Adjc.T
    return Adjc


def func_run_experiment(landFunc, LB, initPos, linkThresh=0.5, tf=1000, nTVars=1000, alpha=0.99, beta=0.9, sigma=0.001, p2explt=1.0, p2msngr=0.0, 
                        randStep_exp=1.0, randStep_msg=0.0001, stepSize=0.002, init_messengers=True, randSeed=-1, bool_use_chain_gradient=True):
    """ Run the experiment with the given parameters and return the results. 
    The function is a Python implementation of the funcEEM_Markov_new function in the MATLAB script.
    args:
        landFunc: The landscape function
        LB: Lower and upper bounds for x and y
        pos: Initial positions of agents
        linkThresh: Link threshold for the spatial network, i.e., the maximum distance between two agents for them to be connected
        tf: Total time frames
        nTVars: Number of time variables for saving data points
        alpha: Alpha parameter: weight on opinion updating
        beta: Beta parameter: weight on pseudo-gradient descending
        sigma: Sigma for noise
        p2explt: Probability to become an exploiter
        p2msngr: Probability to become a messenger
        randStep_exp: Random step for exploiters
        randStep_msg: Random step for messengers
        stepSize: Step size for movement of agents in space
        init_messengers: Flag to initialize messengers based on the given probabilities p2explt and p2msngr, otherwise initialize all as exploiters,
        randSeed: Random seed for reproducibility, if -1, no seed is set
        bool_use_chain_gradient: Flag to use the approximation of the chain rule in pseudo-gradient descend
        
    returns:
        A dictionary containing the arrays of interest (posArr, zsArr, zpArr, etc.) """

    def calculate_walk(Walk, ZS, local_consensus, grad, rndWalkU, randStep, swch, stepSize, time, bool_use_chain_gradient=True):
        """ Implementation of the pseudo-gradient descend: Make a step based on the approximation of gradient and a component of random walk 
        args:
            Walk: Current walk of agents
            ZS: Current noisy measurement of agents
            local_consensus: Neighbors' opinion
            grad: approximation of gradient (based on the difference)
            rndWalkU: Random walk component
            randStep: Random step size
            swch: Switching array
            stepSize: Step size for movement of agents in space
            time: Current time frame
            bool_use_chain_gradient: Flag to use the approximation of the chain rule in pseudo-gradient descend
        returns:
            The updated walk of agents """
        
        if time > 2:
            if(bool_use_chain_gradient):
                # NOTE: here we use (ZS-local_consensus) and multiply it by the gradient as an approximation of the dissonance gradient (ZS-local_consensus)^2 
                # assuming the local_consensus (LC) is constant, then we can say that the gradient G = (ZS-LC)^2~2*(ZS-LC)*d(ZS-LC)/dZS = 2*(ZS-LC)*grad
                Walk = (1-randStep) * (ZS - local_consensus) * (np.sum(Adjc, axis=0) / (1 + np.sum(Adjc, axis=0))) * grad
            else:
                # use this for the actual pseudo-gradient descend, without the approximation of the chain rule
                Walk = (1-randStep) * grad

            Walk += randStep * rndWalkU
            Walk = -stepSize * normc(Walk)
        else:
            Walk = stepSize * rndWalkU
        Walk[:, swch == 0] = stepSize * rndWalkU[:, swch == 0]
        return Walk

    def normc(matrix):
        return matrix / np.linalg.norm(matrix, axis=0, keepdims=True)

    #### Initialize variables
    NPop = initPos.shape[1]

    ### make the initial number of exploiters and messengers
    ## swch shows the state of the agents (0: messenger, 1: exploiter)
    if init_messengers:
        # calculate the expected ratio of exploiters
        ExpNExploiter = p2explt / (p2explt + p2msngr)
        # assign random state to agents based on the expected ratio
        swch = np.random.rand(NPop) < ExpNExploiter
    else:
        swch = np.zeros(NPop).astype(bool)

    randStep = randStep_exp * np.ones(NPop)
    pos = initPos.copy()

    ## calculate the (initial) adjacency matrix based on pos and linkThresh
    Adjc = get_adjacency_matrix(initPos, linkThresh)

    # calculate the initial reading from the landscape function
    ZR = landFunc(initPos[0, :], initPos[1, :])
    # add noise to the reading
    ZS = ZR + sigma * np.random.randn(NPop)
    # initialize the opinion of agents
    ZP = ZS.copy()

    ## intialize pseudo-gradient descend variables
    # initialize the random orientation of agents for random walk in pseudo-gradient descend
    rndOrient = np.random.uniform(-np.pi, np.pi, NPop)
    Walk = 0.0001 * np.random.uniform(-1, 1, (2, NPop))
    grad = np.random.uniform(-1, 1, (2, NPop))
    prevgrad = grad.copy()
    objFArrOld = 0

    ## initialize the arrays for saving the results
    nSkipSave = round(tf / nTVars)
    zpArr = np.full((NPop, nTVars), np.nan)
    zsArr = np.full((NPop, nTVars), np.nan)
    z1Arr = np.full(nTVars, np.nan)
    z1StdArr = np.full(nTVars, np.nan)
    stArr = np.full((NPop, nTVars), np.nan)
    posArr = np.full((2, NPop, nTVars), np.nan)
    debugArr = np.full((NPop, nTVars), np.nan)

    saveCtr = 0
    # set the random seed if it is given
    if randSeed != -1:
        np.random.seed(randSeed)
    else:
        np.random.seed()


    #### loop over the time frames

    for time in range(1, tf + 1):

        prevSwch = swch.copy()

        # receive the (noiseless) reading from the landscape function
        ZR = landFunc(pos[0, :], pos[1, :])
        # add noise to the reading
        ZS = ZR + sigma * np.random.randn(NPop)

        #### iterate the Dichotomous Markov Process (DMP) and update the state of agents
        # make a random array for each agent and check if it is less than the probability to switch
        randArr = np.random.rand(NPop)
        posib2msngr = p2msngr * swch
        posib2explt = p2explt * (~swch.astype(bool))

        expChange = randArr < posib2msngr
        msgChange = randArr < posib2explt

        swch[expChange] = ~swch[expChange]
        swch[msgChange] = ~swch[msgChange]

        # find where the switch happened and use the appropriate random step size: randStep_exp or randStep_msg
        # if swch == 1 (exploiter), use randStep_exp, otherwise use randStep_msg
        randStep[swch == 0] = randStep_msg
        randStep[swch == 1] = randStep_exp

        # reset the gradient and pseudo-gradient descend variables for the agents that switched (they should not have any memory)
        # if the current state is different than the just one before: justSwitched(i) = True
        justSwitched = swch != prevSwch
        grad[:, justSwitched] = np.random.uniform(-1, 1, (2, np.sum(justSwitched)))
        prevgrad[:, justSwitched] = np.random.uniform(-1, 1, (2, np.sum(justSwitched)))


        ####  update the opinions based on the DeGroot model
        # calculate the cardinality of the neighbors of each agent (node_degree: number of neighbors)
        node_degree = np.sum(Adjc, axis=0)
        local_consensus = (ZS + ZP @ Adjc) / (node_degree + 1)

        # do not update the opinion of messengers
        ZP[swch != 1] = ZP[swch != 1]
        # update the opinion of exploiters based on the DeGroot model
        ZP[swch == 1] = alpha * ZP[swch == 1] + (1 - alpha) * local_consensus[swch == 1]

        #### do the pseudo-gradient descend
        ## calculate the objective function value based on the measurement of the landscape function 

        
        if(bool_use_chain_gradient):
            # NOTE: that the actual objective function is something like (ZP-ZS)^2, but we put this into the walk update function where we used (ZS - neigh)
            objFArr = ZS
        else:
            # this is the actual objective function, when not using the approximation of the chain rule
            objFArr = (ZP - ZS)**2


        # difference of the objective function from the last time frame
        diffObjF = objFArr - objFArrOld
        # update the old objective function value
        objFArrOld = objFArr
        # temporary angle for the pseudo-gradient descend
        tmpAng = np.arctan2(Walk[1, :], Walk[0, :])

        # current gradient of the objective function
        curGrad = diffObjF * np.array([np.cos(tmpAng), np.sin(tmpAng)])
        # smooth the gradient using a time-decaying filter with previous gradient
        grad = beta * prevgrad + (1 - beta) * curGrad
        # update the previous gradient
        prevgrad = curGrad.copy()

        # if there is any NaN in the gradient, set it to zero
        if np.isnan(grad).any():
            grad[np.isnan(grad)] = 0

        
        # update the oriention of messengers with a bit of inertia and random walk: for a bit more balistic movement than a pure random walk
        rot_inertia_msg = 0.1
        rndOrient[swch != 1] += rot_inertia_msg * np.random.uniform(-np.pi, np.pi, np.sum(swch != 1))
        rndOrient[swch == 1] = np.random.uniform(-np.pi, np.pi, np.sum(swch == 1))

        # make a unit vector for the random walk orientation
        rndWalkU = np.array([np.cos(rndOrient), np.sin(rndOrient)])

        # calculate the walk, based on the ZS and the local consensus, the gradient, and the random walk component
        Walk = calculate_walk(Walk=Walk, ZS=ZS, local_consensus=local_consensus, grad=grad, rndWalkU=rndWalkU, randStep=randStep, swch=swch, stepSize=stepSize, time=time, bool_use_chain_gradient=bool_use_chain_gradient)

        # if there is any NaN in the walk, set it to a random walk
        if np.isnan(Walk).any():
            indxNaN = np.isnan(Walk)
            Walk[indxNaN] = stepSize * randStep * rndWalkU[indxNaN]

        # if there is any agent that is not connected to any other agent, set its walk to a random walk
        if np.sum(np.sum(Adjc) == 0):
            indxDetach = np.where(np.sum(Adjc, axis=0) == 0)
            Walk[:, indxDetach] = stepSize * rndWalkU[:, indxDetach]

        # normalize the walk to the step size
        nWalk = np.linalg.norm(Walk, axis=0)
        nWalkN = np.minimum(nWalk, stepSize)
        Walk = (Walk / nWalk) * nWalkN

        # take a step based on the walk
        pos += Walk

        # make sure the agents are within the bounds
        pos[0, :] = np.clip(pos[0, :], LB[0, 0], LB[0, 1])
        pos[1, :] = np.clip(pos[1, :], LB[1, 0], LB[1, 1])

        # if an agent is close to the boundary, push it back and change its orientation: x-axis
        punishMargin = 0.05
        if (np.min(pos[0, :]) < (LB[0, 0] + punishMargin * linkThresh)) or (np.max(pos[0, :]) > (LB[0, 1] - punishMargin * linkThresh)):
            indPassed = np.where(pos[0, :] < (LB[0, 0] + punishMargin * linkThresh))[0]
            pos[0, indPassed] += 2 * np.abs(Walk[0, indPassed])
            rndrnd = np.random.uniform(-np.pi, np.pi, NPop)
            rndOrient[indPassed] = rndrnd[indPassed]

            indPassed = np.where(pos[0, :] > (LB[0, 1] - punishMargin * linkThresh))[0]
            pos[0, indPassed] -= 2 * np.abs(Walk[0, indPassed])
            rndOrient[indPassed] = rndrnd[indPassed]
            Walk[:, indPassed] = stepSize * np.array([np.cos(rndOrient[indPassed]), np.sin(rndOrient[indPassed])])

        # if an agent is close to the boundary, push it back and change its orientation: y-axis
        if (np.min(pos[1, :]) < (LB[1, 0] + punishMargin * linkThresh)) or (np.max(pos[1, :]) > (LB[1, 1] - punishMargin * linkThresh)):
            indPassed = np.where(pos[1, :] < (LB[1, 0] + punishMargin * linkThresh))[0]
            pos[1, indPassed] += 2 * np.abs(Walk[1, indPassed])
            rndrnd = np.random.uniform(-np.pi, np.pi, NPop)
            rndOrient[indPassed] = rndrnd[indPassed]

            indPassed = np.where(pos[1, :] > (LB[1, 1] - punishMargin * linkThresh))[0]
            pos[1, indPassed] -= 2 * np.abs(Walk[1, indPassed])
            rndOrient[indPassed] = rndrnd[indPassed]
            Walk[:, indPassed] = stepSize * np.array([np.cos(rndOrient[indPassed]), np.sin(rndOrient[indPassed])])

        # update the adjacency matrix based on the new positions
        Adjc = get_adjacency_matrix(pos, linkThresh)

        # save the results  (if it is the time to save)
        if time % nSkipSave == 0:
            posArr[:, :, saveCtr] = pos
            zsArr[:, saveCtr] = ZS
            zpArr[:, saveCtr] = ZP
            z1Arr[saveCtr] = np.mean(ZS)
            z1StdArr[saveCtr] = np.std(ZS)
            stArr[:, saveCtr] = swch
            debugArr[:, saveCtr] = np.arctan2(Walk[1, :], Walk[0, :])
            saveCtr += 1
            
    return {
        "posArr": posArr,
        "zsArr": zsArr,
        "zpArr": zpArr,
        "z1Arr": z1Arr,
        "z1StdArr": z1StdArr,
        "stArr": stArr,
        # "debugArr": debugArr
    }

def __main__(params=None):

    if(params is None):
        # index of the parameters for probabilities in DMP (according to the paper)

        # Try the following preset parameters

        ## No Messenger
        # i_p2e = 0
        # i_p2m = 48

        ## Fast Switching Messengers (enhanced exploration) (majority exploiters)
        # i_p2e = 6
        # i_p2m = 7

        ## Fast Switching Messengers (majority messengers) (needs longer simulations)
        # i_p2e = 14
        # i_p2m = 3

        ## Specialized few Messengers
        i_p2e = 35
        i_p2m = 44

        ## All Messengers
        # i_p2e = 48
        # i_p2m = 0

        ### Simulation parameters
        tf = int(10e2)
        nTVars = int(10e2) 
        NPop = 100
        linkThresh = 0.4
        ArenaScale = 1.0

        alpha = 0.99
        beta = 0.9
        sigma = 0.001
        randStep_msg = 1.0
        randStep_exp = 0.25 * 0.0001
        stepSize = 0.002

        ### Landscape function
        # cone
        land_func_name = 'cone'
        landFunc_nonScaled = lambda x, y: np.sqrt(x**2 + y**2)

        # # pyramid
        # land_func_name = 'pyramid'
        # landFunc_nonScaled = lambda x, y: (np.abs(x) + np.abs(y))

        # # ramp x
        # land_func_name = 'ramp_x'
        # landFunc_nonScaled = lambda x, y: x

        # # ramp y
        # land_func_name = 'ramp_y'
        # landFunc_nonScaled = lambda x, y: y

        # # ramp x and y
        # land_func_name = 'ramp_x_y'
        # landFunc_nonScaled = lambda x, y: x + y

        # # sine shape x, y
        # land_func_name = 'sine_x_y'
        # landFunc_nonScaled = lambda x, y: np.sin(4*np.pi*x) + np.sin(4*np.pi*y)

        # # sine shape x, y
        # land_func_name = 'sine2_x_y'
        # landFunc_nonScaled = lambda x, y: np.sin(2*np.pi*x) + np.sin(2*np.pi*y)     

        # make a dictionary of params 
        params = {
            'i_p2e': i_p2e,
            'i_p2m': i_p2m,
            'tf': tf,
            'nTVars': nTVars,
            'NPop': NPop,
            'linkThresh': linkThresh,
            'ArenaScale': ArenaScale,
            'alpha': alpha,
            'beta': beta,
            'sigma': sigma,
            'randStep_msg': randStep_msg,
            'randStep_exp': randStep_exp,
            'stepSize': stepSize,
            'land_func_name': land_func_name,
        } 
           
        
    
    else:
        i_p2e = params['i_p2e']
        i_p2m = params['i_p2m']

        ### Simulation parameters
        tf = params['tf']
        nTVars = params['nTVars']
        NPop = params['NPop']
        linkThresh = params['linkThresh']
        ArenaScale = params['ArenaScale']

        alpha = params['alpha']
        beta = params['beta']
        sigma = params['sigma']
        randStep_msg = params['randStep_msg']
        randStep_exp = params['randStep_exp']
        stepSize = params['stepSize']

        ### Landscape function
        land_func_name = params['land_func_name']
        landFunc_nonScaled = params['landFunc_nonScaled']

    # Land bounds and initialization
    landBounds = ArenaScale * np.array([[-1, 1], [-1, 1]])
    initBounds = landBounds
    nMC = 1

    # make the environment landscape function
    nDiscreteSpace = 200    # number of points in each dimension
    x = np.linspace(landBounds[0, 0], landBounds[0, 1], nDiscreteSpace)
    y = np.linspace(landBounds[1, 0], landBounds[1, 1], nDiscreteSpace)
    xx, yy = np.meshgrid(x, y)
    zz = landFunc_nonScaled(xx, yy) # the landscape function (without scaling)
    # scale the landscape function to [0, 1]
    landFunc = lambda x, y: (landFunc_nonScaled(x, y) - np.min(zz)) / (np.max(zz) - np.min(zz)) 
    zz = landFunc(xx, yy)
    # the ground-truth mean of the landscape function
    zEnv = np.mean(zz)

    # parameter space we used for paper
    p2expltArr = np.exp(-np.arange(1, 101, 2) / 5)
    p2msngrArr = p2expltArr[:-1]

    # Single parameter check
    p2explt = p2expltArr[i_p2e]
    p2msngr = p2msngrArr[i_p2m]


    # Filename for saving results
    dateStr = date.today().strftime("%Y_%m_%d")
    # Get the file directory and use it to save the results
    cur_dir = os.path.dirname(os.path.realpath(__file__))

    mainStrng = f"data/single_run_{dateStr}__NPop_{NPop}_Arena_{ArenaScale}__tf_{tf/1000}k__" \
                f"landFunc__{land_func_name}__BasicMarkov__initENumMsngr_sensRang_{linkThresh}__" \
                f"i_p2e_{i_p2e}__i_p2m_{i_p2m}.mat"
    
    mainStrng = os.path.join(cur_dir, mainStrng)
    


    params_strng = mainStrng.replace('.mat', '_Params.mat')
    # Save the parameters to a .mat file
    # savemat(params_strng, {})

    # save the parameters into a pickle file
    with open(params_strng.replace('.mat', '.pkl'), 'wb') as f:
        # save params into a pickle file
        pickle.dump(params, f)

    # # Initialize empty variables
    # dummy variables: here they don't have any use, but in the monte-carlo simulation, they will be used, so I kept them here too.
    # nVar1 = 1
    # nVar2 = 1
    # z1MeanArr = np.full((nVar1, nVar2, nTVars, nMC), np.nan)
    # z1StdArr = np.full_like(z1MeanArr, np.nan)
    # posArr = np.full((nVar1, nVar2, 2, NPop, nTVars, nMC), np.nan)
    # stateArr = np.full((nVar1, nVar2, NPop, nTVars, nMC), np.nan)
    # zpArr = np.full_like(posArr, np.nan)
    # zsArr = np.full_like(posArr, np.nan)

    # Initialize positions randomly
    pos = np.vstack((np.random.uniform(initBounds[0, 0], initBounds[0, 1], NPop),
                    np.random.uniform(initBounds[1, 0], initBounds[1, 1], NPop)))

    # Run the simulation
    iMC = 0

    result = func_run_experiment(landFunc=landFunc, LB=landBounds, initPos=pos, linkThresh=linkThresh, tf=tf, nTVars=nTVars,
                                alpha=alpha, beta=beta, sigma=sigma, p2explt=p2explt, p2msngr=p2msngr,
                                randStep_msg=randStep_msg, randStep_exp=randStep_exp, stepSize=stepSize)



    # savemat(mainStrng, {
    #     'z1MeanArr': z1MeanArr,
    #     'z1StdArr': z1StdArr,
    #     'posArr': posArr,
    #     'stateArr': stateArr,
    #     'zpArr': zpArr,
    #     'zsArr': zsArr,
    # })

    # save the data into a pickle file
    with open(mainStrng.replace('.mat', '.pkl'), 'wb') as f:
        pickle.dump(result, f)

    return result, params



if __name__ == "__main__":
    print("Running one experiment from run_one_config.py", flush=True)
    result, params = __main__()
    print("Finished running one experiment")

    pp.show_animation(result=result, params=params)