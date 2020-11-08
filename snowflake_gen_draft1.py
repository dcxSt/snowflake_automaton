# python snowflake generation version 1

# imports 
import numpy as np
import matplotlib.pyplot as plt
import os

"""
So the grid is hexagonal, and how we're gonna make it is a 
square array, where every other line is displaed to the left 
or to the right a tiny bit.
The first array is slightly to the left and the second slightly 
to the right. Since we start at zero, the arrays with even index 
will be displaced to the left and those with odd things will be 
slightly to the right.
"""

# initializes the grid
def init_grid(beta,n=101):
    # n is the size of the grid
    grid = np.ones((n,n)) * beta
    # assume n is odd, fill middle cell with value 1
    mid = int((n-1)/2)
    grid[mid][mid] = 1
    
    return grid

# plot the grid
def plot_grid(grid,saveas=None):
#     normalized = grid / max(grid.flatten())
    # invert the color 
    n = grid.shape[0]
    max_val = max(grid.flatten())
    inverted_grid = np.zeros((n,n))
    for i,row in enumerate(grid):
        for j,element in enumerate(row):
            inverted_grid[i][j] = max_val - element
    
    # plot it
    plt.figure(figsize=(8,8))
    plt.imshow(inverted_grid, cmap='Greys', interpolation='nearest')
    #plt.colorbar()
    plt.title("display of density")
        
    # save the figure
    if saveas:
        plt.savefig(saveas)
    plt.show()

# convert a gradient thing to a binary solid/liquid
def to_snowflake(grid):
    n = grid.shape[0]
    snowflake_grid = np.zeros((n,n))
    for i,row in enumerate(grid):
        for j,element in enumerate(row):
            if element >= 1: snowflake_grid[i][j] = 1
    return snowflake_grid

# helper for update grid, adds to the receptive cells
def average(grid,alpha):
    n = grid.shape[0]
    averaged_grid = np.zeros((n,n))
    for i,row in enumerate(grid):
        for j,element in enumerate(row):
            averaged_grid[i][j] += (1-alpha*0.5)*element
            
            averaged_grid[(i+1)%n][j] += alpha/12*element
            averaged_grid[i][(j-1)%n] += alpha/12*element
            averaged_grid[i][(j+1)%n] += alpha/12*element
            averaged_grid[(i-1)%n][j%n] += alpha/12*element
            
            if i%2==0:# if the row is even index it's to the left
                averaged_grid[(i+1)%n][(j-1)%n] += alpha/12*element
                averaged_grid[(i-1)%n][(j-1)%n] += alpha/12*element
                
            else:
                averaged_grid[(i+1)%n][(j+1)%n] += alpha/12*element
                averaged_grid[(i-1)%n][(j+1)%n] += alpha/12*element
    return averaged_grid

def get_grid_is_receptive(grid,n):
    grid_is_receptive = np.zeros(())

# update the grid, this is the discrete update function / diff-eq
def timestep(grid,alpha,gamma):
    n = grid.shape[0]
    
    # find the receptive cubes
    grid_is_receptive = np.zeros((n,n))
    for i,row in enumerate(grid):
        for j,element in enumerate(row):
            if element >= 1:
                grid_is_receptive[i][j] = 1
                grid_is_receptive[(i+1)%n][j] = 1
                grid_is_receptive[(i-1)%n][j] = 1
                grid_is_receptive[i][(j-1)%n] = 1
                grid_is_receptive[i][(j+1)%n] = 1
                if i%2==0:
                    grid_is_receptive[(i+1)%n][(j-1)%n] = 1
                    grid_is_receptive[(i-1)%n][(j-1)%n] = 1
                else:
                    grid_is_receptive[(i+1)%n][(j+1)%n] = 1
                    grid_is_receptive[(i-1)%n][(j+1)%n] = 1
                    
    # make receptive and non receptive grids
    grid_receptive = np.zeros((n,n))
    grid_non_receptive = np.zeros((n,n))
    for i,row in enumerate(grid):
        for j,element in enumerate(row):
            if grid_is_receptive[i][j]==1:
                grid_receptive[i][j] = element + gamma
            else:
                grid_non_receptive[i][j] = element
    
    
    grid_non_receptive = average(grid_non_receptive,alpha)
    return grid_receptive + grid_non_receptive

# checks to see if the grid is too big, if it is, returns True, if not it returns false
def snowflake_too_big(grid):
    n = grid.shape[0]
    mid = int((n-1)/2)
    one_ninth = int(n/9)
    if grid[mid][one_ninth] >= 1:
        return True
    return False

# mades frame aka grid bigger by a padded margin of 20, the padding has thickness 20
def pad(grid,beta,pad_size=20):
    n = grid.shape[0]
    bigger_grid = np.ones((n+2*pad_size,n+2*pad_size))*beta
    for i in range(n):
        for j in range(n):
            bigger_grid[i+pad_size][j+pad_size] = grid[i][j]
    return bigger_grid

def make_snowflake_pngs_for_gif(alpha,beta,gamma,iterations,n):
    # open a new directory
    dir_name = "snowflake_pngs_alpha={}_beta={}_gamma={}".format(alpha,beta,gamma)
    os.mkdir(dir_name)
    grid = init_grid(beta,n)
    for i in range(iterations):
        grid = timestep(grid,alpha,gamma)
        if i % 100 == 0:
            # chekc if it's too big, if so stop the loop
            if snowflake_too_big(grid):
                break
            
            # save the figure
            n = grid.shape[0]
            saveas = "{}/snowflake_alpha={}_beta={}_gamma={}_iter={}_n={}.png".format(dir_name,alpha,beta,gamma,i,n)
            plot_grid(grid,saveas)
            np.save("{}/snowflake_alpha={}_beta={}_gamma={}_iter={}_n={}".format(dir_name,alpha,beta,gamma,i,n),grid)
            
            
    # once we're done, save and display one final copy
    n = grid.shape[0]
    saveas = "snowflake_alpha={}_beta={}_gamma={}_iter={}_n={}.png".format(alpha,beta,gamma,i,n)
    plot_grid(grid,saveas)
    np.save("snowflake_alpha={}_beta={}_gamma={}_iter={}_n={}".format(alpha,beta,gamma,i,n),grid)
    return

if __name__=="__main__":
    alpha = 2.43 # this is the diffusion constant
    beta = 0.35 # background level, initial water density in surrounding atmosphere
    gamma = 0.001 # accreation from outside to account for 3rd dimension
    max_iter = 10000 
    n = 450 # size of frame