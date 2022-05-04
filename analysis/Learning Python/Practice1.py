'''
Author: Alex Hopps

'''
#%%
from matplotlib import pyplot as plt
import numpy as np

x_vals = [1, 2, 3, 4]
y_vals = [5, 4, 6, 2]
x = np.linspace(-2*np.pi,2*np.pi,600)
plt.figure()
plt.plot(x_vals, y_vals)
plt.plot(x,np.sin(x))
plt.title("label")

plt.savefig("C:/Users/khopp/OneDrive/Documents/plot2.jpg")
# %%

# %%
