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
import numpy as np
a= np.array([1, 2, 3, 4, 5])
print(a)
print(a.shape)
print(a.dtype)
print(a.ndim)
print(a.size)
print(a.itemsize)

# %%
a = np.array([1,2,3])
print(a[0])
a[0]=5
print(a[0])
b = a * np.array([2,0,2])
print(b)
print(b.sum())

# %%
#functions applied to a numpy array are like matrix math, whereas normal lists just append numbers
l=[1,2,3]
a=np.array([1,2,3])

l2=2*l
a2=2*a
print(l2) #apends an aaray to the end
print(a2) #multiplies each element
# %%
dot=np.dot(a,b)
print(dot)
cross=np.cross(a,b)
print(cross)

# %%
a = np.array([[1,2],[3,4]])
print(a.shape)
print(a)
print(a[0,:]) #the : means all column in the 0 row; you can also put a number there (0) to get the first column
# %%
bool = a > 2
print(a[bool]) #prints the elements of the array that meet the condition 
b = np.where(a>2, a, -1) #replace all elements that don't meet the condition with -1
print(b)
# %%

# %%
