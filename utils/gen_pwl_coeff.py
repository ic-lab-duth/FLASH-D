#%%
from matplotlib import pyplot as plt
import pwlf
import numpy as np
import torch
import math

x = np.linspace(-1, 1, 100)
print(x)
y = np.power(2, x)

# initialize piecewise linear fit with your x and y data
my_pwlf = pwlf.PiecewiseLinFit(x, y)

# my_pwlf.fit(4)

my_pwlf.fit_with_breaks([
#   -2, -1.75, -1.5, -1.25,
  -1.0, -0.75, -0.5, -0.25, 0,
 0.25, 0.50, 0.75, 1.00,
  # 1.25, 1.50, 1.75, 2.00, 
])
# my_pwlf.fit_with_breaks(np.linspace(-8, 8, 17))
# print(np.linspace(-8, 8, 17))

# predict for the determined points
xHat = np.linspace(min(x), max(x), num=10000)
yHat = my_pwlf.predict(xHat)

# plot the results
plt.figure()
plt.plot(x, y, 'o')
plt.plot(xHat, yHat, '-')
plt.show()

#%%
from sympy import Symbol
from sympy.utilities import lambdify
x = Symbol('x')

def get_symbolic_eqn(pwlf_, segment_number):
    if pwlf_.degree < 1:
        raise ValueError('Degree must be at least 1')
    if segment_number < 1 or segment_number > pwlf_.n_segments:
        raise ValueError('segment_number not possible')
    # assemble degree = 1 first
    for line in range(segment_number):
        if line == 0:
            my_eqn = pwlf_.beta[0] + (pwlf_.beta[1])*(x-pwlf_.fit_breaks[0])
        else:
            my_eqn += (pwlf_.beta[line+1])*(x-pwlf_.fit_breaks[line])
    # assemble all other degrees
    if pwlf_.degree > 1:
        for k in range(2, pwlf_.degree + 1):
            for line in range(segment_number):
                beta_index = pwlf_.n_segments*(k-1) + line + 1
                my_eqn += (pwlf_.beta[beta_index])*(x-pwlf_.fit_breaks[line])**k
    return my_eqn.simplify()


eqn_list = []
f_list = []
for i in range(my_pwlf.n_segments):
    eqn_list.append(get_symbolic_eqn(my_pwlf, i + 1))
    print('Equation number: ', i + 1)
    print(eqn_list[-1])
    f_list.append(lambdify(x, eqn_list[-1]))
# %%
