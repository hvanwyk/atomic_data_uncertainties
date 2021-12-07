import numpy as np
import matplotlib.pyplot as plt
import pickle
import corner 
import sys

from scipy.interpolate import CubicSpline, PPoly
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
#
# Load the temperature and rates samples
# 
data = pickle.load(open('../data/be-like-oxygen.pkl','rb'))
T = data['be_like_o_pos']['T']
r = data['be_like_o_pos']['rate_samples'].T

# 
# Replace 0 rate with small number
#
print('Minimum non-zero rate:', np.min(r[r!=0]))
print('Maximum rate:', np.max(r))
r[r==0] = 1e-50

# 
# Convert to logarithmic scale
# 
lT = np.log(T)
lr = np.log(r)

#
# Fit piecewise cubic spline to data 
# 
p = CubicSpline(lT,lr)

#
# Focus on the first interval
# 
n_test = 10000
i1_coef = p.c[:,15,:n_test]

#
# Fit a KDE to the coefficients in the first interval
# 


bandwidths = 10**np.linspace(-3,1,11)
grid = GridSearchCV(KernelDensity(kernel='gaussian'),{'bandwidth':bandwidths})
grid.fit(i1_coef.T)
bw_opt = grid.best_params_['bandwidth']

kde = KernelDensity(kernel='gaussian',bandwidth=bw_opt)
kde.fit(i1_coef.T)
print(kde.score_samples(i1_coef[:,1:10].T))

i1_sample = kde.sample(n_samples=10000)


figure01 = corner.corner(i1_coef.T)

figure02 = corner.corner(i1_sample)

plt.show()

print(sys.getsizeof(kde))

PP = PPoly(i1_coef[:,np.newaxis,:],lT[:2])  # Make piecewise polynomial from this
lt = np.linspace(lT[0],lT[1],101)
plt.plot(lt,PP(lt)[:,:20],'k-',lt,p(lt)[:,:20],'r:', linewidth=0.1)
plt.show()



"""
lT_min, lT_max = np.min(lT), np.max(lT)
lTT = np.linspace(lT_min, lT_max, 101)
plt.plot(lTT,p(lTT)[:,:100], 'k', lT,lr[:,:100],'.r',linewidth=0.1)
plt.show()
"""