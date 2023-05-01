# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 19:14:50 2022

@author: alpno
"""
distribution='uni'
tr_per_block = 150


distrs = dict(
            uni = np.repeat(1, len(FPs)), # uniform distribution
            exp = [8,4,2,1], # exponential FP distribution
            anti = [1,2,4,8], # anti-exponential FP distribution
            exp_ = [8,4,1], # exponential FP distribution (3 FPs)
            anti_ = [1, 4, 8], # anti-exponential FP distribution (3 FPs)
            gauss   = [1,5,1,1], #c.f. Trillenberg, last one is for catch;
            constant = [1], # constant blocks only have a single FP
        )

# if it's a 'famous' distribution, use predefined ratios
if distribution in list(distrs.keys()):
    distribution_ = distrs[distribution]
else:
    distribution_ = distribution

# Take some liberties fitting the precise distributions into the
# given tr_per_block; just treat these as proportions 
distribution = np.array(distribution_, dtype=float)
distribution /= distribution.sum() # proportions

# Assume FPs are in seconds.
FPs = np.array(FPs, dtype = float) # turns catch (None) into nan
np.array(FPs).size == distribution.size

# Make exp as a list of (1) block
full_exp = trs

props = np.round(distribution * tr_per_block).astype(int)
trs = np.repeat(FPs, props).tolist()
np.random.shuffle(trs)


from itertools import chain, repeat

def makeGoNoGoprop(trialPerBlock):
    ls = [int(0.8*trialPerBlock), int(0.2*trialPerBlock)]
    res = list(chain.from_iterable(repeat(j, times = i) for i, j in zip(ls, ['go', 'nogo'])))
    return res

gng = makeGoNoGoprop(trs.size)
gng = np.array(gng).flatten()

full_exp = list(zip(trs, gng))

np.random.shuffle(full_exp)
