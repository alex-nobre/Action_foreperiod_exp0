# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 11:36:17 2023

@author: alpno
"""

x=np.repeat(1, len(FPs))

x

x=np.array(x, dtype=float)

x

x /= x.sum()

x

y=FPs

y

y = np.array(FPs, dtype = float)

y

z= np.round(x * 150).astype(int)

z

alfa=np.repeat(y, z).tolist()

alfa

np.random.shuffle(alfa)

alfa

while np.isnan(alfa[0]):
    np.random.shuffle(alfa)