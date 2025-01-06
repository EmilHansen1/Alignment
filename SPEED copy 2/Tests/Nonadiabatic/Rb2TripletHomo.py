import numpy as np
import scipy as sp
import pandas as pd
from multiprocessing import Pool
from sympy.physics.quantum.cg import wigner_3j
from sympy.physics.wigner import Wigner3j
#from wigners import wigner_3j
from scipy.integrate import quad
from scipy.special import sph_harm
from itertools import starmap

jmax = 40
j0 = 0
jl1 = np.arange(j0, jmax, 2)    # Even Js (0, 2, ...)
jl2 = np.arange(j0+1, jmax, 2)  # Odd Js (1, 3, ...)

# Create lookup lists for Wigner-3j symbols
wgnj1 = np.zeros((jl1.shape[0], 2*jl1[-1]+1, 6))    # Even Js
wgnj2 = np.zeros((jl2.shape[0], 2*jl2[-1]+1, 6))    # Odd Js
for ij in range(len(jl1)):
    for mli in range(-jl1[ij], jl1[ij]+1):
        if (np.abs(mli) <= jl1[ij]-2):
            wgnj1[ij, mli+jl1[ij], 0] = Wigner3j(jl1[ij], jl1[ij]-2, 2, mli, -mli, 0).doit()
        wgnj1[ij, mli+jl1[ij], 1] = Wigner3j(jl1[ij], jl1[ij]-2, 2, mli, -(mli+2), 2).doit()
        wgnj1[ij, mli+jl1[ij], 2] = Wigner3j(jl1[ij], jl1[ij], 2, mli, -mli, 0).doit()
        wgnj1[ij, mli+jl1[ij], 3] = Wigner3j(jl1[ij], jl1[ij], 2, mli, -(mli+2), 2).doit()
        wgnj1[ij, mli+jl1[ij], 4] = Wigner3j(jl1[ij], jl1[ij]+2, 2, mli, -mli, 0).doit()
        wgnj1[ij, mli+jl1[ij], 5] = Wigner3j(jl1[ij], jl1[ij]+2, 2, mli, -(mli+2), 2).doit()
    for mli in range(-jl2[ij], jl2[ij]+1):
        if (np.abs(mli) <= jl2[ij]-2):
            wgnj2[ij, mli+jl2[ij], 0] = Wigner3j(jl2[ij], jl2[ij]-2, 2, mli, -mli, 0).doit()
        wgnj2[ij, mli+jl2[ij], 1] = Wigner3j(jl2[ij], jl2[ij]-2, 2, mli, -(mli+2), 2).doit()
        wgnj2[ij, mli+jl2[ij], 2] = Wigner3j(jl2[ij], jl2[ij], 2, mli, -mli, 0).doit()
        wgnj2[ij, mli+jl2[ij], 3] = Wigner3j(jl2[ij], jl2[ij], 2, mli, -(mli+2), 2).doit()
        wgnj2[ij, mli+jl2[ij], 4] = Wigner3j(jl2[ij], jl2[ij]+2, 2, mli, -mli, 0).doit()
        wgnj2[ij, mli+jl2[ij], 5] = Wigner3j(jl2[ij], jl2[ij]+2, 2, mli, -(mli+2), 2).doit()

size = [np.sum(jl1+1), np.sum(jl1), np.sum(jl2+1), np.sum(jl2)]
cosmat = [np.zeros((size[i],size[i])) for i in range(4)]
sincosmat = [np.zeros((size[i],size[i])) for i in range(4)]
for k in range(2):
    if k == 0:
        jl = jl1
        wgnj = wgnj1
    else:
        jl = jl2
        wgnj = wgnj2
    for i in range(len(jl)):
        sp = np.sum(jl[:i] + 1) if i > 0 else 0
        sc = jl[i] + 1
        cosl = cosmat[2*k]
        sincosl = sincosmat[2*k]
        for mli in range(-jl[i], jl[i] + 1, 2):
            mla = np.abs(mli)
            cosl[sp + (mli+jl[i])//2, sp + (mli+jl[i])//2] = 1/3 + (-1)**mla * (2/3)*(2*jl[i]+1)*wgnj[i,jl[i],2]*wgnj[i,mli+jl[i],2]
            sincosl[sp + (mli+jl[i])//2, sp + (mli+jl[i])//2] = 1/3 - (-1)**mla*(2*jl[i]+1)*(1/3)*wgnj[i,jl[i],2]*wgnj[i,mli+jl[i],2]
            if (jl[i] != jl[-1]):
                cosl[sp + sc + (mli+jl[i+1])//2, sp + (mli+jl[i])//2] = (-1)**mla * (2/3)*np.sqrt((2*jl[i]+1)*(2*jl[i+1]+1))*wgnj[i,jl[i],4]*wgnj[i,mli+jl[i],4]
                sincosl[sp + sc+(mli+jl[i+1])//2, sp + (mli+jl[i])//2] = -(-1)**mla * (1/3)*np.sqrt((2*jl[i]+1)*(2*jl[i+1]+1))*wgnj[i,jl[i],4]*wgnj[i,mli+jl[i],4]
                sincosl[sp + sc+(mli+2+jl[i+1])//2, sp + (mli+jl[i])//2]= (-1)**mla*np.sqrt((2*jl[i]+1)*(2*jl[i+1]+1))*np.sqrt(1/6)*wgnj[i,jl[i],4]*wgnj[i,mli+jl[i],5]
                sincosl[sp + sc+(mli-2+jl[i+1])//2, sp + (mli+jl[i])//2]= (-1)**mla*np.sqrt((2*jl[i]+1)*(2*jl[i+1]+1))*np.sqrt(1/6)*wgnj[i,jl[i],4]*wgnj[i,-mli+jl[i],5]
            if (np.abs(mli) != jl[i]):
                cosl[sp-sc+2+(mli+jl[i-1])//2, sp + (mli+jl[i])//2] = (-1)**mla * (2/3)*np.sqrt((2*jl[i]+1)*(2*jl[i-1]+1))*wgnj[i,jl[i],0]*wgnj[i,mli+jl[i],0]
                sincosl[sp-sc+2+(mli+jl[i-1])//2, sp + (mli+jl[i])//2]=-(-1)**mla * (1/3)*np.sqrt((2*jl[i]+1)*(2*jl[i-1]+1))*wgnj[i,jl[i],0]*wgnj[i,mli+jl[i],0]
            if (mli != jl[i]):
                sincosl[sp + (mli+2+jl[i])//2, sp + (mli+jl[i])//2]= (-1)**mla*(2*jl[i]+1)*np.sqrt(1/6)*wgnj[i,jl[i],2]*wgnj[i,mli+jl[i],3]
            if (mli != -jl[i]):
                sincosl[sp + (mli-2+jl[i])//2, sp + (mli+jl[i])//2]= (-1)**mla*(2*jl[i]+1)*np.sqrt(1/6)*wgnj[i,jl[i],2]*wgnj[i,-mli+jl[i],3]
            if (mli < jl[i]-2):
                sincosl[sp -sc+2+(mli+2+jl[i-1])//2, sp + (mli+jl[i])//2]=(-1)**mla*np.sqrt((2*jl[i]+1)*(2*jl[i-1]+1))*np.sqrt(1/6)*wgnj[i,jl[i],0]*wgnj[i,mli+jl[i],1]
            if (mli > -jl[i]+2):
                sincosl[sp -sc+2+(mli-2+jl[i-1])//2, sp + (mli+jl[i])//2]=(-1)**mla*np.sqrt((2*jl[i]+1)*(2*jl[i-1]+1))*np.sqrt(1/6)*wgnj[i,jl[i],0]*wgnj[i,-mli+jl[i],1]
            
    for i in range(len(jl)):
        sp = np.sum(jl[:i]) if i > 0 else 0
        sc = jl[i]
        cosl = cosmat[2*k+1]
        sincosl = sincosmat[2*k+1]
        for mli in range(-jl[i]+1, jl[i] + 1, 2):
            mla = np.abs(mli)
            cosl[sp + (mli+jl[i]-1)//2, sp + (mli+jl[i]-1)//2] = 1/3 + (-1)**mla * (2/3)*(2*jl[i]+1)*wgnj[i,jl[i],2]*wgnj[i,mli+jl[i],2]
            sincosl[sp + (mli+jl[i]-1)//2, sp + (mli+jl[i]-1)//2] = 1/3 - (-1)**mla*(2*jl[i]+1)*(1/3)*wgnj[i,jl[i],2]*wgnj[i,mli+jl[i],2]
            if (jl[i] != jl[-1]):
                cosl[sp + sc + (mli+jl[i+1]-1)//2, sp + (mli+jl[i]-1)//2] = (-1)**mla * (2/3)*np.sqrt((2*jl[i]+1)*(2*jl[i+1]+1))*wgnj[i,jl[i],4]*wgnj[i,mli+jl[i],4]
                sincosl[sp + sc+(mli+jl[i+1]-1)//2, sp + (mli+jl[i]-1)//2] = -(-1)**mla * (1/3)*np.sqrt((2*jl[i]+1)*(2*jl[i+1]+1))*wgnj[i,jl[i],4]*wgnj[i,mli+jl[i],4]
                sincosl[sp + sc+(mli+1+jl[i+1])//2, sp + (mli+jl[i]-1)//2]= (-1)**mla*np.sqrt((2*jl[i]+1)*(2*jl[i+1]+1))*np.sqrt(1/6)*wgnj[i,jl[i],4]*wgnj[i,mli+jl[i],5]
                sincosl[sp + sc+(mli-3+jl[i+1])//2, sp + (mli+jl[i]-1)//2]= (-1)**mla*np.sqrt((2*jl[i]+1)*(2*jl[i+1]+1))*np.sqrt(1/6)*wgnj[i,jl[i],4]*wgnj[i,-mli+jl[i],5]
            if (np.abs(mli) != jl[i]-1):
                cosl[sp-sc+2+(mli+jl[i-1]-1)//2, sp + (mli+jl[i]-1)//2] = (-1)**mla * (2/3)*np.sqrt((2*jl[i]+1)*(2*jl[i-1]+1))*wgnj[i,jl[i],0]*wgnj[i,mli+jl[i],0]
                sincosl[sp-sc+2+(mli+jl[i-1]-1)//2, sp + (mli+jl[i]-1)//2]=-(-1)**mla * (1/3)*np.sqrt((2*jl[i]+1)*(2*jl[i-1]+1))*wgnj[i,jl[i],0]*wgnj[i,mli+jl[i],0]
            if (mli != jl[i]-1):
                sincosl[sp + (mli+1+jl[i])//2, sp + (mli+jl[i]-1)//2]= (-1)**mla*(2*jl[i]+1)*np.sqrt(1/6)*wgnj[i,jl[i],2]*wgnj[i,mli+jl[i],3]
            if (mli != -jl[i]+1):
                sincosl[sp + (mli-3+jl[i])//2, sp + (mli+jl[i]-1)//2]= (-1)**mla*(2*jl[i]+1)*np.sqrt(1/6)*wgnj[i,jl[i],2]*wgnj[i,-mli+jl[i],3]
            if (mli < jl[i]-3):
                sincosl[sp -sc+2+(mli+1+jl[i-1])//2, sp + (mli+jl[i]-1)//2]=(-1)**mla*np.sqrt((2*jl[i]+1)*(2*jl[i-1]+1))*np.sqrt(1/6)*wgnj[i,jl[i],0]*wgnj[i,mli+jl[i],1]
            if (mli > -jl[i]+3):
                sincosl[sp -sc+2+(mli-3+jl[i-1])//2, sp + (mli+jl[i]-1)//2]=(-1)**mla*np.sqrt((2*jl[i]+1)*(2*jl[i-1]+1))*np.sqrt(1/6)*wgnj[i,jl[i],0]*wgnj[i,-mli+jl[i],1]


Hkin=[np.array([j*(j+1) for j in jl1 for m in range(-j, j + 1, 2)]),
     np.array([j*(j+1) for j in jl1 for m in range(-j+1, j + 1, 2)]),
     np.array([j*(j+1) for j in jl2 for m in range(-j, j + 1, 2)]),
     np.array([j*(j+1) for j in jl2 for m in range(-j+1, j + 1, 2)])]

def IntCosSq(jp, mp, j, m):
    return 0.5*np.pi*(-1)**np.abs(mp) *quad(lambda theta: np.real(sph_harm(-mp, jp, 0, theta) * sph_harm(m, j, 0, theta))*np.abs(np.sin(theta)), 0, np.pi)[0]

def VNa2(theta):
    ac = 29.4225
    cth = np.cos(theta)
    return ac*(-1.0+cth*np.arcsin(cth)/np.sqrt(1.0-cth**2)+ac*np.arcsin(cth)**2)

def IntCosPot(jp, mp, j, m):
    return 2.0*np.pi*(-1)**np.abs(mp) *quad(lambda theta: np.real(sph_harm(-mp, jp, 0, theta) * sph_harm(m, j, 0, theta))*np.abs(np.sin(theta)) * VNa2(theta), 0, np.pi)[0]

cosSqmatF = [np.zeros((np.sum(jl1+1),np.sum(jl1+1))), np.zeros((np.sum(jl1),np.sum(jl1))),
            np.zeros((np.sum(jl2+1),np.sum(jl2+1))), np.zeros((np.sum(jl2),np.sum(jl2)))]
basisF = [np.array([[j,m] for j in jl1 for m in range(-j, j + 1, 2)]),
         np.array([[j,m] for j in jl1 for m in range(-j+1, j + 1, 2)]),
         np.array([[j,m] for j in jl2 for m in range(-j, j + 1, 2)]),
         np.array([[j,m] for j in jl2 for m in range(-j+1, j + 1, 2)])]
for k in range(4):
    basis = basisF[k]
    cosSqmat = cosSqmatF[k]
    for i in range(len(basis)):
        for j in range(i,len(basis)):
            if(i==j):
                cosSqmat[i,j]=0.5
            elif(basis[i,1] == basis[j,1] + 2 or basis[i,1] == basis[j,1] - 2):
                cosSqmat[i,j] = IntCosSq(basis[i,0], basis[i,1], basis[j,0], basis[j,1])
    cosSqmatF[k] = cosSqmat + cosSqmat.T - np.diag(np.diagonal(cosSqmat))

cosVmatF = [np.zeros((np.sum(jl1+1),np.sum(jl1+1))), np.zeros((np.sum(jl1),np.sum(jl1))), 
            np.zeros((np.sum(jl2+1),np.sum(jl2+1))), np.zeros((np.sum(jl2),np.sum(jl2)))]
for k in range(4):
    basis = basisF[k]
    cosVmat = cosVmatF[k]
    for i in range(len(basis)):
        for j in range(i,len(basis)):
            if (basis[i,1] == basis[j,1]):
                cosVmat[i,j]= IntCosPot(basis[i,0], basis[i,1], basis[j,0], basis[j,1])
    cosVmatF[k] = cosVmat + cosVmat.T - np.diag(np.diagonal(cosVmat))

def It(I0, t0, tau, t):
    return I0 * np.exp(-4*np.log(2) * ((t - t0) / tau)**2)
    
def ftD(ti, t, ct, dt, step, idx, Bk, Bp, I0, t0, tau, dal, alp, w0, v0):
    trange = np.arange(ti, t, dt)
    res = np.zeros((int(trange.shape[0]/step), ct.shape[0]), dtype=complex)
    res[0] = ct
    for it, tl in enumerate(trange[1:]):
        if (np.abs(tl - t0) > 1.5 * tau):
            ctc = np.exp(-w0*dt*1j)*ct
        else:
            Il = It(I0, t0, tau, tl)
            Hmat = Bk*np.diag(Hkin[idx])+Bp*cosVmatF[idx]-Il * alp * np.eye(Hkin[idx].shape[0]) - Il * dal * sincosmat[idx]
            w,v = np.linalg.eigh(np.dot(np.conj(v0.T),np.dot(Hmat, v0)))
            ctc = np.dot(v, np.exp(-w*dt*1j) * np.dot(np.conj(v.T), ct))
        if it % step == 0:
            res[(it+1)//step] = ctc.copy()
        ct = ctc.copy()
    return trange[::step],res

def getCosNa2(Blk, Blp,  Il):
    print(Blk, "\t", Blp)
    Hmat = Blk*np.diag(Hkin[0])+Blp*cosVmatF[0]
    w0,v0 = np.linalg.eigh(Hmat)
    Hmat = Blk*np.diag(Hkin[1])+Blp*cosVmatF[1]
    w1,v1 = np.linalg.eigh(Hmat)
    Hmat = Blk*np.diag(Hkin[2])+Blp*cosVmatF[2]
    w2,v2 = np.linalg.eigh(Hmat)
    Hmat = Blk*np.diag(Hkin[3])+Blp*cosVmatF[3]
    w3,v3 = np.linalg.eigh(Hmat)
    w = [w0, w1, w2, w3]
    v = [v0, v1, v2, v3]

    def getCosj(idx, stind):
        ctl = np.zeros(v[idx].shape[0])
        ctl[stind] = 1
        cosMat = np.dot(np.conj(v[idx].T),np.dot(cosSqmatF[idx], v[idx]))
        ctr = ftD(ti, tmax, ctl, dt, step, idx, Blk, Blp, Il, t0, tau, dal, alp, w[idx], v[idx])
        zd1 = np.dot(cosMat, ctr[1].T)
        z = np.real(np.array([np.dot(np.conj(ctr[1][il]), zd1[:,il]) for il in range(ctr[1].shape[0])]))
        return np.real(z)

    thresh = 1e-4
    we = np.sort(np.concatenate((w0, w1)))
    wo = np.sort(np.concatenate((w2, w3)))
    me = 7/12
    mo = 5/12
    z = np.sum(me*np.exp(-we/Ta))+np.sum(mo*np.exp(-wo/Ta))
    prob  = [me * np.exp(-w0/Ta),
            me * np.exp(-w1/Ta),
            mo * np.exp(-w2/Ta),
            mo * np.exp(-w3/Ta)]
    indm = [np.where(pb/z<thresh)[0][0] for pb in prob]
    prob = [prob[i][:indm[i]] for i in range(4)]
    z = np.sum([np.sum(pb) for pb in prob])
    prob = [pb / z for pb in prob]
    lind = [np.array([i * np.ones(indm[i], dtype=int),np.arange(indm[i], dtype=int)]).T for i in range(4)]
    cos2d = [np.array(list(starmap(getCosj, lind[i]))) for i in range(4)]
    return np.sum(np.array([np.dot(prob[i], cos2d[i]) for i in range(4)]), axis=0)

 
I0 = 3.4*1e-2                   # Laser intensity
B = 0.28 * 2 * np.pi *1e-3      # B constant
Bg = 0.317 * 2 * np.pi *1e-3    # Related to V_He
cons = 2 * np.pi / (3 * 1.05)   # Converting intensity units?
Ir = cons * I0
ti = -20.0
t0 = 0.0
tau = 0.79 # 0.67
tmax = 1200.05
dal = 67.82
alp = 0
T = 0.38
Ta = 1.38*T*0.1 / 1.05
sB = 0.05 * 2 * np.pi *1e-3

dt = 0.1
step = 1
    
def getCosIBL(Bl):
    return getCosNa2(Bl, Bg, Ir)

def getCosIB():
    div = 200
    Blt = np.random.normal(B,sB,div)
    print(Blt)
    pool = Pool(processes=4)
    cos2d = np.array(pool.map(getCosIBL, Blt))
    return np.sum(cos2d, axis = 0) / div

if __name__ == '__main__':
    ctr = getCosIB()
    trange = np.arange(ti, tmax, dt)
    trange = trange[::step]
    actr = np.array([trange, ctr]).T
    pf = pd.DataFrame(actr[::])
    pf.to_csv("Rb3DTripletHomoB0.28D0.05I3.4e10WFA.csv")
