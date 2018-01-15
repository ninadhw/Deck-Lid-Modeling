#hyundai code for deck lid
import numpy as np
import scipy.integrate as sp
import matplotlib.pyplot as mpl

r1=80           #mm
r2=75
r3=85
r4=75
d=6
L=932
rg=420.8
rh=377.9
theta_g=23.7
theta_h=37.6
sci_0=20
phi_0=40
beta_0=10.8
W=131.5
Ik=21355
G=75000
k1=12.9
k2=10.5

Ip = ( np.pi * d**4)/32
bk = G*Ip/L

def calculate_angle_41 (a1,a2,a3,a4,theta21):
    #calculate angle theta 41
    K1= 2*a1*a4
    K2= -2*a2*a4
    K3= a1**2 + a2**2 - a3**2 + a4**2 - 2*a1*a2* np.cos( np.deg2rad(theta21))

    a=K1 + K2*np.cos( np.deg2rad(theta21))
    b=K2*np.sin( np.deg2rad(theta21))
    c=K3

    tan= np.roots([(c-a),(2*b),(c+a)])
    atan= np.arctan2(tan,[1,1])
    atan=2*atan
    atan=np.rad2deg(atan)
    return [atan[1]]

ni=180 - (phi_0 + sci_0)


# Theta 41 = calculate_angle_41(r1,r2,r3,r4,ni)
# theta41=np.matrix.transpose(np.rad2deg(ninad))
# Theta 31 = calculate_angle_41(r2,r1,r4,r3,ni)
# theta31=np.matrix.transpose(np.rad2deg(ninad))

def solvern(U, t):
    Mw= W* ( (rg*np.cos(np.deg2rad(theta_g - U[0]))) + (r2*np.cos(np.deg2rad(180 - sci_0 - U[0]))) )
    theta41 = calculate_angle_41 (r1,r2,r3,r4,U[0])
    theta31 = calculate_angle_41 (r2,r1,r4,r3,U[0])
    theta31 = theta31[0] - U[0]
    torsionbar= sci_0 + beta_0 + theta41[0]
    Mb=bk*(np.deg2rad(torsionbar))
    Mt=2*Mb*(r2*np.sin(np.deg2rad(U[0] - theta31)) ) /  (r4*np.sin(np.deg2rad(theta41[0] - theta31)))
    return [U[1], (     (-k1*U[1]) + (-k2*U[1]**2) + (Mw-Mt)     )/Ik]

xs = np.linspace(0,50,500)
Us = sp.odeint(solvern, [120, 0], xs)

U1 = Us[:,1]

mpl.plot(xs,U1)
mpl.show()
print(Us)
