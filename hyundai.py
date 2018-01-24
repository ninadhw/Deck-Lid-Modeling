#hyundai code for deck lid
import numpy as np
from scipy.integrate import ode
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
k1=12.9*0.01
k2=10.5*0.01

Ip = (np.pi * d**4)/32
bk = G*Ip/L

def calculate_angle_41 (a1,a2,a3,a4,theta21):
    #calculate angle theta 41
    
    K1= 2*a1*a4
    K2= -2*a2*a4
    K3= a1**2 + a2**2 - a3**2 + a4**2 - 2*a1*a2* np.cos( np.deg2rad(theta21))

    a=K1 + K2*np.cos( np.deg2rad(theta21))
    b=K2*np.sin( np.deg2rad(theta21))
    c=K3
    print(a,b,c)

    tan= np.roots([(c-a),(2*b),(c+a)])
    atan= np.arctan2(tan,[1,1])
    atan=2*atan
    atan=np.rad2deg(atan)
    return [atan[1]]

##ni=180 - (phi_0 + sci_0)
##theta41 = calculate_angle_41(r1,r2,r3,r4,(180-20-40))
##theta31 = calculate_angle_41(r2,r1,r4,r3,(180-20-40))
##
##print(theta31)
##print(theta41)
th2=[]
th3=[]
th4=[]
tor_mom=[]
def solvern(t , U):
    theta41 = calculate_angle_41 (r1,r2,r3,r4,U[0])
    theta31 = calculate_angle_41 (r2,r1,r4,r3,U[0])
    theta31 = -theta31[0] + U[0]
    th2.append(U[0])
    th3.append(theta31)
    th4.append(theta41)

    Mw= W* ((rg*np.cos(np.deg2rad(theta_g-U[0]))) + (r2*np.cos(np.deg2rad(180-sci_0-U[0]))))
    torsionbar= sci_0 + beta_0 + theta41[0]
    Mb=bk*(np.deg2rad(torsionbar))
    Mt=2*Mb*(r2*np.sin(np.deg2rad(U[0] - theta31))) /(r4*np.sin(np.deg2rad(theta41[0] - theta31)))
    tor_mom.append(Mb)
    
    return [U[1], ((Mt - Mw - k1*U[1] - k2*U[1]*U[1])/Ik) ]


r=ode(solvern, jac=None).set_integrator('dopri5', nsteps=1000)
r.set_initial_value([120,0],0)
t1=100
dt=0.1


i=0
disp=[]
vel=[]
while r.successful() and r.t<t1:
    b=r.integrate(r.t+dt)
    disp.append(b[0])
    vel.append(b[1])
    i=i+1

op_ang = 120*np.ones(len(th2)) - th2
   
tzone=np.linspace(0,t1,(t1/dt))
#mpl.plot(op_ang,th3,'b')
#mpl.show()
#mpl.plot(op_ang,th4,'g')
#mpl.show()
#mpl.plot(np.linspace(0,100,len(op_ang)),op_ang,'g')
mpl.plot(np.linspace(0,100,len(disp)),disp)
#mpl.plot(th2,)
mpl.show()
