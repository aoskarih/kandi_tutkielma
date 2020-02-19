import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import mpl_toolkits.mplot3d.art3d as mpl3d

fig = plt.figure(figsize=(16, 9), edgecolor='w')
plt.rcParams['axes.grid'] = True
plt.rcParams['font.size'] = 26


def rk4(h, t, q0, p0, dq, dp):
    
    k1 = h*dq(t, q0, p0)
    l1 = h*dp(t, q0, p0)

    k2 = h*dq(t+0.5*h, q0+0.5*k1, p0+0.5*l1)
    l2 = h*dp(t+0.5*h, q0+0.5*k1, p0+0.5*l1)

    k3 = h*dq(t+0.5*h, q0+0.5*k2, p0+0.5*l2)
    l3 = h*dp(t+0.5*h, q0+0.5*k2, p0+0.5*l2)

    k4 = h*dq(t+h, q0+k3, p0+l3)
    l4 = h*dp(t+h, q0+k3, p0+l3)

    q1 = q0 + (k1 + 2*k2 + 2*k3 + k4)/6.0
    p1 = p0 + (l1 + 2*l2 + 2*l3 + l4)/6.0

    return [q1, p1]

def rk4_3(h, t, x0, y0, z0, dx, dy, dz):
    
    k1 = h*dx(t, x0, y0, z0)
    l1 = h*dy(t, x0, y0, z0)
    m1 = h*dz(t, x0, y0, z0)

    k2 = h*dx(t+0.5*h, x0+0.5*k1, y0+0.5*l1, z0+0.5*m1)
    l2 = h*dy(t+0.5*h, x0+0.5*k1, y0+0.5*l1, z0+0.5*m1)
    m2 = h*dz(t+0.5*h, x0+0.5*k1, y0+0.5*l1, z0+0.5*m1)

    k3 = h*dx(t+0.5*h, x0+0.5*k2, y0+0.5*l2, z0+0.5*m2)
    l3 = h*dy(t+0.5*h, x0+0.5*k2, y0+0.5*l2, z0+0.5*m2)
    m3 = h*dz(t+0.5*h, x0+0.5*k2, y0+0.5*l2, z0+0.5*m2)

    k4 = h*dx(t+h, x0+k3, y0+l3, z0+m3)
    l4 = h*dy(t+h, x0+k3, y0+l3, z0+m3)
    m4 = h*dz(t+h, x0+k3, y0+l3, z0+m3)

    x1 = x0 + (k1 + 2*k2 + 2*k3 + k4)/6.0
    y1 = y0 + (l1 + 2*l2 + 2*l3 + l4)/6.0
    z1 = z0 + (m1 + 2*m2 + 2*m3 + m4)/6.0

    return [x1, y1, z1]

def euler(h, t, q0, p0, dq, dp):
    
    k1 = dq(t, q0, p0)
    l1 = dp(t, q0, p0)

    q1 = q0 + k1*h
    p1 = p0 + l1*h

    return [q1, p1]


def ham_calculate(method, steps, h, t0, q0, p0, dq, dp):
    
    q = [q0]
    p = [p0]

    t = t0

    for _ in range(steps):
        tmp = method(h, t, q[-1], p[-1], dq, dp)
        q.append(tmp[0])
        p.append(tmp[1])
        t += h

    return q, p

def ham3_calculate(method, steps, h, t0, x0, y0, z0, dx, dy, dz):
    
    x = [x0]
    y = [y0]
    z = [z0]

    t = t0

    for _ in range(steps):
        tmp = method(h, t, x[-1], y[-1], z[-1], dx, dy, dz)
        x.append(tmp[0])
        y.append(tmp[1])
        z.append(tmp[2])
        t += h

    return x, y, z

def lag_calculate(method, steps, h, t0, q0, p0, ddq):
    
    q = [q0]
    p = [p0]

    t = t0

    for _ in range(steps):
        tmp = method(h, t, q[-1], p[-1], ddq)
        q.append(tmp[0])
        p.append(tmp[1])
        t += h

    return q, p


def leapfrog(h, t, q0, p0, ddq):
    
    p12 = p0 + ddq(t, q0)*h*0.5
    q1 = q0 + p12*h
    p1 = p12 + ddq(t, q1)*h*0.5
    
    return [q1, p1]

def main():
    
    def dq(t, q, p): return p
    def dp(t, q, p): return -np.sin(q)
    def ddq(t, q): return -np.sin(q)

    anal_plot = ham_calculate(rk4, 7000, 0.001, 0, 0.0, 1.0, dq, dp)
    e_plot = ham_calculate(euler, 2000, 0.05, 0, 0.0, 1.0, dq, dp)
    rk4_plot = ham_calculate(rk4, 1000, 0.01, 0, 0.0, 1.0, dq, dp)
    lf_plot = lag_calculate(leapfrog, 1000, 0.01, 0, 0.0, 1.0, ddq)

    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    
    ax1.plot(lf_plot[0], lf_plot[1], "-", lw=2.5, label="Leapfrog")
    ax2.plot(rk4_plot[0], rk4_plot[1], "-", lw=2.5, label="Runge-Kutta", color="green")
    ax1.plot(anal_plot[0], anal_plot[1], "--", lw=3.5, label="Analyyttinen", color="red")
    ax2.plot(anal_plot[0], anal_plot[1], "--", lw=3.5, label="Analyyttinen", color="red")

    #lc1 = colorline(lf_plot[0], lf_plot[1], cm="summer")
    #ax1.add_collection(lc1)

    ax1.set_title("")
    ax1.set_xlabel("Paikka")
    ax1.set_ylabel("Liikemäärä")
    ax1.set_xlim(-1.5, 1.5)
    ax1.set_ylim(-1.5, 1.5)
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    

    
    #lc2 = colorline(rk4_plot[0], rk4_plot[1], cm="summer")
    #ax2.add_collection(lc2)

    ax2.set_title("")
    ax2.set_xlabel("Paikka")
    ax2.set_ylabel("Liikemäärä")
    ax2.set_xlim(-1.5, 1.5)
    ax2.set_ylim(-1.5, 1.5)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])

    ax1.legend(loc="upper right", fontsize="small")
    ax2.legend(loc="upper right", fontsize="small")
    

def lotka_volterra():
    a = 2.0/3.0
    b = 4.0/3.0
    d = 1
    g = 1

    def dx(t, x, y): return a*x-b*x*y
    def dy(t, x, y): return d*x*y-g*y

    plt.rcParams['axes.grid'] = False
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20) 

    ax2 = fig.add_subplot(111)
    for i in range(0, 15, 3):
        rk4_plot = ham_calculate(rk4, 4000, 0.005, 0, 1.0, float(i+5)/10, dx, dy)    
        ax2.plot(rk4_plot[0], rk4_plot[1], "-", color=(float(15-i)/15, 0.0, float(i)/15), linewidth=2)

    
    ax2.set_title("")
    ax2.set_xlabel("Saaliseläimet")
    ax2.set_ylabel("Saalistajat")
    ax2.set_xlim(-0.1, 3)
    ax2.set_ylim(-0.1, 2)
    #ax2.set_xticklabels([])
    #ax2.set_yticklabels([0, 1, 2, 3, 4])
    ax2.yaxis.set_ticks(np.arange(0, 2.5, 0.5))

def scroll_attractor():
    def dx(t, x, y, z): return 40*(y-x)
    def dy(t, x, y, z): return (28-40)*x - x*z + 28*y
    def dz(t, x, y, z): return x*y-3*z
    # -0.1, 0.5, -0.6
    scroll_plot = ham3_calculate(rk4_3, 40000, 0.001, 0, -0.1, 0.5, -0.6, dx, dy, dz)

    lc = colorline3(scroll_plot[0], scroll_plot[1], scroll_plot[2])
    ax = fig.add_subplot(111, projection="3d")
    ax.add_collection(lc)
    ax.set_xlim(-40, 40)
    ax.set_ylim(-30, 30)
    ax.set_zlim(0, 40)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    ax._axis3don = False

def lorenz_attractor():
    def dx(t, x, y, z): return 10*(y-x)
    def dy(t, x, y, z): return (28-z)*x - y
    def dz(t, x, y, z): return x*y-8.0/3.0*z
    # -0.1, 0.5, -0.6
    scroll_plot = ham3_calculate(rk4_3, 80000, 0.001, 0, -0.1, 0.5, -0.6, dx, dy, dz)

    lc = colorline3(scroll_plot[0], scroll_plot[1], scroll_plot[2], cm=plt.get_cmap('plasma'))
    ax = fig.add_subplot(111, projection="3d")
    ax.add_collection(lc)
    ax.set_xlim(-40, 40)
    ax.set_ylim(-30, 30)
    ax.set_zlim(0, 40)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    ax._axis3don = False

def euler_ball():
    def dq(t, q, p): return p
    def dp(t, q, p): return -1

    plt.rcParams['axes.grid'] = False
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20) 

    dt = 0.3
    s = 22

    ax2 = fig.add_subplot(111)
    
    e_plot = ham_calculate(euler, s, dt, 0, 0.0, 1.0, dq, dp)
    anal_plot = ham_calculate(rk4, 2500, 0.001, 0, 0.0, 1.0, dq, dp)

    ax2.plot(np.arange(0, (s+0.5)*dt, dt), e_plot[0], "-", color="blue", label="Euler")
    ax2.plot(np.arange(0, (2500+0.5)*0.001, 0.001), anal_plot[0], "-", color="red", label="Analyyttinen")

    ax2.set_title("")
    ax2.set_xlabel("Aika")
    ax2.set_ylabel("Korkeus")
    ax2.set_xlim(0.0, 2.5)
    ax2.set_ylim(0.0, 1)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    plt.legend()
    #ax2.xaxis.set_ticks(np.arange(0, 2.5, 0.3))

def colorline(x, y, cm=plt.get_cmap('plasma'), norm=plt.Normalize(0.0, 1.0), linewidth=2):
    w = np.linspace(0.0, 1.0, len(x))

    if not hasattr(w, "__iter__"):
        w = np.array([w])

    w = np.asarray(w)

    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=w, cmap=cm, norm=norm, linewidth=linewidth)

    return lc

def colorline3(x, y, z, cm=plt.get_cmap('plasma'), norm=plt.Normalize(0.0, 1.0), linewidth=2):
    w = np.linspace(0.0, 1.0, len(x))

    if not hasattr(w, "__iter__"):
        w = np.array([w])

    w = np.asarray(w)

    segments = make_segments3(x, y, z)
    lc = mpl3d.Line3DCollection(segments, array=w, cmap=cm, norm=norm, linewidth=linewidth)

    return lc

def make_segments(x, y):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

def make_segments3(x, y, z):
    points = np.array([x, y, z]).T.reshape(-1, 1, 3)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

if __name__ == "__main__":
    
    main()
    #lorenz_attractor()
    #scroll_attractor()
    #lotka_volterra()
    #euler_ball()
    fig.tight_layout()
    plt.show()
