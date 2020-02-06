import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll
import matplotlib.path as mpath

fig = plt.figure(figsize=(18, 9), edgecolor='w')
plt.rcParams['axes.grid'] = True


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

def euler(h, t, q0, p0, dq, dp):
    
    k1 = h*dq(t, q0, p0)
    l1 = h*dp(t, q0, p0)

    q1 = q0 + k1
    p1 = p0 + l1

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

    e_plot = ham_calculate(euler, 10, 0.01, 0, 1.0, 1.0, dq, dp)
    rk4_plot = ham_calculate(rk4, 100, 0.5, 0, 0.0, 1.0, dq, dp)
    lf_plot = lag_calculate(leapfrog, 100, 0.5, 0, 0.0, 1.0, ddq)

    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    
    lc1 = colorline(lf_plot[0], lf_plot[1])
    ax1.add_collection(lc1)
    ax1.set_title("Leapfrog")
    ax1.set_xlabel("q")
    ax1.set_ylabel("p")
    ax1.set_xlim(-3, 3)
    ax1.set_ylim(-3, 3)

    lc2 = colorline(rk4_plot[0], rk4_plot[1])
    ax2.add_collection(lc2)
    ax2.set_xlabel("q")
    ax2.set_ylabel("p")
    ax2.set_xlim(-3, 3)
    ax2.set_ylim(-3, 3)

def colorline(x, y, cm=plt.get_cmap('plasma'), norm=plt.Normalize(0.0, 1.0), linewidth=1):
    z = np.linspace(0.0, 1.0, len(x))

    if not hasattr(z, "__iter__"):
        z = np.array([z])

    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, cmap=cm, norm=norm, linewidth=linewidth)

    return lc


def make_segments(x, y):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


if __name__ == "__main__":
    
    main()
    plt.show()
