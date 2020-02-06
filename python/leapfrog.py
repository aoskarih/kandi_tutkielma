import numpy as np
import matplotlib.pyplot as plt

dt = 0.1


def leap(x0, v0, f):
    v12 = v0 + f(x0)*dt*0.5
    x1 = x0 + v12*dt
    v1 = v12 + f(x1)*dt*0.5
    #print(f(x0))
    #print(f(x1))
    #print("v0:%f v12:%f v1:%f" % (v0, v12, v1))

    return [x1, v1]


def main(x0, v0):
    
    def pendulum(x): return -np.sin(x)
    
    x = [x0]
    v = [v0]
    
    c = (abs(v0)/3.0 ,0 , (3.0-abs(v0))/3.0)

    for _ in range(100):
        tmp = leap(x0, v0, pendulum)
        x0 = tmp[0]
        v0 = tmp[1]
        #print(tmp)
        x.append(x0)
        v.append(v0)
    
    plt.plot(x, v, "-", color=c, linewidth=2.0)



if __name__ == "__main__":
    for v in range(1, 30):
        main(0.0, v*0.1)
        main(-2*np.pi, v*0.1)
        main(0.0, -v*0.1)
        main(2*np.pi, -v*0.1)
    
    plt.xlabel("q")
    plt.ylabel("p")
    
    plt.xlim(-3.142, 3.142)
    plt.ylim(-3, 3)

    plt.show()
