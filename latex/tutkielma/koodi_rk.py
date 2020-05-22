# Runge-Kutta
def rk4(h, t, q0, p0, dq, dp):
    # lasketaan kulmakertoimet
    k1 = h*dq(t, q0, p0)
    l1 = h*dp(t, q0, p0)

    k2 = h*dq(t+0.5*h, q0+0.5*k1, p0+0.5*l1)
    l2 = h*dp(t+0.5*h, q0+0.5*k1, p0+0.5*l1)

    k3 = h*dq(t+0.5*h, q0+0.5*k2, p0+0.5*l2)
    l3 = h*dp(t+0.5*h, q0+0.5*k2, p0+0.5*l2)

    k4 = h*dq(t+h, q0+k3, p0+l3)
    l4 = h*dp(t+h, q0+k3, p0+l3)
    # paikka ja liikemaara askeleen lopussa
    q1 = q0 + (k1 + 2*k2 + 2*k3 + k4)/6.0
    p1 = p0 + (l1 + 2*l2 + 2*l3 + l4)/6.0
    return [q1, p1]

def calculate(steps, h, t0, q0, p0, dq, dp):
    q = [q0]
    p = [p0]
    t = t0
    # lasketaan arvot jokaiselle askeleelle ja lisataan ne listaan
    for _ in range(steps):
        tmp = rk4(h, t, q[-1], p[-1], dq, dp)
        q.append(tmp[0])
        p.append(tmp[1])
        t += h
    return q, p

def main():
    steps = 50
    h = 0.9
    # maaritellaan paikan ja liikemaaran aikaderivaatat
    def dq(t, q, p): return p
    def dp(t, q, p): return -q
    # lasketaan systeemin ja sen kokonaisenergian kehitys
    rk_plot = calculate(steps, h, 0, 0.0, 1.0, dq, dp)
    tot_e_rk = [(rk_plot[1][i]**2)/2 + (rk_plot[0][i]**2)/2 for i in range(steps)]
