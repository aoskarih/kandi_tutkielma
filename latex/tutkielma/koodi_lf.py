# Loikkakeino
def leapfrog(h, t, q0, p0, ddq):
    # liikemaara askeleen puolessa valissa
    p12 = p0 + ddq(t, q0)*h*0.5
    # paikka askeleen lopussa
    q1 = q0 + p12*h
    # liikemaara askeleen lopussa
    p1 = p12 + ddq(t, q1)*h*0.5
    return [q1, p1]

def calculate(steps, h, t0, q0, p0, ddq):
    q = [q0]
    p = [p0]
    t = t0
    # lasketaan arvot jokaiselle askeleelle ja lisataan ne listaan
    for _ in range(steps):
        tmp = leapfrog(h, t, q[-1], p[-1], ddq)
        q.append(tmp[0])
        p.append(tmp[1])
        t += h
    return q, p

def main():
    steps = 50
    h = 0.9
    # maaritellaan paikan toinen aikaderivaatta
    def ddq(t, q): return -q
    # lasketaan systeemin ja sen kokonaisenergian kehitys
    lf_plot = calculate(steps, h, 0, 0.0, 1.0, ddq)
    tot_e_lf = [(lf_plot[1][i]**2)/2 + (lf_plot[0][i]**2)/2 for i in range(steps)]
