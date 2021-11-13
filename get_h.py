import math as m

def calc_discharge(b, h, m_bank, S_0, **kwargs):
    A = h*(b+h*m_bank)
    P = b+2*h*(m.sqrt(m_bank**2 +1))
    Rh = A/P

    k_st = 20
    for k in kwargs.items():
        if "k" in k[0]:
            try:
                k_st = str(k[1])
                print("Using k_st = " + str(k[1]))
            except:
                print(str(k[1]) + " is not a number (using default value 20).")

    for n in kwargs.items():
        if "n" in n[0]:
            try:
                k_st = str(1/float(n[1]))
                print("Using k_st = " + str(1/float(n[1])))
            except:
                print(str(1/float(n[1])) + " is not a number (using default value 20).")

    for d in kwargs.items():
        if "d" in d[0]:
            try:
                k_st = str(26/pow(d[1],1/6))
                print("Using k_st = " + str(26/pow(d[1],1/6)))
            except:
                print(str(26/pow(d[1],1/6)) + " is not a number (using default value 20).")

    Q = float(k_st)*m.sqrt(S_0) * (Rh**(2/3))*A
    return Q

def interpolate_h(Q, b, m_bank, S, k_st):

    eps=1
    count=100
    n_m=1/k_st
    h=1

    while eps > 10**-3 and count > 0:
        A = h*(b+h*m_bank)
        P = b+2*h*(m.sqrt(m_bank**2 + 1))
        Qk = (A**(5/3)*m.sqrt(S))/(n_m*P**(2/3))
        eps = abs(Q-Qk)/Q

        # Derivate of A
        dA_dh = b + 2*m_bank*h

        # Derivate of P
        dP_dh = 2*m.sqrt(m_bank**2 + 1)

        # Function that should become zero
        F = (n_m*Qk*P**(2/3)) - (A**(5/3) *m.sqrt(S))

        # It's derivate
        dF_dh = ((2/3)*n_m*Qk * P**(-1/3))*dP_dh - ((5/3)* A**(2/3) * m.sqrt(S))*dA_dh

        # Water depth Update
        h = abs((h-F)/dF_dh)
        count -= 1
    return h

if __name__ == '__main__':
    # input parameters
    Q = 15.5        # discharge in (m3/s)
    b = 5.1         # bottom channel width (m)
    m_bank = 2.5    # bank slope
    S_0 = 0.005     # channel slope
    h = 1           # h initial

    Answer = calc_discharge(b, h, m_bank, S_0, k=20)
    print("The result is: %0.2f m3/s" % Answer)

    # call the solver with user-defined channel geometry and discharge
    k_st=20
    h_n = interpolate_h(Q, b, m_bank, S_0, k_st)
    print("The result of the interpolation of h is: %0.2f m3/s" % h_n)

    Answer2 = calc_discharge(b, h_n, m_bank, S_0, k=20)
    print("Finally, the result of Q using the interpolated h is %0.2f m3/s" % Answer2)