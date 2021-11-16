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
    return [Q,float(k_st)]

def interpolate_h(Q, b, m_bank, S_0, **kwargs):

    eps=1
    h = 1
    count=5000

    k_st = 20
    for k in kwargs.items():
        if "k" in k[0]:
            try:
                k_st = str(k[1])
                print("Using k_st = " + str(k[1]))
            except:
                print(str(k[1]) + " is not a number (using default value 20).")
    n_m=1/float(k_st)

    while eps > 10**-3 and count > 0:

        A = h*(b+h*m_bank)
        P = b+2*h*(m.sqrt(m_bank**2 + 1))
        Qk = (A**(5/3)*m.sqrt(S_0))/(n_m*P**(2/3))
        eps = abs(Q-Qk)/Q

        # Derivate of A
        dA_dh = b + 2*m_bank*h

        # Derivate of P
        dP_dh = 2*(m.sqrt(m_bank**2 + 1))

        # Function that should become zero
        F = n_m * Q * P ** (2 / 3) - A ** (5 / 3) * m.sqrt(S_0)

        # It's derivate
        dF_dh = ((2/3) * n_m * Q * P ** (-1 / 3)) * dP_dh - ((5 / 3) * A ** (2 / 3) * m.sqrt(S_0)) * dA_dh

        # Water depth Update
        h = abs((h-F)/dF_dh)
        count -= 1

    return [float(h),float(eps)]

if __name__ == '__main__':

    # input parameters
    Q = 15.5        # discharge in (m3/s)
    b = 5.1         # bottom channel width (m)
    m_bank = 2.5    # bank slope
    S_0 = 0.005     # channel slope
    h = 1           # h initial

    Answer = calc_discharge(b, h, m_bank, S_0, k=20)
    print("The result is: %0.2f m3/s" % Answer[0])

    # call the solver with user-defined channel geometry and discharge
    h_n = interpolate_h(Q, b, m_bank, S_0, k=Answer[1])
    print("The result of the interpolation of h is: %0.2f m and the eps is: %0.4f" % (h_n[0],h_n[1]))

    Answer2 = calc_discharge(b, h_n[0], m_bank, S_0, k=Answer[1])
    print("Finally, the result of Q using the interpolated h is %0.2f m3/s" % Answer2[0])
