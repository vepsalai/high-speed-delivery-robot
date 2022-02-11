import numpy as np

#constants
g=9.81
roo=1.225 # air density to be 1.225 kg/m^3
bat_eff=0.95
mot_eff=0.95


def drive_EV_bw_no_reg(t_data, v_data, ele_data, HW, rrc, m_package, aux=60, C = 0.6, A = 0.2, m_bot_bat=15+1):
    '''Backwards model of EV with NO energy regeneration.
    BLDC motors that the d-bot use are not easy to utilze for energy regeneration.'''
    P_data=[0]
    m = m_bot_bat + m_package
    for i in range(1,len(t_data)):
        dx = v_data[i]*(t_data[i]-t_data[i-1])
        dy = ele_data[i]-ele_data[i-1]
        if dx > 0: # When stopped, energy is not consumed to keep at a stand-still in a slope.
            F_grade = m*g*np.sin(np.arctan(dy/dx))
        else:
            F_grade = 0
        if F_grade < 0:
            F_grade = 0
        F_ma = m*(v_data[i]-v_data[i-1])/(t_data[i]-t_data[i-1]) # F=ma
        if v_data[i] > 0:
            F_rr = rrc*m*g # F_rrc=myy*G
        else:
            F_rr = 0
        if (HW+v_data[i]) >= 0:
            F_d = 0.5*roo*C*A*np.power(np.abs(v_data[i]+HW),2) # F_d
        else: 
            F_d = -0.5*roo*C*A*np.power(np.abs(v_data[i]+HW),2) # F_d
        if F_d < 0:
            F_d = 0
        F = F_ma + F_rr + F_d + F_grade
        P = F*v_data[i]
        if P <= 0: # no regenerative braking!
            P=0
        P_with_loss=P/mot_eff
        P_data.append((P_with_loss+aux)/bat_eff)

    energy = np.trapz(P_data, dx=1)/3600
    return energy

def drive_EV_bw_no_reg_max_power(v_max,a, HW, grade, rrc, m_package, aux=60, C = 0.6, A = 0.2, m_bot_bat=15+1):
    '''Backwards model of EV with NO energy regeneration.
    BLDC motors that the d-bot use are not easy to utilze for energy regeneration.'''
    # Create aggressive cycle
    ts=0.1
    t_data=[0]
    v_data=[0]
    v_max=15/3.6
    a=2
    d=0
    v=0
    t=0
    while d < 20:
        if v < v_max:
            v = v + a*ts
        d = d + v*ts
        t=t+ts
        t_data.append(t)
        v_data.append(v)

    P_data=[0]
    e_data=[0] # Initialize as zero, because we need to start from 2nd time step to get acceleration
    F_d_data=[0]
    F_grade_data=[0]
    F_ma_data=[0]
    F_rr_data=[0]
    eff_loss=[0]
    
    m = m_bot_bat + m_package
    for i in range(1,len(t_data)):
        F_grade = m*g*np.sin(grade)
        F_ma = m*(v_data[i]-v_data[i-1])/(t_data[i]-t_data[i-1]) # F=ma
        if v_data[i] > 0:
            F_rr = rrc*m*g # F_rrc=myy*G
        else:
            F_rr = 0
        if (HW+v_data[i]) >= 0:
            F_d = 0.5*roo*C*A*np.power(np.abs(v_data[i]+HW),2) # F_d
        else: 
            F_d = -0.5*roo*C*A*np.power(np.abs(v_data[i]+HW),2) # F_d
        F = F_ma + F_rr + F_d + F_grade
        P = F*v_data[i]
        if P <= 0: # no regenerative braking!
            P=0
        F_d_data.append(F_d)
        F_grade_data.append(F_grade)
        F_ma_data.append(F_ma)
        F_rr_data.append(F_rr)
        P_with_loss=P/mot_eff
        P_data.append((P_with_loss+aux)/bat_eff)
        eff_loss.append(P_data[-1]- (P+aux))


    energy = np.trapz(P_data, dx=0.1)/3600

    return t_data, v_data, P_data, F_d_data, F_grade_data, F_ma_data, F_rr_data, aux, eff_loss, energy

def drive_EV_bw_no_reg_power_dist(t_data, v_data, ele_data, HW, rrc, m_package, aux=60, C = 0.6, A = 0.2, m_bot_bat=15+1):
    '''Backwards model of EV with NO energy regeneration.
    BLDC motors that the d-bot use are not easy to utilze for energy regeneration.'''
    P_data=[0]
    P_grade=[0]
    P_rr=[0]
    P_ma=[0]
    P_d=[0]
    P_aux=[0]
    P_mot_loss=[0]
    P_bat_loss=[0]

    m = m_bot_bat + m_package
    for i in range(1,len(t_data)):
        dx = v_data[i]*(t_data[i]-t_data[i-1])
        dy = ele_data[i]-ele_data[i-1]
        if dx > 0: # When stopped, energy is not consumed to keep at a stand-still in a slope.
            F_grade = m*g*np.sin(np.arctan(dy/dx))
        else:
            F_grade = 0
        if F_grade < 0:
            F_grade = 0
        F_ma = m*(v_data[i]-v_data[i-1])/(t_data[i]-t_data[i-1]) # F=ma
        if v_data[i] > 0:
            F_rr = rrc*m*g # F_rrc=myy*G
        else:
            F_rr = 0
        if (HW+v_data[i]) >= 0:
            F_d = 0.5*roo*C*A*np.power(np.abs(v_data[i]+HW),2) # F_d
        else: 
            F_d = -0.5*roo*C*A*np.power(np.abs(v_data[i]+HW),2) # F_d
        if F_d < 0:
            F_d = 0
        F = F_ma + F_rr + F_d + F_grade
        P = F*v_data[i]
        if P <= 0: # no regenerative braking!
            P=0
        P_with_loss=P/mot_eff
        P_data.append((P_with_loss+aux)/bat_eff)

        P_d.append(F_d*v_data[i])
        P_rr.append(F_rr*v_data[i])
        P_grade.append(F_grade*v_data[i])
        P_ma.append(F_ma*v_data[i])
        P_aux.append(aux)
        P_mot_loss.append(P*(1-mot_eff))
        P_bat_loss.append((P_with_loss+aux)*(1-bat_eff))

    energy = np.trapz(P_data, dx=1)/3600
    e_d = np.trapz(P_d, dx=1)/3600
    e_rr = np.trapz(P_rr, dx=1)/3600
    e_grade = np.trapz(P_grade, dx=1)/3600
    e_ma = np.trapz(P_ma, dx=1)/3600
    e_aux = np.trapz(P_aux, dx=1)/3600
    e_mot_loss = np.trapz(P_mot_loss, dx=1)/3600
    e_bat_loss = np.trapz(P_bat_loss, dx=1)/3600

    return np.array([energy, e_d, e_rr, e_grade, e_ma, e_aux, e_mot_loss, e_bat_loss])

def drive_EV_bw_reg(t_data, v_data, ele_data, HW, rrc, m_package, aux=60, C = 0.6, A = 0.2, m_bot_bat=15+1):
    '''Backwards model of EV with energy regeneration'''
    P_data=[0]
    m = m_bot_bat + m_package
    for i in range(1,len(t_data)):
        dx = v_data[i]*(t_data[i]-t_data[i-1])
        dy = ele_data[i]-ele_data[i-1]
        F_grade = m*g*np.sin(np.arctan(dy/dx))
        F_ma = m*(v_data[i]-v_data[i-1])/(t_data[i]-t_data[i-1]) # F=ma
        F_rr = rrc*m*g # F_rrc=myy*G
        if (HW+v_data[i]) >= 0:
            F_d = 0.5*roo*C*A*np.abs(v_data[i]+HW)**2 # F_d
        else: 
            F_d = -0.5*roo*C*A*np.abs(v_data[i]+HW)**2 # F_d
        F = F_ma + F_rr + F_d + F_grade
        P = F*v_data[i]
        if P <= 0:
            P_with_loss=P*mot_eff
            P_data.append((P_with_loss+aux)*bat_eff)
        else:
            P_with_loss=P/mot_eff
            P_data.append((P_with_loss+aux)/bat_eff)

    energy = np.trapz(P_data, dx=0.1)/3600
    return energy