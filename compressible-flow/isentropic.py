import numpy as np

def calculate(gamma=1.4, **kwargs):
    
    # parse arguments
    input_variable = list(kwargs.keys())[0]
    input_value = np.array(list(kwargs.values())[0])
    
    # constants
    valid_inputs = ["M", "T_T0", "P_P0", "rho_rho0", "A_As_sub", "A_As_sup", "M_angle", "PM_angle"]
    
    # validate input count 
    assert len(kwargs) == 1,\
        f"1 input expected. Input count: {len(kwargs)}"
    
    # validate input names
    assert input_variable in valid_inputs,\
        f"Valid input variable expected. Variable: {input_variable}"
    
    # calculate mach number based on input
    if input_variable == "M": M = mach_calculate_M(input_value)
    elif input_variable == "T_T0": M = mach_calculate_T_T0(input_value, gamma)
    elif input_variable == "P_P0": M = mach_calculate_P_P0(input_value, gamma)
    elif input_variable == "rho_rho0": M = mach_calculate_rho_rho0(input_value, gamma)
    elif input_variable == "A_As_sub": M = mach_calculate_A_As_sub(input_value, gamma)
    elif input_variable == "A_As_sup": M = mach_calculate_A_As_sup(input_value, gamma)
    elif input_variable == "M_angle": M = mach_calculate_M_angle(input_value, gamma)
    elif input_variable == "PM_angle": M = mach_calculate_PM_angle(input_value, gamma)
    
    # calculate other values based on input
    results = {}
    results["M"] = M
    results["T_T0"] = T_T0_calculate(M, gamma)
    results["T_Ts"] = results["T_T0"]/T_T0_calculate(1, gamma)
    results["P_P0"] = P_P0_calculate(M, gamma)
    results["P_Ps"] = results["P_P0"]/P_P0_calculate(1, gamma)
    results["rho_rho0"] = rho_rho0_calculate(M, gamma)
    results["rho_rhos"] = results["rho_rho0"]/rho_rho0_calculate(1, gamma)
    results["A_As_sub"] = A_As_sub_calculate(M, gamma)
    results["A_As_sup"] = A_As_sup_calculate(M, gamma)
    results["M_angle"] = M_angle_calculate(M, gamma)
    results["PM_angle"] = PM_angle_calculate(M, gamma)
    
    return results
    

def mach_calculate_M(input_value):
    
    # value check
    assert (input_value > 0),\
        f"M must be greater than 0"
        
    # assign value
    M = input_value
    return M
    
def mach_calculate_T_T0(input_value, gamma):
    
    # value check
    assert (0 < input_value < 1),\
        f"T_T0 must be between 0 and 1"
    
    # calculate value
    T_T0 = input_value
    M = np.sqrt(2)*np.sqrt((1 - T_T0)/(T_T0*(gamma - 1)))
    return M
    
def mach_calculate_P_P0(input_value, gamma):
    
    # value check
    assert (0 < input_value < 1),\
        f"P_P0 must be between 0 and 1"
    
    # calculate value
    P_P0 = input_value
    M = np.sqrt((2**(gamma/(gamma - 1)) - 2*P_P0*(2**(gamma/(gamma - 1))/P_P0)**(1/gamma))/\
        (P_P0*(2**(gamma/(gamma - 1))/P_P0)**(1/gamma)*(gamma - 1)))
    return M

def mach_calculate_rho_rho0(input_value, gamma):
    
    # value check
    assert (0 < input_value < 1),\
        f"rho_rho0 must be between 0 and 1"
        
    # calculate value
    rho_rho0 = input_value
    M = np.sqrt((-2**(1 + 1/(gamma - 1)) + rho_rho0*(2**(1/(gamma - 1))/rho_rho0)**gamma)/\
        (2**(1/(gamma - 1))*(gamma - 1)))
    return M

def mach_calculate_A_As_sub(input_value, gamma):
    
    # value check
    assert (input_value > 1),\
        f"A_As_sub must be greater than 1"
        
    return

def mach_calculate_A_As_sup(input_value, gamma):
    
    # value check
    assert (input_value > 1),\
        f"A_As_sup must be greater than 1"
    
    return 

def mach_calculate_M_angle(input_value, gamma):
    
    # value check
    assert (0 < input_value < np.pi/2),\
        f"M_angle must be between 0 and pi/2 radians (0 and 90 degrees)"
    
    return

def mach_calculate_PM_angle(input_value, gamma):
    
    # value check
    assert (0 < input_value < np.deg2rad(130.45)),\
        f"PM_angle must be between 0 and 2.277 radians (0 and 130.45 degrees)"
    
    return


def P_P0_calculate(M, gamma):
    return (1 + ((gamma - 1)/2) * (M ** 2)) ** (-gamma/(gamma-1))

def T_T0_calculate(M, gamma):

    return (1 + ((gamma - 1)/2) * (M ** 2)) ** (-1)

def rho_rho0_calculate(M, gamma):

    return (1 + ((gamma - 1)/2) * (M ** 2)) ** (-1/(gamma-1))

def A_As_sub_calculate(M, gamma):

    return 

def A_As_sup_calculate(M, gamma):

    return 

def M_angle_calculate(M, gamma):

    return 

def PM_angle_calculate(M, gamma):

    return 


if __name__ == "__main__":

    # calculate(gamma=1.4, T_T0=0)
    print(calculate(gamma=1.4, T_T0=0.5))
    # calculate(gamma=1.4, T_T0=1)