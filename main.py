#MATEUS TEIXEIRA MAGALHÃƒES - UFJF - PGMC

import numpy as np
import utils


#---------------ANALYTIC SOLUTIONS------------------
def exact1(x, y, t, Re):
    
    return np.exp(-2 * (np.pi ** 2) * t)\
            * np.sin(np.pi * x)\
            * np.sin(np.pi * y)
             
def exact2(x, y, t, Re):
    
    return 2 * np.pi * np.cos(np.pi * x)\
             * np.cos(np.pi * y)\
             * np.exp((-2 * (np.pi ** 2) * t) / Re)
#--------------------------------------------------



#---------------VELOCITY FUNCTIONS-------------------
def vel_x(i, j, k, dt, h, omega, teta, Re):
    
    x = omega[0] + i *h
    y = omega[0] + j *h
    t = teta[0] + k*dt
    
    return -np.cos(np.pi * x)\
             * np.sin(np.pi * y)\
             * np.exp((-2 * (np.pi ** 2) * t) / Re)

def vel_y(i, j, k, dt, h, omega, teta, Re):
    
    x = omega[0] + i *h
    y = omega[0] + j *h
    t = teta[0] + k*dt
    
    return -np.sin(np.pi * x)\
             * np.cos(np.pi * y)\
             * np.exp((-2 * (np.pi ** 2) * t) / Re)             

def vel_0(i, j, k, dt, h, omega, teta, Re):
    
    return 0
#----------------------------------------------------
    

def p(i, j, dt, h, b):
    
    return (b * dt) / h

def solve_x_line(j, k, w_l, omega, teta, n, Re, h, dt,\
                 EXACT_SOLUTION, VELX, VELY):
    
    epsilon = 1/ Re
    
    sigma = (epsilon * dt) / (h**2)
    
    a, b, c, d = [], [], [], []
    for i in range(1, n - 1):
        
        b1 = VELX(i, j, k, dt, h, omega, teta, Re)
        b2 = VELY(i, j, k, dt, h, omega, teta, Re)
        
        p1 = p(i, j, dt, h, b1)
        p2 = p(i, j, dt, h, b2)
        
        b.append(1 + sigma)
        
        d.append( (1 - sigma) * w_l[i, j]\
                 + ((sigma / 2) - (p2 / 4)) * w_l[i, j+1]\
                 + ((sigma / 2) + (p2 / 4)) * w_l[i, j-1])
        
        if( i < n - 2):
            
            c.append( ((p1 / 4) - (sigma / 2)) )
            
        if( i > 1):
            
            a.append( -((p1 / 4) + (sigma / 2)) )
    
    x_line = np.array(utils.TDMAsolver(a, b, c, d))
    return x_line

def solve_y_collum(i, k, w_l, omega, teta, n, Re, h, dt,\
                   EXACT_SOLUTION, VELX, VELY):
    
    epsilon = 1/ Re
    
    sigma = (epsilon * dt) / (h**2)
    
    a, b, c, d = [], [], [], []
    for j in range(1, n - 1):
        
        b1 = VELX(i, j, k, dt, h, omega, teta, Re)
        b2 = VELY(i, j, k, dt, h, omega, teta, Re)
        
        p1 = p(i, j, dt, h, b1)
        p2 = p(i, j, dt, h, b2)
        
        b.append(1 + sigma)
        
        d.append( (1 - sigma) * w_l[i, j]\
                 + ((sigma / 2) - (p1 / 4)) * w_l[i+1, j]\
                 + ((sigma / 2) + (p1 / 4)) * w_l[i-1, j])
        
        if( j < n - 2):
            
            c.append( ((p2 / 4) - (sigma / 2)) )
            
        if( j > 1):
            
            a.append( -((p2 / 4) + (sigma / 2)) )
    
    y_collum = np.array(utils.TDMAsolver(a, b, c, d))
    return y_collum

def solver(omega, teta, n, Re, EXACT_SOLUTION, VELX, VELY):
    
    h = (omega[-1] - omega[0]) / (n - 1)
    
    dt = h
    
    nt = int((teta[-1] - teta[0]) / dt) + 1
    
    x = np.linspace(omega[0], omega[-1], n)
    y = np.linspace(omega[0], omega[-1], n)
    
    X, Y = np.meshgrid(x, y)
    
    initial = lambda x, y: EXACT_SOLUTION(x, y, 0, Re)
    
    #current solution
    w_c = np.array([[initial(x[i], y[j]) for i in range(n)] for j in range(n)])
    
    #last solution
    w_l = np.array([[initial(x[i], y[j]) for i in range(n)] for j in range(n)])
     
#    utils.plot_3d_surface(X, Y, w_l, 1)
    
#    print(w_c)
#    
#    test = np.array([i+1 for i in range(n-2)])
#    
#    print(test)
#    
#    for j in range(1, n - 1):
#        w_c[j, 1:-1] = test[:]
#    print(w_c)
#    
#    for j in range(1, n - 1):
#        w_c[1:-1, j] = test[:]
#    print(w_c)
    
    k = 1
    w_c, w_l = w_l, w_c
    for k in range(nt):
        
        w_c, w_l = w_l, w_c
        for j in range(1, n - 1):
        
            w_c[1:-1, j] = solve_x_line(j, k, w_l, omega, teta, n, Re, h, dt,\
                          EXACT_SOLUTION, VELX, VELY)[:]
            
        w_c, w_l = w_l, w_c
        for i in range(1, n - 1):
        
            w_c[i, 1:-1] = solve_y_collum(i, k, w_l, omega, teta, n, Re, h, dt,\
                          EXACT_SOLUTION, VELX, VELY)[:]

        exact = np.array([[EXACT_SOLUTION(x[i], y[j], (k+1)*dt, Re)\
                           for i in range(n)]\
                                for j in range(n)])
        
#        utils.plot_3d_numeric_vs_analytic(X, Y, w_c, exact, 'solution/solution.png')
        
#        print("max error: ", utils.max_error(w_c, exact, n, n))
        
    return utils.max_error(w_c, exact, n, n), h
        
def main():
    
    omega = [0, 1]
    
    teta = [0, 0.5]
    
    Re = 1
    
    n_vec = [4, 8, 16, 32, 64]
    
    error_vec = []
    
    h_vec = []
    
    for i in range(len(n_vec)):
        
        n = n_vec[i]
        
        error, h = solver(omega, teta, n, Re, exact1, vel_0, vel_0)
        
        error_vec.append(error)
        
        h_vec.append(h)
        
        if(len(error_vec) > 1):
            
            print("taxa de convergência: (", n_vec[i-1], "->", n_vec[i] ,")",\
                  (np.log(error_vec[-1]) - np.log(error_vec[-2]))\
                  / (np.log(h_vec[-1]) - np.log(h_vec[-2])))
    
    utils.plt.plot(-np.log(h_vec), np.log(error_vec), marker = '.')
    
main()