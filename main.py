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
    
    return np.sin(np.pi * x)\
             * np.cos(np.pi * y)\
             * np.exp((-2 * (np.pi ** 2) * t) / Re)             

def vel_0(i, j, k, dt, h, omega, teta, Re):
    
    return 0
#----------------------------------------------------
    

def p(i, j, dt, h, b):
    
    return (b * dt) / h

def solve_x_line(j, k, w_l, omega, teta, n, Re, h, dt,\
                 EXACT_SOLUTION, VELX, VELY):
    
    epsilon = 1 / Re
    
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
    
    #countour
    d[0] = d[0] - ((-((p1 / 4) + (sigma / 2)))*w_l[0, j])
    d[-1] = d[-1] - (((p1 / 4) - (sigma / 2))*w_l[n-1, j])
    
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
    
    #countour
    d[0] = d[0] - ((-((p2 / 4) + (sigma / 2)))*w_l[i, 0])
    d[-1] = d[-1] - (((p2 / 4) - (sigma / 2))*w_l[i, n-1])
    
    y_collum = np.array(utils.TDMAsolver(a, b, c, d))
    return y_collum

def boudary_condition(omega, teta, k, h, dt, EXACT_SOLUTION, n, Re,\
                      axes = 'x', position = 'start'):
    
    if(axes == 'x'):
        
        if(position == 'start'):
            
            return np.array([EXACT_SOLUTION(omega[0],\
                                            omega[0] + j * h,\
                                            teta[0] + k * dt,\
                                            Re) for j in range(n)])
            
        if(position == 'end'):
        
            return np.array([EXACT_SOLUTION(omega[-1],\
                                            omega[0] + j * h,\
                                            teta[0] + k * dt,\
                                            Re) for j in range(n)])
            
    if(axes == 'y'):
        
        if(position == 'start'):
            
            return np.array([EXACT_SOLUTION(omega[0] + i * h,\
                                            omega[0],\
                                            teta[0] + k * dt,\
                                            Re) for i in range(n)])
            
        if(position == 'end'):  
            
            return np.array([EXACT_SOLUTION(omega[0] + i * h,\
                                            omega[-1],\
                                            teta[0] + k * dt,\
                                            Re) for i in range(n)])

def solver(omega, teta, n, Re, H_FUNC, EXACT_SOLUTION, VELX, VELY,\
           plot_in_time = False, outpu_folder = 'numeric'):
    
    h = (omega[-1] - omega[0]) / (n - 1)
    
    dt = H_FUNC(h)
    
    nt = int((teta[-1] - teta[0]) / dt)
    
    print('n: ' + str(n))
    print('h: ' + str(h))
    print('nt: ' + str(nt))
    print('\n\n\n')
    
    x = np.linspace(omega[0], omega[-1], n)
    y = np.linspace(omega[0], omega[-1], n)
    
    X, Y = np.meshgrid(x, y)
    
    initial = lambda x, y: EXACT_SOLUTION(x, y, 0, Re)
    
    #current solution
    w_c = np.array([[initial(x[i], y[j]) for i in range(n)]\
                     for j in range(n)])
    
    #last solution
    w_l = np.array([[initial(x[i], y[j]) for i in range(n)]\
                     for j in range(n)])
    
    min_z, max_z = np.min(w_c), np.max(w_c)
    
    #plot initial solution
    if(plot_in_time):
            
            utils.plot_3d_surface(X, Y, w_c, 1,\
                                  name = outpu_folder\
                                  + '/numeric_n=' + str(n)\
                                  + '_t=' + str(0)\
                                  + '.png',\
                                  min_z = min_z,\
                                  max_z = max_z)
    
    w_c, w_l = w_l, w_c
    for k in range(nt):
        
        w_c, w_l = w_l, w_c
        for j in range(1, n - 1):
        
            w_c[1:-1, j] = solve_x_line(j, k + 1, w_l, omega,\
                          teta, n, Re, h, dt,\
                          EXACT_SOLUTION, VELX, VELY)[:]
            
        w_c, w_l = w_l, w_c
        for i in range(1, n - 1):
        
            w_c[i, 1:-1] = solve_y_collum(i, k + 1, w_l, omega,\
                          teta, n, Re, h, dt,\
                          EXACT_SOLUTION, VELX, VELY)[:]
        
        w_c[0, :] = boudary_condition(omega, teta, k + 1, h, dt, EXACT_SOLUTION, n, Re,\
                      axes = 'x', position = 'start')[:]
        
        w_c[-1, :] = boudary_condition(omega, teta, k + 1, h, dt, EXACT_SOLUTION, n, Re,\
                      axes = 'x', position = 'end')[:]
        
        w_c[:, 0] = boudary_condition(omega, teta, k + 1, h, dt, EXACT_SOLUTION, n, Re,\
                      axes = 'y', position = 'start')[:]
        
        w_c[:, -1] = boudary_condition(omega, teta, k + 1, h, dt, EXACT_SOLUTION, n, Re,\
                      axes = 'y', position = 'end')[:]
            
        if(plot_in_time):
            
            utils.plot_3d_surface(X, Y, w_c, 1,\
                                  name = outpu_folder\
                                  + '/numeric_n=' + str(n)\
                                  + '_t=' + str(k+1)\
                                  + '.png',\
                                  min_z = min_z,\
                                  max_z = max_z)
    
    #exact solution
    exact = np.array([[EXACT_SOLUTION(x[i], y[j], teta[-1], Re)\
                       for i in range(n)]\
                            for j in range(n)])
        
    return utils.max_error(w_c, exact, n, n), h
        
def gera_taxas_de_convertgencia(omega = [-1, 1], teta = [0, 1], Re = 20,\
                                H_FUNC = lambda h: h, EXACT_SOLUTION = exact2):
    
    n_vec = [4, 8, 16, 32, 64]
    
    error_vec = []
    
    h_vec = []
    
    for i in range(len(n_vec)):
        
        n = n_vec[i]
        
        error, h = solver(omega, teta, n, Re, H_FUNC,\
                          EXACT_SOLUTION, vel_0, vel_0)
        
        error_vec.append(error)
        
        h_vec.append(h)
        
        if(len(error_vec) > 1):
            
            print("taxa de convergência: (", n_vec[i-1], "->", n_vec[i] ,")",\
                  (np.log(error_vec[-1]) - np.log(error_vec[-2]))\
                  / (np.log(h_vec[-1]) - np.log(h_vec[-2])))
    
    utils.plt.plot(-np.log(h_vec), np.log(error_vec),\
                   marker = '.')

def teste_numerica(n = 32, omega = [-1, 1], teta = [0, 1],\
                    EXACT_SOLUTION = exact2, Re = 20, H_FUNC = lambda h: h,\
                    VELX = vel_x, VELY = vel_y):
        
    erro, h = solver(omega, teta, n, Re, H_FUNC, EXACT_SOLUTION, VELX, VELY,\
               plot_in_time= True, outpu_folder='teste_numerica')
    
def teste_analitica(n = 32, omega = [-1, 1], teta = [0, 1],\
                    EXACT_SOLUTION = exact2, Re = 20, H_FUNC = lambda h: h):
    
    h = (omega[-1] - omega[0]) / (n - 1)
        
    dt = H_FUNC(h)
    
    nt = int((teta[-1] - teta[0]) / dt)
    
    #para gerar um plot 3d
    X = [omega[0] + i*h for i in range(0, n)]
    Y = [omega[0] + j*h for j in range(0, n)]
    
    X, Y = np.meshgrid(X, Y)
    
    min_z, max_z = 0, 0
    
    for k in range(nt+1):
        
        print(teta[0] + k * dt)
        
        Z = np.array([[EXACT_SOLUTION(omega[0] + i * h, omega[0] + j * h, teta[0] + k * dt, Re)\
              for i in range(n)] for j in range(n)])
    
        if(k == 0):
            
           min_z, max_z = np.min(Z), np.max(Z)
    
        utils.plot_3d_surface(X, Y, np.array(Z), 1,\
                              name = 'teste_analitica/solution' + str(k) + '.png',\
                              min_z = min_z,\
                              max_z = max_z)
 
gera_taxas_de_convertgencia(omega = [0, 1], teta = [0, 0.5], Re = 1,\
                                H_FUNC = lambda h: h**2, EXACT_SOLUTION = exact1)

#teste_analitica(EXACT_SOLUTION=exact1, omega=[0, 1], teta= [0, 0.5], Re=1)
#
#teste_numerica(EXACT_SOLUTION=exact1, omega=[0, 1], teta= [0, 0.5], Re=1,\
#               VELX=vel_0, VELY=vel_0)