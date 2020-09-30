#MATEUS TEIXEIRA MAGALHÃƒES - UFJF - PGMC

import numpy as np
import os
import shutil

#solver for tri-diagonal matrix
def TDMAsolver(a, b, c, d):
    
    nf = len(d) # number of equations
    
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    
    for it in range(1, nf):
        
        mc = ac[it-1]/bc[it-1]
        
        bc[it] = bc[it] - mc*cc[it-1] 
        
        dc[it] = dc[it] - mc*dc[it-1]
        	    
    xc = bc
    
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc

def save_vtk(nome, dado, nx, ny, delete_folder=False, data_name='data'):
    
    print(nome)
    if not os.path.exists(nome.split('/')[0]):
        os.makedirs(nome.split('/')[0])
    if delete_folder:
        shutil.rmtree(nome.split('/')[0], ignore_errors=True)
        os.makedirs(nome.split('/')[0])
    
    arqvtk = open(nome, "w")
    texto = []
    
    texto.append("# vtk DataFile Version 3.0\n")
    texto.append("vtk output\n")
    texto.append("ASCII\n")
    texto.append("DATASET RECTILINEAR_GRID\n")
    texto.append("DIMENSIONS " + str(nx) + " " + str(ny) + " 1\n")
    
    texto.append("X_COORDINATES " + str(nx) + " double\n")
    for ix in range(nx):
        texto.append(str(ix) + " ")
    texto.append("\n")
    
    texto.append("Y_COORDINATES " + str(ny) + " double\n")
    for jy in range(ny):
        texto.append(str(jy) + " ")
    texto.append("\n")
    
    texto.append("Z_COORDINATES 1 double\n")
    texto.append("0\n")

    texto.append("POINT_DATA " + str((nx)*(ny)) + "\n")
    texto.append("FIELD FieldData 1\n")
    texto.append("" + data_name + " 1 " + str((nx)*(ny)) + " double\n")
    for jy in range(ny):
        for ix in range(nx):
                texto.append(str(dado[ix, jy]) + " ")
    texto.append("\n")
    
    arqvtk.writelines(texto)
    arqvtk.close()

def residuo_and_psi_w_max(w, psi, n, h, omega, Re, p, r):
    
    residuo_psi = np.zeros((n - 2, n - 2))
    residuo_w = np.zeros((n - 2, n - 2))
    max_psi = 0
    max_w = 0
    i_max, j_max = 0, 0
    
    for i in range(1, n-1):
        
        for j in range(1, n-1):
            
            residuo_psi[i-1, j-1] = abs(((psi[r][i+1, j] - (2 * psi[r][i,j]) + psi[r][i-1, j])/(h**2))\
                              + ((psi[r][i, j+1] - (2 * psi[r][i,j]) + psi[r][i, j-1])/(h**2))\
                              + w[p][i,j])
                              
            residuo_w[i-1, j-1] = abs(((1/Re)*((w[p][i, j+1] - (2 * w[p][i,j]) + w[p][i, j-1])/(h**2)))\
                              + ((1/Re)*((w[p][i+1, j] - (2 * w[p][i,j]) + w[p][i-1, j])/(h**2)))\
                              - (((psi[r][i, j+1] - psi[r][i, j-1])/(2*h))*((w[p][i+1, j] - w[p][i-1, j])/(2*h)))\
                              + (((psi[r][i+1, j] - psi[p][i-1, j])/(2*h))*((w[p][i, j+1] - w[p][i, j-1])/(2*h))))
    
            if psi[r][i, j] > max_psi:
                
                max_psi = abs(psi[r][i, j])
                i_max, j_max = i, j
                
            if w[r][i, j] > max_w:
                
                max_w = abs(w[p][i, j])
                i_max, j_max = i, j
    
    max_residuo_psi = np.max(residuo_psi)
    max_residuo_w = np.max(residuo_w)
    max_x = omega[0] + (i_max * h)
    max_y = omega[0] + (j_max * h)
    
    return max_psi, max_w, max_x, max_y, max(max_residuo_psi, max_residuo_w)
    
def solve_x_line(j, w, psi, omega, n, Re, h, dt, sigma, p, q, r, s, variable = 'w'):    
    
    a, b, c, d = [], [], [], []
    for i in range(1, n - 1):

        if variable == 'psi':
            
            b.append(1 + (2*sigma))
            
            c.append(-sigma)
            
            a.append(-sigma)
            
            d.append( ((dt/2)*w[q][i,j])\
                     + ((1-(2*sigma))*psi[s][i,j])\
                     + (sigma*(psi[s][i, j+1] + psi[s][i, j-1])))
            
        if variable == 'w':
            
            px = Re * ( (psi[r][i+1, j] - psi[r][i-1, j]) / 4 )
            py = Re * ( (psi[r][i, j+1] - psi[r][i, j-1]) / 4 )
            
            b.append(1 + (2*sigma))
            
            c.append(sigma * (py - 1))
            
            a.append(-(sigma * (py + 1)))

            d.append( ((1- (2*sigma))*w[q][i,j])\
                     + ((sigma*(px+1))*w[q][i, j+1])\
                     + (sigma*(1-px)*w[q][i, j-1]))
            
    
    if variable == 'psi':
    
        #contorno
        d[0] = d[0] + (sigma * psi[r][0, j])
        d[-1] = d[-1] + (sigma * psi[r][-1, j])     
        
    if variable == 'w':
    
        #contorno
        d[0] = d[0] + ((sigma * (py + 1)) *  w[p][0, j])
        d[-1] = d[-1] - (sigma * (py - 1) * w[p][-1, j])
    
    x_line = np.array(TDMAsolver(a, b, c, d))
    return x_line

def solve_y_collum(i, w, psi, omega, n, Re, h, dt, sigma, p, q, r, s, variable = 'w'):
    
    a, b, c, d = [], [], [], []
    for j in range(1, n - 1):
        
        if variable == 'psi':
            
            b.append(1 + (2*sigma))
            
            c.append(-sigma)
            
            a.append(-sigma)
            
            d.append( ((dt/2)*w[q][i,j])\
                     + ((1-(2*sigma))*psi[s][i,j])\
                     + (sigma*(psi[s][i+1, j] + psi[s][i-1, j])))
            
        if variable == 'w':
            
            px = Re * ( (psi[r][i+1, j] - psi[r][i-1, j]) / 4 )
            py = Re * ( (psi[r][i, j+1] - psi[r][i, j-1]) / 4 )
            
            b.append(1 + (2*sigma))
            
            c.append(-(sigma * (px + 1)))
            
            a.append(sigma * (px - 1))

            d.append( ((1- (2*sigma))*w[q][i,j])\
                     + ((sigma*(1-py))*w[q][i+1, j])\
                     + ((sigma*(1+py))*w[q][i-1, j]))
                
    
    if variable == 'psi':
    
        #contorno
        d[0] = d[0] + (sigma * psi[r][i, 0])
        d[-1] = d[-1] + (sigma * psi[r][j, -1])     
        
    if variable == 'w':
    
        #contorno
        d[0] = d[0] - ((sigma * (px - 1)) *  w[p][i, 0])
        d[-1] = d[-1] + ((sigma * (px + 1)) * w[p][i, -1])
    
    y_collum = np.array(TDMAsolver(a, b, c, d))
    return y_collum


def solver(omega, n, Re, H_FUNC = lambda h: h**2, tol = 1e-6, residuo_it = 1, max_it = 100):
    
    h = (omega[-1] - omega[0]) / (n - 1)
    
    dt = H_FUNC(h)
        
    sigma = dt / (2*(h**2))
    
    print('n: ' + str(n))
    print('h: ' + str(h))
    print('\n\n\n')
        
    #solucao inicial para vorticidade
    w = [np.zeros((n, n)), np.zeros((n, n))]
        
    #solucao inicial para funcao corrente
    psi = [np.zeros((n, n)), np.zeros((n, n))]

#    #contorno em t = 0
#    w[0][:, -1] =  - (2/h)
#    w[1][:, -1] =  - (2/h)

    p, q, = 0, 1
    r, s, = 0, 1
    error = 1
    it = 0
    while error > tol and it <= max_it:
        
        r, s = s, r
        for j in range(1, n - 1):
        
            psi[r][1:-1, j] = solve_x_line(j, w, psi, omega, n, Re, h, dt,\
                                           sigma, p, q, r, s, variable = 'psi')[:]
            
        r, s = s, r
        for i in range(1, n - 1):
        
            psi[r][i, 1:-1] = solve_y_collum(i, w, psi, omega, n, Re, h, dt,\
                                           sigma, p, q, r, s, variable = 'psi')[:]
        
        #aplicando contorno i, 0
        w[p][:, 0] = -2*(psi[r][:, 1]/(h**2))
        
        #aplicando contorno I, j
        w[p][-1, :] = -2*(psi[r][-2, :]/(h**2))
        
        #aplicando contorno 0, j
        w[p][0, :] = -2*(psi[r][1, :]/(h**2))
        
        #aplicando contorno i, J
        w[p][:, -1] = -2*(psi[r][:, -2]/(h**2)) - (2/h)
        
        p, q = q, p
        for j in range(1, n - 1):
        
            w[p][1:-1, j] = solve_x_line(j, w, psi, omega, n, Re, h, dt,\
                                           sigma, p, q, r, s, variable = 'w')[:]
            
        p, q = q, p
        for i in range(1, n - 1):
        
            w[p][i, 1:-1] = solve_y_collum(i, w, psi, omega, n, Re, h, dt,\
                                           sigma, p, q, r, s, variable = 'w')[:]
        
        if it % residuo_it == 0:
            
            max_psi, max_w, x_max, y_max, erro = residuo_and_psi_w_max(w, psi, n, h, omega, Re, p, r)
            print("---------------INTERAÇÃO ", it, "---------------", \
                  "\n\nresiduo: ", erro,\
                  "\nmax psi: ", max_psi,\
                  "\nmax w: ", max_w,\
                  "\nx_max: ", x_max,\
                  "\ny_max: ", y_max)
            
            save_vtk('resultados_cavidade/psi' + str(it) + '.vtk', psi[r], n, n,\
                     delete_folder=it==0, data_name='psi')
            save_vtk('resultados_cavidade/w' + str(it) + '.vtk', w[p], n, n,\
                     delete_folder=it==0, data_name='w')
            
            print("\n---------------------------------------------\n\n\n\n\n");
            
        it = it + 1

def main():

    omega = [0, 1]
    
    n = 32
    
    Re = 1000
    
    H_FUNC = lambda h: ((h**2) / 10)
    
    tol = 1e-6
    
    residuo_it = 25
    
    max_it = 5000
    
    solver(omega, n, Re, H_FUNC = H_FUNC, tol = tol, residuo_it = residuo_it, max_it=max_it)

main()
            
            