# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import numpy as np
from numpy import sqrt, exp
import matplotlib.pyplot as mp
from mpl_toolkits.mplot3d import Axes3D

# +
# 핵심, 1개만 보는것1

a_real = 1.45
b_real = 1.46
c_real = a_real
t1, t2, t3 =1, 0.5, 1.5
kx = np.linspace(0, np.pi, 100)
ky = np.linspace(0, np.pi, 100)
KX, KY = np.meshgrid(kx, ky)

a = [0] #a[0]
a.append(np.array([sqrt(3)/2,-1/2]))  #a[1]
a.append(np.array([0,-1])) #a[2]
a.append(np.array([-sqrt(3)/2,-1/2]))           #a[3]
a = [element * a_real for element in a]
b = [0]
b.append(np.array([a_real,0]))
b.append(np.array([0,b_real]))


n=6
eigenvalues = np.zeros((n, len(kx), len(ky)), dtype=complex)

for i in range(len(kx)):
    for j in range(len(ky)):
        H = np.zeros([n,n], dtype=complex)
        H_t1 = t1*np.diag(np.concatenate(( [exp(1j*(kx[i]*a[c][0]+ky[j]*a[c][1])) for c in [1,2,3]], 
                              [exp(-1j*(kx[i]*a[c][0]+ky[j]*a[c][1])) for c in [1,2]]
                               )), k=1)

        H_t2 = t2*exp(1j*(kx[i]*b[1][0]+ky[j]*b[1][1]))
        H_t3 = t3*exp(1j*(kx[i]*b[2][0]+ky[j]*b[2][1]))

        H[0][3], H[0][5], H[1][5], H[2][4] = H_t3, exp(1j*(kx[i]*a[3][0]+ky[j]*a[3][1])), H_t2, H_t2

        H += H_t1
        H += H.conj().T

        eigenvalues[:, i, j],_ = np.linalg.eigh(H)

fig = mp.figure(figsize=(8,8))

ax = fig.add_subplot(111,projection = '3d')
ax.plot_surface(KX, KY, eigenvalues[2, :, :], cmap= 'cool')
ax.plot_surface(KX, KY, eigenvalues[3, :, :], cmap='spring')
ax.set_xlabel('K_x')
ax.set_ylabel('K_y')
ax.set_zlabel('E')
ax.view_init(30,60)
mp.title(f'$Biphenylene\;t_1 = {t1}, t_2 = {t2}, t_3 = {t3}$')
mp.show()

# -

# %whos

# +
#kx 변화 그림

a_real = 1.45
b_real = 1.46
c_real = a_real
t1, t2, t3 = 1, 0.5, 1.5
kx = np.linspace(0, np.pi, 100)
ky = -0.010577753042389837

a = [0] #a[0]
a.append(np.array([sqrt(3)/2,-1/2]))  #a[1]
a.append(np.array([0,-1])) #a[2]
a.append(np.array([-sqrt(3)/2,-1/2]))           #a[3]
a = [element * a_real for element in a]
b = [0]
b.append(np.array([a_real,0]))
b.append(np.array([0,b_real]))


n=6
eigenvalues = np.zeros((n, len(kx)), dtype=complex)

for i in range(len(kx)):
        H = np.zeros([n,n], dtype=complex)
        H_t1 = t1*np.diag(np.concatenate(( [exp(1j*(kx[i]*a[c][0]+ky*a[c][1])) for c in [1,2,3]], 
                              [exp(-1j*(kx[i]*a[c][0]+ky*a[c][1])) for c in [1,2]]
                               )), k=1)

        H_t2 = t2*exp(1j*(kx[i]*b[1][0]+ky*b[1][1]))
        H_t3 = t3*exp(1j*(kx[i]*b[2][0]+ky*b[2][1]))

        H[0][3], H[0][5], H[1][5], H[2][4] = H_t3, exp(1j*(kx[i]*a[3][0]+ky*a[3][1])), H_t2, H_t2

        H += H_t1
        H += H.conj().T

        eigenvalues[:, i],_ = np.linalg.eigh(H)


fig = mp.figure(figsize=(8,8))

ax = fig.add_subplot(111)
ax.plot(eigenvalues[0, :])
ax.plot(eigenvalues[1, :])
ax.plot(eigenvalues[2, :])
ax.plot(eigenvalues[3, :])
ax.plot(eigenvalues[4, :])
ax.plot(eigenvalues[5, :])
ax.set_xlabel('K_x')
ax.set_ylabel('E')
mp.title(f'$t_1 {t1:.2f}, t_2 {t2}, t_3 {t3} , ky = {ky:.3f}$')
mp.show()


# +
#ky 변화하면서 그림

a_real = 1.45
b_real = 1.46
c_real = a_real
t1, t2, t3 = 1, 1, 2
ky = np.linspace(-np.pi, np.pi, 100)/3
kx = 0

a = [0] #a[0]
a.append(np.array([sqrt(3)/2,-1/2]))  #a[1]
a.append(np.array([0,-1])) #a[2]
a.append(np.array([-sqrt(3)/2,-1/2]))           #a[3]
a = [element * a_real for element in a]
b = [0]
b.append(np.array([(1+sqrt(3))*a_real,0]))
b.append(np.array([0,2*a_real+b_real]))


n=6
eigenvalues = np.zeros((n, len(ky)), dtype=complex)

for i in range(len(ky)):
        H = np.zeros([n,n], dtype=complex)
        H_t1 = t1*np.diag(np.concatenate(( [exp(1j*(kx*a[c][0]+ky[i]*a[c][1])) for c in [1,2,3]], 
                              [exp(-1j*(kx*a[c][0]+ky[i]*a[c][1])) for c in [1,2]]
                               )), k=1)

        H_t2 = t2*exp(1j*(kx*b[1][0]+ky[i]*b[1][1]))
        H_t3 = t3*exp(1j*(kx*b[2][0]+ky[i]*b[2][1]))

        H[0][3], H[0][5], H[1][5], H[2][4] = H_t3, exp(1j*(kx*a[3][0]+ky[i]*a[3][1])), H_t2, H_t2

        H += H_t1
        H += H.conj().T

        eigenvalues[:, i],_ = np.linalg.eigh(H)


fig = mp.figure(figsize=(8,8))

ax = fig.add_subplot(111)
ax.plot(eigenvalues[2, :])
ax.plot(eigenvalues[3, :])
ax.set_xlabel('K_y')
ax.set_ylabel('E')
mp.title('Biphenylene')
mp.show()


# +
#한번에 다그림

a_real = 1.45
b_real = 1.46
c_real = a_real
t1s = np.linspace(0.5, 2, 5)
for t1 in t1s:
    t2s = [0.5, 1, 1.5, 2]
    t3s = [0.5, 1, 1.5, 2]
    kx = np.linspace(-np.pi, np.pi, 100)
    ky = np.linspace(-np.pi, np.pi, 100)
    KX, KY = np.meshgrid(kx, ky)

    a = [0] #a[0]
    a.append(np.array([sqrt(3)/2,-1/2]))  #a[1]
    a.append(np.array([0,-1])) #a[2]
    a.append(np.array([-sqrt(3)/2,-1/2]))           #a[3]
    a = [element * a_real for element in a]
    b = [0]
    b.append(np.array([a_real,0]))
    b.append(np.array([0,-1*b_real]))

    n=6
    eigenvalues = np.zeros((n, len(kx), len(ky)), dtype=complex)
    ts = 1
    fig = mp.figure(figsize=(20,20))
    for t2 in t2s:
        for t3 in t3s:
            for i in range(len(kx)):
                for j in range(len(ky)):
                    H = np.zeros([n,n], dtype=complex)
                    H_t1 = t1*np.diag(np.concatenate(( [exp(1j*(kx[i]*a[c][0]+ky[j]*a[c][1])) for c in [1,2,3]], 
                                          [exp(-1j*(kx[i]*a[c][0]+ky[j]*a[c][1])) for c in [1,2]]
                                           )), k=1)

                    H_t2 = t2*exp(1j*(kx[i]*b[1][0]+ky[j]*b[1][1]))
                    H_t3 = t3*exp(1j*(kx[i]*b[2][0]+ky[j]*b[2][1]))

                    H[0][3], H[0][5], H[1][5], H[2][4] = H_t3, exp(1j*(kx[i]*a[3][0]+ky[j]*a[3][1])), H_t2, H_t2

                    H += H_t1
                    H += H.conj().T

                    eigenvalues[:, i, j],_ = np.linalg.eigh(H)

            ax = fig.add_subplot(len(t2s), len(t2s), ts, projection = '3d')
            ax.plot_surface(KX, KY, eigenvalues[0, :, :], cmap= 'cool')
            ax.plot_surface(KX, KY, eigenvalues[1, :, :], cmap='spring')
            ax.plot_surface(KX, KY, eigenvalues[2, :, :], cmap= 'cool')
            ax.plot_surface(KX, KY, eigenvalues[3, :, :], cmap='spring')
            ax.plot_surface(KX, KY, eigenvalues[4, :, :], cmap= 'cool')
            ax.plot_surface(KX, KY, eigenvalues[5, :, :], cmap='spring')
            ax.set_xlabel('K_x')
            ax.set_ylabel('K_y')a
            ax.set_zlabel('E')
            ax.view_init(30,60)
            mp.title(f'$t_1 = {t1}, t_2 = {t2}, t_3 = {t3}$')
            ts += 1

    mp.show()

# -

## # 밴드끼리 만나는 점 찾기
a_real = 1.45
b_real = 1.46
c_real = a_real
t1s = np.linspace(0.5, 2, 10)
break_outer_loop = False
ax_index = 0
for t1 in t1s:
    print(f't_1 = {t1:.2f}')
    t2s = [0.5, 1, 1.5, 2]
    t3s = [0.5, 1, 1.5, 2]
    kx = np.linspace(-np.pi, np.pi, 100)
    ky = np.linspace(-np.pi, np.pi, 100)
    KX, KY = np.meshgrid(kx, ky)

    a = [0] #a[0]
    a.append(np.array([sqrt(3)/2,-1/2]))  #a[1]
    a.append(np.array([0,-1])) #a[2]
    a.append(np.array([-sqrt(3)/2,-1/2]))           #a[3]
    a = [element * a_real for element in a]
    b = [0]
    b.append(np.array([a_real,0]))
    b.append(np.array([0,-b_real]))

    n=6
    eigenvalues = np.zeros((n, len(kx), len(ky)), dtype=complex)
    ts = 1
    threshold = 0.005
    for t2 in t2s:
        for t3 in t3s:
            for i in range(len(kx)):
                for j in range(len(ky)):
                    H = np.zeros([n,n], dtype=complex)
                    H_t1 = t1*np.diag(np.concatenate(( [exp(1j*(kx[i]*a[c][0]+ky[j]*a[c][1])) for c in [1,2,3]], 
                                          [exp(-1j*(kx[i]*a[c][0]+ky[j]*a[c][1])) for c in [1,2]]
                                           )), k=1)

                    H_t2 = t2*exp(1j*(kx[i]*b[1][0]+ky[j]*b[1][1]))
                    H_t3 = t3*exp(1j*(kx[i]*b[2][0]+ky[j]*b[2][1]))

                    H[0][3], H[0][5], H[1][5], H[2][4] = H_t3, exp(1j*(kx[i]*a[3][0]+ky[j]*a[3][1])), H_t2, H_t2

                    H += H_t1
                    H += H.conj().T

                    eigenvalues[:, i, j],_ = np.linalg.eigh(H)
            for i in range(eigenvalues[0, :, :].shape[0]):
                for j in range(eigenvalues[0, :, :].shape[1]):
                    diff1 = abs(eigenvalues[0, :, :][i,j] - eigenvalues[1, :, :][i, j])
                    diff2 = abs(eigenvalues[2, :, :][i,j] - eigenvalues[3, :, :][i, j])
                    diff3 = abs(eigenvalues[4, :, :][i,j] - eigenvalues[5, :, :][i, j])
                    diff4 = abs(eigenvalues[1, :, :][i,j] - eigenvalues[2, :, :][i, j])
                    diff5 = abs(eigenvalues[3, :, :][i,j] - eigenvalues[4, :, :][i, j])

                    if diff1 <= threshold:
                        print(f"t1 {t1:.2f}, t2 {t2}: Difference = {diff1}, kx = {kx[i]}, ky = {ky[j]}")
                        fig = mp.figure(figsize=(6,6))
                        ax = fig.add_subplot(111 , projection = '3d')
                        ax.plot_surface(KX, KY, eigenvalues[0, :, :], cmap= 'cool')
                        ax.plot_surface(KX, KY, eigenvalues[1, :, :], cmap='spring')
                        ax.set_xlabel('K_x')
                        ax.set_ylabel('K_y')
                        ax.set_zlabel('E')
                        ax.view_init(30,60)
                        mp.title(f'$t_1 = {t1:.2f}, t_2 = {t2}, t_3 = {t3}, 1,2 band $')
                        mp.show()
                        break_outer_loop = True
                        break
                    if diff2 <= threshold:
                        print(f"t1 {t1:.2f}, t2 {t2}: Difference = {diff2}, kx = {kx[i]}, ky = {ky[j]}")
                        fig = mp.figure(figsize=(6,6))
                        ax = fig.add_subplot(111, projection = '3d')
                        ax.plot_surface(KX, KY, eigenvalues[2, :, :], cmap= 'cool')
                        ax.plot_surface(KX, KY, eigenvalues[3, :, :], cmap='spring')
                        ax.set_xlabel('K_x')
                        ax.set_ylabel('K_y')
                        ax.set_zlabel('E')
                        ax.view_init(30,60)
                        mp.title(f'$t_1 = {t1:.2f}, t_2 = {t2}, t_3 = {t3}, 3,4 band $')
                        mp.show()
                        break_outer_loop = True
                        break
                    if diff3 <= threshold:
                        print(f"t1 {t1:.2f}, t2 {t2}: Difference = {diff3}, kx = {kx[i]}, ky = {ky[j]}")
                        fig = mp.figure(figsize=(6,6))
                        ax = fig.add_subplot(111, projection = '3d')
                        ax.plot_surface(KX, KY, eigenvalues[4, :, :], cmap= 'cool')
                        ax.plot_surface(KX, KY, eigenvalues[5, :, :], cmap='spring')
                        ax.set_xlabel('K_x')
                        ax.set_ylabel('K_y')
                        ax.set_zlabel('E')
                        ax.view_init(30,60)
                        mp.title(f'$t_1 = {t1:.2f}, t_2 = {t2}, t_3 = {t3}, 5,6 band $')
                        mp.show()
                        break_outer_loop = True
                        break
                    if diff4 <= threshold:
                        print(f"t1 {t1:.2f}, t2 {t2}: Difference = {diff4}, kx = {kx[i]}, ky = {ky[j]}")
                        fig = mp.figure(figsize=(6,6))
                        ax = fig.add_subplot(111, projection = '3d')
                        ax.plot_surface(KX, KY, eigenvalues[1, :, :], cmap= 'cool')
                        ax.plot_surface(KX, KY, eigenvalues[2, :, :], cmap='spring')
                        ax.set_xlabel('K_x')
                        ax.set_ylabel('K_y')
                        ax.set_zlabel('E')
                        ax.view_init(30,60)
                        mp.title(f'$t_1 = {t1:.2f}, t_2 = {t2}, t_3 = {t3}, 2,3 band $')
                        mp.show()
                        break_outer_loop = True
                        break
                    if diff5 <= threshold:
                        print(f"t1 {t1:.2f}, t2 {t2}: Difference = {diff5}, kx = {kx[i]}, ky = {ky[j]}")
                        fig = mp.figure(figsize=(6,6))
                        ax = fig.add_subplot(111, projection = '3d')
                        ax.plot_surface(KX, KY, eigenvalues[3, :, :], cmap= 'cool')
                        ax.plot_surface(KX, KY, eigenvalues[4, :, :], cmap='spring')
                        ax.set_xlabel('K_x')
                        ax.set_ylabel('K_y')
                        ax.set_zlabel('E')
                        ax.view_init(30,60)
                        mp.title(f'$t_1 = {t1:.2f}, t_2 = {t2}, t_3 = {t3}, 4,5 band $')
                        mp.show()
                        break_outer_loop = True
                        break                    
                if break_outer_loop:
                    break_outer_loop = False
                    ax_index += 1
                    break


# +
condition1, condition2 = '*0', ''
a = findline(KX, KY, condition1, condition2, abs_tol = 1e-7)
condition1 , condition2 = '*(-1)+np.pi', ''
a = np.vstack((a, findline(KX, KY, condition1, condition2, abs_tol = 1e-7,reverse = True)))
condition1, condition2 = '', '*0'
a = np.vstack((a, findline(KX, KY, condition1, condition2, abs_tol = 1e-7, reverse = True, endpoint = True)))

z_points = []

for i in range(len(eigenvalues)):
    z_tmp=[]
    for j,k in a:
        z_tmp.append(eigenvalues[i][j][k])
    z_points.append(z_tmp)
    x_a = range(len(z_points[i]))
    mp.plot(x_a, z_points[i])

mp.xlabel('$\Gamma - X - Y - \Gamma$', fontsize = 16)
mp.xticks([])
mp.ylabel('E', rotation = 0)
mp.title(f'$Biphenylene\;t_1 = {t1}, t_2 = {t2}, t_3 = {t3}$')
mp.show()



# +
def findline(X, Y, sign1, sign2, rel_tol=1e-09, abs_tol=5e-2, reverse = False, endpoint = False):
    from numpy import where, dstack, isclose, flip
    arr = dstack((X, Y))
    result_indices = []
    for i in range(len(X[1])):
        for j in range(len(Y)):
            a = eval(f"{arr[j][i][0]} {sign1}")
            b = eval(f"{arr[j][i][1]} {sign2}")
            if isclose(a, b, rel_tol, abs_tol):
                result_indices.append([i, j])
    if reverse == True:
        result_indices = result_indices[::-1]
    if endpoint == False:
        result_indices = result_indices[:-1]
    return result_indices

def find_nearest_index(arr, target_coord):
    distances = np.linalg.norm(arr - target_coord, axis=2)
    nearest_index = np.unravel_index(np.argmin(distances), distances.shape)
    return nearest_index

def drawline(X, Y, startpoint, endpoint, abs_tol=0.0, reverse = False, endindex = False):
    from numpy import where, dstack, isclose, flip, linspace
    from numpy.linalg import norm
    arr = dstack((X, Y))
    result_indices = []
    incline = (endpoint[1]-startpoint[1])/(endpoint[0]-startpoint[0])
    intercept = incline*startpoint[0]-startpoint[1]
    if endpoint[0] - startpoint[0] == 0:
        incline, intercept = 0
        temp = startpoint[0]
    startpoint = find_nearest_index(arr, startpoint)
    endpoint = find_nearest_index(arr, endpoint)
    print(incline, intercept)
    #line = incline*X+intercept [startpoint ~ endpoint]
    
    for i in range(len(X[1])):
        for j in range(len(Y)):
            a = arr[j][i][0]*incline+intercept
            b = arr[j][i][1]+temp
            print(a)
            if isclose(a, b, atol=abs_tol) and isclose(norm(arr[j][i]- np.array(startpoint))+norm(arr[j][i]-np.array(endpoint)) , norm(np.array(endpoint)-np.array(startpoint)), atol = 3e-2):
                result_indices.append([i, j])
    if reverse == True:
        result_indices = result_indices[::-1]
    if endindex == False:
        result_indices = result_indices[:-1]
    return result_indices



def find_line_segment (X, Y, sign1, sign2, startpoint, endpoint, rel_tol=1e-6, abs_tol=1e-1, reverse = False, endindex = False):
    from numpy import where, dstack, isclose, flip
    arr = dstack((X, Y))
    result_indices, indices = [], []
    for i in range(len(X[1])):
        for j in range(len(Y)):
            a = eval(f"{arr[j][i][0]} {sign1}")
            b = eval(f"{arr[j][i][1]} {sign2}")
            print(type(a))
            if isclose(a, b, atol = abs_tol):
                print(indices[i],[j])
                indices.append([i, j])
    print(indices)
    startpoint = find_nearest_index(arr, startpoint)
    endpoint = find_nearest_index(arr, endpoint)
    start_index = indices.index(startpoint)
    end_index = indices.index(endpoint)
    for a in range(len(indices)):
        if start_index <= a <= end_index or end_index <= a <= start_index:
            result_indices.append(indices[a])
    
    if reverse == True:
        result_indices = result_indices[::-1]
    if endindex == False:
        result_indices = result_indices[:-1]
    return result_indices



# +
a_real = 1.45
b_real = 1.46
c_real = a_real
t1, t2, t3 = 1, 0.5, 1.5
kx = np.linspace(-np.pi,np.pi, 500)
ky = np.linspace(-np.pi,np.pi, 500)
KX, KY = np.meshgrid(kx, ky)

a = [0] #a[0]
a.append(np.array([sqrt(3)/2,-1/2]))  #a[1]
a.append(np.array([0,-1])) #a[2]
a.append(np.array([-sqrt(3)/2,-1/2]))           #a[3]
a = [element * a_real for element in a]
b = [0]
b.append(np.array([a_real,0]))
b.append(np.array([0,-b_real]))
c = [0] #Lattice parameter
c.append(np.array([a_real*(1+np.sqrt(3)), 0]))
c.append(np.array([0, b_real]))


n=6
eigenvalues = np.zeros((n, len(kx), len(ky)), dtype=complex)

for i in range(len(kx)):
    for j in range(len(ky)):
        H = np.zeros([n,n], dtype=complex)
        H_t1 = t1*np.diag(np.concatenate(( [exp(1j*(kx[i]*a[c][0]+ky[j]*a[c][1])) for c in [1,2,3]], 
                                            [exp(-1j*(kx[i]*a[c][0]+ky[j]*a[c][1])) for c in [1,2]]
                                            )), k=1)

        H_t2 = t2*exp(1j*(kx[i]*b[1][0]+ky[j]*b[1][1])/np.linalg.norm(c[1]))
        H_t3 = t3*exp(1j*(kx[i]*b[2][0]+ky[j]*b[2][1])/np.linalg.norm(c[2]))

        H[0][3], H[0][5], H[1][5], H[2][4] = H_t3, exp(1j*(kx[i]*a[3][0]+ky[j]*a[3][1])), H_t2, H_t2
        H += H_t1
        H += H.conj().T

        eigenvalues[:, i, j],_ = np.linalg.eigh(H)
fig = mp.figure(figsize=(8,8))

ax = fig.add_subplot(111,projection = '3d')
for i in range(len(eigenvalues)//2):
    ax.plot_surface(KX, KY, eigenvalues[2*i, :, :], cmap= 'cool')
    ax.plot_surface(KX, KY, eigenvalues[2*i+1,:, :], cmap='spring')

ax.set_xlabel('K_x')
ax.set_ylabel('K_y')
ax.set_zlabel('E')
ax.view_init(30,60)
mp.title(f'$Biphenylene\;t_1 = {t1}, t_2 = {t2}, t_3 = {t3}$')
mp.show()


# +
condition1, condition2 = '*0', ''
startpoint, endpoint = [0,0], [np.pi, 0]
a = find_line_segment(KX, KY, condition1, condition2, startpoint, endpoint, abs_tol = 1e-4)
print('first line')
# startpoint, endpoint = [np.pi,0], [0, np.pi]
# a = np.vstack((a, drawline(KX, KY, startpoint, endpoint, abs_tol = 1e-4)))
# print('second line')
# startpoint, endpoint = [0, np.pi], [0, 0]
# a = np.vstack((a, drawline(KX, KY, startpoint, endpoint, abs_tol = 1e-4, endindex = True)))
# print('last')
z_points = []

for i in range(len(eigenvalues)):
    z_tmp=[]
    for j,k in a:
        z_tmp.append(eigenvalues[i][j][k])
    z_points.append(z_tmp)
    x_a = range(len(z_points[i]))
    mp.plot(x_a, z_points[i])

mp.xlabel('$\Gamma - X - Y - \Gamma$', fontsize = 16)
mp.xticks([])
mp.ylabel('E', rotation = 0)
mp.title(f'$Biphenylene\;t_1 = {t1}, t_2 = {t2}, t_3 = {t3}$')
mp.show()



# +
import numpy as np

def find_nearest_index(arr, target_coord):
    # 배열에서 각 좌표쌍과 입력된 좌표쌍 간의 거리 계산
    distances = np.linalg.norm(arr - target_coord, axis=2)
    # 가장 작은 거리를 가진 좌표쌍의 인덱스 찾기
    nearest_index = np.unravel_index(np.argmin(distances), distances.shape)
    print(np.argmin(distances))
    return nearest_index

# 예제 3차원 배열 (3x3x2)
my_array = np.array([[[1, 2], [3, 4], [5, 6]],
                     [[7, 8], [9, 10], [11, 12]],
                     [[13, 14], [15, 16], [17, 18]]])

# 찾고자 하는 좌표
target_coord = [8, 10]

# 가장 가까운 지점의 인덱스 찾기
nearest_index = find_nearest_index(my_array, target_coord)

print("가장 가까운 지점의 인덱스:", nearest_index)

# -

print(a)

# ##### 
