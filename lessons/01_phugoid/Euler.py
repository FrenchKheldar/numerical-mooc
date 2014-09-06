import numpy
import matplotlib.pyplot as plt

def EulerMethodPhugoid(u,z,N,dt,g,zt):
    for n in range(1,N):
        uprime = g * (1. - z[n-1] / zt)
        u[n] = u[n-1] + dt * uprime
        z[n] = z[n-1] + dt * u[n-1] 
        if n < 10:
            print n,u[n],z[n]
    return (u,z)

T = 100.0
dt = 0.05
N = int(T/dt)+1
#%timeit t = numpy.linspace(0.0, T, N)
#print t

#This is actually the fastest method
 
t = numpy.arange(0.0,T+dt,dt)

# initial conditions
z0 = 100.  #altitude
v  = 10.   #upward velocity resulting from gust
zt = 100.
g  = 9.81

u = numpy.array([z0, v])

# initialize an array to hold the changing angle values
z = numpy.zeros(N)
uz = numpy.zeros(N)
uzprime = numpy.zeros(N)
z[0] = z0  
uz[0] = v

(uz,z) = EulerMethodPhugoid(uz,z,N,dt,g,zt)

plt.figure(figsize=(10,4))   #set plot size
plt.ylim(40,160)             #y-axis plot limits
plt.tick_params(axis='both', labelsize=14) #increase font size for ticks
plt.xlabel('t', fontsize=14) #x label
plt.ylabel('z', fontsize=14) #y label
plt.plot(t,z, 'k-');
plt.show()

z_exact = v*(zt/g)**.5*numpy.sin((g/zt)**.5*t)+\
                    (z0-zt)*numpy.cos((g/zt)**.5*t)+zt

plt.figure(figsize=(10,4))
plt.ylim(40,160)             #y-axis plot limits
plt.tick_params(axis='both', labelsize=14) #increase font size for ticks
plt.xlabel('t', fontsize=14) #x label
plt.ylabel('z', fontsize=14) #y label
plt.plot(t,z)
plt.plot(t, z_exact)
plt.legend(['Numerical Solution','Analytical Solution']);
plt.show()

# time-increment array
dt_values = numpy.array([0.1, 0.05, 0.01, 0.005, 0.001, 0.0001])

# array that will contain solution of each grid
z_values = numpy.empty_like(dt_values, dtype=numpy.ndarray)

for i, dt in enumerate(dt_values):
    N = int(T/dt)+1    # number of time-steps
    ### discretize the time using numpy.linspace() ###
    t = numpy.linspace(0.0, T, N)

    # initial conditions
    u = numpy.array([z0, v])
    uz = numpy.empty_like(t)
    z = numpy.empty_like(t)
    z[0] = z0
    uz[0] = v
    
    # time loop - Euler method
    (uz,z) = EulerMethodPhugoid(uz,z,N,dt,g,zt)
    
    z_values[i] = z.copy()    # store the total elevation calculation grid i

def get_error(z, dt):
    """Returns the error relative to analytical solution using L-1 norm.
    
    Parameters
    ----------
    z : array of float
        numerical solution.
    dt : float
        time increment.
        
    Returns
    -------
    err : float
        L_{1} norm of the error with respect to the exact solution.
    """
    N = len(z)
    t = numpy.linspace(0.0, T, N)
    
    z_exact = v*(zt/g)**.5*numpy.sin((g/zt)**.5*t)+\
                            (z0-zt)*numpy.cos((g/zt)**.5*t)+zt
    
    return dt * numpy.sum(numpy.abs(z-z_exact))

error_values = numpy.empty_like(dt_values)

for i, dt in enumerate(dt_values):
    ### call the function get_error() ###
    error_values[i] = get_error(z_values[i], dt)

plt.figure(figsize=(10, 6))
plt.tick_params(axis='both', labelsize=14) #increase tick font size
plt.grid(True)                         #turn on grid lines
plt.xlabel('$\Delta t$', fontsize=16)  #x label
plt.ylabel('Error', fontsize=16)       #y label
plt.loglog(dt_values, error_values, 'ko-')  #log-log plot
plt.axis('equal')                      #make axes scale equally;
plt.show()
