#Set parameters
function heatsolve(;
t_start = 0,      #start time (s)
t_stop = 100*365*24*60*60,     #stop time (s)
H = 60,            #soil depth (m)
dt = 1800,          #time step (s)
dz = 1,           #height step (m)
D = 7 * 10^(-7))      #soil diffusivity (m^2s^-1)
    #Define initial & boundary condition at surface (z=0)
    N = convert(Int64, (t_stop - t_start) ÷ dt)
    M = convert(Int64, H ÷ dz + 1)
    IC = 273.15 * ones(M)
    BC = 280 * ones(N)

    #Initialize array for temperature
    T = zeros(M, N)
    T[1,:] = BC
    T[:,1] = IC

    #Solve for T
    #using finite difference in space and forward Euler in time
    #and "perfect insulation" neuman boundary condition dT/dz(t,-H)=0 at soil bottom
    ξ = dt*D*(dz^(-2))
    for n in 1:N-1
        for m in 2:M-1
            T[m,n+1] = T[m,n] + ξ*(T[m+1,n]-2*T[m,n]+T[m-1,n])
        end
        T[M,n+1] = T[M,n] + ξ*(-T[M,n]+T[M-1,n])
    end
    return T
end
