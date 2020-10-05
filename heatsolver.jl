
function heatsolve(::Type{MyType};
t_start = 0,      #start time (s)
t_stop = 100*365*24*60*60,     #stop time (s)
H = 60,            #soil depth (m)
dt = 1800,          #time step (s)
dz = 1,           #height step (m)
D = 7 * 10^(-7)      #diffusivity (m^2s^-1), by default that taken by Dawson et al for "soil" (off by factor ~2 from some online sources)
) where MyType<:AbstractFloat
    #Define initial & boundary condition at surface (z=0)
    N = Int64.((t_stop - t_start) ÷ dt)
    M = Int64.(H ÷ dz + 1)
    IC = MyType.(273.15 * ones(M))
    BC = MyType.(280 * ones(N))

    #Initialize array for temperature
    T = MyType.(zeros(M, N))
    T[1,:] = BC
    T[:,1] = IC

    #Solve for T
    #using finite difference in space and forward Euler in time
    #and "perfect insulation" neuman boundary condition dT/dz(t,-H)=0 at soil bottom
    dt = MyType.(dt)
    dz = MyType.(dz)
    D = MyType.(D)
    ξ = dt*D*(dz^(-2))
    for n in 1:N-1
        for m in 2:M-1
            T[m,n+1] = T[m,n] + ξ*(T[m+1,n]-2*T[m,n]+T[m-1,n])
        end
        T[M,n+1] = T[M,n] + ξ*(-T[M,n]+T[M-1,n])
    end
    T = Float64.(T)
    return T
end

function diurnal(t)
    k = 1/(24*60*60) #set frequency 1/(1 day)
    6.85*0.5*(1+sin(2*pi*k*t)) #define a solar forcing with diurnal cycle
end

function diurnal_heatsolve(::Type{MyType};
t_start = 0,      #start time (s)
t_stop = 100*365*24*60*60,     #stop time (s)
H = 60,            #soil depth (m)
dt = 1800,          #time step (s)
dz = 1,           #height step (m)
D = 7 * 10^(-7)      #soil diffusivity (m^2s^-1)
) where MyType<:AbstractFloat
    #Define initial & boundary condition at surface (z=0)
    times = t_start:dt:t_stop
    N = length(times)
    M = Int64.(H ÷ dz + 1)
    IC = MyType.(273.15 * ones(M))
    BC = MyType.(273.15*ones(N) + diurnal.(times))

    #Initialize array for temperature
    T = MyType.(zeros(M, N))
    T[1,:] = BC
    T[:,1] = IC

    #Solve for T
    #using finite difference in space and forward Euler in time
    #and "perfect insulation" neuman boundary condition dT/dz(t,-H)=0 at soil bottom
    dt = MyType.(dt)
    dz = MyType.(dz)
    D = MyType.(D)
    ξ = dt*D*(dz^(-2))
    for n in 1:N-1
        for m in 2:M-1
            T[m,n+1] = T[m,n] + ξ*(T[m+1,n]-2*T[m,n]+T[m-1,n])
        end
        T[M,n+1] = T[M,n] + ξ*(-T[M,n]+T[M-1,n])
    end
    T = Float64.(T)
    return T
end

using Plots

function heatplot(::Type{MyType};
    t_start = 0,      #start time (s)
    t_stop = 100*365*24*60*60,     #stop time (s)
    H = 60,            #soil depth (m)
    dt = 1800,          #time step (s)
    dz = 1,           #height step (m)
    D = 7 * 10^(-7)      #diffusivity, default for soil (m^2s^-1)
    ) where MyType <: AbstractFloat
    T = heatsolve(MyType, t_start=t_start, t_stop=t_stop, H=H, dt = dt, dz=dz, D=D)
    Tweekly = T[:,begin:336:end]
    heatmap(Tweekly)
end

function diurnal_heatplot(::Type{MyType}; dt = 1800) where MyType <: AbstractFloat
    T = diurnal_heatsolve(MyType, dt = dt)
    Tweekly = T[:,begin:336:end]
    heatmap(Tweekly)
end
