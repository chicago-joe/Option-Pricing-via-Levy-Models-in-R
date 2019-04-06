import numpy as np
def main():
    t = np.arange(start=0,stop=0.5,step=0.001)
    x = np.sin(2*np.pi*10*t)
    
    import matplotlib.pyplot as plt
    plt.subplot(2,1,1)
    plt.plot(t,x)
    plt.title('x[n] - original signal')
    plt.xlabel('n')
    plt.ylabel('x[n]')
    
    z = analytic_signal(x)
    
    plt.subplot(2,1,2)
    plt.plot(t,z.real,'k',label='Real(z[n])')
    plt.plot(t,z.imag,'r',label='Imag(z[n])')
    plt.title('Components of Analytic signal')
    plt.xlabel('n')
    plt.ylabel('z_r[n] and z_i[n]')
    plt.legend()

def analytic_signal(x):
    from scipy.fftpack import fft,ifft
    N = len(x)
    X = fft(x,N)
    h = np.zeros(N)
    h[0] = 1
    h[1:N//2] = 2*np.ones(N//2-1)
    h[N//2] = 1
    Z = X*h
    z = ifft(Z,N)
    return z

if __name__ == '__main__':
    main()