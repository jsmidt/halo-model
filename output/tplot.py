from pylab import *


[a, b,c,d,e,f] = loadtxt('test_ukm.dat',unpack=True)
loglog(a,b,label='$10^{11} M_{\odot}$')
loglog(a,c)
loglog(a,d)
loglog(a,e,label='$10^{16} M_{\odot}$')
loglog(a,f)
xlim(0.05,3e3)
ylim(0.001,2.5)
ylabel('u(k|m)',fontsize=16)
xlabel('k  [h/Mpc]',fontsize=16)
legend().draw_frame(0)
show()

[a, b,c,d] = loadtxt('test_powerspec.dat',unpack=True)
loglog(a,b*a**3/2.0/pi**2)
loglog(a,c*a**3/2.0/pi**2)
loglog(a,(b+c)*a**3/2.0/pi**2)
ylim(0.1,1e4)
xlim(0.05,1e3)
show()

[a,b,c,d,e,f] = loadtxt('test_nu_fnu.dat',unpack=True)
loglog(a,sqrt(b))
ylim(0.5,50)
xlim(100,1e15)
show()


exit()

