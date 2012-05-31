from pylab import *

[a, b,c,d,e] = loadtxt('output/test_ukm.dat',unpack=True)
loglog(a,b,label='$10^{11} M_{\odot}$')
loglog(a,c)
loglog(a,d)
loglog(a,e,label='$10^{16} M_{\odot}$')
xlim(0.05,3e3)
ylim(0.001,2.5)
ylabel('u(k|m)',fontsize=16)
xlabel('k  [h/Mpc]',fontsize=16)
title('Compare with Fig. 9 of astro-ph/0206508')
legend().draw_frame(0)
show()


[a, b,c,d] = loadtxt('output/test_powerspec.dat',unpack=True)
loglog(a,b*a**3/2.0/pi**2,label='linear')
loglog(a,c*a**3/2.0/pi**2,label='1h')
loglog(a,d*a**3/2.0/pi**2,label='2h')
loglog(a,(c+d)*a**3/2.0/pi**2,label='total')
ylim(0.1,1e4)
xlim(0.05,1e3)
ylabel('$\Delta^2(k)$',fontsize=16)
xlabel('k  [h/Mpc]',fontsize=16)
title('Compare with Fig. 11 of astro-ph/0206508')
legend(loc='upper left').draw_frame(0)
show()

[a,b,c,d,e] = loadtxt('output/test_nu_fnu.dat',unpack=True)
loglog(a,sqrt(b))
ylim(0.5,30)
xlim(1e2,1e15)
ylabel('$\sigma(m)$',fontsize=16)
xlabel('$M_{\odot}$',fontsize=16)
title('Compare with Fig. 5 of astro-ph/0010468')
show()

[a,b,c,d,e] = loadtxt('output/test_nu_fnu.dat',unpack=True)
loglog(c,d)
ylim(0.001,0.5)
xlim(0.1,30)
ylabel(r'$\nu f(\nu)$',fontsize=16)
xlabel(r'$\nu = (\delta_c/\sigma)^2$',fontsize=16)
title('Compare with Fig. 2 of astro-ph/9901122')
show()

[a,b,c,d,e] = loadtxt('output/test_nu_fnu.dat',unpack=True)
loglog(a,e)
xlim(1e8,1e16)
#ylim(0.1,10)
show()


[a,b,c,d,e] = loadtxt('output/test_nu_fnu.dat',unpack=True)
loglog(a,c)
xlim(1e2,1e16)
#ylim(0.1,10)
show()


exit()

