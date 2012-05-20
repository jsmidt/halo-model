from pylab import *

[a, b] = loadtxt('matterpower0.dat',unpack=True)
loglog(a,b,label="Z = 0")
[a, b] = loadtxt('matterpower01.dat',unpack=True)
loglog(a,b,label="Z = 0.1")
[a, b] = loadtxt('matterpower1.dat',unpack=True)
loglog(a,b,label="Z = 1")
[a, b] = loadtxt('matterpower10.dat',unpack=True)
loglog(a,b,label="Z = 10")
xlim(1e-3,10)
ylim(1,1e5)
legend().draw_frame(0)
xlabel('k')
ylabel('P(k)')

exit()
