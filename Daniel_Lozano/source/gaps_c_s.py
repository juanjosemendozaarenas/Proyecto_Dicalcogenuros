import numpy as np
import matplotlib.pyplot as plt

a=[1.5,2.0,2.5,3.]
Ds=[0.1591,0.5353,1.273,2.125]
Dc=[1.5061,1.5498,2.0919,2.861]

plt.plot(a,Ds,"o-",label=" $ \Delta_s  $")
plt.plot(a,Dc,"^-",label="$ \Delta_c  $")
plt.legend()
plt.grid()
plt.xlim(0,3)
plt.ylim(0,5)
plt.xlabel("$ t_{\\bot} $",size=15)
plt.ylabel("$ \Delta $",size=15)
plt.title("$ U=4,\ \mathrm{Ladder:}32\\times 2$")
plt.savefig("prueba_codigo_ladders.png")
plt.show()
