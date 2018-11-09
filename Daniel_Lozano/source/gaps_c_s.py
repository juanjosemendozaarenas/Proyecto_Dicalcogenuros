import numpy as np
import matplotlib.pyplot as plt

a_4=[1.0,1.5,2.0,2.5,3.]
Ds_4=[0.11039,0.1591,0.5353,1.273,2.125]
Dc_4=[1.3515,1.5061,1.5498,2.0919,2.861]

a_1=[1.0,1.5,2.0,2.5,3.0]
Ds_1=[0.01689,0.0824,0.0969,1.0420,2.0259]
Dc_1=[0.0737,0.1549,0.1720,1.1032,2.081]



plt.plot(a_4,Ds_4,"ks-",label=" $ \Delta_s,\ U=4  $")
plt.plot(a_4,Dc_4,"ks--",alpha=0.3,label="$ \Delta_c, \ U=4 $")
plt.plot(a_1,Ds_1,"b^-",label=" $ \Delta_s,\ U=1 $")
plt.plot(a_1,Dc_1,"b^--",alpha=0.3,label="$ \Delta_c, \ U=1 $")
plt.legend()
plt.grid()
plt.xlim(0,3)
plt.ylim(0,5)
plt.xlabel("$ t_{\\bot} $",size=15)
plt.ylabel("$ \Delta $",size=15)
plt.title("$ \mathrm{Ladder:}32\\times 2$")
plt.savefig("prueba_codigo_ladders.png")
plt.show()
