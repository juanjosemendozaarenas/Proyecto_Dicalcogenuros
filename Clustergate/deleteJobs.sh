#****************************************************************************************************
#************************************* SCRIPT deleteJobs.sh  ****************************************
# Borra una serie de trabajos en un cluster con arquitectura PBS. Los jobs deben tener sus IDs en   *
# orden para que el script funcione. 								    *
# Universidad de los Andes									    *
# Autor :  Juan Sebastian Diaz								   	    *
#                                   							            *
# Date : May 2018		                                                                    *
#                                                                                                   *
#****************************************************************************************************
#****************************************************************************************************


#!/bin/bash

# Lee el numero ID del primer trabajo y lo guarda como "numJob"
echo "Ingrese el ID del primer job"
read numJob

# Lee el numero de jobs a borrar y lo guarda como "n"
echo "Ingrese el numero de jobs que desea borrar"
read n

for ((c=numJob;c<=n+numJob-1;c++))
do
qdel $c.clustermaster.uniandes.edu.co
done



