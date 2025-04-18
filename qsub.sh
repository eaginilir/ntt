#!/bin/sh
#PBS -N qsub
#PBS -e test.e
#PBS -o test.o

/usr/local/bin/pssh -h $PBS_NODEFILE mkdir -p /home/${USER} 1>&2
scp master_ubss1:/home/${USER}/ntt/main /home/${USER} 1>&2
scp -r master_ubss1:/home/${USER}/ntt/files/ /home/${USER}/ 1>&2
/usr/local/bin/pscp -h $PBS_NODEFILE /home/${USER}/main /home/${USER} 1>&2

/home/${USER}/main
rm /home/${USER}/main
scp -r /home/${USER}/files/ master_ubss1:/home/${USER}/ntt/ 2>&1
rm -r /home/${USER}/files/
