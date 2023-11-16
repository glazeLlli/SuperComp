bsub -n 1 -q normal -o /dev/null -e /dev/null -R  "affinity[core(10,same=socket,exclusive=(socket,alljobs)):membind=localonly:distribute=pack(socket=1)]"  "OMP_NUM_THREADS=1" ./seq
