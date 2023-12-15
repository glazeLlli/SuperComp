for p in 1 2 4; do
        mpisubmit.pl -p $p \
            --stdout /dev/null \
            --stderr /dev/null \
            ./mpiPure
done