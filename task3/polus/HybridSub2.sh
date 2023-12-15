for p in 2; do
    for t in 1 2 4 8; do
        mpisubmit.pl -p $p \
            --stdout /dev/null \
            --stderr /dev/null \
            -t $t \
            ./Hybrid80
    done
done