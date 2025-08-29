mkdir -p logs out

for q in {2..72}; do
    echo -e "\n----------- Q: $q/72\n"
    ABI_SCRIPT="s$q.abi"
    # generate abinit script
    cat ddb-gkk.abi | sed s/%QIDX%/$q/ > $ABI_SCRIPT
    # run the abinit script
    time { mpirun -n 4 abinit $ABI_SCRIPT | tee logs/Q$q.log | grep -E "== DATASET|frozen|Perturbation|wfk_read|dfpt_scf"; }
    # remove heavy wavefunctions
    rm out/Q${i}_DS1_1WF*
    rm out/Q${i}_DS1_DEN*
    # remove the script
    rm $ABI_SCRIPT
done

echo "=== Done!"
