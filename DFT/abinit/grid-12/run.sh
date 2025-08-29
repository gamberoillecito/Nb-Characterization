mkdir -p logs out

for q in {2..72}; do
    echo -e "\n----------- Q: $q/72\n"
    ABI_SCRIPT="s$q.abi"
    # generate abinit script
    cat ddb-gkk.abi | sed s/%QIDX%/$q/ > $ABI_SCRIPT
    # run the abinit script
    time { mpirun -n 4 abinit $ABI_SCRIPT | tee logs/Q$q.log | grep -E "== DATASET|frozen part|Perturbation|About to read|dfpt_scf|writing gkk"; }
    # remove heavy wavefunctions
    rm out/Q${q}_DS1_1WF*
    rm out/Q${q}_DS1_DEN*
    # remove the script
    rm $ABI_SCRIPT
done

echo "=== Done!"
