# ais_bandit_bash.sh

for i in {1..30}
do
    echo $i
    for j in greedy proptnl equal
    do
        date
        echo $j
        Rscript ais_bandit_master.r $j $i
    done
done

echo all done
