# ais_pathcompare_bash.sh

for i in {1..30}
do
    echo $i
    date
    # uncomment the line for the path heigh you want to run
    # Rscript ais_pathcompare_master.r ais norm 21 0 200 1 25 $i
    # Rscript ais_pathcompare_master.r ais norm 21 1 200 1 25 $i
    # Rscript ais_pathcompare_master.r ais norm 21 2 200 1 25 $i
    # Rscript ais_pathcompare_master.r ais norm 21 3 200 1 25 $i
    # Rscript ais_pathcompare_master.r ais norm 21 4 200 1 25 $i
done

echo all done