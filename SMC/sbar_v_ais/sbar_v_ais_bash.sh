# sbar_v_ais_bash.sh

for i in {1..500}
do
    echo $i
    date
    # uncomment the line for the param combo(s) you want to run 
    # Rscript sbar_v_ais_master.r ais bimod 6 0 600 1 10 $i
    # Rscript sbar_v_ais_master.r sbar bimod 6 0 600 1 10 $i
    # Rscript sbar_v_ais_master.r ais bimod 11 0 600 1 10 $i
    # Rscript sbar_v_ais_master.r ais bimod 6 0 3000 1 10 $i
    # Rscript sbar_v_ais_master.r ais bimod 6 0 600 1 20 $i

done

echo all done