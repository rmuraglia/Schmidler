# t_distn_bash.sh

for i in {1..500}
do
    echo $i
    date
    # uncomment the line for the param combo(s) you want to run 
    # Rscript t_distn_master.r ais tdistn 6 0 600 0.5 10 $i 
    # Rscript t_distn_master.r seqbar tdistn 6 0 600 0.5 10 $i 
    # Rscript t_distn_master.r crooks tdistn 6 0 600 0.5 10 $i 
    # Rscript t_distn_master.r pcrooks tdistn 6 0 600 0.5 10 $i 
done

echo all done