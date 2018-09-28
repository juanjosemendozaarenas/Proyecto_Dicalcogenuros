cd ..
cd initfiles
ls initi* > init.txt
mv init.txt ../scripts/.
cd ..
cd scripts
for job in $(cat init.txt); do
        echo qsub -v JOB=$job qsub_DMRG.pbs
	qsub -v JOB=$job qsub_DMRG.pbs
done

rm *.txt


