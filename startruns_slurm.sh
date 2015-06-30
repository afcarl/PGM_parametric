source /com/extra/R/3.1.0/load.sh
rm torun*

NUM_THREADS=1
COMMANDFILE=torun.sh
rm -f ${COMMANDFILE}
while read LIST_FILE
do
echo "qx -m 4G -p --nodes=1 -c 1 -t 24:00:00 --dir=/home/michal/molpros/BRCA_study_parametric/2way \
 -i=/home/michal/molpros/BRCA_study_parametric/2way/2d_param_8fold_TvsAN.R \
 -i=/home/michal/molpros/BRCA_study_parametric/2way/data_BRCA.RData \
-o=*_model.RData ${LIST_FILE}" >> ${COMMANDFILE}
done < commandfile.txt
chmod +x torun.sh
./torun.sh

