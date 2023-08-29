#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type END
#SBATCH --mail-user b-barckmann@chu-montpellier.fr
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 6
#SBATCH --mem 64GB



module load bedtools/2.30.0

input_list=("WGcontrol" "9" "11" "13" "15" "17" "19")
prefix=bin100_norm_RFD_cutoff50_bs0.1kb
distance=501 # max distance between minus and plus feature 
inputdir=normalised_RFD/

for x in "${input_list[@]}"; do
# adding the strand for later seperation in plus and minus files 
    awk 'BEGIN { OFS="\t" } {if ($4 > 0 ) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "." "\t" "-"}   else if ($4 < 0  ) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "." "\t" "+"} }' $inputdir/$x\_$prefix\.bedgraph >  $inputdir/$x\_$prefix\_addStrand.bedgraph &&
# merging neigbouring features 
    bedtools merge -i $inputdir/$x\_$prefix\_addStrand.bedgraph -s -c 4,5,6 -o mean,distinct,distinct > $inputdir/$x\_$prefix\_addStrand_merged.bedgraph  &&
# seperation in plus and minus feature files
    awk 'BEGIN { OFS="\t" } {if ($6 == "+" ) {print $0}}' $inputdir/$x\_$prefix\_addStrand_merged.bedgraph > $inputdir/$x\_$prefix\_PLUS.bedgraph &&
    awk 'BEGIN { OFS="\t" } {if ($6 == "-" ) {print $0}}' $inputdir/$x\_$prefix\_addStrand_merged.bedgraph > $inputdir/$x\_$prefix\_MINUS.bedgraph  &&
# selecting minus upstream plus downstream features pairs 
    bedtools closest -iu -D ref -t all -a $inputdir/$x\_$prefix\_MINUS.bedgraph -b $inputdir/$x\_$prefix\_PLUS.bedgraph |  awk -v dist="$distance"  'BEGIN { OFS="\t" } {if ($13 <= dist ) {print $0} }' > $inputdir/ClosestIU_maxdist$distance\_$x\_$prefix\_MINUS.bedgraph  &&
    bedtools closest -id -D ref -t all -a $inputdir/$x\_$prefix\_PLUS.bedgraph  -b $inputdir/$x\_$prefix\_MINUS.bedgraph |  awk -v dist="$distance"  'BEGIN { OFS="\t" } {if ($13 >= -dist ) {print $0} }' > $inputdir/ClosestID_maxdist$distance\_$x\_$prefix\_PLUS.bedgraph  &&
# generate ORI plus and Minus files 
    cat $inputdir/ClosestIU_maxdist$distance\_$x\_$prefix\_MINUS.bedgraph | awk 'BEGIN { OFS="\t" } {if ($13 > 0 ) {print $1 "\t" $3 "\t" $8 "\t" ($4+($10*-1))/2 "\t" $4 "\t" $10 "\t" $6 "\t" $13}}'  > $inputdir/Minus_ORIs_$x\_$prefix\.bedgraph &&
    cat $inputdir/ClosestID_maxdist$distance\_$x\_$prefix\_PLUS.bedgraph | awk 'BEGIN { OFS="\t" } {if ($13 < 0 ) {print $1 "\t" $9 "\t" $2 "\t" (($4*-1)+$10)/2 "\t" $4 "\t" $10 "\t" $6 "\t" $13}}' > $inputdir/Plus_ORIs_$x\_$prefix\.bedgraph  &&  
# combining the plus and minus ORIs and removing all doublons and all with one identical start or end point while keeping the smaller feature
    cat $inputdir/Plus_ORIs_$x\_$prefix\.bedgraph  $inputdir/Minus_ORIs_$x\_$prefix\.bedgraph | sort -k1,1 -k2,2n |  awk '!a[$1 $2 $3]++' |   awk '!seen[$1 $2]++' |sort  -k1,1 -rk2 | awk '!seen[$1 $3]++' | sort -k1,1 -k2,2n > $inputdir/ORI_$x\_$prefix\.bedgraph && 
# to filter the ORIs according to the score of their plus and Minus features
cat $inputdir/ORI_$x\_$prefix\.bedgraph | awk   'BEGIN { OFS="\t" } {if ($4 >= 0.3 ) {print $0} }' >  $inputdir/ORI_$x\_$prefix\_th0.3.bedgraph &&
cat $inputdir/ORI_$x\_$prefix\.bedgraph | awk   'BEGIN { OFS="\t" } {if ($7 == "+" && $5 <= -0.2 && $6 >= 0.2 ) {print $0} else if ($7 == "-" && $5 >= 0.2 && $6 <= -0.2) {print $0} }' >  $inputdir/ORI_$x\_$prefix\_th0.2P+M.bedgraph &&
cat $inputdir/ORI_$x\_$prefix\.bedgraph | awk   'BEGIN { OFS="\t" } {if ($7 == "+" && $5 <= -0.3 && $6 >= 0.3 ) {print $0} else if ($7 == "-" && $5 >= 0.3 && $6 <= -0.3) {print $0} }' >  $inputdir/ORI_$x\_$prefix\_th0.3P+M.bedgraph \
;done