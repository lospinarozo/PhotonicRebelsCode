#######################################################################################
#####################          PHYLOGENETIC ANALYSIS TOOLS         ####################
## Questions about these scripts: ahugall@unimelb.edu.au, ahugall@museum.vic.gov.au
## Unix scripts (mostly awk)
## Fasta format, converted to single line
## Scripts used in conjunction with BioEdit (https://bioedit.software.informer.com/7.2/), Notepad++ and Excel spreadsheets (not included)





#########################################################################
################    GENBANK (NCBI) FASTA DOWNLOAD PARSING
## Simple taxon/locus code downloads and also BLAST searches using reference sequence
## Reference sequence can be nucleotide or amino acid (AA)
## /mnt/c/Users/ahugall/Documents/FASTA/seqdump.txt
## NCBI_parse_list.txt - revise to suit both standard and BLAST downloads
## Option to delete accessions/taxa null (ie _sp. etc) 
## BLASTn and BLASTtn : if aligned section run extra code below; for mtDNA, exclude OV/OU accessions (genome assemblies)
# >AF344616.2:1-206,299-557,635-801,891-977 Anthophora montana long-wavelength rhodopsin gene, partial cds


FAM=Halictidae
LOC=COI

### Single line for sequence (get rid of _ in accession code; remove space after sp. too keep more info)
tr -d '\r' < Andrena_${LOC}.fasta | awk '{if(NR==1) {printf $0; printf "\n"; next}
 else if(/^>/) {printf "\n"; printf $0; printf "\n"} else printf $1} END{printf "\n"}' | sed 's/_//1;s/sp. /sp./1' > ${FAM}_${LOC}_p1.fasta


### Single line for sequence (get rid of _ in accession code; remove space after sp. too keep more info)
tr -d '\r' < ${FAM}_${LOC}_seqdump.fasta | awk '{if(NR==1) {printf $0; printf "\n"; next}
 else if(/^>/) {printf "\n"; printf $0; printf "\n"} else printf $1} END{printf "\n"}' | sed 's/_//1;s/sp. /sp./1' > ${FAM}_${LOC}_p1.fasta


### Delete accession with NCBI_parse_list - TO BE REVISED!!
### ADD filter for spurious chars [;:()@ etc]
## Standard download, second field check
awk 'NR==FNR{x[$1]=$2; next}
 {if(/^>/) {q=$0; if(x[$2]!=-1) a=1} else if(a==1) {printf "%s\n%s\n",q,$1; a=0}}' NCBI_parse_list.txt ${FAM}_${LOC}_p1.fasta > ${FAM}_${LOC}_p2.fasta

### All fields version
awk 'NR==FNR{x[$1]=$2; next}
 {if(/^>/) {a=1; q=$0; for(n=1;n<=NF;n++) if(x[$n]==-1) a=0} else if(a==1) {printf "%s\n%s\n",q,$1}}' NCBI_parse_list.txt ${FAM}_${LOC}_p1.fasta > ${FAM}_${LOC}_p2.fasta


### Cull long sequences (genomic bits etc)
awk '{if(/^>/) q=$0; else if(length($1)<3000) printf "%s\n%s\n",q,$1}' ${FAM}_${LOC}_p2.fasta > ${FAM}_${LOC}_p3.fasta


### Alternate fourth info field (use variable, like family)
awk 'NR==FNR{x[$1]=$2; next}
 {if(/^>/) {i=n=1; while(i<=3) {if(x[$n]<1) {f[i]=$n; i++} n++} printf "%s_%s_%s_%s\n",f[1],f[2],f[3],"'$FAM'"}
 else print $0}' NCBI_parse_list.txt ${FAM}_${LOC}_p3.fasta > ${FAM}_${LOC}_p4b.fasta


grep -c '>' ${FAM}_${LOC}_p1.fasta
grep -c '>' ${FAM}_${LOC}_p2.fasta
grep -c '>' ${FAM}_${LOC}_p3.fasta
grep -c '>' ${FAM}_${LOC}_p4b.fasta



#### ORDER MATCHS BY POSITION
## Processes this type of BLASTtn output format
>NM001328485.1:364-1065_Apis_cerana_Apidae
>NM001328485.1:c1065-364
## Account for compliment: if number starts with 'c' remove and subtract from large dummy number (10000)
## Joins marked by ~

echo $FAM $LOC

paste -sd'\t\n' ${FAM}_${LOC}_p4b.fasta > data.txt

## Five field output: accession, start end, taxon, seq
sed 's/:/\t/1;s/-/\t/1;s/_/\t/1' data.txt > data2.txt

awk '{q=$2;p=$3; if(substr($2,1,1)=="c") {q=10000-gensub("c","","g",$2); p=10000-$3} print $1,q,p,$4,$5}' data2.txt > data3.txt

sort -k1,1 -k2,2n data3.txt > data4.txt

wc -l data.txt
wc -l data2.txt
wc -l data3.txt
wc -l data4.txt

awk '{t[$1]=$4; n[$1]++; b[$1]=b[$1]+$3-$2+1; if(n[$1]==1) a[$1]=$5; else a[$1]=a[$1]"~"$5}
 END{for(x in a) printf "%s_%s~%s~%s\n%s\n",x,t[x],n[x],b[x],a[x]}' data4.txt > ${FAM}_${LOC}_p5.fasta

wc -l ${FAM}_${LOC}_p5.fasta
grep -c '>' ${FAM}_${LOC}_p5.fasta





#########################################################################
################    PICK THE LONGEST
## Max bases per taxon (select fields)
## Assumes single line sequence, case insensitive
## Option to sort by size, taxon, original order (keep NR record)
## Check taxon field combination (accession_genus_species_info)


RUN=Halictidae_COI_p5b

## Set taxon field (versions: standard/p5 type output)
# paste -sd'\t\n' ${RUN}.fasta | tr -d '>' | tr -d '\r' | sed 's/_/ /1' > data.txt
paste -sd'\t\n' ${RUN}.fasta | tr -d '>' | tr -d '\r' | sed 's/_/ /1;s/_/ /2' > data.txt

awk '{s=length(gensub(/[-~?nxXN]/,"","g",$NF)); if(s>b[$2]) {b[$2]=s;a[$2]=$NF;c[$2]=$1}} END{for(x in a) print c[x],x,a[x],b[x]}' data.txt | sort -k2,2 > zzz.txt

awk '{a[$2]++} END{for(x in a) print x,a[x]}' data.txt > temp.txt

wc -l data.txt
wc -l zzz.txt
wc -l temp.txt

awk '{printf ">%s_%s_%s\n%s\n",$1,$2,$4,$3}' zzz.txt > ${RUN}_u.fasta

## Option to limit by length (field 4)
awk '{if($4>100) printf ">%s_%s_%s\n%s\n",$1,$2,$4,$3}' zzz.txt > ${RUN}_u.fasta




################    REMOVE DUPLICATE ACCESSIONS
## SINGLE SEQ LINE FASTA FORMAT

RUN=Scarabaeidae_28Sv3

cat ${RUN}.fasta > temp.fasta
 
## first field only
paste -sd'\t\n' temp.fasta | tr -d '\r' | tr -d '>' | sed 's/_/ /1' | awk '{if(++t[$1]==1) printf ">%s_%s\n%s\n",$1,$2,$3}' > ${RUN}_unique.fasta

grep -c '>' temp.fasta
grep -c '>' ${RUN}_unique.fasta



################    FILTER FASTA FILE BY SEQUENCE LENGTH
## assumes seq on single line
## set to 100 characters (nucleotide)

RUN=apoidea_16S_v3

awk '{if(/^>/) q=$0; else if(length(gensub(/[-~?nxXN]/,"","g",$1))>99) printf "%s\n%s\n",q,$1}' ${RUN}.fasta > ${RUN}_sf.fasta

grep -c '>' ${RUN}.fasta
grep -c '>' ${RUN}_sf.fasta



#########################################################################
################    REFORMAT FASTA ALIGNMENT FOR EXCEL DB
## Rearrage to be label, |seq, label into fields separated by _
## Assume sequence all on one line (out of BioEdit)

RUN=Halictidae_COI_check3

## Temp patch to replace space with _
# sed 's/ /_/g' ${RUN}.fasta | paste -sd'\t\n' | tr -d '>' | tr -d '\r' | awk '{printf "%s\t%s\t|%s\n",NR,$1,$2}' > data.txt
sed 's/ /_/g' ${RUN}.fasta | paste -sd'\t\n' | tr -d '>' | tr -d '\r' | awk '{s=length($2); b=length(gensub(/[-~nNxX?]/,"","g",$2)); printf "%s\t%s\t|%s\t%s\t%s\n",NR,$1,$2,s,b}' > data.txt
awk '{print $2}' data.txt | sed 's/_/\t/g' > lf.txt
paste data.txt lf.txt > ${RUN}_exdb.txt


grep -c '>' ${RUN}.fasta
wc -l ${RUN}_exdb.txt



################    REPLACE TAXON LABELS
## Assume single line fasta
## check tlist file for fields
## substitute ? null base chars

NAME=bee_uce_set2_100SM_TV33_JV50-0

cat Bossert_Sless_nn_tlist.txt > tlist.txt

# cat ${NAME}.txt > gene.txt
paste -sd'\t\n' ${NAME}.fasta | tr -d '>' | tr -d '\r' > gene.txt

awk 'NR==FNR{t[$1]=$2; next} {d=t[$1]; if(d==0) d="xx"; printf ">%s\n%s\n",d,$2}' tlist.txt gene.txt | sed 's/?/n/g' > ${NAME}_nn.fasta




################    FASTA <-> PHYLIP  ########### - TO BE REVISED!!
## Input single label, single line sequence data.txt file


### Single line for sequence (get rid of _ in accession code)
tr -d '\r' < /mnt/c/Users/ahugall/Documents/BEES/FASTA/anthophila/${FAM}/${LOC}/seqdump.txt | awk '{if(NR==1) {printf $0; printf "\n"; next}
 else if(/^>/) {printf "\n"; printf $0; printf "\n"} else printf $1} END{printf "\n"}' | sed 's/_//1' > ${FAM}_${LOC}_p1.fasta

tr -d '\r' <  ${RUN}_mft.fasta | awk '{if(/^>/) printf ("\n%s ",$0); else printf $0}' | sed 1d | sed 's/>//1' > data.txt


## two field single line label, seqs
awk 'END{print NR, length($2)}' data.txt > dim.txt
cat dim.txt data.txt > ${RUN}.phy

awk '{printf ">%s\n%s\n",$1,$2}' data.txt >${NAME}_TAAX.fasta




###################################################################################################
##############################      SINGLE LOCUS SITE TRIMMING       ##############################
#### TRIM MARKLINE generalized n-block version
## nucleotide data only
## Set missing taxa percentage limit (MS), block size (CB)
## load elist, process sample: process loci - site completeness, record markline (each site)
## use  www.txt label and yyy.txt single field output

NAME=bee_uce_set2_100SM

CB=13
MS=33

## Sites as fields version
# cat ${NAME}.txt > gene.txt
paste -sd'\t\n' ${NAME}.fasta | tr -d '>' | tr -d '\r' > gene.txt


### Reformat fasta  to single line for sequence
# tr -d '\r' < ${NAME}.fasta | awk '{if(/^>/) printf ("\n%s ",$0); else printf $0}' | sed 1d | sed 's/>//1;s/://g' > gene.txt


awk '{print $1}' gene.txt > www.txt
awk '{print $2}' gene.txt > xxx.txt

wc -l www.txt
wc -l xxx.txt
awk 'END{print NR,NF,length($1)}' xxx.txt


### Single partition sites as fields version (CB = block size)
## REVISE TO HANDLE MIXED CASE [toupper($s)]
awk 'BEGIN{FS=""; b='$CB'} {l=NF; for(s=1; s<=l; s++) if($s~/[acgtACGT]/) c[s]++}
 END{q=NR*'$MS'/100; for(s=1; s<=l; s=s+b)
 {p=m=0; for(n=1;n<=b;n++) if(c[s+n-1]>q) p++; if(p==b) m=1; for(n=1;n<=b;n++) print m}}' xxx.txt > marklineN.txt

## checks
awk '{m[$1]++} END{for(x in m) print x,m[x]}' marklineN.txt

## trim data using markline
awk 'BEGIN{FS=""} NR==FNR{c[NR]=$1;next} {for(s=1;s<=NF;s++)
 {if(c[s]==1) printf "%s",$s} printf "\n"}' marklineN.txt xxx.txt > aaa.txt


### Generate fasta output with dummy consensus line [input: aaa.txt for trimmed, xxx.txt for original]
## dummy (sites as fields version)
awk 'BEGIN{FS=""} {for(i=1; i<=NF; i++) {ss[i,$i]++; cs[$i]++}}
 END{print ">dummy"; for(i=1; i<=NF; i++) {c="-";n=0; for(x in cs)
 {if(x!~/[-~nNxX?]/ && ss[i,x]>n) {n=ss[i,x]; c=x}} printf c} printf "\n"}' aaa.txt > dummy.fasta


paste -d" " www.txt aaa.txt > zzz.txt
cat zzz.txt > ${NAME}_TV${MS}.txt
awk '{printf ">%s\n%s\n",$1,$2}' zzz.txt > gene.fasta
cat dummy.fasta gene.fasta > ${NAME}_TV${MS}.fasta

### Simple nucleotide p-sites summary
awk '{s=length($2); b=length(gensub(/[-~nNxX?]/,"","g",$2)); print $1,s,b,b/s}' zzz.txt > ${NAME}_TV${MS}_summ.txt



#### GENERATE NEW ELIST USING marklineN.txt

cat ${NAME}_elist.txt > elist.txt

awk 'NR==FNR{mc[NR]=$1;next}
 {b=0; for(s=$3;s<=$4;s++) {if(mc[s]==1) b++}
 if(b>0) {print($1,b,t+1,t+b); t=t+b}}' marklineN.txt elist.txt > ${NAME}_${MS}TCT_elist.txt

### Taxon by Partition p-sites summary
## Use trimmed elist
awk 'NR==FNR{ne++; el[ne]=$2; es[ne]=$3; next} {tb=0;printf $1; for(n=1;n<=ne;n++)
 {q=substr($2,es[n],el[n]); b=length(gensub(/[-~nNxX?]/,"","g",q)); tb=tb+b;\
 printf " %0.2f",b/el[n]} printf " %s %s\n",tb,length($2)}' ${NAME}_${MS}TCT_elist.txt zzz.txt > ${NAME}_summ.txt

## checks
awk 'END{print NR,NF}' ${NAME}_summ.txt
awk '{s[length($2)]++} END{for(x in s) print x,s[x]}' zzz.txt



#### SPECIAL
#### FASTA output of elist dummy.fasta output above
awk 'NR==FNR{s=$1; next}
 {e=substr(s,$3,$2); printf ">%s_dummy\n%s\n",$1,e}' dummy.fasta ${NAME}_${MS}TCT_elist.txt > ${NAME}_${MS}TCT_elist_dummy.fasta






##########################################################
###################    RUN IQTREE      ################### 
## C:\Users\ahugall\Documents\iqtree-1.6.12-Windows

## SIMPLE
bin\iqtree -s andrenidae_COI_align_conszzzz.fasta -m TN93+G -pre andrenidae_COI_align_conszzzz_tn93g -nt 2

## SIMPLE BS
bin\iqtree -s andrenidae_COI_align_conszzzz.fasta -m TN93+G -bb 1000 -pre andrenidae_COI_align_tn93g -nt 4


## PARTITION BS
bin\iqtree -s AND_mat1_10TCT.phy -spp AND_mat1_10TCT_p3pmod.best_scheme.nex -bb 1000 -bnni -pre AND_mat1_10TCT_p3pmBS -nt 4

 
## PARTITION MODEL TESTS
bin\iqtree -s MEGAB_opsin1_blast_mft.fasta -m TESTONLY -pre MEGAB_opsin1_blast_mft_mod -nt 4


bin\iqtree -s AND_mat1_10TCT.phy -spp AND_mat1_p3.txt -m TESTMERGEONLY -pre AND_mat1_10TCT_p3pmod -nt 4


## NON-PARAMETRIC BS
bin\iqtree -s COLUMBI_mat3_10TC.phy -spp COLUMBI_mat3_10TC_mod.best_scheme.nex -b 200 -pre COLUMBI_mat3_10TC_p5mb_npBS -nt 4 -o CUCULIFORMES_n4m7_out





#########################################################################
################    MAFFT ALIGNMENT
## options --keeplength, --mapout, --addlong, addfragments, thread 4

RUN=Halictidae_COI_p5b

mafft --nuc --localpair --thread 4 --maxiterate 100 ${RUN}.fasta > ${RUN}_mft.fasta

# mafft --auto --thread -1 --maxiterate 100 ${RUN}.fasta > ${RUN}_mft.fasta

# mafft --add ANCOS_gpicref.fasta ANCOS_ec16h.fasta > ANCOS_add_mft.fasta

### Single line for sequence (get rid of _ in accession code)
tr -d '\r' < mft.fasta | awk '{if(NR==1) {printf $0; printf "\n"; next}
 else if(/^>/) {printf "\n"; printf $0; printf "\n"} else printf $1} END{printf "\n"}' | sed 's/_//1' > ${RUN}_mft_sl.fasta

 
## Options to reformat
tr -d '\r' <  ${RUN}_mft.fasta | awk '{if(/^>/) printf ("\n%s ",$0); else printf $0}' | sed 1d | sed 's/>//1' > data.txt
# paste -sd'\t\n' ${RUN}_mft.fasta | tr -d '>' | tr -d '\r' > data.txt


## Trim ends back to reference (= first line)
# drop null sequences (limit =50)
# trim ref align of large near end gaps (lg>50 outside 0.05<ref>0.95) - keep number of end bases in target (qb, qe)
awk '{l[NR]=$1;a[NR]=$2}
 END{{b=c=e=lg=qb=qe=0; s=length(a[1]); r=length(gensub("-","","g",a[1])); for(n=1;n<=s;n++)
 {if(substr(a[1],n,1)~"-") {if(c==0) b++; else lg++} else
 {if(lg>50) {if(c<r*0.05) {b=n-1;qb=c} else if(c>r*0.95) {e=n-lg-1; qe=r-c; break}} c++;lg=0} if(c==r) {e=n; break}}}
 for(n=1;n<=NR;n++) {q=substr(a[n],b+1-qb,e-b+qe+qb); if(length(gensub("-","","g",q))>50) print l[n],q}}' data.txt > rtgene.txt









############################################################################
#########   SINGLE TREE ANALYSIS SYSTEM  V1 - UNDER DEVELOPMENT    #########
## Trees treated as rooted

#####  MULTI TIP TAXON MONOPHYLY SCAN
## Need to remove node labels (e.g BS values)
## Select appropriate fields to define taxon unit and output taxon label: accession_genus_species_info
## Creates list of mulit-tip taxon units (glist)

LBL=Halictidae_COI_p5b_mft_cons1503_tn93.treefile

cat ${LBL} > phup.nwk


### Parse newick tree version 2 (option to remove ')
## format onto clade/tip lines; tip field delimiter parse (_); remove node label info
# sed 's/:/\t/g;s/[,;]/\n/g;s/(/(\n/g;s/)/\n)/g' phup.nwk | sed 's/_/ /2' | sed -r 's/\).+\t/)\t/1' | awk 'NF' > temp.txt
sed 's/:/\t/g;s/[,;]/\n/g;s/(/(\n/g;s/)/\n)/g' phup.nwk | sed 's/_/ /g' | sed -r 's/\).+\t/)\t/1' | tr -d "'" | awk 'NF' > tempa.txt

## Select taxon fields (print $2"_"$3 = binom, print $2/$1 = genus gene/supermatrix alignments)
awk '{if(NF>2) print $2,$NF; else print $0}' tempa.txt > tempf.txt
cat tempf.txt > temp.txt


## taxon list, with exclusion list
awk 'BEGIN{q[0]=q["og"]=1} {if(/^[()]/) next; else if(q[$1]!=1) t[$1]++} END{for(x in t) if(t[x]>1) print x,t[x]}' temp.txt > glist.txt



TS=$(awk 'NF{c++} END{print c}' glist.txt)
echo $TS 'taxa'

rm -f ${LBL}_MEC_summT.txt
rm -f ${LBL}_MEC_summTX.txt



for i in $(seq 1 $TS)
do
ET=$(sed -n ${i}p glist.txt | cut -d ' ' -f 1)

### Clade stats with group set stats
## Eight fields: node number, start line, end line, tips, set tips, stem length, node height, stem+height
## Note height values only relevant to ultrametric trees

## ET second field version (bl = last field)
awk '{tr[NR,1]=$1;tr[NR,2]=$NF; if($1=="(") {c[++lb]=NR; for(x in c) cc[x]++}
 else if($1==")") {for(x in c) {cc[x]--; if(cc[x]==0 && rb[x]==0) {rb[x]=NR;csl[x]=$NF}}} else
 {for(x in c) if(rb[x]==0) {ct[x]++; if($1=="'$ET'") ctm[x]++}}}
 END{for(x in c) {nh[x]=0;q=rb[x]; while(tr[q,1]==")")
 {nh[x]=nh[x]+tr[q-1,2];q--} print x,c[x],rb[x],ct[x],ctm[x]+0,csl[x]+0,nh[x],csl[x]+nh[x]}}' temp.txt > temp2.txt


### Most exclusive clade info
### MEC stats: tree, [temp2 info], tree tips, av-height, max-height
sort -k5,5n -k4,4nr temp2.txt | awk '{h=h+$7; if($4>t) {t=$4;m=$7}} END{print '$i',$0,t,h/NR,m}' > temp2S.txt

## MEC summary: tree, total nodes, tips, set tips, node height, stem+height, av-height, max-height, (taxa list up 7)
## Using first 3 fields for taxon label - UPDATE TO NEW INPUT TAXON FIELD FORMAT
 awk 'NR==FNR{t[NR,1]=$1; t[NR,2]=$2"_"$3; next}
 {printf "%s %s %s %s %s %s %s %s %s",$1,"'$ET'",$10,$5,$6,$8,$9,$11,$12; for(n=$3;n<$4;n++)
 {if(t[n,1]!=")" && t[n,1]!="(") if(++q>7) {printf "\n";exit} else printf " %s_%s",t[n,1],t[n,2]} printf "\n"}' tempa.txt temp2S.txt >> ${LBL}_MEC_summT.txt

### MEC discounting null taxon labels (0 label) - WORKING!!
awk 'NR==FNR{t[NR,1]=$1; next}
 {printf "%s %s %s %s %s %s %s %s %s",$1,"'$ET'",$10,$5,$6,$8,$9,$11,$12; for(n=$3;n<$4;n++)
 {if(t[n,1]!=")" && t[n,1]!="(") if(t[n,1]==0) q++} printf " %s\n",q+0}' temp.txt temp2S.txt >> ${LBL}_MEC_summTX.txt

done


### Simple summary
awk 'BEGIN{print "'$LBL'"} {if($4!=$5) {nm++; print $1,$2,$3,$4,$5}}
 END{print NR " multi-tip-otu", nm+0 " not-monophyletic"}' ${LBL}_MEC_summT.txt | tee -a ${LBL}_MEC_summTS.txt

awk 'BEGIN{print "'$LBL'"} {if(($4-$NF)!=$5) {nm++; print $1,$2,$3,$4,$5,$NF}}
 END{print NR " multi-tip-otu", nm+0 " not-monophyletic"}' ${LBL}_MEC_summTX.txt | tee -a ${LBL}_MEC_summTXS.txt






####################################################################################
##############     SPECIES-LEVEL CONSENSUS - NCBI accession data      ##############
## Aligned block of data with trinomial taxon label (genus_species_accession)
## Choose appropriate taxon field combination: accession_genus_species_info
## Sequence needs to be all same case

ASB=Halictidae_COI_p5b_mft

## Reformat accession_genus_species_info fasta  (versions: standard/p5 type output)
# paste -sd'\t\n' ${ASB}.fasta | tr -d '>' | tr -d '\r' | sed 's/[?~]/-/g;s/_/ /1' > asb.txt
paste -sd'\t\n' ${ASB}.fasta | tr -d '>' | tr -d '\r' | sed 's/[?~]/-/g;s/_/ /1;s/_/ /2' > asb.txt
## Convert sequence to same case (option for first position | excel format)
# paste -sd'\t\n' ${ASB}.fasta | tr -d '>' | tr -d '\r' | awk '{if(NF>0) print $1,tolower($2)}' | sed 's/[?~]/-/g;s/_/ /1' > asb.txt
awk '{print NF, length($NF)}' asb.txt



## Get multi-tip taxon unit list - excludes singletons
awk '{t[$2]++} END{for(x in t) print x,t[x]}' asb.txt > glist.txt

awk '{if($2>1) print $0}' glist.txt > glistM.txt

wc -l asb.txt

TS=$(awk 'NF{c++} END{print c}' glist.txt)
echo $TS 'taxa'
rm -f ${ASB}_constats.txt
rm -f ${ASB}_cons.txt



for i in $( seq 1 $TS )
do 
sed -n ${i}p glist.txt > RRF1.txt
RRF1=$(awk '{print $2}' RRF1.txt)
NAME=$(awk '{print $1}' RRF1.txt)
echo $RRF1 $NAME

grep -w $NAME asb.txt > tdata.txt

if test $RRF1 -gt 1; then

## Simple consensus, overlap stats version: sample, n seqs, conbase, overlap, missmatch, x-sites, conseq
## Assumes sequences all the same case; default -, max freq, x for equal with cov2
## check taxon field
awk '{l=length($NF); for(i=1; i<=l; i++) {ss[i,substr($NF,i,1)]++; cs[substr($NF,i,1)]++}}
 END{ov=ms=bc=qx=0; d=""; for(i=1; i<=l; i++) {c="-";n=k=m=0; for(x in cs)
 {if(x!~/[-~nx]/ && ss[i,x]>0) {k++; m=m+ss[i,x]; if(ss[i,x]>n) {n=ss[i,x]; c=x}}}
 if(m>1) {ov++; if(k>1) ms++} if(m==2 && k==2) {c="x";qx++} d=d""c; if(m>0) bc++}
 printf "%s %s %s %s %s %s %s\n","'$NAME'",'$RRF1',bc,ov,ms,qx,d}' tdata.txt > cons.txt

awk '{printf "con%s.%s.%s.%s.%s_%s\t%s\n",$2,$3,$4,$5,$6,$1,$7}' cons.txt >> ${ASB}_cons.txt

awk '{print "'$i'",$1,$2,$3,$4,$5,$6}' cons.txt | tee -a ${ASB}_constats.txt

else
## Trinomial output (drop info field)
awk '{printf "%s_%s\t%s\n",$1,$2,$NF}' tdata.txt >> ${ASB}_cons.txt
fi

done


wc -l ${ASB}_cons.txt
wc -l ${ASB}_constats.txt
awk '{printf ">%s\n%s\n",$1,$2}' ${ASB}_cons.txt > ${ASB}_cons.fasta
sort -k7,7n ${ASB}_constats.txt

## Options for excel formats (all | or |+seq)
# grep 'con' ${ASB}_cons.txt | sed 's/-/|/g' > check.txt
# grep -G '^con' ${ASB}_cons.txt > check.txt


awk '{printf "%s\t%s\t|%s\n",NR,$1,$2}' ${ASB}_cons.txt > data.txt
awk '{print $2}' data.txt | sed 's/_/\t/g' > lf.txt
paste data.txt lf.txt > ${ASB}_exdb.txt



##### Drop conflicting consenses WORKING!!
## limit on x sites
## Assumes 3-fields: accession, taxon, seq
## Replace with longest individual accession

cat ${ASB}_cons.txt | sed 's/_/\t/1' > check.txt
wc -l check.txt

awk '{if($7>50) print $2}' ${ASB}_constats.txt > badcons.txt

grep -w -v -f badcons.txt check.txt > check2.txt

grep -w -f badcons.txt asb.txt > check3.txt
# grep -w -f badcons.txt check.txt > check3c.txt

awk '{b=length(gensub(/[-~nNxX?]/,"","g",$NF)); if(b>l[$2]) {l[$2]=b; r[$2]=$1"\t"$2"\t"$NF}}
 END{for(x in r) print r[x]}' check3.txt > check3b.txt

cat check2.txt check3b.txt | sed 's/\t/_/1' > check4.txt

wc -l badcons.txt
wc -l check.txt
wc -l check3.txt
wc -l check3b.txt
wc -l check4.txt


awk '{printf ">%s\n%s\n",$1,$2}' check4.txt > ${ASB}_cons4.fasta
awk '{printf "%s\t%s\t|%s\n",NR,$1,$2}' check4.txt > data.txt
awk '{print $2}' data.txt | sed 's/_/\t/g' > lf.txt
paste data.txt lf.txt > ${ASB}_exdb4.txt






######################################################################
###############     SUPERMATIX PROCESSOR  version 1    ###############
## DB format: label, loci field seqs (null=0), n-nuke, n-mito, n-bases
## load elist, create dummy seqs (~), process DB file, limit by loci/bases
## DB fields = elist rows-3 (label, n-loci, n-bases)


NAME=BEE_mat3

cat BEE_mat3_elist.txt > elist.txt

## Check/filter DB - all rows have same fields
awk '{s[NF]++} END{for(x in s) print x,s[x]; print NR,NF}' ${NAME}_DB.txt


### Format DB remove \r
## Option: single pipe | before sequence (new excel format)
sed 's/|//g' ${NAME}_DB.txt | tr -d '\r' > DB.txt

wc -l ${NAME}_DB.txt
wc -l DB.txt


### Assemble supermatrix - set taxon data limits - option to convert case
## currently set to include taxa by mimimum total bases (100)
awk 'NR==FNR{ne++; el[ne]=$2; for(n=1;n<=$2;n++) de[ne]=de[ne]"~"; next} {if($NF>100)
 {printf "%s ",$1; for(f=1;f<=ne;f++) {if($(f+1)==0) printf de[f]; else printf $(f+1)} printf "\n"}}' elist.txt DB.txt > SM.txt

### Version to include multiple criteria
awk 'NR==FNR{ne++; el[ne]=$2; for(n=1;n<=$2;n++) de[ne]=de[ne]"~"; next} {if(($(NF-2)+$(NF-1))>1 || $NF>600)
 {printf "%s ",$1; for(f=1;f<=ne;f++) {if($(f+1)==0) printf de[f]; else printf $(f+1)} printf "\n"}}' elist.txt DB.txt > SM.txt

## checks
awk '{s[length($2)]++} END{for(x in s) print x,s[x]; print NR,NF}' SM.txt
# awk '{print $1,length($2)}' SM.txt
## Find bads
awk '{if(length($2)!=35071) print $0}' SM.txt > bads.txt


### Output fasta and phylip
## check missing seq state char (n cannot be used for protein partitions)
# awk '{printf ">%s\n%s\n",$1,$2}' SM.txt > ${NAME}.fasta

### substitute n or x (AA)
sed 's/~/n/g' SM.txt > data.txt
awk 'END{print NR, length($2)}' data.txt > dim.txt
cat dim.txt data.txt > ${NAME}.phy




#### TRIM MARKLINE using elist base group list (fifth field) NEW VERSION
## Use pre defined column 5 in elist file, case insensitive (eb[ne]=$5, or preset eb[ne]=1)
## Set missing taxa percentage limit (MS)
## load elist, process sample: ignore null, process loci - site completeness, record markline (each site)

MS=10

## awk 'NR==FNR{el[++ne]=$2;eb[ne]=$5; next}
## {for(f=1;f<=ne;f++) {if($(f+1)==0) continue; else {te[f]++; for(s=1; s<=el[f]; s++) if(substr($(f+1),s,1)~/[acgtACGT]/) c[f,s]++}}}
## END{for(f=1; f<=ne; f++) {q=te[f]*'$MS'/100; b=eb[f]; for(s=1; s<el[f]; s=s+b)
## {p=m=0; for(n=1;n<=b;n++) if(c[f,s+n-1]>q) p++; if(p==b) m=1; for(n=1;n<=b;n++) print m}}}'  elist.txt DB.txt > marklineN.txt


### version 2 with block remainder site fix (0) - UNDER DEVELEOPMENT
awk 'NR==FNR{el[++ne]=$2;eb[ne]=$5; next}
 {for(f=1;f<=ne;f++) {if($(f+1)==0) continue; else {te[f]++; for(s=1; s<=el[f]; s++) if(substr($(f+1),s,1)~/[acgtACGT]/) c[f,s]++}}}
 END{for(f=1; f<=ne; f++) {q=te[f]*'$MS'/100; b=eb[f]; for(s=1; s<el[f]; s=s+b)
 {p=m=0; for(n=1;n<=b;n++) if(c[f,s+n-1]>q) p++; if(p==b) m=1; for(n=1;n<=b;n++) print m}
 r=el[f]%b; for(n=1;n<=r;n++) print 0}}'  elist.txt DB.txt > marklineN.txt


## checks
awk '{m[$1]++} END{for(x in m) print x,m[x]}' marklineN.txt
awk 'BEGIN{printf ">marklineMS%s\n",'$MS'} {printf $1} END{printf "\n\n"}' marklineN.txt > markline.fasta
awk 'NR==FNR{m[NR]=$1; next} {for(s=$3; s<=$4; s++) c[m[s]]++; print $0,c[0]+0,c[1]+0; c[0]=c[1]=0}' marklineN.txt elist.txt


### Output fasta and phylip
awk '{printf ">%s\n%s\n",$1,$2}' SM.txt > ${NAME}.fasta
cat markline.fasta >> ${NAME}.fasta


##### APPLY MARKLINE
## Revise elist
## Trim datmatrix
## Two-field data
## Simple substr function
## All sites marklineN version, sites as fields

## Generate new elist
awk 'NR==FNR{mc[NR]=$1;next}
 {b=0; for(s=$3;s<=$4;s++) {if(mc[s]==1) b++}
 if(b>0) {print($1,b,t+1,t+b); t=t+b}}' marklineN.txt elist.txt > ${NAME}_${MS}TCT_elist.txt


## Sites as field setup
awk '{print $1}' data.txt > www.txt
awk '{print $2}' data.txt > xxx.txt

wc -l www.txt
wc -l xxx.txt

awk 'END{print NR,length($2)}' data.txt
awk 'END{print NR,length($1)}' xxx.txt
awk 'END{print NR}' marklineN.txt


## Trim the data
awk 'BEGIN{FS=""} NR==FNR{c[NR]=$1;next} {for(s=1;s<=NF;s++)
 {if(c[s]==1) printf $s} printf "\n"}' marklineN.txt xxx.txt > yyy.txt

awk 'END{print NR,length($1)}' yyy.txt

paste -d" " www.txt yyy.txt > ${NAME}_${MS}TCT.txt

## checks
awk '{s[length($2)]++} END{for(x in s) print x,s[x]}' ${NAME}_${MS}TCT.txt


awk 'END{print NR, length($2)}' ${NAME}_${MS}TCT.txt > dim.txt
cat dim.txt ${NAME}_${MS}TCT.txt > ${NAME}_${MS}TCT.phy


### Tagline and consensus dummy (ignore [-~n])
awk 'BEGIN{print ">SMtagline"} {d=$2; n=($2-length($1)-1); if(n>2) {printf "|%s",$1;d=n} for(i=1;i<=d;i++) printf "~"}
 END{printf "\n"}' ${NAME}_${MS}TCT_elist.txt > tagline.fasta

awk '{print length($1)}' tagline.fasta

awk 'BEGIN{FS=""} {for(i=1; i<=NF; i++) {ss[i,$i]++; cs[$i]++}}
 END{print ">dummy"; for(i=1; i<=NF; i++) {c="-";n=0; for(x in cs)
 {if(x!~/[-~nxX]/ && ss[i,x]>n) {n=ss[i,x]; c=x}} printf c} printf "\n"}' yyy.txt > dummyT.fasta

awk '{printf ">%s\n%s\n",$1,$2}' ${NAME}_${MS}TCT.txt > zzz.txt
cat dummyT.fasta tagline.fasta zzz.txt > ${NAME}_${MS}TCT.fasta



#####  TAXON BY PARTITION P-SITES SUMMARY
## Load partition file, process alignmentS

## Option to refer to trimmed data - change input NAME
NAME=BEE_mat3_10TCT
tail -n +2 ${NAME}.phy > data.txt
cat ${NAME}_elist.txt > elist.txt



awk 'NR==FNR{ne++; el[ne]=$2; es[ne]=$3; next} {tb=0;printf $1; for(n=1;n<=ne;n++)
 {q=substr($2,es[n],el[n]); b=length(gensub(/[-~nNxX?]/,"","g",q)); tb=tb+b;\
 printf " %0.2f",b/el[n]} printf " %s %s\n",tb,length($2)}' elist.txt data.txt > ${NAME}_summ.txt

## checks
awk 'END{print NR,NF}' ${NAME}_summ.txt
awk '{s[length($2)]++} END{for(x in s) print x,s[x]}' data.txt




################    GENUS DATA OVERLAP
## Identify species in genus that has no data in common with other congenerics
## Use DB.txt
## Gene fields: 3 to NF-3
## load in arrays: count genus n per locus, count species in genus per locus
## summarize species in multi-species genera for number of genes in common with rest of genus
## for each multi-typic genus, for each species count number of genes score>1
## Option to filter by data limits if(($(NF-2)+$(NF-1))>1), if($NF>100)

## checks
awk '{s[NF]++} END{for(x in s) print x,s[x]; print NR,NF}' DB.txt


sed 's/_/\t/1' DB.txt > DBx.txt

wc -l DB.txt
wc -l DBx.txt

### Multi-taxon genera, version 1 
## six field output: genus, taxa, tot-genes, n-taxa, label, genes-in-common
awk '{if($NF>100) {g[$1]++; s[$1,g[$1]]=$1"_"$2; for(n=3;n<=(NF-3);n++)
 {if($n==0) sg[s[$1,g[$1]],n]=0; else {sg[s[$1,g[$1]],n]=1; gg[$1,n]++}}}}
  END{for(x in g) {if(g[x]>1) {for(n=1;n<=g[x];n++) {p=q=0; for(j=3;j<=(NF-3);j++)
  {if(sg[s[x,n],j]*gg[x,j]>1) q++; if(gg[x,j]>0) p++} print x,g[x],p,n,s[x,n],q}}}}' DBx.txt > DBx_GDC.txt


awk '{if($6<1) print $0}' DBx_GDC.txt > DBx_GDC_summ.txt
wc -l DBx_GDC.txt
wc -l DBx_GDC_summ.txt











############################################################################
###########################        R CODE        ###########################
## Needs packages APE, PHYTOOLS


##########   JEWEL TREE PARSING

##### PRUNE TAXA
## Simple list of taxa to be pruned out
dtax <- read.table("XMAS_mat2b_dlist.txt", as.is = TRUE)
dtaxv <- as.vector(dtax[,1])
length(dtaxv)
td=length(dtaxv)

BSTree <- read.nexus("XMAS_mat2b_bst2f_set3nn.trees")
length(BSTree)
gt=length(BSTree)

plot.phylo(ladderize(BSTree[[11]], right = FALSE), font=1, cex=0.5, no.margin=FALSE)
axisPhylo(side = 1, cex.axis=0.7)

plot.phylo(ladderize(BSTree, right = FALSE), font=1, cex=0.3, no.margin=FALSE)
axisPhylo(side = 1, cex.axis=0.7)
tree <- BSTree


#####
### Code tip labels from list
tree <- BSTree[[10]]
ne=length(tree$edge.length)
nt=length(tree$tip.label)

aa <- matrix(NA, nt, 1)
aa[,1]=0
names(aa) <- tree$tip.label

for (i in 1:td) {
aa[which(names(aa)==dtaxv[i]),1]=1
}

plot.phylo(ladderize(tree, right = FALSE), font=1, cex=0.5, show.tip.label=TRUE, no.margin=FALSE, tip.col=aa+1)
axisPhylo(side = 1, cex.axis=0.7)



#### Process trees - prune first, then add increment
print(c(gt,nt,ne,td))

for (i in 1:gt) {
tree <- BSTree[[i]]

ptree <-drop.tip(tree, dtaxv)
write.tree(ptree, file="XMAS_mat2b_bst2f_set3nn_pt.nwk", append=TRUE)

## note realculation of nt and ne
ne=length(ptree$edge.length)
nt=length(ptree$tip.label)
### Tip increment code, save second lot of trees
for (j in 1:ne) {
if(ptree$edge[j,2]<=nt) ptree$edge.length[j]=ptree$edge.length[j]+0.1
}
write.tree(ptree, file="XMAS_mat2b_bst2f_set3nn_pinct.nwk", append=TRUE)
}


plot.phylo(ladderize(ptree, right = FALSE), font=1, cex=0.6, no.margin=FALSE)

is.binary(ptree)
is.ultrametric(ptree)
sum(ptree$edge.length)

