#!/bin/bash
#PBS -V

while getopts g:p:n:o:a:q:m:e:x:k: option
do
  case "${option}"
    in
      g) geno_dir=${OPTARG};;
      p) p_loc=${OPTARG};;
      n) sname=${OPTARG};;
      o) outdir=${OPTARG};;
      a) acgt=${OPTARG};;
      q) geno_filter=${OPTARG};;
      m) maf_filter=${OPTARG};;
      x) lisa=${OPTARG};;
      e) extract=${OPTARG};;
      k) keep=${OPTARG};;
    esac
done

 if [ $lisa != "yes" ]
 then
  TMPDIR=$lisa
 fi
  if [ $TMPDIR == "" ]
 then 
  exit 110 #FAIL if tmpdir not set!"
 fi
 
 if [ $acgt == "yes" ] 
 then
  acgtflag="--snps-only just-acgt"
 fi
     
 if [ $extract != "xxxx" ] 
 then
  extractflag="--extract $extract"
 fi
     
   if [ $keep != "xxxx" ] 
 then
  keepflag="--keep $keep"
 fi
        
 chr="$PBS_ARRAYID"

 echo $outdir
 ls $geno_dir | grep _chr"$chr".fam | sed 's/.fam//g'   | awk '{print $1".bed",$1".bim",$1".fam"}' > "$outdir"/"$sname"_"$chr".mergelist
 mlist="$outdir"/"$sname"_"$chr".mergelist
 for files in $(ls $geno_dir | grep chr"$chr". ) 
 do
  cp "$geno_dir"/"$files" "$TMPDIR"/.
 done

 cp "$outdir"/"$sname"_"$chr".mergelist "$TMPDIR"/.

 cd "$TMPDIR"
 for files in $(ls | grep .gz)
 do 
  gzip -d $files
 done
 
 
 #Merge data sets
 $p_loc --allow-no-sex --merge-list "$outdir"/"$sname"_"$chr".mergelist --make-bed $acgtflag --out "$TMPDIR"/ZZ"$sname"_"$chr" 
 mv "$TMPDIR"/ZZ"$sname"_"$chr".log "$TMPDIR"/"$sname"_"$chr".log2
 #Quality filter data
 $p_loc --bfile "$TMPDIR"/ZZ"$sname"_"$chr"  --geno $geno_filter --maf $maf_filter  $extractflag $keepflag --make-bed --out "$TMPDIR"/"$sname"_"$chr"
 

 cp "$TMPDIR"/"$sname"_"$chr".* "$outdir"/.
 #rm "$TMPDIR"/*
