
 #List input directory
 indir= #Give full path of best guess genotype data, e.g. 
 /home/mark/tracts/qc1/imputation/dos_qc1.../bgn/
 outdir= #Give full path of where you want the merged data to go 
 plink_path= #Give full path of plink executable bgtype=bgn #arbitrary 
 suffix for filename
 
 mkdir $outdir
 ls "$indir" |  grep .bim$ |  sed 's/.bim//g' |   awk '{print $1".bed",$1".bim",$1".fam"}'  "$outdir"/mergelist_"$studyname"_"$bgtype".txt
 
 #Merge datasets
 $plink_path --merge-list "$outdir"/mergelist_"$studyname"_"$bgtype".txt --allow-no-sex --make-bed --out "$outdir"/"$studyname"_"$bgtype"
 



#Step 1: If there are multiple datasets, merge all genotype datasets into one file

##Important!!: each input genotype dataset should be named e.g. studyname_chromsome, e.g. duke_1 would be chromosome 1 from the duke dataset. This is necessary for a future loop code..

#Set working directory
 wd=/home/maihofer/raj
#Call into it
 cd $wd

#Set Plink binary location
 p_loc=/home/maihofer/katy/plink
 
#Name that merged genotype data will be called
 sname=subfields
#Location of genotype data for all studies to be merged. Data should be split by chromosome, with file naming as above
 starting_geno_dir="$wd"/bfiles2/
#Location where merged genotypes will go (by chromosome)
 outdir_mgeno="$wd"/merged_genotypes
#Location where merged genotypes will go (all chromosomes combined)
 outdir_mgeno2="$wd"/merged_genotypes_allchr
#Location of temporary files folder. Everything in here will be deleted!!
 TMPDIR=

 if [ ! -e "errandout" ]
 then
  mkdir errandout
 fi
 
 if [ ! -e $outdir_mgeno ]
 then
  mkdir $outdir_mgeno
 fi
 if [ ! -e $outdir_mgeno2 ]
 then
  mkdir $outdir_mgeno2
 fi
 if [ ! -e $TMPDIR ]
 then
  mkdir -p $TMPDIR
 fi
 
#only retain markers with missingness < this value
 gtthresh=0.05
#maf: filter merged genotyped data on this observed MAF 
 mafthresh=0.005

#List of SNPS with global maf > 1%, data will be filtered to these.
 filterfile=/home/maihofer/katy/allmaf01_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf

#Merge all datasets, by chromosome
 geno_dir=$starting_geno_dir;
 outdir=$outdir_mgeno
 acgt=yes #By default we'll filter to just acgt SNPs
 geno_filter=$gtthresh
 maf_filter=$mafthresh
 extract=$filterfile
 
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
        
 for chr in {1..22}
 do
  echo "Will output merged genotypes to $outdir"
  ls $geno_dir | grep _chr"$chr".fam | sed 's/.fam//g'   | awk '{print $1".bed",$1".bim",$1".fam"}' > "$outdir"/"$sname"_"$chr".mergelist
  mlist="$outdir"/"$sname"_"$chr".mergelist
  nfcheck= $(wc -l $mlist | awk '{print $1}')
  echo "Will merge $nfcheck files on chromosome $chr. See file $mlist for details"
  
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
  
  echo "Copying files "$TMPDIR"/"$sname"_"$chr" into "$outdir""
  cp "$TMPDIR"/"$sname"_"$chr".* "$outdir"/.
 #rm "$TMPDIR"/*
 done

#Merge all chromosomes into one dataset
 ls $outdir_mgeno/* | grep .bed | sed 's/.bed//g' | awk '{print $1".bed",$1".bim",$1".fam"}' > allchr.mergelist
 $p_loc --allow-no-sex --merge-list "$wd"/allchr.mergelist --make-bed --out "$outdir_mgeno2"/"$sname" 
