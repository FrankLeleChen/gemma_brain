#Producing a mega analysis
#All SNPs with GMAF > 1% - take overall bg data # For each dataset, filter down to superset samples
#Can do this at dosage step without issue... maybe start now.. nah, just take the data then filter
#only need some subset of datasets
#Merge datasets in plink, as dosages, convert gemma
#Gemma - 
#have to convert raj stuff to ricopili ids

 wd=/home/maihofer/raj
#Call into it
 cd $wd

#Study name
 sname=subfields
#Location of genotype data
 starting_geno_dir="$wd"/bfiles2/
#Location of merged genotypes, by chr
 outdir_mgeno="$wd"/merged_genotypes
#Location of merged genotypes, all chr
 outdir_mgeno2="$wd"/merged_genotypes_allchr
  
#Plink binary location
 p_loc=/home/maihofer/katy/plink
 
#List MRI subject pool
 keep="$wd"/mri.subjects
#Working on lisa?
lisa=yes

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
 
#only retain markers genotyped in this % of people
 gtthresh=0.05 #May have to redo to set this
#maf: filter on this MAF
 mafthresh=0.01
#Filter file to these SNPs
 filterfile=/home/maihofer/katy/allmaf01_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
  
#Option 1. Combine all gemma files. REquires sort and join of all datasets.. dont bother this osunds like a shitshow
#Option 2. Join all in PLINK. Convert joined file to GEMMA format in AWK.
#Join all fam files?
#Can use dummy phenotype to insure it worked...

#Step 1: Marker GRMS

#Extract all genotypes
 for dataset in minv pts1 psy3 meg2 betr safr gtpc
 do
  tar xvf /archive/maihofer/dac/"$dataset"__cogbgfile_v1.tar --strip-components=5
  #ls *$dataset* | grep .fini | wc -l
 done
 
#Merge all datasets, into chromosome sized pieces
 qsub -t1-2 -d $wd -lwalltime=00:45:00 scripts/merge_datasets.sh -e errandout/ -o errandout/ -F "-g $starting_geno_dir -p $p_loc -n $sname -o $outdir_mgeno -a yes -q $gtthresh -m $mafthresh -e $filterfile -k $keep -x $lisa" 

#Merge all chromosomes into one dataset
 ls $outdir_mgeno/* | grep .bed | sed 's/.bed//g' | awk '{print $1".bed",$1".bim",$1".fam"}' > allchr.mergelist
 echo "module load plink2; plink --allow-no-sex --merge-list "$wd"/allchr.mergelist --make-bed --out "$outdir_mgeno2"/"$sname" " > scripts/merge_allchr.sh
 qsub -d $wd -lwalltime=01:00:00 -e errandout/ -o errandout/ scripts/merge_allchr.sh
 
 qsub -d $wd -lwalltime=00:05:00 -e errandout/ -o errandout/ scripts/merge_allchr.sh
