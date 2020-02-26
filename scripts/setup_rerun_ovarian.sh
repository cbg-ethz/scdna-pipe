declare -a novaseq=("OKEKIDE-T ODAMILU-T OHAMUHE-T")

for fastq_path in /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_ovarian/*/singlecell_dna/openbis/BSSE*L001_I1_001.fastq.gz;do
  filename=`basename $fastq_path`
  dirname="${fastq_path%/*/*/*}"
  sample_id=`basename $dirname`
  echo $sample_id
  output_dir=$dirname/singlecell_dna/snake_analysis_files
  output_file=$output_dir/config_v1.13.json
  is_novaseq=false
  for item in $novaseq
  do
      if [ "$sample_id" == "$item" ]; then
          is_novaseq=true
          break
      fi
  done
  if $is_novaseq; then
    python /cluster/work/bewi/ngs/projects/tumorProfiler/code/dna-pipeline/scdna-pipe/scripts/create_dnapipeline_config.py -f ${filename} -t ovarian -o ${output_file} --gender female --is_novaseq
  else
    python /cluster/work/bewi/ngs/projects/tumorProfiler/code/dna-pipeline/scdna-pipe/scripts/create_dnapipeline_config.py -f ${filename} -t ovarian -o ${output_file} --gender female
  fi
  cd $output_dir
  snakemake -s /cluster/work/bewi/ngs/projects/tumorProfiler/code/dna-pipeline/scdna-pipe/pipeline/Snakefile --configfile ./config_v1.13.json -j 48 -p -k --unlock
  snakemake -s /cluster/work/bewi/ngs/projects/tumorProfiler/code/dna-pipeline/scdna-pipe/pipeline/Snakefile --configfile ./config_v1.13.json -j 48 -p -k --touch
  bsub -J $sample_id -n 48 -W 23:57 -R fullnode -R "rusage[mem=5000]" snakemake -s /cluster/work/bewi/ngs/projects/tumorProfiler/code/dna-pipeline/scdna-pipe/pipeline/Snakefile --configfile ./config_v1.13.json -j 48 -p -k
done
