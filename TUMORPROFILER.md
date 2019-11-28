# Tumor Profiler workflow

* Get scD data from openBIS (link in email)
* Register sample on LabKey:
  * https://tp-labkey.ethz.ch/labkey/Tumor%20Profiler%20-%20Melanoma/project-begin.view?
  * Overview
  * scDNA
  * Fill in form: First Sequencing Run ID, SOP version, Run Number -> Submit
* go to /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/MIDEKOG-T/singlecell_dna
* prepare folder structure (snake_analysis_files, openbis)
* go to openbis/
* lftp sftp://pedro.ferreira@bsse.ethz.ch@bs-openbis04.ethz.ch:2222
* go to BSSE_BEERENWINKEL_TPROFILER/BEERENWINKEL_TPROFILER/<OPENBIS_EXPERIMENT>/<OPENBIS_DATASET>/original/BSSE_QGF_131047_HNTM7BGXC_1
* `mget *.gz` to get the fastqs
* copy config.json to snake_analysis_files/
  * change sequencing_prefix to Sequencing Run ID from LabKey and sample name

* Running
  * submit the pipeline: `bsub -J MOVAZYQ -n 48 -W 23:57 -R fullnode -R "rusage[mem=5000]" snakemake -s /cluster/work/bewi/ngs/projects/tumorProfiler/code/dna-pipeline/pipeline/Snakefile --configfile ./config.json -j 48 -p -k`

* Upload prep
  * cd to_upload/
  * check derived/summary.txt
  * send figures to team

* Upload
  * ssh leomed from local
  * cd /cluster/work/tumorp/data-drop/scDNA/report
  * mkdir <sample_name>
  * get pwd
  * go to to_upload on euler and run `scp derived/ raw/ pedrof@login.leomed.ethz.ch:<target_path>`
  * go to target_path on leomed and run `touch done.txt`
