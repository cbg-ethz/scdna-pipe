# Get Started with Tumor Profiler 

## Access needed

-   OpenBIS (where the sequencing data is uploaded): email to [vincenzo.capece@id.ethz.ch](mailto:vincenzo.capece@id.ethz.ch)
-   LabKey (where we register the samples): [labkey@nexus.ethz.ch](mailto:labkey@nexus.ethz.ch)
-   LeonhardMed (where the results are uploaded) - contact person: [vipin.sreedharan@nexus.ethz.ch](mailto:vipin.sreedharan@nexus.ethz.ch)
-   This also gives access to [tu-pro.ch](http://tu-pro.ch)
-   Here you need to be added to the group INFK-Raetsch-tp-beerenwinkellab
```
login: ssh -J <NETHZ_USERNAME>@jump.leomed.ethz.ch:22 <NETHZ_USERNAME>@login.leomed.ethz.ch
```
-   Invite to Slack TumorProfiler workspace
    

## Folder structure explained
```
|-- /cluster
  |-- work
    |-- bewi
	  |-- ngs
	    |-- projects
		  |-- tumorProfiler
		    |-- code
			  |-- bulkRNA ← [git: cbg-ethz/single-cell-tumor-profiler]				
			  |-- dna_pipeline ← [git: cbg-ethz/single-cell-tumor-profiler]
			    |-- data ← ref genome, drug list, etc
				  |-- analysis					    
				    |-- trial_melanoma or trial_ovarian
				      |-- MADEGOD-T ← directory name: sample id + suffix
					    |-- singlecell_rna
					    |-- singlecell_dna
						  |-- analysis  
						  |-- fastqs ← ??? empty
						  |-- final_sftp
						  |-- metadata
						  |-- openbis ← fastq.gz files
						  |-- snake_analysis_files
						    |-- config.json
						  |-- snake_temp ← ??? empty
```   

## How to run the pipeline

  

1.  **Register the sc-dna stage on** [**LabKey**](https://tp-labkey.ethz.ch/labkey/Tumor%20Profiler%20-%20Melanoma/project-begin.view):  
<code> Overview → scDNA → fill in form: First Sequencing Run ID, SOP version, Run Number → Submit </code>
 
2.  **Setup the directory structure:**
```
$ cd /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/
$ mkdir <sample_name>
$ cp /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/exampleFolder <sample_name>/.
$ mkdir /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/<sample_name>/singlecell_dna/to_upload
```
---> or run the script in [scripts/setup_directory_structure.sh](https://github.com/cbg-ethz/scdna-pipe/blob/master/scripts/setup_directory_structure.sh "setup_directory_structure.sh")
```
$ ./setup_directory_structure.sh /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/ MADEGOD
```

3.  **Download the sequencing data from OpenBIS** to .../singlecell_dna/openbis
```   
$ lftp sftp://USERNAME@bsse.ethz.ch@bs-openbis04.ethz.ch:2222
$ cd BSSE_BEERENWINKEL_TPROFILER/BEERENWINKEL_TPROFILER/<OPENBIS_EXPERIMENT>/<OPENBIS_DATASET>/original/BSSE_QGF_131047_HNTM7BGXC_1
$ `mget *.gz` to get the fastqs
```

4.  **Create config.json** in .../singlecell_dna/snake_analysis_files/.
---> change sequencing_prefix to Sequencing Run ID from LabKey and sample name or run [script/create_dnapipeline_config.py](https://github.com/cbg-ethz/scdna-pipe/blob/master/scripts/create_dnapipeline_config.py "create_dnapipeline_config.py")
```
$ create_dnapipeline_config.py -s OKEKIDE-T_scD_250c-r1v1.0_r1v1.0-HK5HMDRXX -t ovarian --is_novaseq -o /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_ovarian/OKEKIDE-T/singlecell_dna/snake_analysis_files/config.json
```

5.  **Run the pipeline**
```    
$ module load python/3.6.0 graphviz hdf5/1.8.12 gcc/6.2.0 eth_proxy
$ module load /cluster/work/bewi/modules/cellranger-dna/1.1.0/cellranger-dna-1.1.0/cellranger-dna.modulefile

$ cd path-to/singlecell_dna/snake_analysis_files/
$ bsub -J MADEGOD -n 48 -W 23:57 -R fullnode -R "rusage[mem=5000]" snakemake -s /cluster/work/bewi/ngs/projects/tumorProfiler/code/dna-pipeline/pipeline/Snakefile --configfile ./config.json -j 48 -p -k
```

6.  **Write the results to the to_upload folder**
* $cd to_upload/
* check derived/summary.txt

7.  **Send the results to the team (via email / drive)**
   
8.  **Upload the results to LeonhardMed**
```    
$ ssh -J <NETHZ_USERNAME>@jump.leomed.ethz.ch:22 <NETHZ_USERNAME>@login.leomed.ethz.ch
$ cd /cluster/work/tumorp/data-drop/scDNA/report
$ mkdir <sample_name>
$ get pwd ← target_path
-- go to Euler to .../to_upload
$ scp derived/ raw/ <NETHZ_USERNAME>@login.leomed.ethz.ch:target_path
-- go to target_path on LeoMed
$ touch done.txt (or `done_failed.txt` if failed w/ description)
```
