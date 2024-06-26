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

## Rules run by the pipeline (in order)

merge_files (merge_10x_gzip_files.sh)<br/>
rename_fastqs<br/>
trick_fastqs (cellranger_dna_trick.py)<br/>
move_fastqs<br/>
run_cellranger<br/>
create_raw_symlinks<br/>
copy_cellranger_outputs<br/>
extract_genomic_info<br/>
remove_tenx_artifacts<br/>
create_raw_symlinks<br/>
create_raw_files_list<br/>
create_raw_checksum<br/>
--- breakpoints detection<br/>
detect_breakpoints (sc-dna/bin)<br/>
segment_regions (sc-dna/bin)<br/>
normalise_counts<br/>
---- tree learning<br/>
clustering (secondary_analisys.apply_phenograph)<br/>
create_averaged_region_matrix<br/>
learn_empty_tree<br/>
learn_cluster_trees<br/>
pick_best_tree<br/>
cluster_tree_robustness<br/>
learn_nu_on_cluster_tree<br/>
---- plotting<br/>
plot_cnv_matrix<br/>
cell_assignment<br/>
plot_cluster_cnvs<br/>
plot_cnv_matrix<br/>
full_tree_robustness<br/>
---- process cnvs<br/>
create_bin_gene_region_df<br/>
create_cn_cluster_h5<br/>
---- process cnvs<br/>
visualise_trees<br/>
plot_heatmaps<br/>
create_derived_symlinks<br/>
create_derived_checksum<br/>
all<br/>

## How to run the pipeline

1.  **Register the sc-dna stage on** [**LabKey**](https://tp-labkey.ethz.ch/labkey/Tumor%20Profiler%20-%20Melanoma/project-begin.view):  
<code> Overview → scDNA → fill in form in "Register scDNA Analysis" section: First Sequencing Run ID, SOP version, Run Number → Submit </code>

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
$ python create_dnapipeline_config.py -f BSSE_QGF_131049_HNY5VBGXC_1_MADEGOD_T_scD_250c_r1v1_0_SI-GA-E9_S2_L001_I1_001 -t ovarian --is_novaseq --gender male -o /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_ovarian/OKEKIDE-T/singlecell_dna/snake_analysis_files/config.json
```
 
5.  **Run the pipeline**
```    
$ module load python/3.6.0 graphviz hdf5/1.8.12 gcc/6.2.0 eth_proxy cmake/3.11.4
$ module load /cluster/work/bewi/modules/cellranger-dna/1.1.0/cellranger-dna-1.1.0/cellranger-dna.modulefile
## configure SCICoNE
$ cd /cluster/work/bewi/ngs/projects/tumorProfiler/code/dna-pipeline/scicone
$ mkdir build && cd build                       # Create and enter the build directory
$ cmake ..                                      # Compile the program with cmake
$ make                                          # Build the executables
$ cd /cluster/work/bewi/ngs/projects/tumorProfiler/code/dna-pipeline/scicone/pyscicone
$ pip install -e .
## configure the pipeline
$ cd /cluster/work/bewi/ngs/projects/tumorProfiler/code/dna-pipeline/scdna-pipe/
$ pip install -e .

$ export OMP_NUM_THREADS=24
$ cd path-to/singlecell_dna/snake_analysis_files/
$ bsub -J OHAMUME -n 48 -W 23:57 -R fullnode -R "rusage[mem=5000]" snakemake -s /cluster/work/bewi/ngs/projects/tumorProfiler/code/dna-pipeline/scdna-pipe/pipeline/Snakefile --configfile ./config.json -j 48 -p -k
```

6.  **Write the results to the to_upload folder**
* $cd to_upload/
* check derived/summary.txt

7.  **Send the results to the team (via email / drive)**

8.  **Upload the results to LeonhardMed**
```    
$ ssh leomed
$ cd /cluster/work/tumorp/data-drop/scDNA/report
$ mkdir <sample_name>
$ get pwd ← target_path
-- go to Euler to .../to_upload
$ scp derived/ raw/ <NETHZ_USERNAME>@login.leomed.ethz.ch:target_path
-- go to target_path on LeoMed
$ touch done.txt (or `done_failed.txt` if failed w/ description)
$ chmod 2770 -R target_path
$ bash /cluster/work/tumorp/data-drop/check_data_format.sh target_path
```
