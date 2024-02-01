## General introduction

Microbial community structure, diversity, and population structure are fundamental aspects to understand evolution, niche adaptation and demographic history of bacterial species. As high-throughput sequencing has become cost effective and accessible, sequencing populations of bacteria across the whole metagenome provides unprecedented resolution to investigate within-host evolution, transmission history and population structure. Moreover, analysis of genetic content of microbial communities through metagenomics has become the mainstream methodology. During this one-week course you will learn and apply bioinformatic techniques to perform population genetics and study microbial communities with metagenomic approaches. The goal is to become familiar with the bioinformatic analysis tools and to be able to utilize them in your own research after the course.

Everything is possible thanks to the support of CSC.

## Schedule

The course will be on Monday 5th of February, Wednesday 7th of February and Friday 9th of February in the City Centre Campus.

__Class rooms:__
| Day | Time | Building | Room |  
| --- | --- | --- | --- |  
| Monday 5.2.   | 10-16   | :classical_building: Main Building | __U3041__  |
| Wednesday 7.2.| 10-16   | :classical_building: Main Building | __U3043__  | 
| Friday 9.2.   | 10-16   | :classical_building: Main Building | __U3043__  |

__[Detailed schedule](Schedule.md)__  

## Practicals

Course practical part can be found from here: __[Practicals](Practicals/README.md)__


## Target group
This course is targeted at PhD and MSc students who are interested in performing bacterial genomic analyses

## Learning outcomes
* Foundational skills to work with metagenomic data
* Familiarity and practice with bioinformatics tools
* Perspective and confidence to apply these skills in your own work
* Empower you to ask and answer the questions you have of your own data

## Help with connections

Eduroam or Helsinki UniGuest networks are available

## Organizers and teachers
:woman_technologist: Docent Jenni Hultman, LUKE 
:man_technologist: University lecturer Antti Karkman, University of Helsinki  

# Practicals

__Table of Contents:__
1. [Setting up](#setting-up-the-course-folders)
2. [Interactive use of Puhti](#interactive-use-of-puhti)
3. [QC and trimming for Illumina reads](#qc-and-trimming-for-illumina-reads)
4. [QC and trimming for Nanopore reads](#qc-and-trimming-for-nanopore-reads)
5. [Genome assembly](#genome-assembly)
6. [Assembly QC](#assembly-qc)
7. [Calculate the genome coverage](#calculate-the-genome-coverage)
8. [Genome completeness and contamination](#genome-completeness-and-contamination)
9. [Genome annotation with Prokka](#genome-annotation-with-prokka)
10. [Name the strain](#name-the-strain)
11. [Pangenomics](#pangenomics-with-anvio)
12. [Detection  of secondary  metabolites biosynthesis gene clusters](#detection-of-secondary-metabolites-biosynthesis-gene-clusters)
13. [Comparison of secondary metabolites biosynthesis gene clusters](#comparison-of-secondary-metabolites-biosynthesis-gene-clusters)

## Setting up the course folders
The main course directory is located in `/scratch/project_2005590`.  
There you will set up your own directory where you will perform all the tasks for this course.  

First list all projects you're affiliated with in CSC.
```
csc-workspaces
```

You should see the course project `MBDP_genomics_2024`.
So let's create a folder for you inside the scratch folder, you can find the path in the output from the previous command.

```bash
cd /scratch/project_2005590
mkdir $USER
```

Check with `ls`; which folder did `mkdir $USER` create?  

Navigate to your own folder and copy the course Github repository there. 

```bash
git clone https://github.com/MBDP-bioinformatics-courses/MBDP_Genomics_2024.git
```

This directory (`/scratch/project_2005590/your-user-name/MBDP_Genomics_2024`) is your working directory.  
Every time you log into Puhti, you should use `cd` to navigate to this directory, and **all the scripts are to be run in this folder**.  

The raw data used on this course can be found in `/scratch/project_2005590/RAWDATA/{STRAIN}`. Where {STRAIN} is either `KLB3.1` or `WOD100`.   
Instead of copying the data we will use links to this folder in all of the needed tasks.  
Why don't we want 24 students copying data to their own folders?

Softlinks are made with `ln -s`. When you have chosen the strain you want to work with, make a softlink to all of the reads files under the folder `01_RAW_READS`.  
Remember to change the strain name in the command below.  

```bash
ln -s /scratch/project_2005590/RAWDATA/{STRAIN}/* 01_RAW_READS/
```

## Interactive use of Puhti

Puhti uses a scheduling system called SLURM. Most jobs are sent to the queue,  but smaller jobs can be run interactively.

Interactive session is launched with `sinteractive`   .   
You can specify the resources you need for you interactive work interactively with `sinteractive -i`. Or you can give them as options to `sinteractive`.  
You always need to specify the accounting project (`-A`, `--account`). Otherwise for small jobs you can use the default resources (see below).

| Option | Function | Default | Max |  
| --     | --       | --      | --  |  
| -i, --interactive | Set resources interactively |  |  |  
| -t,  --time | Reservation in minutes or in format d-hh:mm:ss | 24:00:00 | 7-00:00:00 |
| -m, --mem | Memory in Mb       | 2000     | 76000  |  
| -j, --jobname |Job name       | interactive     |   |  
| -c, --cores     | Number of cores       | 1      | 8  |  
| -A, --account     | Accounting project       |       |  |  
| -d, --tmp     | $TMPDIR size (in GiB)      |  32     | 760  |  
| -g, --gpu     | Number of GPUs       | 0     | 0 |  


[__Read more about interactive use of Puhti.__](https://docs.csc.fi/computing/running/interactive-usage/#sinteractive-in-puhti)   


## QC and trimming for Illumina reads
QC for the raw data takes few minutes, depending on the allocation.  
Go to your working directory and make a folder called `01_RAW_READS/FASTQC` for the QC reports of Illumina data.  

QC does not require lot of memory and can be run on the interactive nodes using `sinteractive`.

Activate the biokit environment and open interactive node:

```bash
sinteractive -A project_2005590
module load biokit
```

### Running FastQC
Run `FastQC` for the raw Illumina reads in the 01_RAW_READS folder. What does the `-o` and `-t` flags refer to?

```bash
fastqc path-to-R1-reads --outdir 01_RAW_READS/FASTQC
fastqc path-to-R2-reads --outdir 01_RAW_READS/FASTQC
```

```bash
module load multiqc
multiqc --interactive --outdir 01_RAW_READS/FASTQC/ 01_RAW_READS/FASTQC/*
```

Copy the resulting HTML file (`multiqc_report.html`) to your local machine. You can also copy the individual FastQC output HTML files.  
Have a look at the QC report(s) with your favorite browser.  

After inspecting the output, __what kind of trimming do you think should be done?__

### Running Cutadapt

The adapter sequences that you want to trim are specified with options `-a` and `-A`.  
What is the difference with `-a` and `-A`?  
And what is specified with option `-p` or `-o`?
And how about `-m` and `-j`?  
You can find the answers from Cutadapt [manual](http://cutadapt.readthedocs.io).

Remember to change the strains name and modify the paths to the raw read files.  

```bash
cutadapt \
    -a CTGTCTCTTATACACATCT \
    -A CTGTCTCTTATACACATCT \
    -o 02_TRIMMED_READS/{STRAIN}_1.fastq.gz \
    -p 02_TRIMMED_READS/{STRAIN}_2.fastq.gz \
    path-to-R1-file \
    path-to-R1-file \
    --minimum-length 80 \
    > 00_LOGS/cutadapt.log
```

### Running fastQC on the trimmed reads
You could now check the `cutadapt.log` and answer:

* How many read pairs we had originally?
* How many reads contained adapters?
* How many read pairs were removed because they were too short?
* How many base calls were quality-trimmed?
* Overall, what is the percentage of base pairs that were kept?

Then make a new folder (`FASTQC`) in the `02_TRIMMED_READS` folder for the QC files of the trimmed data and run fastQC and multiQC again as you did before trimming.  
Copy the resulting HTML file to your local machine and see how well the trimming went.  


## QC and trimming for Nanopore reads

The QC for the Nanopore reads can be done with NanoPlot and NanoQC. They are plotting tools for long read sequencing data and alignments. You can read more about them in: [NanoPlot](https://github.com/wdecoster/NanoPlot) and [NanoQC](https://github.com/wdecoster/nanoQC)

This run will require more computing resources, so you can apply for more memory or run as a batch job:

First log out from the computing node (if you're still on one).  

```bash
exit
```

Then open a new interactive task with more memory

```bash
sinteractive -A project_2005590 -m 45000 -c 4
```
NanoPlot and NanoQC are not pre-installed to Puhti so we need to reset the modules and activate the virtual environment. If the environment is already loaded you can skip this step.

Generate graphs for visualization of reads quality and length distribution

```bash
/projappl/project_2005590/nanotools/bin/NanoPlot \
    --outdir 01_RAW_READS/NanoPlot \
    --threads 4 \
    --format png \
    --fastq path-to/your_raw_nanopore_reads.fastq.gz
```

Transfer to your computer and check two plots inside the nanoplot output folder:
Reads quality distribution: `LengthvsQualityScatterPlot_kde.png`
Reads length distribution: `Non_weightedLogTransformed_HistogramReadlength.png`

```bash
/projappl/project_2005590/nanotools/bin/nanoQC \
    --outdir 01_RAW_READS/nanoQC \
    path-to/your_raw_nanopore_reads.fastq.gz
```

Using the Puhti interactive mode, check the file `nanoQC.html` inside the ouput folder of the nanoQC job.

* How is the quality at the beginning and at the end of the reads? How many bases would you cut from these regions?



### Trimming and quality filtering of reads

We'll use a program called [chopper](https://github.com/wdecoster/chopper) for quality filtering and trimming.  

The following command will trim the first XX bases and the last YY bases of each read, exclude reads with a phred score below ZZ and exclude reads with less than XYZ bp.

```bash
gunzip -c path-to/your_raw_nanopore_reads.fastq.gz |\
    /projappl/project_2005590/nanotools/bin/chopper -q ZZ -l XYZ --headcrop XX --tailcrop YY |\
    gzip > 02_TRIMMED_READS/{STRAIN}_nanopore.fastq.gz
```

### Optional - Visualizing the trimmed data

Use NanoPlot to see how the trimming worked out.  

```bash
/projappl/project_2005590/nanotools/bin/NanoPlot ... 
```

## Genome assembly 

Now that you have good trimmed sequences, we can assemble the reads. For assembling you will need more resources than the default.  
Allocate 4 cpus, 40000 Mb of memory (40G) and 2 hours. Remember also the accounting project, `project_2005590`.  

```bash
sinteractive --account --time --mem --cores
```

### Nanopore only assembly with FLye

```bash
/projappl/project_2005590/flye/bin/flye --nano-hq 02_TRIMMED_READS/WOD100_nanopore.fastq.gz --out-dir 03_ASSEMBLIES/flye --threads $SLURM_CPUS_PER_TASK
 ```

### Illumina only assembly with spades

```bash
module purge
module load spades/3.15.0
spades.py -1 02_TRIMMED_READS/WOD100_1.fastq.gz -2 02_TRIMMED_READS/WOD100_2.fastq.gz -o 03_ASSEMBLIES/spades -t $SLURM_CPUS_PER_TASK --isolate
```

### Hybrid assembly with Unicycler

For the hybrid assembly (uses both long- and short-reads) we will use [Unicycler](https://github.com/rrwick/Unicycler). It can assemble short-reads, long-reads or both for hybrid assembly.  
It is a bit out-dated, but might still be a good option for short-read-first hybrid assembly, when the long-read sequencing depth is not optimal. The best option would be to use long-read-first approach. For that there is [Trycycler](https://github.com/rrwick/Trycycler).   

Unicycler has three different modes; conservative, normal and bold.  Conservative will produce the least misassemblies, but is the least likely to produce a cpmplete assembly. Bold is the most likely to produce a complete assembly with the risk of misassemblies. 
Normal is the default, but you can change it if you feel bold or conservative today. Or keep it normal, even if you don't feel totally normal.  

```bash
 /projappl/project_2005590/unicycler/bin/unicycler \
    -l 02_TRIMMED_READS/WOD100_nanopore.fastq.gz \
    -1 02_TRIMMED_READS/WOD100_1.fastq.gz \
    -2 02_TRIMMED_READS/WOD100_2.fastq.gz \
    --keep 0 \
    --mode normal \
    --out 03_ASSEMBLIES/WOD100_unicycler \
    --threads $SLURM_CPUS_PER_TASK
```

After you're done, remember to close the interactive connection and free the resources with `exit`.  
And to make things a bit easier for some of the next steps, copy each of the assembly files to the `03_ASSEMBLIES` folder.  
There's one example, but do it for each of the three assemblies. And remember to nme them so that you know which is which.  

```bash
cp 03_ASSEMBLIES/KLB_flye/assembly.fasta 03_ASSEMBLY/flye_assembly.fasta
```

## Assembly graphs

Each of the assemblers produces also an assembly graph. Download the assembly graphs (`.gfa` or `.fastg`) for each assembly to your local computer and open them with [Bandage](https://rrwick.github.io/Bandage/).  

## Assembly QC

After the assemblies are we will use Quality Assessment Tool for Genome Assemblies, [Quast](http://quast.sourceforge.net/) for (comparing and) evaluating our assemblies.

```bash
module purge
module load quast
quast.py --output-dir 03_ASSEMBLIES/QUAST 03_ASSEMBLIES/*.fasta
```

Now you can move the file `03_ASSEMBLIES/QUAST/report.html` to your computer and compare the different assemblies.  

## Polishing the nanopore assembly (OPTIONAL)

Although Flye does some polishing at the end of the assembly, we can try to polish our nanopore assembly with the nanopore reads using [medaka](https://github.com/nanoporetech/medaka).  
There are also many other polishing methods, you could also use the short-reads to polish the nanopore assembly. For medaka we need to choose a model based on how the basecalling was done.  
We will use the model `r1041_e82_400bps_hac_g632`. 

```bash
module purge
module load medaka
medaka_consensus -i path-to/trimmed_nanopore.fastq.gz -d path-to/flye_assembly.fasta -o 03_ASSEMBLIES/{STRAIN}_polished --threads $SLURM_CPUS_PER_TASK -m r1041_e82_400bps_hac_g632
```

Copy also the polished assembly (`consensus.fasta`) to the `03_ASSEMBLIES` folder. 

## Genome completeness and contamination

Now we have calculated different metrics for our genomes with QUAST and also done some polishing for the nanopore assembly, but we still don't know the "real" quality of our genome.  
One way to get a better sense is to estimate the completeness of the assemblies.

We will use [CheckM2](https://github.com/chklovski/CheckM2) to calculate the completeness and possible contamination in our genome.  
Allocate some resources (>40G memory & 4 threads) and run CheckM2.

Then run CheckM2 with the following command. 

```bash
/projappl/project_2005590/tax_tools/bin/checkm2 predict \
    --input 03_ASSEMBLIES \
    --output-directory 03_ASSEMBLIES/CheckM2 \
    --extension .fasta \
    --ttable 11 \
    --threads $SLURM_CPUS_PER_TASK \
    --database_path /scratch/project_2005590/DB/CheckM2/CheckM2_database/uniref100.KO.1.dmnd 
```

Hopefully now we have come to some conclusion about what is the best assembly and from now on we will continue with only one. 

## Mapping reads to the assembly and calculating the genome coverage (OPTIONAL)

To calculate the genome coverage, all the reads used for the assembly must be mapped to the final genome. As an example we will only map the short reads against the genome.  
For that, we use three programs: Bowtie2 to map the reads; Samtools to sort and make an index of the mapped reads; and bedtools to make the calculation. For long reads the process woul dbe the same, except you would need to use some othere mapping software (e.g. minimap2). 

First build an index from the assembly.

```bash
module load bowtie2
bowtie2-build 03_ASSEMBLIES/KLB3.1_unicycler.fasta 03_ASSEMBLIES/KLB3.1
```

Then map the short reads against the index.

```bash
bowtie2 -x 03_ASSEMBLIES/KLB3.1 -1 02_TRIMMED_READS/KLB3.1_1.fastq.gz -2 02_TRIMMED_READS/KLB3.1_2.fastq.gz --threads $SLURM_CPUS_PER_TASK -S 05_MAPPING/KLB3.1.sam
```

Then we process the resulting alignment file (`.sam`) and produce a sported and compressed format of that (`.bam`)

```bash
module purge
module load samtools
samtools view -Sb 05_MAPPING/KLB3.1.sam |samtools sort > 05_MAPPING/KLB3.1.bam
samtools index 05_MAPPING/KLB3.1.bam
```

Then we can calculate the mean coveragre per each contig with [bamtocov](https://github.com/telatin/bamtocov). We will use one script (`average-coverage.py`) from the program that has can be found from `src` folder. It also needs the location of the actual `bamtocov` that it uses for calulating the coverage.  
```bash
python3 src/average-coverage.py 05_MAPPING/WOD100.bam --bin /projappl/project_2005590/bamtocov/bin/bamtocov
```

Inspect the output from coverage calculation. What is the average coverage over the whole genome? Note that if the genome is in multiple contigs, you will get mean coverage for each of them.  And if some of the contigs are plasmids, their coverage might be different from the rest of the genome. __But why?__

### Long-read mapping with minimap2 (OPTIONAL)

An example how to map the long-reads to the assembly using [minimap2](https://github.com/lh3/minimap2).  

```bash
module load minimap2
module load samtools
minimap2 -ax map-ont path-to/your-assembly.fasta path-to/trimmed-nanopore.fastq,gz |samtools view -Sb |samtools sort > 05_MAPPING/{STRAIN}_nanopore.bam
```

## Genome annotation with Bakta

Now we can annotate our genome assembly using [Bakta](https://github.com/oschwengers/bakta). 

```bash
module purge
/projappl/project_2005590/bakta/bin/bakta \
    --db /scratch/project_2005590/DB/bakta/db/ \
    --prefix {STRAIN} \
    --output 04_ANNOTATION \
    --keep-contig-headers \
    --threads $SLURM_CPUS_PER_TASK \
    path-to/your-assembly
```

Check the files inside the output folder.

## Visualise mapping and annotations with IGV

Download the assembly fasta file, the maping `.bam` (and `.bam.bai`) files and the genome annotation output (`.gff3`) to your own computer.  
We can inspect these with IGV. The next steps will be done together. 


## Annotation and visualization of CRISPR-Cas and Phages (OPTIONAL)

Some genome features are better annotated when considering the genome context of a region involving many genes, instead of looking at only one gene at the time. Two examples of this case are the CRISPR-Cas system and Phages.  

Here we have to examples that are websites where you can upload your assembly.  
The CRISPR-Cas can be annotated using [CRISPRone](https://omics.informatics.indiana.edu/CRISPRone/denovo.php) and Phages can be annotated using [PHASTER](https://phaster.ca/).

Can you find any differences in the annotation of some specific genes when comparing the results of these tools with the Prokka annotation?

## Taxonomic annotation against GTDB

Next thing is to give a name to our strain. In other words, what species did we assemble.  
We will use a tool called [GTDB-tk](https://github.com/ecogenomics/gtdbtk) to give some taxonomy to our genome.

GTDB-tk uses the GTDB database, which is big, but it has been downloaded and can be found from the `DB` folder under the course project.  
We justr need to make an environemntal variable pointing to that, so the tool can find it.  

We will need more memory than previously, so allocate XX G, 6-12 CPUS and 2 hours.

```bash
export GTDBTK_DATA_PATH=/scratch/project_2005590/DB/GTDB/release214/

 /projappl/project_2005590/tax_tools/bin/gtdbtk gtdbtk classify_wf \
    -x fasta \
    --genome_dir 03_ASSEMBLIES \
    --out_dir 04_ANNOTATION/GTDB \
    --cpus $SLURM_CPUS_PER_TASK \
    --tmpdir $LOCAL_SCRATCH
```

The output will contain a file with the most likely taxonomic annotation for your genome. We will go thru the output together.  

## Pangenomics with Anvi'o

### Reference genomes

Before we can do any pangenomics, we need some reference genomes. Go to [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/) and select `Genome`. 
Then search with the taxonomy you got. Either on species level or on genus level. When you get the list of available reference genomes, filter them to include only the ones that have been annotated by NCBI RefSeq and have assembly level at least chromosome, complete would be better. Then select 2-4 genomes for the pangenomic analysis and write their assembly accessions to a file.  

### Setup

We will run the whole anvi'o part interactively. So again, allocate a computing node with enough memory (>40G) and 6 threads for 4 hours.  

Make sure you have at least a copy of each of your genomes in one folder and make a new environmental variable pointing there.

```bash
GENOME_DIR=ABSOLUTE/PATH/TO/GENOME/DIR
```

Then we will polish the contig names for each of the genomes before doing anything with the genomes.  
And also copy few annotated reference genomes to be included in our pangenome.  

```bash
singularity exec --bind $PWD:$PWD,$GENOME_DIR:/genomes $CONTAINERS/anvio_7.sif \
          anvi-script-reformat-fasta --simplify-names -o Oscillatoriales_193.fasta \
          -r reformat_193_report.txt /genomes/GENOME1.fasta

singularity exec --bind $PWD:$PWD,$GENOME_DIR:/genomes $CONTAINERS/anvio_7.sif \
          anvi-script-reformat-fasta --simplify-names -o Oscillatoriales_327_2.fasta \
          -r reformat_327_2_report.txt /genomes/GENOME2.fasta

singularity exec --bind $PWD:$PWD,$GENOME_DIR:/genomes $CONTAINERS/anvio_7.sif \
          anvi-script-reformat-fasta --simplify-names -o Oscillatoriales_328.fasta \
          -r reformat_328_report.txt /genomes/GENOME3.fasta

# copy reference genomes to your own folder
mkdir reference_genomes
cp /scratch/project_2005590/COURSE_FILES/closest_oscillatoriales_genomes/*.gbf ./reference_genomes
```

Although you already annotated your genomes, we'll do it once more, because we changed the names of the contigs.

```bash
module load biokit
for strain in $(ls *.fasta); do prokka --cpus 8 --outdir ./${strain%.fasta}_PROKKA --prefix ${strain%.fasta} $strain; done
```

Now we have both Genbank files from the reference genomes and also from our own genomes.  
Next things is to get the contigs, gene calls and annotations to separate files that anvi'o understands.

```bash
for genome in $(ls */*.gbf)
do
    singularity exec --bind $PWD:$PWD $CONTAINERS/anvio_7.sif \
                                        anvi-script-process-genbank \
                                            -i $genome -O ${genome%.gbf} \
                                            --annotation-source prodigal \
                                            --annotation-version v2.6.3
done
```

The pangenomcis part will be done using a anvi'o workflow (read more from here: )  
And for that we need a file that specifies where the files from previous step are (called `fasta.txt`).  

```bash
echo -e "name\tpath\texternal_gene_calls\tgene_functional_annotation" > fasta.txt
for strain in $(ls */*-contigs.fa)
do
    strain_name=${strain#*/}
    echo -e ${strain_name%-contigs.fa}"\t"$strain"\t"${strain%-contigs.fa}"-external-gene-calls.txt\t"${strain%-contigs.fa}"-external-functions.txt"
done >> fasta.txt
```
In addition to the `fasta.txt` file we need also a configuration file.  
So make a file called `config.json` containing the following.  

```bash
{
    "workflow_name": "pangenomics",
    "config_version": "2",
    "max_threads": "8",
    "project_name": "Oscillatoriales_pangenome",
    "external_genomes": "external-genomes.txt",
    "fasta_txt": "fasta.txt",
    "anvi_gen_contigs_database": {
        "--project-name": "{group}",
        "--description": "",
        "--skip-gene-calling": "",
        "--ignore-internal-stop-codons": true,
        "--skip-mindful-splitting": "",
        "--contigs-fasta": "",
        "--split-length": "",
        "--kmer-size": "",
        "--skip-predict-frame": "",
        "--prodigal-translation-table": "",
        "threads": ""
    },
    "anvi_pan_genome": {
      "threads": "8"
    }
}
```
And then we're ready to run the whole pangenomics workflow.  

```
singularity exec --bind $PWD:$PWD $CONTAINERS/anvio_7.sif anvi-run-workflow -w pangenomics -c config.json
```

When the workflow is ready, we can visualise the results interactively in anvi'o.  
For that we need to connect to Puhti bit differently.  
Log out from the current computing node, open a screen session and a new interactive connection to computing node.
Then take note of the computing node name, the login node number and your port number, you'll need them all in the next part.

### Tunneling the interactive interafce
Although you can install anvi'o on your own computer (and you're free to do so, but we won't have time to help in that), we will run anvi'o in Puhti and tunnel the interactive interface to your local computer.
To be able to to do this, everyone needs to use a different port for tunneling and your port number will be given on the course.

Connecting using a tunnel is a bit tricky and involves several steps, so pay special attention.
Detach from your screen and note on which login node you're on. Then re-attach and note the ID of the computing node your logged in. Then you will also need to remember your port number.

Mini manual for screen:

    screen -S NAME - open a screen and give it a session name NAME
    screen - open new screen without specifying any name
    screen -ls - list all open sessions
    ctrl + a + d - to detach from a session (from inside the screen)
    screen -r NAME - re-attach to a detached session using the name
    screen -rD - re-attach to a attached session
    exit - close the screen and kill all processes running inside the screen (from inside the screen)


Then you can log out and log in again, but this time in a bit different way.
You need to specify your PORT and the NODEID to which you connected and also the NUMBER of the login node you where your screen is running. Also change your username in the command below.
```
ssh -L PORT:NODEID.bullx:PORT USERNAME@puhti-loginNUMBER.csc.fi
```
And in windows using Putty:
In SSH tab select "tunnels". Add:

Source port: PORT
Destination: NODEID.bullx:PORT
Click add and connect to the right login node, login1 or login2.

Then go back to your screen and launch the interactive interface.
Remember to change the PORT.

```
cd 03_PAN
export ANVIOPORT=PORT
singularity exec --bind $PWD:$PWD $CONTAINERS/anvio_7.sif \
                                    anvi-display-pan \
                                        -g Oscillatoriales_pangenome-GENOMES.db \
                                        -p Oscillatoriales_pangenome-PAN.db \
                                        --server-only -P $ANVIOPORT

singularity exec --bind $PWD:$PWD $CONTAINERS/anvio_7.sif \
                                    anvi-get-sequences-for-gene-clusters \
                                        -p Oscillatoriales_pangenome-PAN.db \
                                        -g Oscillatoriales_pangenome-GENOMES.db \
                                        -C default -b SCG \
                                        --concatenate-gene-clusters \
                                        -o single-copy-core-genes.fa                               

singularity exec --bind $PWD:$PWD $CONTAINERS/anvio_7.sif \
                                    anvi-gen-phylogenomic-tree \
                                        -f single-copy-core-genes.fa  \
                                        -o SCG.tre

# study the geosmin phylogeny
singularity exec --bind $PWD:$PWD $CONTAINERS/anvio_7.sif \
                                    anvi-get-sequences-for-gene-clusters \
                                    -p Oscillatoriales_pangenome-PAN.db \
                                    -g Oscillatoriales_pangenome-GENOMES.db \
                                    --report-DNA-sequences \
                                    -C default -b Geosmin -o geosmin.fasta

# This time you need  to align the sequences (use  mafft) and construct the tree from  the aligned sequence file (raxml)
# concatenate the header before alignment
sed -i 's/[|:]/_/g' geosmin.fasta

module load biokit
ginsi geosmin.fasta > geosmin_aln.fasta

module purge
module load raxml
raxmlHPC -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -s geosmin_aln.fasta -n geosmin
```

## Detection of secondary metabolites biosynthesis gene clusters

Biosynthetic genes putatively involved in the synthesis of secondary metabolites can identified using `antiSMASH`

Got to `https://antismash.secondarymetabolites.org/#!/start`. You can load the assembled genome you obtained and turn on all the extra features.

When the analysis is ready, you may be able to answer the following questions:

1. How many secondary metabolites biosynthetic gene clusters (BGC) were detected?
2. Which different types of BGC were detected and what is the difference among these types?
3. Do you think your strain produce all these metabolites? Why?


## Comparison of secondary metabolites biosynthesis gene clusters

Biosynthetic genes clusters can be compared using `BiG-SCAPE` (Biosynthetic Gene Similarity Clustering and Prospecting Engine) and `CORASON` (CORe Analysis of Syntenic Orthologs to prioritize Natural Product-Biosynthetic Gene Cluster)  https://bigscape-corason.secondarymetabolites.org/tutorial/index.html


You can run BiG-SCAPE/CORASON in your folder, using the results obtained from antiSMASH. The .gbk files from the three studied strains and selected reference strains from NCBI are available here: `/scratch/project_2005590/COURSE_FILES/BIGSCAPE/combinedGBK`

You need to load bigscape using a Conda environment:

```bash
export PROJAPPL=/projappl/project_2005590
module load bioconda/3
source activate bigscape
```


```bash
sinteractive -A project_2005590
```

Now you can run BiG-SCAPE/CORASON in your user folder:

```bash
mkdir bigscape
cd bigscape/

python /scratch/project_2005590/COURSE_FILES/BIGSCAPE/BiG-SCAPE/bigscape.py -c 4 -i /scratch/project_2005590/COURSE_FILES/BIGSCAPE/combinedGBK --pfam_dir /scratch/project_2005590/COURSE_FILES/BIGSCAPE/Pfam_database -o bigscape_auto --anchorfile /scratch/project_2005590/COURSE_FILES/BIGSCAPE/anchor_domains.txt --mode auto --hybrids-off --mibig --cutoffs 0.40 0.50 0.60
```

Take a look what it means each parameter used: https://git.wageningenur.nl/medema-group/BiG-SCAPE/-/wikis/home

Once the run is finished, you may transfer the folder to your own computer and observe the results.

Can you find the geosmin and/or 2-methylisoborneol biosynthetic gene clusters?




## Sandbox
Place to store some scratch code while testing.

### checkM
checkM should work from singularity container. Need to pull the right container (tag: 1.1.3--py_0) to course folder and test it once again
```
# needs computing node, otherwise runs out of memory (40G)
singularity exec --bind checkM_test/:/checkM_test /projappl/project_2005590/containers/checkM_1.1.3.sif \
              checkm lineage_wf -x fasta /checkM_test /checkM_test -t 4 --tmpdir /checkM_test
```
### GTDB-tk
Download database before running, needs >200G
```
# download gtdb database
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar -xzf gtdbtk_data.tar.gz

# run gtdbtk
export GTDBTK_DATA_PATH=/scratch/project_2005590 /databases/GTDB/release202/
singularity exec --bind $GTDBTK_DATA_PATH:$GTDBTK_DATA_PATH,$PWD:$PWD  /projappl/project_2005590/containers/gtdbtk_1.7.0.sif \
              gtdbtk classify_wf -x fasta --genome_dir checkM_test/ --out_dir gtdb_test --cpus 4  --tmpdir gtdb_test
```

### basecalling
```

~/projappl/ont-guppy/bin/guppy_basecaller -i fast5_pass/ -s BASECALLED/ -c ~/projappl/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg --device auto --min_qscore 10
cat BASECALLED/pass/*.fastq |gzip > BASECALLED/strain_328_nanopore.fastq.gz
```
