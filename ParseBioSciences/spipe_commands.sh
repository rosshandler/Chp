# spipe is already installed as conda env
conda activate spipe

#split-pipe --version
#split-pipe v0.9.6p

PBS='/data1/ivanir/Chp2022/ParseBS'

# Genome indexing
split-pipe \
--mode mkref \
--genome_name hg38 \
--fasta $PBS/newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
--genes $PBS/newvolume/genomes/Homo_sapiens.GRCh38.93.gtf.gz \
--output_dir $PBS/newvolume/genomes/hg38

# Pipeline running 
split-pipe --mode all --kit WT_mini --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/SLX-21921.DNAA007.H52HVDRX2.s_2.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/SLX-21921.DNAA007.H52HVDRX2.s_2.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/DNAA008 

split-pipe --mode all --kit WT_mini --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/SLX-21921.DNAA008.H52HVDRX2.s_2.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/SLX-21921.DNAA008.H52HVDRX2.s_2.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/DNAA008 

# Combine outputs
split-pipe \
    --mode comb \
    --sublibraries $PBS/newvolume/analysis/DNAA007 $PBS/newvolume/analysis/DNAA008 \
    --output_dir $PBS/newvolume/analysis/Combined
