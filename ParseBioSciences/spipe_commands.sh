# spipe is already installed as conda env
conda activate spipe

# Genome indexing
PBS='/data1/ivanir/Chp2022/ParseBS'

split-pipe \
--mode mkref \
--genome_name hg38 \
--fasta $PBS/newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
--genes $PBS/newvolume/genomes/Homo_sapiens.GRCh38.93.gtf.gz \
--output_dir $PBS/newvolume/genomes/hg38

