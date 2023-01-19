## Generate unspliced/spliced matrices (published)
# conda activate velocyto
velocyto run10x \
  -m /data1/ivanir/resources/hg38_repeats_repeatMasker_allTracks.gtf -@ 8 -v \
  /data2/mlancast/cruk/newfastq/ChoroidPlexusA /data2/genome_builds/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf
  
velocyto run10x \
  -m /data1/ivanir/resources/hg38_repeats_repeatMasker_allTracks.gtf -@ 8 -v \
  /data2/mlancast/cruk/newfastq/ChoroidPlexusB /data2/genome_builds/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf 
  
velocyto run10x \
  -m /data1/ivanir/resources/hg38_repeats_repeatMasker_allTracks.gtf -@ 8 -v \
  /data2/mlancast/cruk/newfastq/ChoroidPlexusC /data2/genome_builds/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf 
  
velocyto run10x \
  -m /data1/ivanir/resources/hg38_repeats_repeatMasker_allTracks.gtf -@ 8 -v \
  /data2/mlancast/cruk/newfastq/ChoroidPlexusD /data2/genome_builds/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf 
  
velocyto run10x \
  -m /data1/ivanir/resources/hg38_repeats_repeatMasker_allTracks.gtf -@ 8 -v \
  /data2/mlancast/cruk/newfastq/ChoroidPlexusE /data2/genome_builds/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf 

## Generate unspliced/spliced matrices (2021 data)
# conda activate velocyto
velocyto run10x \
  -m /data1/ivanir/resources/hg38_repeats_repeatMasker_allTracks.gtf -@ 8 -v \
  /data2/mlancast/cruk/SLX-20683/fastqs/ChoroidPlexusH1D19 /data2/genome_builds/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf
  
velocyto run10x \
  -m /data1/ivanir/resources/hg38_repeats_repeatMasker_allTracks.gtf -@ 8 -v \
  /data2/mlancast/cruk/SLX-20683/fastqs/ChoroidPlexusH1D33 /data2/genome_builds/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf 
  
velocyto run10x \
  -m /data1/ivanir/resources/hg38_repeats_repeatMasker_allTracks.gtf -@ 8 -v \
  /data2/mlancast/cruk/SLX-20683/fastqs/ChoroidPlexusH1D80 /data2/genome_builds/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf 
  
velocyto run10x \
  -m /data1/ivanir/resources/hg38_repeats_repeatMasker_allTracks.gtf -@ 8 -v \
  /data2/mlancast/cruk/SLX-20683/fastqs/ChoroidPlexusH9D63 /data2/genome_builds/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf 
  
## Generate unspliced/spliced matrices (2022 data)
# conda activate velocyto
velocyto run10x \
  -m /data1/ivanir/resources/hg38_repeats_repeatMasker_allTracks.gtf -@ 8 -v \
  /data1/ivanir/Chp2022/10xGen/fastqs/ChPH1D19_A/ /data2/genome_builds/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf
  
velocyto run10x \
  -m /data1/ivanir/resources/hg38_repeats_repeatMasker_allTracks.gtf -@ 8 -v \
  /data1/ivanir/Chp2022/10xGen/fastqs/ChPH1D19_B/ /data2/genome_builds/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf 

