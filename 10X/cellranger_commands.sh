cellranger count \
--id=ChPH1D19 \
--transcriptome=/data2/genome_builds/refdata-cellranger-GRCh38-1.2.0/ \
--fastqs=/data1/ivanir/Chp2022/10xGen/fastqs/ --localcores=40 \
--sample=SLX-21922.SITTA1,SLX-21922.SITTB1 --localmem=48

cellranger count \
--id=ChPH1D19 \
--transcriptome=/data2/genome_builds/refdata-cellranger-GRCh38-1.2.0/ \
--fastqs=/data1/ivanir/Chp2022/10xGen/fastqs/ --localcores=40 \
--sample=SLX-21922.SITTA1 --localmem=48
