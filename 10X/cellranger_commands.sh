cellranger count \
--id=ChPH1D19_A \
--transcriptome=/data2/genome_builds/refdata-cellranger-GRCh38-1.2.0/ \
--fastqs=/data1/ivanir/Chp2022/10xGen/fastqs/ --localcores=40 \
--sample=SITTA1 --localmem=48

cellranger count \
--id=ChPH1D19_B \
--transcriptome=/data2/genome_builds/refdata-cellranger-GRCh38-1.2.0/ \
--fastqs=/data1/ivanir/Chp2022/10xGen/fastqs/ --localcores=40 \
--sample=SITTB1 --localmem=48
