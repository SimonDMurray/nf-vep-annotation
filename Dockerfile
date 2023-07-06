FROM docker.io/ensemblorg/ensembl-vep:latest

USER root

ENV DEBIAN_FRONTEND=noninteractive

RUN cd /opt/vep/.vep && \
    curl -L https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/109/AncestralAllele.pm > AncestralAllele.pm && \
    curl https://ftp.ensembl.org/pub/current_fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz && \
    curl https://ftp.ensembl.org/pub/current_fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz --output homo_sapiens_ancestor_GRCh38.tar.gz && \
    tar xfz homo_sapiens_ancestor_GRCh38.tar.gz && \
    cat homo_sapiens_ancestor_GRCh38/*.fa | bgzip -c > homo_sapiens_ancestor_GRCh38.fa.gz

RUN /opt/vep/src/ensembl-vep/INSTALL.pl -a cf -s homo_sapiens -y GRCh38

RUN cd /opt/vep/.vep && \
    rm -rf homo_sapiens_ancestor_GRCh38/ homo_sapiens_ancestor_GRCh38.tar.gz

RUN cd /opt/vep/.vep && \
    curl -L https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/109/CCDSFilter.pm > CCDSFilter.pm

RUN cd /opt/vep/.vep && \
    curl -L https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/109/CSN.pm > CSN.pm

RUN cd /opt/vep/.vep && \
    curl -L https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/109/Carol.pm > Carol.pm

RUN cd /opt/vep/.vep && \
    curl -L https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/109/Downstream.pm > Downstream.pm

RUN cd /opt/vep/.vep && \
    curl -L https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/109/Draw.pm > Draw.pm

RUN cd /opt/vep/.vep && \
    curl -L https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/109/GXA.pm > GXA.pm

RUN cd /opt/vep/.vep && \
    curl -L https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/109/LOVD.pm > LOVD.pm

RUN cd /opt/vep/.vep && \
    curl -L https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/109/NMD.pm > NMD.pm

RUN cd /opt/vep/.vep && \
    curl -L https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/109/NearestGene.pm > NearestGene.pm
