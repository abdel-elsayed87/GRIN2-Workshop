
#################################################
# GRIN2.0-ALEX library
#######################################################

# get.chrom.length: Retrieve chromosome size data from chr.info txt files available on
# the UCSC genome browser based on the specified genome assembly.

get.chrom.length=function(genome.assembly)   # function support four genome assemblies that include "Human_GRCh38", "Human_GRCh37", "Mouse_HGCm39" and "Mouse_HGCm38"
{
  # retrieve chromosome size data for GRCh38 (hg38) genome build
  if (genome.assembly=="Human_GRCh38")
  {
    chr.size.hg38= circlize::read.chromInfo(species = "hg38")
    chr.size.hg38=as.data.frame(chr.size.hg38)
    chr.size=cbind.data.frame(chrom=chr.size.hg38$chromosome,
                              size=chr.size.hg38$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))
    
    return(chr.size)
  }
  
  # retrieve chromosome size data for GRCh37 (hg19) genome build
  if (genome.assembly=="Human_GRCh37")
  {
    chr.size.hg19= circlize::read.chromInfo(species = "hg19")
    chr.size.hg19=as.data.frame(chr.size.hg19)
    chr.size=cbind.data.frame(chrom=chr.size.hg19$chromosome,
                              size=chr.size.hg19$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))
    
    return(chr.size)
  }
  
  # retrieve chromosome size data for Mouse_HGCm39 (mm39) genome build
  if (genome.assembly=="Mouse_HGCm39")
  {
    chr.size.mm39= circlize::read.chromInfo(species = "mm39")
    chr.size.mm39=as.data.frame(chr.size.mm39)
    chr.size=cbind.data.frame(chrom=chr.size.mm39$chromosome,
                              size=chr.size.mm39$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))
    
    return(chr.size)
  }
  
  # retrieve chromosome size data for Mouse_HGCm38 (mm10) genome build
  if (genome.assembly=="Mouse_HGCm38")
  {
    chr.size.mm38= circlize::read.chromInfo(species = "mm10")
    chr.size.mm38=as.data.frame(chr.size.mm38)
    chr.size=cbind.data.frame(chrom=chr.size.mm38$chromosome,
                              size=chr.size.mm38$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))
    
    return(chr.size)
  }
}

################################################################################
# get.ensembl.annotation: Retrieve gene and regulatory features annotation data from ensembl biomaRt 
# database based on the specified genome assembly.

get.ensembl.annotation=function(genome.assembly)  # one of four genome assemblies that include "Human_GRCh38", "Human_GRCh37", "Mouse_HGCm39" and "Mouse_HGCm38" can be specified
  
{
  # retrieve gene annotation data for human_GRCh38 (hg38) genome assembly from ensembl biomaRt version 104
  if (genome.assembly=="Human_GRCh38")
  {
    message(paste0("Retrieving gene annotation data from Ensembl BioMart: ",date()))
    ensembl_GRCh38 = useEnsembl(biomart="genes",
                                dataset="hsapiens_gene_ensembl",  # specify dataset for homosapiens
                                version = "110")                  # specifying version is critical to get stable search results. If we did not specify version, query will extract data from the most updated version
    
    chromosomes = c(1:22, "X", "Y")                               # specify chromosomes of interest
    hg38_gene_annotation <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                               'start_position','end_position',
                                               'description','external_gene_name',
                                               'gene_biotype', 'strand', 'band'),
                                  filters = 'chromosome_name',  values = chromosomes,
                                  mart = ensembl_GRCh38)         # attributes specify data to be retrieved from biomaRt database
    
    gene=hg38_gene_annotation[,1]
    chrom=hg38_gene_annotation[,2]
    loc.start=hg38_gene_annotation[,3]
    loc.end=hg38_gene_annotation[,4]
    description=hg38_gene_annotation[,5]
    gene.name=hg38_gene_annotation[,6]
    biotype=hg38_gene_annotation[,7]
    chrom.strand=hg38_gene_annotation[,8]
    chrom.band=hg38_gene_annotation[,9]
    
    gene.data=cbind.data.frame(gene=gene,  
                               chrom=chrom,
                               loc.start=loc.start,
                               loc.end=loc.end,
                               description=description,
                               gene.name=gene.name,
                               biotype=biotype,
                               chrom.strand=chrom.strand,
                               chrom.band=chrom.band)
    
    # To retrieve data for computationally predicted regulatory features mapped to GRCh38 (regulatory build)
    message(paste0("Retrieving computationally predicted regulatory features from Ensembl regulatory build: ",date()))
    
    regulatory.hg38 = useEnsembl(biomart="regulation",
                                 dataset="hsapiens_regulatory_feature", # regulatory features includes (promoters, promoter flanking regions, enhancers, CTCF binding sites, TF binding sites and open chromatin regions)
                                 version = '110')                 # specify version for stable results
    
    regulatory.features=c("Promoter", "Promoter Flanking Region", "Enhancer",
                          "CTCF Binding Site", "TF binding site", "Open chromatin") # specify classes of regulatory features
    
    hg38.regulatory= getBM(attributes=c("regulatory_stable_id", "chromosome_name",
                                        "chromosome_start","chromosome_end",
                                        "feature_type_description"),
                           filters =c("regulatory_feature_type_name", "chromosome_name")
                           , values =list(regulatory.features,chromosomes) ,      # We are excluding regulatory features not mapped to chr1:22, X, Y
                           mart =regulatory.hg38)
    
    gene=hg38.regulatory[,1]
    chrom=hg38.regulatory[,2]
    loc.start=hg38.regulatory[,3]
    loc.end=hg38.regulatory[,4]
    description=hg38.regulatory[,5]
    
    regulation.data=cbind.data.frame(gene=gene,
                                     chrom=chrom,
                                     loc.start=loc.start,
                                     loc.end=loc.end,
                                     description=description)
    
    namevector <- c("gene.name", "biotype", "chrom.strand", "chrom.band")
    regulation.data[ , namevector] <- NA
    regulation.data.final=rbind.data.frame(regulation.data, gene.data)
    
    # To retrieve data for experimentally validated regulatory features mapped to GRCh38 (FANTOM project)
    message(paste0("Retrieving experimentally validated regulatory features (FANTOM project): ",date()))
    
    reg.hg38.val = useEnsembl(biomart="regulation",
                              dataset="hsapiens_external_feature", # regulatory features includes (promoters, promoter flanking regions, enhancers, CTCF binding sites, TF binding sites and open chromatin regions)
                              version = '110')                 # specify version for stable results
    
    hg38.reg.val= getBM(attributes=c("display_label", "chromosome_name",
                                     "chromosome_start","chromosome_end",
                                     "feature_type_description", "feature_type_class"),
                        filters =c("external_feature_set_name", "chromosome_name")
                        , values =list("FANTOM predictions",chromosomes) ,      # We are excluding regulatory features not mapped to chr1:22, X, Y
                        mart =reg.hg38.val)
    
    gene=hg38.reg.val[,1]
    chrom=hg38.reg.val[,2]
    loc.start=hg38.reg.val[,3]
    loc.end=hg38.reg.val[,4]
    description=hg38.reg.val[,5]
    
    regulation.data.val=cbind.data.frame(gene=gene,
                                         chrom=chrom,
                                         loc.start=loc.start,
                                         loc.end=loc.end,
                                         description=description)
    
    namevector <- c("gene.name", "biotype", "chrom.strand", "chrom.band")
    regulation.data.val[ , namevector] <- NA
    
    regulation.data.val.final=rbind.data.frame(regulation.data.val, gene.data)
    
    res=list(gene.annotation=gene.data,                    # gene.data represent annotated genes
             reg.annotation.predicted=regulation.data.final,     # predicted regulatory features (ensembl regulatory build)
             reg.annotation.validated=regulation.data.val.final) # experimentally validated regulatory features (FANTOM project)  
    
    return(res)
  }
  
  # retrieve gene annotation data for human_GRCh37 (hg19) genome assembly from ensembl biomaRt version 75
  if (genome.assembly=="Human_GRCh37")
  {
    message(paste0("Retrieving gene annotation data from Ensembl BioMart: ",date()))
    ensemblGRCh37 <- useEnsembl(biomart = 'ensembl',
                                dataset = 'hsapiens_gene_ensembl',
                                version = '75')
    
    chromosomes = c(1:22, "X", "Y")
    hg19_gene_annotation <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                               'start_position','end_position',
                                               'description','external_gene_id',
                                               'gene_biotype', 'strand', 'band'),
                                  filters = 'chromosome_name',  values = chromosomes,
                                  mart = ensemblGRCh37)
    
    gene=hg19_gene_annotation[,1]
    chrom=hg19_gene_annotation[,2]
    loc.start=hg19_gene_annotation[,3]
    loc.end=hg19_gene_annotation[,4]
    description=hg19_gene_annotation[,5]
    gene.name=hg19_gene_annotation[,6]
    biotype=hg19_gene_annotation[,7]
    chrom.strand=hg19_gene_annotation[,8]
    chrom.band=hg19_gene_annotation[,9]
    
    gene.data=cbind.data.frame(gene=gene,
                               chrom=chrom,
                               loc.start=loc.start,
                               loc.end=loc.end,
                               description=description,
                               gene.name=gene.name,
                               biotype=biotype,
                               chrom.strand=chrom.strand,
                               chrom.band=chrom.band)
    
    # To retrieve data for regulatory features mapped to GRCh37
    message(paste0("Retrieving computationally predicted regulatory features from Ensembl regulatory build: ",date()))
    
    regulatory.hg19 = useEnsembl(biomart="regulation",
                                 dataset="hsapiens_regulatory_feature",
                                 version = 'GRCh37')
    
    regulatory.features=c("Promoter", "Promoter Flanking Region", "Enhancer",
                          "CTCF Binding Site",
                          "TF binding site", "Open chromatin")
    
    hg19.regulatory= getBM(attributes=c("regulatory_stable_id", "chromosome_name",
                                        "chromosome_start","chromosome_end",
                                        "feature_type_description"),
                           filters =c("regulatory_feature_type_name", "chromosome_name")
                           , values =list(regulatory.features,chromosomes) ,      # We are excluding regulatory features not mapped to chr1:22, X, Y
                           mart =regulatory.hg19)
    
    gene=hg19.regulatory[,1]
    chrom=hg19.regulatory[,2]
    loc.start=hg19.regulatory[,3]
    loc.end=hg19.regulatory[,4]
    description=hg19.regulatory[,5]
    
    regulation.data=cbind.data.frame(gene=gene,
                                     chrom=chrom,
                                     loc.start=loc.start,
                                     loc.end=loc.end,
                                     description=description)
    
    namevector <- c("gene.name", "biotype", "chrom.strand", "chrom.band")
    regulation.data[ , namevector] <- NA
    regulation.data.final=rbind.data.frame(regulation.data, gene.data)
    
    # To retrieve data for experimentally validated regulatory features mapped to GRCh37 (FANTOM project)
    message(paste0("Retrieving experimentally validated regulatory features (FANTOM project): ",date()))
    
    reg.hg19.val = useEnsembl(biomart="regulation",
                              dataset="hsapiens_external_feature", # regulatory features includes (promoters, promoter flanking regions, enhancers, CTCF binding sites, TF binding sites and open chromatin regions)
                              version = 'GRCh37')                 # specify version for stable results
    
    hg19.reg.val= getBM(attributes=c("display_label", "chromosome_name",
                                     "chromosome_start","chromosome_end",
                                     "feature_type_description", "feature_type_class"),
                        filters =c("external_feature_set_name", "chromosome_name")
                        , values =list("FANTOM predictions",chromosomes) ,      # We are excluding regulatory features not mapped to chr1:22, X, Y
                        mart =reg.hg19.val)
    
    gene=hg19.reg.val[,1]
    chrom=hg19.reg.val[,2]
    loc.start=hg19.reg.val[,3]
    loc.end=hg19.reg.val[,4]
    description=hg19.reg.val[,5]
    
    regulation.data.val=cbind.data.frame(gene=gene,
                                         chrom=chrom,
                                         loc.start=loc.start,
                                         loc.end=loc.end,
                                         description=description)
    
    namevector <- c("gene.name", "biotype", "chrom.strand", "chrom.band")
    regulation.data.val[ , namevector] <- NA
    regulation.data.val.final=rbind.data.frame(regulation.data.val, gene.data)
    
    res=list(gene.annotation=gene.data,                    # gene.data represent annotated genes
             reg.annotation.predicted=regulation.data.final,     # predicted regulatory features (ensembl regulatory build)
             reg.annotation.validated=regulation.data.val.final) # experimentally validated regulatory features (FANTOM project)
    
    return(res)
  }
  
  # retrieve gene annotation data for Mouse_HGCm39 genome assembly from ensembl biomaRt version 104
  if (genome.assembly=="Mouse_HGCm39")
  {
    message(paste0("Retrieving gene annotation data from Ensembl BioMart: ",date()))
    ensembl.HGCm39 = useEnsembl(biomart="genes", 
                                dataset="mmusculus_gene_ensembl",
                                version = "104")
    
    chromosomes = c(1:19, "X", "Y")
    HGCm39.gene_annotation <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                                 'start_position','end_position',
                                                 'description','external_gene_name',
                                                 'gene_biotype', 'strand', 'band'),
                                    filters = 'chromosome_name',  values = chromosomes,
                                    mart = ensembl.HGCm39)
    
    gene=HGCm39.gene_annotation[,1]
    chrom=HGCm39.gene_annotation[,2]
    loc.start=HGCm39.gene_annotation[,3]
    loc.end=HGCm39.gene_annotation[,4]
    description=HGCm39.gene_annotation[,5]
    gene.name=HGCm39.gene_annotation[,6]
    biotype=HGCm39.gene_annotation[,7]
    chrom.strand=HGCm39.gene_annotation[,8]
    chrom.band=HGCm39.gene_annotation[,9]
    
    gene.data=cbind.data.frame(gene=gene,
                               chrom=chrom,
                               loc.start=loc.start,
                               loc.end=loc.end,
                               description=description,
                               gene.name=gene.name,
                               biotype=biotype,
                               chrom.strand=chrom.strand,
                               chrom.band=chrom.band)
    
    # To retrieve data for mouse_regulatory features mapped to mm39
    message(paste0("Retrieving computationally predicted regulatory features from Ensembl regulatory build: ",date()))
    
    regulatory.mm39 = useEnsembl(biomart="regulation",
                                 dataset="mmusculus_regulatory_feature",
                                 version = '104')
    
    regulatory.features=c("Promoter", "Promoter Flanking Region", "Enhancer",
                          "CTCF Binding Site",
                          "TF binding site", "Open chromatin")
    mm39.regulatory= getBM(attributes=c("regulatory_stable_id", "chromosome_name",
                                        "chromosome_start","chromosome_end",
                                        "feature_type_description"),
                           filters =c("regulatory_feature_type_name", "chromosome_name")
                           , values =list(regulatory.features,chromosomes) ,      # We are excluding regulatory features not mapped to chr1:22, X, Y
                           mart =regulatory.mm39)
    
    gene=mm39.regulatory[,1]
    chrom=mm39.regulatory[,2]
    loc.start=mm39.regulatory[,3]
    loc.end=mm39.regulatory[,4]
    description=mm39.regulatory[,5]
    
    regulation.data=cbind.data.frame(gene=gene,
                                     chrom=chrom,
                                     loc.start=loc.start,
                                     loc.end=loc.end,
                                     description=description)
    
    namevector <- c("gene.name", "biotype", "chrom.strand", "chrom.band")
    regulation.data[ , namevector] <- NA
    regulation.data.final=rbind.data.frame(regulation.data, gene.data)
    
    res=list(gene.annotation=gene.data,
             reg.annotation=regulation.data.final)
    return(res)
  }
  
  # retrieve gene annotation data for Mouse_HGCm38 from ensembl biomaRt version 102
  if (genome.assembly=="Mouse_HGCm38")
  {
    message(paste0("Retrieving gene annotation data from Ensembl BioMart: ",date()))
    ensemblGRCm38 <- useEnsembl(biomart = 'genes',
                                dataset = 'mmusculus_gene_ensembl',
                                version = '102')
    
    chromosomes = c(1:19, "X", "Y")
    Mouse_mm10.gene_annotation <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                                     'start_position','end_position',
                                                     'description','external_gene_name',
                                                     'gene_biotype', 'strand', 'band'),
                                        filters = 'chromosome_name',  values = chromosomes,
                                        mart = ensemblGRCm38)
    
    gene=Mouse_mm10.gene_annotation[,1]
    chrom=Mouse_mm10.gene_annotation[,2]
    loc.start=Mouse_mm10.gene_annotation[,3]
    loc.end=Mouse_mm10.gene_annotation[,4]
    description=Mouse_mm10.gene_annotation[,5]
    gene.name=Mouse_mm10.gene_annotation[,6]
    biotype=Mouse_mm10.gene_annotation[,7]
    chrom.strand=Mouse_mm10.gene_annotation[,8]
    chrom.band=Mouse_mm10.gene_annotation[,9]
    
    gene.data=cbind.data.frame(gene=gene,
                               chrom=chrom,
                               loc.start=loc.start,
                               loc.end=loc.end,
                               description=description,
                               gene.name=gene.name,
                               biotype=biotype,
                               chrom.strand=chrom.strand,
                               chrom.band=chrom.band)
    
    # To retrieve data for mouse_regulatory features mapped to mm10
    message(paste0("Retrieving computationally predicted regulatory features from Ensembl regulatory build: ",date()))
    
    regulatory.mm10 = useEnsembl(biomart="regulation",
                                 dataset="mmusculus_regulatory_feature",
                                 version = '102')
    regulatory.features=c("Promoter", "Promoter Flanking Region", "Enhancer",
                          "CTCF Binding Site",
                          "TF binding site", "Open chromatin")
    mm10.regulatory= getBM(attributes=c("regulatory_stable_id", "chromosome_name",
                                        "chromosome_start","chromosome_end",
                                        "feature_type_description"),
                           filters =c("regulatory_feature_type_name", "chromosome_name")
                           , values =list(regulatory.features,chromosomes) ,      # We are excluding regulatory features not mapped to chr1:22, X, Y
                           mart =regulatory.mm10)
    
    gene=mm10.regulatory[,1]
    chrom=mm10.regulatory[,2]
    loc.start=mm10.regulatory[,3]
    loc.end=mm10.regulatory[,4]
    description=mm10.regulatory[,5]
    
    regulation.data=cbind.data.frame(gene=gene,
                                     chrom=chrom,
                                     loc.start=loc.start,
                                     loc.end=loc.end,
                                     description=description)
    
    namevector <- c("gene.name", "biotype", "chrom.strand", "chrom.band")
    regulation.data[ , namevector] <- NA
    regulation.data.final=rbind.data.frame(regulation.data, gene.data)
    
    res=list(gene.annotation=gene.data,
             reg.annotation=regulation.data.final)
    
    return(res)
  }
}


#################################################################################
# order.index.gene.data: This function order and index gene annotation data by chromosome,
# gene start, and end location.

order.index.gene.data=function(gene.data)  # gene annotation data
{
  g=nrow(gene.data)
  gene.ord=order(gene.data[,"chrom"],
                 gene.data[,"loc.start"],
                 gene.data[,"loc.end"])
  
  gene.data=gene.data[gene.ord,]
  new.chrom=which(gene.data[-1,"chrom"]!=gene.data[-g,"chrom"])
  chr.start=c(1,new.chrom+1)
  chr.end=c(new.chrom,g)
  gene.index=cbind.data.frame(chrom=gene.data[chr.start,"chrom"],
                              row.start=chr.start,
                              row.end=chr.end)
  gene.data$gene.row=1:g
  
  res=list(gene.data=gene.data,
           gene.index=gene.index)
  return(res)
}

#################################################################################
# order.index.lsn.data: This function order and index lesion data by lesion type, the
# chromosome on which the lesion is located, and subject.

order.index.lsn.data=function(lsn.data)  # lesion data provided by the user in a GRIN compatible format 
  
{
  l=nrow(lsn.data)
  lsn.ord=order(lsn.data[,"lsn.type"],
                lsn.data[,"chrom"],
                lsn.data[,"ID"])
  lsn.data=lsn.data[lsn.ord,]
  
  lsn.chng=which((lsn.data[-1,"lsn.type"]!=lsn.data[-l,"lsn.type"])|
                   (lsn.data[-1,"chrom"]!=lsn.data[-l,"chrom"])|
                   (lsn.data[-1,"ID"]!=lsn.data[-l,"ID"]))
  lsn.start=c(1,lsn.chng+1)
  lsn.end=c(lsn.chng,l)
  lsn.index=cbind.data.frame(lsn.type=lsn.data[lsn.start,"lsn.type"],
                             chrom=lsn.data[lsn.start,"chrom"],
                             ID=lsn.data[lsn.start,"ID"],
                             row.start=lsn.start,
                             row.end=lsn.end)
  lsn.data$lsn.row=1:l
  
  res=list(lsn.data=lsn.data,
           lsn.index=lsn.index)
  return(res)
}

##################################################################################
# prep.gene.lsn.data: This function prepare gene and lesion data for later computations.

prep.gene.lsn.data=function(lsn.data,       # lesion data file with five columns: "ID" which is patient ID, "chrom" chromosome on which the lesion is located, "loc.start" lesion start position, "loc.end" lesion end position, "lsn.type" lesion category such as gain, loss, mutation, etc...
                            gene.data,      # gene annotation data with four required columns: "gene" has ensembl ID, "chrom" chromosome on which the gene is located, "loc.start" gene start position, "loc.end" which is the gene end position
                            mess.freq=10)   # message frequency: display message every mess.freq^{th} lesion block
{
  saf=options()$stringsAsFactors
  options(stringsAsFactors=F)
  
  # order lesion data by type, chromosome, and subject
  lsn.dset=order.index.lsn.data(lsn.data)
  lsn.data=lsn.dset$lsn.data
  lsn.index=lsn.dset$lsn.index
  
  # order and index gene locus data by chromosome and position
  gene.dset=order.index.gene.data(gene.data)
  gene.data=gene.dset$gene.data
  gene.index=gene.dset$gene.index
  
  # Extract some basic information
  g=nrow(gene.data) # number of genes
  l=nrow(lsn.data)  # number of lesions
  
  # Create gene position data
  message(paste0("Formatting gene position data for counting: ",date()))
  gene.pos.data=rbind.data.frame(cbind.data.frame(ID="",  # gene start data
                                                  lsn.type="",
                                                  lsn.row=NA,
                                                  gene=gene.data[,"gene"],
                                                  gene.row=gene.data[,"gene.row"],
                                                  chrom=gene.data[,"chrom"],
                                                  pos=gene.data[,"loc.start"],
                                                  cty=1),
                                 cbind.data.frame(ID="", # gene end data
                                                  lsn.type="",
                                                  lsn.row=NA,
                                                  gene=gene.data[,"gene"],
                                                  gene.row=gene.data[,"gene.row"],
                                                  chrom=gene.data[,"chrom"],
                                                  pos=gene.data[,"loc.end"],
                                                  cty=4))
  # order gene position data
  ord=order(gene.pos.data[,"chrom"],
            gene.pos.data[,"pos"],
            gene.pos.data[,"cty"])
  gene.pos.data=gene.pos.data[ord,]
  
  # Create lesion position data with one row for each edge of each lesion
  message(paste0("Formatting lesion position data for counting:  ",date()))
  lsn.pos.data=rbind.data.frame(cbind.data.frame(ID=lsn.data[,"ID"],
                                                 lsn.type=lsn.data[,"lsn.type"],
                                                 lsn.row=lsn.data[,"lsn.row"],
                                                 gene="",
                                                 gene.row=NA,
                                                 chrom=lsn.data[,"chrom"],
                                                 pos=lsn.data[,"loc.start"],
                                                 cty=2),
                                cbind.data.frame(ID=lsn.data[,"ID"],
                                                 lsn.type=lsn.data[,"lsn.type"],
                                                 lsn.row=lsn.data[,"lsn.row"],
                                                 gene="",
                                                 gene.row=NA,
                                                 chrom=lsn.data[,"chrom"],
                                                 pos=lsn.data[,"loc.end"],
                                                 cty=3))
  # order lesion position data
  ord=order(lsn.pos.data[,"chrom"],
            lsn.pos.data[,"pos"],
            lsn.pos.data[,"cty"])
  lsn.pos.data=lsn.pos.data[ord,]
  
  # Combine gene & lesion data
  message(paste0("Combining formatted gene and lesion postion data: ",date()))
  gene.lsn.data=rbind.data.frame(gene.pos.data,
                                 lsn.pos.data)
  
  # Order and index gene & lesion data
  ord=order(gene.lsn.data[,"chrom"],
            gene.lsn.data[,"pos"],
            gene.lsn.data[,"cty"])
  gene.lsn.data=gene.lsn.data[ord,]
  m=nrow(gene.lsn.data)
  gene.lsn.data[,"glp.row"]=1:m
  
  # compute vector to order gene.lsn.data by lsn.row and gene.row
  ord=order(gene.lsn.data[,"lsn.row"],
            gene.lsn.data[,"gene.row"],
            gene.lsn.data[,"cty"])
  
  # use that vector to add gene.lsn.data row.start and row.end indices to lsn.data
  lsn.pos=gene.lsn.data[ord[1:(2*l)],]
  lsn.data[,"glp.row.start"]=lsn.pos[2*(1:l)-1,"glp.row"]
  lsn.data[,"glp.row.end"]=lsn.pos[2*(1:l),"glp.row"]
  
  
  # use that vector to add gene.lsn.data row.start and row.end indices to gene.data
  gene.pos=gene.lsn.data[ord[-(1:(2*l))],]
  gene.data[,"glp.row.start"]=gene.pos[2*(1:g)-1,"glp.row"]
  gene.data[,"glp.row.end"]=gene.pos[2*(1:g),"glp.row"]
  
  # Double-check table pointers from lsn.data and gene.data to gene.lsn.data
  
  message(paste0("Verifying structure of combined gene and lesion data: ",date()))
  glp.gene.start=gene.lsn.data[gene.data$glp.row.start,c("gene","chrom","pos")]
  colnames(glp.gene.start)=c("gene","chrom","loc.start")
  ok.glp.gene.start=all(glp.gene.start==gene.data[,c("gene","chrom","loc.start")])
  
  glp.gene.end=gene.lsn.data[gene.data$glp.row.end,c("gene","chrom","pos")]
  colnames(glp.gene.end)=c("gene","chrom","loc.end")
  ok.glp.gene.end=all(glp.gene.end==gene.data[,c("gene","chrom","loc.end")])
  
  glp.lsn.start=gene.lsn.data[lsn.data$glp.row.start,c("ID","chrom","pos","lsn.type")]
  colnames(glp.lsn.start)=c("ID","chrom","loc.start","lsn.type")
  ok.glp.lsn.start=all(glp.lsn.start==lsn.data[,c("ID","chrom","loc.start","lsn.type")])
  
  glp.lsn.end=gene.lsn.data[lsn.data$glp.row.end,c("ID","chrom","pos","lsn.type")]
  colnames(glp.lsn.end)=c("ID","chrom","loc.end","lsn.type")
  ok.glp.lsn.end=all(glp.lsn.end==lsn.data[,c("ID","chrom","loc.end","lsn.type")])
  
  # Double-check table pointers from gene.lsn.data to gene.data and lsn.data
  glp.gene.start=gene.lsn.data[gene.lsn.data$cty==1,c("gene.row","gene","chrom","pos")]
  ok.gene.start=all(glp.gene.start[,c("gene","chrom","pos")]==gene.data[glp.gene.start$gene.row,c("gene","chrom","loc.start")])
  
  glp.gene.end=gene.lsn.data[gene.lsn.data$cty==4,c("gene.row","gene","chrom","pos")]
  ok.gene.end=all(glp.gene.end[,c("gene","chrom","pos")]==gene.data[glp.gene.end$gene.row,c("gene","chrom","loc.end")])
  
  glp.lsn.start=gene.lsn.data[gene.lsn.data$cty==2,c("lsn.row","ID","chrom","pos","lsn.type")]
  ok.lsn.start=all(glp.lsn.start[,c("ID","chrom","pos","lsn.type")]==lsn.data[glp.lsn.start$lsn.row,c("ID","chrom","loc.start","lsn.type")])
  
  glp.lsn.end=gene.lsn.data[gene.lsn.data$cty==3,c("lsn.row","ID","chrom","pos","lsn.type")]
  ok.lsn.end=all(glp.lsn.end[,c("ID","chrom","pos","lsn.type")]==lsn.data[glp.lsn.end$lsn.row,c("ID","chrom","loc.end","lsn.type")])
  
  all.ok=all(c(ok.glp.gene.start,ok.glp.gene.end,
               ok.glp.lsn.start,ok.glp.lsn.end,
               ok.gene.start,ok.gene.end,
               ok.lsn.start,ok.lsn.end))
  
  if (!all.ok)
    stop("Error in constructing and indexing combined lesion and gene data.")
  
  message(paste0("Verified correct construction and indexing of combined lesion and gene data: ",date()))
  
  return(list(lsn.data=lsn.data, # Input lesion data
              gene.data=gene.data, # Input gene annotation data
              gene.lsn.data=gene.lsn.data, # data.frame ordered by gene and lesions start position. Gene start position is coded as 1 in the cty column and gene end position is coded as 4. Lesion start position is coded as 2 in the cty column and lesion end position is coded as 3.
              gene.index=gene.index, # data.frame that shows ordered row start and row end for each chromosome in the gene.lsn.data table
              lsn.index=lsn.index)) # data.frame that shows row start and row end for each lesion in the gene.lsn.data table
  
}

#################################################################################
# find.gene.lsn.overlaps: Use the results of prep.gene.lsn.data function to find
# lesion-gene overlaps.

find.gene.lsn.overlaps=function(gl.data) # A list of five data.frames that represent the output results of the prep.gene.lsn.data function
  
{
  gene.data=gl.data$gene.data
  lsn.data=gl.data$lsn.data
  lsn.index=gl.data$lsn.index
  gene.index=gl.data$gene.index
  gene.lsn.data=gl.data$gene.lsn.data
  m=nrow(gene.lsn.data)
  
  message(paste0("Scanning through combined lesion and gene data to find gene-lesion overlaps: ",date()))
  
  
  gene.row.mtch=NULL  # initialize vector for rows of gene data matched to rows of lesion data
  lsn.row.mtch=NULL   # initialize vector for rows of lesion data matched to rows of gene data
  current.genes=NULL  # initialize vector of genes overlapping this point of the scan
  current.lsns=NULL   # initialize vector of lesions overlapping this point of the scan
  for (i in 1:m)      # loop over rows of gene.lsn.data
  {
    # enter a gene
    if (gene.lsn.data$cty[i]==1)
    {
      # add this gene to the set of current genes
      current.genes=c(current.genes,
                      gene.lsn.data$gene.row[i])
      
      # match this gene to set of current lesions
      lsn.row.mtch=c(lsn.row.mtch,current.lsns)          # add current lesions to lsn.row.mtch
      gene.row.mtch=c(gene.row.mtch,                     # add this gene for each current lesion
                      rep(gene.lsn.data$gene.row[i],
                          length(current.lsns)))
      
    }
    
    # exit a gene
    if (gene.lsn.data$cty[i]==4)
    {
      # drop this gene from the set of current genes
      current.genes=setdiff(current.genes,
                            gene.lsn.data$gene.row[i])
    }
    
    
    if (gene.lsn.data$cty[i]==2)       # enter a lesion
    {
      lsn.row.mtch=c(lsn.row.mtch,
                     rep(gene.lsn.data$lsn.row[i],length(current.genes)))
      gene.row.mtch=c(gene.row.mtch,current.genes)
      current.lsns=c(current.lsns,
                     gene.lsn.data$lsn.row[i])
    }
    
    if (gene.lsn.data$cty[i]==3)
    {
      current.lsns=setdiff(current.lsns,
                           gene.lsn.data$lsn.row[i])
    }
  }
  
  message(paste0("Completed scan of combined gene-lesion data: ",date()))
  
  # Generate the gene-lesion hit data
  gene.lsn.hits=cbind.data.frame(gene.data[gene.row.mtch,
                                           c("gene.row","gene","chrom","loc.start","loc.end")],
                                 lsn.data[lsn.row.mtch,
                                          c("lsn.row","ID","chrom","loc.start","loc.end","lsn.type")])
  
  colnames(gene.lsn.hits)=c("gene.row","gene","gene.chrom","gene.loc.start","gene.loc.end",
                            "lsn.row","ID","lsn.chrom","lsn.loc.start","lsn.loc.end","lsn.type")
  
  res=list(lsn.data=lsn.data, # Input lesion data
           gene.data=gene.data, # Input gene annotation data
           gene.lsn.data=gene.lsn.data, # Data.frame ordered by gene and lesions start position. Gene start position is coded as 1 in the cty column and gene end position is coded as 4. Lesion start position is coded as 2 in the cty column and lesion end position is coded as 3
           gene.lsn.hits=gene.lsn.hits, # Each row represent a gene overlapped by a certain lesion. Gene column shows the overlapped gene and ID column has the patient ID
           gene.index=gene.index, # Data.frame that shows row start and row end for each chromosome in the gene.lsn.data table
           lsn.index=lsn.index) # Data.frame that shows row start and row end for each lesion in the gene.lsn.data table
  
  return(res)
  
}

##################################################################################
# count.hits: The function computes the number of hits affecting each gene by lesion
# type. It also compute the number of subjects with a hit in each annotated gene.

count.hits=function(ov.data) # output results of find.gene.lsn.overlaps function
{
  lsn.data=ov.data$lsn.data
  lsn.index=ov.data$lsn.index
  gene.lsn.hits=ov.data$gene.lsn.hits
  gene.lsn.data=ov.data$gene.lsn.data
  gene.data=ov.data$gene.data
  gene.index=ov.data$gene.index
  
  g=nrow(gene.data)
  
  # Compute the number of hits matrix
  lsn.types=sort(unique(lsn.index[,"lsn.type"]))
  k=length(lsn.types)
  
  nhit.mtx=matrix(0,g,k)
  colnames(nhit.mtx)=lsn.types
  
  nhit.tbl=table(gene.lsn.hits$gene.row,
                 gene.lsn.hits$lsn.type)
  nhit.rows=as.numeric(rownames(nhit.tbl))
  
  for (i in 1:ncol(nhit.tbl))
    nhit.mtx[nhit.rows,colnames(nhit.tbl)[i]]=nhit.tbl[,i]
  
  # Compute the matrix of the number of subjects with a hit
  gene.subj.type=paste0(gene.lsn.hits$gene.row,"_",
                        gene.lsn.hits$ID,"_",
                        gene.lsn.hits$lsn.type)
  dup.gene.subj.type=duplicated(gene.subj.type)
  
  subj.gene.hits=gene.lsn.hits[!dup.gene.subj.type,]
  
  nsubj.mtx=matrix(0,g,k)
  colnames(nsubj.mtx)=lsn.types
  
  nsubj.tbl=table(subj.gene.hits$gene.row,
                  subj.gene.hits$lsn.type)
  
  nsubj.rows=as.numeric(rownames(nsubj.tbl))
  
  for (i in 1:ncol(nsubj.tbl))
    nsubj.mtx[nsubj.rows,colnames(nsubj.tbl)[i]]=nsubj.tbl[,i]
  
  res=list(lsn.data=lsn.data,  # Input lesion data
           lsn.index=lsn.index, # data.frame that shows row start and row end for each lesion in the gene.lsn.data table
           gene.data=gene.data, # Input gene annotation data
           gene.index=gene.index, # data.frame that shows ordered row start and row end for each chromosome in the gene.lsn.data table
           nhit.mtx=nhit.mtx, # A data matrix with number of hits in each gene by lesion type
           nsubj.mtx=nsubj.mtx, # A data matrix with number of affected subjects by lesion type
           gene.lsn.data=gene.lsn.hits, # Each row represent a gene overlapped by a certain lesion. Column "gene" shows the overlapped gene and ID column has the patient ID
           glp.data=gene.lsn.data) # data.frame ordered by gene and lesions start position. Gene start position is coded as 1 in the cty column and gene end position is coded as 4. Lesion start position is coded as 2 in the cty column and lesion end position is coded as 3
  
  return(res)
  
}

#################################################################################
# Following four functions row.prob.subj.hit, row.bern.conv, p.order and pc06.fdr
# will be used by prob.hits function to compute p-values for the probablity of each gene
# to be affected by one or a constellation of multiple types of lesions

# row.prob.subj.hit: Compute the probability that a subject has a hit for each gene.

row.prob.subj.hit=function(P,   # matrix of lesion hit probabilities, rows for genes, columns for lesion types
                           IDs) # vector of subject IDs for each lesion, length must equal ncol(P)
{
  if (length(IDs)!=ncol(P))
    stop("length(IDs) must equal ncol(P).")
  
  g=nrow(P)
  l=ncol(P)
  
  ord=order(IDs)
  IDs=IDs[ord]
  P=matrix(P[,ord],g,l)
  l=length(IDs)
  
  new.ID=which(IDs[-1]!=IDs[-l])
  ID.start=c(1,new.ID+1)
  ID.end=c(new.ID,l)
  
  n=length(ID.start)
  m=nrow(P)
  Pr=matrix(NA,m,n)
  for (i in 1:n)
  {
    pr.mtx=matrix(P[,ID.start[i]:ID.end[i]],m,ID.end[i]-ID.start[i]+1)
    Pr[,i]=1-exp(rowSums(log(1-pr.mtx)))
  }
  
  return(Pr)
}

################################################################################
# row.bern.conv: This function Compute a convolution of Bernoullis for each row of
# a Bernoulli success probability matrix.

row.bern.conv=function(P,          # matrix of lesion hit probabilities, rows for genes, columns for lesion types
                       max.x=NULL) # Maximum number of subjects or maximum number of hits
  
{
  m=nrow(P)
  n=ncol(P)
  if (is.null(max.x))
    max.x=(ncol(P))
  Pr=matrix(0,m,max.x+1)
  Pr[,1]=1
  
  for (i in 1:n)
  {
    P1=Pr*P[,i]
    P0=Pr*(1-P[,i])
    Pr0=P0
    Pr1=cbind(0,P1)
    Pr1[,max.x+1]=Pr1[,max.x+1]+Pr1[,max.x+2]
    Pr=Pr0+Pr1[,-(max.x+2)]
  }
  rs.Pr=rowSums(Pr)
  Pr=Pr/rs.Pr
  return(Pr)
}

#################################################################################
# p.order: This function compute ordered p-values and prepare the data for the
# constellation test which evaluates if the gene is affected by multiple types of lesions.

p.order=function(P) # matrix of lesion hit probabilities, rows for genes, columns for lesion types
{
  k=ncol(P)
  p.mtx=apply(P,1,sort,na.last=T)
  p.mtx=t(p.mtx)
  #n.pvals=rowSums(!is.na(p.mtx))
  
  res=p.mtx
  for (i in 1:k)
  {
    res[,i]=pbeta(p.mtx[,i],i,k-i+1)
  }
  return(res)
}

###############################################################################
# pc06.fdr: Compute FDR with the Pounds & Cheng (2006) estimator of
# the proportion of tests with a true null (pi.hat)
# Reference: https://pubmed.ncbi.nlm.nih.gov/16777905/

pc06.fdr=function(p)
  
{
  pi.hat=min(1,2*mean(p,na.rm=T))
  q=pi.hat*p.adjust(p,method="fdr")
  return(q)
}


#################################################################################
# prob.hits: The function evaluates the probability of a gene to be affected by one
# or a constellation of multiple types of lesions.

prob.hits=function(hit.cnt, # Output results of the count.hits function
                   chr.size=NULL) # A data table showing the size of the 22 autosomes, in addition to X and Y chromosomes in base pairs. data.frame should has two columns "chrom" with the chromosome number and "size" for the size of the chromosome in base pairs
{
  if (is.null(chr.size))
    chr.size=impute.chrom.size(hit.cnt$lsn.data,
                               hit.cnt$gene.data)
  
  lsn.data=hit.cnt$lsn.data
  num.lsn=unique(lsn.data$lsn.type)
  ###################
  # order and index gene.lsn.data
  ord=order(hit.cnt$gene.lsn.data$lsn.type,
            hit.cnt$gene.lsn.data$lsn.chrom,
            hit.cnt$gene.lsn.data$gene.row,
            hit.cnt$gene.lsn.data$ID)
  hit.cnt$gene.lsn.data=hit.cnt$gene.lsn.data[ord,]
  m=nrow(hit.cnt$gene.lsn.data)
  
  new.sect=which((hit.cnt$gene.lsn.data$gene.chrom[-1]!=hit.cnt$gene.lsn.data$gene.chrom[-m])|
                   (hit.cnt$gene.lsn.data$lsn.type[-1]!=hit.cnt$gene.lsn.data$lsn.type[-m])|
                   (hit.cnt$gene.lsn.data$gene.row[-1]!=hit.cnt$gene.lsn.data$gene.row[-m]))
  sect.start=c(1,new.sect+1)
  sect.end=c(new.sect,m)
  gene.lsn.index=cbind.data.frame(lsn.type=hit.cnt$gene.lsn.data$lsn.type[sect.start],
                                  chrom=hit.cnt$gene.lsn.data$gene.chrom[sect.start],
                                  gene.row=hit.cnt$gene.lsn.data$gene.row[sect.start],
                                  row.start=sect.start,
                                  row.end=sect.end,
                                  n.lsns=sect.end-sect.start+1)
  k=nrow(gene.lsn.index)
  
  new.chr=which(gene.lsn.index$chrom[-1]!=gene.lsn.index$chrom[-k])
  chr.start=c(1,new.chr+1)
  chr.end=c(new.chr,k)
  gene.lsn.chr.index=cbind.data.frame(lsn.type=gene.lsn.index$lsn.type[chr.start],
                                      chrom=gene.lsn.index$chrom[chr.start],
                                      row.start=chr.start,
                                      row.end=chr.end,
                                      n.rows=chr.end-chr.start+1)
  
  nr.li=nrow(hit.cnt$lsn.index)
  new.chr=which((hit.cnt$lsn.index$lsn.type[-1]!=hit.cnt$lsn.index$lsn.type[-nr.li])|
                  (hit.cnt$lsn.index$chrom[-1]!=hit.cnt$lsn.index$chrom[-nr.li]))
  chr.start=c(1,new.chr+1)
  chr.end=c(new.chr,nr.li)
  lsn.chr.index=cbind.data.frame(lsn.type=hit.cnt$lsn.index$lsn.type[chr.start],
                                 chrom=hit.cnt$lsn.index$chrom[chr.start],
                                 row.start=chr.start,
                                 row.end=chr.end)
  
  b=nrow(gene.lsn.chr.index)
  g=nrow(hit.cnt$nhit.mtx)
  nlt=ncol(hit.cnt$nhit.mtx)
  
  p.nsubj=p.nhit=matrix(1,g,nlt)
  colnames(p.nsubj)=colnames(p.nhit)=colnames(hit.cnt$nhit.mtx)
  
  for (i in 1:b)
  {
    # find rows for affected genes
    gli.start.row=gene.lsn.chr.index$row.start[i]
    gli.end.row=gene.lsn.chr.index$row.end[i]
    gld.start.row=gene.lsn.index$row.start[gli.start.row]
    gld.end.row=gene.lsn.index$row.end[gli.end.row]
    gld.rows=gld.start.row:gld.end.row
    gene.rows=unique(hit.cnt$gene.lsn.data$gene.row[gld.rows])
    n.genes=length(gene.rows)
    
    # find rows for lesions of this type on this chromosomes
    lsn.chr.mtch=which((lsn.chr.index$lsn.type==gene.lsn.chr.index$lsn.type[i])&
                         (lsn.chr.index$chrom==gene.lsn.chr.index$chrom[i]))
    if (length(lsn.chr.mtch)>0)
    {
      lsn.index.start.row=lsn.chr.index$row.start[lsn.chr.mtch]
      lsn.index.end.row=lsn.chr.index$row.end[lsn.chr.mtch]
      lsn.start.row=hit.cnt$lsn.index$row.start[lsn.index.start.row]
      lsn.end.row=hit.cnt$lsn.index$row.end[lsn.index.end.row]
      lsn.rows=lsn.start.row:lsn.end.row
      n.lsns=length(lsn.rows)
      lsn.type=hit.cnt$lsn.data$lsn.type[lsn.start.row]
      
      message(paste0("Computing p-values for ",
                     n.genes," gene(s) on chromosome ",
                     gene.lsn.chr.index$chrom[i],
                     " affected by ",
                     n.lsns," ",lsn.type,
                     " (data block ",i," of ",b,"): ",date()))
      
      # find chromosome size
      chr.mtch=which(gene.lsn.chr.index$chrom[i]==chr.size$chrom)
      chrom.size=chr.size$size[chr.mtch]
      
      # obtain gene sizes, lesion sizes, and gene hit probabilities
      lsn.size=hit.cnt$lsn.data$loc.end[lsn.rows]-hit.cnt$lsn.data$loc.start[lsn.rows]+1
      gene.size=hit.cnt$gene.data$loc.end[gene.rows]-hit.cnt$gene.data$loc.start[gene.rows]+1
      
      log.pr=log(rep(lsn.size,each=n.genes)+rep(gene.size,times=n.lsns))-log(chrom.size)
      pr.gene.hit=matrix(exp(log.pr),n.genes,n.lsns)
      pr.gene.hit[pr.gene.hit>1]=1
      
      lsn.subj.IDs=hit.cnt$lsn.data$ID[lsn.rows]
      pr.subj=row.prob.subj.hit(pr.gene.hit,lsn.subj.IDs)
      
      max.nsubj=max(hit.cnt$nsubj.mtx[gene.rows,lsn.type])
      max.nhit=max(hit.cnt$nhit.mtx[gene.rows,lsn.type])
      
      pr.nhit=row.bern.conv(pr.gene.hit,max.nhit)
      pr.nsubj=row.bern.conv(pr.subj,max.nsubj)
      
      for (j in 1:n.genes)
      {
        nsubj=hit.cnt$nsubj.mtx[gene.rows[j],lsn.type]
        nhit=hit.cnt$nhit.mtx[gene.rows[j],lsn.type]
        p.nsubj[gene.rows[j],lsn.type]=sum(pr.nsubj[j,(nsubj+1):(max.nsubj+1)])
        p.nhit[gene.rows[j],lsn.type]=sum(pr.nhit[j,(nhit+1):(max.nhit+1)])
      }
    }
    
  }
  
  rownames(p.nhit)=rownames(hit.cnt$nhit.mtx)
  rownames(p.nsubj)=rownames(hit.cnt$nsubj.mtx)
  colnames(hit.cnt$nhit.mtx)=paste0("nhit.",colnames(hit.cnt$nhit.mtx))
  colnames(hit.cnt$nsubj.mtx)=paste0("nsubj.",colnames(hit.cnt$nsubj.mtx))
  colnames(p.nhit)=paste0("p.",colnames(hit.cnt$nhit.mtx))
  colnames(p.nsubj)=paste0("p.",colnames(hit.cnt$nsubj.mtx))
  
  # Compute q-values
  message(paste0("Computing q-values: ",date()))
  q.nhit=p.nhit
  q.nsubj=p.nsubj
  for (i in 1:ncol(q.nhit))
  {
    pi.hat=min(1,2*mean(p.nhit[,i],na.rm=T))
    q.nhit[,i]=pi.hat*p.adjust(p.nhit[,i],method="fdr")
    pi.hat=min(1,2*mean(p.nsubj[,i],na.rm=T))
    q.nsubj[,i]=pi.hat*p.adjust(p.nsubj[,i],method="fdr")
  }
  
  colnames(q.nhit)=paste0("q.",colnames(hit.cnt$nhit.mtx))
  colnames(q.nsubj)=paste0("q.",colnames(hit.cnt$nsubj.mtx))
  
  gd.clms=setdiff(colnames(hit.cnt$gene.data),c("glp.row.start","glp.row.end"))
  lsn.clms=setdiff(colnames(hit.cnt$lsn.data),c("glp.row.start","glp.row.end"))
  
  gd.clms=c("gene.row",setdiff(gd.clms,"gene.row"))
  lsn.clms=c("lsn.row",setdiff(lsn.clms,"lsn.row"))
  
  if (length(num.lsn) >1) {
    # Now get ordered p-values
    message(paste0("Computing p-values for number of lesion types affecting genes: ",date()))
    p.ord.nhit=p.order(p.nhit)
    colnames(p.ord.nhit)=paste0("p",1:ncol(p.nhit),".nhit")
    
    p.ord.nsubj=p.order(p.nsubj)
    colnames(p.ord.nsubj)=paste0("p",1:ncol(p.nsubj),".nsubj")
    
    # q-values of ordered p-values
    q.ord.nhit=p.ord.nhit
    q.ord.nsubj=p.ord.nsubj
    message(paste0("Computing q-values for number of lesion types affecting genes: ",date()))
    for (i in 1:ncol(p.ord.nhit))
    {
      pi.hat=min(1,2*mean(p.ord.nhit[,i],na.rm=T))
      q.ord.nhit[,i]=pi.hat*p.adjust(p.ord.nhit[,i],method="fdr")
      pi.hat=min(1,2*mean(p.ord.nsubj[,i],na.rm=T))
      q.ord.nsubj[,i]=pi.hat*p.adjust(p.ord.nsubj[,i],method="fdr")
    }
    colnames(q.ord.nsubj)=paste0("q",1:ncol(p.nsubj),".nsubj")
    colnames(q.ord.nhit)=paste0("q",1:ncol(p.nhit),".nhit")
    
    gene.res=cbind.data.frame(hit.cnt$gene.data[,gd.clms],
                              hit.cnt$nsubj.mtx,
                              p.nsubj,
                              q.nsubj,
                              p.ord.nsubj,
                              q.ord.nsubj,
                              hit.cnt$nhit.mtx,
                              p.nhit,
                              q.nhit,
                              p.ord.nhit,
                              q.ord.nhit)
  }
  
  # if lesion data has only one type of lesion, skip the constellation test for locus to be affected by multiple types of lesions
  else if (length(num.lsn)==1) {
    
    gene.res=cbind.data.frame(hit.cnt$gene.data[,gd.clms],
                              hit.cnt$nsubj.mtx,
                              p.nsubj,
                              q.nsubj,
                              hit.cnt$nhit.mtx,
                              p.nhit,
                              q.nhit)
  }
  
  res=list(gene.hits=gene.res, # A data table of GRIN results that includes gene annotation, number of subjects and number of hits affecting each locus. in addition, p and FDR adjusted q-values showing the probability of each locus being affected by one or a constellation of multiple types of lesions are also included in the GRIN results.
           lsn.data=hit.cnt$lsn.data[,lsn.clms], # Input lesion data
           gene.data=hit.cnt$gene.data[,gd.clms], # Input gene data
           gene.lsn.data=hit.cnt$gene.lsn.data, # Each row represent a gene overlapped by a certain lesion. Column "gene" shows the overlapped gene and ID column has the patient ID
           chr.size=chr.size, # A data table showing the size of the 22 autosomes, in addition to X and Y chromosomes in base pairs
           gene.index=hit.cnt$gene.index, # Data.frame that shows row start and row end for each chromosome in the gene.lsn.data table
           lsn.index=hit.cnt$lsn.index) # Data.frame that shows row start and row end for each lesion in the gene.lsn.data table
  
  return(res)
}

##############################################################################
# grin.stats: The function run the Genomic Random Interval (GRIN) analysis to
# determine whether a certain locus has an abundance of lesions that is
# statistically significant.

grin.stats=function(lsn.data,            # data.frame with five columns "ID" which is the subject identifier), "chrom" which is the chromosome on which the lesion is located), "loc.start" the lesion start position, "loc.end" the lesion end position, and "lsn.type" which is the lesion category for example gain, mutation, etc..)
                    gene.data=NULL,      # data.frame with four required columns "gene" which is the gene ensembl ID, "chrom" which is the chromosome on which the gene is located, "loc.start" the gene start position, and "loc.end" which is the gene end position
                    chr.size=NULL,       # data.frame with two columns "chrom" which is the chromosome number and "size" which is the size of the chromosome in base pairs)
                    genome.version=NULL) # character string with genome version. Currently, four genome assemblies are supported including "Human_GRCh38", "Human_GRCh37", "Mouse_HGCm39", and "Mouse_HGCm38"
  
{
  if (is.null(genome.version)&&(is.null(gene.data)||is.null(chr.size)))
  {
    genome.version=select.list(c("Human_GRCh38",
                                 "Human_GRCh37",
                                 "Mouse_HGCm39",
                                 "Mouse_HGCm38"))
  }
  
  if (is.character(genome.version))
  {
    ensembl.data=get.ensembl.annotation(genome.version)
    if (is.null(gene.data))
    {
      gene.data=ensembl.data$gene.annotation
    }
    chrom.size=get.chrom.length(genome.version)
    if (is.null(chr.size))
    {
      chr.size=chrom.size
    }
  }
  
  prep.data=prep.gene.lsn.data(lsn.data,
                               gene.data)
  find.overlap=find.gene.lsn.overlaps(prep.data)
  hit.cnt=count.hits(find.overlap)
  hit.pvals=prob.hits(hit.cnt,
                      chr.size)
  return(hit.pvals)
}

##################################################################################
# write.grin.xlsx: This function Write GRIN results to an excel file with multiple
# tabs that include GRIN results, lesion data, gene annotation data, chromosome size,
# gene-lesion overlap and methods paragraph.

write.grin.xlsx=function(grin.result, # output of the grin.stats function
                         output.file) # output file name ".xlsx"
{
  saf=options()$stringsAsFactors
  options(stringsAsFactors = F)
  
  rpt.res=grin.result[c("gene.hits",
                        "gene.lsn.data",
                        "lsn.data",
                        "gene.data",
                        "chr.size")]
  
  ################################
  # Contents of the data sheets
  sheet.int=cbind.data.frame(sheet.name=c("gene.hits",
                                          "gene.lsn.data",
                                          "lsn.data",
                                          "gene.data",
                                          "chr.size"),
                             col.name="entire sheet",
                             meaning=c("GRIN statistical results",
                                       "gene-lesion overlaps",
                                       "input lesion data",
                                       "input gene location data",
                                       "input chromosome size data"))
  
  ######################################
  # Interpretation of gene hit stats
  gh.cols=colnames(rpt.res[["gene.hits"]])
  genehit.int=cbind.data.frame(sheet.name="gene.hits",
                               col.name=gh.cols,
                               meaning=gh.cols)
  rownames(genehit.int)=gh.cols
  rpt.clms=c("gene.row","gene","loc.start","loc.end")
  rpt.defs=c("gene data row index",
             "gene name",
             "locus of left edge of gene",
             "locus of right edge of gene")
  
  rpt.indx=rep(NA,length(rpt.clms))
  for (i in 1:length(rpt.clms))
    rpt.indx[i]=which(genehit.int$col.name==rpt.clms[i])
  
  genehit.int$meaning[rpt.indx]=rpt.defs
  
  nsubj.clms=which(substring(gh.cols,1,5)=="nsubj")
  genehit.int$meaning[nsubj.clms]=paste0("Number of subjects with a ",
                                         substring(gh.cols[nsubj.clms],7)," lesion ",
                                         genehit.int$meaning[nsubj.clms])
  
  p.nsubj.clms=which(substring(gh.cols,1,7)=="p.nsubj")
  genehit.int$meaning[p.nsubj.clms]=paste0("p-value for the number of subjects with a ",
                                           substring(gh.cols[p.nsubj.clms],9)," lesion ",
                                           genehit.int$meaning[p.nsubj.clms])
  
  q.nsubj.clms=which(substring(gh.cols,1,7)=="q.nsubj")
  genehit.int$meaning[q.nsubj.clms]=paste0("FDR estimate for the number of subjects with a ",
                                           substring(gh.cols[q.nsubj.clms],9)," lesion ",
                                           genehit.int$meaning[q.nsubj.clms])
  
  nhit.clms=which(substring(gh.cols,1,4)=="nhit")
  genehit.int$meaning[nhit.clms]=paste0("Number of ",
                                        substring(gh.cols[nhit.clms],6)," lesions ",
                                        genehit.int$meaning[nhit.clms])
  
  p.nhit.clms=which(substring(gh.cols,1,6)=="p.nhit")
  genehit.int$meaning[p.nhit.clms]=paste0("p-value for the number of ",
                                          substring(gh.cols[p.nhit.clms],8)," lesions ",
                                          genehit.int$meaning[p.nhit.clms])
  
  q.nhit.clms=which(substring(gh.cols,1,6)=="q.nhit")
  genehit.int$meaning[q.nhit.clms]=paste0("FDR estimate for the number of ",
                                          substring(gh.cols[q.nhit.clms],8)," lesions ",
                                          genehit.int$meaning[q.nhit.clms])
  
  
  p1.nsubj.clms=paste0("p",1:length(nsubj.clms),".nsubj")
  genehit.int[p1.nsubj.clms,"meaning"]=paste0("p-value for the number of subjects ",
                                              "with any ",1:length(nsubj.clms)," type(s) of lesion ",
                                              "overlapping the gene locus")
  
  q1.nsubj.clms=paste0("q",1:length(nsubj.clms),".nsubj")
  genehit.int[q1.nsubj.clms,"meaning"]=paste0("FDR estimate for the number of subjects ",
                                              "with any ",1:length(nsubj.clms)," type(s) of lesion ",
                                              "overlapping the gene locus")
  
  p1.nhit.clms=paste0("p",1:length(nhit.clms),".nhit")
  genehit.int[p1.nhit.clms,"meaning"]=paste0("p-value for the number of ",
                                             "any ",1:length(nhit.clms)," type(s) of lesion ",
                                             "overlapping the gene locus")
  
  q1.nhit.clms=paste0("q",1:length(nhit.clms),".nhit")
  genehit.int[q1.nhit.clms,"meaning"]=paste0("FDR estimate for the number of",
                                             " any ",1:length(nhit.clms)," type(s) of lesion ",
                                             "overlapping the gene locus")
  
  #######################################
  # Lesion data interpretation
  lsn.clms=colnames(rpt.res[["lsn.data"]])
  lsn.int=cbind.data.frame(sheet.name="lsn.data",
                           col.name=lsn.clms,
                           meaning=lsn.clms)
  
  rpt.clms=c("ID","chrom","loc.start","loc.end","lsn.type")
  int.clms=c(which(lsn.int[,"col.name"]=="ID"),
             which(lsn.int[,"col.name"]=="chrom"),
             which(lsn.int[,"col.name"]=="loc.start"),
             which(lsn.int[,"col.name"]=="loc.end"),
             which(lsn.int[,"col.name"]=="lsn.type"))
  lsn.int[int.clms,"meaning"]=c("Input Subject Identifier",
                                "Input Chromosome",
                                "Input Gene Locus Left Edge",
                                "Input Gene Locus Right Edge",
                                "Input Lesion Type")
  
  ############################################
  # Gene data interpretation
  gene.clms=colnames(rpt.res[["gene.data"]])
  gene.int=cbind.data.frame(sheet.name="gene.data",
                            col.name=gene.clms,
                            meaning="Echoed from Input")
  rpt.clms=c("gene","chrom","loc.start","loc.end",
             "glp.row.start","glp.row.end")
  int.clms=c(which(gene.int[,"col.name"]=="gene"),
             which(gene.int[,"col.name"]=="chrom"),
             which(gene.int[,"col.name"]=="loc.start"),
             which(gene.int[,"col.name"]=="loc.end"))
  gene.int[int.clms,"meaning"]=c("Input Gene Locus Name",
                                 "Input Gene Locus Chromosome",
                                 "Input Gene Locus Left Edge",
                                 "Input Gene Locus Right Edge")
  
  ###############################################
  # Chromosome size data interpretation
  chr.clms=colnames(rpt.res[["chr.size"]])
  chr.int=cbind.data.frame(sheet.name="chr.size",
                           col.name=chr.clms,
                           meaning="Echoed from Input")
  int.clms=c(which(chr.int[,"col.name"]=="chrom"),
             which(chr.int[,"col.name"]=="size"))
  chr.int[int.clms,"meaning"]=c("Input Chromosome",
                                "Input Chromosome Size")
  
  
  rpt.res$interpretation=rbind.data.frame(sheet.int,
                                          genehit.int,
                                          lsn.int,
                                          gene.int,
                                          chr.int)
  
  
  #########################################
  # Write the methods paragraph
  lsn.types=sort(unique(rpt.res[["lsn.data"]][,"lsn.type"]))
  mp=c("The genomic random interval model [ref 1] was used to evaluate ",
       "the statistical significance of the number of subjects with ",
       paste0("each type of lesion (",
              paste0(lsn.types,collapse=","),")"),
       "in each gene.  For each type of lesion, robust false discovery ",
       "estimates were computed from p-values using Storey's q-value [ref 2] ",
       "with the Pounds-Cheng estimator of the proportion of hypothesis ",
       "tests with a true null [ref 3].  Additionally, p-values for the ",
       paste0("number of subjects with any 1 to ",length(lsn.types)," types of lesions "),
       "were computed using the beta distribution derived for order statistics ",
       "of the uniform distribution [ref 4].",
       "",
       "REFERENCES",
       "[ref 1] Pounds S, et al (2013)  A genomic random interval model for statistical analysis of genomic lesion data.  Bioinformatics, 2013 Sep 1;29(17):2088-95 (PMID: 23842812).",
       "[ref 2] Storey J (2002).  A direct approach to false discovery rates.  Journal of the Royal Statistical Society Series B.  64(3): 479-498.  (doi.org/10.1111/1467-9868.00346).",
       "[ref 3] Pounds S and Cheng C (2005)  Robust estimation of the false discovery rate.  Bioinformatics 22(16): 1979-87.  (PMID: 16777905).  ",
       "[ref 4] Casella G and Berger RL (1990)  Statistical Inference.  Wadsworth & Brooks/Cole: Pacific Grove, California.  Example 5.5.1.")
  
  rpt.res$methods.paragraph=cbind.data.frame(methods.paragraph=mp)
  
  m=length(rpt.res)
  for (i in 1:m)
  {
    rpt.res[[i]]=as.data.frame(rpt.res[[i]])
  }
  
  write_xlsx(rpt.res,output.file)
  
  options(stringsAsFactors=saf)
  
  return(invisible())
  
}

############################################################################
# prep.binary.lsn.mtx: The function prepare a lesion matrix with each gene affected 
# by certain lesion type as a row and each patient as a column.

prep.binary.lsn.mtx=function(ov.data,     # output of the find.gene.lsn.overlaps function
                             min.ngrp=0)  # if specified, rows with number of patients affected by this type of lesion that's less than the specified number will be discarded (default is 0; will return all genes affected by this specified type of lesion in at least one patient).
  
{
  gene.lsn=ov.data$gene.lsn.hits
  
  # Order and index data by gene and lesion type
  gene.lsn.type=paste0(gene.lsn$gene,"_",gene.lsn$lsn.type)
  ord=order(gene.lsn$gene,gene.lsn$lsn.type)
  gene.lsn=gene.lsn[ord,]
  m=nrow(gene.lsn)
  new.gene.lsn=which((gene.lsn$gene[-1]!=gene.lsn$gene[-m])|
                       (gene.lsn$lsn.type[-1]!=gene.lsn$lsn.type[-m]))
  
  row.start=c(1,new.gene.lsn+1)
  row.end=c(new.gene.lsn,m)
  gene.lsn.indx=gene.lsn$gene.lsn[row.start]
  gene.lsn.index=cbind.data.frame(gene=gene.lsn$gene[row.start],
                                  lsn.type=gene.lsn$lsn.type[row.start],
                                  row.start=row.start,
                                  row.end=row.end)
  k=nrow(gene.lsn.index)
  uniq.ID=unique(gene.lsn$ID)
  n=length(uniq.ID)
  gene.lsn.mtx=matrix(0,k,n)
  colnames(gene.lsn.mtx)=as.character(uniq.ID)
  rownames(gene.lsn.mtx)=paste0(gene.lsn.indx$gene,"_",
                                gene.lsn.index$lsn.type)
  
  for (i in 1:k)
  {
    rows=(gene.lsn.index$row.start[i]:gene.lsn.index$row.end[i])
    ids=as.character(gene.lsn$ID[rows])
    gene.lsn.mtx[i,ids]=1
  }
  
  n.hit=rowSums(gene.lsn.mtx)
  min.n=pmin(n.hit,n-n.hit)
  
  keep.row=which(min.n>=min.ngrp)
  gene.lsn.mtx=matrix(gene.lsn.mtx[keep.row,],
                      length(keep.row),n)
  colnames(gene.lsn.mtx)=as.character(uniq.ID)
  rownames(gene.lsn.mtx)=paste0(gene.lsn.index$gene,"_",
                                gene.lsn.index$lsn.type)[keep.row]
  
  
  return(gene.lsn.mtx)
}


#############################################################################
# prep.lsn.type.matrix: The function prepare a lesion matrix with all types of lesions
#  affecting one gene as a row and each patient as a column.

prep.lsn.type.matrix=function(ov.data,    # output of the find.gene.lsn.overlaps function
                              min.ngrp=0) # if specified, genes with number of patients affected by all different types of lesions that's less than the specified number will be discarded (default is 0; will return all genes affected by any type of lesions in at least one patient).
{
  gene.lsn=ov.data$gene.lsn.hits
  
  # order and index gene-lesion overlaps by subject ID and gene
  ord=order(gene.lsn$ID,
            gene.lsn$gene)
  gene.lsn=gene.lsn[ord,]
  m=nrow(gene.lsn)
  new.block=which((gene.lsn$ID[-1]!=gene.lsn$ID[-m])|
                    (gene.lsn$gene[-1]!=gene.lsn$gene[-m]))
  row.start=c(1,new.block+1)
  row.end=c(new.block,m)
  ID.gene.index=cbind.data.frame(ID=gene.lsn$ID[row.start],
                                 gene=gene.lsn$gene[row.start],
                                 row.start=row.start,
                                 row.end=row.end)
  
  # initialize the result matrix
  lsn.types=sort(unique(gene.lsn$lsn.type))
  uniq.genes=unique(ID.gene.index$gene)
  uniq.IDs=unique(ID.gene.index$ID)
  
  grp.mtx=matrix("none",
                 length(uniq.genes),
                 length(uniq.IDs))
  rownames(grp.mtx)=uniq.genes
  colnames(grp.mtx)=uniq.IDs
  
  # fill in the matrix
  n.index=nrow(ID.gene.index)
  for (i in 1:n.index)
  {
    gene.lsn.rows=(ID.gene.index$row.start[i]:ID.gene.index$row.end[i])
    block.ID=ID.gene.index$ID[i]
    block.gene=ID.gene.index$gene[i]
    block.lsns=gene.lsn$lsn.type[gene.lsn.rows]
    block.lsns=unique(block.lsns)
    if (length(block.lsns)==1) grp.mtx[block.gene,block.ID]=block.lsns
    else grp.mtx[block.gene,block.ID]="multiple"
  }
  
  n.none=rowSums(grp.mtx== "none")
  lsn.grp.keep=(ncol(grp.mtx)-n.none)>=min.ngrp
  grp.mtx=grp.mtx[lsn.grp.keep,]
  
  return(grp.mtx)
}

#############################################################################
# default.grin.colors: This function specifies default colors for each lesion type
# in GRIN plots.

default.grin.colors=function(lsn.types) # Unique lesion types as specified in the lesion data file
{
  message("Computationally assigning lesion type colors for GRIN plots.")
  uniq.types=sort(unique(lsn.types))
  n.types=length(uniq.types)
  default.colors=c("black", "red","blue",
                    "olivedrab", "purple",
                   "cyan", "brown", "gold", 
                   "orange","steelblue")
  if (length(n.types)>length(default.colors))
    stop(paste0("Too many lesion types for default grin colors; please assign colors manually."))
  
  res=default.colors[1:n.types]
  names(res)=uniq.types
  return(res)
  
}

##################################################################################
# compute.gw.coordinates: This function assign plotting coordinates necessary for
# the genome-wide lesion plot.

compute.gw.coordinates=function(grin.res,    # GRIN results (output of the grin.stats function)
                                scl=1000000) # length of chromosome units in base pairs
  
{
  # Compute new coordinates for chromosomes
  cum.size=cumsum(grin.res$chr.size$size/scl)
  n.chr=nrow(grin.res$chr.size)
  grin.res$chr.size$x.start=c(0,cum.size[-n.chr])
  grin.res$chr.size$x.end=cum.size
  
  grin.res$gene.data$x.start=NA
  grin.res$gene.data$x.end=NA
  grin.res$gene.hits$x.start=NA
  grin.res$gene.hits$x.end=NA
  grin.res$lsn.data$x.start=NA
  grin.res$lsn.data$x.end=NA
  
  gene.data.ord=order(grin.res$gene.data$gene.row)
  grin.res$gene.data=grin.res$gene.data[gene.data.ord,]
  
  gene.hits.ord=order(grin.res$gene.hits$gene.row)
  grin.res$gene.hits=grin.res$gene.hits[gene.hits.ord,]
  
  lsn.ord=order(grin.res$lsn.data$lsn.row)
  grin.res$lsn.data=grin.res$lsn.data[lsn.ord,]
  
  
  ngi=nrow(grin.res$gene.index)
  for (i in 1:ngi)
  {
    chr.mtch=which(grin.res$gene.index$chrom[i]==grin.res$chr.size$chrom)
    chr.start=grin.res$chr.size$x.start[chr.mtch]
    
    gene.rows=(grin.res$gene.index$row.start[i]:grin.res$gene.index$row.end[i])
    grin.res$gene.data$x.start[gene.rows]=chr.start+grin.res$gene.data$loc.start[gene.rows]/scl
    grin.res$gene.data$x.end[gene.rows]=chr.start+grin.res$gene.data$loc.end[gene.rows]/scl
    
    grin.res$gene.hits$x.start[gene.rows]=chr.start+grin.res$gene.hits$loc.start[gene.rows]/scl
    grin.res$gene.hits$x.end[gene.rows]=chr.start+grin.res$gene.hits$loc.end[gene.rows]/scl
  }
  
  nli=nrow(grin.res$lsn.index)
  for (i in 1:nli)
  {
    chr.mtch=which(grin.res$lsn.index$chrom[i]==grin.res$chr.size$chrom)
    chr.start=grin.res$chr.size$x.start[chr.mtch]
    
    lsn.rows=(grin.res$lsn.index$row.start[i]:grin.res$lsn.index$row.end[i])
    grin.res$lsn.data$x.start[lsn.rows]=chr.start+grin.res$lsn.data$loc.start[lsn.rows]/scl
    grin.res$lsn.data$x.end[lsn.rows]=chr.start+grin.res$lsn.data$loc.end[lsn.rows]/scl
  }
  
  return(grin.res) # modified GRIN results to allow adding genome-wide coordinates
}

##################################################################################
# genomewide.lsn.plot:This function return a genome-wide lesion plot (all chromosomes) in the
# middle panel. Panel on the left shows -log10 q value for each locus colored by lesion
# types. Panel on the right has the number of subjects affected by lesions on each locus.

genomewide.lsn.plot=function(grin.res,        # GRIN results (output of the grin.stats function)
                             lsn.colors=NULL, # Lesion colors (If not provided by the user, colors will be automatically assigned using default.grin.colors function).
                             max.log10q=NULL) # Maximum log10 q value for genes in the GRIN results table to be added to the plot. If max.log10q=100, all -log10q values<100, will be adjusted to 100 in the plot
  
{
  if (!is.element("x.start",colnames(grin.res$lsn.data)))
    grin.res=compute.gw.coordinates(grin.res)
  
  grin.res$lsn.data$x.ID=as.numeric(as.factor(grin.res$lsn.data$ID))
  n=max(grin.res$lsn.data$x.ID)
  n.chr=nrow(grin.res$chr.size)
  
  # set up plotting region
  plot(c(-0.2,1.2)*n,c(0,-1.1*grin.res$chr.size$x.end[n.chr]),
       type="n",axes=F,xlab="",ylab="")
  
  # background colors for chromosomes
  rect(0,-grin.res$chr.size$x.start,
       n,-grin.res$chr.size$x.end,
       col=c("lightgray","gray"),
       border=NA)
  
  rect(-0.075*n,-grin.res$chr.size$x.start,
       -0.2*n,-grin.res$chr.size$x.end,
       col=c("lightgray","gray"),
       border=NA)
  
  rect(1.075*n,-grin.res$chr.size$x.start,
       1.2*n,-grin.res$chr.size$x.end,
       col=c("lightgray","gray"),
       border=NA)
  
  lsn.types=sort(unique(grin.res$lsn.index$lsn.type))
  
  if (is.null(lsn.colors))
  {
    lsn.colors=default.grin.colors(lsn.types)
  }
  grin.res$lsn.data$lsn.colors=lsn.colors[grin.res$lsn.data$lsn.type]
  
  grin.res$lsn.data$lsn.size=grin.res$lsn.data$x.end-grin.res$lsn.data$x.start
  ord=order(grin.res$lsn.data$lsn.size,decreasing=T)
  
  rect(grin.res$lsn.data$x.ID-1,
       -grin.res$lsn.data$x.start,
       grin.res$lsn.data$x.ID,
       -grin.res$lsn.data$x.end,
       col=grin.res$lsn.data$lsn.colors,
       border=grin.res$lsn.data$lsn.colors)
  
  text(c(n,0)[(1:n.chr)%%2+1],
       pos=c(4,2)[(1:n.chr)%%2+1],
       -(grin.res$chr.size$x.start+grin.res$chr.size$x.end)/2,
       grin.res$chr.size$chrom,
       cex=0.5)
  
  legend(n/2,-1.025*grin.res$chr.size$x.end[n.chr],
         fill=lsn.colors,cex=0.75,
         legend=names(lsn.colors),
         xjust=0.5,
         ncol=length(lsn.colors),
         border=NA,bty="n")
  
  nsubj.mtx=unlist(grin.res$gene.hits[,paste0("nsubj.",lsn.types)])
  qval.mtx=unlist(grin.res$gene.hits[,paste0("q.nsubj.",lsn.types)])
  nsubj.data=cbind.data.frame(gene=grin.res$gene.hits$gene,
                              x.start=grin.res$gene.hits$x.start,
                              x.end=grin.res$gene.hits$x.end,
                              nsubj=nsubj.mtx,
                              log10q=-log10(qval.mtx),
                              lsn.type=rep(lsn.types,each=nrow(grin.res$gene.hits)))
  nsubj.data=nsubj.data[nsubj.data$nsubj>0,]
  nsubj.data$lsn.colors=lsn.colors[nsubj.data$lsn.type]
  
  ord=order(nsubj.data$nsubj,decreasing=T)
  nsubj.data=nsubj.data[ord,]
  
  segments(1.075*n,
           -(nsubj.data$x.start+nsubj.data$x.end)/2,
           1.075*n+0.125*nsubj.data$nsubj/max(nsubj.data$nsubj)*n,
           col=nsubj.data$lsn.colors)
  
  nsubj.data$log10q[nsubj.data$log10q>max.log10q]=max.log10q
  
  ord=order(nsubj.data$log10q,decreasing=T)
  nsubj.data=nsubj.data[ord,]
  
  segments(-0.075*n,
           -(nsubj.data$x.start+nsubj.data$x.end)/2,
           -0.075*n-0.125*nsubj.data$log10q/max(nsubj.data$log10q)*n,
           col=nsubj.data$lsn.colors)
  
  text(-(0.075+0.20)*n/2,0,
       "-log10(q)",cex=0.75,
       pos=3)
  
  text((1.075+1.2)*n/2,0,
       "Subjects",cex=0.75,
       pos=3)
  
  text(c(-0.075,1.075)*n,
       -grin.res$chr.size$x.end[n.chr],
       0,cex=0.75,pos=1)
  
  text(c(-0.2,1.2)*n,
       -grin.res$chr.size$x.end[n.chr],
       c(max(nsubj.data$log10q),
         max(nsubj.data$nsubj)),
       cex=0.75,pos=1)
}

###################################################################################
# grin.barplt: Function return a stacked bar plot with number of patients affected by all
# different types of lesions in a pre-specified list of genes of interest.

grin.barplt=function(grin.res,        # GRIN results (output of the grin.stats function)
                     count.genes,     # vector with gene names of a list of genes to be added to the bar plot
                     lsn.colors=NULL) # Lesion colors (If not provided by the user, colors will be automatically assigned using default.grin.colors function).
  
  
{
  hits=grin.res$gene.hits
  
  # extract count data for all lesion types from GRIN results table
  count.clms=hits %>% dplyr:: select(starts_with("nsubj"))
  gene.name=hits$gene.name
  count.data=as.data.frame(cbind(gene.name, count.clms))
  count.genes=count.genes
  
  # limit the bar plot to the list of genes of interest
  selected.count=count.data[count.data$gene.name%in%count.genes,]
  names(selected.count) <- sub('^nsubj.', '', names(selected.count))
  count.clms=selected.count[,-1]
  count.clms=stack(count.clms)
  
  nsubj.data=cbind.data.frame(gene.name=selected.count$gene.name,
                              nsubj=count.clms$values,
                              lsn.type=count.clms$ind)
  
  nsubj.data[nsubj.data == 0] <- NA
  nsubj.data=na.omit(nsubj.data)
  
  # assign colors for lesion groups automatically if not provided by the user
  lsn.types=sort(unique(grin.res$lsn.index$lsn.type))
  
  if (is.null(lsn.colors))
  {
    lsn.colors=default.grin.colors(lsn.types)
  }
  
  # summarize the df
  nsubj.data <- nsubj.data %>%
    group_by(gene.name, lsn.type) %>%
    summarize(nsubj = sum(nsubj), .groups = "drop") %>%
    group_by(gene.name) %>%
    arrange(gene.name, lsn.type) %>%
    mutate(label_pos = cumsum(nsubj) - (nsubj / 2))
  
  # generate stacked bar plot for number of patients affected by each type of lesions in the list of genes of interest
  plot=ggplot(data = nsubj.data, 
              aes(x = reorder(gene.name, nsubj, sum)), cex.axis=10) +
    geom_bar(aes(y = nsubj, fill = lsn.type), 
             position = position_stack(reverse = T),
             stat="identity", width = .5)+
    scale_fill_manual(values=c(lsn.colors))+
    geom_text(aes(y = label_pos, label = nsubj),
              color = "white", size = 3.5) +
    labs(x="Gene Name",y="Count", fill = "Lesion Type")+
    coord_flip() + theme_classic()
  
  final.plt=plot+theme(text = element_text(size = 15))
  
  return(final.plt)
  
}

###################################################################################
# grin.gene.plot: This function will directly retrieve transcripts track for a gene of 
# interest from ensembl database, plot lesion data and add GRIN results including number 
# of subjects with each lesion type, -log10p and -log10q values of the gene 
# (regional gene plot). The function can also plot lesions on certain locus specified by 
# the chromosome number, locus start, locus end and the lesion category of interest.

grin.gene.plot=function(grin.res,          # GRIN results (output of the grin.stats function)
                        genome,            # genome assembly (hg19 or hg38)
                        gene=NULL,         # gene name (should be only specified in the regional gene plots)
                        lsn.clrs=NULL,     # Specified colors per lesion types (gene plots when gene name is specified). If not specified, colors will be automatically assigned using default.grin.colors function
                        chrom=NULL,        # chromosome number (should be only specified in the locus plots where plot.start and plot.end for the locus of interest are specified)
                        plot.start=NULL,   # start position of the locus of interest
                        plot.end=NULL,     # end position of the locus of interest
                        lesion.grp=NULL,   # lesion group of interest (should be only specified in locus plots when chrom, plot.start, plot.end are specified)
                        spec.lsn.clr=NULL, # color of the lesion of interest (should be specified when chrom, plot.start, plot.end and lesion.grp are specified)
                        extend.left=NULL,  # specified number will be used to manually align the left side of the gene transcripts track directly retrieved from ensembl database with the gene lesions track
                        extend.right=NULL, # specified number will be used to manually align the right side of the gene transcripts track directly retrieved from ensembl database with the gene lesions track
                        expand=0.0005,     # Controls ratio of the gene locus (start and end position) to the whole plot with default value = 0.0005 (setting expand=0 will only plot the gene locus from the start to the end position without any of the upstream or downstream regions of the gene)
                        hg38.transcripts=NULL, # transcripts data retrieved from annotation hub for hg38 version 110 (should be only specified if genome="hg38")
                        hg19.cytoband=NULL,    # hg19 chromosome bands start and end data retrieved from UCSC genome browser (should be only specified if genome="hg19")
                        hg38.cytoband=NULL)    # hg38 chromosome bands start and end data retrieved from UCSC genome browser (should be only specified if genome="hg38")

{
  # regional gene plot (gene name should be specified)
  if (is.character(gene))
  {
    # Find the requested gene
    gene.data=grin.res[["gene.data"]]
    gene.mtch=which(gene.data[,"gene.name"]==gene)
    
    if (length(gene)!=1)
      stop("Exactly one gene must be specified!")
    if (length(gene.mtch)==0)
      stop(paste0(gene," not found in gene.data."))
    
    if (length(gene.mtch)>1)
      stop(paste0("Multiple matches of ",gene," found in gene data."))
    
    gene.chr=gene.data[gene.mtch,"chrom"]
    gene.start=gene.data[gene.mtch,"loc.start"]
    gene.end=gene.data[gene.mtch,"loc.end"]
    gene.size=(gene.end-gene.start)+1
    
    # Find lesions on the same chromosome as the gene
    lsn.dset=order.index.lsn.data(grin.res[["lsn.data"]])
    lsn.data=lsn.dset$lsn.data
    
    lsn.types=unique(lsn.data$lsn.type)
    if (is.null(lsn.clrs))
      lsn.clrs=default.grin.colors(lsn.types)
    
    lsn.index=lsn.dset$lsn.index
    lsn.ind.mtch=which(lsn.index$chrom==gene.chr)
    if (length(lsn.ind.mtch)==0)
      stop(paste0("No lesions overlap ",gene,"."))
    
    lsn.chr.rows=NULL
    for (i in lsn.ind.mtch)
    {
      blk.rows=(lsn.index$row.start[i]:lsn.index$row.end[i])
      lsn.chr.rows=c(lsn.chr.rows,blk.rows)
    }
    lsn.chr.rows=unlist(lsn.chr.rows)
    
    lsn.chr.data=lsn.data[lsn.chr.rows,]
    
    if (any(lsn.chr.data$chrom!=gene.chr))
      stop(paste0("Error in finding lesions on same chromosome as gene ",gene,"."))
    
    # Find lesions at overlap the gene
    ov.rows=which((lsn.chr.data$loc.start<=gene.end)&(lsn.chr.data$loc.end>=gene.start))
    if (length(ov.rows)==0)
      stop(paste0("No lesions overlap gene ",gene,"."))
    
    lsn.gene=lsn.chr.data[ov.rows,]
    
    lsn.gene$size=lsn.gene$loc.end-
      lsn.gene$loc.start+1
    lsn.gene=lsn.gene[order(lsn.gene$lsn.type, lsn.gene$size), ]
    
    # define plotting data
    x.start=gene.start-expand*gene.size
    x.end=gene.end+expand*gene.size
    
    lsn.gene$index=1:nrow(lsn.gene)
    lsn.gene$subj.num=as.numeric(as.factor(lsn.gene$index))
    lsn.gene$lsn.clr=lsn.clrs[lsn.gene$lsn.type]
    lsn.gene$type.num=as.numeric(as.factor(lsn.gene$lsn.type))
    
    
    n.type=max(lsn.gene$type.num)
    n.subj=max(lsn.gene$subj.num)
    
    lsn.gene$y0=-lsn.gene$subj.num+(lsn.gene$type.num-1)/n.type
    lsn.gene$y1=-lsn.gene$subj.num+lsn.gene$type.num/n.type
    
    lsn.gene$x0=pmax(lsn.gene$loc.start,x.start)
    lsn.gene$x1=pmin(lsn.gene$loc.end,x.end)
    
    # based on how many characters in the gene name, the right side of the lesion data plot should be shifted so that the two panels will be alligned
    if (nchar(gene)==3)
    {
      plot(c(x.start-0.055*(x.end-x.start),x.end+0.13*(x.end-x.start)),
           c(+0.0000001,-1.75)*n.subj,type="n",
           main="",
           xlab="",
           ylab="",axes=F)
    }
    else if (nchar(gene)==4) {
      plot(c(x.start-0.055*(x.end-x.start),x.end+0.145*(x.end-x.start)),
           c(+0.0000001,-1.75)*n.subj,type="n",
           main="",
           xlab="",
           ylab="",axes=F)
    }
    
    else if (nchar(gene)==5) {
      plot(c(x.start-0.055*(x.end-x.start),x.end+0.165*(x.end-x.start)),
           c(+0.0000001,-1.75)*n.subj,type="n",
           main="",
           xlab="",
           ylab="",axes=F)
    }
    
    else if (nchar(gene)==6) {
      plot(c(x.start-0.055*(x.end-x.start),x.end+0.18*(x.end-x.start)),
           c(+0.0000001,-1.75)*n.subj,type="n",
           main="",
           xlab="",
           ylab="",axes=F)
    }
    else if (nchar(gene)==7) {
      plot(c(x.start-0.055*(x.end-x.start),x.end+0.195*(x.end-x.start)),
           c(+0.0000001,-1.75)*n.subj,type="n",
           main="",
           xlab="",
           ylab="",axes=F)
    }
    else if (nchar(gene)==8) {
      plot(c(x.start-0.055*(x.end-x.start),x.end+0.21*(x.end-x.start)),
           c(+0.0000001,-1.75)*n.subj,type="n",
           main="",
           xlab="",
           ylab="",axes=F)
    }
    else if (nchar(gene)>8) {
      plot(c(x.start-0.055*(x.end-x.start),x.end+0.25*(x.end-x.start)),
           c(+0.0000001,-1.75)*n.subj,type="n",
           main="",
           xlab="",
           ylab="",axes=F)
    }
    
    rect(x.start,
         -(1:n.subj),
         x.end,
         -(1:n.subj)+1,
         col=c("snow","gainsboro")[1+(1:n.subj)%%2],
         border=NA)
    
    segments(c(gene.start,gene.end),
             rep(-n.subj,2),
             c(gene.start,gene.end),
             rep(0,2),
             col="darkgray")
    
    rect(lsn.gene$x0,
         lsn.gene$y0,
         lsn.gene$x1,
         lsn.gene$y1,
         col=lsn.gene$lsn.clr,
         border=lsn.gene$lsn.clr)
    
    segments(c(gene.start,gene.end),
             rep(-n.subj,2),
             c(gene.start,gene.end),
             rep(0,2),
             col="darkgray",lty=2)
    
    text(c(gene.start,gene.end),
         -n.subj,
         c(gene.start,gene.end),
         pos=1, cex=0.75)
    
    
    lgd=legend((x.start+x.end)/2,-1.10*n.subj,
               fill=lsn.clrs,
               legend=names(lsn.clrs),
               ncol=length(lsn.clrs),
               xjust=0.4,border=NA,
               cex=0.7,bty="n")
    
    text(lgd$text$x[1]-0.025*diff(range(lgd$text$x)),
         -c(1.25,1.35,1.45)*n.subj,
         c("n","-log10p","-log10q"),pos=2, cex=0.7, font=2)
    
    gene.stats=grin.res[["gene.hits"]]
    stat.mtch=which(gene.stats$gene.name==gene)
    gene.stats=gene.stats[stat.mtch,]
    
    text(lgd$text$x,-1.25*n.subj,
         gene.stats[,paste0("nsubj.",names(lsn.clrs))],
         cex=0.7)
    text(lgd$text$x,-1.35*n.subj,
         round(-log10(gene.stats[,paste0("p.nsubj.",names(lsn.clrs))]),2),
         cex=0.7)
    text(lgd$text$x,-1.45*n.subj,
         round(-log10(gene.stats[,paste0("q.nsubj.",names(lsn.clrs))]),2),
         cex=0.7)
    
    lgd2=legend((x.start+x.end)/2,-1.5*n.subj,
                legend=paste0("const", 1:length(lsn.types), ".typ"),
                ncol=length(lsn.types),
                xjust=0.42,border=NA,
                cex=0.7,bty="n")
    
    text(lgd2$text$x[1]-0.0005*diff(range(lgd2$text$x)),
         -c(1.65,1.75)*n.subj,
         c("const.p", "const.q"),pos=2, cex=0.7, font=2)
    
    text(lgd2$text$x,-1.65*n.subj,
         round(-log10(gene.stats[,paste0("p",1:length(lsn.types), ".nsubj")]),2),
         cex=0.7, pos=4)
    text(lgd2$text$x,-1.75*n.subj,
         round(-log10(gene.stats[,paste0("q",1:length(lsn.types), ".nsubj")]),2),
         cex=0.7, pos=4)
    
    
    grid.echo()
    lesion.plt <- grid.grab()
    lesion.plt <- editGrob(lesion.plt, gp=gpar(fontsize=12))
    
    ## add the track of gene transcripts directly retrived from ensembl database based on the genome assembly
    if (genome=="hg19")
    {
      edb <- EnsDb.Hsapiens.v75      # same version of hg19 genome assembly used to retrieve annotation data in the get.ensembl.annotation function to make sure that coordinates will be consistent between gene annotation and the transcripts track files
      seqlevelsStyle(edb) <- "UCSC"
      txs <- getGeneRegionTrackForGviz(edb, filter = ~ genename == gene)
      
      ## Define the individual tracks:
      ## - Ideogram
      hg19_cytoband=hg19.cytoband
      ideo_track <- IdeogramTrack(genome = "hg19", chromosome = gene.chr, bands = hg19_cytoband)
      ## - Genome axis
      gaxis_track <- GenomeAxisTrack()
      ## - Transcripts
      gene_track <- GeneRegionTrack(txs, showId = TRUE, just.group = "right",
                                    name = "Transcripts", geneSymbol = TRUE, size = 0.5)
      
      grid.newpage()
      
      # specify plotting regions for gene lesions and transcript tracks
      pushViewport(viewport(height=0.78, width=1.2, y=0, just="bottom"))
      grid.draw(lesion.plt)
      popViewport(1)
      
      pushViewport(viewport(height=0.33, width=0.915, y=1, just="top"))
      plotTracks(list(ideo_track, gaxis_track, gene_track), add=TRUE, sizes = c(1,2,6))
      popViewport(1)
      
    }
    
    if (genome=="hg38")
    {
      gtf.v110=hg38.transcripts
      edb <- gtf.V110    # same version of hg38 genome assembly used to retrieve annotation data in the get.ensembl.annotation function to make sure that coordinates will be consistent between gene annotation and the transcripts track files
      seqlevelsStyle(edb) <- "UCSC"
      txs <- getGeneRegionTrackForGviz(edb, filter = ~ genename == gene)
      
      ## Define the individual tracks:
      ## - Ideogram
      hg38_cytoband=hg38.cytoband
      ideo_track <- IdeogramTrack(genome = "hg38", chromosome = gene.chr, bands = hg38_cytoband)
      ## - Genome axis
      gaxis_track <- GenomeAxisTrack()
      ## - Transcripts
      gene_track <- GeneRegionTrack(txs, showId = TRUE, just.group = "right",
                                    name = "Transcripts", geneSymbol = TRUE, size = 0.5)
      
      grid.newpage()
      
      # specify plotting regions for gene lesions and transcript tracks
      pushViewport(viewport(height=0.78, width=1.2, y=0, just="bottom"))
      grid.draw(lesion.plt)
      popViewport(1)
      
      pushViewport(viewport(height=0.33, width=0.915, y=1, just="top"))
      plotTracks(list(ideo_track, gaxis_track, gene_track), add=TRUE, sizes = c(1,2,6))
      popViewport(1)
      
    }
  }
  else {
    
    # if we want to specify a locus of interest instead of a gene
    locus.chr=chrom
    locus.start=plot.start
    locus.end=plot.end
    locus.size=(locus.end-locus.start)+1
    
    # Find lesions on the specified region of the chromosome
    lesion.data=grin.res[["lsn.data"]]
    lsn.type=lesion.grp
    lsn.data=lesion.data[lesion.data$lsn.type==lsn.type,]
    lsn.data=lsn.data[lsn.data$chrom==locus.chr,]
    lsn.dset=order.index.lsn.data(lsn.data)
    lsn.clr=spec.lsn.clr
    
    lsn.index=lsn.dset$lsn.index
    lsn.ind.mtch=which(lsn.index$chrom==locus.chr)
    if (length(lsn.ind.mtch)==0)
      stop(paste0("No lesions overlap"))
    
    lsn.chr.rows=NULL
    for (i in lsn.ind.mtch)
    {
      blk.rows=(lsn.index$row.start[i]:lsn.index$row.end[i])
      lsn.chr.rows=c(lsn.chr.rows,blk.rows)
    }
    lsn.chr.rows=unlist(lsn.chr.rows)
    
    lsn.chr.data=lsn.data[lsn.chr.rows,]
    lsn.chr.data=lsn.chr.data[lsn.chr.data$chrom==locus.chr,]
    
    if (any(lsn.chr.data$chrom!=locus.chr))
      stop(paste0("Error in finding lesions on same chromosome"))
    
    # Find lesions at overlap the locus
    ov.rows=which((lsn.chr.data$loc.start<=locus.end)&(lsn.chr.data$loc.end>=locus.start))
    if (length(ov.rows)==0)
      stop(paste0("No lesions overlap locus"))
    
    lsn.locus=lsn.chr.data[ov.rows,]
    
    lsn.locus$size=lsn.locus$loc.end-lsn.locus$loc.start+1
    lsn.locus=lsn.locus[order(lsn.locus$lsn.type, lsn.locus$size), ]
    
    # define plotting data
    x.start=locus.start-expand*locus.size
    x.end=locus.end+expand*locus.size
    
    lsn.locus$index=1:nrow(lsn.locus)
    lsn.locus$subj.num=as.numeric(as.factor(lsn.locus$index))
    lsn.locus$lsn.clr=lsn.clr
    lsn.locus$type.num=as.numeric(as.factor(lsn.locus$lsn.type))
    n.type=max(lsn.locus$type.num)
    n.subj=max(lsn.locus$subj.num)
    
    lsn.locus$y0=-lsn.locus$subj.num+(lsn.locus$type.num-1)/n.type
    lsn.locus$y1=-lsn.locus$subj.num+lsn.locus$type.num/n.type
    
    lsn.locus$x0=pmax(lsn.locus$loc.start,x.start)
    lsn.locus$x1=pmin(lsn.locus$loc.end,x.end)
    
    plot(c(x.start-0.05*(x.end-x.start),x.end+0.08*(x.end-x.start)),
         c(+0.1,-1.1)*n.subj,type="n",
         main="",
         xlab="",
         ylab="",axes=F)
    
    rect(x.start,
         -(1:n.subj),
         x.end,
         -(1:n.subj)+1,
         col=c("snow","gainsboro")[1+(1:n.subj)%%2],
         border=NA)
    
    segments(c(locus.start,locus.end),
             rep(-n.subj,2),
             c(locus.start,locus.end),
             rep(0,2),
             col="darkgray")
    
    rect(lsn.locus$x0,
         lsn.locus$y0,
         lsn.locus$x1,
         lsn.locus$y1,
         col=lsn.locus$lsn.clr,
         border=lsn.locus$lsn.clr)
    
    segments(c(locus.start,locus.end),
             rep(-n.subj,2),
             c(locus.start,locus.end),
             rep(0,2),
             col="darkgray",lty=2)
    
    text(c(locus.start,locus.end),
         -n.subj,
         c(locus.start,locus.end),
         pos=1, cex=0.85)
    text((locus.start+locus.end)/2,
         0,paste0(lsn.type, " _ ", "chr",locus.chr, ": ", locus.start, " - ", locus.end), pos=3,cex=1)
    
    grid.echo()
    lesion.plt <- grid.grab()
    lesion.plt <- editGrob(lesion.plt, gp=gpar(fontsize=12))
    
    ## add the transcripts track for all genes on the specified region  based on the genome assembly
    if (genome=="hg19")
    {
      edb <- EnsDb.Hsapiens.v75
      seqlevelsStyle(edb) <- "UCSC"
      txs <- getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                       end=plot.end)
      
      ## Define the individual tracks:
      ## - Ideogram
      hg19_cytoband=hg19.cytoband
      ideo_track <- IdeogramTrack(genome = "hg19", chromosome = locus.chr, bands = hg19_cytoband)
      ## - Genome axis
      gaxis_track <- GenomeAxisTrack()
      ## - Transcripts
      gene_track <- GeneRegionTrack(txs, showId = TRUE, just.group = "right",
                                    name = "Transcripts", geneSymbol = TRUE, size = 0.5)
      
      grid.newpage()
      
      # specify plotting regions for the locus lesions and transcript tracks
      pushViewport(viewport(height=0.45, width=1.2, y=0, just="bottom")) 
      grid.draw(lesion.plt)
      popViewport(1)
      
      pushViewport(viewport(height=0.6, width=0.915, y=1, just="top"))
      plotTracks(list(ideo_track, gaxis_track, gene_track), add=TRUE)
      popViewport(1)
      
    }
    
    if (genome=="hg38")
    {
      gtf.v110=hg38.transcripts
      edb <- gtf.V110
      seqlevelsStyle(edb) <- "UCSC"
      txs <- getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                       end=plot.end)
      
      ## Define the individual tracks:
      ## - Ideogram
      hg38_cytoband=hg38.cytoband
      ideo_track <- IdeogramTrack(genome = "hg38", chromosome = locus.chr, bands = hg38_cytoband)
      ## - Genome axis
      gaxis_track <- GenomeAxisTrack()
      ## - Transcripts
      gene_track <- GeneRegionTrack(txs, showId = TRUE, just.group = "right",
                                    name = "Transcripts", geneSymbol = TRUE, size = 0.5)
      
      grid.newpage()
      
      pushViewport(viewport(height=0.45, width=1.2, y=0, just="bottom"))
      grid.draw(lesion.plt)
      popViewport(1)
      
      pushViewport(viewport(height=0.6, width=0.915, y=1, just="top"))
      plotTracks(list(ideo_track, gaxis_track, gene_track), add=TRUE)
      popViewport(1)
    }
  }
}

###################################################################################
# grin.reg.plot: This function will plot lesion data and add GRIN results including number 
# of subjects with each lesion type, -log10p and -log10q values of the gene.
# The function can also plot lesions on certain locus specified by 
# the chromosome number, locus start, and the locus end.

grin.reg.plot=function(grin.res,          # GRIN results for regulatory features (output of the grin.stats function)
                       feature=NULL,      # feature ID from Ensembl regulatory build or FANTOM5 project
                       lsn.clrs=NULL,     # Specified colors per lesion types. If not specified, colors will be automatically assigned using default.grin.colors function
                       chrom=NULL,        # chromosome number (should be only specified in the locus plots where plot.start and plot.end for the locus of interest are specified)
                       plot.start=NULL,   # start position of the locus of interest
                       plot.end=NULL,     # end position of the locus of interest
                       expand=0.0005)     # Controls ratio of the gene locus (start and end position) to the whole plot with default value = 0.0005 (setting expand=0 will only plot the gene locus from the start to the end position without any of the upstream or downstream regions of the gene)
  
{
  # regional regulatory feature plot
  if (is.character(feature))
  {
    # Find the requested regulatory feature
    gene.data=grin.res[["gene.data"]]
    gene.mtch=which(gene.data[,"gene"]==feature)
    
    if (length(feature)!=1)
      stop("Exactly one feature must be specified!")
    if (length(gene.mtch)==0)
      stop(paste0(feature," not found in gene.data."))
    
    if (length(gene.mtch)>1)
      stop(paste0("Multiple matches of ",feature," found in gene data."))
    
    gene.chr=gene.data[gene.mtch,"chrom"]
    gene.start=gene.data[gene.mtch,"loc.start"]
    gene.end=gene.data[gene.mtch,"loc.end"]
    gene.size=(gene.end-gene.start)+1
    
    # Find lesions on the same chromosome as the feature
    lsn.dset=order.index.lsn.data(grin.res[["lsn.data"]])
    lsn.data=lsn.dset$lsn.data
    
    lsn.types=unique(lsn.data$lsn.type)
    if (is.null(lsn.clrs))
      lsn.clrs=default.grin.colors(lsn.types)
    
    lsn.index=lsn.dset$lsn.index
    lsn.ind.mtch=which(lsn.index$chrom==gene.chr)
    if (length(lsn.ind.mtch)==0)
      stop(paste0("No lesions overlap ",feature,"."))
    
    lsn.chr.rows=NULL
    for (i in lsn.ind.mtch)
    {
      blk.rows=(lsn.index$row.start[i]:lsn.index$row.end[i])
      lsn.chr.rows=c(lsn.chr.rows,blk.rows)
    }
    lsn.chr.rows=unlist(lsn.chr.rows)
    
    lsn.chr.data=lsn.data[lsn.chr.rows,]
    
    if (any(lsn.chr.data$chrom!=gene.chr))
      stop(paste0("Error in finding lesions on same chromosome as gene ",feature,"."))
    
    # Find lesions at overlap the regulatory feature
    ov.rows=which((lsn.chr.data$loc.start<=gene.end)&(lsn.chr.data$loc.end>=gene.start))
    if (length(ov.rows)==0)
      stop(paste0("No lesions overlap gene ",feature,"."))
    
    lsn.gene=lsn.chr.data[ov.rows,]
    
    lsn.gene$size=lsn.gene$loc.end-
      lsn.gene$loc.start+1
    lsn.gene=lsn.gene[order(lsn.gene$lsn.type, lsn.gene$size), ]
    
    # define plotting data
    x.start=gene.start-expand*gene.size
    x.end=gene.end+expand*gene.size
    
    lsn.gene$index=1:nrow(lsn.gene)
    lsn.gene$subj.num=as.numeric(as.factor(lsn.gene$index))
    lsn.gene$lsn.clr=lsn.clrs[lsn.gene$lsn.type]
    lsn.gene$type.num=as.numeric(as.factor(lsn.gene$lsn.type))
    
    
    n.type=max(lsn.gene$type.num)
    n.subj=max(lsn.gene$subj.num)
    
    lsn.gene$y0=-lsn.gene$subj.num+(lsn.gene$type.num-1)/n.type
    lsn.gene$y1=-lsn.gene$subj.num+lsn.gene$type.num/n.type
    
    lsn.gene$x0=pmax(lsn.gene$loc.start,x.start)
    lsn.gene$x1=pmin(lsn.gene$loc.end,x.end)
    
    plot(c(x.start-0.055*(x.end-x.start),x.end+0.13*(x.end-x.start)),
         c(+0.05,-1.75)*n.subj,type="n",
         main="",
         xlab="",
         ylab="",axes=F)
    
    rect(x.start,
         -(1:n.subj),
         x.end,
         -(1:n.subj)+1,
         col=c("snow","gainsboro")[1+(1:n.subj)%%2],
         border=NA)
    
    segments(c(gene.start,gene.end),
             rep(-n.subj,2),
             c(gene.start,gene.end),
             rep(0,2),
             col="darkgray")
    
    rect(lsn.gene$x0,
         lsn.gene$y0,
         lsn.gene$x1,
         lsn.gene$y1,
         col=lsn.gene$lsn.clr,
         border=lsn.gene$lsn.clr)
    
    segments(c(gene.start,gene.end),
             rep(-n.subj,2),
             c(gene.start,gene.end),
             rep(0,2),
             col="darkgray",lty=2)
    
    text(c(gene.start,gene.end),
         -n.subj,
         c(gene.start,gene.end),
         pos=1, cex=0.75)
    
    text((gene.start+gene.end)/2,
         0,feature,pos=3,cex=1)
    
    
    lgd=legend((x.start+x.end)/2,-1.10*n.subj,
               fill=lsn.clrs,
               legend=names(lsn.clrs),
               ncol=length(lsn.clrs),
               xjust=0.4,border=NA,
               cex=0.7,bty="n")
    
    text(lgd$text$x[1]-0.025*diff(range(lgd$text$x)),
         -c(1.25,1.35,1.45)*n.subj,
         c("n","-log10p","-log10q"),pos=2, cex=0.7, font=2)
    
    gene.stats=grin.res[["gene.hits"]]
    stat.mtch=which(gene.stats$gene==feature)
    gene.stats=gene.stats[stat.mtch,]
    
    text(lgd$text$x,-1.25*n.subj,
         gene.stats[,paste0("nsubj.",names(lsn.clrs))],
         cex=0.7)
    text(lgd$text$x,-1.35*n.subj,
         round(-log10(gene.stats[,paste0("p.nsubj.",names(lsn.clrs))]),2),
         cex=0.7)
    text(lgd$text$x,-1.45*n.subj,
         round(-log10(gene.stats[,paste0("q.nsubj.",names(lsn.clrs))]),2),
         cex=0.7)
    
    lgd2=legend((x.start+x.end)/2,-1.5*n.subj,
                legend=paste0("const", 1:length(lsn.types), ".typ"),
                ncol=length(lsn.types),
                xjust=0.42,border=NA,
                cex=0.7,bty="n")
    
    text(lgd2$text$x[1]-0.0005*diff(range(lgd2$text$x)),
         -c(1.65,1.75)*n.subj,
         c("const.p", "const.q"),pos=2, cex=0.7, font=2)
    
    text(lgd2$text$x,-1.65*n.subj,
         round(-log10(gene.stats[,paste0("p",1:length(lsn.types), ".nsubj")]),2),
         cex=0.7, pos=4)
    text(lgd2$text$x,-1.75*n.subj,
         round(-log10(gene.stats[,paste0("q",1:length(lsn.types), ".nsubj")]),2),
         cex=0.7, pos=4)
  }
  else {
    
    # To specify a locus of interest for the plotting purpose
    locus.chr=chrom
    locus.start=plot.start
    locus.end=plot.end
    locus.size=(locus.end-locus.start)+1
    
    # Find lesions on the same chromosome
    lsn.dset=order.index.lsn.data(grin.res[["lsn.data"]])
    lsn.data=lsn.dset$lsn.data
    
    lsn.types=unique(lsn.data$lsn.type)
    if (is.null(lsn.clrs))
      lsn.clrs=default.grin.colors(lsn.types)
    
    lsn.index=lsn.dset$lsn.index
    lsn.ind.mtch=which(lsn.index$chrom==locus.chr)
    if (length(lsn.ind.mtch)==0)
      stop(paste0("No lesions overlap ",chrom,"."))
    
    lsn.chr.rows=NULL
    for (i in lsn.ind.mtch)
    {
      blk.rows=(lsn.index$row.start[i]:lsn.index$row.end[i])
      lsn.chr.rows=c(lsn.chr.rows,blk.rows)
    }
    lsn.chr.rows=unlist(lsn.chr.rows)
    
    lsn.chr.data=lsn.data[lsn.chr.rows,]
    
    if (any(lsn.chr.data$chrom!=locus.chr))
      stop(paste0("Error in finding lesions on same chromosome."))
    
    # Find lesions at overlap the locus
    ov.rows=which((lsn.chr.data$loc.start<=locus.end)&(lsn.chr.data$loc.end>=locus.start))
    if (length(ov.rows)==0)
      stop(paste0("No lesions overlap the chromosome."))
    
    lsn.locus=lsn.chr.data[ov.rows,]
    
    lsn.locus$size=lsn.locus$loc.end-lsn.locus$loc.start+1
    lsn.locus=lsn.locus[order(lsn.locus$lsn.type, lsn.locus$size), ]
    
    # define plotting data
    x.start=locus.start-expand*locus.size
    x.end=locus.end+expand*locus.size
    
    lsn.locus$index=1:nrow(lsn.locus)
    lsn.locus$subj.num=as.numeric(as.factor(lsn.locus$index))
    lsn.locus$lsn.clr=lsn.clrs[lsn.locus$lsn.type]
    lsn.locus$type.num=as.numeric(as.factor(lsn.locus$lsn.type))
    n.type=max(lsn.locus$type.num)
    n.subj=max(lsn.locus$subj.num)
    
    lsn.locus$y0=-lsn.locus$subj.num+(lsn.locus$type.num-1)/n.type
    lsn.locus$y1=-lsn.locus$subj.num+lsn.locus$type.num/n.type
    
    lsn.locus$x0=pmax(lsn.locus$loc.start,x.start)
    lsn.locus$x1=pmin(lsn.locus$loc.end,x.end)
    
    plot(c(x.start-0.02*(x.end-x.start),x.end+0.1*(x.end-x.start)),
         c(+0.1,-1.15)*n.subj,type="n",
         main="",
         xlab="",
         ylab="",axes=F)
    
    rect(x.start,
         -(1:n.subj),
         x.end,
         -(1:n.subj)+1,
         col=c("snow","gainsboro")[1+(1:n.subj)%%2],
         border=NA)
    
    segments(c(locus.start,locus.end),
             rep(-n.subj,2),
             c(locus.start,locus.end),
             rep(0,2),
             col="darkgray")
    
    rect(lsn.locus$x0,
         lsn.locus$y0,
         lsn.locus$x1,
         lsn.locus$y1,
         col=lsn.locus$lsn.clr,
         border=lsn.locus$lsn.clr)
    
    segments(c(locus.start,locus.end),
             rep(-n.subj,2),
             c(locus.start,locus.end),
             rep(0,2),
             col="darkgray",lty=2)
    
    text(c(locus.start,locus.end),
         -n.subj,
         c(locus.start,locus.end),
         pos=1, cex=0.85)
    text((locus.start+locus.end)/2,
         0,paste0("chr",locus.chr, ": ", locus.start, " - ", locus.end), pos=3,cex=1)
    
    lgd=legend((x.start+x.end)/2,-1.1*n.subj,
               fill=lsn.clrs,
               legend=names(lsn.clrs),
               ncol=length(lsn.clrs),
               xjust=0.45,border=NA,
               cex=0.75,bty="n")
    
  }
}


#########################################################################
# chrom.lsn.plot: This function plot lesion data of a selected lesion type or all lesion 
# groups located on a certain part or the whole specified chromosome 

chrom.lsn.plot=function(grin.res,          # GRIN results (output of the grin.stats function)
                        genome,            # genome assembly (hg19 or hg38) gene=NULL,         # gene name (should be only specified in the regional gene plots)
                        lsn.clrs=NULL,     # Specified colors per lesion types (gene plots when gene name is specified). If not specified, colors will be automatically assigned using default.grin.colors function
                        chrom=NULL,        # chromosome number (should be only specified in the locus plots where plot.start and plot.end for the locus of interest are specified)
                        plot.start=NULL,   # start position of the locus of interest
                        plot.end=NULL,     # end position of the locus of interest
                        lesion.grp=NULL,   # lesion group of interest (specified to plot just one type of lesions)
                        spec.lsn.clr=NULL, # color of the lesion of interest (should be specified when lesion.grp is specified)
                        extend.left=NULL,  # specified number will be used to manually align the left side of the gene transcripts track directly retrieved from ensembl database with the gene lesions track
                        extend.right=NULL, # specified number will be used to manually align the right side of the gene transcripts track directly retrieved from ensembl database with the gene lesions track
                        expand=0.0005,     # Controls ratio of the gene locus (start and end position) to the whole plot with default value = 0.3 (setting expand=0 will only plot the gene locus from the start to the end position without any of the upstream or downstream regions of the gene)
                        hg38.transcripts=NULL, 
                        hg19.cytoband=NULL,
                        hg38.cytoband=NULL)     

{
  if (is.character(lesion.grp)) {
    # To specify a locus of interest for the plotting purpose
    locus.chr=chrom
    locus.start=plot.start
    locus.end=plot.end
    locus.size=(locus.end-locus.start)+1
    
    # Find lesions on the specified region of the chromosome
    lesion.data=grin.res[["lsn.data"]]
    lsn.type=lesion.grp
    lsn.data=lesion.data[lesion.data$lsn.type==lsn.type,]
    lsn.data=lsn.data[lsn.data$chrom==locus.chr,]
    lsn.dset=order.index.lsn.data(lsn.data)
    lsn.clr=spec.lsn.clr
    
    lsn.index=lsn.dset$lsn.index
    lsn.ind.mtch=which(lsn.index$chrom==locus.chr)
    if (length(lsn.ind.mtch)==0)
      stop(paste0("No lesions overlap"))
    
    lsn.chr.rows=NULL
    for (i in lsn.ind.mtch)
    {
      blk.rows=(lsn.index$row.start[i]:lsn.index$row.end[i])
      lsn.chr.rows=c(lsn.chr.rows,blk.rows)
    }
    lsn.chr.rows=unlist(lsn.chr.rows)
    
    lsn.chr.data=lsn.data[lsn.chr.rows,]
    lsn.chr.data=lsn.chr.data[lsn.chr.data$chrom==locus.chr,]
    
    if (any(lsn.chr.data$chrom!=locus.chr))
      stop(paste0("Error in finding lesions on same chromosome"))
    
    # Find lesions at overlap the locus
    ov.rows=which((lsn.chr.data$loc.start<=locus.end)&(lsn.chr.data$loc.end>=locus.start))
    if (length(ov.rows)==0)
      stop(paste0("No lesions overlap locus"))
    
    lsn.locus=lsn.chr.data[ov.rows,]
    
    lsn.locus$size=lsn.locus$loc.end-lsn.locus$loc.start+1
    lsn.locus=lsn.locus[order(lsn.locus$lsn.type, lsn.locus$size), ]
    
    # define plotting data
    x.start=locus.start-expand*locus.size
    x.end=locus.end+expand*locus.size
    
    lsn.locus$index=1:nrow(lsn.locus)
    lsn.locus$subj.num=as.numeric(as.factor(lsn.locus$index))
    lsn.locus$lsn.clr=lsn.clr
    lsn.locus$type.num=as.numeric(as.factor(lsn.locus$lsn.type))
    n.type=max(lsn.locus$type.num)
    n.subj=max(lsn.locus$subj.num)
    
    lsn.locus$y0=-lsn.locus$subj.num+(lsn.locus$type.num-1)/n.type
    lsn.locus$y1=-lsn.locus$subj.num+lsn.locus$type.num/n.type
    
    lsn.locus$x0=pmax(lsn.locus$loc.start,x.start)
    lsn.locus$x1=pmin(lsn.locus$loc.end,x.end)
    
    plot(c(x.start-0.02*(x.end-x.start),x.end+0.1*(x.end-x.start)),
         c(+0.1,-1.1)*n.subj,type="n",
         main="",
         xlab="",
         ylab="",axes=F)
    
    rect(x.start,
         -(1:n.subj),
         x.end,
         -(1:n.subj)+1,
         col=c("snow","gainsboro")[1+(1:n.subj)%%2],
         border=NA)
    
    segments(c(locus.start,locus.end),
             rep(-n.subj,2),
             c(locus.start,locus.end),
             rep(0,2),
             col="darkgray")
    
    rect(lsn.locus$x0,
         lsn.locus$y0,
         lsn.locus$x1,
         lsn.locus$y1,
         col=lsn.locus$lsn.clr,
         border=lsn.locus$lsn.clr)
    
    segments(c(locus.start,locus.end),
             rep(-n.subj,2),
             c(locus.start,locus.end),
             rep(0,2),
             col="darkgray",lty=2)
    
    text(c(locus.start,locus.end),
         -n.subj,
         c(locus.start,locus.end),
         pos=1, cex=0.85)
    text((locus.start+locus.end)/2,
         0,paste0(lsn.type, " _ ", "chr",locus.chr, ": ", locus.start, " - ", locus.end), pos=3,cex=1)
    
    grid.echo()
    lesion.plt <- grid.grab()
    lesion.plt <- editGrob(lesion.plt, gp=gpar(fontsize=12))
    
    ## add the transcripts track for all genes on the specified region  based on the genome assembly
    if (genome=="hg19")
    {
      edb <- EnsDb.Hsapiens.v75
      seqlevelsStyle(edb) <- "UCSC"
      txs <- getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                       end=plot.end)
      
      ## Define the individual tracks:
      ## - Ideogram
      hg19_cytoband=hg19.cytoband
      ideo_track <- IdeogramTrack(genome = "hg19", chromosome = locus.chr, 
                                  bands = hg19_cytoband, from=plot.start, to=plot.end)
      
      grid.newpage()
      
      # specify plotting regions for the locus lesions and transcript tracks
      pushViewport(viewport(height=1, width=1.25, y=0, just="bottom")) 
      grid.draw(lesion.plt)
      popViewport(1)
      
      pushViewport(viewport(height=0.15, width=0.88, y=1, just="top"))
      plotTracks(ideo_track, from = plot.start, to = plot.end, add=TRUE, showId = FALSE,
                 showBandId = TRUE, cex.bands = 0.4)
      popViewport(1)
      
    }
    
    if (genome=="hg38")
    {
      gtf.v110=hg38.transcripts
      edb <- gtf.V110
      seqlevelsStyle(edb) <- "UCSC"
      txs <- getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                       end=plot.end)
      
      ## Define the individual tracks:
      ## - Ideogram
      hg38_cytoband=hg38.cytoband
      ideo_track <- IdeogramTrack(genome = "hg38", chromosome = locus.chr, 
                                  bands = hg38_cytoband, from=plot.start, to=plot.end)
      
      grid.newpage()
      
      pushViewport(viewport(height=1, width=1.25, y=0, just="bottom"))
      grid.draw(lesion.plt)
      popViewport(1)
      
      pushViewport(viewport(height=0.15, width=0.88, y=1, just="top"))
      plotTracks(ideo_track, from = plot.start, to = plot.end, add=TRUE, showId = FALSE,
                 showBandId = TRUE, cex.bands = 0.4)
      popViewport(1)
    }
  }
  else if (is.null(lesion.grp)) {
    # To specify a locus of interest for the plotting purpose
    locus.chr=chrom
    locus.start=plot.start
    locus.end=plot.end
    locus.size=(locus.end-locus.start)+1
    
    # Find lesions on the same chromosome as the gene
    lsn.dset=order.index.lsn.data(grin.res[["lsn.data"]])
    lsn.data=lsn.dset$lsn.data
    
    lsn.types=unique(lsn.data$lsn.type)
    if (is.null(lsn.clrs))
      lsn.clrs=default.grin.colors(lsn.types)
    
    lsn.index=lsn.dset$lsn.index
    lsn.ind.mtch=which(lsn.index$chrom==locus.chr)
    if (length(lsn.ind.mtch)==0)
      stop(paste0("No lesions overlap ",chrom,"."))
    
    lsn.chr.rows=NULL
    for (i in lsn.ind.mtch)
    {
      blk.rows=(lsn.index$row.start[i]:lsn.index$row.end[i])
      lsn.chr.rows=c(lsn.chr.rows,blk.rows)
    }
    lsn.chr.rows=unlist(lsn.chr.rows)
    
    lsn.chr.data=lsn.data[lsn.chr.rows,]
    
    if (any(lsn.chr.data$chrom!=locus.chr))
      stop(paste0("Error in finding lesions on same chromosome."))
    
    # Find lesions at overlap the locus
    ov.rows=which((lsn.chr.data$loc.start<=locus.end)&(lsn.chr.data$loc.end>=locus.start))
    if (length(ov.rows)==0)
      stop(paste0("No lesions overlap the chromosome."))
    
    lsn.locus=lsn.chr.data[ov.rows,]
    
    lsn.locus$size=lsn.locus$loc.end-lsn.locus$loc.start+1
    lsn.locus=lsn.locus[order(lsn.locus$lsn.type, lsn.locus$size), ]
    
    # define plotting data
    x.start=locus.start-expand*locus.size
    x.end=locus.end+expand*locus.size
    
    lsn.locus$index=1:nrow(lsn.locus)
    lsn.locus$subj.num=as.numeric(as.factor(lsn.locus$index))
    lsn.locus$lsn.clr=lsn.clrs[lsn.locus$lsn.type]
    lsn.locus$type.num=as.numeric(as.factor(lsn.locus$lsn.type))
    n.type=max(lsn.locus$type.num)
    n.subj=max(lsn.locus$subj.num)
    
    lsn.locus$y0=-lsn.locus$subj.num+(lsn.locus$type.num-1)/n.type
    lsn.locus$y1=-lsn.locus$subj.num+lsn.locus$type.num/n.type
    
    lsn.locus$x0=pmax(lsn.locus$loc.start,x.start)
    lsn.locus$x1=pmin(lsn.locus$loc.end,x.end)
    
    plot(c(x.start-0.02*(x.end-x.start),x.end+0.1*(x.end-x.start)),
         c(+0.1,-1.15)*n.subj,type="n",
         main="",
         xlab="",
         ylab="",axes=F)
    
    rect(x.start,
         -(1:n.subj),
         x.end,
         -(1:n.subj)+1,
         col=c("snow","gainsboro")[1+(1:n.subj)%%2],
         border=NA)
    
    segments(c(locus.start,locus.end),
             rep(-n.subj,2),
             c(locus.start,locus.end),
             rep(0,2),
             col="darkgray")
    
    rect(lsn.locus$x0,
         lsn.locus$y0,
         lsn.locus$x1,
         lsn.locus$y1,
         col=lsn.locus$lsn.clr,
         border=lsn.locus$lsn.clr)
    
    segments(c(locus.start,locus.end),
             rep(-n.subj,2),
             c(locus.start,locus.end),
             rep(0,2),
             col="darkgray",lty=2)
    
    text(c(locus.start,locus.end),
         -n.subj,
         c(locus.start,locus.end),
         pos=1, cex=0.85)
    text((locus.start+locus.end)/2,
         0,paste0("chr",locus.chr, ": ", locus.start, " - ", locus.end), pos=3,cex=1)
    
    lgd=legend((x.start+x.end)/2,-1.1*n.subj,
               fill=lsn.clrs,
               legend=names(lsn.clrs),
               ncol=length(lsn.clrs),
               xjust=0.45,border=NA,
               cex=0.75,bty="n")
    
    grid.echo()
    lesion.plt <- grid.grab()
    lesion.plt <- editGrob(lesion.plt, gp=gpar(fontsize=12))
    
    ## add the transcripts track for all genes on the specified region  based on the genome assembly
    if (genome=="hg19")
    {
      edb <- EnsDb.Hsapiens.v75
      seqlevelsStyle(edb) <- "UCSC"
      txs <- getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                       end=plot.end)
      
      ## Define the individual tracks:
      ## - Ideogram
      hg19_cytoband=hg19.cytoband
      ideo_track <- IdeogramTrack(genome = "hg19", chromosome = locus.chr, 
                                  bands = hg19_cytoband, from=plot.start, to=plot.end)
      
      grid.newpage()
      
      # specify plotting regions for the locus lesions and transcript tracks
      pushViewport(viewport(height=1, width=1.25, y=0, just="bottom")) 
      grid.draw(lesion.plt)
      popViewport(1)
      
      pushViewport(viewport(height=0.15, width=0.88, y=1, just="top"))
      plotTracks(ideo_track, from = plot.start, to = plot.end, add=TRUE, showId = FALSE,
                 showBandId = TRUE, cex.bands = 0.4)
      popViewport(1)
      
    }
    
    if (genome=="hg38")
    {
      gtf.v110=hg38.transcripts
      edb <- gtf.V110
      seqlevelsStyle(edb) <- "UCSC"
      txs <- getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                       end=plot.end)
      
      ## Define the individual tracks:
      ## - Ideogram
      hg38_cytoband=hg38.cytoband
      ideo_track <- IdeogramTrack(genome = "hg38", chromosome = locus.chr, 
                                  bands = hg38_cytoband, from=plot.start, to=plot.end)
      
      grid.newpage()
      
      pushViewport(viewport(height=1, width=1.25, y=0, just="bottom"))
      grid.draw(lesion.plt)
      popViewport(1)
      
      pushViewport(viewport(height=0.15, width=0.88, y=1, just="top"))
      plotTracks(ideo_track, from = plot.start, to = plot.end, add=TRUE, showId = FALSE,
                 showBandId = TRUE, cex.bands = 0.4)
      popViewport(1)
    }
  }
}

#########################################################################
# grin.oncoprint.mtx: Function use GRIN results table and prepare the lesion matrix that the user
# can pass to the oncoprint function from complexheat map package to geneare an OncoPrint for a 
# selcted list of genes.

grin.oncoprint.mtx=function(grin.res, # GRIN results (output of the grin.stats function)
                            oncoprint.genes) # vector of ensembl IDs for the selected list of genes
  
{
  selected=unlist(oncoprint.genes)
  selected=as.vector(selected)
  selected.genes= grin.res$gene.lsn.data[grin.res$gene.lsn.data$gene %in% selected,]
  selected.genes=selected.genes[,c(2,7,11)]  # extract patient IDs and lsn type for each gene in the selected genes list
  row.data=paste(selected.genes[,1],
                 selected.genes[,2],
                 selected.genes[,3],
                 sep="_")
  dup.data=duplicated(row.data)
  select.genes=selected.genes[!dup.data,]
  
  ord=order(select.genes$gene,
            select.genes$ID,
            select.genes$lsn.type)
  select.genes=select.genes[ord,]
  
  uniq.genes=unique(select.genes$gene)
  uniq.subj=unique(select.genes$ID)
  n.genes=length(uniq.genes)
  n.subj=length(uniq.subj)
  mtx=matrix("",n.genes,n.subj)   # create a matrix with each gene as a row
  colnames(mtx)=uniq.subj
  rownames(mtx)=uniq.genes
  
  k=nrow(select.genes)
  for (i in 1:k)
  {
    subj.id=select.genes[i,"ID"]
    gene.id=select.genes[i,"gene"]
    mtx[gene.id,subj.id]=paste0(mtx[gene.id,subj.id],
                                select.genes[i,"lsn.type"],";")
  }
  mtx=as.data.frame(mtx)
  mtx<-tibble::rownames_to_column(mtx, "ensembl.ID")
  
  gene.annotation= grin.res$gene.data
  ensembl.annotation=cbind(gene.annotation$gene, gene.annotation$gene.name)
  colnames(ensembl.annotation)=c("ensembl.ID", "gene.name")
  
  mtx.final=merge(ensembl.annotation,mtx,by="ensembl.ID", all.y=TRUE)  # add gene name
  mtx.final=mtx.final[,-1]
  rownames(mtx.final)=mtx.final[,1]
  mtx.final=mtx.final[,-1]
  
  return(mtx.final)
  
}

##############################################################################
# onco.print.props: The function order lesion types based on the average size of each type and assign the
# proportion of the oncoprint rectangle that should be color filled based on the average size of each
# lesion type. Color filled proportion of the oncoprint rectangles can be also specified by the user for 
# each lesion type based on the hgt parameter.

onco.print.props=function(lsn.data,  # data.frame with columns ID (subject identifier), chrom (chromosome on which the lesion is located), loc.start (lesion start position), loc.end (lesion end position), lsn.type (lesion category for example gain, mutation, etc..)
                          clr=NULL,  # Lesion colors (If not provided by the user, colors will be automatically assigned using default.grin.colors function).
                          hgt=NULL)  # manually assign the proportion of the oncoprint rectangle that should be color filled for each lesion group
{
  # get unique lesion types
  results <- list()
  lsn.types=unique(lsn.data$lsn.type)
  
  if (is.null(clr))
  {
    clr=default.grin.colors(sort(lsn.types))
  }
  
  lsn.data$size=lsn.data$loc.end-lsn.data$loc.start+1
  ave.lsn.size=setorder(aggregate(size~lsn.type,data=lsn.data,mean),-size)
  
  # specify proportion of the oncoprint recatangle that should be color filled based on the lesion size
  if (is.null(hgt))
  {
    ave.lsn.size$hgt = length(ave.lsn.size$lsn.type):1
  }
  else{
    hgt=as.data.frame(hgt)
    hgt=tibble::rownames_to_column(hgt, "lsn.type")
    ave.lsn.size=merge(hgt,ave.lsn.size,by="lsn.type", all.y=TRUE)
  }
  
  ave.lsn.size$wdth = rep(1,length(ave.lsn.size$lsn.type))
  ave.lsn.size=ave.lsn.size[order(ave.lsn.size$lsn.type),]
  
  twhc=cbind.data.frame(type=ave.lsn.size$lsn.type,
                        clr=default.grin.colors(sort(lsn.types)),
                        hgt=ave.lsn.size$hgt,
                        wdth=ave.lsn.size$wdth)
  
  res=onco.print.alter.func(twhc)
  res=gsub('\"',"'",res,fixed=T)
  res2=lapply(FUN=eval,X=parse(text=res))
  names(res2) = names(res)
  results$alter_func <- res2
  names(clr) <- ave.lsn.size$lsn.type
  results$col <- sort(clr)
  results$heatmap_legend_param <- oncoprint.legend(ave.lsn.size$lsn.type)
  return(results)
  
}


###################################################################
# onco.print.alter.func: The function define some important characteristics of the oncoprint such as 
# background color, rectangles height and width.

onco.print.alter.func=function(type.wdth.hgt.clr,  # data.frame with columns for lesion type, box size, box color
                               bg.clr="#CCCCCC",   # background color
                               bg.hgt=2,           # background height
                               bg.wdth=0.00005)    # background width
  
{
  res = NULL
  size = nrow(type.wdth.hgt.clr)
  bg.code=paste0("background = function(x,y,w,h){",
                 "grid.rect(x,y,w-unit(",bg.wdth,',"pt"),',
                 "h-unit(",bg.hgt,',"pt"),',
                 "gp = gpar(fill = '",bg.clr,"',col=NA))}")
  names(bg.code) <- "background"
  Rcode=paste0(type.wdth.hgt.clr[,"type"]," = function(x,y,w,h){",
               "grid.rect(x,y,w-unit(",type.wdth.hgt.clr[,"wdth"],',"pt"),',
               "h-unit(h*",(1-(type.wdth.hgt.clr[,"hgt"]*(1/(size+1)))),',"pt"),',
               "gp = gpar(fill ='",type.wdth.hgt.clr[,"clr"],"',col=NA))}")
  
  names(Rcode) = type.wdth.hgt.clr[,"type"]
  Rcode=c(bg.code,Rcode)
  return(Rcode)
}

#######################################################################
# oncoprint.legend: function specify name and colors of the oncoprint legend lesion groups

oncoprint.legend = function(lsn.type){
  res <- list()
  res$title <- 'Lesion Category'
  res$at <- lsn.type
  res$labels <- sub('^(\\w?)', '\\U\\1', lsn.type, perl=T)
  return(res)
}



#########################################################################################
# Associate Lesions with Expression (ALEX) library
########################################################################################

####################################################################
# alex.prep.lsn.expr: This function prepares lesion and expression data matrices
# for the KW.hit.express function that runs kruskal-Wallis test for the association between
# lesion groups and expression level of each gene (it only keep genes with both lesion and 
# expression data with rows ordered by ensembl ID and columns ordered by patient's ID)

alex.prep.lsn.expr=function(expr.mtx,          # Normalized log2 transformed expression data with genes in rows and subjects in columns (first column "ensembl.ID" should be gene ensembl IDs)
                            lsn.data,          # lesion data in GRIN compatible format. Data frame should has five columns that include "ID" with patient ID, "chrom" which is the chromosome on which the lesion is located, "loc.start" which is the lesion start position, "loc.end" the lesion end position and "lsn.type" which is the lesion type for example gain, loss, mutation, fusion, etc...
                            gene.annotation,   # gene annotation data either provided by the user or directly retreived from ensembl biomaRT database using get.ensembl.annotation function included in the GRIN2.0 library if the genome.version is specified. Object should has four columns "gene" which is the ensembl ID of annotated genes, "chrom" which is the chromosome on which the gene is located, "loc.start" which is the gene start position, and "loc.end" the gene end position.
                            min.pts.expr=NULL, # minimum number of patients with expression level of a gene > 0
                            min.pts.lsn=NULL)  # minimum number of patients with any type of lesions in a certain gene
  
{
  gene.lsn=prep.gene.lsn.data(lsn.data, gene.annotation)    # Prepare gene and lesion data for later computations
  gene.lsn.overlap= find.gene.lsn.overlaps(gene.lsn)        # Use the results of prep.gene.lsn.data to find lesion-gene overlaps
  
  message(paste0("Preparing gene-lesion matrix: ",date()))
  lsn.type.matrix=prep.lsn.type.matrix(gene.lsn.overlap)    # prepare lesion type matrix (just one row for each gene with all lesion types included)
  lsn.grp.mtx=as.data.frame(lsn.type.matrix)
  id.hits = colnames(lsn.grp.mtx)
  
  # prepare expression data for ALEX analysis
  expr.mtx=expr.mtx
  rownames(expr.mtx)=expr.mtx[,1]
  expr.mtx=expr.mtx[,-1]
  
  id.expr=colnames(expr.mtx)
  lsn.grp.mtx=lsn.grp.mtx[ ,(names(lsn.grp.mtx) %in% id.expr)]
  
  id.lsn.final=colnames(lsn.grp.mtx)
  expr.mtx=expr.mtx[ ,(names(expr.mtx) %in% id.lsn.final)]
  
  if (is.numeric(min.pts.expr))
  {
    # To keep only genes with expression value > 0 in at least the number of patients specified in min.pts.expr
    zero.expr=rowSums(expr.mtx==0)
    expr.grp.keep=(ncol(expr.mtx)-zero.expr)>=min.pts.expr
    expr.mtx=expr.mtx[expr.grp.keep,]
  }
  ############################
  # To exclude any gene that has no lesions in at least the number of patients specified in min.pts.lsn
  if (is.numeric(min.pts.lsn))
  {
    n.none=rowSums(lsn.grp.mtx== "none")
    lsn.grp.keep=(ncol(lsn.grp.mtx)-n.none)>=min.pts.lsn
    lsn.grp.mtx=lsn.grp.mtx[lsn.grp.keep,]
  }
  # To keep only shared set of genes between expression and lesion matrices
  lsn.genes=rownames(lsn.grp.mtx)
  expr.mtx=expr.mtx[rownames(expr.mtx) %in% lsn.genes,]
  
  expr.genes=rownames(expr.mtx)
  lsn.grp.mtx=lsn.grp.mtx[rownames(lsn.grp.mtx) %in% expr.genes,]
  
  ## order lsn matrix by gene ID first then colnames for patient IDs alphabetically
  lsn.grp.mtx=lsn.grp.mtx[order(rownames(lsn.grp.mtx)),]
  lsn.grp.mtx=lsn.grp.mtx[,order(colnames(lsn.grp.mtx))]
  
  ## order expr matrix by gene ID first then colnames for patient IDs alphabetically
  expr.mtx=expr.mtx[order(rownames(expr.mtx)),]
  expr.mtx=expr.mtx[,order(colnames(expr.mtx))]
  
  # To make sure that both lesion and expression matrices have same patients and genes order
  check.rownames=all(rownames(lsn.grp.mtx)==rownames(expr.mtx))
  if (check.rownames==FALSE)
    stop("Gene-lesion matrix rownames must match expression matrix rownames.")
  
  check.colnames=all(colnames(lsn.grp.mtx)==colnames(expr.mtx))
  if (check.colnames==FALSE)
    stop("Gene-lesion matrix patient IDs must match expression matrix patient IDs.")
  
  # To generate row.mtch file
  lsn.ensembl.id=rownames(lsn.grp.mtx)
  expr.ensembl.id=rownames(expr.mtx)
  
  row.mtch=cbind.data.frame(expr.row=expr.ensembl.id,
                            hit.row=lsn.ensembl.id)
  
  
  res=list(alex.expr=expr.mtx,      # expression data (data table with gene ensembl IDs as row names and each column is a patient)
           alex.lsn=lsn.grp.mtx,    # overlapped gene lesion data table with gene ensembl IDs as row names and each column is a patient
           alex.row.mtch=row.mtch)  # ensembl ID for one row of the alex.lsn table to be paired with one row of the expression data for the Kruskal-Wallis test
  
  return(res)
  
}

##############################################################################################
# KW.hit.express: This function uses Kruskal-Wallis test to evaluate the association
# between lesion groups and expression level of the same corresponding gene.

KW.hit.express=function(alex.data,          # output of the alex.prep.lsn.expr function (list of three data tables "alex.expr" with expression data ready for KW test, "alex.lsn" with lesion data and row.mtch)
                        gene.annotation,    # gene annotation data. First column "gene" should has ensembl IDs
                        min.grp.size=NULL)  # minimum group size to perform the test (there should be at least two groups with number of patients > min.grp.size)
  
{
  expr.matrix=as.matrix(alex.data$alex.expr)
  lsn.matrix=as.matrix(alex.data$alex.lsn)
  row.mtch=alex.data$alex.row.mtch
  
  message(paste0("Computing KW P-value: ",date()))
  p=apply(row.mtch,1,one.KW.pvalue,
          expr.mtx=expr.matrix,
          hit.grps=lsn.matrix,
          min.grp.size=min.grp.size)
  
  message(paste0("Preparing results table: ",date()))
  kw.res=cbind(row.mtch,p.KW=p)
  # To prepare and print KW test results
  
  colnames(kw.res)=c("expr.row", "gene", "p.KW")
  kw.res.annotated=merge(gene.annotation,kw.res,by="gene", all.y=TRUE)
  
  # Compute FDR adjusted q values
  kw.res.annotated$q.KW=p.adjust(kw.res.annotated$p.KW,method="fdr")
  
  # To compute number of subjects, mean, median and stdev by lesion type
  unique.grps=sort(unique(lsn.matrix))
  unique.grps=sort(unique(unique.grps))
  count.by.lesion <- sapply(unique.grps,function(x)rowSums(lsn.matrix==x)) 
  colnames(count.by.lesion) = paste(colnames(count.by.lesion),"n.subjects",sep="_")
  count.by.lesion=as.matrix(count.by.lesion)
  
  mean.by.lsn=row.stats.by.group(expr.matrix, lsn.matrix, mean)
  colnames(mean.by.lsn) = paste(colnames(mean.by.lsn),"mean",sep="_")
  mean.by.lsn=as.matrix(mean.by.lsn)
  
  median.by.lsn=row.stats.by.group(expr.matrix, lsn.matrix, median)
  colnames(median.by.lsn) = paste(colnames(median.by.lsn),"median",sep="_")
  median.by.lsn=as.matrix(median.by.lsn)
  
  sd.by.lsn=row.stats.by.group(expr.matrix, lsn.matrix, sd)
  colnames(sd.by.lsn) = paste(colnames(sd.by.lsn),"sd",sep="_")
  sd.by.lsn=as.matrix(sd.by.lsn)
  
  kw.res.final=cbind(kw.res.annotated, count.by.lesion, mean.by.lsn, median.by.lsn, sd.by.lsn)
  kw.res.final=kw.res.final[!duplicated(kw.res.final[ , "gene.name"]),]
  res=kw.res.final
  
  return(res)
  
}

########################################
# one.KW.pvalue: Function Computes the KW p-value for one row of the lesion data matrix
# paired with one row of the expression data matrix

one.KW.pvalue=function(one.row.mtch,      # one row of the row.mtch matrix (output of the alex.prep.lsn.expr function)
                       expr.mtx,          # expression data matrix (alex.expr, output of the alex.prep.lsn.expr function)
                       hit.grps,          # lesion matrix (alex.lsn, output of the alex.prep.lsn.expr function)
                       min.grp.size=NULL) # minimum group size to perform the test (there should be at least two groups with number of patients > min.grp.size)
  
{
  # extract data for the test
  expr.row=one.row.mtch["expr.row"]
  hit.row=one.row.mtch["hit.row"]
  y=expr.mtx[expr.row,]
  grps=hit.grps[hit.row,]
  
  # identify data to include in the test
  grps=define.grps(grps,min.grp.size)
  exc=(grps=="EXCLUDE")
  grps=grps[!exc]
  y=y[!exc]
  
  # return NA if there is no test to perform
  uniq.grps=unique(grps)
  if (length(uniq.grps)<2) return(NA)
  
  # perform the KW test
  kw.res=kruskal.test(y~grps)
  
  # return the KW p-value
  return(kw.res$p.value)
}


#######################################################
# define.grps: Function defines the lesion groups to be included in the Kruskal-Wallis test.

define.grps=function(grps,               # Lesion groups
                     min.grp.size=NULL)  # minimum number of patients in the lesion group to be included in the KW test.
{
  grp.tbl=table(grps)
  small.grp=(grp.tbl<min.grp.size)
  exclude.grp=grps%in%(names(grp.tbl)[small.grp])
  grps[exclude.grp]="EXCLUDE"
  return(grps)
}

###############################################
# stat.by.group: Function to compute stats for row summary (mean, median and sd)

stat.by.group=function(x,              # vector of expression data
                       g,              # vector of group labels corresponding to x (lesion groups)
                       stat,           # stats to be computed for the expression data by lesion groups (mean, median or standard deviation)
                       all.grps=NULL,  # vector of all the possible group labels
                       ...)            # additional arguments to stat
  
{
  if (is.null(all.grps))               # if all.grps isn't specified, define it from g
    all.grps=sort(unique(g))
  
  all.grps=sort(unique(all.grps))      # ensure grps don't appear twice and are in the same order for each application
  
  res=rep(NA,length(all.grps))         # initialize the result vector
  names(res)=all.grps                  # name the elements of the result vector
  
  for (i in all.grps)                  # loop over the groups
  {
    in.grp=(g%in%i)
    if (any(in.grp))
      res[i]=stat(x[in.grp],...)
  }
  
  return(res)
  
}

###################################################
# row.stats.by.group: Function uses stat.by.group function to compute stats for each row of the
# expression and lesion data matrices and return a data.frame of expression level by lesion groups for
# each specified test

row.stats.by.group=function(X,   # expression data
                            G,   # lesion data
                            stat, # stats to be computed
                            ...)
  
{
  if (any(dim(X)!=dim(G)))
    stop("X and G are of incompatible dimensions.  X and G must have the same dimensions.")
  
  if (any(colnames(X)!=colnames(G)))
    stop("X and G must have column names matched.")
  
  all.grps=sort(unique(G))
  all.grps=sort(unique(all.grps))
  
  k=length(all.grps)
  
  m=nrow(X)
  
  res=matrix(NA,m,k)
  colnames(res)=all.grps
  rownames(res)=rownames(X)
  
  for (i in 1:m)
  {
    res[i,]=stat.by.group(X[i,],G[i,],stat,all.grps,...)
  }
  
  return(res)
}

##################################################################################################
# alex.boxplots: Function return box plots for expression data by lesion groups for selected
# number of genes based on a specified q-value of the kruskal-wallis test results

alex.boxplots=function(alex.data,          # output of the alex.prep.lsn.expr function (list of three data tables "alex.expr" with expression data, "alex.lsn" with overlapped gene lesion data and row.mtch)
                       alex.kw.results,    # ALEX Kruskal-Wallis test results (output of the KW.hit.express function)
                       q,                  # minimum q value for a gene to be included in output PDF file of box plots
                       gene.annotation)    # gene annotation data
  
{
  # To extract genes below specified q.value for boxplots:
  KW.results=alex.kw.results[!is.na(alex.kw.results$q.KW),]
  selected.genes=KW.results[KW.results$q.KW<q,]
  
  selected.IDS=selected.genes$gene
  selected.lsns=alex.data$alex.lsn[rownames(alex.data$alex.lsn) %in% selected.IDS,]
  selected.lsns=t(selected.lsns)
  selected.expr=alex.data$alex.expr[rownames(alex.data$alex.expr) %in% selected.IDS,]
  selected.expr=t(selected.expr)
  
  # To make sure that both lesion and expression matrices have same patients and genes order
  check.rownames=all(rownames(selected.lsns)==rownames(selected.expr))
  if (check.rownames==FALSE)
    stop("Gene-lesion matrix rownames must match expression matrix rownames.")
  
  check.colnames=all(colnames(selected.lsns)==colnames(selected.expr))
  if (check.colnames==FALSE)
    stop("Gene-lesion matrix patient IDs must match expression matrix patient IDs.")
  
  # add ensembl.ID to gene.name column if empty
  #gene.annotation=gene.annotation %>% mutate_all(na_if,"")
  gene.annotation$gene.name <- ifelse(is.na(gene.annotation$gene.name),gene.annotation$gene, gene.annotation$gene.name)
  
  n=ncol(selected.lsns)
  
  for (i in 1:n)
  {
    ensembl.id=colnames(selected.lsns)[i]
    ensembl.df=as.data.frame(ensembl.id)
    colnames(ensembl.df) <- ("gene")
    genes.df.merged=merge(gene.annotation,ensembl.df,by="gene", all.y=TRUE)
    gene.names=as.character(genes.df.merged$gene.name)
    lsn.box=unlist(selected.lsns[,i])
    expr.box=unlist(selected.expr[,i])
    df=as.data.frame(cbind(lsn.box, expr.box))
    df$expr.box=as.numeric(df$expr.box)
    {
      p=ggplot(df, aes(x = fct_reorder(lsn.box, expr.box, .desc =TRUE), y = expr.box)) +
        geom_boxplot(aes(fill = fct_reorder(lsn.box, expr.box,  .desc =TRUE))) +
        geom_jitter(position=position_jitter(0.02)) + # Increasing the number from 0.02 increase the space between the dots but add some extra dots if number of patients is small
        xlab(paste0(gene.names, '_lsn'))+
        ylab(paste0(gene.names, '_Expression'))+
        theme_bw(base_size = 12) +
        scale_fill_discrete(guide = guide_legend(title = "Lesion")) +
        theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"))+
        theme(legend.text=element_text(size=8))
      print(p)
    }
  }
}


###########################################
# alex.waterfall.prep: Function prepares lesion and expression data of a selected gene for the
# alex.waterfall.plot function

alex.waterfall.prep=function(alex.data,             # output of the alex.prep.lsn.expr function (list of three data table "alex.expr" with expression data, "alex.lsn" with lesion data and row.mtch)
                             alex.kw.results,       # ALEX Kruskal-Wallis test results (output of the KW.hit.express function)
                             gene,                  # gene name or ensembl.ID
                             lsn.data)              # lesion data in a GRIN compatible format
  
  
{
  # Find the row of the KW test with results for the gene
  int.row=which((alex.kw.results$gene==gene)|
                  (alex.kw.results$gene.name==gene))
  
  if (length(int.row)!=1) stop("gene must match exactly ONE row of alex.kw.results")
  
  ens.ID=alex.kw.results$gene[int.row]
  gene.ID=alex.kw.results$gene.name[int.row]
  
  gene.lsn.exp=as.data.frame(colnames(alex.data$alex.lsn))
  colnames(gene.lsn.exp)="ID"
  # Add the lesion data for this gene to the gene.lsn.exp data
  rownames(gene.lsn.exp)=gene.lsn.exp$ID
  gene.lsn=paste0(gene.ID,".lsn")
  lsn.row=which(rownames(alex.data$alex.lsn)==ens.ID)
  gene.lsn.exp[,gene.lsn]=""
  lsn.ids=intersect(gene.lsn.exp$ID,colnames(alex.data$alex.lsn))
  gene.lsn.exp[lsn.ids,gene.lsn]=unlist(alex.data$alex.lsn[lsn.row,lsn.ids])
  
  # Add the expression data for this gene to the gene.lsn.exp data
  gene.rna=paste0(gene.ID,".RNA")
  gene.lsn.exp[,gene.rna]=NA
  rna.ids=intersect(gene.lsn.exp$ID,colnames(alex.data$alex.expr))
  rna.row=which(rownames(alex.data$alex.expr)==ens.ID)
  gene.lsn.exp[rna.ids,gene.rna]=unlist(alex.data$alex.expr[rna.row,rna.ids])
  
  # Find lesions that overlap this gene
  lsn.mtch=(lsn.data$chrom==alex.kw.results$chrom[int.row])&
    (lsn.data$loc.end>=alex.kw.results$loc.start[int.row])&
    (lsn.data$loc.start<=alex.kw.results$loc.end[int.row])
  
  lsns=lsn.data[lsn.mtch,]
  
  res=list(gene.lsn.exp=gene.lsn.exp,       # data table with three columns that has patient ID, type of lesions that affect this gene if any, expression level of the gene
           lsns=lsns,                       # all lesions that affect this gene of interest
           stats=alex.kw.results[int.row,], # KW test results for the gene of interest
           gene.ID=gene.ID)                 # gene name
}

###################################################################
# alex.waterfall.plot: Function return a waterfall plot for expression data by lesion groups of
# a selected gene.

alex.waterfall.plot=function(waterfall.prep,   # Output of the alex.waterfall.prep function
                             lsn.data,         # Lesion data in a GRIN compatible format
                             lsn.clrs=NULL,    # Colors assigned for each lesion group. If NULL, the default.grin.colors function will be used to assign lesion colors automatically
                             ex.clrs=rgb(c(0,0,1,1),
                                         rep(0,4),
                                         c(1,1,0,0),
                                         alpha=c(1,0.5,0.5,1)),
                             delta=0.5)
  
{
  unique.grps=sort(unique(lsn.data$lsn.type))
  
  # assign colors for multiple and none lesion groups (colors to be assigned automatically for other lesion groups)
  if (is.null(lsn.clrs))
  {
    lsn.grps.clr=default.grin.colors(unique.grps)
    common.grps.clr=c(none="gray", multiple="violet")
    lsn.clrs=c(lsn.grps.clr, common.grps.clr)
  }
  
  gene.ID=waterfall.prep$gene.ID
  gene.lsn.exp=waterfall.prep$gene.lsn.exp
  
  lsn.clm=paste0(gene.ID,".lsn")
  rna.clm=paste0(gene.ID,".RNA")
  
  gene.lsn.exp.ord=order(gene.lsn.exp[,lsn.clm],
                         gene.lsn.exp[,rna.clm])
  
  gene.lsn.exp=gene.lsn.exp[gene.lsn.exp.ord,]
  
  pt.num=1:nrow(gene.lsn.exp)
  names(pt.num)=gene.lsn.exp$ID
  
  
  ####################################
  # Set up plotting region
  plot(c(-1.1,+1.5),
       c(0.1,-1.1)*nrow(gene.lsn.exp),
       type="n",axes=F,
       xlab="",ylab="")
  
  ###################################
  # DNA lesion plot
  
  # gene locus
  loc.rng=c(waterfall.prep$stats$loc.start,
            waterfall.prep$stats$loc.end)
  loc.lng=diff(loc.rng)
  
  pos.rng=loc.rng+c(-1,1)*delta*loc.lng
  waterfall.prep$lsns$x.start=(waterfall.prep$lsns$loc.start-pos.rng[1])/(diff(pos.rng))-1.1
  waterfall.prep$lsns$x.end=(waterfall.prep$lsns$loc.end-pos.rng[1])/(diff(pos.rng))-1.1
  
  x.locus=(loc.rng-pos.rng[1])/diff(pos.rng)-1.1
  
  # background for lesion plot
  rect(-1.1,-nrow(gene.lsn.exp),
       -0.1,0,
       col=lsn.clrs["none"],
       border=lsn.clrs["none"])
  
  waterfall.prep$lsns$size=waterfall.prep$lsns$loc.end-
    waterfall.prep$lsns$loc.start+1
  ord=rev(order(waterfall.prep$lsns$size))
  waterfall.prep$lsns=waterfall.prep$lsns[ord,]
  rect(pmax(waterfall.prep$lsns$x.start,-1.1),
       -pt.num[waterfall.prep$lsns$ID],
       pmin(waterfall.prep$lsns$x.end,-0.1),
       -pt.num[waterfall.prep$lsns$ID]+1,
       col=lsn.clrs[waterfall.prep$lsns$lsn.type],
       border=lsn.clrs[waterfall.prep$lsns$lsn.type])
  
  segments(x.locus,-nrow(gene.lsn.exp),
           x.locus,0,col="white",lty=3)
  
  text(x.locus,-1.05*nrow(gene.lsn.exp),
       loc.rng,cex=0.75)
  text(-0.55,+0.05*nrow(gene.lsn.exp),
       paste0(gene.ID," DNA Lesions"))
  
  ##############################################
  # RNA expression plot
  
  rna.lbls=pretty(gene.lsn.exp[,rna.clm])
  ok.lbls=(rna.lbls>min(gene.lsn.exp[,rna.clm],na.rm=T))&
    (rna.lbls<max(gene.lsn.exp[,rna.clm],na.rm=T))
  rna.lbls=rna.lbls[ok.lbls]
  
  gene.lsn.exp$x.rna=(gene.lsn.exp[,rna.clm]-min(gene.lsn.exp[,rna.clm],na.rm=T))/diff(range(gene.lsn.exp[,rna.clm],na.rm=T))
  gene.lsn.exp$x.rna=gene.lsn.exp$x.rna+0.1
  
  x.rna.lbls=(rna.lbls-min(gene.lsn.exp[,rna.clm],na.rm=T))/diff(range(gene.lsn.exp[,rna.clm],na.rm=T))+0.1
  
  n.mdn=function(x)
  {
    mdn=median(x,na.rm=T)
    n=sum(!is.na(x))
    res=unlist(c(n=n,mdn=mdn))
    return(res)
  }
  
  lsn.mdn=aggregate(x=gene.lsn.exp$x.rna,
                    by=list(lsn=gene.lsn.exp[,lsn.clm]),
                    FUN=n.mdn)
  rownames(lsn.mdn$x)=lsn.mdn$lsn
  
  x.mdn=lsn.mdn$x[gene.lsn.exp[,lsn.clm],"mdn"]
  
  text(x.rna.lbls,-1.05*nrow(gene.lsn.exp),
       rna.lbls,cex=0.75)
  segments(x.rna.lbls,0,
           x.rna.lbls,-nrow(gene.lsn.exp),
           lty=3)
  
  rect(pmin(gene.lsn.exp$x.rna,x.mdn),-(1:nrow(gene.lsn.exp)),
       pmax(gene.lsn.exp$x.rna,x.mdn),-(1:nrow(gene.lsn.exp))+1,
       col=lsn.clrs[gene.lsn.exp[,lsn.clm]],
       border=NA)
  
  segments(x.mdn,-(1:nrow(gene.lsn.exp)),
           x.mdn,-(1:nrow(gene.lsn.exp))+1,
           col=lsn.clrs[gene.lsn.exp[,lsn.clm]])
  
  text(0.55,0.05*nrow(gene.lsn.exp),
       paste0(gene.ID," RNA Expression"))
  
  #############################
  # Add legend
  
  lsn.inc=names(lsn.clrs)%in%gene.lsn.exp[,lsn.clm]
  legend(1.15,-0.25*nrow(gene.lsn.exp),
         fill=unlist(lsn.clrs[lsn.inc]),
         legend=names(lsn.clrs[lsn.inc]),
         cex=0.75,border=NA,bty="n")
  
}

#############################################################
# top.alex.waterfall.plots: Function return waterfall plots for top significant genes in the KW results table:

top.alex.waterfall.plots=function(out.dir,            # path to the folder to add waterfall plots
                                  alex.data,          # output of the alex.prep.lsn.expr function (list of three data table "alex.expr" with expression data, "alex.lsn" with lesion data and row.mtch)
                                  alex.kw.results,    # ALEX Kruskal-Wallis test results (output of the KW.hit.express function)
                                  q,                  # minimum q value for a gene to be plotted
                                  lsn.data)           # Lesion data in a GRIN compatible format
  
  
{
  # To extract genes below specified q.value for waterfall plots:
  KW.results=alex.kw.results[!is.na(alex.kw.results$q.KW),]
  selected.genes=KW.results[KW.results$q.KW<q,]
  selected.genes=selected.genes[!is.na(selected.genes$gene.name),]
  selected.genes=selected.genes[!duplicated(selected.genes[ , "gene.name"]),]
  top.genes=selected.genes$gene.name
  
  for (i in 1:length(top.genes))
  {
    pdf(paste0(out.dir,top.genes[i],".pdf"),
        height=8,width=10)
    temp=alex.waterfall.prep(alex.data, alex.kw.results,top.genes[i],lsn.data)
    plot.gene=alex.waterfall.plot(temp, lsn.data)
    dev.off()
  }
}

##################################################################################
# compute distance for lesions

dist.lsn=function(lsn.mtx) # subjects in columns, genes in rows
{
  n=ncol(lsn.mtx)
  res=matrix(NA,n,n)
  colnames(res)=rownames(res)=colnames(lsn.mtx)
  diag(res)=0
  for (i in 1:(n-1))
    for (j in (i+1):n)
    {
      res[i,j]=res[j,i]=sum(lsn.mtx[,i]!=lsn.mtx[,j])
    }
  return(res)
}

##############################################################
# alex.pathway: Function compute the distance between subjects in th dataset based on lesions 
# affecting different genes assigned to the pathway of interest and return two panels of lesion and 
# expression data of ordered subjects based on the computed distances.
# Function also return a matrix of lesion and expression data for each subject ordered based
# on the hiearchial clustering analysis (same order of lesion and expression panels of the figure)

alex.pathway=function(alex.data,          # output of the alex.prep.lsn.expr function (list of three data table "alex.expr" with expression data ready for KW test, "alex.lsn" with overlapped gene lesion data and row.mtch)
                      lsn.data,           # lesion data in a GRIN compatible format
                      pathways,           # data.table with three columns "gene.name" with gene symbols, "ensembl.id" with gene ensembl ID and "pathway" that has the pathway name
                      selected.pathway)   # pathway of interest
  
{
  expr=alex.data$alex.expr
  lsn=alex.data$alex.lsn
  selected.genes=pathways[pathways$pathway==selected.pathway,]
  pathway.ensembl=selected.genes$ensembl.id
  
  path.expr=expr[(rownames(expr)%in%pathway.ensembl),]
  path.lsn=lsn[(rownames(lsn)%in%pathway.ensembl),]
  
  # extract lesion and expression data for genes assigned to the pathway of interest
  path.expr=as.matrix(path.expr)
  path.lsn=as.matrix(path.lsn)
  
  # assign colors to lesion groups
  unique.grps=sort(unique(lsn.data$lsn.type))
  lsn.grps.clr=default.grin.colors(unique.grps)
  common.grps.clr=c(none="gray", multiple="violet")
  lsn.clrs=c(lsn.grps.clr, common.grps.clr)
  
  lsn.grps=names(lsn.clrs)
  clrs=as.character(lsn.clrs)
  path.lsn.clr=path.lsn
  
  # to replace lesion groups with colors in the lsn matrix
  path.lsn.clr[path.lsn.clr %in% lsn.grps] <- clrs[match(path.lsn.clr, lsn.grps, nomatch = 0)]
  path.genes=selected.genes$ensembl.id
  path.gene.names=selected.genes$gene.name
  names(path.gene.names)=path.genes
  
  # to replace ensembl IDs with gene name
  rownames(path.expr)=path.gene.names[rownames(path.expr)]
  rownames(path.lsn)=path.gene.names[rownames(path.lsn)]
  
  # compute the distance between each two genes based on the lesion data using dist.lsn function
  dist.lsn.genes=dist.lsn(t(path.lsn))
  hcl.lsn.genes=hclust(as.dist(dist.lsn.genes),"complete")
  
  # compute distance between subjects based on lesions affecting pathway genes using dist.lsn function
  path.lsn.dist=dist.lsn(path.lsn)
  path.lsn.dist=as.matrix(path.lsn.dist)
  
  path.ls.hcl=hclust(as.dist(path.lsn.dist),method="ward.D2")
  sum.none=rowSums(path.lsn.clr=="gray")
  subj.ord=path.ls.hcl$order
  subj.labels=path.ls.hcl$labels
  
  #####################################
  # side-by-side heatmap
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(c(0,+1.1*ncol(path.lsn.clr)),
       c(0,-(nrow(path.lsn.clr)+3+nrow(path.expr))),
       type="n",xlab="",ylab="",
       axes=F)
  
  for (i in 1:nrow(path.lsn.clr))
    rect(0:(ncol(path.lsn.clr)-1),-(i-1),
         1:ncol(path.lsn.clr),-i,
         col=path.lsn.clr[hcl.lsn.genes$order[i],subj.ord],
         border=NA)
  text(ncol(path.lsn),-(1:nrow(path.lsn))+0.5,
       rownames(path.expr)[hcl.lsn.genes$order],
       pos=4,cex=0.75)
  
  # Add legend
  lsn.inc=names(lsn.clrs)%in%path.lsn
  legend("topright", inset=c(-0.25,0.01),
         fill=unlist(lsn.clrs[lsn.inc]),
         legend=names(lsn.clrs[lsn.inc]),
         cex=0.72,border=NA,bty="n")
  
  for (i in 1:nrow(path.expr))
  {
    y=path.expr[hcl.lsn.genes$order[i],]
    z=(y-mean(y))/sd(y)
    clr=rgb((z>0),0,(z<0),alpha=sqrt(1-exp(-abs(z))))
    rect(0:(ncol(path.lsn.clr)-1),-nrow(path.lsn.clr)-3-(i-1),
         1:ncol(path.lsn.clr),-nrow(path.lsn.clr)-3-i,
         col=clr[subj.ord],
         border=NA)
  }
  
  legend("bottomright", inset=c(-0.25,0.2),
         legend=c("Z.expr<0",0,"Z.expr>0"),
         fill = c("blue", "white", "red"),
         cex=0.72,border=NA,bty="n")
  
  text(ncol(path.expr),-nrow(path.lsn.clr)-3-1:nrow(path.expr)+0.5,
       rownames(path.expr)[hcl.lsn.genes$order],cex=0.75,pos=4)
  
  # Extract ordered lesion and expression data for the pathway genes
  pts.labels=as.data.frame(subj.labels)
  pts.labels<-tibble::rownames_to_column(pts.labels, "subj.ord")
  pts.order=as.data.frame(subj.ord)
  pts.order$index=1:nrow(pts.order)
  ordered.subj.final=merge(pts.labels,pts.order,by="subj.ord", all.y=TRUE)
  ordered.subj.final=ordered.subj.final[order(ordered.subj.final$index),]
  path.lsn.df=as.data.frame(path.lsn)
  ordered.lsn.data=path.lsn.df[ordered.subj.final$subj.labels]
  rownames(ordered.lsn.data) = paste(rownames(ordered.lsn.data),"_lsn")
  path.expr.df=as.data.frame(path.expr)
  ordered.expr.data=path.expr.df[ordered.subj.final$subj.labels]
  rownames(ordered.expr.data) = paste(rownames(ordered.expr.data),"_expr")
  ordered.path.data=rbind(ordered.lsn.data, ordered.expr.data)
  
  return(ordered.path.data)
  
}

##################################################################################
# grin.assoc.lsn.outome: function run association analysis between the binary lesion matrix 
# (output of prep.binary.lsn.mtx function) and treatment outcomes including MRD, EFS and OS

grin.assoc.lsn.outcome=function(lsn.mtx,            # output of prep.binary.lsn.mtx function, each lesion affecting a gene is represented in a separate row and patient will be coded as 1 if he/she is affected with the lesion and 0 otherwise
                                clin.data,          # clinical data table. First column "ID" has patient IDs and each row is a patient clinical data
                                annotation.data,    # gene annotation data file
                                clinvars,           # clinical variables (survival variables such as EFS and OS should be first coded as survival objects using Surv function and added as new columns to the clinical data file, binary variables such as MRD should be coded as 0, 1)
                                covariate=NULL)     # covariates that the model will adjust for if any
{
  lsn.mtx=t(lsn.mtx)
  lsn.df=as.data.frame(lsn.mtx)
  
  clin.data=clin.data
  clin.ids=clin.data$ID
  
  lsn.df=lsn.df[(rownames(lsn.df) %in% clin.ids),]
  lsn.ids=rownames(lsn.df)
  clin.data=clin.data[clin.data$ID %in% lsn.ids,]
  lsn.df=lsn.df[order(rownames(lsn.df)),]
  clin.data=clin.data[order(clin.data$ID),]
  # To make sure that both lesion and clinical data have same patient IDs
  check.rownames=all(rownames(lsn.df)==clin.data$ID)
  if (check.rownames==FALSE)
    stop("Gene-lesion matrix rownames must match patient IDs in the clinical data.")
  
  merged.data=cbind(lsn.df, clin.data)
  
  res<-data.frame(matrix(NA, length(lsn.df),))
  
  for (i in 1:length(clinvars))
  {
    var.name=clinvars[i]
    thisvar<- merged.data[, clinvars[i]]
    
    # If the variable is a survival object, run COXPH models
    if (is.Surv(thisvar))
    {
      message(paste0("Running COXPH models for association with ", var.name, ":",date()))
      surv.time=thisvar[,1]
      surv.censor=thisvar[,2]
      surv.data=cbind(merged.data, surv.time, surv.censor)
      surv.data=surv.data[!is.na(surv.data$surv.time),]
      surv.censor=surv.data$surv.censor
      surv.lsn.clms=surv.data[ ,grepl("ENSG", names(surv.data))]
      
      surv.data.count=cbind(surv.censor, surv.lsn.clms)
      
      surv.event= surv.data.count[surv.data.count$surv.censor==1,]
      surv.without.event= surv.data.count[surv.data.count$surv.censor==0,]
      
      surv.event=surv.event[,-1]
      surv.without.event=surv.without.event[,-1]
      
      # count the number of event free and patients with events who are affected or not affected by genomic lesions
      surv.event.with.lsn=as.data.frame(colSums(surv.event== 1))
      surv.event.without.lsn=as.data.frame(colSums(surv.event== 0))
      surv.no.event.with.lsn=as.data.frame(colSums(surv.without.event== 1))
      surv.no.event.without.lsn=as.data.frame(colSums(surv.without.event== 0))
      
      surv.count=cbind(surv.event.with.lsn, surv.event.without.lsn,
                       surv.no.event.with.lsn, surv.no.event.without.lsn)
      
      colnames(surv.count)=c(paste0(var.name,".event.with.lsn"), paste0(var.name,".event.without.lsn"),
                             paste0(var.name,".no.event.with.lsn"), paste0(var.name,".no.event.without.lsn"))
      
      # Run COXPH models for association between genomic lesions and the survival object without adjustment for any covariate
      if (is.null(covariate)) {
        surv.model<-lapply(surv.lsn.clms,function(x) coxph(Surv(surv.time, surv.censor) ~ x , data=surv.data))
        
        ## To extract coefficients:
        surv.Coeff<-lapply(surv.model,function(f) summary(f)$coefficients[,1])
        ## For univariate mosels, number of columns should be 1
        surv.coeff.mat <- do.call(rbind,lapply(surv.Coeff,matrix,ncol=1,byrow=TRUE))
        
        ## To extract hazard.ratio:
        surv.hazard<-lapply(surv.model,function(f) summary(f)$coefficients[,2])
        surv.hazard.mat = as.matrix(surv.hazard)
        
        ## To extract model error:
        surv.error = lapply(surv.model,function(f) summary(f)$coefficients[,3])
        surv.error.mat <- do.call(rbind,lapply(surv.error,matrix,ncol=1,byrow=TRUE))
        
        ## To calculate lower95% CI:
        surv.lower95 = exp(surv.coeff.mat[,1] - 1.96*surv.error.mat[,1])
        surv.lower.mat = as.matrix(surv.lower95)
        
        ## To calculate upper95% CI:
        surv.upper95 = exp(surv.coeff.mat[,1] + 1.96*surv.error.mat[,1])
        surv.upper.mat = as.matrix(surv.upper95)
        
        ## To extract p-value:
        surv.Pvalue <-lapply(surv.model,function(f) summary(f)$coefficients[,5])
        surv.Pvalue.mat = as.matrix(surv.Pvalue)
        surv.pvalue.mat=do.call(rbind,lapply(surv.Pvalue,matrix,ncol=1,byrow=TRUE))
        pvalue.surv= surv.pvalue.mat[,1]
        
        # Compute FDR with the Pounds & Cheng (2006) estimator of the proportion of tests with a true null (pi.hat
        pi.hat=min(1,2*mean(pvalue.surv,na.rm=T))
        q.surv=pi.hat*p.adjust(pvalue.surv,method="fdr")
        
        ## To prepare final results without covariates adjustment:
        surv.results=cbind(surv.hazard.mat, surv.lower.mat, surv.upper.mat, pvalue.surv, q.surv)
        surv.lsn=as.data.frame(colnames(lsn.mtx))
        colnames(surv.lsn)="Gene_lsn"
        results.surv.final=cbind(surv.lsn, surv.results)
        colnames(results.surv.final)=c("Gene_lsn",paste0(var.name, ".HR") ,
                                       paste0(var.name, ".lower95"), paste0(var.name, ".upper95"), 
                                       paste0(var.name, ".p-value"), paste0(var.name, ".q-value.adj"))
        thisres=cbind(results.surv.final, surv.count)
        
      } else {
        
        # Run COXPH models for association between genomic lesions and survival object while adjusting for some covariates
        covariates=covariate
        covariates=surv.data[,covariates]
        
        surv.model.adj<-lapply(surv.lsn.clms,function(x) coxph(Surv(surv.time, surv.censor) ~ x+covariates , data=surv.data))
        
        ## To extract coefficients:
        surv.Coeff.adj<-lapply(surv.model.adj,function(f) summary(f)$coefficients[,1])
        ## Number of columns depend on how many variables we have in the model
        surv.coeff.mat.adj <- do.call(rbind,lapply(surv.Coeff.adj,matrix,byrow=TRUE))
        num.grps=length(surv.coeff.mat.adj)/ncol(lsn.mtx)
        surv.coeff.mat.adj <- do.call(rbind,lapply(surv.Coeff.adj,matrix,ncol=num.grps,byrow=TRUE))
        
        ## To extract hazard.ratio:
        surv.hazard.adj<-lapply(surv.model.adj,function(f) summary(f)$coefficients[,2])
        surv.hazard.mat.adj <- do.call(rbind,lapply(surv.hazard.adj,matrix,ncol=num.grps,byrow=TRUE))
        
        ## To extract model error:
        surv.error.adj = lapply(surv.model.adj,function(f) summary(f)$coefficients[,3])
        surv.error.mat.adj <- do.call(rbind,lapply(surv.error.adj,matrix,ncol=num.grps,byrow=TRUE))
        
        ## To calculate lower95% CI:
        surv.lower95.adj = exp(surv.coeff.mat.adj[,1] - 1.96*surv.error.mat.adj[,1])
        surv.lower.mat.adj = as.matrix(surv.lower95.adj)
        
        ## To calculate upper95% CI:
        surv.upper95.adj = exp(surv.coeff.mat.adj[,1] + 1.96*surv.error.mat.adj[,1])
        surv.upper.mat.adj = as.matrix(surv.upper95.adj)
        
        ## To extract p-value:
        surv.Pvalue.adj <-lapply(surv.model.adj,function(f) summary(f)$coefficients[,5])
        surv.Pvalue.mat.adj = as.matrix(surv.Pvalue.adj)
        surv.pvalue.mat.adj=do.call(rbind,lapply(surv.Pvalue.adj,matrix,ncol=num.grps,byrow=TRUE))
        pvalue.surv.adj=surv.pvalue.mat.adj[,1]
        
        # Compute FDR with the Pounds & Cheng (2006) estimator of the proportion of tests with a true null (pi.hat
        pi.hat=min(1,2*mean(pvalue.surv.adj,na.rm=T))
        q.surv.adj=pi.hat*p.adjust(pvalue.surv.adj,method="fdr")
        
        ## To prepare final results (covariates adjusted):
        surv.results.adj=cbind(surv.hazard.mat.adj[,1], surv.lower.mat.adj, surv.upper.mat.adj,
                               pvalue.surv.adj, q.surv.adj)
        surv.lsn=as.data.frame(colnames(lsn.mtx))
        colnames(surv.lsn)="Gene_lsn"
        results.surv.final.adj=cbind(surv.lsn, surv.results.adj)
        colnames(results.surv.final.adj)=c("Gene_lsn", paste0(var.name, ".HR.adj"), paste0(var.name, ".lower95.adj"),
                                           paste0(var.name, ".upper95.adj"), paste0(var.name, ".p-value.adj"),
                                           paste0(var.name, ".q-value.adj"))
        thisres=cbind(results.surv.final.adj, surv.count)
      }
    }
    
    # if the variable is numeric such as MRD coded as 0 for MRD negative and 1 for MRD positive patients, run logistic regression models
    else if (is.numeric(thisvar))
    {
      message(paste0("Running logistic regression models for association with ", var.name, ":",date()))
      num.variable=thisvar
      num.data=cbind(merged.data, num.variable)
      num.data=num.data[!is.na(num.data$num.variable),]
      
      num.lsn.clms=num.data[ ,grepl("ENSG", names(num.data))]
      num.variable=num.data$num.variable
      
      num.data.count=cbind(num.variable, num.lsn.clms)
      
      num.event= num.data.count[num.data.count$num.variable==1,]
      num.without.event= num.data.count[num.data.count$num.variable==0,]
      
      num.event=num.event[,-1]
      num.without.event=num.without.event[,-1]
      
      # count the number of event free and patients with events who are affected or not affected by genomic lesions
      num.event.with.lsn=as.data.frame(colSums(num.event== 1))
      num.event.without.lsn=as.data.frame(colSums(num.event== 0))
      num.no.event.with.lsn=as.data.frame(colSums(num.without.event== 1))
      num.no.event.without.lsn=as.data.frame(colSums(num.without.event== 0))
      
      num.count=cbind(num.event.with.lsn, num.event.without.lsn,
                      num.no.event.with.lsn, num.no.event.without.lsn)
      
      colnames(num.count)=c(paste0(var.name,".event.with.lsn"), paste0(var.name,".event.without.lsn"),
                            paste0(var.name,".no.event.with.lsn"), paste0(var.name,".no.event.without.lsn"))
      
      
      # Run logistic regression models for association between genomic lesions and binary variables without adjustment for any covariate
      if (is.null(covariate)) {
        log.fit<-lapply(num.lsn.clms,function(x) glm(num.variable  ~ x ,data = num.data, family = "binomial"))
        # To extract model coefficients
        coeff = lapply(log.fit,function(f) summary(f)$coefficients[,1])
        
        ## for univariate models, number of columns "ncol" should be 2 (column for intercept and column for model coefficient):
        coeff.mat <- do.call(rbind,lapply(coeff,matrix,ncol=2,byrow=TRUE))
        mod_coeff = lapply(coeff, head, 2)
        
        ## To extract model error:
        error = lapply(log.fit,function(f) summary(f)$coefficients[,2])
        error.mat <- do.call(rbind,lapply(error,matrix,ncol=2,byrow=TRUE))
        
        ## To calculate odds ratio and 95% CI:
        ## N.B) Here, we are using the standard errors not the robust SE
        odds.ratio= exp(coeff.mat[,2])
        odds.mat=as.matrix(odds.ratio)
        lower.confint = exp(coeff.mat[,2] - 1.96*error.mat[,2])
        lower.mat=as.matrix(lower.confint)
        higher.confint = exp(coeff.mat[,2] + 1.96*error.mat[,2])
        higher.mat=as.matrix(higher.confint)
        
        ## To extract P-value:
        Pvalue<-lapply(log.fit,function(f) summary(f)$coefficients[,4])
        pvalue.mat=do.call(rbind,lapply(Pvalue,matrix,ncol=2,byrow=TRUE))
        pvalue.num= pvalue.mat[,2]
        
        # Compute FDR with the Pounds & Cheng (2006) estimator of the proportion of tests with a true null (pi.hat
        pi.hat=min(1,2*mean(pvalue.num,na.rm=T))
        q.num=pi.hat*p.adjust(pvalue.num,method="fdr")
        
        # prepare final results without adjustment for any covariate
        results.logistic=cbind(odds.mat, lower.mat, higher.mat, pvalue.num, q.num)
        logistic.lsn=as.data.frame(colnames(lsn.mtx))
        colnames(logistic.lsn)="Gene_lsn"
        results.logistic.final=cbind(logistic.lsn, results.logistic)
        colnames(results.logistic.final)=c("Gene_lsn", paste0(var.name, ".odds.ratio"),
                                           paste0(var.name, ".lower95"), paste0(var.name, ".upper95"),
                                           paste0(var.name, ".p-value"), paste0(var.name, ".q-value"))
        num.results.count.added=cbind(results.logistic.final, num.count)
        thisres=num.results.count.added
        
      } else {
        
        # Run logistic regression models for association between genomic lesions and binary variable with covariate adjustment
        covariates=covariate
        covariates=num.data[,covariates]
        log.fit.adj<-lapply(num.lsn.clms,function(x) glm(data = num.data, num.variable  ~ x + covariates , family = "binomial"))
        ## To extract coefficients
        coeff.adj = lapply(log.fit.adj,function(f) summary(f)$coefficients[,1])
        coeff.mat.adj <- do.call(rbind,lapply(coeff.adj,matrix,byrow=TRUE))
        num.grps=round(length(coeff.mat.adj)/ncol(lsn.mtx))
        coeff.mat.adj <- do.call(rbind,lapply(coeff.adj,matrix,ncol=num.grps,byrow=TRUE))
        
        ## To extract model error:
        error.adj = lapply(log.fit.adj,function(f) summary(f)$coefficients[,2])
        error.mat.adj <- do.call(rbind,lapply(error.adj,matrix,ncol=num.grps,byrow=TRUE))
        
        ## To calculate odds ratio and 95% CI:
        ## N.B) Here, we are using the standard errors not the robust SE
        odds.ratio.adj= exp(coeff.mat.adj[,2])
        odds.mat.adj=as.matrix(odds.ratio.adj)
        lower.confint.adj = exp(coeff.mat.adj[,2] - 1.96*error.mat.adj[,2])
        lower.mat.adj=as.matrix(lower.confint.adj)
        higher.confint.adj = exp(coeff.mat.adj[,2] + 1.96*error.mat.adj[,2])
        higher.mat.adj=as.matrix(higher.confint.adj)
        
        ## To extract P-value:
        Pvalue.adj<-lapply(log.fit.adj,function(f) summary(f)$coefficients[,4])
        pvalue.mat.adj=do.call(rbind,lapply(Pvalue.adj,matrix,ncol=num.grps,byrow=TRUE))
        pvalue.num.adj= pvalue.mat.adj[,2]
        
        # Compute FDR with the Pounds & Cheng (2006) estimator of the proportion of tests with a true null (pi.hat
        pi.hat=min(1,2*mean(pvalue.num.adj,na.rm=T))
        q.num.adj=pi.hat*p.adjust(pvalue.num.adj,method="fdr")
        
        # prepare final results while adjsuting for covariates
        results.logistic.adj=cbind(odds.mat.adj, lower.mat.adj, higher.mat.adj, pvalue.num.adj, q.num.adj)
        logistic.lsn=as.data.frame(colnames(lsn.mtx))
        colnames(logistic.lsn)="Gene_lsn"
        results.logistic.final.adj=cbind(logistic.lsn, results.logistic.adj)
        colnames(results.logistic.final.adj)=c("Gene_lsn", paste0(var.name, ".odds.ratio.adj"),
                                               paste0(var.name, ".lower95.adj"), paste0(var.name, ".upper95.adj"),
                                               paste0(var.name, ".p-value.adj"), paste0(var.name, ".q-value.adj"))
        num.results.count.added.adj=cbind(results.logistic.final.adj, num.count)
        thisres=num.results.count.added.adj
        
      }  
    } 
    res<-cbind.data.frame(res, thisres)
    res$matrix.NA..length.lsn.df....<-NULL
    
    gene.list=str_split_fixed(res$Gene_lsn, "_", 2)
    colnames(gene.list)=c("gene","lsn")
    annotation.merged=merge(annotation.data,gene.list,by="gene", all.y=TRUE)
    res.final=cbind(annotation.merged, res)
    res.final$lsn<-NULL
  }
  
  return(res.final)
}


########################################################
# Lesion boundaries evaluation
###################################################

# grin.lsn.boundaries: The function evaluates Copy number variations that include gain and deletions
# as boundaries and return a table of ordered lesion boundaries by start and end positions on each
# chromosome.

grin.lsn.boundaries=function(lsn.data, # Lesion data file that should be limited to include gain OR deletions (If gains are splitted to gain and amplification based on the log2Ratio value of the CNV segmentation file, the two categories can be included, same for homozygous and heterozygous deletions)
                             chrom.size) # Chromosize size table
{
  lsn.data=lsn.data
  chr.size=chrom.size
  chroms=chr.size$chrom
  lsn.bound = data.frame()
  {
    for (i in chroms) {
      chr=i
      chr.lsns=lsn.data[lsn.data$chrom==chr,]
      ord=order(chr.lsns$loc.start,
                chr.lsns$loc.end)
      
      # specify all unique start and end positions for included CNVs
      uniq.start=unique(chr.lsns$loc.start)
      uniq.end=unique(chr.lsns$loc.end)
      chrom.size=chr.size$size[chr.size$chrom==chr]
      
      # boundary will be the region between each unique start and end positions 
      # large size lesions will be splitted into multiple boundaries based on other smaller size lesions that affect the same region if any
      all.start=unique(c(1,uniq.start,uniq.end+1))
      all.end=unique(c(uniq.start-1,uniq.end,chrom.size))
      all.end=all.end[!all.end ==0]
      all.start=sort(all.start)
      all.end=sort(all.end)
      chr.lsn.loci=cbind.data.frame(gene=paste0("chr",chr,"_",all.start,"_",all.end),
                                    chrom=chr,
                                    loc.start=all.start,
                                    loc.end=all.end)
      chr.lsn.bound = chr.lsn.loci
      lsn.bound=rbind(lsn.bound, chr.lsn.loci)
      
    }
    lsn.bound$diff=(lsn.bound$loc.end - lsn.bound$loc.start)
    lsn.bound=lsn.bound[!lsn.bound$diff<=0,]
    
    return(lsn.bound)
    
  }
}

############################################################################################
# genomewide.log10q.plot: The function return a genome-wide lesion plot based on -log(10) q-value
# of each of the evaluated lesion boundaries. The plot is lesion type specific.

genomewide.log10q.plot=function(grin.res,        # GRIN results (output of the grin.stats function)
                                lsn.grps,        # selected lesion groups to be added to the plot
                                lsn.colors=NULL, # Lesion colors
                                max.log10q=NULL) # Maximum log10 q value to be added to the plot
  
{
  if (!is.element("x.start",colnames(grin.res$lsn.data)))
    grin.res=compute.gw.coordinates(grin.res)
  
  grin.res$lsn.data$x.ID=as.numeric(as.factor(grin.res$lsn.data$ID))
  n=max(grin.res$lsn.data$x.ID)
  n.chr=nrow(grin.res$chr.size)
  
  # set up plotting region
  plot(c(-0.05,1.2)*n,c(0,-1.1*grin.res$chr.size$x.end[n.chr]),
       type="n",axes=F,xlab="",ylab="")
  
  # background colors for chromosomes
  rect(0*n,-grin.res$chr.size$x.start,
       1*n,-grin.res$chr.size$x.end,
       col=c("lightgray","gray"),
       border=NA)
  
  lsn.types=lsn.grps
  selected.clms = list()
  
  for (i in 1:length(lsn.types)) {
    temp=grin.res$gene.hits[grepl(lsn.types[i],colnames(grin.res$gene.hits))]
    selected.clms[[i]] <- temp
  }
  
  final.data = do.call(cbind, selected.clms)
  
  if (is.null(lsn.colors))
  {
    lsn.colors=default.grin.colors(lsn.types)
  }
  grin.res$lsn.data$lsn.colors=lsn.colors[grin.res$lsn.data$lsn.type]
  
  grin.res$lsn.data$lsn.size=grin.res$lsn.data$x.end-grin.res$lsn.data$x.start
  ord=order(grin.res$lsn.data$lsn.size,decreasing=T)
  
  nsubj.mtx=unlist(grin.res$gene.hits[,paste0("nsubj.",lsn.types)])
  qval.mtx=unlist(grin.res$gene.hits[,paste0("q.nsubj.",lsn.types)])
  nsubj.data=cbind.data.frame(gene=grin.res$gene.hits$gene,
                              x.start=grin.res$gene.hits$x.start,
                              x.end=grin.res$gene.hits$x.end,
                              nsubj=nsubj.mtx,
                              log10q=-log10(qval.mtx),
                              lsn.type=rep(lsn.types,each=nrow(grin.res$gene.hits)))
  nsubj.data=nsubj.data[nsubj.data$nsubj>0,]
  nsubj.data$lsn.colors=lsn.colors[nsubj.data$lsn.type]
  
  ord=order(nsubj.data$nsubj,decreasing=T)
  nsubj.data=nsubj.data[ord,]
  
  nsubj.data$log10q[nsubj.data$log10q>max.log10q]=max.log10q
  
  ord=order(nsubj.data$log10q,decreasing=T)
  nsubj.data=nsubj.data[ord,]
  
  segments(0*n,
           -(nsubj.data$x.start+nsubj.data$x.end)/2,
           0*n+1*nsubj.data$log10q/max(nsubj.data$log10q)*n,
           col=nsubj.data$lsn.colors)
  
  text(c(n,0)[(1:n.chr)%%2+1],
       pos=c(4,2)[(1:n.chr)%%2+1],
       -(grin.res$chr.size$x.start+grin.res$chr.size$x.end)/2,
       grin.res$chr.size$chrom,
       cex=0.75)
  
  legend(n/2,-1.05*grin.res$chr.size$x.end[n.chr],
         fill=lsn.colors,cex=0.9,
         legend=names(lsn.colors),
         xjust=0.5,
         ncol=length(lsn.colors),
         border=NA,bty="n")
  
  text(1*n/2,0,
       "-log10(q)",cex=0.95,
       pos=3)
  
  text(c(0,0.25, 0.5, 0.75, 1)*n,
       -grin.res$chr.size$x.end[n.chr],
       c(0,max.log10q/4, max.log10q/2, round(max.log10q/1.3333333333, 1), max.log10q),
       cex=0.75,pos=1)
}


###########################################################
# prep.grin.gsda: Function prepare a modified binary lesion matrix where all patients without
# any type of lesions "none" are coded as 0 and patients with any or multiple types of lesions
# are coded as 1. Function also return a selected group og gene-sets retrieved from the 
# Molecular Signature Database (MSigDB) that has a certain gene size range based on specified
# small.gset and big.gset arguments. Function will also prepare the expression data file for 
# the grin.gsda function if provided by the user.

prep.grin.gsda=function(lsn.bin.mtx,       # lesion binary matrix (output of the prep.lsn.type.matrix function) 
                        lsn.data,          # lesion data
                        expr.data=NULL,    # expression data, first column "gene"
                        big.gset=NULL,     # largest number of genes per gene-set
                        small.gset=NULL)   # smallest number of genes per gene-set
  
{
  lsn.mtx=lsn.bin.mtx
  
  # generate the binary lesion matrix by coding all patients without any type of lesions "none" and patients with any or multiple types of lesions as 1  
  lsn.mtx[lsn.mtx=="none"]<-0       
  lsn.types=unique(lsn.data$lsn.type)
  lsn.mtx[lsn.mtx%in%lsn.types]<-1
  lsn.mtx[lsn.mtx=="multiple"]<-1
  lsn.df=as.data.frame(lsn.mtx)
  
  omic.clms=ncol(lsn.df)
  lsn.df <- lsn.df %>% mutate_at(1:omic.clms, as.numeric)
  lsn.mtx.final=data.matrix(lsn.df)
  
  # obtain gene-set data
  message(paste0("Obtaining gene-set data from Molecular Signature Database (MSigDB): ",date()))
  msigdb.data=msigdbr(species = "Homo sapiens")
  msigdb.data=as.data.frame(msigdb.data)
  
  #######################################
  # Prepare GSDA vset data
  
  vset.data=cbind.data.frame(vID=msigdb.data$ensembl_gene,
                             Vgene.name=msigdb.data$human_gene_symbol,
                             vset=msigdb.data$gs_name)
  
  
  # Keep only gene-sets with reasonable number of genes that range between numbers specified by the user
  vset.tbl=table(vset.data$vset)
  big.vset=which(vset.tbl>big.gset)
  small.vset=which(vset.tbl<small.gset)
  drop.vset=(vset.data$vset%in%(names(vset.tbl)[c(big.vset,small.vset)]))
  vset.data=vset.data[!drop.vset,]
  ord=order(vset.data$vset)
  vset.data=vset.data[ord,]
  
  if (!is.null(expr.data)) {
    
    # prepare expression data for the grin.gsda function
    rownames(expr.data) <- expr.data[,1]
    expr.data[,1] <- NULL
    expr.data=data.matrix(expr.data)
    
    res=list(grin.gsda.omic=lsn.mtx.final,
             grin.gsda.expr=expr.data,
             grin.gsda.vset=vset.data)
  }
  else{
    
    res=list(grin.gsda.omic=lsn.mtx.final,
             grin.gsda.vset=vset.data)
  } 
  return(res)
}


###########################################################
# grin.gsda: Function run gene-set distance analysis (GSDA) for association between lesions
# that affect group of genes that are part of gene-sets downloaded from Molecular 
# Signature Database (MSigDB) or provided by the user and clinical variables of interest. 
# Function can also run GSDA analysis using gene expression data if provided by the user.

grin.gsda=function(grin.gsda.data,    # list of three data.frames with lesion binary matrix, expression data and gene-sets data of the size specified by the user in the prep.grin.gsda function (output of the prep.grin.gsda function)
                   clin.data,         # data.frame of clinical data matrix with subjects as rows and variables as columns, must include a column named "ID" with identifiers matching column names of omic.data
                   gene.annotation,   # gene annotation data file
                   lesion= FALSE,     # should be specified as TRUE if the lesion binary matrix will be evaluated in the GSDA analysis
                   expression=FALSE,  # should be specified as TRUE if the expression data will be evaluated in the GSDA analysis
                   clin.vars,         # vector name or numeric index of endpoint/treatment variable(s) in clin.data
                   path=NULL,         # pathways provided by the user (other than the GSDA gene sets directly retrieved using msigdb package). data.table should has three columns "vID" with ensembl IDs, "Vgene.name" with gene names, and "vset" with the pathway name
                   omic.dist,         # distance metric for omic data, may be "oe" (overall Euclidean), "me" (marginal Euclidean), "om" (overall Manhattan), or "mm" (marginal Manhattan)
                   clin.dist,         # distance metric for clinical data, may be "oe", "me", "om", or "mm" as defined above, or "ct" (categorical) or "st" (survival time)
                   mess.freq=200)     # print a progress message every mess.freq gene-sets)  
  
{
  if (lesion==TRUE) {
    # prepare lesion data for analysis
    t0=proc.time()
    
    omic.data=grin.gsda.data$grin.gsda.omic
    if (is.data.frame(path)) {
      vset.data=path
    } else {
      vset.data=grin.gsda.data$grin.gsda.vset
    }
    clin.data=clin.data
    ensembl.annotation=gene.annotation
    
    ensembl.annotation=cbind.data.frame(ensembl.ID=hg19.gene.annotation$gene,
                                        Gene=hg19.gene.annotation$gene.name)
    
    
    gsda.data=prep.gsda(omic.data,
                        clin.data,
                        vset.data)
    
    omic.data=gsda.data$omic.data
    clin.data=gsda.data$clin.data
    vset.data=gsda.data$vset.data
    vset.index=gsda.data$vset.index
    t1=proc.time()
    prep.time=(t1-t0)[3]
    comp.time=(t1-t1)[3]
    
    n=nrow(clin.data)          # number of subjects
    nset=nrow(vset.index)      # number of gene-sets
    p.vset=rep(0,nset)         # initialize p-value vector
    vset.list=rep("NA",nset)   # initialize vector of variable set lists 
    dCor.stats=rep(0,nset)    # initialize vector of marginal distance correlation statistics
    comp.time=rep(NA,nset)     # initialize vecto of compute time for each variable set
    
    # Loop over gene-sets
    for (i in 1:nset)
    {
      if (i%%mess.freq==1)
        message(paste("Computing distance correlation for gene-set ",
                      i," of ",nset,": ",date()))
      # Get variables in this variable set
      t0=proc.time()
      vset.ind=(vset.index$start.index[i]:vset.index$end.index[i]) # rows with variables in this variable set
      vset.vIDs=vset.data$vID[vset.ind]                            # ids of variables in this variable set
      genes.df=as.data.frame(vset.vIDs)
      colnames(genes.df) <- ("ensembl.ID")
      genes.df.merged=merge(ensembl.annotation,genes.df,by="ensembl.ID", all.y=TRUE)
      gene.names=as.character(genes.df.merged$Gene)
      vset.list[i]=paste(gene.names,collapse=" | ")
      
      # Extract data matrix for this variable set
      omic.data=as.data.frame(omic.data)
      vset.mtx=omic.data[vset.vIDs, ]
      vset.mtx=tibble::rownames_to_column(vset.mtx, "Gene")
      vset.mtx=subset(vset.mtx, grepl('^ENSG', Gene))
      rownames(vset.mtx) <- vset.mtx[,1]
      vset.mtx[,1] <- NULL
      vset.mtx=data.matrix(vset.mtx)
      
      
      # Compute distances for this variable set
      vset.mtx=t(vset.mtx)                                         # transpose the matrix for use in dist.corr
      t1=proc.time()                                               # check clock to calculate computing times
      prep.time=prep.time+(t1-t0)[3]                               # calculate data prep time for this variable set
      
      if(ncol(vset.mtx) > 0)
      {     
        temp.res=dist.corr(vset.mtx,                                 # compute the distance correlation stats for this variable set
                           clin.data[,clin.vars],
                           omic.dist,clin.dist)
        
        p.vset[i]=temp.res$p.odCor                                   # get the overall distance correlation p-value
        dCor.stats[i]=temp.res$odCor                               # get the overall distance correlation statistic
      }
      if(ncol(vset.mtx) == 0)
      {
        p.vset[i]=1
        dCor.stats[i]=0        
      }
      t2=proc.time()                                               # get time for calculating compute time
      comp.time[i]=(t2-t1)[3]     
    } # end loop over gene-sets
    
    ################
    # result object: a data.frame with one row per gene-set and the following colums:
    res.lsn=cbind.data.frame(vset=vset.index$value,  # name of variable set (gene-set)
                             vIDs=vset.list,         # character string with list of variables in the variable set (gene-set)
                             dCor=dCor.stats,      # distance statistic for the variable set
                             p.vset=p.vset,          # p-value
                             comp.time=comp.time)    # computing time
    res=res.lsn
  }
  
  if (expression==TRUE) {
    # prepare data for analysis
    t0=proc.time()
    
    omic.data=grin.gsda.data$grin.gsda.expr
    
    if (is.data.frame(path)) {
      vset.data=path
    } else {
      vset.data=grin.gsda.data$grin.gsda.vset
    }
    
    clin.data=clin.data
    ensembl.annotation=gene.annotation
    ensembl.annotation=cbind.data.frame(ensembl.ID=hg19.gene.annotation$gene,
                                        Gene=hg19.gene.annotation$gene.name)
    
    gsda.data=prep.gsda(omic.data,
                        clin.data,
                        vset.data)
    
    omic.data=gsda.data$omic.data
    clin.data=gsda.data$clin.data
    vset.data=gsda.data$vset.data
    vset.index=gsda.data$vset.index
    t1=proc.time()
    prep.time=(t1-t0)[3]
    comp.time=(t1-t1)[3]
    
    n=nrow(clin.data)          # number of subjects
    nset=nrow(vset.index)      # number of gene-sets
    p.vset=rep(0,nset)         # initialize p-value vector
    vset.list=rep("NA",nset)   # initialize vector of variable set lists 
    dCor.stats=rep(0,nset)    # initialize vector of marginal distance correlation statistics
    comp.time=rep(NA,nset)     # initialize vecto of compute time for each variable set
    
    # Loop over gene-sets
    for (i in 1:nset)
    {
      if (i%%mess.freq==1)
        message(paste("Computing distance correlation for gene-set ",
                      i," of ",nset,": ",date()))
      
      # Get variables in this variable set
      t0=proc.time()
      vset.ind=(vset.index$start.index[i]:vset.index$end.index[i]) # rows with variables in this variable set
      vset.vIDs=vset.data$vID[vset.ind]                            # ids of variables in this variable set
      genes.df=as.data.frame(vset.vIDs)
      colnames(genes.df) <- ("ensembl.ID")
      genes.df.merged=merge(ensembl.annotation,genes.df,by="ensembl.ID", all.y=TRUE)
      gene.names=as.character(genes.df.merged$Gene)
      vset.list[i]=paste(gene.names,collapse=" | ")
      
      # Extract data matrix for this variable set
      omic.data=as.data.frame(omic.data)
      vset.mtx=omic.data[vset.vIDs, ]
      vset.mtx=tibble::rownames_to_column(vset.mtx, "Gene")
      vset.mtx=subset(vset.mtx, grepl('^ENSG', Gene))
      rownames(vset.mtx) <- vset.mtx[,1]
      vset.mtx[,1] <- NULL
      vset.mtx=data.matrix(vset.mtx)
      
      # Compute distances for this variable set
      vset.mtx=t(vset.mtx)                                         # transpose the matrix for use in dist.corr
      t1=proc.time()                                               # check clock to calculate computing times
      prep.time=prep.time+(t1-t0)[3]                               # calculate data prep time for this variable set
      
      if(ncol(vset.mtx) > 0)
      {
        temp.res=dist.corr(vset.mtx,                                 # compute the distance correlation stats for this variable set
                           clin.data[,clin.vars],
                           omic.dist,clin.dist)
        
        p.vset[i]=temp.res$p.odCor                                   # get the overall distance correlation p-value
        dCor.stats[i]=temp.res$odCor                                # get the overall distance correlation statistic
      }
      if(ncol(vset.mtx) == 0)
      {
        p.vset[i]=1
        dCor.stats[i]=0
      }
      t2=proc.time()                                         # get time for calculating compute time
      comp.time[i]=(t2-t1)[3]     
    } # end loop over gene-sets
    
    ################
    # result object: a data.frame with one row per gene-set and the following colums:
    res.expr=cbind.data.frame(vset=vset.index$value,  # name of variable set (gene-set)
                              vIDs=vset.list,         # character string with list of variables in the variable set (gene-set)
                              dCor=dCor.stats,      # distance statistic for the variable set
                              p.vset=p.vset,          # p-value
                              comp.time=comp.time)    # computing time
    
    res=res.expr
  }
  
  if (expression==TRUE && lesion==TRUE) {
    res=list(gsda.lesion.results=res.lsn,
             gsda.expression.results=res.expr)
  }
  
  return(res)
}
