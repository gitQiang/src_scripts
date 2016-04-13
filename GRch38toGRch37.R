## mapping use http://genome.ucsc.edu/cgi-bin/hgLiftOver

### get the mapped variants from Clinvar itself
Clinvar <- read.delim("variant_summary.txt")
gr37 <- Clinvar[Clinvar[,"Assembly"]=="GRCh37",]
genedx <- read.delim("clinvar_result.txt")
remap <- genedx[!(genedx[,"Name"] %in% gr37[,"Name"]),]
qwt(remap,file="toGRch37_482.txt",flag=2)

## try to mapping the rest 5 using liftover
beda <- remap[,c("Chromosome","Location")]
bedb <- sapply(1:dim(beda)[1],function(i) {
         if(any(grepl("-",beda[i,2]))){ unlist(strsplit(beda[i,2]," - "))        ;
         }else{ rep(beda[i,2],2); }
})
bedb <- t(bedb)
c <- paste("chr",beda[,1],":",bedb[,1],"-",bedb[,2],sep="")
subs <- c(c=="chr:-" | c=="chrChromosome:Location-Location")
c <- c[!subs]
qwt(c,file="GeneDx_5.bed")

### get ref and alt from variant names to run wannovar
genedx37 <- gr37[gr37[,"Name"] %in% genedx[,"Name"], ]
wanno <- genedx37[,c("Chromosome","Start","Stop","ReferenceAllele","AlternateAllele")]
qwt(wanno,file="Genedx37_23609.txt")

anno5 <- read.table("GeneDx_5_gr37.txt")
colnames(anno5) <- colnames(wanno)
wannoa <- rbind(wanno,anno5)
qwt(wannoa,file="Genedx_GR37.txt")


### get ref and alt for 975 P or VLP variants
Pvlp <- read.delim("Table3Sup.txt")
tmp <- Pvlp[,4]
tmp <- gsub("+","",tmp)
tmp <- gsub("c.","",tmp)
tmp <- gsub("_","",tmp)
tmp <- gsub("-","",tmp)

refalt <- sapply(1:length(tmp), function(i){
    if(grepl("del",tmp[i])){
        c(gsub("del","",gsub("\\d","",tmp[i])),"-")
    }else if(grepl("ins",tmp[i])){
        c("-",gsub("ins","",gsub("\\d","",tmp[i])))
    }else if(grepl("dup",tmp[i])){
        onea <- gsub("dup","",gsub("\\d","",tmp[i]))
        c(onea,paste(onea,onea,sep=""))
    }else{
        unlist(strsplit(gsub("\\d","",tmp[i]),">"))
    }
    })
refalt <- t(refalt)
#spchar <- substr(refalt[1,1],1,1)
#refalt[,1] <- gsub(spchar,"",refalt[,1])

wann  <- cbind(Pvlp[,2],Pvlp[,3],Pvlp[,3],refalt)
qwt(wann,file="GeneDx_PVLP.txt")

## fixed version by manual double check

tmp <- read.delim("GeneDx_PVLP_fixed.txt",header=FALSE)
tmp[tmp[,4]=="",4] <- "-"
tmp[tmp[,5]=="",5] <- "-"
qwt(tmp,file="GeneDx_PVLP_fixed_alt.txt")
tmp[tmp=="-"] <- "*"
qwt(tmp,file="GeneDx_PVLP_fixed_alt1.txt")
tmp[tmp=="-"] <- "."
qwt(tmp,file="GeneDx_PVLP_fixed_alt2.txt")

