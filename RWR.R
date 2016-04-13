# RWR.R: the RWR algorithm
# InterfaceTest: the interface for this algorithm
# Descriptions:
# The random walk with restart algorithm on protein-protein interaction (PPI) network with a mutated gene list and candidate gene list
# geneList: a file name of candidate gene list, human gene symbol
# mutatedGenes: a file name of mutated gene list, human gene symbol
# network: the PPI network used, STRING (STRING database) and FI (Functional interaction network, Guanming, Xin Feng and Lincoln Stein, A human functional protein interaction network and its application to cancer data analysis, Genome Biology, 2010)
# outfile: subname in output files
InterfaceTest <- function(geneList=NULL,mutatedGenes,network="STRING",outfile){
	
	
	if(network=="STRING"){
		file <- "networks/STRING.interaction600.txt"
		W <- STRINGnetwork(file)
		}else if(network=="FI") {
			file <- "networks/FIs_043009.txt"	
			W  <- FInetwork(file)
			}else if(network=="PrePPI"){
				file <- "networks/PrePPI.txt"	
				W  <- FInetwork(file)
			}else{print("The network is not support yet!")}
	
	
	# normalize the matrix	
	nodes <- rownames(W)
	n.nodes <- length(nodes)
	normal.vector <- colSums(W)

	normal.vector[normal.vector==0] <- 1
	normal.matrix <- matrix(normal.vector,n.nodes,n.nodes,byrow=TRUE)
	W <- W/normal.matrix
	
	# read the ID mapping table
	idmapfile <- "networks/ID mapping.txt"
	maptable <- as.matrix(read.delim(idmapfile,sep=","))
	
	# get the gene sets from genelist and mutatedGenes
	if(!is.null(geneList)){
		geneInterest <- as.matrix(read.table(geneList))
	}
	mutatedgene <- as.matrix(read.table(mutatedGenes))
	
	
	netMap <- matrix(0,n.nodes,2)
	netMap[,1] <- nodes
	
	# id mapping for genelist and mutatedGenes: Hmuan Symbol
	if(network=="STRING"){
		# mapping to ensemble protein ID
		if(!is.null(geneList)){
			geneInterestMap <- geneListMap(geneInterest,maptable,C=c(3,2))
		}
		mutatedgeneMap <- geneListMap(mutatedgene,maptable,C=c(3,2))
		
		tmp <- match(nodes,maptable[,2],nomatch=-1)
		tmpsubs <- tmp>0
		netMap[tmpsubs,2] <- maptable[tmp[tmpsubs],3]
		}else if(network=="FI" | network=="PrePPI") {
		# mapping to Uniprot ID
		if(!is.null(geneList)){
			geneInterestMap <- geneListMap(geneInterest,maptable,C=c(3,4))
		}
		mutatedgeneMap <- geneListMap(mutatedgene,maptable,C=c(3,4))
		
		tmp <- match(nodes,maptable[,4],nomatch=-1)
		tmpsubs <- tmp>0
		netMap[tmpsubs,2] <- maptable[tmp[tmpsubs],3]	
			}
			
	rm(tmp)


	# run the RWR algorithm
	subs <- nodes %in% mutatedgeneMap[,2]
	p0 <- 1:n.nodes
	p0[subs] <- 1/sum(subs)
	p0[!subs] <- 0
	
	d <- p0
	alpha <- 0.9
	p1 <- RWR(W,p0,alpha,d)
	
	
	X <- sort(p1,decreasing=TRUE,index.return=TRUE)
	XTable <- cbind(netMap[X$ix,c(2,1)],p1[X$ix],1:n.nodes)
	write.table(XTable,file=paste(network,outfile,"all.txt",sep="_"),quote=FALSE,row.names=FALSE,col.names=FALSE)
	
	if(!is.null(geneList)){
		subs <- XTable[,2] %in% geneInterestMap[,2]
		resultTable <-  XTable[subs,]
		write.table(resultTable,file=paste(network,outfile,"based.txt",sep="_"),quote=FALSE,row.names=FALSE,col.names=FALSE)
	}
	

}


RWR <- function(W,p0,alpha,d){
	# RWR: random walk with restart 
	
	# parameters:
	# W: the edge weights matrix
	# p0 : the start point for iteration
	# alpha: the restart probability
	# d: starting vector

	p1 <- (1-alpha)*W%*%p0 + alpha*d;
	iter <- 1
	while (sum(abs(p0-p1))>10^(-10)){
    		p0=p1;
    		p1 <- (1-alpha)*W%*%p0 + alpha*d;
			iter <- iter + 1
	}
	
	print(iter)
	p1

}

STRINGnetwork <- function(file){
	
	#file <- "STRING.protein.links.v9.1.txt"
	#links.text <- as.matrix(read.table(file, fill=T, as.is=T, col.names=1:max(count.fields(file)),skip=1))
	#links.text <- links.text[as.numeric(links.text[,3])>=600,]
	
	#file <- "networks/STRING.interaction600.txt"
	links.text <- as.matrix(read.table(file))
	
	net.text <- cbind(substr(links.text[,1],6,20),substr(links.text[,2],6,20),links.text[,3])
	
	net <- read_net(net.text)
	
	W <- net$matrix
	W <- W/max(W)
	
	W
	
}

FInetwork <- function(file){
	
	#file <- "networks/FIs_043009.txt"	
	net.text <- as.matrix(read.table(file))
	net.text <-  cbind(net.text,1)
	net.text <-  rbind(net.text,cbind(net.text[,2],net.text[,1],net.text[,3]))

	net <- read_net(net.text)
	
	W <- net$matrix
	
	W
}

read_net <- function(net.text){
	
	net.node <- unique(union(net.text[,1],net.text[,2]))
	net.node <- net.node[net.node != ""]
	net.size <- length(net.node)
	net.edge <- cbind(as.character(net.text[,1]), as.character(net.text[,2]))
	net.edge <- net.edge[net.edge[,2] != "", ]
	net.edge <- net.edge[net.edge[,1] != "", ]
	net.matrix <- matrix(0, net.size, net.size, dimnames=list(net.node, net.node))
	net.matrix[net.edge] <- as.numeric(net.text[,3])
	list(size=net.size, node=net.node, matrix=net.matrix)

}

geneListMap <- function(glist,maptable,C=c(3,2)){
	
	tmp  <- maptable[,C[1]] %in% glist
	Map <- maptable[tmp,C]
	
	Map
	
}

testData <- function(){
	
	source("RWR.R")
	geneList <- "genelist.txt"
	mutatedGenes <- "TCGAmutatedgenes.txt"
	network <- "STRING"
	InterfaceTest(geneList,mutatedGenes,network,outfile="TCGA")
	
	geneList <- "genelist.txt"
	mutatedGenes <- "TCGAmutatedgenes.txt"
	network <- "FI"
	InterfaceTest(geneList,mutatedGenes,network,outfile="TCGA")
	
	geneList <- "genelist.txt"
	mutatedGenes <- "GWASmutatedgenes.txt"
	network <- "STRING"
	InterfaceTest(geneList,mutatedGenes,network,outfile="GWAS")
	
	geneList <- "genelist.txt"
	mutatedGenes <- "GWASmutatedgenes.txt"
	network <- "FI"
	InterfaceTest(geneList,mutatedGenes,network,outfile="GWAS")
	
}

testDataFIandPREPPI <- function(){
	
	#file <- "PREPPI human.c600.sm.hc"
	#W <- PREPPInetwork(file)
	file <- "functional network/FIs_043009.txt"	
	W <- FInetwork(file)
	
	nodes <- rownames(W)
	n.nodes <- length(nodes)
	normal.vector <- colSums(W)

	normal.vector[normal.vector==0] <- 1
	normal.matrix <- matrix(normal.vector,n.nodes,n.nodes,byrow=TRUE)
	W <- W/normal.matrix
	
			
	#TCGAfile <- "TCGAmutatedgenes.txt"
	TCGAfile <- "GWASmutatedgenes.txt"
	genelist <- "genelist.txt"
	nodeTCGA <- as.matrix(read.table(TCGAfile))
	nodeGene <- as.matrix(read.table(genelist))[,2]
	
	subs <- nodes %in% nodeTCGA
	p0 <- 1:n.nodes
	p0[subs] <- 1/sum(subs)
	#p0[subs] <- 1
	p0[!subs] <- 0
	
	d <- p0
	
	alpha <- 0.9
	p1 <- RWR(W,p0,alpha,d)
	
	X <- sort(p1,decreasing=TRUE,index.return=TRUE)
	XTable <- cbind(nodes[X$ix],p1[X$ix],1:n.nodes)
	
	nodesOrder <-  XTable[,1]
	subs <- nodesOrder %in% nodeGene
	mapgene <- nodesOrder[subs]

	
	nodeGenemap <- as.matrix(read.table(genelist))
	genesymbol <- nodeGenemap[match(mapgene,nodeGenemap[,2]),1]
	
	resultTable <-  cbind(genesymbol,mapgene,XTable[subs,c(2,3)])
	
	write.table(XTable,file="FInetworkGWASall",quote=FALSE,row.names=FALSE,col.names=FALSE)
	write.table(resultTable,file="FInetworkGWASbased",quote=FALSE,row.names=FALSE,col.names=FALSE)
	
	
	resultTable
}

tmpPlot <- function(){
	source("RWR.R")
	a <- read.table("STRINGnetworkall")
	b <- read.table("STRINGnetworkbased")
	hist(as.numeric(a[,2]),main="The whole network probability distribution",xlab="Probability")

	hist(as.numeric(b[,4]),main="Our gene list probability distribution",xlab="Probability")
	
	
	a <- read.table("FInetworkall")
	b <- read.table("FInetworkbased")
	hist(as.numeric(a[,2]),main="The whole network probability distribution",xlab="Probability")

	hist(as.numeric(b[,3]),main="Our gene list probability distribution",xlab="Probability")
	
	
	a <- read.table("STRINGnetworkGWASall")
	b <- read.table("STRINGnetworkGWASbased")
	hist(as.numeric(a[,2]),main="The whole network probability distribution",xlab="Probability")

	hist(as.numeric(b[,4]),main="Our gene list probability distribution",xlab="Probability")
	
	
	a <- read.table("FInetworkGWASall")
	b <- read.table("FInetworkGWASbased")
	hist(as.numeric(a[,2]),main="The whole network probability distribution",xlab="Probability")

	hist(as.numeric(b[,3]),main="Our gene list probability distribution",xlab="Probability")
	
}