CCBA_ssGSEA_project_dataset.v1 <- function(
input.ds,
output.ds,
gene.set.databases,
gene.set.selection  = "ALL",   # "ALL" or list with names of gene sets
sample.norm.type    = "rank",  # "rank", "log" or "log.rank"
weight              = 0.25,
statistic           = "area.under.RES",
output.score.type   = "ES",    # "ES" or "NES"
save_norm_dataset   = NULL,
nperm               = 200,     # number of random permutations for NES case
combine.mode        = "combine.off",  # "combine.off" do not combine *_UP and *_DN versions in 
min.overlap         = 1,
gene.names.in.desc  = F,      # in Protein, RNAi Ataris or hairpin gct files the gene symbols are in the descs column
correl.type         = "rank") # "rank", "z.score", "symm.rank"
{ 

# Sample normalization
if (sample.norm.type == "rank") {
	for (j in 1:Ns) {  # column rank normalization 
		m[,j] <- rank(m[,j], ties.method = "average")
	}
	m <- 10000*m/Ng
} else if (sample.norm.type == "log.rank") {
	for (j in 1:Ns) {  # column rank normalization 
		m[,j] <- rank(m[,j], ties.method = "average")
	}
	m <- log(10000*m/Ng + exp(1))
} else if (sample.norm.type == "log") {
	m[m < 1] <- 1
	m <- log(m + exp(1))
}

# Loop over gene sets
score.matrix <- score.matrix.2 <- matrix(0, nrow=N.gs, ncol=Ns)
print(paste("Size score.matrix:", dim(score.matrix)))
print(paste("Size score.matrix.2:", dim(score.matrix.2)))                
for (gs.i in 1:N.gs) {
	#browser()
	gene.set <- gs[gs.i, 1:size.G[gs.i]]
	gene.overlap <- intersect(gene.set, gene.names)
	print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", length(gene.overlap)))
	if (length(gene.overlap) < min.overlap) { 
		score.matrix[gs.i, ] <- rep(NA, Ns)
		print(paste("Size score.matrix:", dim(score.matrix)))                        
		next
	} else {
		gene.set.locs <- match(gene.overlap, gene.set)
		gene.names.locs <- match(gene.overlap, gene.names)
		msig <- m[gene.names.locs,]
		msig.names <- gene.names[gene.names.locs]
		if (output.score.type == "ES") {
			OPAM <- CCBA_ssGSEA.Projection.v1(data.array = m, gene.names = gene.names, n.cols = Ns, 
					n.rows = Ng, weight = weight, statistic = statistic,
					gene.set = gene.overlap, nperm = 1, correl.type = correl.type)
			score.matrix[gs.i,] <- as.matrix(t(OPAM$ES.vector))
			print(paste("Size score.matrix:", dim(score.matrix)))                                
		} else if (output.score.type == "NES") {
			OPAM <- CCBA_ssGSEA.Projection.v1(data.array = m, gene.names = gene.names, n.cols = Ns, 
					n.rows = Ng, weight = weight, statistic = statistic,
					gene.set = gene.overlap, nperm = nperm, correl.type = correl.type)
			score.matrix[gs.i,] <- as.matrix(t(OPAM$NES.vector))
			print(paste("Size score.matrix:", dim(score.matrix)))                                
		}
	}
}


CCBA_ssGSEA.Projection.v1 <- function(
data.array,
gene.names,
n.cols,
n.rows,
weight = 0,
statistic    = "Kolmogorov-Smirnov",
gene.set,
nperm = 200,
correl.type  = "rank")
{

ES.vector <- vector(length=n.cols)
NES.vector <- vector(length=n.cols)
p.val.vector <- vector(length=n.cols)
correl.vector <- vector(length=n.rows, mode="numeric")

# Compute ES score for signatures in each sample

#   print("Computing GSEA.....")
phi <- array(0, c(n.cols, nperm))
for (sample.index in 1:n.cols) {
	gene.list <- order(data.array[, sample.index], decreasing=T)            
	gene.set2 <- match(gene.set, gene.names)
	
	if (weight == 0) {
		correl.vector <- rep(1, n.rows)
	} else if (weight > 0) {
		if (correl.type == "rank") {
			correl.vector <- data.array[gene.list, sample.index]
		} else if (correl.type == "symm.rank") {
			correl.vector <- data.array[gene.list, sample.index]
			correl.vector <- ifelse(correl.vector > correl.vector[ceiling(n.rows/2)], 
					correl.vector,
					correl.vector + correl.vector - correl.vector[ceiling(n.rows/2)]) 
		} else if (correl.type == "z.score") {
			x <- data.array[gene.list, sample.index]
			correl.vector <- (x - mean(x))/sd(x)
		}
	}
	### Olga's Additions ###
#		ptm.new = proc.time()
	tag.indicator <- sign(match(gene.list, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
	no.tag.indicator <- 1 - tag.indicator 
	N <- length(gene.list) 
	Nh <- length(gene.set2) 
	Nm <-  N - Nh 
	orig.correl.vector <- correl.vector
	if (weight == 0) correl.vector <- rep(1, N)   # unweighted case
	ind = which(tag.indicator==1)
	correl.vector <- abs(correl.vector[ind])^weight
	
	
	sum.correl = sum(correl.vector)
	up = correl.vector/sum.correl     # "up" represents the peaks in the mountain plot
	gaps = (c(ind-1, N) - c(0, ind))  # gaps between ranked pathway genes
	down = gaps/Nm
	
	RES = cumsum(c(up,up[Nh])-down)
	valleys = RES[1:Nh]-up
	
	max.ES = max(RES)
	min.ES = min(valleys)
	
	if( statistic == "Kolmogorov-Smirnov" ){
		if( max.ES > -min.ES ){
			ES <- signif(max.ES, digits=5)
			arg.ES <- which.max(RES)
		} else{
			ES <- signif(min.ES, digits=5)
			arg.ES <- which.min(RES)
		}
	}
	
	if( statistic == "area.under.RES"){
		if( max.ES > -min.ES ){
			arg.ES <- which.max(RES)
		} else{
			arg.ES <- which.min(RES)
		}
		gaps = gaps+1
		RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES[1:Nh]) - c(valleys,0) ) * (gaps)
		ES = sum(RES)
	}
	GSEA.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator)
#		new.time <<- new.time + (proc.time() - ptm.new)
	### End Olga's Additions ###
	#GSEA.results <- GSEA.EnrichmentScore5(gene.list=gene.list, gene.set=gene.set2,
	#		statistic = statistic, alpha = weight, correl.vector = correl.vector)
	ES.vector[sample.index] <- GSEA.results$ES
	
	if (nperm == 0) {
		NES.vector[sample.index] <- ES.vector[sample.index]
		p.val.vector[sample.index] <- 1
	} else {
		for (r in 1:nperm) {
			reshuffled.gene.labels <- sample(1:n.rows)
			if (weight == 0) {
				correl.vector <- rep(1, n.rows)
			} else if (weight > 0) {
				correl.vector <- data.array[reshuffled.gene.labels, sample.index]
			} 
#				GSEA.results <- GSEA.EnrichmentScore5(gene.list=reshuffled.gene.labels, gene.set=gene.set2,
#						statistic = statistic, alpha = weight, correl.vector = correl.vector)
			### Olga's Additions ###
			tag.indicator <- sign(match(reshuffled.gene.labels, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
			no.tag.indicator <- 1 - tag.indicator 
			N <- length(reshuffled.gene.labels) 
			Nh <- length(gene.set2) 
			Nm <-  N - Nh 
#   orig.correl.vector <- correl.vector
			if (weight == 0) correl.vector <- rep(1, N)   # unweighted case
			ind <- which(tag.indicator==1)
			correl.vector <- abs(correl.vector[ind])^weight   
			
			sum.correl <- sum(correl.vector)
			up = correl.vector/sum.correl
			gaps = (c(ind-1, N) - c(0, ind))
			down = gaps/Nm
			
			RES = cumsum(c(up,up[Nh])-down)
			valleys = RES[1:Nh]-up
			
			max.ES = max(RES)
			min.ES = min(valleys)
			
			if( statistic == "Kolmogorov-Smirnov" ){
				if( max.ES > -min.ES ){
					ES <- signif(max.ES, digits=5)
					arg.ES <- which.max(RES)
				} else{
					ES <- signif(min.ES, digits=5)
					arg.ES <- which.min(RES)
				}
			}
			
			if( statistic == "area.under.RES"){
				if( max.ES > -min.ES ){
					arg.ES <- which.max(RES)
				} else{
					arg.ES <- which.min(RES)
				}
				gaps = gaps+1
				RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES[1:Nh]) - c(valleys,0) ) * (gaps)
				ES = sum(RES)
			}
			
			GSEA.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator)
			### End Olga's Additions ###
			phi[sample.index, r] <- GSEA.results$ES
		}
		if (ES.vector[sample.index] >= 0) {
			pos.phi <- phi[sample.index, phi[sample.index, ] >= 0]
			if (length(pos.phi) == 0) pos.phi <- 0.5
			pos.m <- mean(pos.phi)
			NES.vector[sample.index] <- ES.vector[sample.index]/pos.m
			s <- sum(pos.phi >= ES.vector[sample.index])/length(pos.phi)
			p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
		} else {
			neg.phi <-  phi[sample.index, phi[sample.index, ] < 0]
			if (length(neg.phi) == 0) neg.phi <- 0.5 
			neg.m <- mean(neg.phi)
			NES.vector[sample.index] <- ES.vector[sample.index]/abs(neg.m)
			s <- sum(neg.phi <= ES.vector[sample.index])/length(neg.phi)
			p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
		}
	}
return(list(ES.vector = ES.vector, NES.vector =  NES.vector, p.val.vector = p.val.vector))
}
