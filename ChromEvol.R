#ChromEvol by N. Cusimano
	#cusimano@bio.lmu.de
	
	
# 1. PREPARING INPUT DATA 
	###ce.input reads a tree in newick format from file or in phylo format from R and a table from simple text file with two colums, 1st with species names, 2nd with haploid chromosome numbers or from fasta file. 
	#Writes: 
		# a tree and the table input file with chromosome numbers that have matching species names and that can be directly used as input for chromEvol
		# the params.txt file with which you start the analysis. It uses the current working directory


#2. READ DATA 
#After analysis read the results obtained by fitting a specific model to the data with 
	###read.cE just by giving the path of the folder of the respective model
		#value is an object of class "ChromEvol'
			#list of 3....

#  	3. PLOT RESULTS
#the following functions allow you to get a good idea of your results in a short time, with
	###plot.anc: haploid numbers at nodes from ML reconstruction or best numbers obtained from Bayesian inference together with pie charts giving the probabilities of the different numbers are plotted
	###plot.EVENTS: plots the numbers of all expected events along branches OR all events inferred with a PP > 0.5 along branches
	###plot.noNo: plots node numbers

#	4. SUMMARIZE PARAMETERS
#with
	###sum.Models you get a summary of all relevant parameters of either all models, or only the best ones, i.e. those that differ max. 2 units from that one with the lowest AIC value. Input is the path of the folder(_outDir) of the analysis that have been run with the "_mainType All_Models" option
	
	
######################################
#ce.input 
#writes tree, table and parameters file in the current working directory	
	
ce.input <- function(phy,table, missing=NA, filename="ChromEvol_analysis", folder=NA, percentage=TRUE, mainType="All_Models", fix.root=FALSE, maxChrNum = "-10", minChrNum = "1", branchMul = "1", rootAt = "N1", numberSim=NA, ladderize=TRUE){
	if (class(phy) != "phylo")
		phy <- read.tree(phy)
	if (!is.null(phy$posterior) && is.null(phy$node.label))  
			phy$node.label <- round(phy$posterior, digits=2)  
	if (ladderize)
		phy <- ladderize(phy)
	phy <- fixNodes(phy)
	if (!is.matrix(table)){
		tableS <- scan(table, what ="list")
		if (length(grep(">", tableS))>1){
			table <- read.cE.tab(table)
			print(table)}
		else
			table <- read.table(table, stringsAsFactors = FALSE)
			
		}
	table2 <- match.table2tip(phy, table)
	table2[,1] <- phy$tip.label
		#drop tips not in table
	if (length(which(is.na(table2[,2])))!= 0){
		if (is.na(missing)){
			phy2 <- drop.tip(phy, which(is.na(table2[,2])),1)	
			table3 <- table2[-which(is.na(table2[,2])),]
		}
	
		else{
			table2[which(is.na(table2[,2])),2] <- "X"
			table3 <- table2
			phy2 <- phy
		}
	}
	else {
		table3 <- table2
		phy2 <- phy
		}
	
	#get species with multiple counts
	
	if (percentage){
		multi <- grep("_", table3[,2])
		if(length(multi) > 0){
			perc <- grep("=", table3[,2])
			if (length(perc) == 0){
				multi2 <- strsplit(table3[multi,2], "_")
				for (i in 1:length(multi2)){
					p <- round(1/length(multi2[[i]]), digits=2)
					if (length(multi2[[i]])*p != 1)
						p1 <- (1-length(multi2[[i]])*p)+p
					else
						p1 <- p
					ll <- paste(multi2[[i]][1], "=", p1,  sep="")
					for (j in 2:(length(multi2[[i]]))){
						ll <- paste(ll,"_",multi2[[i]][j], "=", p,  sep="")
					}	
					multi2[[i]] <- ll 
					table3[multi[i],2] <- multi2[[i]]
				}
			}
		}
	}
	if (is.na(numberSim)){
		if(mainType == "All_Models")
			numberSim <- 0
		else
			numberSim <- 10000
	}
	if (!is.na(folder))
		setwd(folder)
	
	
	write.table(paste(">", table3[1,1], "\n", table3[1,2]), file=paste(filename, "_input_table.txt", sep=""), quote=FALSE, row.names = FALSE, col.names = FALSE)
	for (i in 2:length(table3[,1]))
		write.table(paste(">", table3[i,1], "\n", table3[i,2]), file=paste(filename, "_input_table.txt", sep=""), quote=FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)
	fnp <- paste(filename, "_params.txt", sep="")
	write.tree(phy2, file=paste(filename, "_with_bootstraps.tree", sep=""))
	phy2$node.label <- NULL
	write.tree(phy2, file=paste(filename, "_input.tree", sep=""))
	write.param(filename, mainType, maxChrNum, minChrNum, branchMul, rootAt, numberSim, fix.root)
	if (mainType != "All_Models")
		specify.model(mainType, filename)
	cat("\tStructure of the input table\n")
	system(paste("open", fnp))
	list(table3, phy2)
	cat("number of tips of the original tree:", length(phy$tip.label), "\n",#
		"number of tips of the ChromEvol input tree:", length(phy2$tip.label), "\n",#
		length(table[,1])-length(table3[,1]), "species have been dropped from the table\n",#
		length(phy$tip.label)-length(phy2$tip.label), "tips have been dropped from tree:\n")
	cat(phy$tip.label[which(is.na(table2[,2]))], "\n",sep="\n")
	if(!is.logical(fix.root))
		write.fix.root(fix.root)
	setwd("~/R")
	
	invisible(list(table3, phy2))
	}
		
####################
write.param <- function(filename, mainType,  maxChrNum, minChrNum, branchMul, rootAt, numberSim=numberSim, fix.root){
	wd <- getwd()
	if (mainType == "All_Models")
		write(paste("_mainType", mainType), file=paste(filename, "_params.txt", sep=""))
	else{
		if (length(grep("fix_", mainType)) !=1)
			write("_mainType Optimize_Model", file=paste(filename, "_params.txt", sep=""))
		else
			write("_mainType Run_Fix_Param", file=paste(filename, "_params.txt", sep=""))
	}
	write(paste("_outDir", wd), file=paste(filename, "_params.txt", sep=""), append=TRUE)
	write(paste("_dataFile ",wd, "/", filename, "_input_table.txt", sep=""), file=paste(filename, "_params.txt", sep=""), append=TRUE)
	write(paste("_treeFile ",wd, "/", filename, "_input.tree", sep=""), file=paste(filename, "_params.txt", sep=""), append=TRUE)
	write(paste("_logFile ", filename, "_log1.txt", sep=""), file=paste(filename, "_params.txt", sep=""), append=TRUE)
	if (!is.logical(fix.root)){
		write(paste("_freqFile ",  wd, "/freq.txt", sep="" ), file=paste(filename, "_params.txt", sep=""), append=TRUE)
		write("_rootFreqType FIXED", file=paste(filename, "_params.txt", sep=""), append=TRUE)
		}
	write(paste("_maxChrNum", maxChrNum), file=paste(filename, "_params.txt", sep=""), append=TRUE)
	write(paste("_minChrNum", minChrNum), file=paste(filename, "_params.txt", sep=""), append=TRUE)
	write(paste("_rootAt", rootAt), file=paste(filename, "_params.txt", sep=""), append=TRUE)
	write(paste("_branchMul", branchMul), file=paste(filename, "_params.txt", sep=""), append=TRUE)
	write(paste("_simulationsNum", numberSim), file=paste(filename, "_params.txt", sep=""), append=TRUE)
	#	write("", file=paste(filename, "_params.txt", sep=""), append=TRUE)
	}
#####################
specify.model <- function(model, filename){
	model <- gsub("fix_","", model)
	if (model == "cr")
			write("_gainConstR\n_lossConstR\n_duplConstR\n_demiPloidyR -999\n_gainLinearR -999\n_lossLinearR -999", file=paste(filename, "_params.txt", sep=""), append=TRUE)
	if (model == "crd")
		write("_gainConstR\n_lossConstR\n_duplConstR\n_demiPloidyR -2\n_gainLinearR -999\n_lossLinearR -999", file=paste(filename, "_params.txt", sep=""), append=TRUE)
	if (model == "crde")
		write("_gainConstR\n_lossConstR\n_duplConstR\n_demiPloidyR\n_gainLinearR -999\n_lossLinearR -999", file=paste(filename, "_params.txt", sep=""), append=TRUE)
	if (model == "crnd")
		write("_gainConstR\n_lossConstR\n_duplConstR -999\n_demiPloidyR -999\n_gainLinearR -999\n_lossLinearR -999", file=paste(filename, "_params.txt", sep=""), append=TRUE)
	if (model == "lr")
		write("_gainConstR\n_lossConstR\n_gainLinearR\n_lossLinearR\n_duplConstR\n_demiPloidyR -999", file=paste(filename, "_params.txt", sep=""), append=TRUE)
	if (model == "lrd")
		write("_gainConstR\n_lossConstR\n_gainLinearR\n_lossLinearR\n_duplConstR\n_demiPloidyR -2", file=paste(filename, "_params.txt", sep=""), append=TRUE)
	if (model == "lrde")
		write("_gainConstR\n_lossConstR\n_gainLinearR\n_lossLinearR\n_duplConstR\n_demiPloidyR", file=paste(filename, "_params.txt", sep=""), append=TRUE)
	if (model == "lrnd")
		write("_gainConstR\n_lossConstR\n_gainLinearR\n_lossLinearR\n_duplConstR -999\n_demiPloidyR -999", file=paste(filename, "_params.txt", sep=""), append=TRUE)
	
}

write.fix.root <- function(x){
	if (length(x)==1)
		write(paste("F[", x, "]=1", sep=""), file="freq.txt")
	else{
		for (i in 1:length(x))
			write(paste("F[", x[i], "]=", sep=""), file="freq.txt", append=TRUE)
		system("open freq.txt")
		cat("	ADD THE PROBABILITIES TO THE FIXED ROOT NODE NUMBERS IN THE freq.txt FILE (ALREADY OPENED)\n")
		}
	}


#####################

read.ce <- function(path){
	tree <- read.mlAnc.tree(paste(path, "/mlAncestors.tree", sep=""))
	ppTR <- read.PP.tree(paste(path, "/posteriorAncestors.tree", sep=""))
	tree$PP <- ppTR$node.label
	
	anc <- read.ancProb.tab(paste(path, "/ancestorsProbs.txt", sep=""))
	exp <- read.expect.tab(paste(path, "/expectations.txt", sep=""))
	param <- read.param(path, optimize=TRUE)
	cE <- list(POSTERIOR_PROBABILITIES_of_CN=anc, EVENTS_along_edges=exp, TREE_with_ML_reconstruction_as_node.label=tree, Parameters=param)
	class(cE) <- "ChromEvol"
	cE
	}
	
	
	
##########	
	
read.mlAnc.tree <- function(filename){
	phy <- read.tree(file=filename)
	phy$node.number <- gsub("-.*[[:punct:]]", "",phy$node.label)
	phy$node.number <- gsub("[[:punct:]]", "", phy$node.number)
	phy$node.label <- gsub("[[:punct:]]N.*-", "",phy$node.label)
	phy$node.label <- gsub("[[:punct:]]", "",phy$node.label)
	phy
	}

########

read.PP.tree <- function(filename){
	phy <- read.tree(file=filename)
	phy$node.label <- gsub("[[:punct:]]N.*_", "",phy$node.label)
	phy$node.label <- gsub("[]]", "",phy$node.label)
	phy$node.label  <- gsub("-", ":",phy$node.label)
	phy$node.label  <- gsub("//", "-",phy$node.label)
	phy$node.label  <- gsub("-.*:0$", "",phy$node.label)
	phy$node.label  <- gsub("-.*:0.0[0-9]$", "",phy$node.label)
	phy
	}	
	
read.ancProb.tab <- function(filename){
	anc.tab <- read.table(filename, header=TRUE, stringsAsFactors = FALSE)
	anc.tab2 <- as.matrix(anc.tab)
	anc.tab3 <- apply(anc.tab2, 2, as.numeric)
	dimnames(anc.tab3)[[1]] <- anc.tab2[,1]
	anc.tab3 <- anc.tab3[,-1]
	dimnames(anc.tab3)[[2]] <- gsub("X", "",dimnames(anc.tab3)[[2]])
	ancNodes <- anc.tab3[grep("^N[0-9]",dimnames(anc.tab3)[[1]]),]
	ancTips <- anc.tab3[-(grep("^N[0-9]",dimnames(anc.tab3)[[1]])),]
	foo <- function(x) as.numeric(dimnames(ancNodes)[[2]][which(x == max(x))])
	pp <- apply(ancNodes, 1,foo)
	anc <- list(Nodes=ancNodes, Tips=ancTips, Nodes_best_CN=pp)
	anc
	}
	
	
####################
#####################

#######
read.expect.tab <- function(filename){
	exp <- scan(filename, what="list", sep="\n")
##
	text <- c("#Nodes with GAIN events with expectation above 0.5", "#Nodes with LOSS events with expectation above 0.5","#Nodes with duplication events with expectation above 0.5", "#Nodes with demi-duplication events with expectation above 0.5", "#+++++++++++++++++++++++++++++")
	events <- list()
	EVENTS <- list()
	EVENTS2 <- list()
	for (j in 1:4){
		events[[j]] <- exp[(which(exp[1:length(exp)] == text[j])+1):(which(exp[1:length(exp)] == text[5])[j]-1)]

		events[[j]] <- strsplit(events[[j]], split=": ")
 	
	 	ev <- NULL
		for (i in 1:length(events[[j]]))
			ev <- rbind(ev,events[[j]][[i]])
		EVENTS[[j]] <- ev	
		}
####
	write.table(exp[(which(exp[1:length(exp)] == "#ALL EVENTS EXPECTATIONS PER NODE")+2):(which(exp[1:length(exp)] == "#Expected number of events from root to leaf")-2)], sep="/t",file="exp.tab")
	exp.tab <- read.table("exp.tab", stringsAsFactors = FALSE, skip=1, quote="",sep="\t")
	exp.tab[,1] <- gsub("\"", "",exp.tab[,1])
	exp.tab[,5] <- as.numeric(gsub("\"", "",exp.tab[,5]))
	exp.tab2 <- exp.tab
	foo <- function(x)
		round(x, digits=3)
	exp.tab2[2:5] <- lapply(exp.tab[2:5], foo)
	#for (i in 2:5)
		#exp.tab2[which(exp.tab2[,i]==0),i] <- NA
	exp.tab2[,1] <- gsub("[0-9].*/t", "", exp.tab2[,1])

#Read heuristic estimations for events not accounted for in the simulation
	he <- exp[(which(exp[1:length(exp)] == "#EVENTS NOT ACCOUNTED FOR IN THE SIMULATIONS: ")+2):(which(exp[1:length(exp)] == "#ALL EVENTS EXPECTATIONS PER NODE")-2)]
	if(length(grep("HEURISTIC ESTIMATION", he)) > 0){
	he2 <- strsplit(he, split="\t")
	ev <- c("Gains = ", "Losses = ", "Duplications = ", "Demi-dupl = ")
	HeurEst <- data.frame()
	for (i in grep("Gains", he2)){
		for(j in 3:6)
			HeurEst[i,j-2] <- as.numeric(gsub(ev[j-2], "", he2[[i]][j]))
	}
	HeurEst <- HeurEst[grep("Gains", he2),]
	}
	else HeurEst <- data.frame(rbind(rep(NA,4)))
	###########

	list(Events_with_high_PP=EVENTS, All_events=exp.tab2, Heuristic_Estimation=HeurEst)
	}
#####	


####################
#####################

plot.ce <- function(chromEvol, cex.tips=0.6, offset = 0.004,...){
	plot(chromEvol[[3]], cex=cex.tips, no.margin=TRUE, label.offset=offset, ...)
	}



####################
## #####################
#ChromEvol: object of class "ChromEvol"

#noNo: nodeNumbers: indicates which kind of inferred ancestral numbers should be plotted at nodes:  "ML": the inferred ML numbers are plotted; "PP": the numbers inferred under Bayesian optimization with the highest posterior probability;  "2bPP": the 2 numbers with the highest and second-highest PP will be plotted together with the respective PP; if the best number has a PP>0.9 the second-best will be omitted; NA: no number will be plotted

#pie: if TRUE, pie charts of posterior probabilites of the chromosome numbers are plotted at nodes

#tip.pie: pie charts with PP of tip numbers are plotted between the terminal branches of of the tree and the species names

#inferred.tip: if TRUE inferred chromosome number for the tips with the highest probability wil be plotted next to the tiplabels

#label.offset: see plot.phylo

#adj.tip.pie:

#cut.off:

#cex.tip, cex.pie, cex.node.No: a numeric value giving the factor scaling of the tip and node labels and pie charts (Character EXpansion). The default is to take the current value from the graphical parameters.

#add: if FALSE the tree of the analysis will be plotted together with the node labels and pie charts;  if TRUE pie charts and/or node labels will be added to an exitent plot.

#title: plots a title above the phylogeny

#x,y, position of legend, if y=NA  the legend will be placed by clicking on the figure with locator()

#...: further graphical parameters to be passed to plot.phylo

plot.anc <- function(ChromEvol, node.pie=TRUE, node.No="PP", info=FALSE, tip.pie=TRUE, inferred.tip=FALSE, chroCOL=NULL, cex.node.pie=0.8, cut.off=0.1, cex.tip=0.8, cex.tip.pie=0.3, cut.tipNo=FALSE, adj.tip.pie=0.5,cex.node.No=0.6, cex.legend=0.8, add=FALSE, label.offset=NA, title=NULL, cex.title=2, frame="c", x=0, y=NA, only.legend=FALSE, ...){
	if (class(ChromEvol) == "ChromEvol"){
		ChromEvol.anc <- ChromEvol[[1]]
		anc <- ChromEvol.anc[[1]]
		tip <- ChromEvol.anc[[2]]
		phy <- ChromEvol[[3]]
		pp <- ChromEvol.anc[[3]]
		tipNO <- gsub(".*-", "", phy$tip.label)
		multiples <- grep("=", tipNO)
		tipNO <- as.numeric(gsub("=.*", "", tipNO))
		if(cut.tipNo)
			phy$tip.label <- gsub("-.*", "", phy$tip.label)
		 foo <- function(x)
			as.numeric(dimnames(tip)[[2]][which(x == max(x))])
		bestTip <- apply(tip,1,foo)

		chroNO <- sort(as.numeric(unique(c(unique(dimnames(anc)[[2]][which(anc > 0.1, arr.ind=TRUE)[,2]]), unique(dimnames(tip)[[2]][which(tip > 0.1, arr.ind=TRUE)[,2]])))))
		#unique(which(anc > cut.off, arr.ind=TRUE)[,2])
		COL <- rep("black", length(anc[1,]))
		if (max(as.numeric(dimnames(anc)[[2]])) != length(anc[1,]))
			d <- max(as.numeric(dimnames(anc)[[2]])) -length(anc[1,])	
		else d <- 0
		if (is.null(chroCOL)){
			chroCOL <- c("grey", rainbow(length(chroNO)-1, start=.9, end=.8))
			if (length(chroCOL) > 11)			
				chroCOL[c(11,12)] <- c("green3","darkolivegreen3")
			COL[chroNO-d] <- chroCOL
			}
		else
			COL <- chroCOL
	if (is.na(label.offset)){
		if (tip.pie || inferred.tip)
			label.offset <- 0.004##
		else label.offset <- 0
		}
	if (!only.legend){	
		if (!add)
			plot(phy, cex=cex.tip,  label.offset=label.offset, ...)
		if (node.pie)
			nodelabels(pie=anc, piecol=COL, cex=cex.node.pie)
		if (tip.pie)
			tiplabels(pie=tip, cex= cex.tip.pie, piecol=COL, adj=adj.tip.pie, font=2)
		
				
		if (node.No == "ML"){
			nodelabels(phy$node.label, cex=cex.node.No, font=2,  frame=frame, col="black", bg="white")
			if (info)
				mtext("Numbers at nodes: ML reconstruction", adj=0, line=-1, cex=0.8)
			}
		if (node.No =="PP"){
			nodelabels(pp, cex=cex.node.No, font=2,  frame=frame, col="black", bg="white")
			if (info)
				mtext("Numbers at nodes:  ChrNo with highest PP from Bayesian optimization", line=-1, adj=0, cex=0.8)
				}
		if (node.No =="2bPP"){
			nodelabels(phy$PP, cex=cex.node.No, font=2,  frame="r", col="black", bg="white")
			if (info)
				mtext("Numbers at nodes: 2 ChrNos with highest PP from Bayesian optimization", line=-1, adj=0, cex=0.8)
				}
		if (inferred.tip){
			tipCOL <- vector(mode="character",length=length(tipNO))
			for (i in 1:length(tipNO))
				tipCOL[i] <- COL[tipNO[i]]
			tipNO[multiples] <- paste(tipNO[multiples], "+", sep="")
			if (isTRUE(tip.pie))
				tipsymbols(phy, text=tipNO, cex=cex.tip, col=tipCOL, differ=FALSE, symbol.offset=label.offset)
			else				
				tiplabels(tipNO, frame="n", cex= 0.6, col=tipCOL, adj=-0.5)
		}#,
	
	}
	if (node.pie){
		coll <- NULL
		maxcoll=NULL
		if(!is.null(ChromEvol$collapse)){
			for(i in 1:length(ChromEvol$collapse[,1])){
				coll <- cbind(coll, which(phy$tip.label==dimnames(ChromEvol$collapse)[[1]][i]))
				maxcoll <- cbind(maxcoll, which.max(ChromEvol$collapse[i,]))
				}
			
			tiplabels(tip=coll, pie=ChromEvol$collapse, piecol=COL, cex=cex.node.pie)
			tiplabels(tip=coll, text=maxcoll, cex=cex.node.No, font=2,  frame=frame, col="black", bg="white")
	}
	}
if (node.pie || tip.pie || only.legend){
		if (is.na(y) && is.numeric(x)){
			cat("Bring the plotting device to the front and click on the plot where the upper left corner of the legend should be placed")
			xy <- locator(1)
			x <- xy[[1]]
			y <- xy[[2]]
			}
		else{
			if (!is.numeric(x))
				y <- NULL
			}
			legend(x, y,legend=chroNO, fill=COL[chroNO-d], title="n =", cex=cex.legend)
			}
	mtext(title, line=-2, adj=0, cex=cex.title, font=2)
			}
	else stop()
	invisible(COL)
	}
	
########################
###plots node numbering as used by ChromEvol

plot.noNO <- function(ChromEvol, ...){
	plot(ChromEvol[[3]], cex=0.6, no.margin=TRUE, ...)
	nodelabels(ChromEvol[[3]]$node.number, frame="n", cex=0.8)

	}
	
########################
#### PLOTTING EVENTS #####


plot.EVENTS <- function(ChromEvol, events="best", frame=FALSE, info=TRUE, cex.tip=0.6, cex.event=0.6, cex.node.No=0.4, cex.legend=0.8, chromNo="ML", add=FALSE, position="above", title=NULL, x=0, y=NA, only.legend=FALSE, ...){
	if (class(ChromEvol) == "ChromEvol"){
		ChromEvol.exp <- ChromEvol[[2]]
		exp.table <- ChromEvol.exp[[2]]
		EVENTS.tab <- ChromEvol.exp[[1]]		
		phy <- ChromEvol[[3]]
		pp <- ChromEvol[[1]][[3]]
		col.EVENTS <- c( rainbow(4)[4], rainbow(4)[1], rainbow(4)[3],"darkorange")
		
		if (!add)
			plot(phy, cex=cex.tip, ...)
		else
			chromNo <- FALSE
		events.tab <- sum.Events(ChromEvol)
		text.col <- c(rep("black",3),"white")
		EVENTS.tab2 <- list()
		edges2 <- list()
		legend.text <- c("Chromosome gains", "Chromosome losses", "Duplications", "Demiduplications")
		legend.text2 <- NULL
		sumEvents <- 0
		col.eventsInframes <- c("white", rep("black",3))
		if (!is.logical(chromNo))
			if (chromNo =="ML"){
				nodelabels(ChromEvol[[3]]$node.label, cex=cex.node.No, font=2,  adj=c(1.1,-0.3),frame="n", col="black", bg="white")
				if (info)
					mtext("Numbers at nodes: ML reconstruction", adj=0, line=-1, cex=0.8)
			}
			if (chromNo =="PP"){
				nodelabels(pp, cex= cex.node.No, adj=c(1.1,-0.3), font=2,  frame="n", col="black")
				if (info)
				mtext("Numbers at nodes:  ChrNo with highest PP from Bayesian optimization", line=-1, adj=0, cex=0.8)
			}
		if (position == "above")
			ad <- -0.3
		else
			ad <- 1.5	
		if (events =="best"){
			for(i in 1:4){
				b <- exp.table[,i+1]
				b[which(exp.table[,i+1]<=0.5)] <- NA
				if (!only.legend)
					edgelabels(round(b, digits=1),  adj = c(3-i, ad), frame="n", cex=cex.event, font=2, col=col.EVENTS[i])
				legend.text[i] <- paste(legend.text[i], ": ", events.tab[,3][i])
				}
			titleL <- paste("Events inferred with an exp. > 0.5:", sum(events.tab[-5,3]))	
			}
			if (events == "SIMbest"){	
				for (i in 1:4){
					if (length(EVENTS.tab[[i]][,1])>1){
						w <- NULL
						EVENTS.tab[[i]] <- cbind(EVENTS.tab[[i]], NA)
						for (j in 1:(length(EVENTS.tab[[i]][,1])-1))
							EVENTS.tab[[i]][j,3] <- which(exp.table[,1]==EVENTS.tab[[i]][j,1])
						for (h in length(EVENTS.tab[[i]][,1])){
							lel <- length(phy$edge.length)
							lmiss <- length(h:lel)
							xx <- 1:lel
							yy <- EVENTS.tab[[i]][-h,3]
							zz <- xx[-as.numeric(yy)]
							EVENTS.tab2[[i]] <- rbind(EVENTS.tab[[i]][-h,], cbind(rep(NA, lmiss), rep(NA, lmiss), zz))
						}
					
						ord <- order(as.numeric(EVENTS.tab2[[i]][,3]), EVENTS.tab2[[i]][,2])
						EVENTS.tab2[[i]][ord,2]
						if (!only.legend){
							if(!frame)
								edgelabels(round(as.numeric(EVENTS.tab2[[i]][ord,2]), digits=1),  adj = c(3-i, ad), frame="n", cex=cex.event, font=2, col=col.EVENTS[i])
							else
								edgelabels(round(as.numeric(EVENTS.tab2[[i]][ord,2]), digits=1),  adj = c(3-i, ad), frame="r", col=col.eventsInframes[i], cex=cex.event, font=2, bg=col.EVENTS[i])
						}
						else EVENTS.tab2[[i]] <- NA		
					}	
					a <- events.tab[,4][i]
					legend.text[i] <- paste(legend.text[i], ": ", a)
				sumEvents <- sumEvents + a
				}
			titleL <- paste("Events inferred by sim. with exp. >0.5:", sumEvents)
			}
		
		if (events == "all"){
			for (i in 1:4){
				if (!only.legend)
					edgelabels(round(exp.table[,i+1], digits=1), adj = c(3-i, ad), cex= cex.event, font=2,frame="n", col=col.EVENTS[i])
				legend.text[i] <- paste(legend.text[i], ": ", events.tab[i,1])
				titleL <- paste("All inferred events:", round(sum(events.tab[-5,1]), digits=1))	
			}
		}
		if (is.na(y)){
			cat("Bring the Quartz window to the front and click on the plot where the upper left corner of the legend should be placed\n")
			xy <- locator(1)
			x <- xy[[1]]
			y <- xy[[2]]
			}
			
		legend(x, y, legend.text, fill=col.EVENTS, title=titleL, cex=cex.legend)
		mtext(title, line=-1)
	}
	else stop()
	events.tab
	}

#################################################################################################

################
## summarizes parameters


sum.Models <- function(path, best=TRUE, node=1, write.table=TRUE){
	modSum <- read.modSum(path)
	filename <- "/summary_ALL_models"
	if (best){
		filename <- "/summary_BEST_models"
		models <- modSum[which(modSum[,5]<=2),1]
		modSum <- modSum[which(modSum[,5]<=2),]
		}	
	else
		models <- modSum[,1]
	foo <- function(x){
		y <- paste(path, "/", x, sep="")
		y
		}
	models <- unlist(lapply(models,foo))
	
	param <- read.param(models)
	events <- read.results(models, node=node)
	PP <- get.PP(models, node=node)
	modSum <- cbind(modSum, param,PP, events)
		if (write.table)
		write.table(modSum, file=paste(path, filename, sep=""), quote=FALSE, sep="\t", row.names=FALSE)
		cat("\nThe table below has been written in the file \"", filename, "\" in the folder \"", path, "\"\n\n", sep="")
	system((paste("open", paste(path, filename, sep=""))))
	return(modSum)

	}
#################################################################################################

################
#summarizes parameters of a ChromEvol object

sum.Parameters <- function(ChromEvol, node=1, write.table=TRUE){
		#root ML
		
		CN_ML <- as.numeric(ChromEvol[[3]]$node.label[node])
		events <- sum.Events(ChromEvol)
		pp <- vector()
		pp[1] <- gsub("-.*","", ChromEvol[[3]]$PP[1])
		pp[2] <- gsub(".*-","", ChromEvol[[3]]$PP[1])
		names(pp) <- c("CN_PP_root_BEST", "CN_PP_root_2ndB")

		modSum <- list(c(ChromEvol[[4]] , pp, "CN_ML_root"=CN_ML), events)
		cat("EVENTS\t",  file="sumParameters")
		write.table(as.data.frame(events),file="sumParameters",sep="\t",quote = FALSE,append=TRUE)
		cat(paste("\n", names(modSum[[1]])[1], "\t", as.numeric(modSum[[1]][1]),"\n", sep=""), file="sumParameters", append=TRUE)
		for (i in 2:length(modSum[[1]]))
			cat(paste(names(modSum[[1]])[i],  as.vector(modSum[[1]][i]),"\n", sep="\t"), file="sumParameters", append=TRUE)
	system((paste("open", "sumParameters")))
	modSum
	}

#################################################################################################

################


###################
## reads file models_summary.txt and calculates expectations of models and differences in AIC values

read.modSum <- function(path, best=TRUE){
	modSum <- read.table(paste(path, "/models_summary.txt", sep=""), header=TRUE, stringsAsFactors = FALSE)
	x <- modSum[,3]
	min <- min(x)
	y <- NULL
	d <- NULL
	for (i in 1:length(x)){
		y[i] <- exp((min-x[i])/2)
		d[i] <- x[i]-min
		}
	modSum <- cbind(modSum, exp=round(y, digits=2), diff=round(d, digits=1))
	}
	


############################################
## reads rates of model parameters from files in specific model folder
# ...: path of model folders
read.param <- function(..., optimize=FALSE){
	models <- list(...)
	models <- unlist(models)
	SCAN <- list()
	param <- vector("list", 6)
	text <- c("LOSS_CONST","GAIN_CONST", "DUPL", "HALF_DUPL","LOSS_LINEAR", "GAIN_LINEAR")
	for (i in 1:length(models)){
		SCAN[[i]] <- scan(paste( models[[i]], "/chromEvol.res", sep=""), what="list", sep="\t")
		for (j in 1:6){
			#param[[j]] <- vector()
			if (length(grep(text[j], SCAN[[i]]))>0) 
				param[[j]][i]  <- round(as.numeric(SCAN[[i]][grep(text[j], SCAN[[i]])[1]+1]), digits=2) 
			else 
				param[[j]][i] <- NA
			}
		}
	text2 <- c("lossC", "gainC", "dupl", "demi", "lossL","gainL")
	names(param) <- paste("rate_", text2, sep="")

	if (optimize){
		ll <- gsub("LogLikelihood = ", "", SCAN[[1]][grep("LogLikelihood", SCAN[[1]])])
		AIC <- gsub("AIC.*= ", "", SCAN[[1]][grep("AIC", SCAN[[1]])])
		param <- c(AIC=as.numeric(AIC), param)
		param <- c(logLikelihood=as.numeric(ll), param)
		}
		param
	}	
	
############################################
#reads No of events with probabilities >0.5, and the chromosome number at root node inferred in an ML framework

read.results <- function(..., node=node){
	models <- list(...)
	models <- unlist(models)
	CE <- vector("list",length(models))
	ML <- vector("list", 1)
	events <- vector("list", 6)
	for (i in 1:length(models)){
		CE[[i]] <- read.cE(paste(models[[i]], sep=""))
		
		se <- sum.Events(CE[[i]])
		#root ML
		events[[1]][[i]] <- as.numeric(CE[[i]][[3]]$node.label[node])
		#events rates
		for (j in 1:4)
			events[[j+1]][i] <- se[j,3]			
		events[[6]][i] <- se[5,3]
		}
	names(events) <- c("CN_ML","No_GAINS>0.5", "No_LOSS>0.5", "No_DUPL>0.5", "No_DEMI>0.5", "totNoEvents")
	events
	}

########################################

#get numbers of the chromosome number at ROOT node inferred in a Bayesian framework, giving the two numbers with the best PP and the respective PP value

get.PP <- function(..., node=1){
	models <- list(...)
	models <- unlist(models)
	#pp.trees <- vector("list",length(models))
	pp <- vector("list", 2)
	for (i in 1:length(models)){
		pp.tree <- read.tree(paste(models[[i]], "/posteriorAncestors.tree", sep=""))
		pp.tree$node.label[node] <- gsub(paste("[[]N", node, "_", sep=""), "",pp.tree$node.label[node])
		pp.tree$node.label[node] <- gsub("[]]", "",pp.tree$node.label[node])
		pp[[1]] <- c(pp[[1]], gsub("//.*", "", pp.tree$node.label[node]))
		pp[[2]] <- c(pp[[2]], gsub(".*//", "", pp.tree$node.label[node]))
		}
	names(pp) <- c("CN_PP_root_BEST", "CN_PP_root_2ndB")
	pp
	
}
#get numbers of the chromosome number at a SPECIFIC node inferred in a Bayesian framework, giving the two numbers with the best PP and the respective PP value

get.node.info <- function(path, nodes){
	cE <- read.cE(path)
	pp <- vector("list", 3)
	pp.tree <- read.tree(paste(path, "/posteriorAncestors.tree", sep=""))
	if (all(nodes > length(pp.tree$tip.label)))
		nodes <- nodes-length(pp.tree$tip.label)

	for (i in nodes){
		#PP 
		pp.tree$node.label[i] <- gsub(paste("[[]N", as.character(i), "_", sep=""), "",pp.tree$node.label[i])
		pp.tree$node.label[i] <- gsub("[]]", "",pp.tree$node.label[i])
		pp[[1]] <- c(pp[[1]], gsub("//.*", "", pp.tree$node.label[i]))
		pp[[2]] <- c(pp[[2]], gsub(".*//", "", pp.tree$node.label[i]))
		#ML
		pp[[3]] <- as.numeric(cE[[3]]$node.label[i])
	}
	pp <- as.data.frame(pp)
	dimnames(pp)[[1]] <- as.character(nodes)
	dimnames(pp)[[2]] <- c("Best CN - PP", "2nd best CN - PP", "CN-ML")
	pp
}

###
#give you all chromosome numbers inferred for a certain node with the respective posterior probabilites if it is greater than cutoff.
#node numbers can be according to chromEvol or phylo node numbering

PP.nodes <- function(ChromEvol, nodes, cutoff=0.1, percentage = 0.8){
	pp <- round(ChromEvol[[1]][[1]], digits=3)
	if (all(nodes > length(pp[,1])))
		nodes <- nodes-length(ChromEvol[[3]]$tip.label)
	xx <- vector("list",length(nodes))
	for(i in 1:length(nodes)){
		p <- pp[nodes[i],which(pp[nodes[i],] > cutoff)]	
		while (sum(p) < percentage){
			cutoff <- cutoff-0.01
			p <- pp[nodes[i],which(pp[nodes[i],] > cutoff)]
		}
		names(p) <- names(which(pp[nodes[i],] > cutoff))
		xx[[i]] <- c(round(p, digits=2), SUM = sum(p), ML=round(as.integer(ChromEvol[[3]]$node.label[nodes[i]]), digits=0))
		names(xx) <- paste(nodes, " - ", nodes+length(ChromEvol[[3]]$tip.label))

}	
	xx
	}
####################################
sum.Events <- function(ChromEvol){
		ChromEvol.exp <- ChromEvol[[2]]
		exp.table <- ChromEvol.exp[[2]]
		EVENTS.tab <- ChromEvol.exp[[1]]
		
		ALL <- NULL
		sm05 <- NULL
		big05 <- NULL
		NoTable <-NULL
		HE <- NULL
		for (i in 1:4){
		#####
			#exp.table[which(exp.table[,i+1] < 0.5), i+1] <- NA
			ALL <- c(ALL, sum(round(exp.table[,i+1], digits=1),na.rm=TRUE))
			sm05 <- c(sm05, sum(round(exp.table[which(exp.table[,i+1] <= 0.5), i+1], digits=1)))
			big05 <- c(big05, sum(round(exp.table[which(exp.table[,i+1] > 0.5), i+1], digits=1)))
			NoTable <- c(NoTable, sum(round(as.numeric(EVENTS.tab[[i]][-length(EVENTS.tab[[i]][,1]), 2]), digits=1)))
			HE <- c(HE, sum(ChromEvol.exp[[3]][,i]))
	}
	data.frame(All=c(ALL, sum(ALL)), SmallerThan0.5=c(sm05, sum(sm05)), BiggerThan0.5=c(big05, sum(big05)), OnlySimBigger0.5=c(NoTable, sum(NoTable)), OnlyHeur=c(HE, sum(HE)), row.names=c("Gains", "Losses", "Dupl.", "Demi.", "Sum"))
}

#######################################
#Extracts a clade from a phylogeny in a ChromEvol object, together with all node, branch and tip information. 
#NOTE : in the events table columns "OnlySimBigger05" and "OnlyHeur" are not adapted to the extracted clade and should therefore not be considered (information still from the original, complete phylogeny)

extract.clade.ce <- function(ChromEvol, node){
	ceN <- ChromEvol
	nodes <- c(node, descendants(ChromEvol[[3]], node, type="i"))
	N <- nodes - length(ChromEvol[[3]]$tip.label)
	E <- which.edge(ChromEvol[[3]], descendants(ChromEvol[[3]], node))
	T <- descendants(ChromEvol[[3]], node)

	ceN[[1]][[1]] <- ChromEvol[[1]][[1]][N,]
	ceN[[1]][[2]] <- ChromEvol[[1]][[2]][T,]
	ceN[[1]][[3]] <- ChromEvol[[1]][[3]][N]

	ceN[[2]][[2]] <- ChromEvol[[2]][[2]][E,]

	ceN[[3]] <- extract.clade(ChromEvol[[3]], node)
ceN
}
#############################################
collapse.clade.ce <- function(ChromEvol, node, name="collapsed", numbers=""){
	
	ceN <- ChromEvol
	#find all tips descending from this node
	tips <- descendants(ChromEvol[[3]],node)
	#T <- tips[-length(tips)]
	T <- tips[-1]
	#find all nodes descending from node
	nodes <- c(node,descendants(ChromEvol[[3]],node, type="i"))
	
	N <- nodes - length(ChromEvol[[3]]$tip.label)
	E <- which.edge(ChromEvol[[3]], T)
	E <- c(E[1]-1, E)

	ceN[[1]][[1]] <- ChromEvol[[1]][[1]][-N,]
	ceN[[1]][[2]] <- ChromEvol[[1]][[2]][-T,]
	ceN[[1]][[3]] <- ChromEvol[[1]][[3]][-N]

	ceN[[2]][[2]] <- ChromEvol[[2]][[2]][-E,]
	if (numbers != "")
		numbers <- paste("-",numbers, sep="")
	name <- paste(name, " (", length(tips), " tips)", sep="")
	ChromEvol[[3]]$tip.label[tips[1]] <- name
	ceN[[3]] <- drop.tip(ChromEvol[[3]], T)
	
	ceN$collapse <- rbind(ceN$collapse, ChromEvol[[1]][[1]][node-length(ChromEvol[[3]]$tip.label),])
	dimnames(ceN$collapse)[[1]][length(ceN$collapse[,1])] <- name
	ceN
}
##############
### Use this function to order information about species according to the tips of a tree for plotting purposes. The first column of 'table' or rownames must contain tip labels as used in 'tree', the following column the information to be plotted
#Author: N.Cusimano, 26.05.07 
match.table2tip<-function(tree, table, nomatch=NA)
{	if (dim(table)[2] == 1)	
		table <- cbind(dimnames(table)[[1]], table[,1])
	d <- duplicated(c(tree$tip.label, table[,1]))
	dd <- d[(length(tree$tip.label)+1):length(d)]
		table <- table[dd,]
	if (length(tree$tip.label) > length(table[[1]]))
	{	table2<-NULL
		for (i in 1:length(table[1,]))
			table2 <- cbind(table2,table[,i][match(tree$tip.label, table[,1], nomatch=nomatch)])
		
		}
	else
	{
		m<-tree$tip.label[match(table[[1]], tree$tip.label, nomatch=NA)]
		
		table1<-data.frame(table[[1]][!is.na(m)])
		
		for (i in 2:length(table))
			table1[i]<-table[[i]][!is.na(m)]
			
		table2<-table1
		
		for (i in 1:length(table1))
			table2[i]<-table1[[i]][match(tree$tip.label, table1[[1]], nomatch=NA)]
	}
	names(table2)<-names(table)
	return(table2)
	}
	

	
####calculate tree length
	tree.length <- function(phy){
		xx <- numeric(length(phy$tip.label) + phy$Nnode)
    	for (i in 1:dim(phy$edge)[1]) xx[phy$edge[i, 2]] <- 
    		xx[phy$edge[i,1]] + phy$edge.length[i]
    	x <- max(xx[1:length(phy$tip.label)])
    	x
    	}
    	
 ###make table with chromosome numbers
 
 make.input.table <- function(phy, file="Chromosome_numbers.txt", sort=TRUE){
 	if (sort)
	 	
	 	write.table(cbind(sort(phy$tip.label), rep(NA, length(phy$tip.label))), quote=FALSE, col.names=FALSE, sep="\t",row.names=FALSE, file=file)
	 else
	 	write.table(cbind(phy$tip.label, rep(NA, length(phy$tip.label))), quote=FALSE, col.names=FALSE, sep="\t",row.names=FALSE, file=file)
 }
 
 ############## read. ChromEvol table ###########
# table: path of the chromEvol input table
 
 read.cE.tab <- function(table, geiger=FALSE){
 	table <- read.table(table, sep="\n", stringsAsFactor=FALSE)
 	table2 <- data.frame()
 	table2 <- table[grep(">", table[,1]),1]
 	table2 <- gsub(">", "", table2)
 	table2 <- cbind(table2, table[grep(">", table[,1], invert=TRUE),1])
 	table2 <- gsub(" ", "", table2)
 	table2 <- as.data.frame(table2, stringsAsFactors=FALSE)
 	
 	if (geiger){
 		table2[,2] <- gsub("=.*", "", table2[,2])
 		table3 <- as.factor(table2[,2])	
 		names(table3) <- table2[,1]
 		table2 <- table3
 		}
 	table2
 	}