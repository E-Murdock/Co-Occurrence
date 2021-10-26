#pull nodes and edges

##setwd("YOUR FILE LOCATION")##

nodes = read.csv("./data/synthase 5 by 70 node.csv")
edges = read.csv("./data/synthase 5 by 70 edge.csv")
edges = edges$name
edges = gsub(","," ", edges)
edges = strsplit(edges, " ")
edges = data.frame(matrix(unlist(edges), nrow=length(edges), byrow=TRUE))

#list of connections
connect = vector(mode = "list", length = nrow(nodes))
for (i in 1:nrow(nodes)) {
  x = edges[which(edges$X1 %in% nodes[i,2]),2]
  y = edges[which(edges$X2 %in% nodes[i,2]),1]
  z = c(nodes[i,2],x,y)
  connect[[i]] = z
}

#pull singlets and cluster members
singlets = nodes[which(lengths(connect) < 2),]
clustered = nodes[which(lengths(connect) > 1),]

#process singlets
singlets = sub(" .*$", "", singlets$Description)

#process clustered
for (i in 1:nrow(clustered)) {
  clustered[i,1]= sub(" .*$", "", clustered[i,1])
}
remove = rev(which(lengths(connect) < 2))
for (i in remove) {
  connect[[i]] = NULL
}

#create bins
for (i in 1:length(connect)) {
  for (j in 1:length(connect)){
    if (connect[[i]][1] != 0) {
      if (connect[[j]][1] !=0) {
        x = is.element(connect[[i]],connect[[j]])
        if (is.element("TRUE", x) == "TRUE") {
          connect[[i]] = c(connect[[i]],connect[[j]])
          connect[[i]] = unique(connect[[i]])
          if (i != j) {
            connect[[j]] = 0
          }
        }
      }
    }
  }
}

#isolating bins, get IDs
SIDs = vector(mode = "list", length = 0)
for (i in which(lengths(connect) > 1)) {
  SIDs = append(SIDs, connect[i])
}
for (i in 1:length(SIDs)) {
  SIDs[[i]] = clustered[which(clustered$name %in% SIDs[[i]]),1]
  SIDs[[i]] = strsplit(SIDs[[i]], " ")
  SIDs[[i]] = unlist(SIDs[[i]])
  SIDs[[i]] = unique(SIDs[[i]])
}
SIDs = c(SIDs,singlets)

#synthase geneID dataframe
for (i in 1:length(SIDs)) {
  length(SIDs[[i]])=max(lengths(SIDs))
}
SIDs = data.frame(matrix(unlist(SIDs), nrow=length(SIDs), byrow=TRUE))

#gene to genome
Sgenetogenome = read.csv("./data/Sgene to genome.csv")
SIDs = as.data.frame(t(SIDs))
for (i in 1:length(SIDs)) {
  Scells=which(Sgenetogenome$Gene.ID %in% SIDs[,i])
  length(Scells)=nrow(SIDs)
  SIDs[,i]=Sgenetogenome[Scells,4]
}
Sgroups=paste("synthase.bin",1:length(SIDs),sep=".")
colnames(SIDs)=Sgroups

#pull nodes and edges
nodes = read.csv("./data/5 by 85 node.csv")
edges = read.csv("./data/5 by 85 edge.csv")
edges = edges$name
edges = gsub(","," ", edges)
edges = strsplit(edges, " ")
edges = data.frame(matrix(unlist(edges), nrow=length(edges), byrow=TRUE))

#list of connections
connect = vector(mode = "list", length = nrow(nodes))
for (i in 1:nrow(nodes)) {
  x = edges[which(edges$X1 %in% nodes[i,3]),2]
  y = edges[which(edges$X2 %in% nodes[i,3]),1]
  z = c(nodes[i,3],x,y)
  connect[[i]] = z
}

#pull singlets and cluster members
singlets = nodes[which(lengths(connect) < 2),]
repnodes = singlets[which(singlets$Number.of.IDs.in.Rep.Node > 1),]
singlets = singlets[which(singlets$Number.of.IDs.in.Rep.Node < 2),]
clustered = which(lengths(connect) > 1)
clustered = nodes[clustered,]

#process singlets
singlets = sub(" .*$", "", singlets$Description)

#process repnodes
repIDs = vector(mode = "list", length = nrow(repnodes))
for (i in 1:nrow(repnodes)) {
  repnodes[i,1]=gsub("\\|"," ", repnodes[i,1])
  repnodes[i,1]=gsub("(\\S+[a-zA-Z]\\S+|[^0-9])"," ",repnodes[i,1])
  repnodes[i,1]=gsub("\\s+"," ",repnodes[i,1])
  repnodes[i,1]=trimws(repnodes[i,1])
  repIDs[[i]] = repnodes[i,1]
  repIDs[[i]]=strsplit(repIDs[[i]]," ")
  repIDs[[i]]=repIDs[[i]][[1]]
}

#process clustered
for (i in 1:nrow(clustered)) {
  clustered[i,1]=gsub("\\|"," ", clustered[i,1])
  clustered[i,1]=gsub("(\\S+[a-zA-Z]\\S+|[^0-9])"," ",clustered[i,1])
  clustered[i,1]=gsub("\\s+"," ",clustered[i,1])
  clustered[i,1]=trimws(clustered[i,1])
}
remove = rev(which(lengths(connect) < 2))
for (i in remove) {
  connect[[i]] = NULL
}

#create bins
for (i in 1:length(connect)) {
  for (j in 1:length(connect)){
    if (connect[[i]][1] != 0) {
      if (connect[[j]][1] !=0) {
        x = is.element(connect[[i]],connect[[j]])
        if (is.element("TRUE", x) == "TRUE") {
          connect[[i]] = c(connect[[i]],connect[[j]])
          connect[[i]] = unique(connect[[i]])
          if (i != j) {
            connect[[j]] = 0
          }
        }
      }
    }
  }
}

#isolating bins, get IDs
clusters = vector(mode = "list", length = 0)
for (i in which(lengths(connect) > 1)) {
  clusters = append(clusters, connect[i])
}
for (i in 1:length(clusters)) {
  clusters[[i]] = clustered[which(clustered$name %in% clusters[[i]]),1]
  clusters[[i]] = strsplit(clusters[[i]], " ")
  clusters[[i]] = unlist(clusters[[i]])
  clusters[[i]] = unique(clusters[[i]])
}
clusters = c(clusters,repIDs,singlets)

#Ligase geneID dataframe
for (i in 1:length(clusters)) {
  length(clusters[[i]])=max(lengths(clusters))
}
clusters = data.frame(matrix(unlist(clusters), nrow=length(clusters), byrow=TRUE))

#convert to genome ID
Lgenetogenome = read.csv("./data/Lgene to genome.csv")
clusters = as.data.frame(t(clusters))
for (i in 1:length(clusters)) {
  Lcells=which(Lgenetogenome$Gene.ID %in% clusters[,i])
  length(Lcells)=nrow(clusters)
  clusters[,i]=Lgenetogenome[Lcells,5]
}
Lgroups=paste("ligase.bin",1:length(clusters),sep=".")
colnames(clusters)=Lgroups

#combining
TotalIDs=c(clusters,SIDs)
for (i in 1:length(TotalIDs)) {
  TotalIDs[[i]]=TotalIDs[[i]][!is.na(TotalIDs[[i]])]
}

#presence absence matrix
GenomeIDs=read.csv("./data/GenomeIDs.csv")
PAM=data.frame(matrix(ncol=nrow(GenomeIDs),nrow=length(TotalIDs)))
colnames(PAM)=GenomeIDs$IMG.Genome.ID
rownames(PAM)=names(TotalIDs)
for (i in 1:nrow(PAM)) {
  Pcells=which(GenomeIDs[,7] %in% TotalIDs[[i]])
  PAM[i,Pcells]=1
}
PAM[is.na(PAM)] = 0

##if jaccard is unnecessary, and only sufficiency matters use this option
##this will still take several hours
#sufficiency measure
sufficiency = vector(mode = "list", length  = 0)
for (i in 1:ncol(clusters)) {
  for (j in 1:ncol(SIDs)) {
    ligase.presence = which(PAM[i,] == 1)
    synthase.presence = which(PAM[j+ncol(clusters),] == 1)
    co.occurrence = which(synthase.presence %in% ligase.presence)
    synthase.in.ligase = length(co.occurrence)
    synthase.out.ligase = sum(PAM[j+ncol(clusters),]) - synthase.in.ligase
    synthase.percent = synthase.in.ligase/sum(PAM[j+ncol(clusters),])*100
    sufficiency = c(sufficiency,list(c(Lgroups[i],Sgroups[j],sum(PAM[i,]),
                                 sum(PAM[j+ncol(clusters),]),synthase.in.ligase,
                                 synthase.out.ligase,synthase.percent)))
  }
}
sufficiency = data.frame(matrix(unlist(sufficiency),
                            nrow=(ncol(clusters)*ncol(SIDs)), byrow=TRUE))
colnames(sufficiency) = c("ligase.bin","synthase.bin","ligase occurence",
                       "synthase occurence","synthase in", "synthase out",
                       "synthase percent")

##running jaccards takes a long time, 5 by 85 jaccard data is in the github already
##if you want to run it on different nodes and edges use this option

##only need to run this once
#install.packages("jaccard")

#jaccard w/ ratios
#library(jaccard)
#jaccards = vector(mode = "list", length  = 0)
#for (i in 1:ncol(clusters)) {
#  for (j in 1:ncol(SIDs)) {
#    j.value = jaccard(as.numeric(PAM[i,]),as.numeric(PAM[j+ncol(clusters),]))
#    j.sig = jaccard.test(as.numeric(PAM[i,]),as.numeric(PAM[j+ncol(clusters),]),
#                         method = "exact")
#    ligase.presence = which(PAM[i,] == 1)
#    synthase.presence = which(PAM[j+ncol(clusters),] == 1)
#    co.occurrence = which(synthase.presence %in% ligase.presence)
#    synthase.in.ligase = length(co.occurrence)
#    synthase.out.ligase = sum(PAM[j+ncol(clusters),]) - synthase.in.ligase
#    synthase.percent = synthase.in.ligase/sum(PAM[j+ncol(clusters),])*100
#    jaccards = c(jaccards,list(c(Lgroups[i],Sgroups[j],j.value,sum(PAM[i,]),
#                                 sum(PAM[j+ncol(clusters),]),synthase.in.ligase,
#                                 synthase.out.ligase,synthase.percent,
#                                 j.sig[[1]][1],j.sig[[2]][1],j.sig[[3]][1])))
#  }
#}
#jaccards = data.frame(matrix(unlist(jaccards),
#                             nrow=(ncol(clusters)*ncol(SIDs)), byrow=TRUE))
#colnames(jaccards) = c("ligase.bin","synthase.bin","jaccard","ligase occurence",
#                       "synthase occurence","synthase in", "synthase out",
#                       "synthase percent","centered","p.value","expectation")
