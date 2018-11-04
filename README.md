# Cat predation on vertebrates

*Marco Girardello* (marco.girardello@gmail.com) 

R code developed as part of the manuscript "Dog days are over...cat ones still ongoing! Domestic cats are a most serious threat to native biodiversity in a Mediterranean country "Mori et al. (in review) . A copy of the main script (AnalysesCats.R) is shown below. 

```r
#########################
# Cat predation analyses 
########################

# load required libraries
library(tidyverse);library(stringi)
library(FD);library(ape)
library(magrittr);library(data.table)
library(network);library(sna)
library(ggnetwork)

# load data files
cats<-read.csv("/mnt/data1tb/Dropbox/gatti/data/gatti/Gatti.csv")

# trait datasets
mamtraits<-read.csv("/mnt/data1tb/Dropbox/HistFunc/traitdata/wilman14/traitfinalbats.csv")
birdtraits<-read.csv("/mnt/data1tb/Dropbox/FDTD/birdsFDfiles/traitbirdsALL.csv")
othertraits<-read.csv("/mnt/data1tb/Dropbox/gatti/data/traits/othergroups.csv")

####################
# Data preparation
####################

# prepare trait data for analyses
mamtraits %>%
  dplyr::select(binomial,Diet.Inv,Diet.Vend,Diet.Vect,Diet.Vfish,Diet.Vunk,Diet.Scav,Diet.Fruit,
                Diet.Nect,Diet.Seed,Diet.PlantO,BodyMass.Value) %>%
  rename(Species=binomial) ->mamtraits1

birdtraits %>%
  dplyr::select(species,Diet.Inv,Diet.Vend,Diet.Vect,Diet.Vfish,Diet.Vunk,Diet.Scav,Diet.Fruit,
                Diet.Nect,Diet.Seed,Diet.PlantO,BodyMass.Value) %>%
  rename(Species=species) %>%
  bind_rows(mamtraits1) ->intermd

othertraits %>%
  dplyr::select(-c(Classe))  %>%
  mutate(Species=paste(stri_split_fixed(Species," ",simplify=T)[,1],stri_split_fixed(Species," ",simplify=T)[,2],sep="_")) %>%
  bind_rows(intermd) ->alltraits

write.csv(alltraits,file="/mnt/data1tb/Dropbox/gatti/data/traits/alltraits.csv",row.names=FALSE)

# prepare cat data for analyses
cats %>%
  # trim white spaces
  mutate(Species=trimws(Species),Classe=trimws(Classe),TOPONIMO=trimws(TOPONIMO)) %>%
  # change Squamata for Reptilia
  mutate(Classe=ifelse(Classe=="Squamata","Reptilia",Classe)) %>%
  # exclude Suncus etruscus ALBINO and Gallus gallus
  filter(Species!="Suncus etruscus ALBINO" & Species!="Gallus gallus") %>%
  # change names to match the style of the trait dataframes
  mutate(Species=paste(stri_split_fixed(Species," ",simplify=T)[,1],stri_split_fixed(Species," ",simplify=T)[,2],sep="_"))  %>%
  # fix a few names 
  filter(Species!="Erithachus_rubecula") %>%
  mutate(Species=ifelse(Species=="Passer_italiae","Passer_domesticus",Species),
         Species=ifelse(Species=="Sylvia_subalpina","Sylvia_cantillans",Species),
         Species=ifelse(Species=="Sylvia_subalpina","Sylvia_cantillans",Species)) %>%
  dplyr::select(TOPONIMO,Species,Classe,Italian.IUCN,International.IUCN) %>%
  # extract unique site and species combinations
  group_by(TOPONIMO) %>%
  distinct(Species,Classe,TOPONIMO,Italian.IUCN,International.IUCN) %>%
  as.data.frame() ->cats1

# combine cat dataset with trait information
cats1 %>%
  # extract species column only
  dplyr::select(Species) %>%
  # unique list of species
  distinct(Species) %>%
  # merge with mammal traits
  left_join(alltraits) %>%
  # exclude NAs
  filter(!is.na(Diet.Inv))  %>%
  # take the log of body mass
  mutate(BodyMass.Value=log(BodyMass.Value)) ->traitdf

# Missing species from TRAIT ANALYSIS (FIX!!!)
anti_join(edge_df1 %>% dplyr::select(group2) %>% 
            rename(Species=group2) %>%
            mutate(Species=stri_replace_all_fixed(Species," ","_")),traitdf %>%
            dplyr::select(Species)) 

####################
# Data analysis
####################

#--- Barchart of impacts group according to IUCN threat categories
cats1 %>%
  filter(Italian.IUCN!="ND" & Classe!="" & Italian.IUCN!="INVALID SPECIES" & International.IUCN!="ND") %>%
  dplyr::select(-c(TOPONIMO,Species))   %>%
  gather(IUCN,iucn.cat,-c(Classe))  %>%
  mutate(IUCN=ifelse(IUCN=="International.IUCN","International IUCN","Italian IUCN")) %>%
  group_by(Classe,IUCN,iucn.cat) %>%
  summarise(count=n()) %>%
  group_by(IUCN) %>%
  mutate(freq = count / sum(count)) %>%
  mutate(iucn.cat=factor(iucn.cat,levels=c("VU","NT","LC","DD","NE","Alien"))) %>%
  ggplot(.,aes(x=Classe,y=freq,fill=iucn.cat))+geom_bar(stat="identity",colour="black")+
  theme_bw()+scale_fill_brewer(palette='RdBu')+
  theme(axis.text.x = element_text(size=17,color="black",angle=40,hjust=1),
        axis.text.y=element_text(size=17,color="black"),
        strip.text =element_text(size=19,face="bold"),axis.title =element_text(size=20),
        legend.text = element_text(size=15),legend.title=element_text(size=15))+
  facet_wrap(~IUCN)+xlab("Class")+ylab("Proportion of species killed")+
  labs(fill="IUCN Category") ->barcat

ggsave(barcat,filename="/mnt/data1tb/Dropbox/gatti/figures/barcat.png",width=12,height=9,dpi=400)

# --- Impacts on the functional structure of vertebrate communities

# Principal coordinate analysis
row.names(traitdf)<-as.character(traitdf$Species)
pcoares<-ape::pcoa(gowdis(traitdf[2:12]),correction="cailliez")

# insert information of animal Class
data.frame(pcoares$vectors[,1:2],Species=row.names(pcoares$vectors)) %>% 
  left_join(cats1 %>% 
              dplyr::select(Species,Classe)  %>%
              distinct(Species,Classe)) %$%
  as.data.table(.) ->results

# calculate convex hulls and create pretty figure
results[, .SD[chull(Axis.1,Axis.2)], by =Classe] ->hulls  

# create plot
ggplot(data=results,aes(x=Axis.1,y=Axis.2,colour=Classe))+geom_point(size=3)+theme_bw()+
  geom_polygon(data = hulls,aes(fill=Classe),alpha = 0.2)+
  theme(axis.text =element_text(size=19,colour="black"),axis.title = 
  element_text(size=19,colour="black"),legend.title=element_blank(),legend.text=element_text(size=19))+
  ylab("Dim 2")+xlab("Dim 1")->finalp  

# save plot
ggsave(finalp,file="/mnt/data1tb/Dropbox/gatti/figures/Fvolumes.png",width=9,height=8,
       dpi=400)

# --- Multivariate analysis of variance
res.man <- manova(cbind(Axis.1, Axis.2) ~Classe, data = results1)
sink("/mnt/data1tb/Dropbox/gatti/results/manova.txt")
summary(res.man)
sink()

# --- Create network interaction matrix
cats1 %>%
  mutate(counter=1) %>%
  group_by(Classe,Species) %>%
  summarise(connection_strength=sum(counter)) %>%
  filter(Classe!="") %>%
  mutate(group1="Felis catus") %>%
  group_by(Classe)  %>%
  mutate(connection_strength=connection_strength/sum(connection_strength)) %>%
  rename(group2=Species) %>%
  as.data.frame() %>%
  filter(group2!="Sorex_sp." & group2!="Pipistrellus_sp.") %>%
  mutate(group2=stri_replace_all_fixed(group2,"_"," ")) ->edge_df

# --- Separate network plot for each animal class

classe.l<-unique(edge_df$Classe)
res<-NULL

for (i in 1:length(classe.l)){
  tmp<-subset(edge_df,Classe==classe.l[i])
  # create network
  n<-network(tmp[,c(2,4)],directed = FALSE)
  # insert node labels
  n %v% "species" <- tmp$group2
  # connection strength
  set.edge.attribute(n, "connection_strength",tmp$connection_strength)
  # save into a dataframe format
  tmp.net.df<-ggnetwork(n)
  tmp.net.df$Classe<-unique(tmp$Classe)
  res<-rbind(tmp.net.df,res)
  
}

res$Classe<-factor(res$Classe,levels=c("Amphibia","Aves","Mammalia","Reptilia"))

# create plot
netcat<-ggplot(data=res, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black",aes(size=connection_strength))+
  geom_nodes(color = "grey", size = 8)+geom_nodelabel_repel(aes(label = vertex.names),
  fontface = "italic", box.padding = unit(1, "lines"),size=3)+theme_bw()+
  theme(legend.position="none")+facet_wrap(~Classe)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"),
  strip.text=element_text(face="bold",size=16),axis.text = element_blank())+
  ylab("")+xlab("")

# save plot
ggsave(netcat,filename = "/mnt/data1tb/Dropbox/gatti/figures/network.png",
       width=11,height=10,dpi=400)

#---- Network plots for all classes pooled together

cats1 %>%
  filter(Species!="_") %>%
  mutate(counter=1) %>%
  group_by(Species) %>%
  summarise(connection_strength=sum(counter))  %>%
  mutate(group1="Felis catus") %>%
  mutate(connection_strength=connection_strength/sum(connection_strength)) %>%
  rename(group2=Species) %>%
  as.data.frame() %>%
  filter(group2!="Sorex_sp." & group2!="Pipistrellus_sp.") %>%
  mutate(group2=stri_replace_all_fixed(group2,"_"," ")) ->edge_df1

# create network
n<-network(edge_df1[,c(1,3)],directed = FALSE)
# insert node labels
n %v% "species" <- tmp$group2
# connection strength
set.edge.attribute(n, "connection_strength",edge_df1$connection_strength)
# save into a dataframe format
res1<-ggnetwork(n)

# create network plot
netcat1<-ggplot(data=res1, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black",aes(size=connection_strength))+
  geom_nodes(color = "grey", size = 8)+geom_nodelabel_repel(aes(label = vertex.names),
  fontface = "italic", box.padding = unit(1, "lines"),size=3)+theme_bw()+
  theme(legend.position="none")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"),
  strip.text=element_text(face="bold",size=16),axis.text = element_blank())+
  ylab("")+xlab("")

ggsave(netcat1,filename = "/mnt/data1tb/Dropbox/gatti/figures/network1.png",
       width=10,height=10,dpi=400)
```