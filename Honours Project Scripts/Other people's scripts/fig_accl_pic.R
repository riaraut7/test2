# Let's do the acclimation vs. adaptation figure again using PICs

# Make tree
library(stringr)
library(V.PhyloMaker)
library(ggtree)


# Generate a phylogeny and plot
traits_sub = traits %>% subset(!is.na(bio10))

traits_sub %>% 
  select(species_binomial_TNRS, Family = family) %>% #making two new columns; binomial name and family. #oh wait you already have both? 
  unique() ->a #storing the unique species into a dataframe called 'a'

a$species_binomial_TNRS %>%  #taking the binomial name and splitting it into genus and species 
  str_split_fixed(" ", 2) %>% 
  as.data.frame() %>% 
  rename(Genus = V1, Species = V2) %>% #naming your new split columns appropriately 
  bind_cols(a) -> b #storing this dataset as 'b' 

b %>% rename(genus = Genus, species = species_binomial_TNRS, family = Family) %>% 
  select(-Species) -> b #from b, take out everything (curves?) that don't have? an appropriate species name? 

b = b %>% select(species, genus, family) #only select the species, genus, family column 

tree <- phylo.maker(b, tree = GBOTB.extended, output.tree = T, r = 1) #make a tree 

tree3 <- tree$scenario.3

#from the same dataset, you're selecting your columns of interest, grouping by species, and getting average 
#traits of interest (ex. r_max) : 
traits_sub %>% 
  select(species_binomial_TNRS,topt,e,eh,r_max,breadth, bio10, Rolling_Mean_C) %>% 
  group_by(species_binomial_TNRS) %>% 
  summarize(topt = mean(topt), 
            e = mean(e), 
            eh = mean(eh), 
            r_max = mean(r_max),
            breadth = mean(breadth),
            Thome = mean(bio10),
            Tgrowth = mean(Rolling_Mean_C)) %>% 
  as.data.frame() -> x #save this dataframe as x 
rownames(x) = x$species_binomial_TNRS
x = x %>% select(-species_binomial_TNRS)
rownames(x) = rownames(x) %>% gsub(" ", "_",.) #replacing spaces with underscores For consisteny's sake, I think g
x = as.matrix(x)
##* now Josef has a dataset such that there's a column for species name, and another column 
##* for each average value of interest (in your case, this would have to be mean(Amax), etc )

# Next, how to do a PIC regression?
#library(ape)
#cent.tree<-read.tree("Centrarchidae.tre")
tree3
#buccal.length<-setNames(obj[,"buccal.length"],
#                        rownames(obj))
topt = setNames(x[,"topt"],rownames(x))
ea = setNames(x[,"e"], rownames(x))
ed = setNames(x[,"eh"], rownames(x))
r_max = setNames(x[,"r_max"], rownames(x))
breadth = setNames(x[,"breadth"], rownames(x))
Thome = setNames(x[,"Thome"], rownames(x))
Tgrowth = setNames(x[,"Tgrowth"], rownames(x))

pic.topt<-pic(topt,tree3) # Note that this automatically scales variables; can be turned off if desired
pic.ea<-pic(ea,tree3)
pic.ed<-pic(ed,tree3)
pic.r_max<-pic(r_max,tree3)
pic.breadth<-pic(breadth,tree3)
pic.Tahome<-pic(Thome,tree3)
pic.Tgrowth<-pic(Tgrowth,tree3)

df = data.frame(pic.topt,pic.ea,pic.ed,pic.r_max,pic.breadth,pic.Tgrowth,pic.Tahome)
df = df %>% 
  pivot_longer(cols = 1:5, names_to="parameter",values_to="param_value") %>% 
  pivot_longer(cols = c("pic.Tgrowth","pic.Tahome"), names_to = "T_type", values_to = "T_value")

models = df %>% group_by(parameter, T_type) %>% do(model = lm(param_value ~ T_value, data = .))
summary(models$model[[1]])

significance = data.frame(parameter = rep(c("pic.topt","pic.r_max","pic.ea","pic.ed","pic.breadth"),2),
                          T_type = c(rep("pic.Tahome",5),rep("pic.Tgrowth",5)),
                          sig = as.factor(c(1,0,0,0,0,1,1,1,0,0)))

df = merge(df, significance, by=c("parameter","T_type"))


p_all = ggplot(df, aes(x = T_value, y = param_value)) +
  #geom_abline(slope=1,lty=2) +
  geom_point(size = 1, fill = "lightblue",pch=21) +
  geom_smooth(method = "lm", formula = y~x+0,color="black", aes(lty=sig),linewidth=0.75) +
  scale_linetype_manual(values = c(2,1), guide="none") +
  my_theme +
  xlab(NULL) +
  ylab(NULL) +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.direction = "vertical")+
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  #guides(color = guide_legend( 
  #  override.aes=list(shape = 21))) +
  facet_grid(rows=vars(parameter), 
             cols=vars(T_type), 
             scales = "free",
             switch = "both",
             labeller = as_labeller(c(pic.breadth = "PIC~~Omega",#~(phantom()^o*C)",
                                      pic.Tahome = "PIC~~T[home]",#~(phantom()^o*C)",
                                      pic.Tgrowth = "PIC~~T[growth]",#~(phantom()^o*C)",
                                      pic.ea = "PIC~~E[A]",#~(eV)",
                                      pic.ed = "PIC~~E[D]",#~(eV)",
                                      pic.r_max = "PIC~~A[max]",#~(μmol~m^-2~s^-1)",
                                      pic.topt = "PIC~~T[opt]"),#~(phantom()^o*C)"),
                                    default = label_parsed))


svg("figures/fig_accl_adap_pic.svg", width=4, height = 6.5)
p_all
dev.off()
