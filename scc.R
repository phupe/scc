#######################################################################################
# This file is part of geniac.
#
# Copyright Institut Curie 2023.
#
# This software is a computer program whose purpose is to perform
# statistics on squamous cell carcinoma data
#
# You can use, modify and/ or redistribute the software under the terms
# of license (see the LICENSE file for more details).
#
# The software is distributed in the hope that it will be useful,
# but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
# Users are therefore encouraged to test the software's suitability as regards
# their requirements in conditions enabling the security of their systems and/or data.
#
# The fact that you are presently reading this means that you have had knowledge
# of the license and that you accept its terms.
#######################################################################################


require(ggplot2)
require(scales)


### GLOBOCAN 2020 ##
### Number of new cancer cases and death per year
### data downloaded from https://gco.iarc.fr/tomorrow in march 2023
### source: https://gco.iarc.fr/tomorrow
### reference: Sung H, et al., Global cancer statistics 2020: GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers in 185 countries. CA Cancer J Clin. 2021

dir.in <- 'globocan'

### Organs with squamous cell carcinoma (SCC)
head.neck.organs <-  c('Hypopharynx', 'Larynx', 'Lip, oral cavity', 'Oropharynx', 'Nasopharynx')
hn <- "Head and neck"
lung <- "Lung"
trachea <- "Trachea, bronchus and lung"
squamous.organs <- c('Anus', 'Cervix uteri', hn, 'Oesophagus', 'Penis', trachea, 'Vulva')

### Percentage of squamous cancer for each organ (with references)
squamous.infos <- read.table(paste(dir.in, 'scc-infos.tsv', sep='/'), sep='\t', header=TRUE, stringsAsFactors=FALSE)
years <- c(2020, 2040)

### STEP 1 - Read all the GLOBOCAN data for World, European union and Argentina
### extract the information for the SCCs

all.cancer <- NULL
for (f in list.files(dir.in, pattern='*.csv')) {
	zone <- gsub('.*-', '', f)
	zone <- gsub('\\.csv', '', zone)
	indice <- gsub('.*2040-', '', f)
	indice <- gsub('-.*', '', indice)
	values <- read.csv(paste(dir.in, f, sep='/'), header=TRUE, stringsAsFactors=FALSE)
	values$zone <- zone
	ind <- which(values$Cancer.label == 'All cancers' & (values$Year == 2020 | values$Year == 2040))
	info.all.cancer <- values[ind, c('Cancer.label', 'Year', 'Prediction')]
    info.all.cancer$zone <- zone 
	if(indice == 'mort') {
		values$indice <- 'Death'
    	info.all.cancer$indice <- 'Death'
	}
	if(indice == 'inc') {
		values$indice <- 'New'
    	info.all.cancer$indice <- 'New'
	}
	if(zone == 'world') {
    	info.all.cancer$zone <- 'World'
	}
	if(zone == 'eu') {
    	info.all.cancer$zone <- 'European Union'
	}
	if(zone == 'argentina') {
    	info.all.cancer$zone <- 'Argentina'
	}
	all.cancer <- rbind(all.cancer, info.all.cancer)
	ind.hn <- which(values$Cancer.label %in% head.neck.organs)
    values$Cancer.label[ind.hn] <- hn
	ind.years <- which(values$Year %in% years)
	values <- values[ind.years,]
	ind.squamous <- which(values$Cancer.label %in% squamous.organs)
	values <- values[ind.squamous,]
	values$number <- values$Change.in.number.of.cases
	ind.2020 <- which(values$Change.in.number.of.cases == 0)
	values$number[ind.2020] <- values$number[ind.2020] + values$Prediction[ind.2020]
	values <- values[, c('Population', 'Cancer.label', 'Year', 'zone', 'indice', 'number')]
	values$Year <- as.factor(values$Year)
	values$Year <- factor(values$Year, levels=rev(levels(values$Year)))
	values.agg <- aggregate(number ~ Cancer.label + Year + zone + indice, values, sum)
	values.merge <- merge(values.agg, squamous.infos)
	values.merge$scc <- as.integer(values.merge$percentscc * values.merge$number)
	values.merge$Cancer.label[which(values.merge$Cancer.label == trachea)] <- lung
    values.merge$Cancer.label <- factor(values.merge$Cancer.label, levels=c(hn, lung, 'Oesophagus', 'Cervix uteri', 'Penis', 'Anus', 'Vulva'))
	assign(paste(zone, indice, sep='.'), values.merge)


}

### STEP 2 - concatenation of information for World, European union and Argentina

eu <- rbind(eu.mort, eu.inc)
world <- rbind(world.mort, world.inc)
argentina <- rbind(argentina.mort, argentina.inc)


### STEP 3 - generate the plots with new cases and death for SCCs
# color-blind-friendly palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot.epidemio <- function(data, title, variable) {
    size.axis <- 20
    size.strip <- 18
    size.title <- 25
    size.legend.title <- 20
    size.legend.text <- 20
    gg.zone <- ggplot(data=data, aes(x=indice, y=scc, fill=Year)) +
      geom_bar(stat="identity") +
      facet_grid(as.formula(paste('~', variable, sep=''))) +
      labs(title=title) +
	  scale_y_continuous(labels = label_number(suffix = "K", scale = 1e-3)) +
      scale_fill_manual(values=cbPalette) +
      theme(axis.text=element_text(size=size.axis),
    		axis.title.x=element_blank(),
    		axis.title.y=element_blank(),
    		strip.text.x=element_text(size=size.strip),
    		plot.title=element_text(size=size.title),
    		legend.title=element_text(size=size.legend.title),
    		legend.text=element_text(size=size.legend.text),
    		)
    return(gg.zone)
}

######################
### SCC statistics ###
######################

### Plots by regions and organs
gg.world <- plot.epidemio(world, 'Squamous cell carcinoma by tumor location / World', 'Cancer.label')
gg.world
ggsave('images/world.png', width=15, height=10)
gg.eu <- plot.epidemio(eu, 'Squamous cell carcinoma  by tumor location / European Union', 'Cancer.label')
gg.eu
ggsave('images/eu.png', width=15, height=10)
gg.argentina <- plot.epidemio(argentina, 'Squamous cell carcinoma by tumor location / Argentina', 'Cancer.label')
gg.argentina
ggsave('images/argentina.png', width=15, height=10)

### Plots by regions over all SCCs
world.agg <- aggregate(scc ~ Year + indice, world, sum)
world.agg$zone <- "World" 
eu.agg <- aggregate(scc ~ Year + indice, eu, sum)
eu.agg$zone <- "European Union"
argentina.agg <- aggregate(scc ~ Year + indice, argentina, sum)
argentina.agg$zone <- "Argentina"
scc.agg <- rbind(world.agg, eu.agg, argentina.agg)

gg.scc <- plot.epidemio(scc.agg, 'Squamous cell carcinoma over all tumor locations', 'zone')
gg.scc
ggsave('images/scc.png', width=15, height=10)
names(all.cancer)[which(names(all.cancer) == "Prediction")] <- "all.cancers"
ratio.scc <- merge(all.cancer[,-1], scc.agg)
ratio.scc$ratio <- round(100*ratio.scc$scc/ratio.scc$all.cancers, 1)
ratio.scc.new.eu <- ratio.scc$ratio[which(ratio.scc$zone == "European Union" & ratio.scc$indice == "New" & ratio.scc$Year == 2020)]
ratio.scc$all.cancers <- round(ratio.scc$all.cancers/10^3, 0)
ratio.scc$scc <- round(ratio.scc$scc/10^3, 0)

# Print the statistics over all SCCs 
print(ratio.scc)

##############################
### All cancers statistics ###
##############################

all.cancer <- all.cancer[order(all.cancer$zone, all.cancer$indice, all.cancer$Year),]

###################
### Cancer cost ###
###################

### cost in million euros
### Hofmarcher et al., The cost of cancer in Europe 2018. Eur J Cancer. 2020.
### source https://doi.org/10.1016/j.ejca.2020.01.011
cancer.cost <-  read.csv(paste('hofmarcher', 'hofmarcher-2020.csv', sep='/'), header=TRUE, stringsAsFactors=FALSE)
cancer.cost.eu <- cancer.cost[which(cancer.cost$euro == 'Y'),]

total.cost.eu <- round(sum(cancer.cost.eu$total)/10^3, 1)
healthcare.cost.eu <- round(sum(cancer.cost.eu$direct.healthcare)/10^3, 1)
drugs.cost.eu <- round(sum(cancer.cost.eu$direct.drugs)/10^3, 1)
informal.cost.eu <- round(sum(cancer.cost.eu$informal.care)/10^3, 1)
productivity.cost.eu <- round((sum(cancer.cost.eu$indirect.mortality) + sum(cancer.cost.eu$indirect.morbidity))/10^3, 1)


# Print the cost for World, European Union and Argentina
all.scc.cost <- data.frame(cost=c('Total', 'Healthcare', 'Drug', 'Informal care', 'Productivity loss'), all.cancers=c(total.cost.eu, healthcare.cost.eu - drugs.cost.eu, drugs.cost.eu, informal.cost.eu, productivity.cost.eu))
all.scc.cost$scc <- round(all.scc.cost$all.cancers*ratio.scc.new.eu/100, 1)
write.table(as.data.frame(t(all.scc.cost)), file='hofmarcher/scc-cost.csv', sep='\t', col.names=FALSE, quote=FALSE)
cat('Total \t All cancers:', total.cost.eu, 'BE', '\t SCC: ', round(total.cost.eu*ratio.scc.new.eu/100, 1), 'BE\n')
cat('Healthcare \t All cancers:', healthcare.cost.eu - drugs.cost.eu, 'BE', '- SCC: ', round((healthcare.cost.eu - drugs.cost.eu)*ratio.scc.new.eu/100, 1), 'BE\n')
cat('Drug \t All cancers:', drugs.cost.eu, 'BE', '\t SCC: ', round(drugs.cost.eu*ratio.scc.new.eu/100, 1), 'BE\n')
cat('Informal care \t All cancers:', informal.cost.eu, 'BE', '\t SCC: ', round(informal.cost.eu*ratio.scc.new.eu/100, 1), 'BE\n')
cat('Productivity loss \t All cancers:', productivity.cost.eu, 'BE', '\t SCC: ', round(productivity.cost.eu*ratio.scc.new.eu/100, 1), 'BE\n')

#################################################
### Cost of immunotherapy
### First-line immunotherapy (platinum-sensitive disease, CPS≥1), Pembrolizumab 200 mg (flat dose) q21 until disease progression or unacceptable toxicity
### 7,600 euro/cycle
### Median PFS (CPS≥1): 3.2 months at least half of patients are treated for at least 12-13 weeks: at least 4 administrations
### 30,400 euro for 3 months of treatment (per-patient cost)
### One year of treatment (52 weeks/3-week cycles = 17 cycles) of a single patient: 129,200 euro
#################################################


# cost for immunotherapy during 12 weeks (euros) 
cost.immuno.12weeks <- 30400

# percentage of patients who will respond to immunotherapy
immuno.response <- 0.2 

### Patient eligible for immuno:
### Cost of immunotherapy
# 50% of patients with Head and Neck cancer will receive immunotherapy
# 80% of patients with squamous lung cancer will receive immunotherapy
# 40% of patients with cervix cancer will receive immunotherapy

df.immuno.12weeks <- data.frame(Cancer.label=c('Head and neck', 'Lung', 'Cervix uteri'), immuno.eligibility=c(0.5, 0.8, 0.4))
df.immuno.12weeks <- merge(df.immuno.12weeks, rbind(world,eu, argentina))
ind.new <- which(df.immuno.12weeks$indice == 'New' & df.immuno.12weeks$Year == '2020')
df.immuno.12weeks <- df.immuno.12weeks[ind.new,]
df.immuno.12weeks$immuno.cost <- df.immuno.12weeks$scc * df.immuno.12weeks$immuno.eligibility * cost.immuno.12weeks
df.immuno.12weeks$immuno.cost.no.response <- df.immuno.12weeks$immuno.cost * (1- immuno.response)


plot.immuno.cost <- function(data, title){
	#
    gg.immuno <- ggplot(data=data,
    						aes(x=Year,y=immuno.cost.no.response, fill=Cancer.label)) +
    					geom_bar(stat="identity") +
          				labs(title=title) +
    	  				scale_y_continuous('Cost (billion euros)', labels = label_number(suffix = "B", scale = 1e-9)) +
          				scale_fill_manual(values=cbPalette) +
    				    labs(fill="SCC") +
						theme(axis.title.x=element_blank(), axis.text.x=element_blank(), text=element_text(size=35), axis.text = element_text(size = 45))
	return(gg.immuno)

}

cost.immuno.agg <- aggregate(immuno.cost.no.response ~ zone  ,df.immuno.12weeks, sum)
cat('cost for immunotherapy (patient with no response)\n')

gg.immuno.world <- plot.immuno.cost(df.immuno.12weeks[which(df.immuno.12weeks$zone == 'world'),], 'World')
gg.immuno.world
ggsave('images/world-immuno.png', width=7, height=10)
gg.immuno.eu <- plot.immuno.cost(df.immuno.12weeks[which(df.immuno.12weeks$zone == 'eu'),], 'Europe')
gg.immuno.eu
ggsave('images/eu-immuno.png', width=7, height=10)
gg.immuno.argentina <- plot.immuno.cost(df.immuno.12weeks[which(df.immuno.12weeks$zone == 'argentina'),], 'Argentina')
gg.immuno.argentina
ggsave('images/argentina-immuno.png', width=7, height=10)

