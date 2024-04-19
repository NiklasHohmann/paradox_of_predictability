#### Investigating P-matrices, evolvabilities and conditional evolvabilities  
library(evolvability)
library(nvctr)

#read custom codes

source("pool_cov.R") # Written by Gene Hunt
source("mee312674-sup-0001-appendixs1.r") # From Grabowski and Porto (2017) MEE

#read in data

load('data/Bjorklund_2017.RData')
load('data/Brombacher.et.al_2017_Globoconella_puncticulata.Rdata')
load('data/Brombacher.et.al_2017_Truncorotalia_crassaformis.Rdata')
load('data/bryo.RData')
load('data/Guralnick_et.al.2022_Peromyscus_maniculatus_females.Rdata')
load('data/Guralnick_et.al.2022_Peromyscus_maniculatus_males.Rdata')
load('data/Firmat.et.al_2014.Rdata')
load('data/Geary_conj.RData')
load('data/Geary_dipro.RData')
load('data/Geladi.RData')
load('data/Poseidonamicus_LH.Rdata')
load('data/Mattioli_2021.RData')
load('data/Pimiento.Balk.2015.Rdata')
load('data/Theriot_et_al_2006.Rdata')
load('data/Ulaski_et_al_2022.Rdata')
load('data/Bert_data.RData')
load('data/Voje.et.al.2022.Rdata')
load('data/Waller_2017.RData')
load('data/Witts_2020.Rdata')
load('data/Yamaguchi.et.al_2017.Rdata')

#gather data into one object

full_data = list(
		Arubb,
		Brombacher.et.al_2017_Globoconella_puncticulata,
		Brombacher.et.al_2017_Truncorotalia_crassaformis,
		comb[[1]],
		comb[[2]],
		Firmat.et.al_2014,
		flycatcher_Male_and_Female[[1]],
		flycatcher_Male_and_Female[[2]],
		Guralnick_et.al.2022_Peromyscus_maniculatus_females,
		Guralnick_et.al.2022_Peromyscus_maniculatus_males,
		Hunt_2007[[1]],
		Hunt_2007[[2]],
		Hunt_2007[[3]],
		Hunt_2007[[4]],
		Hunt_2007[[5]],
		Hunt_2007[[6]],
		Hunt_2007[[7]],
		Hunt_2007[[8]],
		Hunt_2007[[9]],
		Mattioli_2021[[1]],
		Mattioli_2021[[2]],
		mvTS,
		mvTS.1,
		mvTS.2,
		mvTS.3,
		mvTS.4,
		Pimiento.Balk.2015,
		Rguat,
		Theriot_2006,
		Voje.et.al.2022,
		Waller_2017[[1]],
		Waller_2017[[2]],
		Waller_2017[[3]],
		Waller_2017[[4]],
		Yamaguchi.et.al_2017,
		comb2[[1]],
		comb2[[2]]		
)

## Creating data frames to store results
results_summary<-matrix(data=NA,nrow= 37, ncol=13)
colnames(results_summary)<-c("taxon", "sex", "e_exponent", "e_exponent_SE", "e_intercept", "e_intercept_SE", "e_r_squared", "e_adj_r_squared", "e_p_value", "time_lapse", "time_step", "data_type", "reference")

results_summary_pooled_P<-results_summary
results_summary_magnitude<-results_summary
results_summary_pooled_P_magnitude<-results_summary

#Nested loop: run from line 85 to 209

for (j in c(c(1:4, 6:10, 13, 15, 17, 20:24, 27:37))){

test_dat <- full_data[[j]]
sample_identifier = test_dat $nsamp
out_matrix_evolvabilities<-matrix(data=NA,nrow=sample_identifier-1, ncol=11)
out_matrix_pooled_P<-matrix(data=NA,nrow=sample_identifier-1, ncol=11)


# Estimating the pooled P across all samples:
pooled_P<-pool_cov(test_dat$S, test_dat$nn, return.list = FALSE)

for (i in 1:(sample_identifier-1)){
  
  #Compute P matrix in i-1 time 
  P_matrix<-test_dat$S[[i]]
  ## Using Grabowski and Porto 2017 to assess whether the P matrix for a given sample/population should be replaced by the pooled P.
  # The maximum inaccuracy (squared deviation) from the mean was set to 0.05). 
  pool.test.cor<-(mean(c(cov2cor(as.matrix(full_data[[j]]$S[[i]])))))^2 # Estimated average integration (R2) in P based on correlation matrix
  pool.test.n<-full_data[[j]]$nvar # dimensionality of the P matrix (number of traits)
  if (howmany(0.05,pool.test.cor, pool.test.n)[2,1] > full_data[[j]]$nn[i]) P_matrix<- pooled_P else P_matrix<-P_matrix

   # Calculate the vector that defines the observed divergence between sample i an i-1
  first<- as.numeric(test_dat $M[i,])
  second<-as.numeric(test_dat $M[i+1,])
  unit_length<-unname(nvctr::unit(second-first))
  
  # The evolvability in the direction of divergence
  out_matrix_evolvabilities[i,1]<-obs.evolvability<-t(unit_length)%*%as.matrix(P_matrix)%*%unit_length
  out_matrix_pooled_P[i,1]<-obs.evolvability<-t(unit_length)%*%as.matrix(pooled_P)%*%unit_length
  
  out_matrix_evolvabilities[i,2]<-abs(mean(second-first)) # Absolute multivariate trait change 
  out_matrix_evolvabilities[i,3]<-sqrt(sum((second-first)^2)) # Magnitude of change (length of the response vector)

  out_matrix_pooled_P[i,2]<-abs(mean(second-first)) # Absolute multivariate trait change
  out_matrix_pooled_P[i,3]<-sqrt(sum((second-first)^2)) # Magnitude of change (length of the response vector)
  
}

##### Individually estimated P per population/sample (as long as they have been estimated based on an N that provides sufficient accuracy)

## Divergence measured as mean absolute multivariate trait change 

results_summary[j,1] <- full_data[[j]] $taxon
results_summary[j,2] <- full_data[[j]] $sex

#evolvability:
results_summary[j,3] <-  lm(log(out_matrix_evolvabilities[,2])~log(out_matrix_evolvabilities[,1]))$ coefficients[2]
results_summary[j,4] <-  summary(lm(log(out_matrix_evolvabilities[,2])~log(out_matrix_evolvabilities[,1])))$ coefficients[2,2]
results_summary[j,5] <-  lm(log(out_matrix_evolvabilities[,2])~log(out_matrix_evolvabilities[,1]))$ coefficients[1]
results_summary[j,6] <-  summary(lm(log(out_matrix_evolvabilities[,2])~log(out_matrix_evolvabilities[,1])))$coefficients[2,1]
results_summary[j,7] <-  summary(lm(log(out_matrix_evolvabilities[,2])~log(out_matrix_evolvabilities[,1])))$ r.squared
results_summary[j,8] <-  summary(lm(log(out_matrix_evolvabilities[,2])~log(out_matrix_evolvabilities[,1])))$ adj.r.squared
results_summary[j,9] <-  summary(lm(log(out_matrix_evolvabilities[,2])~log(out_matrix_evolvabilities[,1])))$ coefficients[2,4]

# Saving info about the time series: 
results_summary[j,10] <-  ifelse(full_data[[j]] $time.units == "Myr", as.numeric(max(full_data[[j]] $tt)*10^6), as.numeric(max(full_data[[j]] $tt)))
results_summary[j,11] <- length(full_data[[j]] $nn) 
results_summary[j,12] <- ifelse(full_data[[j]] $time.units == "Myr", "paleo", "neo")
results_summary[j,13] <- full_data[[j]] $reference


## Divergence measured as magnitude of trait change (length of response vector) 

results_summary_magnitude[j,1] <- full_data[[j]] $taxon
results_summary_magnitude[j,2] <- full_data[[j]] $sex

# evolvability:
results_summary_magnitude[j,3] <-  lm(log(out_matrix_evolvabilities[,3])~log(out_matrix_evolvabilities[,1]))$ coefficients[2]
results_summary_magnitude[j,4] <-  summary(lm(log(out_matrix_evolvabilities[,3])~log(out_matrix_evolvabilities[,1])))$ coefficients[2,2]
results_summary_magnitude[j,5] <-  lm(log(out_matrix_evolvabilities[,3])~log(out_matrix_evolvabilities[,1]))$ coefficients[1]
results_summary_magnitude[j,6] <-  summary(lm(log(out_matrix_evolvabilities[,3])~log(out_matrix_evolvabilities[,1])))$coefficients[2,1]
results_summary_magnitude[j,7] <-  summary(lm(log(out_matrix_evolvabilities[,3])~log(out_matrix_evolvabilities[,1])))$ r.squared
results_summary_magnitude[j,8] <-  summary(lm(log(out_matrix_evolvabilities[,3])~log(out_matrix_evolvabilities[,1])))$ adj.r.squared
results_summary_magnitude[j,9] <-  summary(lm(log(out_matrix_evolvabilities[,3])~log(out_matrix_evolvabilities[,1])))$ coefficients[2,4]

# Saving info about the time series: 
results_summary_magnitude[j,10] <-  ifelse(full_data[[j]] $time.units == "Myr", as.numeric(max(full_data[[j]] $tt)*10^6), as.numeric(max(full_data[[j]] $tt)))
results_summary_magnitude[j,11] <- length(full_data[[j]] $nn) 
results_summary_magnitude[j,12] <- ifelse(full_data[[j]] $time.units == "Myr", "paleo", "neo")
results_summary_magnitude[j,13] <- full_data[[j]] $reference

write.csv(out_matrix_evolvabilities, paste0("data/output/",paste(full_data[[j]] $reference, full_data[[j]] $taxon, full_data[[j]] $sex, sep = "_"), ".csv"))

##### Pooled P #####

## Divergence measured as mean absolute multivariate trait change 

results_summary_pooled_P[j,1] <- full_data[[j]] $taxon
results_summary_pooled_P[j,2] <- full_data[[j]] $sex

# evolvability:
results_summary_pooled_P[j,3] <-  lm(log(out_matrix_pooled_P[,2])~log(out_matrix_pooled_P[,1]))$ coefficients[2]
results_summary_pooled_P[j,4] <-  summary(lm(log(out_matrix_pooled_P[,2])~log(out_matrix_pooled_P[,1])))$ coefficients[2,2]
results_summary_pooled_P[j,5] <-  lm(log(out_matrix_pooled_P[,2])~log(out_matrix_pooled_P[,1]))$ coefficients[1]
results_summary_pooled_P[j,6] <-  summary(lm(log(out_matrix_pooled_P[,2])~log(out_matrix_pooled_P[,1])))$coefficients[2,1]
results_summary_pooled_P[j,7] <-  summary(lm(log(out_matrix_pooled_P[,2])~log(out_matrix_pooled_P[,1])))$ r.squared
results_summary_pooled_P[j,8] <-  summary(lm(log(out_matrix_pooled_P[,2])~log(out_matrix_pooled_P[,1])))$ adj.r.squared
results_summary_pooled_P[j,9] <-  summary(lm(log(out_matrix_pooled_P[,2])~log(out_matrix_pooled_P[,1])))$ coefficients[2,4]

results_summary_pooled_P[j,10] <-  ifelse(full_data[[j]] $time.units == "Myr", as.numeric(max(full_data[[j]] $tt)*10^6), as.numeric(max(full_data[[j]] $tt)))
results_summary_pooled_P[j,11] <- length(full_data[[j]] $nn) 
results_summary_pooled_P[j,12] <- ifelse(full_data[[j]] $time.units == "Myr", "paleo", "neo")
results_summary_pooled_P[j,13] <- full_data[[j]] $reference 

## Divergence measured as magnitude of trait change (length of response vector) 

results_summary_pooled_P_magnitude[j,1] <- full_data[[j]] $taxon
results_summary_pooled_P_magnitude[j,2] <- full_data[[j]] $sex

# evolvability:
results_summary_pooled_P_magnitude[j,3] <-  lm(log(out_matrix_pooled_P[,3])~log(out_matrix_pooled_P[,1]))$ coefficients[2]
results_summary_pooled_P_magnitude[j,4] <-  summary(lm(log(out_matrix_pooled_P[,3])~log(out_matrix_pooled_P[,1])))$ coefficients[2,2]
results_summary_pooled_P_magnitude[j,5] <-  lm(log(out_matrix_pooled_P[,3])~log(out_matrix_pooled_P[,1]))$ coefficients[1]
results_summary_pooled_P_magnitude[j,6] <-  summary(lm(log(out_matrix_pooled_P[,3])~log(out_matrix_pooled_P[,1])))$coefficients[2,1]
results_summary_pooled_P_magnitude[j,7] <-  summary(lm(log(out_matrix_pooled_P[,3])~log(out_matrix_pooled_P[,1])))$ r.squared
results_summary_pooled_P_magnitude[j,8] <-  summary(lm(log(out_matrix_pooled_P[,3])~log(out_matrix_pooled_P[,1])))$ adj.r.squared
results_summary_pooled_P_magnitude[j,9] <-  summary(lm(log(out_matrix_pooled_P[,3])~log(out_matrix_pooled_P[,1])))$ coefficients[2,4]

results_summary_pooled_P_magnitude[j,10] <-  ifelse(full_data[[j]] $time.units == "Myr", as.numeric(max(full_data[[j]] $tt)*10^6), as.numeric(max(full_data[[j]] $tt)))
results_summary_pooled_P_magnitude[j,11] <- length(full_data[[j]] $nn) 
results_summary_pooled_P_magnitude[j,12] <- ifelse(full_data[[j]] $time.units == "Myr", "paleo", "neo")
results_summary_pooled_P_magnitude[j,13] <- full_data[[j]] $reference 


}

#Table S2

na.omit(results_summary_magnitude[, c(1:4, 7,12)])

# Figures

library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(tidyverse)
library(scales)

#Figure 1b

#setwd("~/Dropbox/ESEB_Timeseries_analyses_personal/analyses/var_div_dat")#this folder should have data-wise outputs from the loop above

variance_dat <- read.csv("data/output/Björklund_2017_Ficedula_albicollis_male.csv")

Figure_1b <- ggplot(data = as.data.frame(variance_dat), aes(x = log(V1), y = log(V2))) +
  geom_point(size = 1, alpha = 0.1) +
  theme_classic()	+
  geom_line(stat = "smooth", method = "lm", se = FALSE, alpha = .5) +
  xlab("log(variance)") + ylab("log(divergence)") +
  annotate(geom = "text", x = -7.2, y = -6.6, label = expression("= slope"), hjust = 0, size= 2.1)+
  annotate(geom = "text", x = -7.5, y = -8.3, label = expression("Bjorklund et al. 2017"), hjust = 0, size= 2.1)+
  annotate(geom = "segment", x = log(min(as.data.frame(variance_dat)$V1)), xend = log(min(as.data.frame(variance_dat)$V1))+0.5, y = log(min(as.data.frame(variance_dat)$V2))+1.2, , yend = log(min(as.data.frame(variance_dat)$V2))+1.2, arrow = arrow(type = "closed", length = unit(0.02, "npc")), linewidth = .3, lty = "dashed")+
  annotate(geom = "segment", x = log(min(as.data.frame(variance_dat)$V1))+0.5, xend = log(min(as.data.frame(variance_dat)$V1))+0.5, y = log(min(as.data.frame(variance_dat)$V2))+1.25 , yend = log(min(as.data.frame(variance_dat)$V2))+(1.25 + 0.5* 2.4828), arrow = arrow(type = "closed", length = unit(0.02, "npc")), linewidth = .3, lty = "dashed")+
  theme(		axis.text=element_blank(),
  				axis.title=element_text(size=7),
  				axis.ticks = element_blank())

#Here I first correct several entries that were incorrect

results_summary_magnitude_corrected <- as.data.frame(results_summary_magnitude) %>%
	mutate(	time_lapse = ifelse(reference=="Theriot_et_al_2006", as.numeric(time_lapse)*1000, time_lapse),
					data_type = ifelse(reference=="Theriot_et_al_2006", "paleo", data_type),
					data_type = ifelse(reference=="Guralnick.et.al.2022", "neo", data_type),
					data_type = ifelse(reference=="VanBocxlaer_2013", "paleo", data_type),
					data_dype_2 = ifelse (data_type == "neo", "extant data", "fossil data"))

#Run a linear model weigted by the inverse of SE to obtain estimates presented in Figure 1c

summary(lm( e_exponent ~ log10(as.numeric(time_lapse)), weight = 1/as.numeric(e_exponent_SE), data = subset(results_summary_magnitude_corrected)))

#Figure 1c

Figure_1c <- ggplot(data = na.omit(results_summary_magnitude_corrected), aes(x = as.numeric(time_lapse), y = as.numeric(e_exponent), shape = data_dype_2)) +
	geom_point(size = 3, alpha = 0.5) +
	geom_errorbar(aes(ymin = as.numeric(e_exponent) - as.numeric(e_exponent_SE)*1.96, ymax = as.numeric(e_exponent) + as.numeric(e_exponent_SE)*1.96), width = .1, alpha = .5) +
	geom_hline(yintercept = 0, linetype = "dashed") +
	geom_hline(yintercept = 1, linetype = "dotted") +
	geom_abline(slope = -0.02, intercept= 0.65691, color = "grey", lty = 1)+
	theme_classic() +
	xlab("time lapse (years)") + ylab("slope")+
	annotate(geom = "text", x = 0.3, y = -2.5, label = expression("regression: slope = -0.02×log"[10]*"(time lapse)+0.64"), hjust = 0, size= 3)+
	annotate(geom = "text", x = 0.3, y = -3, label = expression("p-value: 0.51, r"^2*": 0.02"), hjust = 0, size= 3)+
	scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) +
	theme(	legend.position = c(0.2,0.9),
				legend.title = element_blank())

#Assenmble panels and write out

row1 <- plot_grid(NULL, NULL, Figure_1b, nrow = 3, rel_heights = c(4,1,4),
  label_size = 12, label_fontface = 'bold', label_fontfamily = 'Helvetica',  labels = c("a", " ", "b"))

plot_grid(row1, Figure_1c, ncol = 2, rel_widths = c(1,2.5),
  label_size = 12, label_fontface = 'bold', label_fontfamily = 'Helvetica',  labels = c(" ", "c"))

ggsave("figs/time_vs_exponent_ver3.pdf", device = "pdf", width = 14, height = 10, units = "cm")

#Rest of the figure was created using keynote

#Figure S2

variance_dat2 <- read.csv("data/output/Brombacher.et.al_2017_Globoconella_puncticulata_unknown.csv")
variance_dat3 <- read.csv("data/output/Brombacher.et.al_2017_Truncorotalia_crassaformis_unknown.csv")

summary(lm(log(V2)~log(V1), data = as.data.frame(variance_dat2)))
summary(lm(log(V2)~log(V1), data = as.data.frame(variance_dat3)))

brom_1 = ggplot(data = as.data.frame(variance_dat2), aes(x = V1, y = V2)) +
  geom_point(size = 1, alpha = 0.1) +
  theme_classic()	+
  geom_line(stat = "smooth", method = "lm", se = FALSE, alpha = .5) +
  xlab("variance (%)") + ylab("divergence (%)") + 
  scale_x_log10() +
  scale_y_log10()

brom_2 = ggplot(data = as.data.frame(variance_dat3), aes(x = V1, y = V2)) +
  geom_point(size = 1, alpha = 0.1) +
  theme_classic()	+
  geom_line(stat = "smooth", method = "lm", se = FALSE, alpha = .5) +
  xlab("variance (%)") + ylab("divergence (%)") + 
  scale_x_log10() +
  scale_y_log10()

plot_grid(brom_1, brom_2, ncol = 2, rel_widths = c(1,1),
  label_size = 10, label_fontface = 'bold', label_fontfamily = 'Helvetica',  labels = c("Globoconella puncticulata", "Truncorotalia crassaformis"))

ggsave("figs/Brombacher_plot.pdf", device = "pdf", width = 15, height = 15/2, units = "cm")

#Figure S3

ggplot(data = as.data.frame(results_summary_magnitude_corrected), aes(x = as.numeric(e_exponent), fill = data_dype_2, linetype = data_dype_2)) +
  geom_histogram(bins = 15, alpha = .3) +
  scale_fill_manual(name = "data_type2", values = c("black", "grey50"))+
  scale_color_manual(name = "data_type2", values = c("black", "grey50"))+
  geom_density(alpha = .1)+
  theme_classic()	+
 xlab("slope") + ylab("count") +
  theme(	legend.title = element_blank(),
  				legend.position = "bottom")

ggsave("figs/hist.pdf", device = "pdf", width = 12, height = 12, units = "cm")

#Figure S4

#Prepare data such that stepwise estimates of variance and divergence are stored in one dataframe

#setwd("data/var_div_dat")#this folder should have data-wise outputs from the loop above

df <-
  paste0("data/output/" , list.files(path = "data/output/", pattern = "*.csv")) %>% 
  map_df(~read_csv(.))
df

names(df) <- c("dummy", "variance", "abs_divergence", "mag_divergence", "data_type", "a", "b" ,"c", "d", "e", "f", "g")

#estimate the relationship, paleontological and neontological data separately

summary(lm(log10(abs_divergence)~ log10(variance), data = subset(df, data_type == "paleo")))
summary(lm(log10(abs_divergence)~ log10(variance), data = subset(df, data_type == "neo")))

#Create a dummy variable for figure legend

df <- df %>%
	mutate(data_type2 = ifelse(data_type=="paleo", "fossil: b = 0.54±0.05 (27%)", "extant: b = 0.43±0.04 (32%)"))

ggplot(data=df, aes(x = variance*100, y = abs_divergence*100, fill = data_type2 , color= data_type2, lty = data_type2)) +
	geom_point(alpha = 0.6, stroke = .5) +
	scale_shape_manual(values = c(1,19))+
	geom_smooth(method = "lm", se = F,  linewidth = .5, alpha = 0.5)+
	scale_fill_manual(name = "data_type2", values = c("black", "grey50"))+
	scale_color_manual(name = "data_type2", values = c("black", "grey50"))+
	scale_linetype_manual(values = c(1, 5))+
	theme_classic() +
	scale_x_log10() +
    scale_y_log10()+
    xlab("variance (%)") + ylab("divergence (%)")+
    theme(	legend.title = element_blank(),
    			legend.position = c(0.3, 0.95),
    			legend.background = element_rect(fill = "transparent", colour = "transparent"))

ggsave('figs/test.pdf',h = 100, w = 100, units = 'mm', scale = 1, device ="pdf")

#Figure S5

#Prepare data

results_summary <-as.data.frame(results_summary)
results_summary_pooled_P <-as.data.frame(results_summary_pooled_P)
results_summary_magnitude <-as.data.frame(results_summary_magnitude)
results_summary_pooled_P_magnitude <-as.data.frame(results_summary_pooled_P_magnitude)

results_summary <- as.data.frame(results_summary) %>%
	mutate(	time_lapse = ifelse(reference=="Theriot_et_al_2006", as.numeric(time_lapse)*1000, time_lapse),
					data_type = ifelse(reference=="Theriot_et_al_2006", "paleo", data_type),
					data_type = ifelse(reference=="Guralnick.et.al.2022", "neo", data_type),
					data_type = ifelse(reference=="VanBocxlaer_2013", "paleo", data_type),
					data_dype_2 = ifelse (data_type == "neo", "extant data", "fossil data"))

results_summary_pooled_P <- as.data.frame(results_summary_pooled_P) %>%
	mutate(	time_lapse = ifelse(reference=="Theriot_et_al_2006", as.numeric(time_lapse)*1000, time_lapse),
					data_type = ifelse(reference=="Theriot_et_al_2006", "paleo", data_type),
					data_type = ifelse(reference=="Guralnick.et.al.2022", "neo", data_type),
					data_type = ifelse(reference=="VanBocxlaer_2013", "paleo", data_type),
					data_dype_2 = ifelse (data_type == "neo", "extant data", "fossil data"))

results_summary_pooled_P_magnitude <- as.data.frame(results_summary_pooled_P_magnitude) %>%
	mutate(	time_lapse = ifelse(reference=="Theriot_et_al_2006", as.numeric(time_lapse)*1000, time_lapse),
					data_type = ifelse(reference=="Theriot_et_al_2006", "paleo", data_type),
					data_type = ifelse(reference=="Guralnick.et.al.2022", "neo", data_type),
					data_type = ifelse(reference=="VanBocxlaer_2013", "paleo", data_type),
					data_dype_2 = ifelse (data_type == "neo", "extant data", "fossil data"))

mag_pool <- ggplot(data = na.omit(results_summary_pooled_P_magnitude)[-21,], aes(x = as.numeric(time_lapse), y = as.numeric(e_exponent), shape = data_dype_2)) +
	geom_point(size = 3, alpha = 0.5) +
	geom_errorbar(aes(ymin = as.numeric(e_exponent) - as.numeric(e_exponent_SE)*1.96, ymax = as.numeric(e_exponent) + as.numeric(e_exponent_SE)*1.96), width = .1, alpha = .5) +
	geom_hline(yintercept = 0, linetype = "dashed") +
	geom_hline(yintercept = 1, linetype = "dotted") +
	geom_abline(slope = -0.005744, intercept= 0.645567, color = "grey", lty = 1)+
	theme_classic() +
	xlab("time lapse (years)") + ylab("slope")+
	ylim(-4, 6.5) +
	annotate(geom = "text", x = 1, y = -2.5, label = expression("regression: slope = -0.0057×log"[10]*"(time lapse)+0.65"), hjust = 0, size= 3)+
	annotate(geom = "text", x = 1, y = -3, label = expression("p-value: 0.89, r"^2*": 0.001"), hjust = 0, size= 3)+
	annotate(geom = "text", x = 1, y = 6, label = "slope: magnitude ~ pooled P", hjust = 0, size= 3)+
	scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) +
	theme(	legend.position = "none")

abs_raw <- ggplot(data = na.omit(results_summary), aes(x = as.numeric(time_lapse), y = as.numeric(e_exponent), shape = data_dype_2)) +
	geom_point(size = 3, alpha = 0.5) +
	geom_errorbar(aes(ymin = as.numeric(e_exponent) - as.numeric(e_exponent_SE)*1.96, ymax = as.numeric(e_exponent) + as.numeric(e_exponent_SE)*1.96), width = .1, alpha = .5) +
	geom_hline(yintercept = 0, linetype = "dashed") +
	geom_hline(yintercept = 1, linetype = "dotted") +
	geom_abline(slope = -0.14559, intercept= 1.55743, color = "grey", lty = 1)+
	theme_classic() +
	xlab("time lapse (years)") + ylab("slope")+
	ylim(-4, 6.5) +
	annotate(geom = "text", x = 1, y = -2.5, label = expression("regression: slope = -0.15×log"[10]*"(time lapse)+1.56"), hjust = 0, size= 3)+
	annotate(geom = "text", x = 1, y = -3, label = expression("p-value: 0.01, r"^2*": 0.21"), hjust = 0, size= 3)+
	annotate(geom = "text", x = 1, y = 6, label = "slope: absolute change ~ raw P", hjust = 0, size= 3)+
	scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) +
	theme(	legend.position = "none")

abs_pool <- ggplot(data = na.omit(results_summary_pooled_P)[-21,], aes(x = as.numeric(time_lapse), y = as.numeric(e_exponent), shape = data_dype_2)) +
	geom_point(size = 3, alpha = 0.5) +
	geom_errorbar(aes(ymin = as.numeric(e_exponent) - as.numeric(e_exponent_SE)*1.96, ymax = as.numeric(e_exponent) + as.numeric(e_exponent_SE)*1.96), width = .1, alpha = .5) +
	geom_hline(yintercept = 0, linetype = "dashed") +
	geom_hline(yintercept = 1, linetype = "dotted") +
	geom_abline(slope = -0.01097, intercept= 1.41632, color = "grey", lty = 1)+
	theme_classic() +
	xlab("time lapse (years)") + ylab("slope")+
	ylim(-4, 6.5) +
	annotate(geom = "text", x = 1, y = -2.5, label = expression("regression: slope = -0.01×log"[10]*"(time lapse)+1.42"), hjust = 0, size= 3)+
	annotate(geom = "text", x = 1, y = -3, label = expression("p-value: 0.87, r"^2*": 0.001"), hjust = 0, size= 3)+
	annotate(geom = "text", x = 1, y = 6, label = "slope: absolute change ~ pooled P", hjust = 0, size= 3)+
	scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) +
	theme(	legend.position = "none")

cowplot::plot_grid(mag_pool, abs_raw, abs_pool , nrow = 1,
  label_size = 12, label_fontface = 'bold', label_fontfamily = 'Helvetica',  labels = c("a", "b", "c"))

ggsave("figs/sensitivity_revision.pdf", device = "pdf", width = 30, height = 10, units = "cm")

#Statistical tests

summary(lm( e_exponent ~ log10(as.numeric(time_lapse)), weight = 1/as.numeric(e_exponent_SE), data = subset(results_summary_magnitude_corrected)))
summary(lm( e_exponent ~ log10(as.numeric(time_lapse)), weight = 1/as.numeric(e_exponent_SE), data = subset(results_summary_pooled_P_magnitude)))
summary(lm( e_exponent ~ log10(as.numeric(time_lapse)), weight = 1/as.numeric(e_exponent_SE), data = subset(results_summary)))
summary(lm( e_exponent ~ log10(as.numeric(time_lapse)), weight = 1/as.numeric(e_exponent_SE), data = subset(results_summary_pooled_P)))
