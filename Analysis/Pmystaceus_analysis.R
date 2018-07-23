#----------------------------------------------------------------------------------#
# Whiting, Noble & Qi. Diemetic signal in Phyrnocephalus mystaceus
# Visual modeling code and statistical analysis
# Date: 2017-04-10 09:43:39 AEST
# Contact: Daniel Noble - daniel.noble@unsw.edu.au (UNSW)
#----------------------------------------------------------------------------------#
# Clean working space
	rm(list = ls())

# Load packages
	library(pavo)
	library(plotrix)
	library(plyr)
	library(MCMCglmm)
	library(Hmisc)

## Load functions. Note for Martin: Make sure you set the working dir in R before running code setwd("~/Desktop/Dropbox/Phrynocephalus_mystaceus/")
	source("./Analysis/func.R")

# Morphology data for lizards
	morph <- read.csv("./Data/morph.csv")
	adults <- na.omit(morph[morph$Sex %in% c("f", "m"), c("ID", "Sex","SVL", "right.flap.height", "right.flap.length", "HW", "HL", "Tail", "Mass")])

## Import spec data for lizards
	specs <- getspec(where = "Data/mystaceus_all_spec_files", lim = c(300,700), ext = 'txt')
	specs <- procspec(specs, fixneg = "addmin")

##----------------------------- BODY REGIONS ANALYSIS------------------------------#
## Extract the body regions of interest for the analysis.
## Body region vector
	regions    <- c("mouth", "flap", "dorsum")

## Subset data and store in a list
	subsetdat <- list()

	for(i in 1:length(regions)){
		subsetdat [[i]] <- specs[,c("wl", grep(regions[i], colnames(specs), value = TRUE))]	
	}
	names(subsetdat) <- regions

# Body regions - Smooth spectral curves, fix negative reflectance and prepare for plotting and visual modeling. 

	regions_specs_sm <- list()

	for(i in 1:length(regions)){
		regions_specs_sm[[i]]   <- procspec(subsetdat[[regions[i]]], opt = "smooth" , span = 0.2, fixneg = 'addmin')
	}
	names(regions_specs_sm) <- regions

	sand <- getspec(where = "Data/Tukai_sand/", lim = c(300,700), ext = 'txt')

# Background - Smooth spectral curves and prepare for plotting and visual modeling
	 sand_sm <- procspec(sand, opt = "smooth" , span = 0.2, fixneg = 'addmin')
	sand_avg <- aggspec(sand_sm)

#--------------------------- INDIVIDUAL SPECTRAL ANALYSIS -------------------------------#

## Condense spectral curves across individuals - i.e. average spectra within individuals for body regions

	ind_agg <- list()

	ind <- lapply(regions_specs_sm, function(x) {substr(names(x[2:length(x)]), 1, 7)})

	for(i in 1:length(regions_specs_sm)){
		ind_agg[[i]] <- aggspec(regions_specs_sm[[i]], by = ind[[i]])
	}
	names(ind_agg) <- regions

	lapply(ind_agg, function(x) unique(colnames(x)))

# Add individual background - sand

	ind_agg <- lapply(ind_agg, function(x) cbind(x, sand = sand_avg[,2]))
	lapply(ind_agg, function(x) colnames(x))

#Irradiance for one measurement

	    IRA <- getspec(where="Data/mystaceus_irradiance_FullSunHorizontal")
	 IRA_sm <- procspec(IRA, opt = "smooth" , span = 0.2, fixneg="addmin")
	IRA_avg <- aggspec(IRA_sm)


##-------------------------------VISUAL MODELING-----------------------------------------#

# Create a dataframe containing the visual sensitivities of cones from lizards - taken from pg 1336 of Chan et al. 2009. Behavioural Ecology, 20: 1334-1342

# Lizard visual system
	liz_vis <- sensmodel(c(360, 440, 493, 571)) 

	visual_ind_liz    <- lapply(ind_agg, function(x) vismodel(x, visual = liz_vis, illum= IRA_avg[,2], bkg = sand_avg[,2], relative = FALSE, achromatic = "l", vonkries = TRUE, qcatch = "fi"))

	coldistance_liz   <- lapply(visual_ind_liz, function(x) coldist(x, achro = TRUE, noise = "neural", n = c(1, 1, 3.5, 6),  weber = 0.1))

	# Here, we only want the individual contrasts against the "sand" background. So, subset these data.
	ind_vs_bkg_liz    <- lapply(coldistance_liz, function(x) x[x$patch2 == "sand",])

	lapply(ind_vs_bkg_liz, head)
	lapply(ind_vs_bkg_liz, tail)

## Now that we have the ds and dl contrasts against the background using a lizard visual system we can test whether "greater contrast against a background based on a lizard visual system" is associated with morphology and performance controlling for sex differences. We are making an important assumption that ds between leaves and lizard patches across lizards can be detected by lizards. In other words, the difference between ds for each lizard are discriminable. Pretty big assumption, but seems to be used in the colour literature (See Chan et al. 2009. BE, 20:1334-1342 as example. They implicitly assume that higher ds values of females are more discriminable as they increase based on measurements across the season.)

# Bird visual system
	visual_ind_bird  <- lapply(ind_agg, function(x) vismodel(x, visual = "avg.uv", illum= IRA_avg[,2], bkg = sand_avg[,2], achromatic = "l", relative = FALSE, vonkries = TRUE, qcatch = "fi"))

	coldistance_bird <- lapply(visual_ind_bird, function(x) coldist(x, achro = TRUE, noise = "neural", n = c(1,2,3,3), weber = 0.1))

	# Here, we only want the individual contrasts against the "sand" background. So, subset these data.
	ind_vs_bkg_bird  <- lapply(coldistance_bird, function(x) x[x$patch2 == "sand",])

	lapply(ind_vs_bkg_bird, head)

# Snake visual system. 
	snakeSensModel <- sensmodel(c(360, 482, 554)) # Taken from Sillman et al. 1997.
	
	visual_ind_snake  <- lapply(ind_agg, function(x) vismodel(x, visual = snakeSensModel, illum= IRA_avg[,2], bkg = sand_avg[,2], achromatic = "l",relative = FALSE, vonkries = TRUE, qcatch = "fi"))

	coldistance_snake <- lapply(visual_ind_snake, function(x) coldist(x, achro = TRUE, noise = "neural", n = c(1,1.6,7.3), weber = 0.1))

	# Here, we only want the individual contrasts against the "sand" background. So, subset these data.
	ind_vs_bkg_snake  <- lapply(coldistance_snake, function(x) x[x$patch2 == "sand",])

	lapply(ind_vs_bkg_snake, head)

# Add ind ID and Sex to the data frames so that they can be subset and grouped as well as matched with IDs

	ind_Sex_col <- lapply(ind_agg, function(x) substr(colnames(x[-1]), 1, 1))

	aggregate <- list()

	for(i in 1:3){
		aggregate[[i]] <- aggspec(ind_agg[[i]], by = ind_Sex_col[[i]])
	}

	error <- list()
	for(i in 1:3){
		error[[i]] <- aggspec(ind_agg[[i]], by = ind_Sex_col[[i]], function(x) sd(x)/sqrt(length(x)))
	}

	names(aggregate) <- names(ind_agg)
	names(error) <- names(ind_agg)

#Create sex vector
	sex  <- lapply(ind_vs_bkg_bird, function(x) substr(x[,1],1,1))
	ind_ID_col <- lapply(ind_vs_bkg_bird, function(x) gsub("[fjm].liz","", x[,1]))

	sex2  <- lapply(ind_vs_bkg_snake, function(x) substr(x[,1],1,1))
	ind_ID_col2 <- lapply(ind_vs_bkg_snake, function(x) gsub("[fjm].liz","", x[,1]))

# Create a new data frame for each region and add sex
	birdJND <- list()

	for(i in 1:3){
		birdJND[[i]] <- cbind(ind_vs_bkg_bird[[i]], sex = sex[[i]], ID = ind_ID_col[[i]])	
	}
	names(birdJND) <- names(ind_vs_bkg_bird)

	snakeJND <- list()

	for(i in 1:3){
		snakeJND[[i]] <- cbind(ind_vs_bkg_snake[[i]], sex = sex2[[i]], ID = ind_ID_col2[[i]])	
	}
	names(snakeJND) <- names(ind_vs_bkg_snake)

# Match the morphology data with the ID column so that chromatic and achromatic contrasts match.
	ind_vs_bkg_bird2 <- lapply(birdJND, function(x) merge(x, adults, by = "ID"))
	ind_vs_bkg_snake2 <- lapply(snakeJND, function(x) merge(x, adults, by = "ID"))

# Test for differences between the sexes
	savemodel = FALSE
	
	# Model birds
	modelsbirdMass <- lapply(ind_vs_bkg_bird2, function(x) MCMCglmm(c(dS, dL) ~ sex + Mass, rcov = ~us(trait):units, family = c("gaussian", "gaussian"), nitt = 1000000, thin = 100, data = x))
	
	if(savemodel == TRUE){
		saveRDS(modelsbirdMass, file = "./output/models_sex_MR_dSdL_birds")
	}

	lapply(modelsbirdMass, function(x) summary(x))
	lapply(modelsbirdMass, function(x) plot(x))
	lapply(modelsbirdMass, function(x) autocorr(x$VCV))
	lapply(modelsbirdMass, function(x) autocorr(x$Sol))

	# Model snakes
	modelsSnakeMass <- lapply(ind_vs_bkg_snake2, function(x) MCMCglmm(c(dS, dL) ~ sex + Mass, rcov = ~us(trait):units, family = c("gaussian", "gaussian"), nitt = 1000000, thin = 100, data = x))
	
	if(savemodel == TRUE){
		saveRDS(modelsSnakeMass, file = "./output/models_sex_MR_dSdL_snakes")
	}
	lapply(modelsSnakeMass, function(x) summary(x))
	lapply(modelsSnakeMass, function(x) plot(x))	
	lapply(modelsSnakeMass, function(x) autocorr(x$VCV))
	lapply(modelsSnakeMass, function(x) autocorr(x$Sol))

# Average JNDs for each region and produce se's
	avgbirdJND_sum <- lapply(birdJND, function(x) ddply(x, .(sex), summarise, meandS = mean(dS), sedS = std.error(dS), meandL = mean(dL), sedL = std.error(dL), n = length(dS))) 

	avgsnakeJND_sum <- lapply(snakeJND, function(x) ddply(x, .(sex), summarise, meandS = mean(dS), sedS = std.error(dS), meandL = mean(dL), sedL = std.error(dL), n = length(dS))) 

# Make sure to source the JND plot function
	avgbirdJND_sum <- lapply(avgbirdJND_sum, function(x){rownames(x) <- c("f", "j", "m"); x})

	error_bird <- list()
	error_bird <- lapply(avgbirdJND_sum, function(x) x[,c(3,5)])

	avgsnakeJND_sum <- lapply(avgsnakeJND_sum, function(x){rownames(x) <- c("f", "j", "m"); x})

	error_snake <- list()
	error_snake <- lapply(avgsnakeJND_sum, function(x) x[,c(3,5)])

# Check out sexual dimorphism between flaps
	# Multi-response model. Best because it accounts for covariance between two traits. Also, impact of SVL on both traits in a single analysis
	mod <- MCMCglmm(c(log(right.flap.height), log(right.flap.length)) ~ Sex + log(SVL), family = c("gaussian", "gaussian"), nitt = 1000000, thin = 100, data = adults, rcov = ~us(trait):units)

	summary(mod)
	plot(mod)
	autocorr(mod$VCV)
	autocorr(mod$Sol)

	#mod.height <- glm(log(right.flap.height) ~ Sex + log(SVL), family = "gaussian", data = adults)
	#mod.length <- glm(log(right.flap.length) ~ Sex + log(SVL), family = "gaussian", data = adults)

# Flaps flaring analysis
	matrix <- matrix(c(21,4, 12.5, 12.5), ncol = 2, nrow = 2)
	colnames(matrix) <- c("obs", "exp")
	rownames(matrix) <- c("flare", "noFlare")

	chisq.test(matrix)
	fisher.test(matrix)

#####################################################################################
##---------------------------------- Figures--------------------------------------###
#####################################################################################
	## Figure 2 - Proportion behaviours in different experiments
		## Male and females tethering and enclosure trials
			datProp <- read.csv(file = "Data/Pmystaceus_graphs.csv", header = TRUE)

		# Proportions
			Props <- as.matrix(datProp[, 2:7])
			rownames(Props) <- datProp[,1]

		## Sample sizes
			N <- as.matrix(datProp[, 8:ncol(datProp)])
			rownames(N) <- datProp[,1]
		
		# Here we can conduct a fisher's exact test on counts		
			fisher.test(N)

			name <- substr(colnames(Props), 1, 9)
			name <- gsub("Tethering", "Tether", name)

			pdf(file = "./Figures/Figure2.pdf", width = 7.08, height = 6.30)
				JNDBarplot(data = Props, error = 0, ylab = "Proportions", ylim = c(0, 1.8),  pos = 18, col = c("brown",  "blue", "white"), names.arg = name,  fontsize = 1.5) -> bp.out
				text(x = bp.out, y = Props+0.05, N)
				arrows(x0 = bp.out[1,1], x1 = bp.out[3,2], y0=1.2, y1 = 1.2, length = 0)
				arrows(x0 = bp.out[1,3], x1 = bp.out[3,4], y0=1.2, y1 = 1.2, length = 0)
				arrows(x0 = bp.out[1,5], x1 = bp.out[3,5], y0=1.2, y1 = 1.2, length = 0)
				arrows(x0 = bp.out[1,6], x1 = bp.out[3,6], y0=1.2, y1 = 1.2, length = 0)
				text(paste0("\\MA",":","\\MA"), vfont = c("sans serif", "bold"), xpd = TRUE, x = bp.out[3,2]/2+0.75, y =1.27, cex = 2)
				text(paste0("\\VE",":","\\VE"), vfont = c("sans serif", "bold"), xpd = TRUE, x = bp.out[3,4]-2.85, y =1.27, cex = 2)
				text(paste0("\\VE",":","\\MA"), vfont = c("sans serif", "bold"), xpd = TRUE, x = bp.out[2,5], y =1.27, cex = 2)
				text(paste0("\\MA",":","\\VE"), vfont = c("sans serif", "bold"), xpd = TRUE, x = bp.out[2,6], y =1.27, cex = 2)
				box()
			dev.off()

	## Figure 3 - dS and dL graphs for each sex and body region
		pdf(file = "./Figures/Figure3.pdf", height = 5.973568, width = 8.86)
			reg <- c("flap", "mouth", "dorsum")
			par(mfrow = c(2, 3),  cex.lab = 1.2, mgp = c(1.8,0.5,0), mar = c(4,3,1,0.4))
			for(i in reg){
				JNDBarplot(data = as.matrix(avgbirdJND_sum[[i]][,c(2,4)]), error = as.matrix(error_bird[[i]]), ylab = "", ylim = c(0, 10), name = capitalize(i), pos = 9, col = c("brown", "white", "blue"), names.arg = c("Chromatic", "Achromatic"), fontsize = 1.5, las = 1)
				box()
			}

			for(i in reg){
				JNDBarplot(data = as.matrix(avgsnakeJND_sum[[i]][,c(2,4)]), error = as.matrix(error_snake[[i]]), ylab = "", ylim = c(0, 10), name = capitalize(i), pos = 9, col = c("brown", "white", "blue"), names.arg = c("Chromatic", "Achromatic"), fontsize = 1.5, las = 1)
				box()
			}

			mtext("Just Noticeable Differences (JNDs)", side = 2, outer = TRUE, adj = 0.5, padj = 1.5)
		dev.off()

	## Figure 4 - Spectral reflectance curves
		pdf(file="./Figures/Figure3.pdf", height = 4.09, width = 10.74 )
			par(mfrow = c(1,3),mar = c(4,3.5,1,0.3), mgp = c(2,0.5,0))
			region <- c("flap", "mouth", "dorsum")
			ylabel <- c("Reflectance (%)", "", "")
			xlabel <- c("", "Wavelength (nm)", "")

			#re-order if needed
			aggregate_order <- aggregate[region]
			    error_order <- error[region]
			for(i in 1:3){
				plotcol(aggregate_order[[i]], error_order[[i]], region = capitalize(region[i]), x0 = 300, x1 = 320, y0 = 48, y1 = 48, ylab = 	ylabel [i], xlab = xlabel[i], ylim = c(0,50), cex.lab = 1.5)
			}
		dev.off()

	## Figure 5 - Bird flaring trials tethering
			BirdProp <- read.csv(file = "Data/graphs_birdfield.csv")[1:4,]

		# Behaviour proportions
			props <- as.matrix(BirdProp[,2:4])
			rownames(props) <- firstup(as.character(BirdProp[,1]))

		# N proportions
			N <- as.matrix(BirdProp[,5:ncol(BirdProp)])
			rownames(N) <- firstup(as.character(BirdProp[,1]))

			pdf(file = "./Figures/Figure5.pdf", height = 6.10, width = 6.5)
				JNDBarplot(data = props, error = 0, ylab = "Proportions", ylim = c(0, 1.6), pos = 1.8, col = c("brown",  "blue", "white", "gray"), names.arg = c("Males", "Females", "Juveniles"),  fontsize = 1.5) -> bp.out
				text(x = bp.out, y = props+0.05, N)
				box()
			dev.off()

	## Figure 6 - Noosing trials field
			NooseProp <- read.csv(file = "Data/noosing_fig6.csv")

		# Behaviour proportions

			props <- as.matrix(NooseProp[,2])
			rownames(props) <- NooseProp[,1]

		# N proportions
			N <- NooseProp[,3]
			rownames(N) <- NooseProp[,1]

			pdf(file = "./Figures/Figure6.pdf", height = 6.10, width = 6.5)
				JNDBarplot(data = props, ylab = "Proportions", error = 0, ylim = c(0, 1.8),  pos = 1.8, name = "", col = c("brown",  "blue", "white"), names.arg = c("Flaps flared"),  fontsize = 1.5, space = 0.10) -> bp.out
				text(x = bp.out, y = props/2, paste0(props*100, "%"))
				text(x = bp.out, y = props+0.05, N)
				box()
			dev.off()
