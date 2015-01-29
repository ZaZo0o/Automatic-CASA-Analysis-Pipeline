######################################
# Automatized Casa Analysis Pipeline
# v0.10
# 2015-01-29
#
# Isabelle Stévant
# University of Geneva


# MIT Licence
# http://opensource.org/licenses/MIT

# Copyright (c) 2015

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

######################################



# /////// ChangeLog /////// #

##### v0.10 #####
# Implement (partial) unit support:
#    - split unit from concentration values
#    - convert mM values to µM
#    - convert concentration values to numerics (proper ordering on the x axis of the graphs)
#####



##### Input File Structure:

# analysis_date_time  ---------------------------  Date of the experiment
# animal_id  ------------------------------------  Molecules (and condition, i.e. if starts with "C_" it is a control)
# lab_tech_name  --------------------------------  ID of the donnor
# motile_mean_alh  ------------------------------  ...
# motile_mean_vcl  ------------------------------  ...
# motile_percent_of_total  ----------------------  ...
# progressive_percent_of_total  -----------------  ...
# sort_a_percent_of_total  ----------------------  ...
# study_number  ---------------------------------  Concentration




list.of.packages <- c("ggplot2", "tcltk", "gridExtra", "stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

require("ggplot2")
require("tcltk")
require("gridExtra")
require("stringr")



done <- tclVar(0)

mydialog <- function(){

	xvar <- tclVar("")
	load_file_name <- tclVar("")
	save_file_location <- tclVar("")
	status <- tclVar(" ")
	tt <- tktoplevel()
	tkwm.title(tt,"Automatized Casa Analysis")
	x.entry <- tkentry(tt, textvariable=xvar)
	entry.load_file_name <-tkentry(tt, width = "50", textvariable = load_file_name)
	entry.save_file_location <-tkentry(tt, width = "50", textvariable = save_file_location)

	submit <- function() {
		tclvalue(status) <- " "
		project_name <- tclvalue(xvar)
		casa_mesures <- read.table(name,header=TRUE, sep="\t")
		assign("project_name", project_name, envir = .GlobalEnv)
		assign("casa_mesures", casa_mesures, envir = .GlobalEnv)
		launchAnalysis(project_name, casa_mesures, result_path)
		tclvalue(status) <- "Done!"
	}


	getfile <- function() {
	    name <- tclvalue(tkgetOpenFile(
		    filetypes = "{{CASA files} {.tab .txt}} {{All files} *}")
		)
		if (name == "") return;
	   	assign("name", name, envir = .GlobalEnv)
	   	tclvalue(load_file_name) <- as.character(name)
	}

	savefile <- function() {
	    name <- tclvalue(tkchooseDirectory())
		if (name == "") return;
	   	assign("result_path", name, envir = .GlobalEnv)
	   	tclvalue(save_file_location) <- as.character(name)
	}

	launchAnalysis <- function(project_name, casa_mesures, result_path){


		result_table <- c(
			"molecule", 
			"n", 
			"conditions",
			"percent_motile", 
			"se_percent_motile", 
			"Paired_t_test_motile",
			"percent_mean_ALH", 
			"se_percent_mean_ALH", 
			"Paired_t_test_ALH",
			"percent_mean_VCL", 
			"se_percent_mean_VCL",
			"Paired_t_test_VCL",
			"percent_mean_progressive", 
			"se_percent_mean_progressive",
			"Paired_t_test_progressive",
			"percent_mean_hyperactive", 
			"se_percent_mean_hyperactive",
			"Paired_t_test_hyperactive"
		)

		getCoef <- function(controlVar){
			100/controlVar
		}

		getPercentCtrl <- function(expVar,controlVar){
			if (controlVar == 0){
				return(NA)
			} else {
				(100*expVar)/controlVar
			}
		}		

		getStdError <- function(expVar){
			sd(expVar, na.rm=TRUE)/sqrt(length(expVar))
		}

		getConcentrationAsNumeric <- function(rawConcentration){
			concentrationValue <- vector()
			for (concentration in rawConcentration) {
				unit <- gsub("[[:digit:][:punct:]]","",concentration)
				if (unit=="mm"){
					value <- as.numeric(gsub("[[:alpha:]]","",concentration))*1000
				} else if (unit=="um"){
					value <- as.numeric(gsub("[[:alpha:]]","",concentration))
				} else {

					value <- "Unknown Unit"
				}
				print(unit)
				concentrationValue <- c(concentrationValue, value)
			}
			return(concentrationValue)
		}


		# Replace NA values by 0
		casa_mesures[is.na(casa_mesures)] <- 0

		# Replace empty string values by 0
		casa_mesures[casa_mesures==""] <- 0

		# Remove space and transform string to lower case to unify concentration values
		casa_mesures$study_number <- tolower(gsub(" ", "", casa_mesures$study_number))
		casa_mesures$analysis_date_time <- str_split_fixed(casa_mesures$analysis_date_time, fixed(" "), 2)[, 1]

		casa_mesures$study_number <- getConcentrationAsNumeric(casa_mesures$study_number)

		expId <- paste(casa_mesures$analysis_date_time, casa_mesures$lab_tech_name, sep="_")

		casa_mesures <- cbind(casa_mesures, expId)

		# Get the list of all the molecules tested and transform to lower case
		molecule_list <- unique(tolower(gsub("C_", "", casa_mesures$animal_id)))


		# For each molecule in the list
		for (molecule in molecule_list) {

			# Isolate the molecule test
			molecule_test <- casa_mesures[grep(paste(molecule,"$", sep=""), tolower(casa_mesures$animal_id)), ]

			# Get replicate list (differentiable with the date of the experiment)
			experiment_list <- unique(molecule_test$expId)

			# Get the number of replicates
			replicates <- length(experiment_list)


			condition_table <- c(
				"experiment",
				"molecule",
				"concentration",
				"percent_motile",
				"percent_alh",
				"percent_vcl",
				"percent_progressive",
				"percent_hyperactive"
			)

			# For each replicate
			for (experiment in experiment_list){

				# Isolate the replicate
				replicate_exp <- molecule_test[grep(paste(experiment,"$",sep=""),molecule_test$expId),]

				# Get the control
				control <- replicate_exp[grep("^C_", replicate_exp$animal_id), ]


				# Get the list of conditions, i.e. the concentrations
				condition_list <- unique(tolower(replicate_exp[grep("C_", replicate_exp$animal_id, invert=T), "study_number"]))


				for (concentration in condition_list) {

					condition <- replicate_exp[grep(paste(paste("^", concentration, sep=""),"$", sep=""), replicate_exp$study_number),]

					percent_motile <- getPercentCtrl(
						replicate_exp[grep(paste(paste("^", concentration, sep=""),"$", sep=""), replicate_exp$study_number),"motile_percent_of_total"],
						control$motile_percent_of_total
					)


					percent_alh <- getPercentCtrl(
						replicate_exp[grep(paste(paste("^", concentration, sep=""),"$", sep=""), replicate_exp$study_number),"motile_mean_alh"],
						control$motile_mean_alh
					)


					percent_vcl <- getPercentCtrl(
						replicate_exp[grep(paste(paste("^", concentration, sep=""),"$", sep=""), replicate_exp$study_number),"motile_mean_vcl"],
						control$motile_mean_vcl
					)


					percent_progressive <- getPercentCtrl(
						replicate_exp[grep(paste(paste("^", concentration, sep=""),"$", sep=""), replicate_exp$study_number),"progressive_percent_of_total"],
						control$progressive_percent_of_total
					)


					percent_hyperactive <- getPercentCtrl(
						replicate_exp[grep(paste(paste("^", concentration, sep=""),"$", sep=""), replicate_exp$study_number),"sort_a_percent_of_total"],
						control$sort_a_percent_of_total
					)
					print(replicate_exp[grep(paste(paste("^", concentration, sep=""),"$", sep=""), replicate_exp$study_number),"sort_a_percent_of_total"])

					condition_replicate <- cbind(
						experiment,
						molecule,
						concentration,
						percent_motile,
						percent_alh,
						percent_vcl,
						percent_progressive,
						percent_hyperactive
					)

					condition_table <- as.data.frame(rbind(
						condition_table,
						condition_replicate
					))

				} # End for concentration

			} # End for experiment

			condition_table <- as.data.frame(condition_table[-1,])


			# condition_table <- as.data.frame(sapply(condition_table,as.numeric))


			control <- molecule_test[grep("^C_", molecule_test$animal_id), ]


			condition_list <- unique(tolower(molecule_test[grep("C_", molecule_test$animal_id, invert=T), "study_number"]))

			result_table <- rbind(result_table, c(molecule, replicates, 0, 100, 0, 0, 100, 0, 0, 100, 0, 0, 100, 0, 0, 100, 0, 0),deparse.level = 0)


			for (concentration in condition_list) {

				condition <- as.data.frame(condition_table[grep(paste(paste("^", concentration, sep=""),"$", sep=""), condition_table$concentration),])

				mean_percent_mobile_ctrl <- mean(as.numeric(as.character(condition$percent_motile)), na.rm=TRUE)
				mean_percent_alh_ctrl <- mean(as.numeric(as.character(condition$percent_alh)), na.rm=TRUE)
				mean_percent_vcl_ctrl <- mean(as.numeric(as.character(condition$percent_vcl)), na.rm=TRUE)
				mean_percent_progressive_ctrl <- mean(as.numeric(as.character(condition$percent_progressive)), na.rm=TRUE)
				mean_percent_hyperactive_ctrl <- mean(as.numeric(as.character(condition$percent_hyperactive)), na.rm=TRUE)

				se_percent_mobile_ctrl <- getStdError(as.numeric(as.character(condition$percent_motile)))
				se_percent_alh_ctrl <- getStdError(as.numeric(as.character(condition$percent_alh)))
				se_percent_vcl_ctrl <- getStdError(as.numeric(as.character(condition$percent_vcl)))
				se_percent_progressive_ctrl <- getStdError(as.numeric(as.character(condition$percent_progressive)))
				se_percent_hyperactive_ctrl <- getStdError(as.numeric(as.character(condition$percent_hyperactive)))


				t.test.p.value <- function(a,b) {
					obj<-try(t.test(a,b), silent=TRUE)
					if (is(obj, "try-error")) return(NA) else return(paste("*", obj$p.value, sep=""))
				}

				pt.test.p.value <- function(a,b) {
					obj<-try(t.test(a,b,paired=TRUE), silent=TRUE)
					if (is(obj, "try-error")) return(t.test.p.value(a,b)) else return(obj$p.value)
				}

				if(replicates>1){

					pvalue_motile <- pt.test.p.value(
						control$motile_percent_of_total,
						molecule_test[grep(paste(paste("^", concentration, sep=""),"$", sep=""), molecule_test$study_number), "motile_percent_of_total" ]
					)

					pvalue_alh <- pt.test.p.value(
						control$motile_mean_alh,
						molecule_test[grep(paste(paste("^", concentration, sep=""),"$", sep=""), molecule_test$study_number), "motile_mean_alh" ]
					)

					pvalue_vcl <- pt.test.p.value(
						control$motile_mean_vcl,
						molecule_test[grep(paste(paste("^", concentration, sep=""),"$", sep=""), molecule_test$study_number), "motile_mean_vcl" ]
					)

					pvalue_progressive <- pt.test.p.value(
						control$progressive_percent_of_total,
						molecule_test[grep(paste(paste("^", concentration, sep=""),"$", sep=""), molecule_test$study_number), "progressive_percent_of_total" ]
					)

					pvalue_hyperactive <- pt.test.p.value(
						control$sort_a_percent_of_total,
						molecule_test[grep(paste(paste("^", concentration, sep=""),"$", sep=""), molecule_test$study_number), "sort_a_percent_of_total" ]
					)

				} else {

					pvalue_motile <- 0

					pvalue_alh <- 0

					pvalue_vcl <- 0

					pvalue_progressive <- 0

					pvalue_hyperactive <- 0
				}

				result_table <- rbind(
						result_table,cbind(
							molecule, 
							replicates, 
							concentration, 
							mean_percent_mobile_ctrl,  
							se_percent_mobile_ctrl, 
							pvalue_motile,
							mean_percent_alh_ctrl, 
							se_percent_alh_ctrl, 
							pvalue_alh,
							mean_percent_vcl_ctrl, 
							se_percent_vcl_ctrl,
							pvalue_vcl,
							mean_percent_progressive_ctrl,
							se_percent_progressive_ctrl,
							pvalue_progressive,
							mean_percent_hyperactive_ctrl,
							se_percent_hyperactive_ctrl,
							pvalue_hyperactive
						)
					)



			} # End for concentration

		} # End for molecule


		colnames(result_table)<- result_table[1,]
		result_table<- as.data.frame(result_table[-1,])
		# result_table[is.na(result_table)] <- 0

		result_table$percent_motile <- as.numeric(as.character(result_table$percent_motile))
		result_table$se_percent_motile <- as.numeric(as.character(result_table$se_percent_motile))

		result_table$percent_mean_ALH <- as.numeric(as.character(result_table$percent_mean_ALH))
		result_table$se_percent_mean_ALH <- as.numeric(as.character(result_table$se_percent_mean_ALH))

		result_table$percent_mean_VCL <- as.numeric(as.character(result_table$percent_mean_VCL))
		result_table$se_percent_mean_VCL <- as.numeric(as.character(result_table$se_percent_mean_VCL))

		result_table$percent_mean_progressive <- as.numeric(as.character(result_table$percent_mean_progressive))
		result_table$se_percent_mean_progressive <- as.numeric(as.character(result_table$se_percent_mean_progressive))

		result_table$percent_mean_hyperactive <- as.numeric(as.character(result_table$percent_mean_hyperactive))
		result_table$se_percent_mean_hyperactive <- as.numeric(as.character(result_table$se_percent_mean_hyperactive))



		# result_table[is.na(result_table)] <- 0


		write.table(result_table, file=paste(result_path,paste(paste("stats",project_name, sep="_"),"txt",sep="."), sep="/"), row.names=F)

		molecule_list2 <- unique(result_table$molecule)

		for (molecule in molecule_list2){
			
			experiment2 <- result_table[grep(paste(molecule,"$", sep=""), result_table$molecule), ]

			experiment2$percent_motile <- as.numeric(as.character(experiment2$percent_motile))
			experiment2$se_percent_motile <- as.numeric(as.character(experiment2$se_percent_motile))
			experiment2$Paired_t_test_motile <- strtrim(experiment2$Paired_t_test_motile,5)
			experiment2$Paired_t_test_motile[experiment2$condition=="0.00nM"] <- ""
			experiment2$Paired_t_test_motile[experiment2$Paired_t_test_motile==0] <- ""

			experiment2$percent_mean_ALH <- as.numeric(as.character(experiment2$percent_mean_ALH))
			experiment2$se_percent_mean_ALH <- as.numeric(as.character(experiment2$se_percent_mean_ALH))
			experiment2$Paired_t_test_ALH <- strtrim(experiment2$Paired_t_test_ALH,5)
			experiment2$Paired_t_test_ALH[experiment2$condition=="0.00nM"] <- ""
			experiment2$Paired_t_test_ALH[experiment2$Paired_t_test_ALH==0] <- ""

			experiment2$percent_mean_VCL <- as.numeric(as.character(experiment2$percent_mean_VCL))
			experiment2$se_percent_mean_VCL <- as.numeric(as.character(experiment2$se_percent_mean_VCL))
			experiment2$Paired_t_test_VCL <- strtrim(experiment2$Paired_t_test_VCL,5)
			experiment2$Paired_t_test_VCL[experiment2$condition=="0.00nM"] <- ""
			experiment2$Paired_t_test_VCL[experiment2$Paired_t_test_VCL==0] <- ""

			experiment2$percent_mean_progressive <- as.numeric(as.character(experiment2$percent_mean_progressive))
			experiment2$se_percent_mean_progressive <- as.numeric(as.character(experiment2$se_percent_mean_progressive))
			experiment2$Paired_t_test_progressive <- strtrim(experiment2$Paired_t_test_progressive,5)
			experiment2$Paired_t_test_progressive[experiment2$condition=="0.00nM"] <- ""
			experiment2$Paired_t_test_progressive[experiment2$Paired_t_test_progressive==0] <- ""

			experiment2$percent_mean_hyperactive <- as.numeric(as.character(experiment2$percent_mean_hyperactive))
			experiment2$se_percent_mean_hyperactive <- as.numeric(as.character(experiment2$se_percent_mean_hyperactive))
			experiment2$Paired_t_test_hyperactive <- strtrim(experiment2$Paired_t_test_hyperactive,5)
			experiment2$Paired_t_test_hyperactive[experiment2$condition=="0.00nM"] <- ""
			experiment2$Paired_t_test_hyperactive[experiment2$Paired_t_test_hyperactive==0] <- ""



			pd <- position_dodge(.5)

			n <- experiment2[experiment2[,"conditions"]==0,"n"]

			my_grob = grobTree(textGrob(paste("n=",n,sep=""), x=0.85,  y=0.95, hjust=0, gp=gpar(fontsize=10)))

			p_mobile <- ggplot(experiment2, aes(x=conditions, y=percent_motile, color=molecule, ymin=percent_motile-se_percent_motile, ymax=percent_motile+se_percent_motile), position=pd) + 
				geom_errorbar(aes(ymin=percent_motile-se_percent_motile, ymax=percent_motile+se_percent_motile), width=.1, position=pd) +
				geom_line(aes(group=molecule), position=pd) +
				geom_point(position=pd, size=3, shape=21, fill="white") + 
				theme(legend.position="bottom") +
				annotation_custom(my_grob)+
				geom_text(aes(label=Paired_t_test_motile), hjust=-0.5, vjust=1.2, colour = "black", size=3)+
				scale_x_discrete(name=expression(paste("Concentration (", mu, "M)"))) +

				scale_y_continuous(name="% of Motile Sperm Control")

			p_ALH <- ggplot(experiment2, aes(x=conditions, y=percent_mean_ALH, color=molecule, ymin=percent_mean_ALH-se_percent_mean_ALH, ymax=percent_mean_ALH+se_percent_mean_ALH), position=pd) + 
				geom_errorbar(aes(ymin=percent_mean_ALH-se_percent_mean_ALH, ymax=percent_mean_ALH+se_percent_mean_ALH), width=.1, position=pd) +
				geom_line(aes(group=molecule), position=pd) +
				geom_point(position=pd, size=3, shape=21, fill="white") + 
				theme(legend.position="bottom") +
				annotation_custom(my_grob)+
				geom_text(aes(label=Paired_t_test_ALH), hjust=-0.5, vjust=1.2, colour = "black", size=3)+
				scale_x_discrete(name=expression(paste("Concentration (", mu, "M)"))) +

    			scale_y_continuous(name="% of ALH Control")

			p_VCL <- ggplot(experiment2, aes(x=conditions, y=percent_mean_VCL, color=molecule, ymin=percent_mean_VCL-se_percent_mean_VCL, ymax=percent_mean_VCL+se_percent_mean_VCL), position=pd) + 
				geom_errorbar(aes(ymin=percent_mean_VCL-se_percent_mean_VCL, ymax=percent_mean_VCL+se_percent_mean_VCL), width=.1, position=pd) +
				geom_line(aes(group=molecule), position=pd) +
				geom_point(position=pd, size=3, shape=21, fill="white") + 
				theme(legend.position="bottom") +
				annotation_custom(my_grob)+
				geom_text(aes(label=Paired_t_test_VCL), hjust=-0.5, vjust=1.2, colour = "black", size=3)+
				scale_x_discrete(name=expression(paste("Concentration (", mu, "M)"))) +

    			scale_y_continuous(name="% of VCL Control")

			p_prog <- ggplot(experiment2, aes(x=conditions, y=percent_mean_progressive, color=molecule, ymin=percent_mean_progressive-se_percent_mean_progressive, ymax=percent_mean_progressive+se_percent_mean_progressive), position=pd) + 
				geom_errorbar(aes(ymin=percent_mean_progressive-se_percent_mean_progressive, ymax=percent_mean_progressive+se_percent_mean_progressive), width=.1, position=pd) +
				geom_line(aes(group=molecule), position=pd) +
				geom_point(position=pd, size=3, shape=21, fill="white") + 
				theme(legend.position="bottom") +
				annotation_custom(my_grob)+
				geom_text(aes(label=Paired_t_test_progressive), hjust=-0.5, vjust=1.2, colour = "black", size=3)+
				scale_x_discrete(name=expression(paste("Concentration (", mu, "M)"))) +

    			scale_y_continuous(name="% of Progressive Sperm Control")

			p_hyp <- ggplot(experiment2, aes(x=conditions, y=percent_mean_hyperactive, color=molecule, ymin=percent_mean_hyperactive-se_percent_mean_hyperactive, ymax=percent_mean_hyperactive+se_percent_mean_hyperactive), position=pd) + 
				geom_errorbar(aes(ymin=percent_mean_hyperactive-se_percent_mean_hyperactive, ymax=percent_mean_hyperactive+se_percent_mean_hyperactive), width=.1, position=pd) +
				geom_line(aes(group=molecule), position=pd) +
				geom_point(position=pd, size=3, shape=21, fill="white") + 
				theme(legend.position="bottom") +
				annotation_custom(my_grob)+
				geom_text(aes(label=Paired_t_test_hyperactive), hjust=-0.2, vjust=1.2, colour = "black", size=3)+
				scale_x_discrete(name=expression(paste("Concentration (", mu, "M)"))) +
    			scale_y_continuous(name="% of Hyperactivity Sperm Control")

		# pdf(paste(result_path,paste(paste("graphs",paste(project_name, molecule, sep="_"), sep="_"),"pdf",sep="."), sep="/"))

		plots <- arrangeGrob(
			p_mobile,
			# p_ALH,
			p_VCL,
			p_prog,
			p_hyp,
			ncol=2,
			nrow=2
		)

		ggsave(plots, file=paste(result_path,paste(paste("graphs",paste(project_name, molecule, sep="_"), sep="_"),"pdf",sep="."), sep="/"),  width = 210, height = 210, units = "mm")

		# dev.off()

		}

	}


	choose_file.but1 <- tkbutton(tt, text = "Browse...", command = getfile)
	choose_file.but2 <- tkbutton(tt, text = "Browse...", command = savefile)

	run_status <- tklabel(tt,textvariable=status)

	submit.but <- tkbutton(tt, text="Launch Analysis", command=submit)
	quit.but <- tkbutton(tt, text = "Quit",
		command = function() { 
	    	tclvalue(done) <- 1
			tkdestroy(tt)
	    }
	)

	tkgrid(tklabel(tt,text="Give a name to this project:"), x.entry, pady= 10, padx= 10)
	tkgrid(tklabel(tt,text="Choose a file to analyse:"),entry.load_file_name, choose_file.but1, pady = 10)
	tkgrid(tklabel(tt,text="Choose result path:"),entry.save_file_location, choose_file.but2, pady = 10)
	tkgrid(submit.but, run_status, quit.but, pady= 10, padx= 10)

	tkfocus(tt)
	tkwait.variable(done)
	tkdestroy(tt)


}

mydialog()

