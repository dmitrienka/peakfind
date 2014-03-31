#!/usr/bin/Rscript

options(warn=-1)
Sys.setenv(LANG="EN")
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("pastecs"))

## Functions
	## peak.plot (plots peaklines):
		peak.plot <- function(peak.vector)
	           	  {mapply(function(x)
	                      {abline(v = x, lwd=0.2, col="red")
	                       return(NULL)},
	                     peak.vector)
	              return(NULL)}


start.x <- function(file){
    start <-  grep("start_X", readLines(file), perl=TRUE, value=TRUE)
    if (length(start) == 0)
        FALSE
    else
        as.numeric(gsub("^[^0-9]+([0-9.]+).*$", "\\1", start, perl=TRUE))}

is.shoulder <- function(peak.vector, pattern.deg, pattern.exp, delta)
{
    logik.vector <- mapply(function(peak)
       {
           delta.nh.exp <- pattern.exp[{pattern.deg < {peak + delta}} & {pattern.deg > {peak - delta}}]
           delta.nh.deg <- pattern.deg[{pattern.deg < {peak + delta}} & {pattern.deg > {peak - delta}}]
           if (max(delta.nh.exp) == delta.nh.exp[[1]] | max(delta.nh.exp) == delta.nh.exp[[length(delta.nh.exp)]])
           {return(TRUE)} else {return(FALSE)}
       }, peak.vector)
    return(logik.vector)
}


peak.refine <-  function(peak.vector, pattern.deg, pattern.exp, delta)
{
    refined.matrix <- mapply(function(peak)
      {
          delta.nh.exp <- pattern.exp[{pattern.deg < {peak + delta}} & {pattern.deg > {peak - delta}}]
          delta.nh.deg <- pattern.deg[{pattern.deg < {peak + delta}} & {pattern.deg > {peak - delta}}]
          if (max(delta.nh.exp) == delta.nh.exp[[1]] | max(delta.nh.exp) == delta.nh.exp[[length(delta.nh.exp)]])
          {return(c(peak, delta.nh.exp[[which(delta.nh.deg == peak)]]))}
          else
              {return(c(delta.nh.deg[which(delta.nh.exp == max(delta.nh.exp))][[1]], max(delta.nh.exp)))}

      }, peak.vector)
    refined.deg <- refined.matrix[1,]
    refined.int <- refined.matrix[2,]
    refined.int.rel <- refined.int / max(refined.int) * 100
    return(data.frame(refined.deg, refined.int, refined.int.rel))
}

go.to.dot <- function(dirty.vector, pattern.deg)
{
    mapply(function(dirty.dot)
       {
           pattern.deg[which(abs(pattern.deg - dirty.dot) == min(abs(pattern.deg - dirty.dot)))][[1]]
       }, dirty.vector)
}

read.hkl <- function(file.inp){
    b <- sapply(grep("hkl_m_d_th2",
                     readLines(file.inp),
                     perl=TRUE, value=TRUE),
                function(hkl.string){
                    a <- strsplit(hkl.string, " ")[[1]]
                    return(c(as.integer(a[2]),
                             as.integer(a[3]),
                             as.integer(a[4]),
                             as.integer(a[5]),
                             as.numeric(a[6]),
                             as.numeric(a[7]),
                             as.numeric(strsplit(a[10], "_")[[1]][1])))})
    lp.factor <-  grep("LP_Factor", readLines(file.inp), perl=TRUE, value=TRUE)

    if (length(lp.factor) == 0)
        lp.factor <-  27.3
    else
        lp.factor <- as.numeric(gsub("^[^0-9]+([0-9.]+).*$", "\\1", lp.factor, perl=TRUE))

    lp.scale <- function (th2, i){
        (1 + cos(lp.factor)^2 * cos(th2/180*pi)^2) / (sin(th2/360*pi)^2 * cos(th2/360*pi)) * i}

    return(data.frame(h = b[1,],
                      k = b[2,],
                      l = b[3,],
                      m = b[4,],
                      d = b[5,],
                      th2 = b[6,],
                      I = lp.scale(b[6,], b[7,]), row.names=NULL))
}

append.hkl <- function(peaks.frame, hkl.frame)
{
    hkl.matrix <- mapply(function(peak.dot)
                        {
                            abs.diff <- abs(hkl.frame[[6]] - peak.dot)
                            my.row.number <-  which(abs.diff == min(abs.diff))[[1]]
                            return(c(hkl.frame[[1]][my.row.number],
                                     hkl.frame[[2]][my.row.number],
                                     hkl.frame[[3]][my.row.number]))
                        },
                            peaks.frame[[1]])
    lambda <-  1.540596
    h <- hkl.matrix[1,]
    k <- hkl.matrix[2,]
    l <- hkl.matrix[3,]
    degree <- peaks.frame[[1]]
    d <- lambda/(2 * sin(degree/360*pi))
    d <- as.numeric(sprintf("%.5f", d))
    falses <- rep("false", length(degree))
    dots <- rep(".", length(degree))
    intencity <- peaks.frame[[2]]
    return(data.frame(d, degree, intencity, h, k, l, falses, dots))
}


pattern.print <- function(pattern, n)
{
    if (length(pattern) <= n)
    {
        cat(paste(pattern))
        cat("\n")
    }
    else
    {
        cat(paste(pattern[1:n]))
        cat("\n")
        pattern <- pattern[-(1:n)]
        pattern.print(pattern, n)
    }
}


is.overlaped <- function(degrees, d)  ## ААААА! Быдлокод!!!!
{
    l <- length(degrees)
    mult.flag <- character(l)
    for (i in 2:(l-1))
    {
        if (degrees[i] - degrees[(i-1)] < d |
            degrees[(i+1)] - degrees[i] < d)
        {
            mult.flag[i] <- "true"
        }
        else
        {
            mult.flag[i] <- "false"
        }
    }
    if (degrees[2]-degrees[1] < d)
    {mult.flag[1] <- "true"}else{mult.flag[1] <- "false"}
    if (degrees[l]-degrees[l-1] < d)
    {mult.flag[l] <- "true"}else{mult.flag[l] <- "false"}
    return(mult.flag)
}




option_list <- list(
    make_option(c("-o", "--output"), default="icdd.cif",
                help="Name of output file, [default \"%default\"]"),
    make_option(c("-g", "--graph"), default="check.eps",
                help="Name of the graph file, [default \"%default\"]"),
    make_option(c("-c", "--crop"), type="double", default=60,
                help="Max 2theta, [default %default]"),
    make_option(c("-d", "--delta"), type="double", default=0.037,
                help="Delta for peak position correction [default %default]"),
    make_option(c("-p", "--peak"),  default="empty",
                help="Add peaks manualy. Use notation like this: -p \"1.1 5.987 10\"  "),
    make_option(c("-n", "--DisableGraph"), action="store_true", default=FALSE,
                help="Do not plot the graph"),
    make_option(c("-s", "--shoulders"), action="store_true", default=FALSE,
                help="Remove peak's shoulders"),
    make_option(c("-x", "--startx"),  type="double", default=0,
                help="Start X"),
    make_option(c("-t", "--threshold"), type="double", default=0.35,
                help="About deleting weak peaks [default %default]"),
    make_option(c("-b", "--BkgSub"), action="store_true", default=FALSE,
                help="Write data with background subtracted")			
    )

parser    <- OptionParser(usage = "peakfind.r [options] data.file pro.file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt       <- arguments$options

if (length(arguments$args) != 2)
{
    print("Incorrect number of required positional arguments\n\n")
    print_help(parser)
    stop()
}

## Reading data

mask <- max(start.x(arguments$args[[2]]), opt$startx)

topas.data <- read.table(arguments$args[[1]], sep=",", skip=2, col.names=c("X","Yobs","Ycalc","Ydiff","Bkg"))
topas.data <- topas.data[topas.data$X > mask,]
hkl.data <- read.hkl(arguments$args[[2]])
bkg <- topas.data$Bkg[topas.data$X < opt$crop]
work.area.degree <- topas.data$X[topas.data$X < opt$crop]
work.area.calc <- topas.data$Ycalc[topas.data$X < opt$crop]
work.area.exp <- topas.data$Yobs[topas.data$X < opt$crop]
no.bkg.calc <- work.area.calc - bkg
no.bkg.exp <- work.area.exp - bkg

## Finding peaks


extrema.logic <- extract.turnpoints(turnpoints(work.area.calc), no.tp=FALSE, pit=FALSE, peak=TRUE)
peak.list <- work.area.degree[which(extrema.logic)]
print(paste("Found ", length(peak.list), " peaks"))

if (opt$peak == "empty")
{    print(paste("No manual peaks were added"))} else
{

    peaks.manual <-  as.numeric(strsplit(opt$peak, "[ ]+")[[1]])
    peaks.manual <- go.to.dot(peaks.manual, work.area.degree)

    print(paste(length(peaks.manual), " peaks were added manually"))

    peak.list <- sort(c(peak.list, peaks.manual))
}

## Refining peak positions

print(paste("Refining peaks positions..."))
result <- peak.refine(peak.list, work.area.degree, no.bkg.exp, opt$delta)
print(paste("Maximal differens: ", max(result$refined.deg - peak.list)))

print("Removing low peaks...")
old.length <- length(result$refined.deg)
result <- subset(result, refined.int.rel > opt$threshold)
print(paste( old.length - length(result$refined.deg)  ," weak peaks were removed"))

print("Counting peaks...")
new.crop <-  tail(result$refined.deg, n=1)
if (length(result$refined.deg) > 180)
    {
        result <-  result[1:180,]
        new.crop <-  result$refined.deg[180]
    }

print("Check for collisions...")
l.unique <- length(unique(result$refined.deg))
l.origin <- length(result$refined.deg)
if (l.unique == l.origin)
{    print("All right, no problems with delta value")} else
{    print(paste("WARNING: peak collisions found. Check ", opt$output, " and try to reduce the --delta value"))}

if (opt$shoulders)
{
print("Removing shoulders...")
old.length <- length(result$refined.deg)
result <- result[which(!is.shoulder(result$refined.deg, work.area.degree, no.bkg.exp, opt$delta)), ]
print(paste( old.length - length(result$refined.deg)  ," shoulders  were removed"))
}


## Writing data


#result$refined.int.rel <- as.numeric(sprintf("%.2f", result$refined.int.rel))
result <- append.hkl(result, hkl.data)

my.hkls <- hkl.data[(hkl.data[[6]] < new.crop), ]  # crop from peakcounting
my.hkls <- my.hkls[(my.hkls[[7]] > 0.0000000001), ]
peaks.calc <- data.frame(d = my.hkls[[5]],
                         degree = my.hkls[[6]],
                         intensity = my.hkls[[7]],
                         h = my.hkls[[1]],
                         k = my.hkls[[2]],
                         l = my.hkls[[3]],
                         mult = is.overlaped(my.hkls[[6]], 0.01),
                         dots = rep(".", length(my.hkls[[1]])))

##Just Writing. Kill me, pls.


sink(file=opt$output)
cat("
#################
## DQActivityFile=C:\\ICDDDQActivity\\2013\\Mar\\DQACTIVITY_20130331222759.log#################

data_AUDIT

_audit_creation_date   '2013/03/31 22:32:18'


_audit_creation_method   'DataQuacker CIF Creation: Version 9.b'


_icdd_PDF_vno   9
_icdd_PDF_granteeid   11-03

data_Phase1


#
# Data from PDFEAPC CellPar Panel
#
#CHANGE_ME
#
# cellvalues
#

_space_group_name_H-M_alt  #CHANGE_ME
_cell_length_a  #CHANGE_ME
_cell_length_b  #CHANGE_ME
_cell_length_c  #CHANGE_ME
_cell_angle_alpha  #CHANGE_ME
_cell_angle_beta  #CHANGE_ME
_cell_angle_gamma  #CHANGE_ME
_cell_formula_units_Z  #CHANGE_ME
_cell_volume  #CHANGE_ME
_exptl_crystal_density_meas  ?
_exptl_crystal_density_diffrn  ?
_chemical_formula_weight  ?
_space_group_crystal_system  #CHANGE_ME
_icdd_PDF_space_group_origin  NONE


#
#  experimental conditions of cell related measurements


#

_exptl_crystal_density_meas_temp  ?



#
# Data from PDFEAPC Name Panel_Phase1
#

loop_
_chemical_name_systematic
#CHANGE_ME


#
# Data from PDFEAPC Formula Panel:Phase1
#

_chemical_formula_sum  ?
_chemical_formula_iupac  '#CHANGE_ME'
_chemical_formula_structural  '#CHANGE_ME'
_chemical_formula_analytical  ?



#
# Data from DataQUACkER Instrumentation Panel:Phase1
#

_pd_instr_cons_illum_flag   no
_pd_meas_scan_method   cont
_diffrn_radiation_monochromator   Ge(111)

##################################
_icdd_PDF_intensity_esd  ?
_icdd_PDF_RIR  ?
_pd_meas_2theta_range_min  6
_icdd_PDF_goniometer_radius  217.5
_pd_calib_std_external_name  Si
_pd_instr_geometry  'diffractometer:bragg brentano'


_pd_meas_units_of_intensity   Diffractometer:diffractometer:

#
# Data from PDFEAPC ConditionsPanel:Phase1
#

_pd_prep_temperature  ?


_pd_prep_pressure  ?


_diffrn_ambient_temperature  ?


_diffrn_ambient_pressure  ?


_pd_prep_conditions   ?

#
# Data from PDFEAPC Chemical Properties:Phase1
#

_chemical_melting_point_gt  #CHANGE_ME
_chemical_melting_point_lt  #CHANGE_ME


_chemical_temperature_decomposition  ?


_chemical_temperature_sublimation  ?


_pd_char_colour   colourless
_chemical_compound_source
;
#CHANGE_ME
;

_icdd_pdf_comments_pp   'stable white powder'
#
# FileType OUTPUT
#

loop_
_icdd_PDF_filetype
   O


#
# SubFile OUTPUT
#

loop_
_icdd_PDF_subfiles
   PHR


#
# Data from DataQUACKER null:Phase1
#

loop_
_icdd_PDF_elemental_analysis_type
   Other:chemical


_icdd_PDF_elemental_analysis_details   ?
_icdd_PDF_elemental_analysis_percenttype   Wt
loop_
      _icdd_PDF_elemental_analysis_element
      _icdd_PDF_elemental_analysis_percent
C   #CHANGE_ME
H   #CHANGE_ME
N   #CHANGE_ME

##########################################

#
# Data from DataQUACkER Reference Panel_Phase1
#

loop_
      _icdd_PDF_reference_code
      _journal_coden_astm
      _journal_name_full
      _icdd_PDF_reference_paper_title
      _journal_volume
      _journal_issue
      _journal_year
      _journal_page_first
      _journal_page_last
      _publ_author_name
SM   .   .	.	.	.	.	.	.

##########################################

#
# Data from PDFEAPC Comment Panel:Phase1
#

loop_
      _icdd_PDF_comment_type
      _icdd_PDF_comment
BI   '#CHANGE_ME'
RC
;
The wide peak at 22.6-22.8 degrees 2 theta corresponds to the Kapton thin filmus
ed as a disposable sample holder in the transmission measurements.The peak is ex
cluded from peak lists.

;


##########################################

data_AtomicCoordinates_Phase1

#
# Data from PDFEAPC null
#

_refine_ls_R_factor_all  ?
_refine_ls_wR_factor_all  ?
_pd_proc_ls_prof_R_factor  ?
_pd_proc_ls_prof_wR_factor  ?
_refine_ls_goodness_of_fit_all  ?
_pd_proc_ls_prof_wR_expected  ?
_refine_ls_R_I_factor  ?


#
# Data from PDFEAPC CellPar Panel
#

#
# cellvalues
#
_space_group_name_H-M_alt  #CHANGE_ME
_cell_length_a  #CHANGE_ME
_cell_length_b  #CHANGE_ME
_cell_length_c  #CHANGE_ME
_cell_angle_alpha  #CHANGE_ME
_cell_angle_beta  #CHANGE_ME
_cell_angle_gamma  #CHANGE_ME
_cell_formula_units_Z  #CHANGE_ME
_cell_volume  #CHANGE_ME
_exptl_crystal_density_meas  ?
_exptl_crystal_density_diffrn  ?
_chemical_formula_weight  ?
_space_group_crystal_system  #CHANGE_ME
_icdd_PDF_space_group_origin  NONE

#
# Data from DataQUACkER Symmetry Operators_Phase1
#

_cell_measurement_temperature  ?


_cell_measurement_pressure  ?


_icdd_PDF_global_displacement_parameter_type   NONE
#
# Data from DataQUACKER Coordinate Panel
#

loop_
      _local_atom_site_seq
      _atom_site_label
      _atom_site_type_symbol
      _atom_site_fract_x
      _atom_site_fract_y
      _atom_site_fract_z
      _atom_site_occupancy
      _atom_site_Wyckoff_symbol
      _atom_site_U_iso_or_equiv
 .  .  .  .  .  .  .  .  .

##########################################
#
# Data from DataQUACKER Coordinate Panel
#

loop_
      _atom_site_aniso_label
      _atom_site_aniso_B_11
      _atom_site_aniso_B_22
      _atom_site_aniso_Beta_33
      _atom_site_aniso_B_12
      _atom_site_aniso_B_13
      _atom_site_aniso_B_23
 .  .  .  .  .  .  .

##########################################

data_DIList/PowderData:DILIST-I

#
# Data from DataQUACkER DIList Panel
#

loop_
      _icdd_PDF_comment_type
      _icdd_PDF_comment
DQ   .
D2   .
#
#
#


##################################
_diffrn_radiation_type   CuKa1
_diffrn_radiation_wavelength   1.54056
_diffrn_radiation_probe   x-ray

##################################
_icdd_PDF_dilist_intensity_type   0
loop_
      _refln_d_spacing
      _refln_intensity_meas
      _refln_index_h
      _refln_index_k
      _refln_index_l
      _icdd_pdf_refln_multiple_flag
      _icdd_pdf_dilist_ied
")
    sink()

write.table(result[-2], file=opt$output, row.names=FALSE, col.names=FALSE, quote=FALSE, append = TRUE)

sink(file=opt$output, append = TRUE)
cat("
##########################################
#
# #########
#

data_DIList/PowderData:DILIST-II

#
# Data from DataQUACkER DIList Panel
#

loop_
      _icdd_PDF_comment_type
      _icdd_PDF_comment
DQ   .
D2   Pawley
#
#
#


##################################
_diffrn_radiation_type   CuKa1
_diffrn_radiation_wavelength   1.54056
_diffrn_radiation_probe   x-ray

##################################
_icdd_PDF_dilist_intensity_type   1
loop_
      _refln_d_spacing
      _icdd_pdf_dilist_intensity_int
      _refln_index_h
      _refln_index_k
      _refln_index_l
      _icdd_pdf_refln_multiple_flag
      _icdd_pdf_dilist_ied
")
sink()

    write.table(peaks.calc[-2], file=opt$output, row.names=FALSE, col.names=FALSE, quote=FALSE, append = TRUE)

sink(file=opt$output, append = TRUE)
cat("
##########################################
#
# #########
#

data_DIList/PowderData:RAWDATA

#
# Data from DataQUACkER POWDERPATTERN Panel
")
sink()


sink(file=opt$output, append = TRUE)
    cat(paste("_pd_meas_2theta_range_max   ", topas.data[[1]][[length(topas.data[[1]])]], "\n"))
    cat(paste("_pd_meas_2theta_range_min", topas.data[[1]][[1]], "\n"))
    cat(paste("_pd_meas_2theta_range_inc   ", topas.data[[1]][[2]]-topas.data[[1]][[1]], "\n"))
cat("
##################################
_diffrn_radiation_type   CuKa1
_diffrn_radiation_wavelength   1.54056
_diffrn_radiation_probe   x-ray

##################################
_icdd_PDF_dilist_intensity_type   0
loop_
      _pd_meas_counts_total
")
full.exp <- topas.data[[2]]
if (opt$BkgSub)
{
full.exp <- round(full.exp) - round(topas.data$Bkg)
full.exp <- full.exp - min(full.exp) + 100
work.area.exp <- no.bkg.exp - min(no.bkg.exp) + 100
}
pattern.print(full.exp, 10)
cat("##########################################
#
# #########
#


data_REFERENCEPATT:RAWDATA

#
# Data from DataQUACkER POWDERPATTERN Panel
#


##################################
_diffrn_radiation_type   ?
_diffrn_radiation_wavelength   ?
_diffrn_radiation_probe   ?

##################################
_icdd_PDF_dilist_intensity_type   ?
#
# #########
#
")
sink()


## About Graphs

if (opt$DisableGraph)
{print(paste("DisableGraph option is active, so that's all!"))} else
{
print("Ploting graph...")
setEPS()
postscript(opt$graph, horizontal=FALSE, width=(76/2.5), height=(10/2.5), paper="a4")
par(las=1, mar=c(3, 5, 1, 1), lwd = 0.3, mgp = c(40, 1, 0))
plot.new()
plot(work.area.degree, work.area.exp, type="l", lwd=0.3)
lines(work.area.degree, bkg, type="l", lwd=0.3)
peak.plot(result$degree)
dev.off()
print(paste("Graph is ploted to ", opt$graph))
}
