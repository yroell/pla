###################################
# User input (change lines 13-16 as needed):
#   1. path = location where file to read is placed and where results will be stored
#   2. filename = name of file to work with
#       a) must be csv
#       b) column with results needs to be called: Result
#       c) column with dilution information needs to be called: Dilution_Factor
#       d) column with sample id needs to be called: Sample_ID
#       e) column with antigen needs to be called: Antigen
#       f) column with antibody needs to be called Antibody
#   3. lab = unique identifier for results in output files
#   4. standard_id = name of standard within file
path = "/path/to/working_directory/"
filename = "File_with_Results"
lab = "Name_of_lab"
standard_id = "WHO-IS"
###################################

###################################
# install and load necessary packages and load functions
list.of.packages = c("stringr", "dplyr")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)
library(stringr)
library(dplyr)

changeSciNot = function(n) {
  output = format(n, scientific = TRUE, digits = 3) # Transforms the number into scientific notation even if small
  output = sub("e", "*10^", output) # Replace e with 10^
  output = sub("\\+0?", "", output) # Remove + symbol and leading zeros on exponent, if > 1
  output = sub("-0?", "-", output) # Leaves - symbol but removes leading zeros on exponent, if < 1
  output
}

# create folder system, read in data, and clean up data
setwd(path)
dir.create(paste0("./", lab))
dir.create(paste0("./", lab, "/", lab, "_figures"))
dir.create(paste0("./", lab, "/", lab, "_figures/not_parallel"))
x = read.csv(paste0(filename, ".csv"))
x$Sample_ID = str_squish(x$Sample_ID)
x$Result = as.numeric(gsub("%|#VALUE!", "", x$Result))
x = x[!is.na(x$Result), ]
x = x[x$Result > 0, ]
x$Antibody = str_squish(gsub("%", "perc_", x$Antibody))
x$Antibody = gsub("Total \\(IgM/IgG/IgA\\)", "Total", x$Antibody)

# create output files
writefile = paste0("./", lab, "/", lab, "_not_plotted_samples.txt")
cat(paste0("Samples not plotted from ", lab, " lab:\n"),
    file = writefile, append = F)
logfile = paste0("./", lab, "/", lab, "_log_file.txt")
cat(paste0("Log file for ", lab, " lab:\n"),
    file = logfile, append = F)
df_output = data.frame(
  Sample_ID = c(), Antigen = c(), Anitbody = c(),
  Relative_Slope = c(), 
  Unknown_BAU = c(), Unknown_CI = c()
)

# get dose in correct format based on data type of input
if (class(x$Dilution_Factor) == "character") {
  # split dilution character to get numeric value for dose
  x$den = as.numeric(str_split(x$Dilution_Factor, pattern = "/", simplify = TRUE)[, 2])
  x$dose = 1 / x$den
}
if (class(x$Dilution_Factor) == "integer") {
  x$dose = 1 / x$Dilution_Factor
}
if (class(x$Dilution_Factor) == "numeric") {
  x$dose = x$Dilution_Factor
}

# log transform dose and result
x$dose_log = log10(x$dose)
x$result_log = log10(x$Result)

# get all unique combinations of sample_id, antigen, and antibody
just_samples = x %>%
  filter(Sample_ID != standard_id)
u_id = unique(just_samples$Sample_ID)
u_gen = unique(just_samples$Antigen)
u_body = unique(just_samples$Antibody)
u_comb = expand.grid(u_id, u_gen, u_body)
cat(paste0("The name of the standard is: ", standard_id, "\n"),
    file = logfile, append = T)
cat(paste0("There are ", nrow(u_comb), " unique combinations\n"),
    file = logfile, append = T)

# loop through unique combinations
for (i in 1:nrow(u_comb)) {
  cat(paste("\n", i, "begin processing unique combination:",
            u_comb[i, 1], u_comb[i, 2], u_comb[i, 3], "\n"),
      file = logfile, append = T)
  # split data into standard based on sample_id, antigen, antibody
  standard = x %>%
    filter(Sample_ID == standard_id,
           Antigen == u_comb[i, 2],
           Antibody == u_comb[i, 3])
  standard_4_not_parallel = standard
  standard_dilution_count = standard %>%
    group_by(dose) %>%
    summarize(n = n()) %>%
    filter(n >= 3)
  standard = x %>%
    filter(Sample_ID == standard_id,
           Antigen == u_comb[i, 2],
           Antibody == u_comb[i, 3],
           dose %in% standard_dilution_count$dose) %>%
    mutate(group = "standard")
  cat("   get data for standard\n",
      file = logfile, append = T)
  
  # if no data in sample for unique combination, report in file and go to next
  if (nrow(standard) == 0) {
    cat(paste0("COMBINATION NOT IN STANDARD - ",
               i,". Sample_ID = ", u_comb[i, 1],
               "; Antigen = ", u_comb[i, 2],
               "; Antibody = ", u_comb[i, 3], "\n"),
        file = writefile, append = T)
    df_output = rbind(
      df_output,
      data.frame(
        Sample_ID = u_comb[i, 1], Antigen = u_comb[i, 2], Antibody = u_comb[i, 3],
        Relative_Slope = NA,
        Unknown_BAU = NA, Unknown_CI = NA))
    cat("   combination not in standard\n",
        file = logfile, append = T)
    next
  }
  
  # get slope and intercept for standard based on log transformed result
  who_lm = summary(lm(result_log ~ dose_log, standard))
  who_slope = who_lm$coefficients[2, 1]
  who_inter = who_lm$coefficients[1, 1]
  who_r2 = who_lm$r.squared
  cat(paste("   standard slope:", who_slope, "\n", 
            "   standard intercept:", who_inter, "\n",
            "   standard R2:", who_r2, "\n"),
      file = logfile, append = T)
  cat(paste("   number of original dilutions in standard:", nrow(standard_dilution_count), "\n"),
      file = logfile, append = T)
  
  # remove top and bottom dilutions if needed
  top3 = standard %>%
    filter(dose_log %in% sort(unique(standard$dose_log), decreasing = T)[1:3])
  bottom3 = standard %>%
    filter(dose_log %in% sort(unique(standard$dose_log), decreasing = F)[1:3])
  top3_lm = summary(lm(result_log ~ dose_log, top3))
  bottom3_lm = summary(lm(result_log ~ dose_log, bottom3))
  top3_slope = top3_lm$coefficients[2, 1]
  bottom3_slope = bottom3_lm$coefficients[2, 1]
  
  while (((top3_slope <= who_slope - who_slope * 0.2 | top3_slope >= who_slope + who_slope * 0.2) | 
          (bottom3_slope <= who_slope - who_slope * 0.2 | bottom3_slope >= who_slope + who_slope * 0.2)) &
         (length(unique(standard$dose_log)) >= 5) & (who_r2 < 0.95)) {
    if ((top3_slope <= who_slope - who_slope * 0.2 | top3_slope >= who_slope + who_slope * 0.2) &
        (bottom3_slope <= who_slope - who_slope * 0.2 | bottom3_slope >= who_slope + who_slope * 0.2)) {
      cat("   removing upper and lower point from standard\n",
          file = logfile, append = T)
      standard = standard %>%
        filter(dose_log != sort(unique(standard$dose_log), decreasing = T)[1]) %>%
        filter(dose_log != sort(unique(standard$dose_log), decreasing = F)[1])
      who_lm = summary(lm(result_log ~ dose_log, standard))
      who_slope = who_lm$coefficients[2, 1]
      who_inter = who_lm$coefficients[1, 1]
      who_r2 = who_lm$r.squared
      top3 = standard %>%
        filter(dose_log %in% sort(unique(standard$dose_log), decreasing = T)[1:3])
      bottom3 = standard %>%
        filter(dose_log %in% sort(unique(standard$dose_log), decreasing = T)[1:3])
      top3_lm = summary(lm(result_log ~ dose_log, top3))
      top3_slope = top3_lm$coefficients[2, 1]
      bottom3_lm = summary(lm(result_log ~ dose_log, bottom3))
      bottom3_slope = bottom3_lm$coefficients[2, 1]
      next
    } else if (top3_slope <= who_slope - who_slope * 0.2 | top3_slope >= who_slope + who_slope * 0.2) {
      cat("   removing upper point from standard\n",
          file = logfile, append = T)
      standard = standard %>%
        filter(dose_log != sort(unique(standard$dose_log), decreasing = T)[1])
      who_lm = summary(lm(result_log ~ dose_log, standard))
      who_slope = who_lm$coefficients[2, 1]
      who_inter = who_lm$coefficients[1, 1]
      who_r2 = who_lm$r.squared
      top3 = standard %>%
        filter(dose_log %in% sort(unique(standard$dose_log), decreasing = T)[1:3])
      top3_lm = summary(lm(result_log ~ dose_log, top3))
      top3_slope = top3_lm$coefficients[2, 1]
      next
    } else if (bottom3_slope <= who_slope - who_slope * 0.2 | bottom3_slope >= who_slope + who_slope * 0.2) {
      cat("   removing lower point from standard\n",
          file = logfile, append = T)
      standard = standard %>%
        filter(dose_log != sort(unique(standard$dose_log), decreasing = F)[1])
      who_lm = summary(lm(result_log ~ dose_log, standard))
      who_slope = who_lm$coefficients[2, 1]
      who_inter = who_lm$coefficients[1, 1]
      who_r2 = who_lm$r.squared
      bottom3 = standard %>%
        filter(dose_log %in% sort(unique(standard$dose_log), decreasing = F)[1:3])
      bottom3_lm = summary(lm(result_log ~ dose_log, bottom3))
      bottom3_slope = bottom3_lm$coefficients[2, 1]
      next
    } else {
      break
    }
  }
  cat(paste("   number of dilutions in standard after removing upper or lower points:", length(unique(standard$dose)), "\n"),
      file = logfile, append = T)
  cat(paste("   standard slope after removing points:", who_slope, "\n",
            "   standard intercept after removing points:", who_inter, "\n",
            "   standard R2 after removing points:", who_r2, "\n"),
      file = logfile, append = T)
  
  # split data into sample based on sample_id, antigen, antibody
  cat("   get data for sample\n",
      file = logfile, append = T)
  samples = x %>%
    filter(Sample_ID == u_comb[i, 1],
           Antigen == u_comb[i, 2],
           Antibody == u_comb[i, 3])
  samples_4_not_parallel = samples
  num_dilutions = length(unique(samples$dose))
  samples_dilution_count = samples %>%
    group_by(dose) %>%
    summarize(n = n()) %>%
    filter(n >= 3)
  samples = x %>%
    filter(Sample_ID == u_comb[i, 1],
           Antigen == u_comb[i, 2],
           Antibody == u_comb[i, 3],
           dose %in% samples_dilution_count$dose,
           dose >= min(standard$dose),
           dose <= max(standard$dose)) %>%
    mutate(group = "samples")
  cat(paste("   number of original dilutions in sample:", nrow(samples_dilution_count), "\n"),
      file = logfile, append = T)
  cat(paste("   number of dilutions in sample after matching with standard:", length(unique(samples$dose)), "\n"),
      file = logfile, append = T)
  
  # if no data in sample for unique combination, report in file and go to next
  if (nrow(samples) == 0) {
    cat(paste0("COMBINATION NOT IN SAMPLE - ",
               i, ". Sample_ID = ", u_comb[i, 1],
               "; Antigen = ", u_comb[i, 2],
               "; Antibody = ", u_comb[i, 3],"\n"),
        file = writefile,append = T)
    df_output = rbind(
      df_output,
      data.frame(
        Sample_ID = u_comb[i, 1], Antigen = u_comb[i, 2], Antibody = u_comb[i, 3],
        Relative_Slope = NA,
        Unknown_BAU = NA, Unknown_CI = NA))
    cat("   combination not in sample\n",
        file = logfile, append = T)
    next
  }
  
  # if sample has less than three dilutions, report in file and go to next
  if (length(unique(samples$dose_log)) < 3) {
    cat(paste0("COMBINATION HAS LESS THAN THREE DILUTION - ",
               i, ". Sample_ID = ", u_comb[i, 1], 
               "; Antigen = ", u_comb[i, 2],
               "; Antibody = ", u_comb[i, 3], "\n"),
        file = writefile, append = T)
    df_output = rbind(
      df_output,
      data.frame(
        Sample_ID = u_comb[i, 1], Antigen = u_comb[i, 2], Antibody = u_comb[i, 3],
        Relative_Slope = NA,
        Unknown_BAU = NA, Unknown_CI = NA))
    cat("   sample has less than three dilutions\n",
        file = logfile, append = T)
    next
  }
  
  # combine standard and sample
  together = rbind(standard, samples)
  together_4_not_parallel = rbind(standard_4_not_parallel, samples_4_not_parallel)
  
  # get slope and intercept for sample based on log transformed result
  sample_lm = summary(lm(result_log ~ dose_log, samples))
  sample_slope = sample_lm$coefficients[2, 1]
  sample_inter = sample_lm$coefficients[1, 1]
  sample_r2 = sample_lm$r.squared
  rel_slope = who_slope / sample_slope
  cat(paste("   sample slope:", sample_slope, "\n",
            "   sample intercept:", sample_inter, "\n",
            "   sample R2:", sample_r2, "\n",
            "   relative slope:", rel_slope, "\n"),
      file = logfile, append = T)
  
  # if sample is within 20% of slope, plot
  if (rel_slope >= 0.8 & rel_slope <= 1.2) {
    cat("   relative slope is within 20%, proceed to calculating BAU and plotting\n",
        file = logfile, append = T)
    # calculate mean and CI for relative potency using bootstrapping
    co = c()
    set.seed(1)
    for (j in 1:100) {
      b_df_standard = data.frame(result_log = c(), dose_log = c())
      b_df_samples = data.frame(result_log = c(), dose_log = c())
      for (k in unique(standard$dose_log)) {
        b_all_standard = subset(standard, dose_log == k)
        b_standard = sample(b_all_standard$result_log, 3, replace = T)
        b_sub_standard = data.frame(result_log = b_standard, dose_log = k)
        b_df_standard = rbind(b_df_standard, b_sub_standard)
      }
      for (k in unique(samples$dose_log)) {
        b_all_samples = subset(samples, dose_log == k)
        b_samples = sample(b_all_samples$result_log, 3, replace = T)
        b_sub_samples = data.frame(result_log = b_samples, dose_log = k)
        b_df_samples = rbind(b_df_samples, b_sub_samples)
      }
      b_who_lm = summary(lm(result_log ~ dose_log, b_df_standard))
      b_sample_lm = summary(lm(result_log ~ dose_log, b_df_samples))
      b_who_slope = b_who_lm$coefficients[2, 1]
      b_sample_slope = b_sample_lm$coefficients[2, 1]
      b_who_inter = b_who_lm$coefficients[1, 1]
      b_sample_inter = b_sample_lm$coefficients[1, 1]
      
      interval_result = seq(
        from = min(b_df_standard$result_log),
        to = max(b_df_standard$result_log),
        length.out = 100
      )
      
      who_potency = (interval_result - b_who_inter) / b_who_slope
      sample_potency = (interval_result - b_sample_inter) / b_sample_slope
      log_rel_potency = who_potency - sample_potency
      mean_rel_pot = mean(10 ^ log_rel_potency) * 1000
      co[j] = mean_rel_pot
    }
    m = round(mean(co))
    e = round(qnorm(0.975) * sd(co) / sqrt(100))
    cil = m - e
    ciu = m + e
    cat(paste("   BAU:", m, "\n",
              "   CI:", e, "\n"),
        file = logfile, append = T)
    df_output = rbind(
      df_output,
      data.frame(
        Sample_ID = u_comb[i, 1], Antigen = u_comb[i, 2], Antibody = u_comb[i, 3],
        Relative_Slope = rel_slope,
        Unknown_BAU = m, Unknown_CI = e))
    
    # create plot and save
    cat("   create plot\n",
        file = logfile, append = T)
    pdf(paste0("./",lab,"/",lab,"_figures/",samples$Sample_ID[1],"_",samples$Antigen[1],"_",samples$Antibody[1],".pdf"))
    plot(standard$dose_log, standard$result_log,
         col = "red", pch = 1,
         ylim = c(min(together$result_log), max(together$result_log)),
         ylab = "Result (log transformed)",
         xaxt = "n", xlab = "Dilution",
         main = paste0(lab," lab\n","Sample_ID = ",samples$Sample_ID[1],"\nAntigen = ",samples$Antigen[1],"; Antibody = ",samples$Antibody[1]))
    points(samples$dose_log, samples$result_log,
           col = "blue", pch = 1)
    abline(b = who_slope, a = who_inter, col = "red")
    abline(b = sample_slope, a = sample_inter, col = "blue")
    legend("bottomright",
           legend = c(standard_id, samples$Sample_ID[1]),
           col = c("red", "blue"), pch = 1, bty = "y")
    legend("topleft",
           legend = c(paste0("BAU = ", m), paste0("Conf. Int. = [", cil, ", ", ciu, "]")),
           bty = "n")
    axis(1, at = sort(unique(samples$dose_log), decreasing = F),
         labels = changeSciNot(sort(unique(samples$dose), decreasing = F)))
    dev.off()
    next
    
  } else {
    cat("   relative slope is not within 20%, remove points and recalculate relative slope\n",
        file = logfile,append = T)
    while (length(unique(samples$dose_log)) >= 3) {
      # drop dilutions from standard and sample until parallel
      standard_top = standard %>%
        filter(dose_log != max(samples$dose_log))
      samples_top = samples %>%
        filter(dose_log != max(samples$dose_log))
      standard_bottom = standard %>%
        filter(dose_log != min(samples$dose_log))
      samples_bottom = samples %>%
        filter(dose_log != min(samples$dose_log))
      
      # if data has less than two dilutions after removing points, report in file and go to next
      if (length(unique(samples_top$dose_log)) == 2) {
        cat(paste0("COMBINATION NOT PARALLEL - ",
                   i, ". Sample_ID = ", samples$Sample_ID[1],
                   "; Antigen = ", samples$Antigen[1],
                   "; Antibody = ", samples$Antibody[1], "\n"),
            file = writefile, append = T)
        df_output = rbind(
          df_output,
          data.frame(
            Sample_ID = u_comb[i, 1], Antigen = u_comb[i, 2], Antibody = u_comb[i, 3],
            Relative_Slope = NA,
            Unknown_BAU = NA, Unknown_CI = NA))
        
        # create plot for unparallel and save
        pdf(paste0("./",lab,"/",lab,"_figures/not_parallel/",samples$Sample_ID[1],"_",samples$Antigen[1],"_",samples$Antibody[1],".pdf"))
        plot(standard_4_not_parallel$dose_log, standard_4_not_parallel$result_log,
             col = "red", pch = 1,
             ylim = c(min(together_4_not_parallel$result_log), max(together_4_not_parallel$result_log)),
             ylab = "Result (NOT PARALLEL)",
             xaxt = "n", xlab = "Dilution",
             main = paste0(lab," lab\n","Sample_ID = ",samples$Sample_ID[1]," NOT PARALLEL\nAntigen = ",samples$Antigen[1],"; Antibody = ",samples$Antibody[1]))
        points(samples_4_not_parallel$dose_log, samples_4_not_parallel$result_log,
               col = "blue", pch = 1)
        legend("bottomright",
               legend = c(standard_id, samples_4_not_parallel$Sample_ID[1]),
               col = c("red", "blue"), pch = 1, bty = "y")
        axis(1, at = sort(unique(samples_4_not_parallel$dose_log), decreasing = F),
             labels = changeSciNot(sort(unique(samples_4_not_parallel$dose), decreasing = F)))
        dev.off()
        cat(paste("   combination not parallel\n"),
            file = logfile, append = T)
        break
      }
      
      # remove top or bottom point based on which returns a better relative slope
      who_top_lm = summary(lm(result_log ~ dose_log, standard_top))
      who_bottom_lm = summary(lm(result_log ~ dose_log, standard_bottom))
      sample_top_lm = summary(lm(result_log ~ dose_log, samples_top))
      sample_bottom_lm = summary(lm(result_log ~ dose_log, samples_bottom))
      rel_slope_top = who_top_lm$coefficients[2, 1] / sample_top_lm$coefficients[2, 1]
      rel_slope_bottom = who_bottom_lm$coefficients[2, 1] / sample_bottom_lm$coefficients[2, 1]
      
      if (abs(rel_slope_top - 1) <= abs(rel_slope_bottom - 1)) {
        standard = standard_top
        samples = samples_top
        cat("   removing upper point from standard and sample\n",
            file = logfile, append = T)
      } else {
        standard = standard_bottom
        samples = samples_bottom
        cat("   removing lower point from standard and sample\n",
            file = logfile, append = T)
      }
      
      # combine standard and sample
      together = rbind(standard, samples)
      
      # get slope and intercept for WHO and sample based on log transformed result
      who_lm = summary(lm(result_log ~ dose_log, standard))
      sample_lm = summary(lm(result_log ~ dose_log, samples))
      who_slope = who_lm$coefficients[2, 1]
      sample_slope = sample_lm$coefficients[2, 1]
      who_inter = who_lm$coefficients[1, 1]
      sample_inter = sample_lm$coefficients[1, 1]
      who_r2 = who_lm$r.squared
      sample_r2 = sample_lm$r.squared
      rel_slope = who_slope / sample_slope
      cat(paste("   standard slope after removing points:", who_slope, "\n",
                "   standard intercept after removing points:", who_inter, "\n",
                "   standard R2 after removing points:", who_r2, "\n",
                "   sample slope after removing points:", sample_slope, "\n",
                "   sample intercept after removing points:", sample_inter, "\n",
                "   sample R2 after removing points:", sample_r2, "\n",
                "   relative slope after removing points:", rel_slope, "\n"),
          file = logfile, append = T)
      
      # if sample is within 20% of slope, plot
      if (rel_slope >= 0.8 & rel_slope <= 1.2) {
        cat("   relative slope is within 20% after removing points, proceed to calculating BAU and plotting\n",
            file = logfile, append = T)        
        # calculate mean and CI for relative potency using bootstrapping
        co = c()
        set.seed(1)
        for (j in 1:100) {
          b_df_standard = data.frame(result_log = c(), dose_log = c())
          b_df_samples = data.frame(result_log = c(), dose_log = c())
          for (k in unique(standard$dose_log)) {
            b_all_standard = subset(standard, dose_log == k)
            b_standard = sample(b_all_standard$result_log, 3, replace = T)
            b_sub_standard = data.frame(result_log = b_standard, dose_log = k)
            b_df_standard = rbind(b_df_standard, b_sub_standard)
          }
          for (k in unique(samples$dose_log)) {
            b_all_samples = subset(samples, dose_log == k)
            b_samples = sample(b_all_samples$result_log, 3, replace = T)
            b_sub_samples = data.frame(result_log = b_samples, dose_log = k)
            b_df_samples = rbind(b_df_samples, b_sub_samples)
          }
          b_who_lm = summary(lm(result_log ~ dose_log, b_df_standard))
          b_sample_lm = summary(lm(result_log ~ dose_log, b_df_samples))
          b_who_slope = b_who_lm$coefficients[2, 1]
          b_sample_slope = b_sample_lm$coefficients[2, 1]
          b_who_inter = b_who_lm$coefficients[1, 1]
          b_sample_inter = b_sample_lm$coefficients[1, 1]
          
          interval_result = seq(
            from = min(b_df_standard$result_log),
            to = max(b_df_standard$result_log),
            length.out = 100
          )
          
          who_potency = (interval_result - b_who_inter) / b_who_slope
          sample_potency = (interval_result - b_sample_inter) / b_sample_slope
          log_rel_potency = who_potency - sample_potency
          mean_rel_pot = mean(10 ^ log_rel_potency) * 1000
          co[j] = mean_rel_pot
        }
        m = round(mean(co))
        e = round(qnorm(0.975) * sd(co) / sqrt(100))
        cil = m - e
        ciu = m + e
        cat(paste("   BAU:", m, "\n",
                  "   CI:", e, "\n"),
            file = logfile, append = T)
        df_output = rbind(
          df_output,
          data.frame(
            Sample_ID = u_comb[i, 1], Antigen = u_comb[i, 2], Antibody = u_comb[i, 3],
            Relative_Slope = rel_slope,
            Unknown_BAU = m, Unknown_CI = e))
        
        # create plot and save
        cat("   create plot\n",
            file = logfile, append = T)
        pdf(paste0("./",lab,"/",lab,"_figures/",samples$Sample_ID[1],"_",samples$Antigen[1],"_",samples$Antibody[1],".pdf"))
        plot(standard$dose_log, standard$result_log,
             col = "red", pch = 1,
             ylim = c(min(together$result_log), max(together$result_log)),
             ylab = paste0("Result (log transformed and dropped ", num_dilutions - length(unique(samples$dose_log)), " dilutions)"),
             xaxt = "n", xlab = "Dilution",
             main = paste0(lab," lab\n","Sample_ID = ",samples$Sample_ID[1],"\nAntigen = ",samples$Antigen[1],"; Antibody = ",samples$Antibody[1]))
        points(samples$dose_log, samples$result_log,
               col = "blue", pch = 1)
        abline(b = who_slope, a = who_inter, col = "red")
        abline(b = sample_slope, a = sample_inter, col = "blue")
        legend("bottomright",
               legend = c(standard_id, samples$Sample_ID[1]),
               col = c("red", "blue"), pch = 1, bty = "y")
        legend("topleft",
               legend = c(paste0("BAU = ", m), paste0("Conf. Int. = [", cil, ", ", ciu, "]")),
               bty = "n")
        axis(1, at = sort(unique(samples$dose_log), decreasing = F),
             labels = changeSciNot(sort(unique(samples$dose), decreasing = F)))
        dev.off()
        break
        
      } else {
        cat("   not possible to make combination parallel\n",
            file = logfile, append = T)
        next
      }
    }
  }
}
colnames(df_output) = c("Sample_ID", "Antigen", "Antibody",
                        "Relative_Slope",
                        "Unknown_BAU", "Unknown_CI")
write.csv(df_output, paste0("./", lab, "/", lab, "_BAU_table.csv"))
cat("\nsave BAUs\nend script\n",
    file = logfile, append = T)
###################################

###################################
# Uncomment following code if all bau tables produced should be combined into one table

# setwd(path)
# bau_files = list.files(pattern = "BAU_table.csv$", recursive = T)
# all_bau = data.frame(Lab_Name = c(),
#                      X = c(),
#                      Sample_ID = c(),
#                      Antigen = c(),
#                      Antibody = c(),
#                      Relative_Slope = c(),
#                      Unknown_BAU = c(),
#                      Unknown_CI = c())
# for (i in bau_files) {
#   one_bau_table = read.csv(i)
#   one_bau_table = cbind(str_split(i, "/")[[1]][1], one_bau_table)
#   all_bau = rbind(all_bau, one_bau_table)
# }
# colnames(all_bau) = c("Lab_Name", "Unique_Combo",
#                       "Sample_ID",
#                       "Antigen", "Antibody",
#                       "Relative_Slope",
#                       "Unknown_BAU", "Unknown_CI")
# write.csv(all_bau, "pla_all_bau.csv")
