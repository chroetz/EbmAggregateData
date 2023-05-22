rootDir <- "/p"
varNames <- c("hurs", "huss", "pr", "prsn", "ps", "rlds", "rsds", "sfcwind", "tas", "tasmax", "tasmin")
for (varName in varNames) {
	jobName <- paste0("stats_gswp3-w5w5_year_country_", varName, "_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
	cat("Starting SLURM job", jobName, "\n")
	clcom <- paste0(
	  "sbatch",
	  " --qos=standby",
	  " --time=24:00:00",
	  " --ntasks=1",
	  " --cpus-per-task=4", # RAM: ~ cpus-per-task * 4 GB
	  " --job-name=", jobName,
	  " --output=", jobName, "_%j.out",
	  " --error=", jobName, "_%j.err",
	  " --mail-type=END",
	  " --wrap=\"",
	  paste(
		"Rscript stats_gswp3-w5w5_year_runCountry.R", varName),
	  "\"")
	cat(clcom, "\n")
	system(clcom)
}