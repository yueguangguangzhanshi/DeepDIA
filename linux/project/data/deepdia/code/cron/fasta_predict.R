f="/home/saitoasuka/shiny-server/project/data/deepdia/code/loop_predict.R"
cmd <- cron_rscript(f)
cron_add(command = cmd, frequency = 'hourly', id = 'fasta predict', description = 'fasta predict', tags = c('lab', 'xyz'))
