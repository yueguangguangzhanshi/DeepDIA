library(cronR)
f="/home/saitoasuka/shiny-server/project/data/deepdia/code/loop_train.R"
cmd <- cron_rscript(f)
cron_add(command = cmd, frequency = 'hourly', id = 'train', description = 'train', tags = c('lab', 'xyz'))
