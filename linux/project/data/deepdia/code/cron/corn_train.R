library(cronR)
f="/srv/shiny-server/project/data/deepdia/code/loop_train.R"
cmd <- cron_rscript(f)
cron_add(command = cmd, frequency = 'minutely', id = 'train', description = 'train', tags = c('lab', 'xyz'))
cron_add(command = cmd, frequency = 'hourly', id = 'train', description = 'train', tags = c('lab', 'xyz'))

# cron_ls()
# cron_clear(ask=FALSE)
