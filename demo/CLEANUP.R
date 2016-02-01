rm(DATA_ENV, ROOT, SEED, REFRESH)
ours <- unlist(eapply(.GlobalEnv, FUN = function(x) is(x, "stanreg") | is(x, "loo")))
ours <- names(ours[ours])
rm(list = ours)
rm(ours)

