fd@keep.names <- FALSE
m <- melt(fd)
ll <- lm(value ~ Population*Stim.Condition, m)
w.var <- attr(effects(ll), 'assign')

lrout2 <- zlm.SingleCellAssay(value ~ Population*Stim.Condition, fd)
