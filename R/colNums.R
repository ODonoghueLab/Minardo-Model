cols_matFifty = list(col_clus = 1, col_x = 2, col_y = 3, col_dir = 4)
cols_matFifty_v2 = list(clus = 1, x = 2, y = 3, dir = 4, startTp=5, endTp = 6)


cols_clusPlotObjs = list(col_x0 = 1, col_y0 =2, col_x1 = 3, col_y1 = 4, col_col = 5)

cols_grayLines_v2 = list(col_x0=1, col_x1=2, col_isSemiCirc=3, col_col=4)


cols_timeReg <- list(cluster=1, zScore=2, tpStart=3, tpEnd=4, dir=5)

Col_labels <- list(label=1, color=2)

Col_events <- list(clus=1, x=2, y=3, dir=4, startTp=5, endTp=6, order=7, combinedDatasetNum=8)

Color_multiomics <- list(list(incr="#ff0000", decr="#f3a956", bar="#ffc2b5"), list(incr="#712cdf", decr="#b5a1de", bar="#c4b1ea"), list(incr="#4baba4", decr="#8cc7c7", bar="#7ae6c3"))


colors_orderedEvents = list(incr = "#e6aeae", decr=rgb(0,0,0,alpha=0))

colors_orderedEventsNew = list(incr_phos = "#e6aeae", incr_prot = "#89d8d0", incr_rnaSeq="#d7aeff", decr=rgb(0,0,0,alpha=0)) # incr color #F78A72
#
color_events = list(up="#ea3424", down="#C99999")
color_eventDist = list(up=rgb(1, 0, 0, 0.5), down=rgb(0, 0, 1, 0.5))

# colors_orderedEvents <- vector(mode="list", length=3)
# names(colors_orderedEvents) <- c("incr", "decr")
# colors_orderedEvents[[1]] <- "#989898" # incr - dark gray
# colors_orderedEvents[[2]] <- "#DCDCDC" # decr - light gray

# red increased - #ea3424

# background boxes - #ededed


eventOrderTest = list(param="t-test", nonParam="wilcox")

################## OTHER CONSTANTS

cols_rectPoints = list(x0=1, x1=2)
cols_missingStats = list(numNa=1, numTotal=2, percentNa=3)
