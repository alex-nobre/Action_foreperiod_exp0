textwidth = 17.5 # cm

# sizes (line, points)
lsz = 0.75
psz = 1.5

# color/shape scales, for consistency across graphs
shpscale <- scale_shape_manual(values=c(21,24))
colcat <- scale_color_manual(values=c('#A35C00', '#00523D'))
fillcat<- scale_fill_manual(values=c('#A35C00','#00523D')) 

colbetas  <- scale_color_manual(values=c( '#0C1618', '#154D66', '#8E1529', '#490838'))
fillbetas <- scale_fill_manual(values=c( '#0C1618', '#154D66', '#8E1529', '#490838'))

# group scales into a list
myscales <- list(shpscale, colcat, fillcat)
# define a 'theme'
mytheme	 <- theme_classic() + theme(
  legend.position = c(1,1), legend.justification = c(1,1), legend.direction = 'horizontal', 
  panel.grid.major.y=element_line(linetype='dotted',color='black',linewidth=0.25*lsz))