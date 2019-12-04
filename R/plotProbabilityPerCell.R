# prob_plot <- ggplot(data=r_sub[r_sub$stat=='mean',], aes(x=cells, y=probs, colour=gene))+
#   # geom_point(size=1) +
#   geom_line(aes(linetype=cell_type), size=2) +
#   # geom_errorbar(aes(ymin=meanminstd_prob, ymax=meanplusstd_prob),width=.1,
#                 # position=position_dodge(0.05)) +
#   scale_linetype_manual(labels = c("iPS","Liver","S.I."), values=c( "solid", "dotted", "dashed")) +
#   scale_color_manual(labels = c("Oncogenic mutations","TP53","BRAF.V600E"), values = c('#F8766D','#619CFF','#00BA38')) + #green: #00BA38, blue: #619CFF
#   # geom_smooth(method='loess', formula = y ~ exp(x), se=F, span=0.15, size=2) +
#   # geom_abline(aes(slope=1e-3, intercept = 0)) +
#   labs(x='Number of cells', y=ylab_probs, 
#        color='Type of mutations',
#        # shape='Cell types',
#        linetype='Cell types') 

# if (any(vline_prob!=F)){
#   for (i in vline_prob){
#     prob_plot <- prob_plot + geom_vline(xintercept=as.numeric(i), linetype='dotted') 
#       # geom_text(aes( as.numeric(i), 0, label = format(as.numeric(i), scientific = T), vjust = -1), size = 8, colour="black")
#   }
# }
# if (any(hline_prob!=F)){
#   for (i in hline_prob){
#     prob_plot <- prob_plot +  geom_hline(yintercept=as.numeric(i), linetype='dotted') 
#   }
# }


# if (log_y_mut){
#   if (log_x_mut){
#     mut_plot <- mut_plot + scale_y_log10(breaks=breaks_y, labels=labels_y) + scale_x_log10(breaks=breaks_x,labels=labels_x)
#     # mut_plot <- mut_plot + scale_y_log10() + scale_x_log10(breaks=breaks_x,labels=labels_x)
#     # mut_plot <- mut_plot + scale_y_log10() + scale_x_log10()
#   } else {
#     mut_plot <- mut_plot + scale_y_log10()
#   }
# } else {
#   if (log_x_mut){
#     mut_plot <- mut_plot + scale_x_log10()
#   }
# }
# 
# if (log_y_prob){
#   if (log_x_prob){
#     prob_plot <- prob_plot + scale_y_log10() + scale_x_log10()
#   } else {
#     prob_plot <- prob_plot + scale_y_log10()
#   }
# } else {
#   if (log_x_prob){
#     prob_plot <- prob_plot + scale_x_log10()
#   }
# }

# prob_plot <- prob_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                              panel.background = element_blank(), axis.line = element_line(colour = "black"),
#                              text = element_text(size=30))

