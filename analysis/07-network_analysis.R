
# load packages and read data ---------------------------------------------
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)

gene_module_assignment <- read.csv("data/processed/network/gene_module_assignment.csv")
module_summary <- read.csv("data/processed/network/module_summary.csv")
cell_eigengenes <- read.csv("data/processed/network/cell_eigengenes.csv")

ring_genes <- read.csv("data/processed/network/ring_genes_0610.csv")

# module eigengene + trajectories plot ------------------------------------

name1=paste0("Module ",1:8)
names(name1)=paste0("eigen",1:8)

name2=paste0("Module ",9:16)
names(name2)=paste0("eigen",9:16)


p1 <- cell_eigengenes %>%
  pivot_longer(cols = eigen1:eigen8, names_to="module", values_to="expr") %>%
  group_by(module) %>%
  mutate(expr_scaled=(expr-mean(expr))/sd(expr)) %>%
  mutate(expr_scaled=ifelse(expr_scaled>6, 6,expr_scaled)) %>%
  arrange(expr_scaled) %>%
  ggplot(aes(V1, V2, colour = expr_scaled)) +
  geom_point() +
  facet_grid(rows = vars(module),labeller = labeller(module=name1),
             switch = "both")+
  scale_colour_viridis_c() +
  labs(x = "", y = "",
       colour = "Expr") +
  ggtitle("Module eigengenes") +
  theme_void() +
  theme(plot.title = element_text(face="bold",vjust =3,hjust=0.6),
        legend.position = 'bottom',
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size=5),
        legend.margin =margin(r=10,l=5,t=5,b=5)) +
  coord_fixed()


p2 <- cell_eigengenes %>%
  pivot_longer(cols = eigen1:eigen8, names_to="module", values_to="expr") %>%
  filter(traj%in%paste0("Traj",1:3)) %>%
  group_by(traj) %>%
  mutate(expr_scaled=(expr-mean(expr))/sd(expr)) %>%
  ggplot(aes(time,expr_scaled)) +
  geom_point(color="grey")+
  geom_hline(yintercept=0,size=0.5)+
  geom_smooth(aes(col=traj),method = "gam",formula = y ~ s(x),size=1.4) +
  scale_x_continuous(breaks=c(0,0.1,0.2))+
  xlab("pseudotime")+
  ylab("expression")+
  facet_grid(rows=vars(module), cols = vars(traj),
             labeller = labeller(module=name1, traj=c(Traj1="PPP", Traj2="CC", Traj3="PSE"))) +
  scale_color_brewer(palette = "Dark2",name = "Trajectory", labels = c("PPP", "CC", "PSE"))+
  ggtitle("Trajectories")+
  theme_minimal()+
  theme(plot.title = element_text(face="bold"),
        strip.text.x = element_text(size = 10),
        legend.position = 'bottom',
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size=5),
        legend.margin =margin(r=10,l=5,t=5,b=5))


p3 <- cell_eigengenes %>%
  pivot_longer(cols = eigen9:eigen16, names_to="module", values_to="expr") %>%
  group_by(module) %>%
  mutate(expr_scaled=(expr-mean(expr))/sd(expr)) %>%
  mutate(expr_scaled=ifelse(expr_scaled>6, 6,expr_scaled)) %>%
  arrange(expr_scaled) %>%
  ggplot(aes(V1, V2, colour = expr_scaled)) +
  geom_point() +
  facet_grid(rows = vars(module),labeller = labeller(module=name2),
             switch = "both")+
  scale_colour_viridis_c() +
  labs(x = "", y = "",
       colour = "Expr") +
  ggtitle("Module eigengenes") +
  theme_void() +
  theme(plot.title = element_text(face="bold",vjust =3,hjust=0.6),
        legend.position = 'bottom',
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size=5),
        legend.margin =margin(r=10,l=5,t=5,b=5)) +
  coord_fixed()


p4 <- cell_eigengenes %>%
  pivot_longer(cols = eigen9:eigen16, names_to="module", values_to="expr") %>%
  filter(traj%in%paste0("Traj",1:3)) %>%
  group_by(traj) %>%
  mutate(expr_scaled=(expr-mean(expr))/sd(expr)) %>%
  ggplot(aes(time,expr_scaled)) +
  geom_point(color="grey")+
  geom_hline(yintercept=0,size=0.5)+
  geom_smooth(aes(col=traj),method = "gam",formula = y ~ s(x),size=1.4) +
  scale_x_continuous(breaks=c(0,0.1,0.2))+
  xlab("pseudotime")+
  ylab("expression")+
  facet_grid(rows=vars(module), cols = vars(traj),
             labeller = labeller(module=name2, traj=c(Traj1="PPP", Traj2="CC", Traj3="PSE"))) +
  scale_color_brewer(palette = "Dark2",name = "Trajectory", labels = c("PPP", "CC", "PSE"))+
  ggtitle("Trajectories")+
  theme_minimal()+
  theme(plot.title = element_text(face="bold"),
        strip.text.x = element_text(size = 10),
        legend.position = 'bottom',
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size=5),
        legend.margin =margin(r=10,l=5,t=5,b=5))


png("module eigengenes and trajectories 1_8.png",
    width = 2480, height = 3508, units = 'px',
    pointsize = 12, res = 300)
p1+p2+plot_layout(ncol = 2, widths = c(1,1.2))
dev.off()

png("module eigengenes and trajectories 9_16.png",
    width = 2480, height = 3508, units = 'px',
    pointsize = 12, res = 300)
p3+p4+plot_layout(ncol = 2, widths = c(1,1.2))
dev.off()


# sub-module eigenegenes of Module 1 --------------------------------------

name=paste0("sub-module ",1:15)
names(name)=paste0("eigen1.",1:15)

png("Module 1 sub-module eigengenes.png",
    width = 2480, height = 3508, units = 'px',
    pointsize = 12, res = 300)

cell_eigengenes %>%
  pivot_longer(cols = eigen1.1:eigen1.15, names_to="module", values_to="expr") %>%
  group_by(module) %>%
  mutate(expr_scaled=(expr-mean(expr))/sd(expr)) %>%
  mutate(expr_scaled=ifelse(expr_scaled>6, 6,expr_scaled)) %>%
  arrange(expr_scaled) %>%
  mutate(module=factor(module,levels=c(paste0("eigen1.",1:9),paste0("eigen1.",10:15)))) %>%
  ggplot(aes(V1, V2, colour = expr_scaled)) +
  geom_point() +
  facet_wrap(vars(module),ncol = 3, nrow=5, shrink = FALSE,
             labeller = labeller(module=name))+
  scale_colour_viridis_c() +
  labs(x = "", y = "",
       colour = "Expr") +
  ggtitle("Module 1: sub-module eigengenes")  +
  theme_test()+
  theme(plot.title = element_text(face="bold",vjust =2),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size=5),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  coord_fixed()

dev.off()


# p-values ----------------------------------------------------------------

ring_genes <- ring_genes %>%
  left_join(gene_module_assignment, by="gene")
write.csv(ring_genes,"data/processed/network/ring_genes_0610.csv")

message("Number of ring genes in Module 1 is ",sum(ring_genes$module=="M1",na.rm=TRUE))
message("Number of ring genes in sub-module 1 is ",sum(ring_genes$sub_module=="M1.1",na.rm=TRUE))

# 8 ring genes fall in Module 1
#total genes in network
N=sum(module_summary$number_genes);
# module 1 size
m=module_summary[1,"number_genes"];
n=N-m;
# ring genes in network
k=nrow(ring_genes);
# observed in Module1
x=sum(ring_genes$module=="M1",na.rm=TRUE);

# Probability for 8 or more out of 10 ring genes observed in Module 1
p=phyper(x-1,m,n,k,lower.tail = FALSE)


# 7 ring genes fall in sub-module 1 of Module 1
#total genes in Module 1
N1=module_summary[1,"number_genes"];
# module 1 sub-module 1
m1=sum(gene_module_assignment$sub_module=="M1.1",na.rm = TRUE);
n1=N1-m1;
# ring genes in Module 1
k1=sum(ring_genes$module=="M1",na.rm=TRUE);
# observed ring genes in module 1 sub 3
x1=sum(ring_genes$sub_module=="M1.1",na.rm=TRUE);

# Probability for 7 out of 8 ring genes observed in sub-module 3
p1=phyper(x1-1,m1,n1,k1,lower.tail = FALSE)

