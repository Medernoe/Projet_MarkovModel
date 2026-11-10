library(seqinr)
library(ggplot2)
library(tidyr)
library(HMM)

setwd('M2_BIMS/Markov/Projet')
source("fonctions.R")

# ========= Ouverture Data =========

# --- Fichier --- 
# fichier FASTA
f <- read.fasta("E_coli/GCF_000005845.2_ASM584v2_genomic.fna")
seq <- f$NC_000913.3
n <- length(seq)

# fichier GFF
gff <- read.delim("E_coli/genomic.gff", comment.char="#", header=FALSE)
colnames(gff) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")



# --- Sequence etat cache ---
# Etat caches : Annotation codant non codant 
# Etat observer : acide nucleique 
# On connait tout du genome, celui ci est utiliser a des fins de test pour des modelisations
# Les annotations sont en partie cache et servent au calcul de precision 

# Isoler les régions codantes (CDS)
cds_gff <- gff[gff$feature == "CDS", ]

# Sequence de région codante / non codante 
# 1 non codants
# 2 codants 
etats <- rep("1", n)
# Remplissage de l'état "2" en fonction des coordonnées du GFF
for (i in 1:nrow(cds_gff)) {
  start_pos <- cds_gff$start[i]
  end_pos <- cds_gff$end[i]
  if (start_pos <= n && end_pos <= n) {
    etats[start_pos:end_pos] <- "2"
  }
}

# pourcentages de séquences codante 
(sum(etats == "2") / length(etats)) * 100


# ========= Partie A =========

# --- calcul Log-vraisemble, AIC, BIC ---
k_max <- 7

# tableau de résultat
resultats_modeles <- data.frame(
  Ordre = 0:k_max,
  Log_Vraisemblance = numeric(k_max + 1),
  AIC = numeric(k_max + 1),
  BIC = numeric(k_max + 1)
)

# Calcul des métriques pour chaque ordre (AIC, BIC, Log-Vraisemblance)
for (i in 0:k_max) {
  # R indexe à partir de 1, donc on utilise i + 1 pour l'index
  idx <- i + 1 
  
  # Calcul de la log-vraisemblance
  log_v <- log_vraisemblance(seq, ordre = i)
  
  resultats_modeles$Log_Vraisemblance[idx] <- log_v
  resultats_modeles$AIC[idx] <- AIC(log_vrais = log_v, ordre = i)
  resultats_modeles$BIC[idx] <- BIC(log_vrais = log_v, ordre = i, seq = seq)
}

print(resultats_modeles)


# --- Comparaison par Rapport de Vraisemblance ---

# Tableau de résultats : modèles M(i) vs M(i+1)
rapports_tests <- data.frame(
  Comparaison = character(k_max),
  Statistique = numeric(k_max),
  Quantile = numeric(k_max),
  DDL = numeric(k_max),
  P_Value = numeric(k_max)
)

# Boucle pour tester les modèles successifs
for (i in 0:(k_max - 1)) {
  idx <- i + 1
  
  # Récupération des log-vraisemblances depuis le tableau précédent
  L0 <- resultats_modeles$Log_Vraisemblance[idx]     # Modèle ordre i
  L1 <- resultats_modeles$Log_Vraisemblance[idx + 1] # Modèle ordre i+1
  
  # Lancement du test
  test <- rapport_vraisemblance(
    log_vrais0 = L0, 
    log_vrais1 = L1, 
    ordre0 = i, 
    ordre1 = i + 1
  )
  
  rapports_tests$Comparaison[idx] <- paste("Ordre", i, "vs", i + 1)
  rapports_tests$Statistique[idx] <- test$statistique
  rapports_tests$Quantile[idx] <- test$quantile
  rapports_tests$DDL[idx] <- test$ddl
  rapports_tests$P_Value[idx] <- test$pvalue
}

print(rapports_tests)

# --- Plot ---
# plot Log-Vraisemblance
plot_logV <- ggplot(resultats_modeles, aes(x = Ordre, y = Log_Vraisemblance)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue", size = 3) +
  theme_minimal() +
  labs(title = "Évolution de la Log-Vraisemblance en fonction de l'ordre k du modèle",
       x = "Ordre du modèle de Markov",
       y = "Log-Vraisemblance")

print(plot_logV)


# plot comparatif : AIC vs BIC
plot_criteres <- ggplot(resultats_modeles, aes(x = Ordre)) +
  # AIC
  geom_line(aes(y = AIC, color = "AIC"), linewidth = 1) +
  geom_point(aes(y = AIC, color = "AIC"), size = 3) +
  
  # BIC
  geom_line(aes(y = BIC, color = "BIC"), linewidth = 1) +
  geom_point(aes(y = BIC, color = "BIC"), size = 3) +
  
  # Lignes verticales 
  geom_vline(xintercept = 6, linetype = "dashed", color = "#2ca02c", alpha = 0.5) + 
  geom_vline(xintercept = 5, linetype = "dashed", color = "#d62728", alpha = 0.5) +
  
  # Légende
  theme_minimal() +
  labs(title = "Comparaison des critères AIC et BIC",
       x = "Ordre du modèle de Markov",
       y = "Valeur du critère",
       color = "Critère") +
  scale_color_manual(values = c("AIC" = "#2ca02c", "BIC" = "#d62728")) 

print(plot_criteres)


# plot comparatif : Statistique de test du rapport de Vraisemblance
rapports_tests$Ordre_Test <- 0:(k_max-1)

plot_rapport <- ggplot(rapports_tests, aes(x = Comparaison, y = Statistique, group = 1)) +
  geom_line(color = "purple", linewidth = 1) +
  geom_point(color = "purple", size = 3) +
  theme_minimal() +
  labs(title = "Statistique du Rapport de Vraisemblance",
       subtitle = "Comparaison d'un modèle d'ordre i avec l'ordre i+1",
       x = "Ordre (i) du modèle de base",
       y = "Valeur de la statistique (Chi-2)")

print(plot_rapport)










# ========= Partie B =========
# Modeliser une sequence d'etats caches a partir d'observation pour la prediction 

# --- Data hold out (Version Aléatoire Sécurisée) --- 
taille <- 500000
n <- length(seq)

#  indice Train set
set.seed(1234)
start_idx_train <- sample(1:(n - taille + 1), 1) 
idx_train <- start_idx_train:(start_idx_train + taille - 1)

# indice Validation set
# Il est crucial d'exclure la zone utilisée par le Train set pour éviter tout chevauchement entre Train et Validation.
zone_a_exclure <- (start_idx_train - taille):(start_idx_train + taille)
indices_possibles_valid <- setdiff(1:(n - taille + 1), zone_a_exclure)

set.seed(5678)
start_idx_valid <- sample(indices_possibles_valid, 1)
idx_valid <- start_idx_valid:(start_idx_valid + taille - 1)

# Création des jeux de donnes 
seq_train <- seq[idx_train]
etats_train <- etats[idx_train] 

seq_valid <- seq[idx_valid]
etats_valid <- etats[idx_valid]


# --- EM ---
# estimer les parametres 
# Définition : états, alphabet, matrice de transition et matrice d'émission
etats_caches <- c("1", "2")
alphabet <- c("a", "c", "g", "t")

table(etats)/n*100 # etats 1 = non codant / 2 = codant
#1        2 
#13.73976 86.26024 

table(seq)/n*100 # pourcentage de lettre 
#       a        c        g        t 
#24.61930 25.42394 25.36677 24.59000 


# loi initial 
# Proche du pourcentage de codant / non codant 
mu = c(0.86, 0.14)
names(mu) <- etats_caches


# transition
#L'algorithme EM a besoin de conditions initiales asymétriques pour ne pas rester bloqué dans un optimum local symétrique.
# les proba. des matrices ont ete mise de facon subjectives 

transition <- matrix(c(0.6, 0.4, 
                       0.3, 0.7), 
                     nrow = 2, 
                     byrow = TRUE)
rownames(transition) <- etats_caches
colnames(transition) <- etats_caches

# Emision 
emission <- matrix(c(0.26, 0.24, 0.24, 0.26,
                     0.2, 0.3, 0.3, 0.2), 
                   nrow = 2,
                   byrow = TRUE)
rownames(emission) <- etats_caches
colnames(emission) <- alphabet

# Initialisation du modèle HMM (M1-M0)
hmm_init <- initHMM(
  States = etats_caches,
  Symbols = alphabet,
  startProbs = mu,
  transProbs = transition,
  emissionProbs = emission
)

# Algo baumWelch (EM)
hmm_fit <- baumWelch(hmm_init, 
                     observation = seq_train, 
                     maxIterations = 100)

print(hmm_fit$hmm)


# --- Viterbi ---

# Verifier les resultats et l'efficacite de la re implementation de l'algo de viterbi par rapport au package HMM
# print(a)
# Segmentation_package <- HMM::viterbi(hmm_fit$hmm, seq_train)
# z=Sys.time()
# print(z-a)

# a=Sys.time()
# Segmentation_maison <- viterbi(seq_train, hmm_fit$hmm$startProbs, hmm_fit$hmm$transProbs, hmm_fit$hmm$emissionProbs,log = TRUE)
# z=Sys.time()
# print(z-a)

# print(identical(as.character(Segmentation_package), as.character(Segmentation_maison$chemin)))


# Viterbi
# Estime la succession d'etats caches la plus vraisemblable 
# on applique sur seq valid car on cherche a voir la precision de nos prediction des estimateurs via EM 
# Si M1-M0 suffit pour modeliser une region codante / non codantes alors on aura un bon score d'accuracy et surtout de mcc
res_viterbi <- viterbi(seq_valid, hmm_fit$hmm$startProbs, hmm_fit$hmm$transProbs, hmm_fit$hmm$emissionProbs, log = TRUE)
etats_predits <- res_viterbi$chemins

TP <- sum(etats_predits == "2" & etats_valid == "2")
TN <- sum(etats_predits == "1" & etats_valid == "1")
FP <- sum(etats_predits == "2" & etats_valid == "1")
FN <- sum(etats_predits == "1" & etats_valid == "2")

# Score de precision
accuracy <- (TP + TN) / length(etats_valid)
print(accuracy)
# 0.655252

# Score MCC 
# Score binaryaccuracy non pertinent du au desequilibre de classe (predire uniquement codant donnerai un tres bon score >86%)
# Score entre -1 et 1, avec proche de 0 = prediction aleatoire
mcc <- (as.numeric(TP) * TN - as.numeric(FP) * FN) / sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) * as.numeric(TN + FP) * as.numeric(TN + FN))
print(mcc)
#  0.1069814


# plot de visualisation de la segmentation
plot_segmentation <- function(sequence_caracteres, taille = 2000) {
  
  # Convertir la séquence de caractères en nombres entiers
  vecteur_binaire <- as.integer(sequence_caracteres)
  vecteur_tronque <- head(vecteur_binaire, n = taille)
  df <- data.frame(
    Index = 1:length(vecteur_tronque),
    Valeur = vecteur_tronque
  )
  
  # plot
  ggplot(df, aes(x = Index, y = Valeur)) +
    geom_step(color = "red", direction = "hv") +
    theme_minimal() +                           
    labs(x = "Position", y = "Valeur")                       
}
seg_predite<-plot_segmentation(etats_predits)

seg_valid<-plot_segmentation(etats_valid)



# ========= Export des Figures =========

# Création d'un sous-dossier "figures" s'il n'existe pas déjà
dir.create("figures", showWarnings = FALSE)

# Évolution de la Log-Vraisemblance
ggsave(filename = "figures/1_Log_Vraisemblance.png", 
       plot = plot_logV, 
       width = 8, height = 6, dpi = 300, bg = "white")

# Comparaison AIC vs BIC
ggsave(filename = "figures/2_Criteres_AIC_BIC.png", 
       plot = plot_criteres, 
       width = 8, height = 6, dpi = 300, bg = "white")

# Rapport de Vraisemblance
ggsave(filename = "figures/3_Rapport_Vraisemblance.png", 
       plot = plot_rapport, 
       width = 8, height = 6, dpi = 300, bg = "white")

# Visualisation de la segmentation (Viterbi vs Réalité)
# Format plus "large" (width = 10, height = 4) adapté aux séquences
ggsave(filename = "figures/4a_Segmentation_Predite.png", 
       plot = seg_predite, 
       width = 10, height = 4, dpi = 300, bg = "white")

ggsave(filename = "figures/4b_Segmentation_Reelle.png", 
       plot = seg_valid, 
       width = 10, height = 4, dpi = 300, bg = "white")