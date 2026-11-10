library(seqinr)

# ========= Estimateur =========
# Estimateur matrice de transition Mk
# Avec Estimateur Pab_g_est = Nabg/Nab. pour tout a,b,g € A
estimateur_transition <- function(sequence, ordre = 1, alphabet = c('a','c','g','t')){

  # cas ordre 0
  if(ordre == 0){
    #compte les occurences de mots de 1 lettre
    comp1 <- count(sequence,1) 
    mat = comp1/sum(comp1)
    return(mat)
  }
  
  # Comptage des mots de taille (ordre + 1) et (ordre)
  compkk <- count(sequence, ordre + 1)
  compk <- count(sequence, ordre)
  
  # Récupération des mots de départ (k-mers) pour les colonnes de la matrice
  ab <- names(compk)
  
  # Initialisation de la matrice de transition
  mat <- matrix(0, nrow = length(alphabet), ncol = length(ab))
  rownames(mat) <- alphabet
  colnames(mat) <- ab
  
  # Remplissage de la matrice
  for (combi in ab) {
    
    # Extraire les comptages associés à ces mots de k+1 lettres
    abg <- paste0(combi, alphabet)
    Nabg <- compkk[abg]
    
    # Remplacer les NA (les combinaisons non observées) par 0
    Nabg[is.na(Nabg)] <- 0
    
    # Nombre total d'observations partant de la combinaison 'combi'
    Nab <- sum(Nabg)
    
    # Calcul des probabilités de transition (Estimateur de maximum de vraisemblance)
    if (Nab > 0) {
      mat[,combi] <- Nabg / Nab
    }
  }
  
  return(mat)
}


# ========= Generation =========
generer_sequence_HMM <- function(n, mu, transition, emission) {
  etats <- c()
  obs <- c()
  etats_possibles <- rownames(transition)
  alphabet <- colnames(emission)
  
  # Initialisation
  etats[1] <- sample(etats_possibles, 1, prob = mu)
  obs[1] <- sample(alphabet, 1, prob = emission[etats[1], ])
  
  # Récurrence
  for (t in 2:n) {
    etats[t] <- sample(etats_possibles, 1, prob = transition[etats[t-1], ])
    obs[t] <- sample(alphabet, 1, prob = emission[etats[t], ])
  }
  return(obs)
}


# ========= log_vraisemblance =========
# log(Pi) * Somme(Nabg*log(Pabg_est))
log_vraisemblance <- function(sequence, ordre = 1, proba_init = NULL, alphabet = c('a','c','g','t')){
  
  # estimation des paramètre à partir d'une séquence 
  matrice_transition <- estimateur_transition(sequence = sequence, ordre = ordre)
  
  # cas ordre 0
  if(ordre == 0){
    #calcule de la log vraisemblance 
    log_vrais <- sum(count(sequence, 1) * log(matrice_transition))
    return(log_vrais)
  }
  
  # Comptage des mots de taille (k + 1)
  compk <- count(sequence, ordre + 1)
  
  # Log proba de la condition initiale
  if (!is.null(proba_init)) {
    log_vrais <- log(proba_init[sequence[1]])
  } else {
    log_vrais <- 0 
  }
  
  # Récupération des mots observés dans la séquence
  mots <- names(compk)
  
  # Somme
  for (mot in mots) {
    
    # Séparation de l'état de départ (k-mer) et de l'état d'arrivée (1 lettre)
    kmer_depart <- substr(mot, 1, ordre)
    lettre_arrivee <- substr(mot, ordre + 1, ordre + 1)
    
    # Extraction de la probabilité depuis la matrice
    prob <- matrice_transition[lettre_arrivee, kmer_depart]
    
    # On vérifie que la probabilité est > 0 pour éviter l'erreur log(0) = -Inf
    if (!is.na(prob) && prob > 0) {
      log_vrais <- log_vrais + compk[[mot]] * log(prob)
    }
  }
  return(log_vrais)
}


# ========= Comparaison de modèles =========

# test du rapport de vraisemblance
rapport_vraisemblance <- function(log_vrais0, log_vrais1, 
                                  ordre0, ordre1, 
                                  taille_alphabet = 4, alpha = 0.05) {
  
  # Calcul du nombre de paramètres pour chaque modèle : A^k * (A - 1)
  K0 <- (taille_alphabet^ordre0) * (taille_alphabet - 1)
  K1 <- (taille_alphabet^ordre1) * (taille_alphabet - 1)
  
  # Degrés de liberté 
  ddl <- K1 - K0
  
  # Statistique de test
  stat <- -2 * (log_vrais0 - log_vrais1)
  
  # Quantile ordre alpha
  q <- qchisq(1 - alpha, ddl)
  
  # p-value
  p <- pchisq(stat, ddl, lower.tail = FALSE)
  
  return(list(
    statistique = stat, 
    ddl = ddl, 
    quantile = q, 
    pvalue = p
  ))
}


# AIC 
AIC <- function(log_vrais, ordre, taille_alphabet = 4) {
  
  # Calcul du nombre de parametres 
  K <- (taille_alphabet^ordre) * (taille_alphabet - 1)
  
  # AIC = −2L+2∗K
  return(-2*log_vrais + 2*K)
}


# BIC 
BIC <- function(log_vrais, ordre, seq, taille_alphabet = 4) {
  
  # Calcul du nombre de parametres 
  K <- (taille_alphabet^ordre) * (taille_alphabet - 1)
  
  # BIC=−2L+Klog(n)
  return(-2*log_vrais + K * log(length(seq)))
}


# ========= Viterbi  =========

viterbi <- function(Y, mu, transition, emission, log = TRUE, alphabet = c('a','c','g','t'), etat = c('1','2')) {
  n <- length(Y)
  K <- length(etat)
  
  # Initialisation des structures de résultats
  result <- list(v = list(), to = list())
  
  # Vecteur pour stocker les probabilités courantes à l'instant t
  vt <- c()
  
  if (log) {
    # --- version log ---
    # nécéssaire pour n grand car les proba tendent vers 0, on se retrouve avec un arrondie à 0 pour les probas trop grandes 
    
    # Transformation log 
    log_mu <- log(mu)
    log_trans <- log(transition)
    log_emiss <- log(emission)
    
    # initialisation
    t <- 1 
    emiss_t <- log_emiss[, colnames(log_emiss) == Y[t]]
    
    # vt prend la valeur initiale + émission
    vt <- log_mu + emiss_t
    
    result$v[[t]] <- vt
    result$to[[t]] <- rep(NA, K)
    
    # récurrence 
    while(t < n){
      t <- t + 1
      emiss_t <- log_emiss[, colnames(log_emiss) == Y[t]]
      
      nv_vt <- c()
      nv_to <- c()
      
      for (k in 1:K) {
        chemins <- vt + log_trans[, k]
        
        # On calcule le max et l'index du max
        nv_vt[k] <- emiss_t[k] + max(chemins)
        nv_to[k] <- which.max(chemins)
      }
      
      vt <- nv_vt
      result$v[[t]] <- vt
      result$to[[t]] <- nv_to
    }
    
  } else {
    # --- version classique ---
    # version sans log 
    
    # initialisation 
    t <- 1 
    emiss_t <- emission[, colnames(emission) == Y[t]]
    vt <- mu * emiss_t
    result$v[[t]] <- vt 
    result$to[[t]] <- rep(NA, K)
    
    # récurrence 
    while(t < n){
      t <- t + 1
      emiss_t <- emission[, colnames(emission) == Y[t]]
      
      nv_vt <- c()
      nv_to <- c()
      
      for (k in 1:K) {
        chemins <- vt * transition[, k]
        
        nv_vt[k] <- emiss_t[k] * max(chemins)
        nv_to[k] <- which.max(chemins)
      }
      
      vt <- nv_vt
      result$v[[t]] <- vt
      result$to[[t]] <- nv_to
    }
  }
  
  # --- remonter ---
  chemin_indices <- c()
  chemin_labels <- c()
  
  # Dernier état
  chemin_indices[n] <- which.max(result$v[[n]])
  chemin_labels[n] <- etat[chemin_indices[n]]
  
  # Remontée
  for (i in (n-1):1) {
    chemin_indices[i] <- result$to[[i+1]][chemin_indices[i+1]]
    chemin_labels[i] <- etat[chemin_indices[i]]
  }
  result$chemins <- chemin_labels
  
  return(result)
}