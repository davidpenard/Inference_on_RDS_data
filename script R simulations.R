rm(list = ls())
# Packages
library(igraph)
library(tidyverse)
library(dbplyr)

### Fonctions

SimulGraph_pk = function(p0k,N0_hat,p1k,N1_hat)
{
  d_sim0 = sample(x=1:length(p0k),size=N0_hat,prob=p0k,replace=T) #distribution de degrés (z=0)
  d_sim1 = sample(x=1:length(p1k),size=N1_hat,prob=p1k,replace=T) #distribution de degrés (z=1)
  d_sim = c(d_sim0,d_sim1)
  z_sim = c(rep(0,N0_hat),rep(1,N1_hat))
  # si la somme des demi-arêtes n'est pas paire
  if(sum(d_sim) %% 2 != 0) d_sim[1]=d_sim[1]+1
  nbaretes = sum(d_sim)/2
  shuffle<-sample(x=rep(1:length(d_sim),d_sim),size=nbaretes*2, replace=FALSE)
  aretes<-matrix(data=shuffle,nrow=nbaretes,ncol=2 )
  graph0 <- graph.data.frame(aretes,directed=FALSE)
  return(list(graph0=graph0,aretes=aretes,z_sim=z_sim ))
}


RDS = function(NbCoupons, # nombre total de coupons qu'on peut distribuer
               CouponsADistribuer,# nombre de coupons qu'on donne à chaque personne)
               NbVaguesMax = 1e10,
               exploRDS = 1000, 
               Graph){
  ### Initialisation
  graph0 = Graph$graph0 
  aretes = Graph$aretes
  sommets=as.numeric(V(graph0)$name)
  N = length(sommets)
  components <- decompose.graph(graph0,mode="weak")
  length(components)
  comp.sizes <- unlist(lapply(components,vcount)) # vcount() 
  comp.sizes.d <- table(comp.sizes)
  
  couponsrestants=NbCoupons
  aretes_rds = matrix(ncol=2,byrow = T)
  
  # la proba de tirage détermine comment on distribue les coupons : 
  # quelle est la proba de tirer les individus 1 (jamais mentionnés) ou 2 (mentionnés avant)
  # (les codes 3 et 4 sont toujours à 0 car on ne leur redonne pas de coupons)
  # idée : si probatirage=c(0,1,0,0), on ne distribue des coupons qu'aux individus déjà mentionnés
  # on ne redonne pas de coupon à des gens déjà tirés ou ayant déjà un coupon
  probatirage=c(1,1,0,0)# on ne tire que les individus qui n'ont pas encore été tirés
  etats=rep(1,N) 
  # codage : 1=non mentionné, 2=mentionné mais pas de coupon, 3=dispose d'un coupon, 4=interrogé
  
  # Echantillonnage #
  v=sample.int(N,1) #v est le premier individu interrogé, choisi au hasard
  etats[v]=4
  if(length(adjacent_vertices(graph0, sommets == sommets[v], mode = "all"))>0){
    # adjacent_vertices() récupère les points connectés avec "v" (igraph) 
    
    listevoisins=which(sommets%in%sommets[as.numeric(names(table(adjacent_vertices(graph0, sommets == sommets[v], mode = "all"))))])
    etats[listevoisins[etats[listevoisins]==1]]=2 # les voisins sont mentionnés
    
    # On note les arètes entre v et ses voisins 
    for(voisin in listevoisins) aretes_rds = rbind(aretes_rds,c(sommets[v],sommets[voisin]))
    coup=as.numeric(sample(as.character(listevoisins),min(CouponsADistribuer,length(listevoisins),couponsrestants),replace=FALSE))
    #on tire "CouponsADistribuer" voisins parmi tout les voisins de l'individu échantillonné 
    # le coup du as.numeric as.character est pour éviter que R ne tire parmi les entiers entre 1 et listevoisins lorsque cette liste ne contient qu'un seul entier
    etats[coup]=3
    couponsrestants=NbCoupons-min(CouponsADistribuer,length(listevoisins),couponsrestants)
  }
  t=1
  
  ### Propagation ###
  while((couponsrestants>0 &sum(etats>=2)<exploRDS & sum(etats<3)>0))
  {
    #s'il n'y a plus de coupons en circulation, on interroge une nouvelle personne prise au hasard
    if(sum(etats==3)==0)
    {
      # print("aucun coupon en circulation !")
      v=sample.int(N,1,replace=FALSE,prob=(etats<3))
    }
    # s'il reste des coupons, l'un des individus avec un coupon est interrogé
    if(sum(etats==3)>0)
    {
      v=sample.int(N,1,replace=FALSE,prob=(etats==3))
    }
    etats[v]=4
    #ses voisins sont mentionnés
    voisins1=0
    voisins2=0
    
    if(length(adjacent_vertices(graph0, sommets == sommets[v], mode = "all"))>0)
    {
      listevoisins=which(sommets%in%sommets[as.numeric(names(table(adjacent_vertices(graph0, sommets == sommets[v], mode = "all"))))])
      # On note les arètes entre v et ses voisins 
      for(voisin in listevoisins) aretes_rds = rbind(aretes_rds,c(sommets[v],sommets[voisin]))
      voisins1=sum(etats[listevoisins]==1)
      voisins2=sum(etats[listevoisins]==2)
    }  
    if(voisins1+voisins2>0)
    {
      etats[listevoisins[etats[listevoisins]==1]]=2
      if(couponsrestants>0)
      {
        coup=as.numeric(sample(as.character(listevoisins),size=min(CouponsADistribuer,couponsrestants,voisins1*(probatirage[1]>0)+voisins2*(probatirage[2]>0)),replace=FALSE,prob=probatirage[etats[listevoisins]]))
        etats[coup]=3
        couponsrestants=couponsrestants-length(coup)
      } 
    }
    #print(t);
    t=t+1
  }
  #S_i, d_i, z_i #
  S = as.numeric(etats[order(sommets)]>=3)
  # print(table(etats))
  # print(etats[order(sommets)])
  #Matrice d'adjacence
  rdsGraph = graph.data.frame(aretes_rds[-1,],directed=FALSE)
  # print(aretes_rds)
  # print(rdsGraph)
  # vertex.attributes(rdsGraph,index = S)
  # AdjMat = matrix(data = as.numeric(unlist(map(1:N, function(i) 
  #   1:N%in%as.numeric(voisins[[i]][-1])))), nrow = N)   
  # diag(AdjMat) = rep(0,N)
  # rownames(AdjMat) = colnames(AdjMat) = 1:N
  return(list(rdsGraph=rdsGraph, S=S))
}

F_c = function(C,x) unlist(map(x,.f = function(u) sum(C[C<=u])/sum(C)))
# "Fonction de répartition de C"

count_vert_z = function(graphei,z,Z_data){
  return(unlist(lapply(decompose.graph(graphei,mode="weak"),
                       FUN = function(graphComp) sum(Z_data[as.numeric(names(V(graphComp)))]==z)
                       )))}

d_comp = function(graphe1,graphe2,z){
  # Les différentes composantes ordonnées par taille décroissantes
  # On travail avec des proportions
  comp1 = count_vert_z(graphe1,z,z_True) %>% 
    .[order(.,decreasing = F)]
  # print(comp1)
  comp2 = count_vert_z(graphe2,z,z_sim) %>% 
    .[order(.,decreasing = F)]
  # print(comp2)
  comp12 = c(comp1,comp2) %>% 
    .[order(.,decreasing = F)]
  # On renvoie la différence entre les 2 répartions
  return(sum(unlist(map(1:(length(comp12)-1),function(i) 
    abs(F_c(comp1,comp12[i])-F_c(comp2,comp12[i]))*(comp12[i+1]-comp12[i]))
  )))
}

RecSim = function(d_crit, K,iteration){
  if(d_crit <= 0 ) return(1)
  else return(rbinom(1,1,prob = exp(-K*d_crit*log(iteration))) )
}

# Graphe Simulé

N = N_True= 1000 # taille de la population simulée
PropCov1= 1/5

# deg_dist = read.csv(file = "work/pk_table.csv")$pk
deg_dist = c(0.754900676, 0.140258284, 0.056241795, 0.036657892, 0.008764305, 0.002101549, 0, 0.001075499)
True_pop = SimulGraph_pk(deg_dist,N*(1-PropCov1),deg_dist,N*PropCov1)
N0_True = round(N*(1-PropCov1))
N1_True = round(N*PropCov1)
d_True = unlist(map(V(True_pop$graph0), function(vertex)
  length(table(adjacent_vertices(True_pop$graph0,V(True_pop$graph0) == vertex, mode = "all")))))
### Covariables (ex : genre)
z_True = True_pop$z_sim
p_geom = 1/sum((1:max(d_True))*deg_dist[1:max(d_True)])
pk_True = unlist(map(1:N, function(k) ifelse(sum(d_True==k)==0,0, sum(d_True==k)/length(d_True))))
p0k_True = unlist(map(1:max(d_True), function(k) sum(d_True==k&z_True==0)/length(d_True[z_True==0])))
p1k_True = unlist(map(1:max(d_True), function(k) sum(d_True==k&z_True==1)/length(d_True[z_True==1])))

#paramètres du RDS
NbCoupons= 1000 # nombre total de coupons qu'on peut distribuer
# S'arreter après 150 indiv interrogés
CouponsADistribuer= 3 # nombre de coupons qu'on donne à chaque personne
NbVaguesMax = 50
exploRDS = 245
NbModel = 25# M1
# importance du critère d_TV par rapport au critère d_C
alpha = 1/2





### Utilisation des données de l'étude
data_reseau = read.csv(file = "work/donnees_aretes.csv",header = F,sep = " ")
data_individus = read.csv(file = "work/donnees_individus.csv",header = T,sep = ",")
aretes_data = data_reseau[,-3]
aretes_data = matrix(aretes_data[!is.na(aretes_data)],ncol=2)
deg = data.frame(row.names = 1:max(aretes_data), unlist(map(1:max(aretes_data), function(i) sum(aretes_data == i))))
RDS_Iris = graph.data.frame(aretes_data[,],directed = F)
### Variable échantillonnage (interrogé/coupon ou juste mentionné)
S_RDS_Iris = rep(0,max(data_individus$N..ego))
for( u in 1:max(data_individus$N..ego)){
  if(u%in%data_individus$N..ego)
    if(data_individus$questionnaire[which(data_individus$N..ego==u)] %in% c("1","2"))
      S_RDS_Iris[u] = 1
}
### Variable genre (H=0, F=1)
z_RDS = rep(NA,max(data_individus$N..ego))
for( u in 1:max(data_individus$N..ego)){
  if(u%in%data_individus$N..ego)
    z_RDS[u] = ifelse(data_individus$sexe[which(data_individus$N..ego==u)]=="F",1,0)
}
### Séparation des différentes composantes en 3 graphes (Paris, Aulnay + Slammeurs)
id_indiv_chemsex = data_individus$N..ego[which(data_individus$Chemsex==1)]
id_indiv_Aulnay =  data_individus$N..ego[which(data_individus$CAARUD=="Aulnay")]
id_indiv_Paris = data_individus$N..ego[which(!(data_individus$N..ego%in%c(id_indiv_chemsex,id_indiv_Aulnay))&data_individus$questionnaire==1)]

comp_graphe = decompose.graph(RDS_Iris,mode="weak")

graphe_Paris = comp_graphe[unlist(map(comp_graphe,function(comp) any(as.numeric(names(V(comp)))%in%id_indiv_Paris)))] %>%
  graph.union
id_indiv_Paris2=as.numeric(names(V(graphe_Paris)))
id_indiv_Paris2[order(id_indiv_Paris2)]
z_True = z_RDS


### Simulations

#Mesure des grandeurs à la fin de chaque simulation
mesure_EMV = c() # estimation de p_geom
mesure_d_TV = c() # distance en variation totale sur la distribution de degrés
mesure_N = c() # Estimation de N
repetitions = 100
# list_K_d_C = rep(c(1,0.1,0.01),each=repetitions)
list_K_d_C = rep(1,repetitions)
simulation=1
list_seed = array(sample.int(10000000,repetitions),dim=c(6*repetitions))
for(simulation in 1:length(list_K_d_C)){
  set.seed(list_seed[simulation])
  print(paste("start new simulation : K =",list_K_d_C[simulation]," simulation : ",simulation%%repetitions))
  list_N0_hat_t=list_N1_hat_t=liste_d_C0_t=liste_d_C1_t=liste_d_TV_t=NULL  # Première exploration par RDS
  ## Sur données simulées
  # ResRDS = RDS(NbCoupons = NbCoupons,CouponsADistribuer = CouponsADistribuer,NbVaguesMax = NbVaguesMax,exploRDS=exploRDS,Graph = True_pop)
  # rdsGraph = ResRDS$rdsGraph
  # S = ResRDS$S # les individus mentionnés
  
  # Sur données réelles
  rdsGraph = graphe_Paris
  exploRDS = 245 #length(V(rdsGraph)) 
  S = as.numeric(S_RDS_Iris==1&(1:399)%in%names(V(rdsGraph)))
sum(S)
  
  Y_data = get.adjacency(rdsGraph,)
  Y_data = matrix(ifelse(Y_data>0,1,0),nrow = nrow(Y_data),ncol = ncol(Y_data),
                  dimnames = list(colnames(Y_data),colnames(Y_data))) #ordre graphe
  diag(Y_data)=rep(0,nrow(Y_data))
  nodenames = as.numeric(colnames(Y_data))
  # d_data = unlist(map(1:nrow(Y_data), function(i) sum(Y_data[nodenames==i,]))) %>% 
  #   .[which(S[nodenames]==1)]# ordre graphe
  d_data = unlist(map(1:length(S), function(i){
    if(i%in%nodenames) return(sum(Y_data[nodenames==i,]))
    else return(0)}))
  # d_data = d_data[S==1]
  d_data_z0 = ifelse((1-z_True)*S==1,d_data,NA)
  d_data_z1 = ifelse(z_True*S==1,d_data,NA)
  
  # z_True = z_True[S==1]
  # S= S[S==1]
  # distribution de degrés observée 
  pk_obs = unlist(map(1:max(d_data), function(k) sum(d_data==k)/sum(S)))
  # probabilité d'échantillonnage initiale proportionnelle à d_i
  Pi0_k = (1:N)/sum(d_data[z_True[S==1]==0]) *sum(S==1&z_True==0) # ordre des k (1:N) #sum(S) = nombre de personnes interrogées
  Pi1_k = (1:N)/sum(d_data[z_True[S==1]==1]) *sum(S==1&z_True==1)
  # N_hat = round(sum(1/Pi_k[d_data[which(Pi_k[d_data]!=0)]])) # ignorer les pi_k = 0
  # N0_hat = round(length(V(rdsGraph))*mean(1-z_True[S==1]))*2
  # N1_hat = round(length(V(rdsGraph))*mean(z_True[S==1]))*2
  N0_hat = 800
  N1_hat = 200
  print(paste("N0 : ",N0_hat,"N1 : ",N1_hat))
  
  pk =  rep(0,max(d_data))
  for(k in 1:length(pk)){
    pk[k] =  ifelse(Pi0_k[k]!=0, sum((d_data==k)/Pi0_k[k])/sum(S),0)}
  pk = pk/sum(pk)
  d_TV_courant = 100
  d_C0_courant = d_C1_courant = 100
  
  # #### Mesures des erreurs par rapport au vrai graphe
  mesure_d_EMV = c(abs(1/mean((1:max(d_data))*pk)-p_geom)) # estimateur du max de vraissemblance
  mesure_d_N = c(abs(N_True-(N0_hat+N1_hat)))   # différence N_hat - N
  mesure_d_TV = c(sum(abs(pk - pk_True[1:max(d_data)])))
  
  
  ### Boucle (3)
  
  for(iteration in 1:50){
    ## Simulation de M graphes suivant la distribution de degrés pk
    list_d_TV = c()
    list_d_C0 = c()
    list_d_C1 = c()
    table_Pi0_k_model = table_Pi1_k_model = matrix(ncol = max(d_data))
    for(gph in 1:NbModel){#NbModel
      # Création du graphe (fonction à créer)
      resSimul = SimulGraph_pk(pk,N0_hat,pk,N1_hat)
      #exploration par RDS
      simRDS = RDS(NbCoupons=NbCoupons, CouponsADistribuer=CouponsADistribuer, NbVaguesMax=NbVaguesMax,exploRDS=exploRDS, Graph=resSimul)
      # extraire les vecteurs d_sim et S_sim et estimation de Pi_k
      d_sim = unlist(map(1:(N0_hat+N1_hat), function(k) sum(resSimul$aretes==k)))
      S_sim = simRDS$S
      # sum(S_sim)
      z_sim = resSimul$z_sim
      Pi0_k_model = unlist(map(1:max(d_data), function(k) ifelse(sum(d_sim==k&z_sim==0)==0,0,sum(d_sim==k&S_sim==1&z_sim==0)/sum(d_sim==k&z_sim==0))))
      Pi1_k_model = unlist(map(1:max(d_data), function(k) ifelse(sum(d_sim==k&z_sim==1)==0,0,sum(d_sim==k&S_sim==1&z_sim==1)/sum(d_sim==k&z_sim==1))))
      
      # calcul de pk (la distribution de degrés) du graphe simulé
      pk_model = unlist(map(1:max(d_data), function(k) sum(d_sim==k&S_sim==1)/sum(S_sim)))
      # p1k_model = unlist(map(1:max(d_data), function(k) sum(d_sim==k&S_sim==1&z_sim==1)/sum(S_sim[z_sim==1])))
      
      #calcul du d_TV pour le modèle
      if(length(pk_obs)<=length(pk_model)){d_TV = sum(abs(pk_model -c(pk_obs,rep(0,length(pk_model)-length(pk_obs)))))
      }else d_TV = sum(abs(pk_obs -c(pk_model,rep(0,length(pk_obs)-length(pk_model)))))
      list_d_TV = c(list_d_TV,d_TV)
      
      # d_comp distance entre les distribution de taille de composantes
      d_C0 = d_comp(rdsGraph,simRDS$rdsGraph,0)
      list_d_C0 = c(list_d_C0,d_C0)
      d_C1 = d_comp(rdsGraph,simRDS$rdsGraph,1)
      list_d_C1 = c(list_d_C1,d_C1)
      
      table_Pi0_k_model = rbind(table_Pi0_k_model,Pi0_k_model)
      table_Pi1_k_model = rbind(table_Pi1_k_model,Pi1_k_model) 
    }
    list_crit = alpha*list_d_TV + (1-alpha)*(list_d_C0+list_d_C1) # critère de sélection des graphes 
    table_Pi0_k_model = table_Pi0_k_model[-1,] # on retire la première ligne de NA 
    table_Pi1_k_model = table_Pi1_k_model[-1,]
    
    
    # on prend les meilleurs 10%
    in_seuil_crit = list_crit <= list_crit[order(list_crit)][round(0.2*NbModel)] 
    # Moyenne des meilleurs Pi_k estimés
    Pi0_k = as.numeric(colSums(table_Pi0_k_model[in_seuil_crit,])/sum(in_seuil_crit))  
    Pi1_k = as.numeric(colSums(table_Pi1_k_model[in_seuil_crit,])/sum(in_seuil_crit))
    
    print(Pi0_k)
    print(Pi1_k)
    #estimation des proportions des degrés (provisoire)
    # print(pk)
    d_TV_nouveau = mean(list_d_TV[in_seuil_crit])
    pk_nouveau = rep(0,max(d_data))
    for(k in 1:length(pk)){
      pk_nouveau[k] =  ifelse(Pi0_k[k]!=0, sum((d_data==k)/Pi0_k[k])/N0_hat,0)
    }
    #recuit simulé
    # z=0
    res_RecSim = RecSim(d_crit=d_TV_nouveau-d_TV_courant,K = 1 ,iteration = iteration )
    d_TV_courant = res_RecSim*d_TV_nouveau + (1-res_RecSim)*d_TV_courant
    pk = res_RecSim*pk_nouveau + (1-res_RecSim)*pk
    pk= pk/sum(pk)
    
    # critère de tailles des composantes 
    # z=0
    d_C0_nouveau = mean(list_d_C0[in_seuil_crit])
    N0_hat_nouveau = round(sum(1/Pi0_k[d_data_z0[which(Pi0_k[d_data_z0]!=0)]])) # ignorer les pi_k = 0
    #recuit simulé
    res_RecSim = RecSim(d_crit=d_C0_nouveau-d_C0_courant,K = list_K_d_C[simulation],iteration = iteration )
    print(paste("delta_d_C0 :",d_C0_nouveau-d_C0_courant))
    N0_hat = res_RecSim*N0_hat_nouveau + (1-res_RecSim)*N0_hat
    d_C0_courant = res_RecSim*d_C0_nouveau + (1-res_RecSim)*d_C0_courant
    # z=1
    d_C1_nouveau = mean(list_d_C1[in_seuil_crit])
    N1_hat_nouveau = round(sum(1/Pi1_k[d_data_z1[which(Pi1_k[d_data_z1]!=0)]])) # ignorer les pi_k = 0
    #recuit simulé
    res_RecSim = RecSim(d_crit=d_C1_nouveau-d_C1_courant,K = list_K_d_C[simulation],iteration = iteration )
    print(paste("delta_d_C1 :",d_C1_nouveau-d_C1_courant))
    N1_hat = res_RecSim*N1_hat_nouveau + (1-res_RecSim)*N1_hat
    d_C1_courant = res_RecSim*d_C1_nouveau + (1-res_RecSim)*d_C1_courant
    print(paste("N0:",N0_hat, "  N1:",N1_hat,"  N:",N0_hat+N1_hat))
    if(N0_hat+N1_hat<exploRDS){
      print("Erreur : N_hat est trop petit ")
      break}
    
    if(N0_hat+N1_hat>5*N){
      print("Erreur : N_hat est trop grand ")
      prop1=N0_hat/(N0_hat+N1_hat)
      N0_hat = round((N0_hat+N1_hat)*prop1)
      N1_hat = round((N0_hat+N1_hat)*(1-prop1))
      break}
    
    list_N0_hat_t = c(list_N0_hat_t, N0_hat)
    list_N1_hat_t = c(list_N1_hat_t, N1_hat)
    liste_d_C0_t = c(liste_d_C0_t,d_C0_courant)
    liste_d_C1_t = c(liste_d_C1_t,d_C1_courant)
    liste_d_TV_t = c(liste_d_TV_t, d_TV_courant)
simulation_table = data.frame(list_N0_hat_t=list_N0_hat_t,list_N1_hat_t=list_N1_hat_t,
                              liste_d_C0_t=liste_d_C0_t,liste_d_C1_t=liste_d_C1_t,
                              liste_d_TV_t=liste_d_TV_t)
write.csv(simulation_table, 
          file = paste0("resultats sur données/Resultat_données_sim",simulation%%repetitions,".csv"))
  }
  mesure_EMV = c(mesure_EMV,1/sum((1:max(d_data))*pk)) # estimation de p_geom par EMV
  mesure_d_TV = c(mesure_d_TV,d_TV_courant) # distance en variation totale sur la distribution de degrés
  mesure_N = c(mesure_N,N0_hat+N1_hat) # Estimation de N
}


resultats_estim = data.frame(EMV = mesure_EMV,d_TV = mesure_d_TV,  N = mesure_N, K = list_K_d_C)
  write_csv(x = resultats_estim,file = "Resultats calibrage K//resultats_estim_serveur_calibrage.csv",append = FALSE)


