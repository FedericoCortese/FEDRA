load("FEDRA_varsxclust_all.RData") # 27 variabili, ovvero più o meno tutti gli score principali dei task fatti
load("FEDRA_varsxclust_k7.RData") # 7 variabili selezionate per rappresentare teoricamente costrutti ben distinti ed essenziali
load("FEDRA_varsxclust_k12.RData") # 12 variabili, che rappresentano costrutti non così distinti ma comunque diversi e potenzialmente utili

# In ogni file c'è una variabile IDK che rappresenta il soggetto e una variabile GRP che rappresenta il gruppo. 
# I soggetti con GRP = EXC vanno esclusi.


