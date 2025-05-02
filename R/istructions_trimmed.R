# Inizializzazione:
#   
#   Parti con k centroidi (scelti a caso o con k-means++).
# 
# Assegnamento e trimming (1a iterazione):
#   
#   Assegni tutte le osservazioni ai centroidi (anche quelle che poi scarterai).
# 
# Calcoli le distanze di tutte le osservazioni dai rispettivi centroidi.
# 
# Ordini le distanze e scarti il trim_lev più distante (es. il 10% più lontano), indipendentemente dal cluster.
# 
# Le osservazioni trimmate non vengono usate per aggiornare i centroidi.
# 
# Aggiornamento:
#   
#   Ricalcoli i centroidi solo con le osservazioni non trimmate.
# 
# Iterazioni successive:
#   
#   Riutilizzi tutte le osservazioni (incluse quelle che erano outlier prima).
# 
# Le riassegni a un cluster.
# 
# Ricalcoli di nuovo le distanze.
# 
# Scarti un nuovo set di outlier: potrebbe includere altri punti, o persino riaccettare quelli scartati nella prima iterazione se ora sono più vicini.