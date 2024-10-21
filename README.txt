This file serves as supportive documentation to Shelby McCahon's neonicotinoid body condition analysis. 

Date created: 2024-10-21

-------------------------

Objective: Determine if there are any relationships between shorebird body condition and neonicotinoid concentrations/detections. We are assessing detections (yes or no) as an explanatory variable due to the rapid metabolism of neonicotinoids. 

Body condition response variables: body mass, size corrected body mass, fat, pectoral muscle size, pectoral muscle score, scaled mass index, fattening index (metabolite data)

Explanatory variables of interest: neonicotinoid detection (yes or no), neonicotinoid concentration (ug/L)

Controlling variables: sex, year, Julian date (if non-migratory species), date into migration season (if migratory species), season (spring/fall)

Species/group and sample size: 
1. Lesser Yellowlegs (n = 54)
2. Least Sandpiper (n = 29)
3. Calidris sp. (Least Sandpipers, Semipalmated Sandpipers, n = 49)
4. Killdeer (n = )
5. Pectoral Sandpiper (n = )
6. Semipalmated Sandpipier (n = )
7. Tringa sp. (Lesser Yellowlegs and Willet, n = )
8. Wilson's Phalarope (n = )
9. Multiple Species (those with >10 birds/species, n = )
10. Multiple Spcies (those with >15 birds/species, n = )
11. Rare Species (those with <10 birds/species, n = )

-------------------------

Analysis Overview:

I used a two-stage AICc model selection approach to determine if neonicotinoid concentrations/detections further explain variation in body condition beyond the controlling variables. The first stage AICc model selection includes all non-correlated variables aside from neonicotinoids along with all two-way interactions among them. The top model(s) (within 2 delta AICc) from this first stage became the "informed null model." The second stage AICc model selection included the informed null model, all possible combinations and interactions of neonicotinoid concentrations/detections with the informed null model variables, and all possible combinations and interactions of neonicotinoid concentrations/detections with the global model variables. 

-------------------------

