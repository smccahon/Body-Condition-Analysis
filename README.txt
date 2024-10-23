This file serves as supportive documentation to Shelby McCahon's neonicotinoid body condition analysis. 

Date created: 2024-10-21

-------------------------

Objective: Determine if there are any relationships between shorebird body condition and neonicotinoid concentrations/detections. We are assessing detections (yes or no) as an explanatory variable due to the rapid metabolism of neonicotinoids. 

Body condition response variables: body mass, size corrected body mass, fat, pectoral muscle size, pectoral muscle score, scaled mass index, fattening index (metabolite data)

Explanatory variables of interest: neonicotinoid detection (yes or no), neonicotinoid concentration (ug/L)

Controlling variables: sex, year, Julian date (if non-migratory species), date into migration season (if migratory species), season (spring/fall)

Species/group and sample size: 
1. Lesser Yellowlegs (n = 54; 23 with detections)
2. Least Sandpiper (n = 29; 4 with detections)
3. Semipalmated Sandpiper (n = 20; 7 with detections)
4. Small Sandpipers (Least Sandpipers, Semipalmated Sandpipers, n = 49; 11 with detections)
5. Killdeer (n = 14; 4 with detections)
6. Pectoral Sandpiper (n = 16; 2 with detections <- did not do analysis because of this)
7. Tringa sp. (Lesser Yellowlegs and Willet, n = )
8. Wilson's Phalarope (n = )
9. Multiple Species (those with >10 birds/species, n = )
10. Multiple Species (those with >15 birds/species, n = )
11. Rare Species (those with <10 birds/species, n = )
12. Calidris sp. (Least Sandpiper, Pectoral Sandpiper, Semipalmated Sandpiper, n = 63; 13 with detections)

-------------------------

Analysis Overview:

I used a two-stage AICc model selection approach to determine if neonicotinoid concentrations/detections further explain variation in body condition beyond the controlling variables. The first stage AICc model selection includes all non-correlated variables aside from neonicotinoids along with all two-way interactions among them. The top model(s) (within 2 delta AICc) from this first stage became the "informed null model." The second stage AICc model selection included the informed null model, all possible combinations and interactions of neonicotinoid concentrations/detections with the informed null model variables, and all possible combinations and interactions of neonicotinoid concentrations/detections with the global model variables. 

-------------------------

Conclusions (Mass ~ Neonicotinoid):

1. Lesser Yellowlegs (Mass ~ Neonicotinoid): Neonicotinoids (raw concentrations, log transformations, and detections) do not significantly impact body mass. Year is the only variable that had a significant impact on Lesser Yellowlegs body mass. Lesser Yellowlegs in 2023 had significantly lower body mass compared to birds in 2021. 
2. Least Sandpipers (Mass ~ Neonicotinoid): Birds with neonicotinoid detections (and higher concentrations) had significantly higher body mass (Mann-Whitney, W = 90, p = 0.01). However, only 4 birds had detections. This result is most likely due to small sample size.
3. Semipalmated Sandpipers (Mass ~ Neonicotinoid): There were no variables (including neonicotinoid concentrations/detections) that had any significant impact on body mass. 
4. Small pipers (Mass ~ Neonicotinoid): Conflicting results across models. Models with neonicotinoid detections suggest that birds with detections had significantly lower body mass. However, when graphing the relationship independently of the model, boxplot shows that birds WITH detections have significantly higher body mass. For both models that used neonicotinoids as a continuous variable (log transformation and no log transformation), model suggests birds with higher neonic concentrations have higher body mass (positive relationship). 
5. Killdeer (Mass ~ Neonicotinoid): Neonicotinoid detection does not have any influence on Killdeer body mass. I did not analyze concentrations as a continuous variable due to extreme violations in homoscedasticity. 
6. Calidris sp. (Mass ~ Neonicotinoid): 

