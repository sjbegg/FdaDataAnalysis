# FdaDataAnalysis


Exploratory analysis of adverse event reports from the FDA database (https://open.fda.gov/drug/event/)

The matlab code in this repository carries out a preliminary investigation of adverse events reported by country, and whether there is correlation either between individual events, or between countries. 

## Code


Matlab code is provided in the file 'AdverseEventsInDifferentCountries.m'.


## Output 

The folder 'Visualisations' contains the following output:

* CountriesThatAreTopTenReporters.png: A plot showing the frequency with which a given country appears as one of the 'top 10' reporting countries for that event. Only 25 countries create this set. 

* ReactionsThatAreTopTenReported.png: A plot showing the frequency with which a given reaction is one of the 'top 10' reported across all countries. 72 reactions create this set. 

* HeatMapOfCorrelationBetweenReactions.png: A heatmap of the Pearson correlation coefficients for each pair of adverse events reported across all countries.

* HeatMapOfReactionsWithPearsonPoint8.png: A heatmap of 0 / 1 values indicating whether a given pair of reactions has a Pearson correlation coefficient greater than 0.8, indicating very strong correlation across all reporting countries. 

* HeatMapOfCorrelationBetweenCountries.png: A heatmap of the Pearson correlation coefficient for each pair of countries across all reported adverse events. 

* HeatMapOfCountriesWithPearsonPoint8.png: A heatmap of 0 / 1 values indicating whether a given pair of countries has a Pearson correlation coefficient greater than 0.8, indicating very strong correlation across all events. 

* HierarchicalClustering.png: A plot showing hierarchical clustering (Euclidean distance metric) + dendrogram. 
