% Exploring FDA adverse drug effect data
% SJ Dunn

%//////////////////////////////////////////////////////////////
% Are different adverse events reported in different countries?
%//////////////////////////////////////////////////////////////

% Preliminary analysis for 100 reactions reveals that 

% 25 countries that span the top ten reporting countries for each reaction
% 14 countries span the top five
% 11 countries span the top three
% The USA is the top reporting country for each reaction

% 72 reactions span the top 10 reported reactions for each country.
% Together with the above, this could suggest that some countries 'do
% better' at filing the reports, which could skew the results

% Correlation analysis (using Pearson's coefficient):
% A small number of reactions are well-correlated

% There are many instances of pairs of countries reporting the same events,
% but cluster analysis does not reveal significant groupings of countries
% that might suggest that events are linked to geographical regions or
% demographics

% If we reduce the number of reactions that we look for the results don't
% show a change in trend

clear all 
close all
clc

% Set a limit the number of reactions we search for
limit = 100;

% API calls
api = 'https://api.fda.gov/drug/event.json?api_key=SHmwbhZldi8xHYddZTHfVI6iVd5sVhKEbgXam9ij&';

% Set query to count the top patient reactions
query = ['count=patient.reaction.reactionmeddrapt.exact&limit=', num2str(limit)];

% Set up the url to make the query
url = [api query];

% Pull out the results of the search as a struct
numInstancesReactions = webread(url);


% Loop over the different reactions to count the number of countries that
% report them

% Do first case separately to get order of countries to compare against
query = ['search=patient.reaction.reactionmeddrapt:"', numInstancesReactions.results(1).term,...
    '"&count=occurcountry'];
url = [api query];
numByCountries = webread(url);

% Pull out the list of country references
orderOfCountries = {numByCountries.results(:).term};

% Sort countries into alphabetical order 
[sortedCountries, I] = sort(orderOfCountries);

counts = [numByCountries.results.count];
sortedCounts = counts(I);   % Sort the counts according to alphabetical order

% Initialise array with first results to store the counts per country
countPerCountry = sortedCounts';

% Now do the remaining reactions, and sort the results into the order
% determined above
for reaction = 2:size(numInstancesReactions.results,1)
        
    query = ['search=patient.reaction.reactionmeddrapt:"', numInstancesReactions.results(reaction).term,...
        '"&count=occurcountry'];
    url = [api query];
    numByCountries = webread(url);
    
    orderOfCountries = {numByCountries.results(:).term};
    counts = [numByCountries.results.count];
    
    for i = 1:length(orderOfCountries)        
        indexC = strfind(sortedCountries, orderOfCountries{i});
        index = find(not(cellfun('isempty', indexC)));
        sortedCounts(index) = counts(i);
    end
    
    % Rows are countries, cols are reactions
    % Matrix order: sortedCountries x listOfEvents
    countPerCountry = [countPerCountry sortedCounts'];
           
end 

%                       //////////////////////////
%                           Now analyse the data
%                       //////////////////////////

%//////////////////////////////////////////////////////////////////////
% Frequency of a given country being a top 10 reporter of the reaction
%//////////////////////////////////////////////////////////////////////

% Sort the matrix by column so that each event is ordered by the most
% reports 
[colsOrderedByCountry, indices] = sort(countPerCountry, 'descend');

topN = 10;  % Set to examine the top n reporting countries

% Top ten countries per reaction (use the indices to access the country
% names)
topNCountries = indices(1:topN,:);
topNCountries = topNCountries(:);

% Transform country codes to actual names
[~,countryCodesToNames] = xlsread('CountryCodes.xlsx');
countryCodesToNames(1,:) = []; % Remove headers

xLabels = sortedCountries(unique(topNCountries)); % Get the codes
newLabels = cell(size(xLabels));    
for i = 1:length(xLabels)
      index = find(strcmp(countryCodesToNames(:,2),xLabels{i}));
      newLabels{i} = countryCodesToNames{index,1};      % Pull out full name
end

figure
[y, ~] = hist(topNCountries, 100);  % Pull out frequency for each country (bins=100 due to indices being from 1-100)
y(y==0) = [];
[y, index] = sort(y, 'descend');
x = 1:length(y);
bar(x,y)
ylabel('Frequency')
axis([0 max(x)+1 0 max(y)+10])
set(gca, 'XTick', x);
set(gca, 'XTickLabel', newLabels(index));
set(gca, 'XTickLabelRotation', 90);
title(['Frequency of Being a Top ', num2str(topN), ' Reporter (', num2str(limit), ' Reactions)'])

%//////////////////////////////////////////////////////////////////////
% Frequency of a given reaction being in the top 10 reported by country
%//////////////////////////////////////////////////////////////////////

% Sort the matrix by rows so that each country is ordered by the most
% reported events
[rowsOrderedByEvents, eventIndices] = sort(countPerCountry, 2, 'descend');

% Top n events per country 
topNEventsPerCountry = eventIndices(:,1:topN);
topNEventsPerCountry = topNEventsPerCountry(:);

% Transform event to actual description
% Pull out the list of event references
listOfEvents = {numInstancesReactions.results(:).term};
xLabels = listOfEvents(unique(topNEventsPerCountry)); 

figure
[y, ~] = hist(topNEventsPerCountry, 100);  % Pull out frequency for each country (bins=100 due to indices being from 1-100)
y(y==0) = [];
[y, index] = sort(y, 'descend');
x = 1:length(y);
bar(x,y)
ylabel('Frequency')
set(gca, 'XTick', x);
axis([0 max(x)+1 0 max(y)+10])
set(gca, 'XTickLabel', xLabels(index));
set(gca, 'XTickLabelRotation', 90);
title(['Frequency of Being a Top ', num2str(topN), ' Reported Reaction'])


%//////////////////////////////////////////////////////////////////////
% What is the extent of (linear) correlation between any two reactions 
% across the reporting countries
%//////////////////////////////////////////////////////////////////////

% Normalise the data so that we account for different extents of reporting
% between countries
normalisedCountsPerCountry = bsxfun(@rdivide, countPerCountry, max(countPerCountry,[],2));

% Returns correlation between columns
[corrBetweenReactions, pValuesReactions] = corr(normalisedCountsPerCountry, 'type', 'Pearson');

% Set to zero any correlations with p value greater than 0.01
pIndices = find(pValuesReactions>0.01);
corrBetweenReactions(pIndices) = 0;
% Make sure we haven't screwed up the main diagonal
corrBetweenReactions(1:(size(corrBetweenReactions,1)+1):end) = 1;

% Plot heatmap of correlation
figure
imagesc(corrBetweenReactions);
colorbar;
title(['Correlation Between Top ', num2str(limit), ' Reactions'])
xlabel('Reactions')
ylabel('Reactions')

% Plot only those with high correlation (plotting bools)
correlationThreshold = 0.8;
figure
imagesc(corrBetweenReactions>correlationThreshold);
colorbar;
title(['Those Reactions with \rho > ', num2str(correlationThreshold)])
xlabel('Reactions')
ylabel('Reactions')

% Pull out the pairs of highly correlated reactions (ignoring self-pairs)
[whichReactionsRow, whichReactionsCol] = find(triu(corrBetweenReactions)>correlationThreshold & triu(corrBetweenReactions)<1);
correlatedReactions = {};
% Store list of highly correlated reactions
if ~isempty(whichReactionsRow)
    for i = 1:length(whichReactionsRow)
        correlatedReactions{i,1} = [cell2mat(listOfEvents(whichReactionsRow(i))), ' and ' , cell2mat(listOfEvents(whichReactionsCol(i)))];
    end

    xlswrite(['ReactionsWithCorrelationCoefficientGreaterThanPoint', num2str(correlationThreshold*10)], correlatedReactions);
end

%////////////////////////////////////////////////////////////
% What is the extent of correlation between any two countries?
%////////////////////////////////////////////////////////////

[corrBetweenCountries, pValuesCountries] = corr(normalisedCountsPerCountry', 'type', 'Pearson');

% Set to zero any correlations with p value greater than 0.01
pIndices = find(pValuesCountries>0.01);
corrBetweenCountries(pIndices) = 0;
corrBetweenCountries(1:(size(corrBetweenCountries,1)+1):end) = 1;

% Plot heatmap of correlation
figure
imagesc(corrBetweenCountries);
colorbar;
title(['Correlation Between Countries (', num2str(limit), ' Reactions)'])
xlabel('Countries')
ylabel('Countries')

% Plot only those with high correlation (plotting bools)
correlationThreshold = 0.8;
figure
imagesc(corrBetweenCountries>correlationThreshold);
colorbar;
title(['Those countries with \rho > ', num2str(correlationThreshold)])
xlabel('Countries')
ylabel('Countries')

% Pull out the pairs of highly correlated countries (ignoring self-pairs)
[whichCountriesRow, whichCountriesCol] = find(triu(corrBetweenCountries)>correlationThreshold & triu(corrBetweenCountries)<1);
correlatedCountries = {};
if ~isempty(whichCountriesRow)
    % Store list of highly correlated countries 
    for i = 1:length(whichCountriesRow)

        rowIndex = strcmp(countryCodesToNames(:,2), sortedCountries{whichCountriesRow(i)});
        colIndex = strcmp(countryCodesToNames(:,2), sortedCountries{whichCountriesCol(i)});

        correlatedCountries{i,1} = [countryCodesToNames{rowIndex,1}, ' and ' , countryCodesToNames{colIndex,1}];
    end

    xlswrite(['CountriesWithCorrelationCoefficientGreaterThanPoint', num2str(correlationThreshold*10)], correlatedCountries);
end


%////////////////////////////////////////////////////////////
% Hierarchical Clustering 
%////////////////////////////////////////////////////////////

% Cluster along the columns so that we're clustering countries
% Matlab function that implements euclidean distance metric
clusterData = clustergram(countPerCountry, 'Standardize', 'row', 'Cluster', 'column','Colormap',redbluecmap);
set(clusterData,'RowLabels',sortedCountries)
colorbar
addXLabel(clusterData, 'Reaction')
addYLabel(clusterData, 'Country')
addTitle(clusterData, ['Hierarchical Clustering (', num2str(limit), ' Reactions)'])



