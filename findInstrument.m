function arrayInstr = findInstrument(pks)

numInstruments = length(database);
numHarmonics = 8;

% arrayInstr = ["Piano", "Guitar", "Flute", "Pipe Organ", "Clarinet", "Saxophone", "Trumpet", "French Horn", "Tuba"]

% obtain common frequencies and compare their locations
pksT = zeros(1:numHarmonics);
% locsT = zeros(1,length(locs));

% figure out which harmonics should be set to 0
% for rowNum = 1:1:length(pks)
%     for j = 1:1:length(pks)
%         if (abs(locsT(j)-database(rowNum)(j)) < 20)
%         
%         end
%     end
% end

% take first numHarmonics peaks of pks and locs using loop above
if length(pks) >= numHarmonics
    pksT = pks(1:numHarmonics);
% fill with zeros if overtones aren't picked up/present
else
    pksT = pks(1:length(pks));
    pksT(length(pks)+1:numHarmonics) = zeros;
%     for j = length(pks)+1:1:numHarmonics
%         pksT(j) = 0;
%     end
end

pksD = database(1,1:8)';
% locsT = locs(1:8);
% locsD = database(1+numInstruments,1:8);

% array to hold mse of instruments
% MSElist = zeros(1:8);
diff = zeros(1:8);

% find the mean squared error for each instrument
MSElist = zeros(1:numHarmonics);
for rowNum = 1:1:numInstruments
    pksD = database(rowNum,1:8);
    for l = 1:1:8
        diff(l) = pksD(l) - pksT(l);
    end
%     locsD = database(1+numInstruments,1:8);
    D = abs(diff).^2;
    MSE = sum(D(:))/numel(pksD);
    MSElist(rowNum) = MSE;
end

% find the instrument with least mean squared error
minNum = min(MSElist);
minIdx = 0;
for m = 1:1:numInstruments
    if (MSElist(m) == minNum)
        minIdx = m;
    end
end
p = minIdx;

% print output
if (p==1)
        arrayInstr = 'Piano';
    elseif (p==2)
        arrayInstr = 'Guitar';
    elseif (p==3)
        arrayInstr = 'Flute';
    elseif (p==4)
        arrayInstr = 'Pipe Organ';
    elseif (p==5)
        arrayInstr = 'Clarinet';
    elseif (p==6)
        arrayInstr = 'Saxophone';
    elseif (p==7)
        arrayInstr = 'Trumpet';
    elseif (p==8)
        arrayInstr = 'French Horn';
    else
        arrayInstr = 'Tuba';
end

end