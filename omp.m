clear all;
warning off all;
%Variables
N=10:5:100;%measrement level
d=256; %dimension
n=2:2:90; %sparsty lrvrl
numTrials = 500;
%counter
numN = size(N, 2);
numd = size(d, 2);
numn = size(n, 2);
%Data Collection
numCorr = zeros(numN, numd, numn);
mostSpars = zeros(numN);
id = 1;
for iN=1:numN
in=1;
done=1;
keepGo=1;
while in <= numn && keepGo
for trial = 1:numTrials
tN = N(1, iN);
td = d(1, id);
tn = n(1, in);
%Set Matrix
Phi = randn(tN, td); %randn Normally distributed pseudorandom numbers
Phi = sign(Phi);
I = zeros(1,1);
%Set signal
z = randperm(td); % randperm Random permutation
supp=z(1:tn);% print the 2 starting column
v = zeros(td, 1);%v is the matrix of order 256*1 having zero elements in it
for t = 1:tn%print the values in sequence from 1 to value of tn get by previous for loop
v(z(t))=1;
end
x = Phi * v;
r = x;
%Run OMP
while length(I)-1 < tn
u = Phi' * r;
absu = abs(u);
[b, ix] = sort(absu, 'descend');
bestInd = ix(1);
bestVal = b(1);
%Update I
I(length(I)+1) = bestInd;
%Update the residual
PhiSubI = Phi(:, I(2));
for c=3:length(I)
if ~ismember(I(2:c-1),I(c))
PhiSubI(:,c-1) = Phi(:,I(c));
end
end
y = lscov(PhiSubI, x);
r = x - PhiSubI * y;
end
if ismember(supp, I)
numCorr(iN, id, in) = numCorr(iN, id, in) +1;
end
end
if numCorr(iN, id, in) / numTrials > 0.98 && done
mostSpars(iN) = tn;
else
done=0;
end
if numCorr(iN, id, in) <= 0.01
keepGo=0;
end
in = in +1;
end % n
end % N