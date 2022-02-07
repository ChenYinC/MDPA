function AddCircAtomicwithCent(numA, rad, centX, centY, cirAtomType, InitDist, Temp)
global C
global x y AtomSpacing
global nAtoms
global AtomType Vx Vy Mass0 Mass1

Types = zeros(1, numA+1);
Types(2:length(Types)) = cirAtomType;
Types(1) = ~(cirAtomType);
X0 = zeros(1, numA+1);
Y0 = zeros(1, numA+1);
X0(1) = centX;
Y0(1) = centY;
theta = linspace(0, 2*pi, numA);
for i=1:numA
    X0(i+1) = rad*AtomSpacing*cos(theta(i));
    Y0(i+1) = rad*AtomSpacing*sin(theta(i));
end

VX0 = zeros(1, numA+1);
VY0 = zeros(1, numA+1);
xp = X0;
yp = Y0;
numAtoms = numA+1;

x(nAtoms + 1:nAtoms+numAtoms) = xp + (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist;
y(nAtoms + 1:nAtoms+numAtoms) = yp + (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist;

iTypes = Types == 1;
Mass = ones(1,numAtoms)*Mass0;
Mass(iTypes) = Mass1;

AtomType(nAtoms + 1:nAtoms + numAtoms) = Types;


if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;
    Vy(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb * Temp ./ Mass);

    Vx(nAtoms + 1:nAtoms + numAtoms) = std0 .* randn(1, numAtoms);
    Vy(nAtoms + 1:nAtoms + numAtoms) = std0 .* randn(1, numAtoms);
end

Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms)) + VX0;
Vy(nAtoms + 1:nAtoms + numAtoms) = Vy(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vy(nAtoms + 1:nAtoms + numAtoms)) + VY0;

nAtoms = nAtoms + numAtoms;

end