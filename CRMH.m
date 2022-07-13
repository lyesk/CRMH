%% 0D Collisional-Radiative Model for Hydrogen

ne=5e11; %[cm^-3]
Te=4;

%% AMJUEL H.4
%Contains processes that involve collisions with electrons
%Electron impact ionization and dissociation
datasets=["H.4_2.1.5", "H.4_2.2.5g", "H.4_2.2.9", "H.4_2.2.10", "H.4_2.2.11", "H.4_2.2.12", "H.4_2.2.14"];

for i=1:length(datasets)
    inputFolder='/home/lyes/Documents/PhD/CRM_H/Eirene_data/mat/';
    inputFile=[inputFolder char(datasets(i)) '.mat'];
    reactionRates=importdata(inputFile);
    R(i)=interp1(reactionRates(:,1),reactionRates(:,2),Te);
end

H4_Matrix=ne*[-R(1),    2*R(2)+R(4)   , 0,    R(6)+2*R(7) ; ... 
                0  , -(R(2)+R(3)+R(4)), 0,         0      ; ...
               R(1),        R(4)      , 0,    2*R(5)+R(6) ; ...
                0  ,        R(3)      , 0, -(R(5)+R(6)+R(7))];

%% Recombination and association at the walls
%Everything that touches the walls becomes H2

c=3e8;
kB=1.38e-23;
m_proton=1.67e-27;
alpha=0.5;
RAID_Length=1.5; %[m]
RAID_Radius=0.2; %[m]

%If switching between H and D, change the mass here
m_atom=1; %In units of proton mass
m_molecule=2;

Wall_Matrix=zeros(4);

%Ions: they follow the magnetic field and recombine on axial end of RAID

vBohm_atom=sqrt(Te/(m_atom*1e9))*c; %Bohm velocity
vBohm_molecule=sqrt(Te/(m_molecule*1e9))*c;
v_eff_atom=vBohm_atom*alpha; %Factor alpha=1/2 from Tonello2021
v_eff_molecule=vBohm_molecule*alpha;

Wall_Recombination_Rate_atom=v_eff_atom/RAID_Length; %[1/s]
Wall_Recombination_Rate_molecule=v_eff_molecule/RAID_Length;

Wall_Matrix(3,3)=-Wall_Recombination_Rate_atom;
Wall_Matrix(2,3)=Wall_Recombination_Rate_atom/2; %Factor 1/2 from H^+ -> (1/2)*H_2

Wall_Matrix(4,4)=-Wall_Recombination_Rate_molecule;
Wall_Matrix(2,4)=Wall_Recombination_Rate_molecule;

%Atoms: they follow ballistic motion at thermal velocity and associate on radial surface of RAID
T_H=20000; %H atoms assumed to be at room temperature
vth=sqrt(kB*T_H/(m_atom*m_proton));
Wall_Association_Rate=vth/RAID_Radius;

Wall_Matrix(1,1)=-Wall_Association_Rate;
Wall_Matrix(2,1)=Wall_Association_Rate/2; %Factor 1/2 from H -> (1/2)*H_2

%% Solve for the steady-state densities
M=H4_Matrix+Wall_Matrix;
poplevels=null(M);
poplevels=poplevels/sum(poplevels);

%Calculate absolute densities by setting (nH^+ + nH_2^+)=ne
total_density=ne*1e6/(poplevels(3)+poplevels(4));
densities=poplevels*total_density;
