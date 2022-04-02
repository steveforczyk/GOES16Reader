function MakeOpticalSizeMovieNameFromGOESDataFile(numfiles)
% This routine will make a movie name for the ABI-L2_CMI Cloud Top Phase
% the first name is a list of GOESFileNames
% Written By: Stephen Forczyk
% Created: April 1,2022
% Revised: ----
% Classification: Unclassified
global TempMovieName iaddLightning MapFormFactor numncfiles;



nframes=num2str(numfiles);
if(iaddLightning<1)
    TempMovieName=strcat('ParticleSize-',MapFormFactor,'-Fr',nframes,'.mp4');
else
    TempMovieName=strcat('ParticleSizewLightning-',MapFormFactor,'-Fr',nframes,'.mp4');
end
end

