function [hklIntensitiesMatrix] = hkl_intensities(etaVector,sVector,hklRandomIntensitiesMatrix,PreferentialOrientationMatrix)
% HKL_INTENSITIES Returns the tridimensional array formed by the calculated hkl intensities
%   [hklIntensitiesMatrix] = hkl_intensities(etaVector,sVector,hklRandomIntensitiesMatrix,PreferentialOrientationMatrix)
% Input 
%   etaVector: row vector formed by the list of values for eta   
%   sVector : column vector formed by the list of values for s 
%   hklRandomIntensitiesMatrix : matrix formed by the calculated hkl random intensities.Each hklRandomIntensitiesMatrix column corresponds to a hkl plane, and each line to a s value.
%   PreferentialOrientationMatrix : matrix formed by the preferential orientation coefficients. Each PreferentialOrientationMatrix column corresponds to an eta value, and each line to a hkl plane.
% Output
%   hklIntensitiesMatrix : tridimensional array formed by the calculated hkl intensities. Each hklIntensitiesMatrix column corresponds to an eta value, each line to a s value and the third dimension corresponds to the hkl planes.
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
hklRandomIntensitiesMatrix = repmat(hklRandomIntensitiesMatrix,[1 size(etaVector,2) 1]);
PreferentialOrientationMatrix = repmat(PreferentialOrientationMatrix,[size(sVector,1) 1 1]); 

hklIntensitiesMatrix = hklRandomIntensitiesMatrix.*PreferentialOrientationMatrix;
end

