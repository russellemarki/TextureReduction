function[TDI_scl] = grabTDIScale(CS, maxL)
% Scale the TDI max possible distance in an ODF

%new
%TricTric 0.6*log(maxL) + 4/3
%OrthTric (1/3)*log(maxL) + 0.5
%HexTric  (1/6)*log(maxL) + 0.4
%CubTric  (2/15)*log(maxL) + 0.15
%OrthOrth 0.3*log(maxL) + 0.52
%CubOrth  (2/15)*log(maxL) + 0.15
%HexOrth  (1/7)*log(maxL) + 0.15

if ~isnumeric(CS)
    CS = CS.id;
end

% if CS == 2
%     TDI_scl = 1.0/(0.6*log(maxL) + 4/3);
% elseif CS == 16
%     TDI_scl = 1.0/((1/3)*log(maxL) + 0.5);
% elseif CS== 40
%     TDI_scl = 1.0/((1/6)*log(maxL) + 0.4);
% elseif CS == 45
%     TDI_scl = 1.0/((2/15)*log(maxL) + 0.15);
% else
%     warning('No TDI scale for symmetry')
%     TDI_scl = 1.0;
% end


if CS == 2
    TDI_scl = 1.0/(1.9215*log(maxL + 2.1429) - 1.4333);
elseif CS == 16
    TDI_scl = 1.0/(0.4804*log(maxL + 0.7646) + 0.1382);
elseif CS== 40
    TDI_scl = 1.0/(0.1601*log(maxL + 0.2728) + 0.2045);
elseif CS == 45
    TDI_scl = 1.0/(0.1069*log(maxL + 2.0585) - 0.0776);
else
    warning('No TDI scale for symmetry')
    TDI_scl = 1.0;
end
end