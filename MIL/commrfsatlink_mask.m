function commrfsatlink_mask(selLvlBO, selTemp, selDoppler,...
    selImbal, selNoPh, selDcOffset, selCorrDoppler, ...
    selCorrIQ, altitude,frequency, antDiam, samplesPerSymbol, filterSpan,selPreDistortion, selADC)
% COMMRFSATLINK_MASK Sets up workspace variables for the RF Link example.
%

% Copyright 1996-2018 The MathWorks, Inc.

    paramRFSatLink.Altitude = altitude;
    paramRFSatLink.Frequency = frequency;
    paramRFSatLink.SamplesPerSymbol = samplesPerSymbol;
    paramRFSatLink.FilterSpan = filterSpan;
    paramRFSatLink.preDistortion = selPreDistortion;
    rfSatSourceSampleTime = 1e-5;
    rfSatSourceSamplesPerSymbol = 400;
    paramRFSatLink.resetBER = 0;
    
    % Update the input gain and output gain for HPA block (memoryless
    % nonlinearity)
    [paramRFSatLink.GindB,paramRFSatLink.GoutdB] = updateHPA(selLvlBO);

    % Icon code change on RX Thermal Noise block
    tmpTemp = [0 20 290 500];
    paramRFSatLink.RXTemp = tmpTemp(selTemp);
    
    % Update the doppler impairment parameter
    dopOffset = [0 0.7 3];
    paramRFSatLink.DoppOffset = dopOffset(selDoppler);
    
    % Imbalances: Amp, Phase, I DC Off, Q DC Off
    tmp = [3 20 1e-8 5e-8];
    paramRFSatLink.IQImbal = zeros(1,4);
    if selImbal > 1
        paramRFSatLink.IQImbal(selImbal-1) = tmp(selImbal-1);
    end
    
    tmpNoise = [-100 -55 -48];
    rfSatPhaseNoise = tmpNoise(selNoPh);
    
    % Enable/disable DC blocker
    paramRFSatLink.DCBlock = logical(selDcOffset); 
    
    % Enable/disable carrier synchronizer
    paramRFSatLink.CarrierSync = logical(selCorrDoppler);
    
    % Enable/disable I/Q imbalance compensation
    paramRFSatLink.IQComp = logical(selCorrIQ);
    
    % Set antenna gains
    neff= 0.55; % middle of the road efficiency
    % gain = sqrt(neff)*pi*antDiam*(frequency .* 1e6)/3e8; % Proakis pp 316 
    paramRFSatLink.TXAntGain_Downlink = 10^((1.7) / 10);
    paramRFSatLink.RXAntGain_Downlink = 10^((11 -(10+4+1+4)) / 10);
    
    paramRFSatLink.TXAntGain_Uplink = 10^((11) / 10);
    paramRFSatLink.RXAntGain_Uplink = 10^((1.7 -(10+4+1+4)) / 10);

    paramRFSatLink.TXAntGain_ISL = 10^((1.7) / 10);
    paramRFSatLink.RXAntGain_ISL = 10^((1.7 -(10+4+1+4)) / 10);


    
    % Set ADC parameters
    switch selADC
      case 'No ADC'
        paramRFSatLink.ADC = false;
        paramRFSatLink.ADC_NumBits = 16;
      case 'High resolution (16 bits)'
        paramRFSatLink.ADC = true;
        paramRFSatLink.ADC_NumBits = 16;
      case 'Medium resolution (8 bits)'
        paramRFSatLink.ADC = true;
        paramRFSatLink.ADC_NumBits = 8;
      case 'Low resolution (3 bits)'
        paramRFSatLink.ADC = true;
        paramRFSatLink.ADC_NumBits = 3;
    end
    paramRFSatLink.ADC_VFS = 1.1;

    assignin('base', 'paramRFSatLink', paramRFSatLink);
    assignin('base', 'rfSatSourceSampleTime', rfSatSourceSampleTime);
    assignin('base', 'rfSatSourceSamplesPerSymbol', rfSatSourceSamplesPerSymbol);
    assignin('base', 'rfSatPhaseNoise', rfSatPhaseNoise)
    
%*********************************************************************
% Function Name:     updateHPA
% Description:       update nonlinear amplifier input and output gains      
%********************************************************************
function [GindB, GoutdB] = updateHPA(selLvlBO)

    % Update the saturation level parameters
    valsBO = [-30 -7 -1];       % values for backoff
    rrcComp = 20*log10(.38);    % compensation for RRC filter P2P power
    gainLin = 18;               % fixed HPA linear gain
    alpha = 20*log10(2.1587);   % difference between linear gain and
                                % small signal gain
    sps = 8;                    % Sample per symbol
    rctGain = 10*log10(sps);    % Raised Cosine Filter Gain
    
    gainIP = valsBO - rrcComp - rctGain;
    GindB = gainIP(selLvlBO);
    gainOP = -valsBO + rrcComp - alpha + gainLin;
    GoutdB = gainOP(selLvlBO);
% EOF
