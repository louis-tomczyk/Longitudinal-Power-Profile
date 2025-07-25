function [Eadc,Axis]   = adc(Eelecfilt,Axis,rx)

    % codes compatibility
    ADCsps  = rx.ADC.Nsps;
    Nsps    = Axis.Nsps;
    Nsymb   = Axis.Nsymb;
    Tsymb   = Axis.Tsymb;

    % new axis
    Nsamp_adc       = ADCsps*Nsymb; 
    te              = Tsymb/ADCsps; 
    Be              = 1/te;
    
    f_rec           = (-Nsamp_adc/2+1:Nsamp_adc/2)/Nsamp_adc*Be;
    time_rec        = (0:(Nsamp_adc-1))*te;

    Axis.ADC.freq   = f_rec;
    Axis.ADC.time   = time_rec;
    Axis.ADC.df     = abs(f_rec(2)-f_rec(1));
    Axis.ADC.dt     = time_rec(2)-time_rec(1);
    Axis.ADC.sps    = ADCsps;
    Axis.ADC.Nsamp  = Nsamp_adc;

    % downsampling
    if size(Eelecfilt.field,2) == 1
        Xout    = Eelecfilt.field(:,1);
        Rx_ADC  = downsample(Xout, floorNsps/ADCsps));
    else
        Xout    = Eelecfilt.field(:,1);
        Yout    = Eelecfilt.field(:,2);

        Rx_ADC_X= downsample(Xout, floor(Nsps/ADCsps));
        Rx_ADC_Y= downsample(Yout, floor(Nsps/ADCsps));
        Rx_ADC  = [Rx_ADC_X,Rx_ADC_Y];
    end

    Eadc        = Eelecfilt;
    Eadc.field  = Rx_ADC;

end