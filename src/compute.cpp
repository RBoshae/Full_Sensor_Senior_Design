void 
calculate_cp(int tempCnt, 
             double &CPMethane, 
             double &CPCO2, 
             double &CPEthane, 
             double &thermConMethane, 
             double &thermConEthane, 
             double &thermConCO2) 
{
    if (tempCnt == -20) {
        CPMethane = methane1 * gasvalue;
        CPCO2 = CO2neg20 * gasvalue;
        CPEthane = ethaneneg20 * gasvalue;
        thermConMethane = thermalConductMethaneneg20;
        thermConEthane = thermalConductEthaneneg20; 
        thermConCO2 = thermalConductCO2neg20; 
    }
    if (tempCnt == 0) {
        CPMethane = methane2 * gasvalue;
        CPCO2 = CO20 * gasvalue;
        CPEthane = ethane0 * gasvalue;
        thermConMethane = thermalConductMethane0;
        thermConEthane = thermalConductEthane0; 
        thermConCO2 = thermalConductCO2zero; 
    }
    if (tempCnt == 20) {
        CPMethane = methane3 * gasvalue;
        CPCO2 = CO220 * gasvalue;
        CPEthane = ethane20 * gasvalue; 
        thermConMethane = thermalConductMethane20;
        thermConEthane = thermalConductEthane20; 
        thermConCO2 = thermalConductCO220; 
    }
    if (tempCnt == 40) {
        CPMethane = methane4 * gasvalue; 
        CPCO2 = CO240 * gasvalue;
        CPEthane = ethane40 * gasvalue;
        thermConMethane = thermalConductMethane40;
        thermConEthane = thermalConductEthane40; 
        thermConCO2 = thermalConductCO240; 
    }
    if (tempCnt == 60) {
        CPMethane = methane5 * gasvalue;
        CPCO2 = CO260 * gasvalue;
        CPEthane = ethane60 * gasvalue; 
        thermConMethane = thermalConductMethane60;
        thermConEthane = thermalConductEthane60; 
        thermConCO2 = thermalConductCO260; 
    }
    if (tempCnt == 80) {
        CPMethane = methane6 * gasvalue; 
        CPCO2 = CO280 * gasvalue;
        CPEthane = ethane80 * gasvalue; 
        thermConMethane = thermalConductMethane80;
        thermConEthane = thermalConductEthane80; 
        thermConCO2 = thermalConductCO280; 
    }
}



int count = 1;

for(int i = -20; i <= 80; i = i + 20) {
    calculate_cp(i, CPMethane, CPCO2, CPEthane, thermConMethane, thermConEthane, thermConCO2);

    for(int pressure = 500; pressure <= 3000; pressure = pressure + 500) {

        for(int j = 0; j < v.size(); ++j) {
            obtain_combinations(j, CombinationFossil, CombinationAnaerobic, CombinationLandfill);

            MethaneRatio = texasPipleineMethane*percentTexas + 
                           rockyPipelineMethane*percentRocky + 
                           peruvianLNGMethane*percentPeruvian + 
                           highAssociated1Methane*percentHigh1 + 
                           highAssociated2Methane*percentHigh2;

            //Here the percent variable are the percentages the permutations should be generating

            EthaneRatio = texasPipleineEthane*percentTexas + 
                          rockyPipelineEthane*percentRocky + 
                          peruvianLNGEthane*percentPeruvian + 
                          highAssociated1Ethane*percentHigh1 + 
                          highAssociated2Ethane*percentHigh2;


            CO2Ratio = texasPipleineCO2*percentTexas + 
                       rockyPipelineCO2*percentRocky + 
                       peruvianLNGCO2*percentPeruvian + 
                       highAssociated1CO2*percentHigh1 + 
                       highAssociated2CO2percentHigh2;

            file << endl << " --------------------------------------------------------" << endl;

            file << "(" << count << ")" << endl;
            count++;

            temp = i + 273.15;

            file << "At " << temp - 273.15 << " degrees, the Cp of Methane is: " << CPMethane << endl;
            file << "At " << temp << " degrees, the Cp of CO2 is: " << CPCO2 << endl;
            file << "Pressure: " << pressure << endl;

            CO2Ratio = 100 - MethaneRatio - EthaneRatio;
            file << "The ratio is " << MethaneRatio << "% Methane and " 
                << CO2Ratio  << "% CO2 and " << EthaneRatio  << "% Ethane." << endl;

            CpMixture = (CPMethane * MethaneRatio) + 
                        (CPCO2 * CO2Ratio) + 
                        (CPEthane*EthaneRatio); // important we want this
            
            file << "The Cp of the mixture is : " << CpMixture << endl; 
            CpMixture = CpMixture/100 ;
            CvMixture = CpMixture - gasvalue; // important we want this
            MethaneRatio = MethaneRatio / 100;
            CO2Ratio = CO2Ratio / 100;
            EthaneRatio = EthaneRatio / 100;
            YouRatio = (CpMixture / CvMixture);

            avgmolecularweight = MethaneRatio * molarmassmethane + 
                                 CO2Ratio * molarmassCO2 + 
                                 EthaneRatio * molarmassethane;

            specificgravity = avgmolecularweight / airmolecularweight; // important we want this
            file << "The specific gravity of the mixture is: " << specificgravity << endl;
            z=1;
            mixtureDensity = ((avgmolecularweight * pressure) / 
                             (z*gasDensityvalue* temp)) * sidensityconvert; // important we want this
            file << "The molar density is : " << mixtureDensity << endl;
            ThermalConducitivyMixture = ((thermConMethane * MethaneRatio) + 
                                         (thermConEthane * EthaneRatio) + 
                                         ( thermConCO2 * CO2Ratio))/1000; // important we want this

            soundVelocity = sqrt((YouRatio*z*idealgasvalue*temp)/(avgmolecularweight));
            file << "The Thermal Conductivity of the Mixture is : " << ThermalConducitivyMixture; // important we want this
            file << endl << " --------------------------------------------------------" << endl;

            file1 << i << ", " << pressure << ", " << count << ", " 
                  << CombinationFossil << ", " << CombinationAnaerobic << ", " 
                  << CombinationLandfill;
            file1 << ", , " << MethaneRatio * 100 << ", " << CO2Ratio * 100 
                  << ", " << EthaneRatio * 100 << ", , " << ThermalConducitivyMixture 
                  << ", " << CpMixture << ", " << CvMixture 
                  << ", " << mixtureDensity << ", " << specificgravity 
                  << ", " << soundVelocity << endl; 
            }
        }
    }
}







