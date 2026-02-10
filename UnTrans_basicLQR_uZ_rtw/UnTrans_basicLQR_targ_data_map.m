    function targMap = targDataMap(),

    ;%***********************
    ;% Create Parameter Map *
    ;%***********************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 2;
        sectIdxOffset = 0;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc paramMap
        ;%
        paramMap.nSections           = nTotSects;
        paramMap.sectIdxOffset       = sectIdxOffset;
            paramMap.sections(nTotSects) = dumSection; %prealloc
        paramMap.nTotData            = -1;

        ;%
        ;% Auto data (UnTrans_basicLQR_P)
        ;%
            section.nData     = 68;
            section.data(68)  = dumData; %prealloc

                    ;% UnTrans_basicLQR_P.SFunction3_P1_Size
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

                    ;% UnTrans_basicLQR_P.Ts
                    section.data(2).logicalSrcIdx = 1;
                    section.data(2).dtTransOffset = 2;

                    ;% UnTrans_basicLQR_P.Translation_gain
                    section.data(3).logicalSrcIdx = 2;
                    section.data(3).dtTransOffset = 3;

                    ;% UnTrans_basicLQR_P.Rotation_gain
                    section.data(4).logicalSrcIdx = 3;
                    section.data(4).dtTransOffset = 4;

                    ;% UnTrans_basicLQR_P.CoulombViscousFriction2_gain
                    section.data(5).logicalSrcIdx = 4;
                    section.data(5).dtTransOffset = 5;

                    ;% UnTrans_basicLQR_P.CoulombViscousFriction1_gain
                    section.data(6).logicalSrcIdx = 5;
                    section.data(6).dtTransOffset = 6;

                    ;% UnTrans_basicLQR_P.CoulombViscousFriction2_offset
                    section.data(7).logicalSrcIdx = 6;
                    section.data(7).dtTransOffset = 7;

                    ;% UnTrans_basicLQR_P.CoulombViscousFriction1_offset
                    section.data(8).logicalSrcIdx = 7;
                    section.data(8).dtTransOffset = 8;

                    ;% UnTrans_basicLQR_P.NoCtrl_Value
                    section.data(9).logicalSrcIdx = 8;
                    section.data(9).dtTransOffset = 9;

                    ;% UnTrans_basicLQR_P.deg2rad_Gain
                    section.data(10).logicalSrcIdx = 9;
                    section.data(10).dtTransOffset = 11;

                    ;% UnTrans_basicLQR_P.Uprightanglecorrection_Value
                    section.data(11).logicalSrcIdx = 10;
                    section.data(11).dtTransOffset = 12;

                    ;% UnTrans_basicLQR_P.u_Value
                    section.data(12).logicalSrcIdx = 11;
                    section.data(12).dtTransOffset = 13;

                    ;% UnTrans_basicLQR_P.m2rad_Gain
                    section.data(13).logicalSrcIdx = 12;
                    section.data(13).dtTransOffset = 14;

                    ;% UnTrans_basicLQR_P.LQR_Gain
                    section.data(14).logicalSrcIdx = 13;
                    section.data(14).dtTransOffset = 15;

                    ;% UnTrans_basicLQR_P.Normal_Value
                    section.data(15).logicalSrcIdx = 14;
                    section.data(15).dtTransOffset = 27;

                    ;% UnTrans_basicLQR_P.Switch_Threshold
                    section.data(16).logicalSrcIdx = 15;
                    section.data(16).dtTransOffset = 28;

                    ;% UnTrans_basicLQR_P.Mu2_Gain
                    section.data(17).logicalSrcIdx = 16;
                    section.data(17).dtTransOffset = 29;

                    ;% UnTrans_basicLQR_P.UnitDelay_InitialCondition
                    section.data(18).logicalSrcIdx = 17;
                    section.data(18).dtTransOffset = 30;

                    ;% UnTrans_basicLQR_P.UnitDelay1_InitialCondition
                    section.data(19).logicalSrcIdx = 18;
                    section.data(19).dtTransOffset = 31;

                    ;% UnTrans_basicLQR_P.raw2g_Gain
                    section.data(20).logicalSrcIdx = 19;
                    section.data(20).dtTransOffset = 32;

                    ;% UnTrans_basicLQR_P.Gain4_Gain
                    section.data(21).logicalSrcIdx = 20;
                    section.data(21).dtTransOffset = 33;

                    ;% UnTrans_basicLQR_P.raw2dps_Gain
                    section.data(22).logicalSrcIdx = 21;
                    section.data(22).dtTransOffset = 34;

                    ;% UnTrans_basicLQR_P.UnitConversion_Gain
                    section.data(23).logicalSrcIdx = 22;
                    section.data(23).dtTransOffset = 35;

                    ;% UnTrans_basicLQR_P.Relay_OnVal
                    section.data(24).logicalSrcIdx = 23;
                    section.data(24).dtTransOffset = 36;

                    ;% UnTrans_basicLQR_P.Relay_OffVal
                    section.data(25).logicalSrcIdx = 24;
                    section.data(25).dtTransOffset = 37;

                    ;% UnTrans_basicLQR_P.Relay_YOn
                    section.data(26).logicalSrcIdx = 25;
                    section.data(26).dtTransOffset = 38;

                    ;% UnTrans_basicLQR_P.Relay_YOff
                    section.data(27).logicalSrcIdx = 26;
                    section.data(27).dtTransOffset = 39;

                    ;% UnTrans_basicLQR_P.Gain3_Gain
                    section.data(28).logicalSrcIdx = 27;
                    section.data(28).dtTransOffset = 40;

                    ;% UnTrans_basicLQR_P.Gain1_Gain
                    section.data(29).logicalSrcIdx = 28;
                    section.data(29).dtTransOffset = 41;

                    ;% UnTrans_basicLQR_P.Pulse2SI1_Gain
                    section.data(30).logicalSrcIdx = 29;
                    section.data(30).dtTransOffset = 42;

                    ;% UnTrans_basicLQR_P.Memory2_InitialCondition
                    section.data(31).logicalSrcIdx = 30;
                    section.data(31).dtTransOffset = 43;

                    ;% UnTrans_basicLQR_P.Saturation_UpperSat
                    section.data(32).logicalSrcIdx = 31;
                    section.data(32).dtTransOffset = 44;

                    ;% UnTrans_basicLQR_P.Saturation_LowerSat
                    section.data(33).logicalSrcIdx = 32;
                    section.data(33).dtTransOffset = 45;

                    ;% UnTrans_basicLQR_P.DiscreteFilter_NumCoef
                    section.data(34).logicalSrcIdx = 33;
                    section.data(34).dtTransOffset = 46;

                    ;% UnTrans_basicLQR_P.DiscreteFilter_DenCoef
                    section.data(35).logicalSrcIdx = 34;
                    section.data(35).dtTransOffset = 48;

                    ;% UnTrans_basicLQR_P.DiscreteFilter_InitialStates
                    section.data(36).logicalSrcIdx = 35;
                    section.data(36).dtTransOffset = 50;

                    ;% UnTrans_basicLQR_P.Pulse2SI2_Gain
                    section.data(37).logicalSrcIdx = 36;
                    section.data(37).dtTransOffset = 51;

                    ;% UnTrans_basicLQR_P.Memory2_InitialCondition_i
                    section.data(38).logicalSrcIdx = 37;
                    section.data(38).dtTransOffset = 52;

                    ;% UnTrans_basicLQR_P.Saturation_UpperSat_i
                    section.data(39).logicalSrcIdx = 38;
                    section.data(39).dtTransOffset = 53;

                    ;% UnTrans_basicLQR_P.Saturation_LowerSat_f
                    section.data(40).logicalSrcIdx = 39;
                    section.data(40).dtTransOffset = 54;

                    ;% UnTrans_basicLQR_P.DiscreteFilter1_NumCoef
                    section.data(41).logicalSrcIdx = 40;
                    section.data(41).dtTransOffset = 55;

                    ;% UnTrans_basicLQR_P.DiscreteFilter1_DenCoef
                    section.data(42).logicalSrcIdx = 41;
                    section.data(42).dtTransOffset = 57;

                    ;% UnTrans_basicLQR_P.DiscreteFilter1_InitialStates
                    section.data(43).logicalSrcIdx = 42;
                    section.data(43).dtTransOffset = 59;

                    ;% UnTrans_basicLQR_P.Stepm_Value
                    section.data(44).logicalSrcIdx = 43;
                    section.data(44).dtTransOffset = 60;

                    ;% UnTrans_basicLQR_P.Stepdeg_Value
                    section.data(45).logicalSrcIdx = 44;
                    section.data(45).dtTransOffset = 61;

                    ;% UnTrans_basicLQR_P.Switch_Threshold_h
                    section.data(46).logicalSrcIdx = 45;
                    section.data(46).dtTransOffset = 62;

                    ;% UnTrans_basicLQR_P.Gain_Gain
                    section.data(47).logicalSrcIdx = 46;
                    section.data(47).dtTransOffset = 63;

                    ;% UnTrans_basicLQR_P.Saturation1_UpperSat
                    section.data(48).logicalSrcIdx = 47;
                    section.data(48).dtTransOffset = 64;

                    ;% UnTrans_basicLQR_P.Saturation1_LowerSat
                    section.data(49).logicalSrcIdx = 48;
                    section.data(49).dtTransOffset = 65;

                    ;% UnTrans_basicLQR_P.Saturation3_UpperSat
                    section.data(50).logicalSrcIdx = 49;
                    section.data(50).dtTransOffset = 66;

                    ;% UnTrans_basicLQR_P.Saturation3_LowerSat
                    section.data(51).logicalSrcIdx = 50;
                    section.data(51).dtTransOffset = 67;

                    ;% UnTrans_basicLQR_P.Reset_Value
                    section.data(52).logicalSrcIdx = 51;
                    section.data(52).dtTransOffset = 68;

                    ;% UnTrans_basicLQR_P.Constant1_Value
                    section.data(53).logicalSrcIdx = 52;
                    section.data(53).dtTransOffset = 69;

                    ;% UnTrans_basicLQR_P.Gain_Gain_n
                    section.data(54).logicalSrcIdx = 53;
                    section.data(54).dtTransOffset = 70;

                    ;% UnTrans_basicLQR_P.Gain2_Gain
                    section.data(55).logicalSrcIdx = 54;
                    section.data(55).dtTransOffset = 71;

                    ;% UnTrans_basicLQR_P.Sumi_gainval
                    section.data(56).logicalSrcIdx = 55;
                    section.data(56).dtTransOffset = 72;

                    ;% UnTrans_basicLQR_P.Sumi_IC
                    section.data(57).logicalSrcIdx = 56;
                    section.data(57).dtTransOffset = 73;

                    ;% UnTrans_basicLQR_P.Constant4_Value
                    section.data(58).logicalSrcIdx = 57;
                    section.data(58).dtTransOffset = 74;

                    ;% UnTrans_basicLQR_P.Sumi1_gainval
                    section.data(59).logicalSrcIdx = 58;
                    section.data(59).dtTransOffset = 75;

                    ;% UnTrans_basicLQR_P.Sumi1_IC
                    section.data(60).logicalSrcIdx = 59;
                    section.data(60).dtTransOffset = 76;

                    ;% UnTrans_basicLQR_P.uDLookupTable_tableData
                    section.data(61).logicalSrcIdx = 60;
                    section.data(61).dtTransOffset = 77;

                    ;% UnTrans_basicLQR_P.uDLookupTable_bp01Data
                    section.data(62).logicalSrcIdx = 61;
                    section.data(62).dtTransOffset = 85;

                    ;% UnTrans_basicLQR_P.Memory1_InitialCondition
                    section.data(63).logicalSrcIdx = 62;
                    section.data(63).dtTransOffset = 93;

                    ;% UnTrans_basicLQR_P.Constant1_Value_f
                    section.data(64).logicalSrcIdx = 63;
                    section.data(64).dtTransOffset = 94;

                    ;% UnTrans_basicLQR_P.Memory_InitialCondition
                    section.data(65).logicalSrcIdx = 64;
                    section.data(65).dtTransOffset = 95;

                    ;% UnTrans_basicLQR_P.Mu1_Gain
                    section.data(66).logicalSrcIdx = 65;
                    section.data(66).dtTransOffset = 96;

                    ;% UnTrans_basicLQR_P.Mu2_Gain_p
                    section.data(67).logicalSrcIdx = 66;
                    section.data(67).dtTransOffset = 97;

                    ;% UnTrans_basicLQR_P.Constant2_Value
                    section.data(68).logicalSrcIdx = 67;
                    section.data(68).dtTransOffset = 98;

            nTotData = nTotData + section.nData;
            paramMap.sections(1) = section;
            clear section

            section.nData     = 1;
            section.data(1)  = dumData; %prealloc

                    ;% UnTrans_basicLQR_P.ResetEncodersWatchdogandIntegra
                    section.data(1).logicalSrcIdx = 68;
                    section.data(1).dtTransOffset = 0;

            nTotData = nTotData + section.nData;
            paramMap.sections(2) = section;
            clear section


            ;%
            ;% Non-auto Data (parameter)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        paramMap.nTotData = nTotData;



    ;%**************************
    ;% Create Block Output Map *
    ;%**************************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 1;
        sectIdxOffset = 0;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc sigMap
        ;%
        sigMap.nSections           = nTotSects;
        sigMap.sectIdxOffset       = sectIdxOffset;
            sigMap.sections(nTotSects) = dumSection; %prealloc
        sigMap.nTotData            = -1;

        ;%
        ;% Auto data (UnTrans_basicLQR_B)
        ;%
            section.nData     = 24;
            section.data(24)  = dumData; %prealloc

                    ;% UnTrans_basicLQR_B.Encoders
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

                    ;% UnTrans_basicLQR_B.IMU
                    section.data(2).logicalSrcIdx = 1;
                    section.data(2).dtTransOffset = 2;

                    ;% UnTrans_basicLQR_B.PWMstatus
                    section.data(3).logicalSrcIdx = 2;
                    section.data(3).dtTransOffset = 8;

                    ;% UnTrans_basicLQR_B.IMUstatus
                    section.data(4).logicalSrcIdx = 3;
                    section.data(4).dtTransOffset = 9;

                    ;% UnTrans_basicLQR_B.theta
                    section.data(5).logicalSrcIdx = 4;
                    section.data(5).dtTransOffset = 10;

                    ;% UnTrans_basicLQR_B.DiscreteFilter
                    section.data(6).logicalSrcIdx = 5;
                    section.data(6).dtTransOffset = 11;

                    ;% UnTrans_basicLQR_B.fi
                    section.data(7).logicalSrcIdx = 6;
                    section.data(7).dtTransOffset = 12;

                    ;% UnTrans_basicLQR_B.DiscreteFilter1
                    section.data(8).logicalSrcIdx = 7;
                    section.data(8).dtTransOffset = 13;

                    ;% UnTrans_basicLQR_B.Sum
                    section.data(9).logicalSrcIdx = 8;
                    section.data(9).dtTransOffset = 14;

                    ;% UnTrans_basicLQR_B.Sum_i
                    section.data(10).logicalSrcIdx = 9;
                    section.data(10).dtTransOffset = 15;

                    ;% UnTrans_basicLQR_B.ResetEncodersWatchdogandIntegra
                    section.data(11).logicalSrcIdx = 10;
                    section.data(11).dtTransOffset = 16;

                    ;% UnTrans_basicLQR_B.PWMperiodsec
                    section.data(12).logicalSrcIdx = 11;
                    section.data(12).dtTransOffset = 17;

                    ;% UnTrans_basicLQR_B.Gain
                    section.data(13).logicalSrcIdx = 12;
                    section.data(13).dtTransOffset = 18;

                    ;% UnTrans_basicLQR_B.Gain2
                    section.data(14).logicalSrcIdx = 13;
                    section.data(14).dtTransOffset = 19;

                    ;% UnTrans_basicLQR_B.Sumi
                    section.data(15).logicalSrcIdx = 14;
                    section.data(15).dtTransOffset = 20;

                    ;% UnTrans_basicLQR_B.SFunction
                    section.data(16).logicalSrcIdx = 15;
                    section.data(16).dtTransOffset = 21;

                    ;% UnTrans_basicLQR_B.uDLookupTable
                    section.data(17).logicalSrcIdx = 16;
                    section.data(17).dtTransOffset = 22;

                    ;% UnTrans_basicLQR_B.SFunction2
                    section.data(18).logicalSrcIdx = 17;
                    section.data(18).dtTransOffset = 23;

                    ;% UnTrans_basicLQR_B.Constant2
                    section.data(19).logicalSrcIdx = 18;
                    section.data(19).dtTransOffset = 24;

                    ;% UnTrans_basicLQR_B.SFunction3
                    section.data(20).logicalSrcIdx = 19;
                    section.data(20).dtTransOffset = 25;

                    ;% UnTrans_basicLQR_B.MinMax
                    section.data(21).logicalSrcIdx = 20;
                    section.data(21).dtTransOffset = 34;

                    ;% UnTrans_basicLQR_B.MinMax1
                    section.data(22).logicalSrcIdx = 21;
                    section.data(22).dtTransOffset = 35;

                    ;% UnTrans_basicLQR_B.psi
                    section.data(23).logicalSrcIdx = 22;
                    section.data(23).dtTransOffset = 36;

                    ;% UnTrans_basicLQR_B.gyroOut
                    section.data(24).logicalSrcIdx = 23;
                    section.data(24).dtTransOffset = 37;

            nTotData = nTotData + section.nData;
            sigMap.sections(1) = section;
            clear section


            ;%
            ;% Non-auto Data (signal)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        sigMap.nTotData = nTotData;



    ;%*******************
    ;% Create DWork Map *
    ;%*******************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 4;
        sectIdxOffset = 1;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc dworkMap
        ;%
        dworkMap.nSections           = nTotSects;
        dworkMap.sectIdxOffset       = sectIdxOffset;
            dworkMap.sections(nTotSects) = dumSection; %prealloc
        dworkMap.nTotData            = -1;

        ;%
        ;% Auto data (UnTrans_basicLQR_DW)
        ;%
            section.nData     = 17;
            section.data(17)  = dumData; %prealloc

                    ;% UnTrans_basicLQR_DW.DiscreteFilter_states
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

                    ;% UnTrans_basicLQR_DW.DiscreteFilter1_states
                    section.data(2).logicalSrcIdx = 1;
                    section.data(2).dtTransOffset = 1;

                    ;% UnTrans_basicLQR_DW.Sumi_DSTATE
                    section.data(3).logicalSrcIdx = 2;
                    section.data(3).dtTransOffset = 2;

                    ;% UnTrans_basicLQR_DW.Sumi1_DSTATE
                    section.data(4).logicalSrcIdx = 3;
                    section.data(4).dtTransOffset = 3;

                    ;% UnTrans_basicLQR_DW.UnitDelay_DSTATE
                    section.data(5).logicalSrcIdx = 4;
                    section.data(5).dtTransOffset = 4;

                    ;% UnTrans_basicLQR_DW.UnitDelay1_DSTATE
                    section.data(6).logicalSrcIdx = 5;
                    section.data(6).dtTransOffset = 5;

                    ;% UnTrans_basicLQR_DW.Memory2_PreviousInput
                    section.data(7).logicalSrcIdx = 6;
                    section.data(7).dtTransOffset = 6;

                    ;% UnTrans_basicLQR_DW.Memory2_PreviousInput_j
                    section.data(8).logicalSrcIdx = 8;
                    section.data(8).dtTransOffset = 7;

                    ;% UnTrans_basicLQR_DW.Memory1_PreviousInput
                    section.data(9).logicalSrcIdx = 10;
                    section.data(9).dtTransOffset = 8;

                    ;% UnTrans_basicLQR_DW.Memory_PreviousInput
                    section.data(10).logicalSrcIdx = 11;
                    section.data(10).dtTransOffset = 9;

                    ;% UnTrans_basicLQR_DW.H
                    section.data(11).logicalSrcIdx = 12;
                    section.data(11).dtTransOffset = 10;

                    ;% UnTrans_basicLQR_DW.Q
                    section.data(12).logicalSrcIdx = 13;
                    section.data(12).dtTransOffset = 13;

                    ;% UnTrans_basicLQR_DW.R
                    section.data(13).logicalSrcIdx = 14;
                    section.data(13).dtTransOffset = 22;

                    ;% UnTrans_basicLQR_DW.P
                    section.data(14).logicalSrcIdx = 15;
                    section.data(14).dtTransOffset = 23;

                    ;% UnTrans_basicLQR_DW.I
                    section.data(15).logicalSrcIdx = 16;
                    section.data(15).dtTransOffset = 32;

                    ;% UnTrans_basicLQR_DW.x_est
                    section.data(16).logicalSrcIdx = 17;
                    section.data(16).dtTransOffset = 41;

                    ;% UnTrans_basicLQR_DW.q
                    section.data(17).logicalSrcIdx = 18;
                    section.data(17).dtTransOffset = 44;

            nTotData = nTotData + section.nData;
            dworkMap.sections(1) = section;
            clear section

            section.nData     = 2;
            section.data(2)  = dumData; %prealloc

                    ;% UnTrans_basicLQR_DW.Ctrl_PWORK.LoggedData
                    section.data(1).logicalSrcIdx = 19;
                    section.data(1).dtTransOffset = 0;

                    ;% UnTrans_basicLQR_DW.State_PWORK.LoggedData
                    section.data(2).logicalSrcIdx = 20;
                    section.data(2).dtTransOffset = 2;

            nTotData = nTotData + section.nData;
            dworkMap.sections(2) = section;
            clear section

            section.nData     = 1;
            section.data(1)  = dumData; %prealloc

                    ;% UnTrans_basicLQR_DW.MinMaxPeriod_SubsysRanBC
                    section.data(1).logicalSrcIdx = 21;
                    section.data(1).dtTransOffset = 0;

            nTotData = nTotData + section.nData;
            dworkMap.sections(3) = section;
            clear section

            section.nData     = 2;
            section.data(2)  = dumData; %prealloc

                    ;% UnTrans_basicLQR_DW.Relay_Mode
                    section.data(1).logicalSrcIdx = 22;
                    section.data(1).dtTransOffset = 0;

                    ;% UnTrans_basicLQR_DW.P_not_empty
                    section.data(2).logicalSrcIdx = 23;
                    section.data(2).dtTransOffset = 1;

            nTotData = nTotData + section.nData;
            dworkMap.sections(4) = section;
            clear section


            ;%
            ;% Non-auto Data (dwork)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        dworkMap.nTotData = nTotData;



    ;%
    ;% Add individual maps to base struct.
    ;%

    targMap.paramMap  = paramMap;
    targMap.signalMap = sigMap;
    targMap.dworkMap  = dworkMap;

    ;%
    ;% Add checksums to base struct.
    ;%


    targMap.checksum0 = 3106538297;
    targMap.checksum1 = 1166197426;
    targMap.checksum2 = 1437391692;
    targMap.checksum3 = 2965872559;

