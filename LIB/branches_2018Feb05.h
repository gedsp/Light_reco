/*
Library with the variable definition to read from the charge-light file
******************************************************************************
*Tree    :anatree   : analysis tree                                          *
*Entries :      335 : Total =       939430623 bytes  File  Size =  260494788 *
*        :          : Tree compression factor =   3.61                       *
******************************************************************************
*Br    0 :Run       : run/I                                                  *
*Entries :      335 : Total  Size=       1879 bytes  File Size  =        110 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  12.85     *
*............................................................................*
*Br    1 :Subrun    : subrun/I                                               *
*Entries :      335 : Total  Size=       1894 bytes  File Size  =        109 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  12.99     *
*............................................................................*
*Br    2 :EventNumberInRun : event/I                                         *
*Entries :      335 : Total  Size=       1911 bytes  File Size  =        624 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   2.29     *
*............................................................................*
*Br    3 :EventTimeSeconds : evttime_seconds/I                               *
*Entries :      335 : Total  Size=       1941 bytes  File Size  =        522 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   2.73     *
*............................................................................*
*Br    4 :EventTimeNanoseconds : evttime_nanoseconds/I                       *
*Entries :      335 : Total  Size=       1961 bytes  File Size  =       1430 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br    5 :IsData    : isdata/B                                               *
*Entries :      335 : Total  Size=        883 bytes  File Size  =         99 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   4.15     *
*............................................................................*
*Br    6 :NumberOfHits : no_hits/I                                           *
*Entries :      335 : Total  Size=       1909 bytes  File Size  =        919 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.55     *
*............................................................................*
*Br    7 :Hit_TPC   : hit_tpc[no_hits]/S                                     *
*Entries :      335 : Total  Size=    1148227 bytes  File Size  =      10988 *
*Baskets :       40 : Basket Size=      32000 bytes  Compression= 104.39     *
*............................................................................*
*Br    8 :Hit_View  : hit_view[no_hits]/S                                    *
*Entries :      335 : Total  Size=    1148271 bytes  File Size  =      14160 *
*Baskets :       40 : Basket Size=      32000 bytes  Compression=  81.01     *
*............................................................................*
*Br    9 :Hit_Channel : hit_channel[no_hits]/S                               *
*Entries :      335 : Total  Size=    1148403 bytes  File Size  =     379173 *
*Baskets :       40 : Basket Size=      32000 bytes  Compression=   3.03     *
*............................................................................*
*Br   10 :Hit_PeakTime : hit_peakT[no_hits]/F                                *
*Entries :      335 : Total  Size=    2294802 bytes  File Size  =    2044542 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   1.12     *
*............................................................................*
*Br   11 :Hit_ChargeSummedADC : hit_chargesum[no_hits]/F                     *
*Entries :      335 : Total  Size=    2295360 bytes  File Size  =    2102675 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   1.09     *
*............................................................................*
*Br   12 :Hit_ChargeIntegral : hit_chargeintegral[no_hits]/F                 *
*Entries :      335 : Total  Size=    2295297 bytes  File Size  =    2094223 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   1.10     *
*............................................................................*
*Br   13 :Hit_PeakHeight : hit_ph[no_hits]/F                                 *
*Entries :      335 : Total  Size=    2294949 bytes  File Size  =    2087165 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   1.10     *
*............................................................................*
*Br   14 :Hit_StartTime : hit_startT[no_hits]/F                              *
*Entries :      335 : Total  Size=    2294883 bytes  File Size  =     831876 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   2.76     *
*............................................................................*
*Br   15 :Hit_EndTime : hit_endT[no_hits]/F                                  *
*Entries :      335 : Total  Size=    2294721 bytes  File Size  =     847385 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   2.71     *
*............................................................................*
*Br   16 :Hit_Width : hit_rms[no_hits]/F                                     *
*Entries :      335 : Total  Size=    2294562 bytes  File Size  =    1153240 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   1.99     *
*............................................................................*
*Br   17 :Hit_GoodnessOfFit : hit_goodnessOfFit[no_hits]/F                   *
*Entries :      335 : Total  Size=    2295216 bytes  File Size  =    1286315 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   1.78     *
*............................................................................*
*Br   18 :Hit_Multiplicity : hit_multiplicity[no_hits]/S                     *
*Entries :      335 : Total  Size=    1148623 bytes  File Size  =     155870 *
*Baskets :       40 : Basket Size=      32000 bytes  Compression=   7.36     *
*............................................................................*
*Br   19 :Hit_TrackID : hit_trkid[no_hits]/S                                 *
*Entries :      335 : Total  Size=    1148397 bytes  File Size  =     159706 *
*Baskets :       40 : Basket Size=      32000 bytes  Compression=   7.18     *
*............................................................................*
*Br   20 :Hit_ClusterID : hit_clusterid[no_hits]/S                           *
*Entries :      335 : Total  Size=    1148491 bytes  File Size  =     300368 *
*Baskets :       40 : Basket Size=      32000 bytes  Compression=   3.82     *
*............................................................................*
*Br   21 :NumberOfClusters : nclusters/S                                     *
*Entries :      335 : Total  Size=       1249 bytes  File Size  =        578 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.31     *
*............................................................................*
*Br   22 :ClusterID : clusterId[nclusters]/S                                 *
*Entries :      335 : Total  Size=      42394 bytes  File Size  =       3708 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=  11.28     *
*............................................................................*
*Br   23 :Cluster_NumberOfHits : cluster_NHits[nclusters]/S                  *
*Entries :      335 : Total  Size=      42439 bytes  File Size  =      17877 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   2.34     *
*............................................................................*
*Br   24 :Cluster_View : clusterView[nclusters]/S                            *
*Entries :      335 : Total  Size=      42409 bytes  File Size  =       3157 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=  13.25     *
*............................................................................*
*Br   25 :Cluster_ChargeIntegral : cluster_Integral[nclusters]/F             *
*Entries :      335 : Total  Size=      82863 bytes  File Size  =      72775 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.13     *
*............................................................................*
*Br   26 :Cluster_ChargeIntegralAveragePerHit :                              *
*         | cluster_IntegralAverage[nclusters]/F                             *
*Entries :      335 : Total  Size=      82936 bytes  File Size  =      71957 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.14     *
*............................................................................*
*Br   27 :Cluster_StartChannel : cluster_StartWire[nclusters]/S              *
*Entries :      335 : Total  Size=      42451 bytes  File Size  =      29917 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.40     *
*............................................................................*
*Br   28 :Cluster_StartTick : cluster_StartTick[nclusters]/S                 *
*Entries :      335 : Total  Size=      42442 bytes  File Size  =      33490 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.25     *
*............................................................................*
*Br   29 :Cluster_EndChannel : cluster_EndWire[nclusters]/S                  *
*Entries :      335 : Total  Size=      42439 bytes  File Size  =      30365 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.38     *
*............................................................................*
*Br   30 :Cluster_EndTick : cluster_EndTick[nclusters]/S                     *
*Entries :      335 : Total  Size=      42430 bytes  File Size  =      33505 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   1.25     *
*............................................................................*
*Br   31 :Cluster_StartCharge : cluster_StartCharge[nclusters]/F             *
*Entries :      335 : Total  Size=      82860 bytes  File Size  =      73142 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.12     *
*............................................................................*
*Br   32 :Cluster_StartAngle : cluster_StartAngle[nclusters]/F               *
*Entries :      335 : Total  Size=      82853 bytes  File Size  =      74941 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.10     *
*............................................................................*
*Br   33 :Cluster_EndCharge : cluster_EndCharge[nclusters]/F                 *
*Entries :      335 : Total  Size=      82846 bytes  File Size  =      74401 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.11     *
*............................................................................*
*Br   34 :Cluster_EndAngle : cluster_EndAngle[nclusters]/F                   *
*Entries :      335 : Total  Size=      82839 bytes  File Size  =      74399 *
*Baskets :        3 : Basket Size=      32000 bytes  Compression=   1.11     *
*............................................................................*
*Br   35 :NumberOfTracks_pmtrack : NumberOfTracks_pmtrack/S                  *
*Entries :      335 : Total  Size=       1300 bytes  File Size  =        437 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.74     *
*............................................................................*
*Br   36 :TrackID_pmtrack : TrackID_pmtrack[NumberOfTracks_pmtrack]/S        *
*Entries :      335 : Total  Size=       8099 bytes  File Size  =       1773 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   4.20     *
*............................................................................*
*Br   37 :Track_NumberOfHits_pmtrack :                                       *
*         | Track_NumberOfHits_pmtrack[NumberOfTracks_pmtrack]/S             *
*Entries :      335 : Total  Size=       8154 bytes  File Size  =       4323 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.72     *
*............................................................................*
*Br   38 :Track_Length_pmtrack :                                             *
*         | Track_Length_pmtrack[NumberOfTracks_pmtrack]/F                   *
*Entries :      335 : Total  Size=      14143 bytes  File Size  =      11926 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.13     *
*............................................................................*
*Br   39 :Track_StartPoint_X_pmtrack :                                       *
*         | Track_StartPoint_X_pmtrack[NumberOfTracks_pmtrack]/F             *
*Entries :      335 : Total  Size=      14173 bytes  File Size  =      12048 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.12     *
*............................................................................*
*Br   40 :Track_StartPoint_Y_pmtrack :                                       *
*         | Track_StartPoint_Y_pmtrack[NumberOfTracks_pmtrack]/F             *
*Entries :      335 : Total  Size=      14173 bytes  File Size  =      12056 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.12     *
*............................................................................*
*Br   41 :Track_StartPoint_Z_pmtrack :                                       *
*         | Track_StartPoint_Z_pmtrack[NumberOfTracks_pmtrack]/F             *
*Entries :      335 : Total  Size=      14173 bytes  File Size  =      11733 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.15     *
*............................................................................*
*Br   42 :Track_StartPoint_DistanceToBoundary_pmtrack :                      *
*         | Track_StartPoint_DistanceToBoundary_pmtrack[NumberOfTracks_pmtrack]/F*
*Entries :      335 : Total  Size=      14258 bytes  File Size  =      11909 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.13     *
*............................................................................*
*Br   43 :Track_EndPoint_X_pmtrack :                                         *
*         | Track_EndPoint_X_pmtrack[NumberOfTracks_pmtrack]/F               *
*Entries :      335 : Total  Size=      14163 bytes  File Size  =      12033 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.12     *
*............................................................................*
*Br   44 :Track_EndPoint_Y_pmtrack :                                         *
*         | Track_EndPoint_Y_pmtrack[NumberOfTracks_pmtrack]/F               *
*Entries :      335 : Total  Size=      14163 bytes  File Size  =      12063 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.12     *
*............................................................................*
*Br   45 :Track_EndPoint_Z_pmtrack :                                         *
*         | Track_EndPoint_Z_pmtrack[NumberOfTracks_pmtrack]/F               *
*Entries :      335 : Total  Size=      14163 bytes  File Size  =      11656 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.15     *
*............................................................................*
*Br   46 :Track_EndPoint_DistanceToBoundary_pmtrack :                        *
*         | Track_EndPoint_DistanceToBoundary_pmtrack[NumberOfTracks_pmtrack]/F*
*Entries :      335 : Total  Size=      14248 bytes  File Size  =      11924 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.13     *
*............................................................................*
*Br   47 :Track_StartDirection_Theta_pmtrack :                               *
*         | Track_StartDirection_Theta_pmtrack[NumberOfTracks_pmtrack]/F     *
*Entries :      335 : Total  Size=      14213 bytes  File Size  =      11792 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.14     *
*............................................................................*
*Br   48 :Track_StartDirection_Phi_pmtrack :                                 *
*         | Track_StartDirection_Phi_pmtrack[NumberOfTracks_pmtrack]/F       *
*Entries :      335 : Total  Size=      14203 bytes  File Size  =      12052 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.12     *
*............................................................................*
*Br   49 :Track_StartDirection_X_pmtrack :                                   *
*         | Track_StartDirection_X_pmtrack[NumberOfTracks_pmtrack]/F         *
*Entries :      335 : Total  Size=      14193 bytes  File Size  =      11901 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.13     *
*............................................................................*
*Br   50 :Track_StartDirection_Y_pmtrack :                                   *
*         | Track_StartDirection_Y_pmtrack[NumberOfTracks_pmtrack]/F         *
*Entries :      335 : Total  Size=      14193 bytes  File Size  =      12121 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.11     *
*............................................................................*
*Br   51 :Track_StartDirection_Z_pmtrack :                                   *
*         | Track_StartDirection_Z_pmtrack[NumberOfTracks_pmtrack]/F         *
*Entries :      335 : Total  Size=      14193 bytes  File Size  =      11905 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.13     *
*............................................................................*
*Br   52 :Track_EndDirection_Theta_pmtrack :                                 *
*         | Track_EndDirection_Theta_pmtrack[NumberOfTracks_pmtrack]/F       *
*Entries :      335 : Total  Size=      14203 bytes  File Size  =      11792 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.14     *
*............................................................................*
*Br   53 :Track_EndDirection_Phi_pmtrack :                                   *
*         | Track_EndDirection_Phi_pmtrack[NumberOfTracks_pmtrack]/F         *
*Entries :      335 : Total  Size=      14193 bytes  File Size  =      12078 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.11     *
*............................................................................*
*Br   54 :Track_EndDirection_X_pmtrack :                                     *
*         | Track_EndDirection_X_pmtrack[NumberOfTracks_pmtrack]/F           *
*Entries :      335 : Total  Size=      14183 bytes  File Size  =      11889 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.13     *
*............................................................................*
*Br   55 :Track_EndDirection_Y_pmtrack :                                     *
*         | Track_EndDirection_Y_pmtrack[NumberOfTracks_pmtrack]/F           *
*Entries :      335 : Total  Size=      14183 bytes  File Size  =      12141 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.11     *
*............................................................................*
*Br   56 :Track_EndDirection_Z_pmtrack :                                     *
*         | Track_EndDirection_Z_pmtrack[NumberOfTracks_pmtrack]/F           *
*Entries :      335 : Total  Size=      14183 bytes  File Size  =      11932 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.13     *
*............................................................................*
*Br   57 :Track_PitchInViews_pmtrack :                                       *
*         | Track_PitchInViews_pmtrack[NumberOfTracks_pmtrack][2]/F          *
*Entries :      335 : Total  Size=      26195 bytes  File Size  =      22773 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.12     *
*............................................................................*
*Br   58 :Track_NumberOfHitsPerView_pmtrack :                                *
*         | Track_NumberOfHitsPerView_pmtrack[NumberOfTracks_pmtrack][2]/S   *
*Entries :      335 : Total  Size=      14203 bytes  File Size  =       6912 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.95     *
*............................................................................*
*Br   59 :Track_Hit_X_pmtrack : Track_Hit_X_pmtrack[no_hits]/F               *
*Entries :      335 : Total  Size=    2295378 bytes  File Size  =     665005 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   3.45     *
*............................................................................*
*Br   60 :Track_Hit_Y_pmtrack : Track_Hit_Y_pmtrack[no_hits]/F               *
*Entries :      335 : Total  Size=    2295378 bytes  File Size  =     662944 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   3.46     *
*............................................................................*
*Br   61 :Track_Hit_Z_pmtrack : Track_Hit_Z_pmtrack[no_hits]/F               *
*Entries :      335 : Total  Size=    2295378 bytes  File Size  =     640023 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   3.58     *
*............................................................................*
*Br   62 :Track_Hit_dx_LocalTrackDirection_pmtrack :                         *
*         | Track_Hit_dx_LocalTrackDirection_pmtrack[no_hits]/F              *
*Entries :      335 : Total  Size=    2297079 bytes  File Size  =     178841 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=  12.83     *
*............................................................................*
*Br   63 :Track_Hit_dx_3DPosition_pmtrack :                                  *
*         | Track_Hit_dx_3DPosition_pmtrack[no_hits]/F                       *
*Entries :      335 : Total  Size=    2296350 bytes  File Size  =     628297 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   3.65     *
*............................................................................*
*Br   64 :Track_Hit_TPC_pmtrack : Track_Hit_TPC_pmtrack[no_hits]/S           *
*Entries :      335 : Total  Size=    1148843 bytes  File Size  =      14731 *
*Baskets :       40 : Basket Size=      32000 bytes  Compression=  77.90     *
*............................................................................*
*Br   65 :Track_Hit_View_pmtrack : Track_Hit_View_pmtrack[no_hits]/S         *
*Entries :      335 : Total  Size=    1148887 bytes  File Size  =      59021 *
*Baskets :       40 : Basket Size=      32000 bytes  Compression=  19.44     *
*............................................................................*
*Br   66 :Track_Hit_Channel_pmtrack : Track_Hit_Channel_pmtrack[no_hits]/S   *
*Entries :      335 : Total  Size=    1149019 bytes  File Size  =     249863 *
*Baskets :       40 : Basket Size=      32000 bytes  Compression=   4.59     *
*............................................................................*
*Br   67 :Track_Hit_PeakTime_pmtrack : Track_Hit_PeakTime_pmtrack[no_hits]/F *
*Entries :      335 : Total  Size=    2295945 bytes  File Size  =     640648 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   3.58     *
*............................................................................*
*Br   68 :Track_Hit_ChargeSummedADC_pmtrack :                                *
*         | Track_Hit_ChargeSummedADC_pmtrack[no_hits]/F                     *
*Entries :      335 : Total  Size=    2296512 bytes  File Size  =     653879 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   3.51     *
*............................................................................*
*Br   69 :Track_Hit_ChargeIntegral_pmtrack :                                 *
*         | Track_Hit_ChargeIntegral_pmtrack[no_hits]/F                      *
*Entries :      335 : Total  Size=    2296431 bytes  File Size  =     652839 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   3.51     *
*............................................................................*
*Br   70 :Track_Hit_PeakHeight_pmtrack :                                     *
*         | Track_Hit_PeakHeight_pmtrack[no_hits]/F                          *
*Entries :      335 : Total  Size=    2296107 bytes  File Size  =      50322 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=  45.59     *
*............................................................................*
*Br   71 :Track_Hit_StartTime_pmtrack :                                      *
*         | Track_Hit_StartTime_pmtrack[no_hits]/F                           *
*Entries :      335 : Total  Size=    2296026 bytes  File Size  =     342948 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   6.69     *
*............................................................................*
*Br   72 :Track_Hit_EndTime_pmtrack : Track_Hit_EndTime_pmtrack[no_hits]/F   *
*Entries :      335 : Total  Size=    2295864 bytes  File Size  =     355140 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   6.46     *
*............................................................................*
*Br   73 :Track_Hit_Width_pmtrack : Track_Hit_Width_pmtrack[no_hits]/F       *
*Entries :      335 : Total  Size=    2295702 bytes  File Size  =     543811 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   4.22     *
*............................................................................*
*Br   74 :Track_Hit_GoodnessOfFit_pmtrack :                                  *
*         | Track_Hit_GoodnessOfFit_pmtrack[no_hits]/F                       *
*Entries :      335 : Total  Size=    2296350 bytes  File Size  =     621454 *
*Baskets :       77 : Basket Size=      32000 bytes  Compression=   3.69     *
*............................................................................*
*Br   75 :Track_Hit_Multiplicity_pmtrack :                                   *
*         | Track_Hit_Multiplicity_pmtrack[no_hits]/S                        *
*Entries :      335 : Total  Size=    1149239 bytes  File Size  =      76665 *
*Baskets :       40 : Basket Size=      32000 bytes  Compression=  14.97     *
*............................................................................*
*Br   76 :light_run : light_run/I                                            *
*Entries :      335 : Total  Size=       1909 bytes  File Size  =        116 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  12.23     *
*............................................................................*
*Br   77 :light_event : light_event/I                                        *
*Entries :      335 : Total  Size=       1919 bytes  File Size  =        632 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   2.25     *
*............................................................................*
*Br   78 :light_nchannels : light_nchannels/I                                *
*Entries :      335 : Total  Size=       1939 bytes  File Size  =        122 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  11.68     *
*............................................................................*
*Br   79 :light_TimeSample : light_time_sample/I                             *
*Entries :      335 : Total  Size=       1947 bytes  File Size  =        123 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  11.59     *
*............................................................................*
*Br   80 :light_nsamples : light_nsamples/I                                  *
*Entries :      335 : Total  Size=       1934 bytes  File Size  =        121 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  11.77     *
*............................................................................*
*Br   81 :light_adc_value_0 : light_adc_value_0[light_nsamples]/S            *
*Entries :      335 : Total  Size=  175676810 bytes  File Size  =   44495272 *
*Baskets :      335 : Basket Size=      32000 bytes  Compression=   3.95     *
*............................................................................*
*Br   82 :light_adc_value_1 : light_adc_value_1[light_nsamples]/S            *
*Entries :      335 : Total  Size=  175676810 bytes  File Size  =   46181014 *
*Baskets :      335 : Basket Size=      32000 bytes  Compression=   3.80     *
*............................................................................*
*Br   83 :light_adc_value_2 : light_adc_value_2[light_nsamples]/S            *
*Entries :      335 : Total  Size=  175676810 bytes  File Size  =   52309483 *
*Baskets :      335 : Basket Size=      32000 bytes  Compression=   3.36     *
*............................................................................*
*Br   84 :light_adc_value_3 : light_adc_value_3[light_nsamples]/S            *
*Entries :      335 : Total  Size=  175676810 bytes  File Size  =   49544157 *
*Baskets :      335 : Basket Size=      32000 bytes  Compression=   3.55     *
*............................................................................*
*Br   85 :light_adc_value_4 : light_adc_value_4[light_nsamples]/S            *
*Entries :      335 : Total  Size=  175676810 bytes  File Size  =   46551178 *
*Baskets :      335 : Basket Size=      32000 bytes  Compression=   3.77     *
*............................................................................*
*Br   86 :light_PCTimeTag : light_PCTimeTag[3]/I                             *
*Entries :      335 : Total  Size=       4625 bytes  File Size  =       1688 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   2.43     *
*............................................................................*
*Br   87 :light_crt_daq_match : light_crt_daq_match/I                        *
*Entries :      335 : Total  Size=       1959 bytes  File Size  =        122 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  11.71     *
*............................................................................*
*Br   88 :light_crt_reco : light_crt_reco/I                                  *
*Entries :      335 : Total  Size=       1934 bytes  File Size  =        117 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  12.17     *
*............................................................................*
*Br   89 :crt_adc   : light_crt_adc[4][32]/I                                 *
*Entries :      335 : Total  Size=     172496 bytes  File Size  =       1405 *
*Baskets :        6 : Basket Size=      32000 bytes  Compression= 122.41     *
*............................................................................*
*Br   90 :crt_track_pos0 : crt_track_pos0[3]/F                               *
*Entries :      335 : Total  Size=       4620 bytes  File Size  =        134 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  30.63     *
*............................................................................*
*Br   91 :crt_track_pos1 : crt_track_pos1[3]/F                               *
*Entries :      335 : Total  Size=       4620 bytes  File Size  =        134 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  30.63     *
*............................................................................*
*Br   92 :crt_ToF   : crt_ToF/F                                              *
*Entries :      670 : Total  Size=       3239 bytes  File Size  =        121 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  22.79     *
*............................................................................*
*Br   93 :crt_isFV  : crt_isFV/I                                             *
*Entries :        0 : Total  Size=        486 bytes  One basket in memory    *
*Baskets :        0 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br   94 :crt_track_lenFV : crt_track_lenFV/F                                *
*Entries :      335 : Total  Size=       1939 bytes  File Size  =        118 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  12.08     *
*............................................................................*
*Br   95 :crt_point_in : crt_point_in[3]/F                                   *
*Entries :      335 : Total  Size=       4610 bytes  File Size  =        132 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  31.08     *
*............................................................................*
*Br   96 :crt_point_out : crt_point_out[3]/F                                 *
*Entries :      335 : Total  Size=       4615 bytes  File Size  =        133 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  30.85     *
*............................................................................*
*Br   97 :crt_pmt_dist : crt_pmt_dist[5]/F                                   *
*Entries :      335 : Total  Size=       7290 bytes  File Size  =        142 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  47.76     *
*............................................................................*

*/


	const int kmaxnhits=100000;
	const int kmaxntracks=1000;

	int _runcharge;
	int _subruncharge;
	int _event;
    	int _eventtime_seconds;
    	int _eventtime_nanoseconds;

    	//Charge variables:
	TBranch * _b2_runcharge;	  	t2->SetBranchAddress("Run"       		, &_runcharge       		, &_b2_runcharge       );
	TBranch * _b2_subruncharge;	  	t2->SetBranchAddress("Subrun"       		, &_subruncharge       		, &_b2_subruncharge       );
	TBranch * _b2_event;	  		t2->SetBranchAddress("EventNumberInRun"       		, &_event       		, &_b2_event       );
	TBranch * _b2_eventtime_seconds;	t2->SetBranchAddress("EventTimeSeconds"		, &_eventtime_seconds		, &_b2_eventtime_seconds       );
	TBranch * _b2_eventtime_nanoseconds;	t2->SetBranchAddress("EventTimeNanoseconds"	, &_eventtime_nanoseconds       , &_b2_eventtime_nanoseconds       );

	int _no_hits;
	float _Hit_ChargeIntegral[kmaxnhits];
	short _hit_trkid[kmaxnhits];
	short _NumberOfTracks_pmtrack;
	short _nclusters;

	float _Track_PitchInViews_pmtrack[kmaxntracks][2];
	short _Track_NumberOfHitsPerView_pmtrack[kmaxntracks][2];

	float _Track_TPC_pmtrack[kmaxnhits];
	float _Track_HitX_pmtrack[kmaxnhits];
	float _Track_HitY_pmtrack[kmaxnhits];
	float _Track_HitZ_pmtrack[kmaxnhits];
	float _Track_StartX_pmtrack[kmaxntracks];
	float _Track_StartY_pmtrack[kmaxntracks];
	float _Track_StartZ_pmtrack[kmaxntracks];

	float _Track_EndX_pmtrack[kmaxntracks];
	float _Track_EndY_pmtrack[kmaxntracks];
	float _Track_EndZ_pmtrack[kmaxntracks];
	float _Track_StartDirection_Theta_pmtrack[kmaxntracks];
	float _Track_Length_pmtrack[kmaxntracks];
	float _Track_StartDirection_Phi_pmtrack[kmaxntracks];
	float _Hit_PeakTime[kmaxnhits];
	short _Hit_Channel[kmaxnhits];
	short _Hit_View[kmaxnhits];
	short _Hit_TPC[kmaxnhits];

	TBranch * _b2_no_hits;			t2->SetBranchAddress("NumberOfHits"		, &_no_hits			, &_b2_no_hits       );
	TBranch * _b2_Hit_ChargeIntegral;	t2->SetBranchAddress("Hit_ChargeIntegral"	, _Hit_ChargeIntegral		, &_b2_Hit_ChargeIntegral      );
	TBranch * _b2_hit_trkid;		t2->SetBranchAddress("Hit_TrackID"		, _hit_trkid			, &_b2_hit_trkid      );
	TBranch * _b2_NumberOfTracks_pmtrack;	t2->SetBranchAddress("NumberOfTracks_pmtrack"	, &_NumberOfTracks_pmtrack	, &_b2_NumberOfTracks_pmtrack      );
	TBranch * _b2_nclusters;		t2->SetBranchAddress("NumberOfClusters"		, &_nclusters			, &_b2_nclusters       );
	TBranch * _b2_Track_PitchInViews_pmtrack;t2->SetBranchAddress("Track_PitchInViews_pmtrack"		, _Track_PitchInViews_pmtrack			, &_b2_Track_PitchInViews_pmtrack       );

	TBranch * _b2_Hit_TPC;			t2->SetBranchAddress("Hit_TPC"		, _Hit_TPC		, &_b2_Hit_TPC		);
	TBranch * _b2_Hit_View;			t2->SetBranchAddress("Hit_View"		, _Hit_View		, &_b2_Hit_View		);
	TBranch * _b2_Hit_Channel;		t2->SetBranchAddress("Hit_Channel"	, _Hit_Channel		, &_b2_Hit_Channel	);
	TBranch * _b2_Hit_PeakTime;		t2->SetBranchAddress("Hit_PeakTime"	, _Hit_PeakTime		, &_b2_Hit_PeakTime	);

	TBranch * _b2_Track_NumberOfHitsPerView_pmtrack;t2->SetBranchAddress("Track_NumberOfHitsPerView_pmtrack", _Track_NumberOfHitsPerView_pmtrack, &_b2_Track_NumberOfHitsPerView_pmtrack);

	
	TBranch * _b2_Track_Hit_ChargeSummedADC_pmtrack ;float _Track_Hit_ChargeSummedADC_pmtrack[kmaxnhits]; t2->SetBranchAddress("Track_Hit_ChargeSummedADC_pmtrack"		, _Track_Hit_ChargeSummedADC_pmtrack			, &_b2_Track_Hit_ChargeSummedADC_pmtrack       );
	TBranch * _b2_Track_Hit_ChargeIntegral_pmtrack ;float _Track_Hit_ChargeIntegral_pmtrack[kmaxnhits];t2->SetBranchAddress("Track_Hit_ChargeIntegral_pmtrack"		, _Track_Hit_ChargeIntegral_pmtrack			, &_b2_Track_Hit_ChargeIntegral_pmtrack       );


	TBranch * _b2_Track_Hit_dx_LocalTrackDirection_pmtrack ;float _Track_Hit_dx_LocalTrackDirection_pmtrack[kmaxnhits];t2->SetBranchAddress("Track_Hit_dx_LocalTrackDirection_pmtrack"		, _Track_Hit_dx_LocalTrackDirection_pmtrack			, &_b2_Track_Hit_ChargeIntegral_pmtrack       );
	TBranch * _b2_Track_Hit_dx_3DPosition_pmtrack ;float _Track_Hit_dx_3DPosition_pmtrack[kmaxnhits];t2->SetBranchAddress("Track_Hit_dx_3DPosition_pmtrack"		, _Track_Hit_dx_3DPosition_pmtrack			, &_b2_Track_Hit_dx_3DPosition_pmtrack       );


	TBranch * _b2_Track_TPC_pmtrack		;t2->SetBranchAddress("Track_Hit_TPC_pmtrack"		, _Track_TPC_pmtrack			, &_b2_Track_TPC_pmtrack       );
	TBranch * _b2_Track_HitX_pmtrack	;t2->SetBranchAddress("Track_Hit_X_pmtrack"		, _Track_HitX_pmtrack			, &_b2_Track_HitX_pmtrack       );
	TBranch * _b2_Track_HitY_pmtrack	;t2->SetBranchAddress("Track_Hit_Y_pmtrack"		, _Track_HitY_pmtrack			, &_b2_Track_HitY_pmtrack       );
	TBranch * _b2_Track_HitZ_pmtrack	;t2->SetBranchAddress("Track_Hit_Z_pmtrack"		, _Track_HitZ_pmtrack			, &_b2_Track_HitZ_pmtrack       );

	TBranch * _b2_Track_StartX_pmtrack	;t2->SetBranchAddress("Track_StartPoint_X_pmtrack"		, _Track_StartX_pmtrack			, &_b2_Track_StartX_pmtrack       );
	TBranch * _b2_Track_StartY_pmtrack	;t2->SetBranchAddress("Track_StartPoint_Y_pmtrack"		, _Track_StartY_pmtrack			, &_b2_Track_StartY_pmtrack       );
	TBranch * _b2_Track_StartZ_pmtrack	;t2->SetBranchAddress("Track_StartPoint_Z_pmtrack"		, _Track_StartZ_pmtrack			, &_b2_Track_StartZ_pmtrack       );

	TBranch * _b2_Track_EndX_pmtrack	;t2->SetBranchAddress("Track_EndPoint_X_pmtrack"		, _Track_EndX_pmtrack			, &_b2_Track_EndX_pmtrack       );
	TBranch * _b2_Track_EndY_pmtrack	;t2->SetBranchAddress("Track_EndPoint_Y_pmtrack"		, _Track_EndY_pmtrack			, &_b2_Track_EndY_pmtrack       );
	TBranch * _b2_Track_EndZ_pmtrack	;t2->SetBranchAddress("Track_EndPoint_Z_pmtrack"		, _Track_EndZ_pmtrack			, &_b2_Track_EndZ_pmtrack       );


	TBranch * _b2_Track_StartDirection_Theta_pmtrack	;t2->SetBranchAddress("Track_StartDirection_Theta_pmtrack"		, _Track_StartDirection_Theta_pmtrack			, &_b2_Track_StartDirection_Theta_pmtrack      );
	TBranch * _b2_Track_Length_pmtrack	;t2->SetBranchAddress("Track_Length_pmtrack"		, _Track_Length_pmtrack			, &_b2_Track_Length_pmtrack       );
	TBranch * _b2_Track_StartDirection_Phi_pmtrack	;	t2->SetBranchAddress("Track_StartDirection_Phi_pmtrack"		, _Track_StartDirection_Phi_pmtrack			, &_b2_Track_StartDirection_Phi_pmtrack       );

	//Light Variables
	int _PCTimeTag[3];
	int _nsamples;
	int _time_stamp;
    	int _nevent;
    	int _nchannels;
  	int _crt_daq_match;
  	int _crt_reco;
    	int _time_sample;
	int _run;
	short _adc_value[5][300000];

	TBranch * _b2_run;		t2->SetBranchAddress("light_run"		, &_run			, &_b2_run       	);
	TBranch * _b2_light_event;	t2->SetBranchAddress("light_event"		, &_nevent    		, &_b2_light_event	);
	TBranch * _b2_time_sample;	t2->SetBranchAddress("light_TimeSample"		, &_time_sample		, &_b2_time_sample   	);	
	TBranch * _b2_nchannels;	t2->SetBranchAddress("light_nchannels"		, &_nchannels		, &_b2_nchannels   	);
	TBranch * _b2_PCTimeTag;	t2->SetBranchAddress("light_PCTimeTag"  	, _PCTimeTag  		, &_b2_PCTimeTag	);

	TBranch * _b2_nsamples;		t2->SetBranchAddress("light_nsamples"   	, &_nsamples  		,  &_b2_nsamples	);	
	TBranch * _b2_adc_value_0; 	t2->SetBranchAddress("light_adc_value_0"	, _adc_value[0]  	, &_b2_adc_value_0 	);
	TBranch * _b2_adc_value_1; 	t2->SetBranchAddress("light_adc_value_1"	, _adc_value[1]  	, &_b2_adc_value_1 	);
	TBranch * _b2_adc_value_2; 	t2->SetBranchAddress("light_adc_value_2"	, _adc_value[2]  	, &_b2_adc_value_2 	);
	TBranch * _b2_adc_value_3;  	t2->SetBranchAddress("light_adc_value_3"	, _adc_value[3]  	, &_b2_adc_value_3 	);
	TBranch * _b2_adc_value_4;  	t2->SetBranchAddress("light_adc_value_4"	, _adc_value[4]  	, &_b2_adc_value_4 	);

	TBranch * _b2_crt_daq_match;	t2->SetBranchAddress("light_crt_daq_match", &_crt_daq_match	, &_b2_crt_daq_match		);
	TBranch * _b2_crt_reco;		t2->SetBranchAddress("light_crt_reco"	, &_crt_reco		, &_b2_crt_reco		);
	TBranch * _b2_crt_adc      	; int	_crt_adc[4][32]      ; t2->SetBranchAddress("crt_adc"       , _crt_adc          , &_b2_crt_adc       );
	TBranch * _b2_crt_track_pos0    ; float	_crt_track_pos0[3]     ; t2->SetBranchAddress("crt_track_pos0"       , _crt_track_pos0          , &_b2_crt_track_pos0       );
	TBranch * _b2_crt_track_pos1    ; float	_crt_track_pos1[3]     ; t2->SetBranchAddress("crt_track_pos1"       , _crt_track_pos1          , &_b2_crt_track_pos1       );
	TBranch * _b2_crt_ToF      	; float	_crt_ToF     ; t2->SetBranchAddress("crt_ToF"       , &_crt_ToF          , &_b2_crt_ToF       );


	TBranch * _b2_crt_isFV      	; int	_crt_isFV     		; t2->SetBranchAddress("crt_isFV"        , &_crt_isFV         , &_b2_crt_isFV        	);
	TBranch * _b2_crt_track_lenFV   ; float	_crt_track_lenFV     	; t2->SetBranchAddress("crt_track_lenFV" , &_crt_track_lenFV  , &_b2_crt_track_lenFV	);
	TBranch * _b2_crt_point_in      ; float	_crt_point_in[3]     	; t2->SetBranchAddress("crt_point_in"    , _crt_point_in      , &_b2_crt_point_in       );
	TBranch * _b2_crt_point_out     ; float _crt_point_out[3]     	; t2->SetBranchAddress("crt_point_out"    , _crt_point_out      , &_b2_crt_point_out    );
	TBranch * _b2_crt_pmt_dist     ; float _crt_pmt_dist[3]     	; t2->SetBranchAddress("crt_pmt_dist"    , _crt_pmt_dist      , &_b2_crt_pmt_dist       );





