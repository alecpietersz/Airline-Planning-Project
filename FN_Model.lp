\ Model FN
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 45.42746255909827 wLIRA-LIMF - 35.36613299523423 wLIRA-LIBD
   - 45.42746255909827 wLIMF-LIRA - 60.38950860166594 wLIMF-LIBD
   - 35.36613299523423 wLIBD-LIRA - 60.38950860166594 wLIBD-LIMF
   - 50.47495839899808 xLIRA-LIMF - 39.29570332803803 xLIRA-LIBD
   - 50.47495839899808 xLIMF-LIRA - 67.09945400185104 xLIMF-LIBD
   - 39.29570332803803 xLIBD-LIRA - 67.09945400185104 xLIBD-LIMF
   + 5400 zLIRA-LIRA-Regional_Turboprop
   + 49476 zLIRA-LIMF-Regional_Turboprop
   + 35139.41454545454 zLIRA-LIBD-Regional_Turboprop
   + 49476 zLIMF-LIRA-Regional_Turboprop
   + 13500 zLIMF-LIMF-Regional_Turboprop
   + 103428.5454545455 zLIMF-LIBD-Regional_Turboprop
   + 35139.41454545454 zLIBD-LIRA-Regional_Turboprop
   + 103428.5454545455 zLIBD-LIMF-Regional_Turboprop
   + 13500 zLIBD-LIBD-Regional_Turboprop + 15000 ACRegional_Turboprop
Subject To
 qnorLIRA-LIRA: wLIRA-LIRA + xLIRA-LIRA <= 0
 qnorLIRA-LIMF: wLIRA-LIMF + xLIRA-LIMF <= 0
 qnorLIRA-LIBD: wLIRA-LIBD + xLIRA-LIBD <= 0
 qnorLIMF-LIRA: wLIMF-LIRA + xLIMF-LIRA <= 0
 qnorLIMF-LIMF: wLIMF-LIMF + xLIMF-LIMF <= 0
 qnorLIMF-LIBD: wLIMF-LIBD + xLIMF-LIBD <= 0
 qnorLIBD-LIRA: wLIBD-LIRA + xLIBD-LIRA <= 0
 qnorLIBD-LIMF: wLIBD-LIMF + xLIBD-LIMF <= 0
 qnorLIBD-LIBD: wLIBD-LIBD + xLIBD-LIBD <= 0
 qtransLIRA-LIRA: wLIRA-LIRA <= 0
 qtransLIRA-LIMF: wLIRA-LIMF <= 0
 qtransLIRA-LIBD: wLIRA-LIBD <= 0
 qtransLIMF-LIRA: wLIMF-LIRA <= 0
 qtransLIMF-LIMF: wLIMF-LIMF <= 0
 qtransLIMF-LIBD: wLIMF-LIBD <= 0
 qtransLIBD-LIRA: wLIBD-LIRA <= 0
 qtransLIBD-LIMF: wLIBD-LIMF <= 0
 qtransLIBD-LIBD: wLIBD-LIBD <= 0
 CapacityLIRA-LIRA: 2 wLIRA-LIRA + wLIRA-LIMF + wLIRA-LIBD + wLIMF-LIRA
   + wLIBD-LIRA + xLIRA-LIRA - 36 zLIRA-LIRA-Regional_Turboprop <= 0
 CapacityLIRA-LIMF: wLIRA-LIMF + wLIMF-LIMF + wLIBD-LIMF + xLIRA-LIMF
   - 36 zLIRA-LIMF-Regional_Turboprop <= 0
 CapacityLIRA-LIBD: wLIRA-LIBD + wLIMF-LIBD + wLIBD-LIBD + xLIRA-LIBD
   - 36 zLIRA-LIBD-Regional_Turboprop <= 0
 CapacityLIMF-LIRA: wLIMF-LIRA + wLIMF-LIMF + wLIMF-LIBD + xLIMF-LIRA
   - 36 zLIMF-LIRA-Regional_Turboprop <= 0
 CapacityLIMF-LIMF: xLIMF-LIMF - 36 zLIMF-LIMF-Regional_Turboprop <= 0
 CapacityLIMF-LIBD: xLIMF-LIBD - 36 zLIMF-LIBD-Regional_Turboprop <= 0
 CapacityLIBD-LIRA: wLIBD-LIRA + wLIBD-LIMF + wLIBD-LIBD + xLIBD-LIRA
   - 36 zLIBD-LIRA-Regional_Turboprop <= 0
 CapacityLIBD-LIMF: xLIBD-LIMF - 36 zLIBD-LIMF-Regional_Turboprop <= 0
 CapacityLIBD-LIBD: xLIBD-LIBD - 36 zLIBD-LIBD-Regional_Turboprop <= 0
 ContinuityLIRARegional_Turboprop: zLIRA-LIMF-Regional_Turboprop
   + zLIRA-LIBD-Regional_Turboprop - zLIMF-LIRA-Regional_Turboprop
   - zLIBD-LIRA-Regional_Turboprop = 0
 ContinuityLIMFRegional_Turboprop: - zLIRA-LIMF-Regional_Turboprop
   + zLIMF-LIRA-Regional_Turboprop + zLIMF-LIBD-Regional_Turboprop
   - zLIBD-LIMF-Regional_Turboprop = 0
 ContinuityLIBDRegional_Turboprop: - zLIRA-LIBD-Regional_Turboprop
   - zLIMF-LIBD-Regional_Turboprop + zLIBD-LIRA-Regional_Turboprop
   + zLIBD-LIMF-Regional_Turboprop = 0
 ACProductivityRegional_Turboprop:
   0.4166666666666667 zLIRA-LIRA-Regional_Turboprop
   + 1.416666666666667 zLIRA-LIMF-Regional_Turboprop
   + 1.058484848484849 zLIRA-LIBD-Regional_Turboprop
   + 1.416666666666667 zLIMF-LIRA-Regional_Turboprop
   + 0.4166666666666667 zLIMF-LIMF-Regional_Turboprop
   + 1.989393939393939 zLIMF-LIBD-Regional_Turboprop
   + 1.058484848484849 zLIBD-LIRA-Regional_Turboprop
   + 1.989393939393939 zLIBD-LIMF-Regional_Turboprop
   + 0.4166666666666667 zLIBD-LIBD-Regional_Turboprop
   - 10 ACRegional_Turboprop <= 0
 runwayLIRA-LIRARegional_Turboprop: zLIRA-LIRA-Regional_Turboprop <= 10000
 runwayLIRA-LIMFRegional_Turboprop: zLIRA-LIMF-Regional_Turboprop <= 10000
 runwayLIRA-LIBDRegional_Turboprop: zLIRA-LIBD-Regional_Turboprop <= 10000
 runwayLIMF-LIRARegional_Turboprop: zLIMF-LIRA-Regional_Turboprop <= 10000
 runwayLIMF-LIMFRegional_Turboprop: zLIMF-LIMF-Regional_Turboprop <= 10000
 runwayLIMF-LIBDRegional_Turboprop: zLIMF-LIBD-Regional_Turboprop <= 10000
 runwayLIBD-LIRARegional_Turboprop: zLIBD-LIRA-Regional_Turboprop <= 10000
 runwayLIBD-LIMFRegional_Turboprop: zLIBD-LIMF-Regional_Turboprop <= 10000
 runwayLIBD-LIBDRegional_Turboprop: zLIBD-LIBD-Regional_Turboprop <= 10000
 runwayLIRA-LIRARegional_Turboprop: zLIRA-LIRA-Regional_Turboprop <= 10000
 runwayLIRA-LIMFRegional_Turboprop: zLIRA-LIMF-Regional_Turboprop <= 10000
 runwayLIRA-LIBDRegional_Turboprop: zLIRA-LIBD-Regional_Turboprop <= 10000
 runwayLIMF-LIRARegional_Turboprop: zLIMF-LIRA-Regional_Turboprop <= 10000
 runwayLIMF-LIMFRegional_Turboprop: zLIMF-LIMF-Regional_Turboprop <= 10000
 runwayLIMF-LIBDRegional_Turboprop: zLIMF-LIBD-Regional_Turboprop <= 10000
 runwayLIBD-LIRARegional_Turboprop: zLIBD-LIRA-Regional_Turboprop <= 10000
 runwayLIBD-LIMFRegional_Turboprop: zLIBD-LIMF-Regional_Turboprop <= 10000
 runwayLIBD-LIBDRegional_Turboprop: zLIBD-LIBD-Regional_Turboprop <= 10000
Bounds
Generals
 wLIRA-LIRA wLIRA-LIMF wLIRA-LIBD wLIMF-LIRA wLIMF-LIMF wLIMF-LIBD
 wLIBD-LIRA wLIBD-LIMF wLIBD-LIBD xLIRA-LIRA xLIRA-LIMF xLIRA-LIBD
 xLIMF-LIRA xLIMF-LIMF xLIMF-LIBD xLIBD-LIRA xLIBD-LIMF xLIBD-LIBD
 zLIRA-LIRA-Regional_Turboprop zLIRA-LIMF-Regional_Turboprop
 zLIRA-LIBD-Regional_Turboprop zLIMF-LIRA-Regional_Turboprop
 zLIMF-LIMF-Regional_Turboprop zLIMF-LIBD-Regional_Turboprop
 zLIBD-LIRA-Regional_Turboprop zLIBD-LIMF-Regional_Turboprop
 zLIBD-LIBD-Regional_Turboprop ACRegional_Turboprop
End
