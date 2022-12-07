\ Model FN
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 45.42860250817204 wLIRA-LIMF - 35.41533474797972 wLIRA-LIBD
   - 45.42860250817204 wLIMF-LIRA - 60.41186043801749 wLIMF-LIBD
   - 35.41533474797972 wLIBD-LIRA - 60.41186043801749 wLIBD-LIMF
   - 50.47622500908005 xLIRA-LIMF - 39.35037194219969 xLIRA-LIBD
   - 50.47622500908005 xLIMF-LIRA - 67.12428937557499 xLIMF-LIBD
   - 39.35037194219969 xLIBD-LIRA - 67.12428937557499 xLIBD-LIMF
   + 120 zLIRA-LIRA-Regional_Turboprop
   + 1099.504110553851 zLIRA-LIMF-Regional_Turboprop
   + 782.3647664559321 zLIRA-LIBD-Regional_Turboprop
   + 1099.504110553851 zLIMF-LIRA-Regional_Turboprop
   + 300 zLIMF-LIMF-Regional_Turboprop
   + 2299.530685660759 zLIMF-LIBD-Regional_Turboprop
   + 782.3647664559321 zLIBD-LIRA-Regional_Turboprop
   + 2299.530685660759 zLIBD-LIMF-Regional_Turboprop
   + 300 zLIBD-LIBD-Regional_Turboprop + 1500 ACRegional_Turboprop
Subject To
 qnorLIRA-LIRA: wLIRA-LIRA + xLIRA-LIRA <= 0
 qnorLIRA-LIMF: wLIRA-LIMF + xLIRA-LIMF <= 110.8216688451043
 qnorLIRA-LIBD: wLIRA-LIBD + xLIRA-LIBD <= 215.6576325987535
 qnorLIMF-LIRA: wLIMF-LIRA + xLIMF-LIRA <= 110.8216688451043
 qnorLIMF-LIMF: wLIMF-LIMF + xLIMF-LIMF <= 0
 qnorLIMF-LIBD: wLIMF-LIBD + xLIMF-LIBD <= 66.70536156686033
 qnorLIBD-LIRA: wLIBD-LIRA + xLIBD-LIRA <= 215.6576325987535
 qnorLIBD-LIMF: wLIBD-LIMF + xLIBD-LIMF <= 66.70536156686033
 qnorLIBD-LIBD: wLIBD-LIBD + xLIBD-LIBD <= 0
 qtransLIRA-LIRA: wLIRA-LIRA <= 0
 qtransLIRA-LIMF: wLIRA-LIMF <= 0
 qtransLIRA-LIBD: wLIRA-LIBD <= 0
 qtransLIMF-LIRA: wLIMF-LIRA <= 0
 qtransLIMF-LIMF: wLIMF-LIMF <= 0
 qtransLIMF-LIBD: wLIMF-LIBD <= 66.70536156686033
 qtransLIBD-LIRA: wLIBD-LIRA <= 0
 qtransLIBD-LIMF: wLIBD-LIMF <= 66.70536156686033
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
 ACProductivityRegional_Turboprop: 0.625 zLIRA-LIRA-Regional_Turboprop
   + 1.416708763676692 zLIRA-LIMF-Regional_Turboprop
   + 1.060158759069529 zLIRA-LIBD-Regional_Turboprop
   + 1.625042097010025 zLIMF-LIRA-Regional_Turboprop
   + 0.4166666666666667 zLIMF-LIMF-Regional_Turboprop
   + 1.99027423670399 zLIMF-LIBD-Regional_Turboprop
   + 1.268492092402862 zLIBD-LIRA-Regional_Turboprop
   + 1.99027423670399 zLIBD-LIMF-Regional_Turboprop
   + 0.4166666666666667 zLIBD-LIBD-Regional_Turboprop
   - 70 ACRegional_Turboprop <= 0
 runwayLIRA-LIRARegional_Turboprop: zLIRA-LIRA-Regional_Turboprop <= 10000
 runwayLIRA-LIMFRegional_Turboprop: zLIRA-LIMF-Regional_Turboprop <= 10000
 runwayLIRA-LIBDRegional_Turboprop: zLIRA-LIBD-Regional_Turboprop <= 10000
 runwayLIMF-LIRARegional_Turboprop: zLIMF-LIRA-Regional_Turboprop <= 10000
 runwayLIMF-LIMFRegional_Turboprop: zLIMF-LIMF-Regional_Turboprop <= 10000
 runwayLIMF-LIBDRegional_Turboprop: zLIMF-LIBD-Regional_Turboprop <= 10000
 runwayLIBD-LIRARegional_Turboprop: zLIBD-LIRA-Regional_Turboprop <= 10000
 runwayLIBD-LIMFRegional_Turboprop: zLIBD-LIMF-Regional_Turboprop <= 10000
 runwayLIBD-LIBDRegional_Turboprop: zLIBD-LIBD-Regional_Turboprop <= 10000
 rangeLIRA-LIRARegional_Turboprop: zLIRA-LIRA-Regional_Turboprop <= 10000
 rangeLIRA-LIMFRegional_Turboprop: zLIRA-LIMF-Regional_Turboprop <= 10000
 rangeLIRA-LIBDRegional_Turboprop: zLIRA-LIBD-Regional_Turboprop <= 10000
 rangeLIMF-LIRARegional_Turboprop: zLIMF-LIRA-Regional_Turboprop <= 10000
 rangeLIMF-LIMFRegional_Turboprop: zLIMF-LIMF-Regional_Turboprop <= 10000
 rangeLIMF-LIBDRegional_Turboprop: zLIMF-LIBD-Regional_Turboprop <= 10000
 rangeLIBD-LIRARegional_Turboprop: zLIBD-LIRA-Regional_Turboprop <= 10000
 rangeLIBD-LIMFRegional_Turboprop: zLIBD-LIMF-Regional_Turboprop <= 10000
 rangeLIBD-LIBDRegional_Turboprop: zLIBD-LIBD-Regional_Turboprop <= 10000
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
