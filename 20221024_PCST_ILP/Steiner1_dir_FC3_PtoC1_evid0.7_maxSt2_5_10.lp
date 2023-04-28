\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 0.021 x_Dbf4_Cdc5 + 0.027 x_Cdc5_Ipl1 + 0.115 x_Cdc5_Ndc80 + 0.241 x_Ndc80_Ipl1 + 0.076 x_Dbf4_Rad53 + 0.172 x_Rad53_Cdc5 + 0.26 x_Rad53_Hht1 + 0.273 x_Hht1_Ipl1 + 0.154 x_Rad9_Cdc28 + 0.202 x_Cdc28_Cdc5 + 0.208 x_Cdc28_Slt2 + 0.194 x_Slt2_Spa2 + 0.023 x_Rad9_Hta1 + 0.003 x_Hta1_Hht1 + 0.023 x_Rad9_Hta2 + 0.011 x_Hta2_Hht1 + 0.017 x_Rad9_Rad53 + 0.022 x_Rad9_Rfa1 + 0.259 x_Rfa1_Hht1 + 0.237 x_Rad9_Rps31 + 0.2 x_Rps31_Nop1 + 0.001 x_Nop1_Nsr1 + 0.185 x_Rps31_Nop56 + 0.001 x_Nop56_Nsr1 + 0.161 x_Rps31_Snu13 + 0.038 x_Snu13_Nsr1 + 0.001 x_Snu13_Snu66 + 0.246 x_Rpo21_Glc7 + 0.037 x_Glc7_Cef1 + 0.081 x_Cef1_Snu66 + 0.189 x_Glc7_Hsh155 + 0.001 x_Hsh155_Snu66 + 0.175 x_Glc7_Ipl1 + 0.234 x_Glc7_Smd1 + 0.016 x_Smd1_Snu66 + 0.232 x_Glc7_Smd2 + 0.001 x_Smd2_Snu66 + 0.214 x_Rpo21_Hog1 + 0.214 x_Hog1_Slt2 + 0.298 x_Rpo21_Npl3 + 0.022 x_Npl3_Cef1 + 0.193 x_Npl3_Mud1 + 0.001 x_Mud1_Snu66 + 0.133 x_Npl3_Nam8 + 0.001 x_Nam8_Snu66 + 0.185 x_Npl3_Snp1 + 0.028 x_Snp1_Snu66 + 0.183 x_Npl3_Sto1 + 0.001 x_Sto1_Snu66 + 0.00900000000000001 x_Rpo21_Sgv1 + 0.197 x_Sgv1_Sto1 + 0.23 x_Rpo21_Spt16 + 0.047 x_Spt16_Hht1 + 0.005 x_Rpo21_Ssu72 + 0.001 x_Ssu72_Glc7 + 0.014 x_Rpo21_Sth1 + 0.042 x_Sth1_Hht1 + 0.024 x_Rpo21_Taf14 + 0.212 x_Taf14_Hht1 + 0.001 x_Rpo21_Spt6 + 0.239 x_Spt6_Hht1 + 0.004 x_Rpo21_Set2 + 0.084 x_Set2_Hht1 + 0.082 x_Stb1_Cdc28 + 0.08 x_Ptp3_Hog1 + 0.196 x_Ptp3_Ptp2 + 0.26 x_Ptp2_Slt2 - 1.74838578814271 y_Dbf4 - 0.383515987743196 y_Cdc5 - 0 y_Ndc80 - 0 y_Rad53 - 0 y_Hht1 - 1.52061735460832 y_Rad9 - 0.121462293216658 y_Cdc28 - 0 y_Slt2 - 0 y_Hta1 - 0 y_Hta2 - 0.059240798414649 y_Rfa1 - 0.227967035078496 y_Rps31 - 0 y_Nop1 - 0 y_Nop56 - 0 y_Snu13 - 1.8162863958549 y_Rpo21 - 0 y_Glc7 - 0.070726960384864 y_Cef1 - 0 y_Hsh155 - 0 y_Smd1 - 0 y_Smd2 - 1.81347223589972 y_Hog1 - 0.063200435312774 y_Npl3 - 0 y_Mud1 - 0 y_Nam8 - 0 y_Snp1 - 0 y_Sto1 - 0.158874401057077 y_Sgv1 - 0.0822614010998221 y_Spt16 - 0 y_Ssu72 - 0 y_Sth1 - 0 y_Taf14 - 0.216056193029292 y_Spt6 - 0.0493386562549077 y_Set2 - 1.84450577831956 y_Stb1 - 2.02585031541065 y_Ptp3 - 0.173434899484957 y_Ptp2 - 2.22879578185143 y_Ipl1 - 4.20790455734154 y_Spa2 - 1.91624538979193 y_Nsr1 - 1.12613426133923 y_Snu66

subject to
x_Dbf4_Cdc5 + x_Dbf4_Rad53 >= 1 \edges from terminal
x_Rad9_Rad53 + x_Rad9_Cdc28 + x_Rad9_Hta1 + x_Rad9_Hta2 + x_Rad9_Rfa1 + x_Rad9_Rps31 >= 1 \edges from terminal
x_Rpo21_Glc7 + x_Rpo21_Hog1 + x_Rpo21_Npl3 + x_Rpo21_Sgv1 + x_Rpo21_Spt16 + x_Rpo21_Ssu72 + x_Rpo21_Sth1 + x_Rpo21_Taf14 + x_Rpo21_Spt6 + x_Rpo21_Set2 >= 1 \edges from terminal
x_Stb1_Cdc28 >= 1 \edges from terminal
x_Ptp3_Hog1 + x_Ptp3_Ptp2 >= 1 \edges from terminal
x_Cdc5_Ipl1 + x_Ndc80_Ipl1 + x_Hht1_Ipl1 + x_Glc7_Ipl1 >= 1 \edges to terminal
x_Slt2_Spa2 >= 1 \edges to terminal
x_Nop1_Nsr1 + x_Nop56_Nsr1 + x_Snu13_Nsr1 >= 1 \edges to terminal
x_Snu13_Snu66 + x_Cef1_Snu66 + x_Hsh155_Snu66 + x_Smd1_Snu66 + x_Smd2_Snu66 + x_Mud1_Snu66 + x_Nam8_Snu66 + x_Snp1_Snu66 + x_Sto1_Snu66 >= 1 \edges to terminal
y_Dbf4 >= 1  \terminal
y_Rad9 >= 1  \terminal
y_Rpo21 >= 1  \terminal
y_Stb1 >= 1  \terminal
y_Ptp3 >= 1  \terminal
y_Ipl1 >= 1  \terminal
y_Spa2 >= 1  \terminal
y_Nsr1 >= 1  \terminal
y_Snu66 >= 1  \terminal
y_Cdc5 - x_Cdc5_Ndc80 - x_Cdc5_Ipl1 <= 0 \Steiner node out
y_Ndc80 - x_Ndc80_Ipl1 <= 0 \Steiner node out
y_Rad53 - x_Rad53_Cdc5 - x_Rad53_Hht1 <= 0 \Steiner node out
y_Hht1 - x_Hht1_Ipl1 <= 0 \Steiner node out
y_Cdc28 - x_Cdc28_Cdc5 - x_Cdc28_Slt2 <= 0 \Steiner node out
y_Slt2 - x_Slt2_Spa2 <= 0 \Steiner node out
y_Hta1 - x_Hta1_Hht1 <= 0 \Steiner node out
y_Hta2 - x_Hta2_Hht1 <= 0 \Steiner node out
y_Rfa1 - x_Rfa1_Hht1 <= 0 \Steiner node out
y_Rps31 - x_Rps31_Nop1 - x_Rps31_Nop56 - x_Rps31_Snu13 <= 0 \Steiner node out
y_Nop1 - x_Nop1_Nsr1 <= 0 \Steiner node out
y_Nop56 - x_Nop56_Nsr1 <= 0 \Steiner node out
y_Snu13 - x_Snu13_Nsr1 - x_Snu13_Snu66 <= 0 \Steiner node out
y_Glc7 - x_Glc7_Cef1 - x_Glc7_Hsh155 - x_Glc7_Smd1 - x_Glc7_Smd2 - x_Glc7_Ipl1 <= 0 \Steiner node out
y_Cef1 - x_Cef1_Snu66 <= 0 \Steiner node out
y_Hsh155 - x_Hsh155_Snu66 <= 0 \Steiner node out
y_Smd1 - x_Smd1_Snu66 <= 0 \Steiner node out
y_Smd2 - x_Smd2_Snu66 <= 0 \Steiner node out
y_Hog1 - x_Hog1_Slt2 <= 0 \Steiner node out
y_Npl3 - x_Npl3_Cef1 - x_Npl3_Mud1 - x_Npl3_Nam8 - x_Npl3_Snp1 - x_Npl3_Sto1 <= 0 \Steiner node out
y_Mud1 - x_Mud1_Snu66 <= 0 \Steiner node out
y_Nam8 - x_Nam8_Snu66 <= 0 \Steiner node out
y_Snp1 - x_Snp1_Snu66 <= 0 \Steiner node out
y_Sto1 - x_Sto1_Snu66 <= 0 \Steiner node out
y_Sgv1 - x_Sgv1_Sto1 <= 0 \Steiner node out
y_Spt16 - x_Spt16_Hht1 <= 0 \Steiner node out
y_Ssu72 - x_Ssu72_Glc7 <= 0 \Steiner node out
y_Sth1 - x_Sth1_Hht1 <= 0 \Steiner node out
y_Taf14 - x_Taf14_Hht1 <= 0 \Steiner node out
y_Spt6 - x_Spt6_Hht1 <= 0 \Steiner node out
y_Set2 - x_Set2_Hht1 <= 0 \Steiner node out
y_Ptp2 - x_Ptp2_Slt2 <= 0 \Steiner node out
y_Cdc5 - x_Dbf4_Cdc5 - x_Rad53_Cdc5 - x_Cdc28_Cdc5 <= 0 \Steiner node in
y_Ndc80 - x_Cdc5_Ndc80 <= 0 \Steiner node in
y_Rad53 - x_Dbf4_Rad53 - x_Rad9_Rad53 <= 0 \Steiner node in
y_Hht1 - x_Rad53_Hht1 - x_Hta1_Hht1 - x_Hta2_Hht1 - x_Rfa1_Hht1 - x_Spt16_Hht1 - x_Sth1_Hht1 - x_Taf14_Hht1 - x_Spt6_Hht1 - x_Set2_Hht1 <= 0 \Steiner node in
y_Cdc28 - x_Rad9_Cdc28 - x_Stb1_Cdc28 <= 0 \Steiner node in
y_Slt2 - x_Cdc28_Slt2 - x_Hog1_Slt2 - x_Ptp2_Slt2 <= 0 \Steiner node in
y_Hta1 - x_Rad9_Hta1 <= 0 \Steiner node in
y_Hta2 - x_Rad9_Hta2 <= 0 \Steiner node in
y_Rfa1 - x_Rad9_Rfa1 <= 0 \Steiner node in
y_Rps31 - x_Rad9_Rps31 <= 0 \Steiner node in
y_Nop1 - x_Rps31_Nop1 <= 0 \Steiner node in
y_Nop56 - x_Rps31_Nop56 <= 0 \Steiner node in
y_Snu13 - x_Rps31_Snu13 <= 0 \Steiner node in
y_Glc7 - x_Rpo21_Glc7 - x_Ssu72_Glc7 <= 0 \Steiner node in
y_Cef1 - x_Glc7_Cef1 - x_Npl3_Cef1 <= 0 \Steiner node in
y_Hsh155 - x_Glc7_Hsh155 <= 0 \Steiner node in
y_Smd1 - x_Glc7_Smd1 <= 0 \Steiner node in
y_Smd2 - x_Glc7_Smd2 <= 0 \Steiner node in
y_Hog1 - x_Rpo21_Hog1 - x_Ptp3_Hog1 <= 0 \Steiner node in
y_Npl3 - x_Rpo21_Npl3 <= 0 \Steiner node in
y_Mud1 - x_Npl3_Mud1 <= 0 \Steiner node in
y_Nam8 - x_Npl3_Nam8 <= 0 \Steiner node in
y_Snp1 - x_Npl3_Snp1 <= 0 \Steiner node in
y_Sto1 - x_Npl3_Sto1 - x_Sgv1_Sto1 <= 0 \Steiner node in
y_Sgv1 - x_Rpo21_Sgv1 <= 0 \Steiner node in
y_Spt16 - x_Rpo21_Spt16 <= 0 \Steiner node in
y_Ssu72 - x_Rpo21_Ssu72 <= 0 \Steiner node in
y_Sth1 - x_Rpo21_Sth1 <= 0 \Steiner node in
y_Taf14 - x_Rpo21_Taf14 <= 0 \Steiner node in
y_Spt6 - x_Rpo21_Spt6 <= 0 \Steiner node in
y_Set2 - x_Rpo21_Set2 <= 0 \Steiner node in
y_Ptp2 - x_Ptp3_Ptp2 <= 0 \Steiner node in
y_Dbf4 - x_Dbf4_Cdc5 >= 0 \select node if connected (outgoing)
y_Dbf4 - x_Dbf4_Rad53 >= 0 \select node if connected (outgoing)
y_Cdc5 - x_Cdc5_Ndc80 >= 0 \select node if connected (outgoing)
y_Cdc5 - x_Cdc5_Ipl1 >= 0 \select node if connected (outgoing)
y_Ndc80 - x_Ndc80_Ipl1 >= 0 \select node if connected (outgoing)
y_Rad53 - x_Rad53_Cdc5 >= 0 \select node if connected (outgoing)
y_Rad53 - x_Rad53_Hht1 >= 0 \select node if connected (outgoing)
y_Hht1 - x_Hht1_Ipl1 >= 0 \select node if connected (outgoing)
y_Rad9 - x_Rad9_Rad53 >= 0 \select node if connected (outgoing)
y_Rad9 - x_Rad9_Cdc28 >= 0 \select node if connected (outgoing)
y_Rad9 - x_Rad9_Hta1 >= 0 \select node if connected (outgoing)
y_Rad9 - x_Rad9_Hta2 >= 0 \select node if connected (outgoing)
y_Rad9 - x_Rad9_Rfa1 >= 0 \select node if connected (outgoing)
y_Rad9 - x_Rad9_Rps31 >= 0 \select node if connected (outgoing)
y_Cdc28 - x_Cdc28_Cdc5 >= 0 \select node if connected (outgoing)
y_Cdc28 - x_Cdc28_Slt2 >= 0 \select node if connected (outgoing)
y_Slt2 - x_Slt2_Spa2 >= 0 \select node if connected (outgoing)
y_Hta1 - x_Hta1_Hht1 >= 0 \select node if connected (outgoing)
y_Hta2 - x_Hta2_Hht1 >= 0 \select node if connected (outgoing)
y_Rfa1 - x_Rfa1_Hht1 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Nop1 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Nop56 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Snu13 >= 0 \select node if connected (outgoing)
y_Nop1 - x_Nop1_Nsr1 >= 0 \select node if connected (outgoing)
y_Nop56 - x_Nop56_Nsr1 >= 0 \select node if connected (outgoing)
y_Snu13 - x_Snu13_Nsr1 >= 0 \select node if connected (outgoing)
y_Snu13 - x_Snu13_Snu66 >= 0 \select node if connected (outgoing)
y_Rpo21 - x_Rpo21_Glc7 >= 0 \select node if connected (outgoing)
y_Rpo21 - x_Rpo21_Hog1 >= 0 \select node if connected (outgoing)
y_Rpo21 - x_Rpo21_Npl3 >= 0 \select node if connected (outgoing)
y_Rpo21 - x_Rpo21_Sgv1 >= 0 \select node if connected (outgoing)
y_Rpo21 - x_Rpo21_Spt16 >= 0 \select node if connected (outgoing)
y_Rpo21 - x_Rpo21_Ssu72 >= 0 \select node if connected (outgoing)
y_Rpo21 - x_Rpo21_Sth1 >= 0 \select node if connected (outgoing)
y_Rpo21 - x_Rpo21_Taf14 >= 0 \select node if connected (outgoing)
y_Rpo21 - x_Rpo21_Spt6 >= 0 \select node if connected (outgoing)
y_Rpo21 - x_Rpo21_Set2 >= 0 \select node if connected (outgoing)
y_Glc7 - x_Glc7_Cef1 >= 0 \select node if connected (outgoing)
y_Glc7 - x_Glc7_Hsh155 >= 0 \select node if connected (outgoing)
y_Glc7 - x_Glc7_Smd1 >= 0 \select node if connected (outgoing)
y_Glc7 - x_Glc7_Smd2 >= 0 \select node if connected (outgoing)
y_Glc7 - x_Glc7_Ipl1 >= 0 \select node if connected (outgoing)
y_Cef1 - x_Cef1_Snu66 >= 0 \select node if connected (outgoing)
y_Hsh155 - x_Hsh155_Snu66 >= 0 \select node if connected (outgoing)
y_Smd1 - x_Smd1_Snu66 >= 0 \select node if connected (outgoing)
y_Smd2 - x_Smd2_Snu66 >= 0 \select node if connected (outgoing)
y_Hog1 - x_Hog1_Slt2 >= 0 \select node if connected (outgoing)
y_Npl3 - x_Npl3_Cef1 >= 0 \select node if connected (outgoing)
y_Npl3 - x_Npl3_Mud1 >= 0 \select node if connected (outgoing)
y_Npl3 - x_Npl3_Nam8 >= 0 \select node if connected (outgoing)
y_Npl3 - x_Npl3_Snp1 >= 0 \select node if connected (outgoing)
y_Npl3 - x_Npl3_Sto1 >= 0 \select node if connected (outgoing)
y_Mud1 - x_Mud1_Snu66 >= 0 \select node if connected (outgoing)
y_Nam8 - x_Nam8_Snu66 >= 0 \select node if connected (outgoing)
y_Snp1 - x_Snp1_Snu66 >= 0 \select node if connected (outgoing)
y_Sto1 - x_Sto1_Snu66 >= 0 \select node if connected (outgoing)
y_Sgv1 - x_Sgv1_Sto1 >= 0 \select node if connected (outgoing)
y_Spt16 - x_Spt16_Hht1 >= 0 \select node if connected (outgoing)
y_Ssu72 - x_Ssu72_Glc7 >= 0 \select node if connected (outgoing)
y_Sth1 - x_Sth1_Hht1 >= 0 \select node if connected (outgoing)
y_Taf14 - x_Taf14_Hht1 >= 0 \select node if connected (outgoing)
y_Spt6 - x_Spt6_Hht1 >= 0 \select node if connected (outgoing)
y_Set2 - x_Set2_Hht1 >= 0 \select node if connected (outgoing)
y_Stb1 - x_Stb1_Cdc28 >= 0 \select node if connected (outgoing)
y_Ptp3 - x_Ptp3_Hog1 >= 0 \select node if connected (outgoing)
y_Ptp3 - x_Ptp3_Ptp2 >= 0 \select node if connected (outgoing)
y_Ptp2 - x_Ptp2_Slt2 >= 0 \select node if connected (outgoing)
y_Cdc5 - x_Dbf4_Cdc5 >= 0 \select node if connected (incoming)
y_Cdc5 - x_Rad53_Cdc5 >= 0 \select node if connected (incoming)
y_Cdc5 - x_Cdc28_Cdc5 >= 0 \select node if connected (incoming)
y_Ndc80 - x_Cdc5_Ndc80 >= 0 \select node if connected (incoming)
y_Rad53 - x_Dbf4_Rad53 >= 0 \select node if connected (incoming)
y_Rad53 - x_Rad9_Rad53 >= 0 \select node if connected (incoming)
y_Hht1 - x_Rad53_Hht1 >= 0 \select node if connected (incoming)
y_Hht1 - x_Hta1_Hht1 >= 0 \select node if connected (incoming)
y_Hht1 - x_Hta2_Hht1 >= 0 \select node if connected (incoming)
y_Hht1 - x_Rfa1_Hht1 >= 0 \select node if connected (incoming)
y_Hht1 - x_Spt16_Hht1 >= 0 \select node if connected (incoming)
y_Hht1 - x_Sth1_Hht1 >= 0 \select node if connected (incoming)
y_Hht1 - x_Taf14_Hht1 >= 0 \select node if connected (incoming)
y_Hht1 - x_Spt6_Hht1 >= 0 \select node if connected (incoming)
y_Hht1 - x_Set2_Hht1 >= 0 \select node if connected (incoming)
y_Cdc28 - x_Rad9_Cdc28 >= 0 \select node if connected (incoming)
y_Cdc28 - x_Stb1_Cdc28 >= 0 \select node if connected (incoming)
y_Slt2 - x_Cdc28_Slt2 >= 0 \select node if connected (incoming)
y_Slt2 - x_Hog1_Slt2 >= 0 \select node if connected (incoming)
y_Slt2 - x_Ptp2_Slt2 >= 0 \select node if connected (incoming)
y_Hta1 - x_Rad9_Hta1 >= 0 \select node if connected (incoming)
y_Hta2 - x_Rad9_Hta2 >= 0 \select node if connected (incoming)
y_Rfa1 - x_Rad9_Rfa1 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rad9_Rps31 >= 0 \select node if connected (incoming)
y_Nop1 - x_Rps31_Nop1 >= 0 \select node if connected (incoming)
y_Nop56 - x_Rps31_Nop56 >= 0 \select node if connected (incoming)
y_Snu13 - x_Rps31_Snu13 >= 0 \select node if connected (incoming)
y_Glc7 - x_Rpo21_Glc7 >= 0 \select node if connected (incoming)
y_Glc7 - x_Ssu72_Glc7 >= 0 \select node if connected (incoming)
y_Cef1 - x_Glc7_Cef1 >= 0 \select node if connected (incoming)
y_Cef1 - x_Npl3_Cef1 >= 0 \select node if connected (incoming)
y_Hsh155 - x_Glc7_Hsh155 >= 0 \select node if connected (incoming)
y_Smd1 - x_Glc7_Smd1 >= 0 \select node if connected (incoming)
y_Smd2 - x_Glc7_Smd2 >= 0 \select node if connected (incoming)
y_Hog1 - x_Rpo21_Hog1 >= 0 \select node if connected (incoming)
y_Hog1 - x_Ptp3_Hog1 >= 0 \select node if connected (incoming)
y_Npl3 - x_Rpo21_Npl3 >= 0 \select node if connected (incoming)
y_Mud1 - x_Npl3_Mud1 >= 0 \select node if connected (incoming)
y_Nam8 - x_Npl3_Nam8 >= 0 \select node if connected (incoming)
y_Snp1 - x_Npl3_Snp1 >= 0 \select node if connected (incoming)
y_Sto1 - x_Npl3_Sto1 >= 0 \select node if connected (incoming)
y_Sto1 - x_Sgv1_Sto1 >= 0 \select node if connected (incoming)
y_Sgv1 - x_Rpo21_Sgv1 >= 0 \select node if connected (incoming)
y_Spt16 - x_Rpo21_Spt16 >= 0 \select node if connected (incoming)
y_Ssu72 - x_Rpo21_Ssu72 >= 0 \select node if connected (incoming)
y_Sth1 - x_Rpo21_Sth1 >= 0 \select node if connected (incoming)
y_Taf14 - x_Rpo21_Taf14 >= 0 \select node if connected (incoming)
y_Spt6 - x_Rpo21_Spt6 >= 0 \select node if connected (incoming)
y_Set2 - x_Rpo21_Set2 >= 0 \select node if connected (incoming)
y_Ptp2 - x_Ptp3_Ptp2 >= 0 \select node if connected (incoming)
y_Ipl1 - x_Cdc5_Ipl1 >= 0 \select node if connected (incoming)
y_Ipl1 - x_Ndc80_Ipl1 >= 0 \select node if connected (incoming)
y_Ipl1 - x_Hht1_Ipl1 >= 0 \select node if connected (incoming)
y_Ipl1 - x_Glc7_Ipl1 >= 0 \select node if connected (incoming)
y_Spa2 - x_Slt2_Spa2 >= 0 \select node if connected (incoming)
y_Nsr1 - x_Nop1_Nsr1 >= 0 \select node if connected (incoming)
y_Nsr1 - x_Nop56_Nsr1 >= 0 \select node if connected (incoming)
y_Nsr1 - x_Snu13_Nsr1 >= 0 \select node if connected (incoming)
y_Snu66 - x_Snu13_Snu66 >= 0 \select node if connected (incoming)
y_Snu66 - x_Cef1_Snu66 >= 0 \select node if connected (incoming)
y_Snu66 - x_Hsh155_Snu66 >= 0 \select node if connected (incoming)
y_Snu66 - x_Smd1_Snu66 >= 0 \select node if connected (incoming)
y_Snu66 - x_Smd2_Snu66 >= 0 \select node if connected (incoming)
y_Snu66 - x_Mud1_Snu66 >= 0 \select node if connected (incoming)
y_Snu66 - x_Nam8_Snu66 >= 0 \select node if connected (incoming)
y_Snu66 - x_Snp1_Snu66 >= 0 \select node if connected (incoming)
y_Snu66 - x_Sto1_Snu66 >= 0 \select node if connected (incoming)

binary
y_Dbf4
y_Cdc5
y_Ndc80
y_Rad53
y_Hht1
y_Rad9
y_Cdc28
y_Slt2
y_Hta1
y_Hta2
y_Rfa1
y_Rps31
y_Nop1
y_Nop56
y_Snu13
y_Rpo21
y_Glc7
y_Cef1
y_Hsh155
y_Smd1
y_Smd2
y_Hog1
y_Npl3
y_Mud1
y_Nam8
y_Snp1
y_Sto1
y_Sgv1
y_Spt16
y_Ssu72
y_Sth1
y_Taf14
y_Spt6
y_Set2
y_Stb1
y_Ptp3
y_Ptp2
y_Ipl1
y_Spa2
y_Nsr1
y_Snu66
x_Dbf4_Cdc5
x_Dbf4_Rad53
x_Cdc5_Ndc80
x_Cdc5_Ipl1
x_Ndc80_Ipl1
x_Rad53_Cdc5
x_Rad53_Hht1
x_Hht1_Ipl1
x_Rad9_Rad53
x_Rad9_Cdc28
x_Rad9_Hta1
x_Rad9_Hta2
x_Rad9_Rfa1
x_Rad9_Rps31
x_Cdc28_Cdc5
x_Cdc28_Slt2
x_Slt2_Spa2
x_Hta1_Hht1
x_Hta2_Hht1
x_Rfa1_Hht1
x_Rps31_Nop1
x_Rps31_Nop56
x_Rps31_Snu13
x_Nop1_Nsr1
x_Nop56_Nsr1
x_Snu13_Nsr1
x_Snu13_Snu66
x_Rpo21_Glc7
x_Rpo21_Hog1
x_Rpo21_Npl3
x_Rpo21_Sgv1
x_Rpo21_Spt16
x_Rpo21_Ssu72
x_Rpo21_Sth1
x_Rpo21_Taf14
x_Rpo21_Spt6
x_Rpo21_Set2
x_Glc7_Cef1
x_Glc7_Hsh155
x_Glc7_Smd1
x_Glc7_Smd2
x_Glc7_Ipl1
x_Cef1_Snu66
x_Hsh155_Snu66
x_Smd1_Snu66
x_Smd2_Snu66
x_Hog1_Slt2
x_Npl3_Cef1
x_Npl3_Mud1
x_Npl3_Nam8
x_Npl3_Snp1
x_Npl3_Sto1
x_Mud1_Snu66
x_Nam8_Snu66
x_Snp1_Snu66
x_Sto1_Snu66
x_Sgv1_Sto1
x_Spt16_Hht1
x_Ssu72_Glc7
x_Sth1_Hht1
x_Taf14_Hht1
x_Spt6_Hht1
x_Set2_Hht1
x_Stb1_Cdc28
x_Ptp3_Hog1
x_Ptp3_Ptp2
x_Ptp2_Slt2

End