\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 0.224 x_Bem2_Bud14 + 0.232 x_Hog1_Kss1 + 0.152 x_Kss1_Fus3 + 0.221 x_Fus3_Pbs2 + 0.212 x_Kss1_Pbs2 + 0.175 x_Kss1_Ste11 + 0.128 x_Ste11_Pbs2 + 0.035 x_Hog1_Mih1 + 0.023 x_Mih1_Pbs2 + 0.001 x_Hog1_Pbs2 + 0.186 x_Hog1_Ptc1 + 0.024 x_Ptc1_Nbp2 + 0.205 x_Nbp2_Pbs2 + 0.203 x_Ptc1_Pbs2 + 0.095 x_Hog1_Ptp2 + 0.239 x_Ptp2_Ptc1 + 0.214 x_Hog1_Rpo21 + 0.246 x_Rpo21_Glc7 + 0.005 x_Glc7_Bud14 + 0.298 x_Rpo21_Npl3 + 0.00600000000000001 x_Npl3_Sky1 + 0.249 x_Hog1_Sho1 + 0.083 x_Sho1_Pbs2 + 0.269 x_Sho1_Ste11 + 0.001 x_Hog1_Sko1 + 0.235 x_Sko1_Fus3 + 0.268 x_Sko1_Kss1 + 0.214 x_Hog1_Slt2 + 0.183 x_Slt2_Ptc1 + 0.1 x_Hog1_Ste7 + 0.001 x_Ste7_Fus3 + 0.004 x_Ste7_Kss1 + 0.039 x_Ste7_Ras2 + 0.178 x_Ras2_Gpa2 + 0.005 x_Ste7_Ste11 + 0.029 x_Hog1_Rck2 + 0.215 x_Rim15_Pho85 + 0.039 x_Pho85_Mih1 + 0.182 x_Rim15_Sch9 + 0.174 x_Sch9_Ras2 + 0.268 x_Rim15_Tpk1 + 0.115 x_Tpk1_Ras2 + 0.00700000000000001 x_Tif2_Tif4631 + 0.259 x_Tif4631_Npl3 + 0.266 x_Tif4631_Yra1 + 0.004 x_Yra1_Sub2 + 0.124 x_Tif2_Tif4632 + 0.266 x_Tif4632_Yra1 - 3.6868019743945 y_Bem2 - 1.81347223589972 y_Hog1 - 0 y_Kss1 - 0 y_Fus3 - 0.10305136110519 y_Ste11 - 0.116628340365083 y_Mih1 - 0.0854152687766444 y_Ptc1 - 0.663781033634114 y_Nbp2 - 0.173434899484957 y_Ptp2 - 1.8162863958549 y_Rpo21 - 0 y_Glc7 - 0.063200435312774 y_Npl3 - 0 y_Sho1 - 1.50119873871931 y_Sko1 - 0 y_Slt2 - 0 y_Ste7 - 0.080694252350959 y_Ras2 - 1.00036416841318 y_Rim15 - 0 y_Pho85 - 1.51300508365901 y_Sch9 - 0 y_Tpk1 - 2.893087867955 y_Tif2 - 0.367065514001302 y_Tif4631 - 0.110683243874969 y_Yra1 - 0.517617071461765 y_Tif4632 - 1.87884256040446 y_Bud14 - 1.7465565764114 y_Pbs2 - 2.43145961134509 y_Sky1 - 1.01894679922407 y_Gpa2 - 3.58764096905305 y_Rck2 - 1.32099603495184 y_Sub2

subject to
x_Bem2_Bud14 >= 1 \edges from terminal
x_Hog1_Kss1 + x_Hog1_Mih1 + x_Hog1_Ptc1 + x_Hog1_Ptp2 + x_Hog1_Rpo21 + x_Hog1_Sho1 + x_Hog1_Sko1 + x_Hog1_Slt2 + x_Hog1_Ste7 + x_Hog1_Pbs2 + x_Hog1_Rck2 >= 1 \edges from terminal
x_Rim15_Pho85 + x_Rim15_Sch9 + x_Rim15_Tpk1 >= 1 \edges from terminal
x_Tif2_Tif4631 + x_Tif2_Tif4632 >= 1 \edges from terminal
x_Bem2_Bud14 + x_Glc7_Bud14 >= 1 \edges to terminal
x_Hog1_Pbs2 + x_Kss1_Pbs2 + x_Fus3_Pbs2 + x_Ste11_Pbs2 + x_Mih1_Pbs2 + x_Ptc1_Pbs2 + x_Nbp2_Pbs2 + x_Sho1_Pbs2 >= 1 \edges to terminal
x_Npl3_Sky1 >= 1 \edges to terminal
x_Ras2_Gpa2 >= 1 \edges to terminal
x_Hog1_Rck2 >= 1 \edges to terminal
x_Yra1_Sub2 >= 1 \edges to terminal
y_Bem2 >= 1  \terminal
y_Hog1 >= 1  \terminal
y_Rim15 >= 1  \terminal
y_Tif2 >= 1  \terminal
y_Bud14 >= 1  \terminal
y_Pbs2 >= 1  \terminal
y_Sky1 >= 1  \terminal
y_Gpa2 >= 1  \terminal
y_Rck2 >= 1  \terminal
y_Sub2 >= 1  \terminal
y_Kss1 - x_Kss1_Fus3 - x_Kss1_Ste11 - x_Kss1_Pbs2 <= 0 \Steiner node out
y_Fus3 - x_Fus3_Pbs2 <= 0 \Steiner node out
y_Ste11 - x_Ste11_Pbs2 <= 0 \Steiner node out
y_Mih1 - x_Mih1_Pbs2 <= 0 \Steiner node out
y_Ptc1 - x_Ptc1_Nbp2 - x_Ptc1_Pbs2 <= 0 \Steiner node out
y_Nbp2 - x_Nbp2_Pbs2 <= 0 \Steiner node out
y_Ptp2 - x_Ptp2_Ptc1 <= 0 \Steiner node out
y_Rpo21 - x_Rpo21_Glc7 - x_Rpo21_Npl3 <= 0 \Steiner node out
y_Glc7 - x_Glc7_Bud14 <= 0 \Steiner node out
y_Npl3 - x_Npl3_Sky1 <= 0 \Steiner node out
y_Sho1 - x_Sho1_Ste11 - x_Sho1_Pbs2 <= 0 \Steiner node out
y_Sko1 - x_Sko1_Kss1 - x_Sko1_Fus3 <= 0 \Steiner node out
y_Slt2 - x_Slt2_Ptc1 <= 0 \Steiner node out
y_Ste7 - x_Ste7_Kss1 - x_Ste7_Fus3 - x_Ste7_Ste11 - x_Ste7_Ras2 <= 0 \Steiner node out
y_Ras2 - x_Ras2_Gpa2 <= 0 \Steiner node out
y_Pho85 - x_Pho85_Mih1 <= 0 \Steiner node out
y_Sch9 - x_Sch9_Ras2 <= 0 \Steiner node out
y_Tpk1 - x_Tpk1_Ras2 <= 0 \Steiner node out
y_Tif4631 - x_Tif4631_Npl3 - x_Tif4631_Yra1 <= 0 \Steiner node out
y_Yra1 - x_Yra1_Sub2 <= 0 \Steiner node out
y_Tif4632 - x_Tif4632_Yra1 <= 0 \Steiner node out
y_Kss1 - x_Hog1_Kss1 - x_Sko1_Kss1 - x_Ste7_Kss1 <= 0 \Steiner node in
y_Fus3 - x_Kss1_Fus3 - x_Sko1_Fus3 - x_Ste7_Fus3 <= 0 \Steiner node in
y_Ste11 - x_Kss1_Ste11 - x_Sho1_Ste11 - x_Ste7_Ste11 <= 0 \Steiner node in
y_Mih1 - x_Hog1_Mih1 - x_Pho85_Mih1 <= 0 \Steiner node in
y_Ptc1 - x_Hog1_Ptc1 - x_Ptp2_Ptc1 - x_Slt2_Ptc1 <= 0 \Steiner node in
y_Nbp2 - x_Ptc1_Nbp2 <= 0 \Steiner node in
y_Ptp2 - x_Hog1_Ptp2 <= 0 \Steiner node in
y_Rpo21 - x_Hog1_Rpo21 <= 0 \Steiner node in
y_Glc7 - x_Rpo21_Glc7 <= 0 \Steiner node in
y_Npl3 - x_Rpo21_Npl3 - x_Tif4631_Npl3 <= 0 \Steiner node in
y_Sho1 - x_Hog1_Sho1 <= 0 \Steiner node in
y_Sko1 - x_Hog1_Sko1 <= 0 \Steiner node in
y_Slt2 - x_Hog1_Slt2 <= 0 \Steiner node in
y_Ste7 - x_Hog1_Ste7 <= 0 \Steiner node in
y_Ras2 - x_Ste7_Ras2 - x_Sch9_Ras2 - x_Tpk1_Ras2 <= 0 \Steiner node in
y_Pho85 - x_Rim15_Pho85 <= 0 \Steiner node in
y_Sch9 - x_Rim15_Sch9 <= 0 \Steiner node in
y_Tpk1 - x_Rim15_Tpk1 <= 0 \Steiner node in
y_Tif4631 - x_Tif2_Tif4631 <= 0 \Steiner node in
y_Yra1 - x_Tif4631_Yra1 - x_Tif4632_Yra1 <= 0 \Steiner node in
y_Tif4632 - x_Tif2_Tif4632 <= 0 \Steiner node in
y_Bem2 - x_Bem2_Bud14 >= 0 \select node if connected (outgoing)
y_Hog1 - x_Hog1_Kss1 >= 0 \select node if connected (outgoing)
y_Hog1 - x_Hog1_Mih1 >= 0 \select node if connected (outgoing)
y_Hog1 - x_Hog1_Ptc1 >= 0 \select node if connected (outgoing)
y_Hog1 - x_Hog1_Ptp2 >= 0 \select node if connected (outgoing)
y_Hog1 - x_Hog1_Rpo21 >= 0 \select node if connected (outgoing)
y_Hog1 - x_Hog1_Sho1 >= 0 \select node if connected (outgoing)
y_Hog1 - x_Hog1_Sko1 >= 0 \select node if connected (outgoing)
y_Hog1 - x_Hog1_Slt2 >= 0 \select node if connected (outgoing)
y_Hog1 - x_Hog1_Ste7 >= 0 \select node if connected (outgoing)
y_Hog1 - x_Hog1_Pbs2 >= 0 \select node if connected (outgoing)
y_Hog1 - x_Hog1_Rck2 >= 0 \select node if connected (outgoing)
y_Kss1 - x_Kss1_Fus3 >= 0 \select node if connected (outgoing)
y_Kss1 - x_Kss1_Ste11 >= 0 \select node if connected (outgoing)
y_Kss1 - x_Kss1_Pbs2 >= 0 \select node if connected (outgoing)
y_Fus3 - x_Fus3_Pbs2 >= 0 \select node if connected (outgoing)
y_Ste11 - x_Ste11_Pbs2 >= 0 \select node if connected (outgoing)
y_Mih1 - x_Mih1_Pbs2 >= 0 \select node if connected (outgoing)
y_Ptc1 - x_Ptc1_Nbp2 >= 0 \select node if connected (outgoing)
y_Ptc1 - x_Ptc1_Pbs2 >= 0 \select node if connected (outgoing)
y_Nbp2 - x_Nbp2_Pbs2 >= 0 \select node if connected (outgoing)
y_Ptp2 - x_Ptp2_Ptc1 >= 0 \select node if connected (outgoing)
y_Rpo21 - x_Rpo21_Glc7 >= 0 \select node if connected (outgoing)
y_Rpo21 - x_Rpo21_Npl3 >= 0 \select node if connected (outgoing)
y_Glc7 - x_Glc7_Bud14 >= 0 \select node if connected (outgoing)
y_Npl3 - x_Npl3_Sky1 >= 0 \select node if connected (outgoing)
y_Sho1 - x_Sho1_Ste11 >= 0 \select node if connected (outgoing)
y_Sho1 - x_Sho1_Pbs2 >= 0 \select node if connected (outgoing)
y_Sko1 - x_Sko1_Kss1 >= 0 \select node if connected (outgoing)
y_Sko1 - x_Sko1_Fus3 >= 0 \select node if connected (outgoing)
y_Slt2 - x_Slt2_Ptc1 >= 0 \select node if connected (outgoing)
y_Ste7 - x_Ste7_Kss1 >= 0 \select node if connected (outgoing)
y_Ste7 - x_Ste7_Fus3 >= 0 \select node if connected (outgoing)
y_Ste7 - x_Ste7_Ste11 >= 0 \select node if connected (outgoing)
y_Ste7 - x_Ste7_Ras2 >= 0 \select node if connected (outgoing)
y_Ras2 - x_Ras2_Gpa2 >= 0 \select node if connected (outgoing)
y_Rim15 - x_Rim15_Pho85 >= 0 \select node if connected (outgoing)
y_Rim15 - x_Rim15_Sch9 >= 0 \select node if connected (outgoing)
y_Rim15 - x_Rim15_Tpk1 >= 0 \select node if connected (outgoing)
y_Pho85 - x_Pho85_Mih1 >= 0 \select node if connected (outgoing)
y_Sch9 - x_Sch9_Ras2 >= 0 \select node if connected (outgoing)
y_Tpk1 - x_Tpk1_Ras2 >= 0 \select node if connected (outgoing)
y_Tif2 - x_Tif2_Tif4631 >= 0 \select node if connected (outgoing)
y_Tif2 - x_Tif2_Tif4632 >= 0 \select node if connected (outgoing)
y_Tif4631 - x_Tif4631_Npl3 >= 0 \select node if connected (outgoing)
y_Tif4631 - x_Tif4631_Yra1 >= 0 \select node if connected (outgoing)
y_Yra1 - x_Yra1_Sub2 >= 0 \select node if connected (outgoing)
y_Tif4632 - x_Tif4632_Yra1 >= 0 \select node if connected (outgoing)
y_Kss1 - x_Hog1_Kss1 >= 0 \select node if connected (incoming)
y_Kss1 - x_Sko1_Kss1 >= 0 \select node if connected (incoming)
y_Kss1 - x_Ste7_Kss1 >= 0 \select node if connected (incoming)
y_Fus3 - x_Kss1_Fus3 >= 0 \select node if connected (incoming)
y_Fus3 - x_Sko1_Fus3 >= 0 \select node if connected (incoming)
y_Fus3 - x_Ste7_Fus3 >= 0 \select node if connected (incoming)
y_Ste11 - x_Kss1_Ste11 >= 0 \select node if connected (incoming)
y_Ste11 - x_Sho1_Ste11 >= 0 \select node if connected (incoming)
y_Ste11 - x_Ste7_Ste11 >= 0 \select node if connected (incoming)
y_Mih1 - x_Hog1_Mih1 >= 0 \select node if connected (incoming)
y_Mih1 - x_Pho85_Mih1 >= 0 \select node if connected (incoming)
y_Ptc1 - x_Hog1_Ptc1 >= 0 \select node if connected (incoming)
y_Ptc1 - x_Ptp2_Ptc1 >= 0 \select node if connected (incoming)
y_Ptc1 - x_Slt2_Ptc1 >= 0 \select node if connected (incoming)
y_Nbp2 - x_Ptc1_Nbp2 >= 0 \select node if connected (incoming)
y_Ptp2 - x_Hog1_Ptp2 >= 0 \select node if connected (incoming)
y_Rpo21 - x_Hog1_Rpo21 >= 0 \select node if connected (incoming)
y_Glc7 - x_Rpo21_Glc7 >= 0 \select node if connected (incoming)
y_Npl3 - x_Rpo21_Npl3 >= 0 \select node if connected (incoming)
y_Npl3 - x_Tif4631_Npl3 >= 0 \select node if connected (incoming)
y_Sho1 - x_Hog1_Sho1 >= 0 \select node if connected (incoming)
y_Sko1 - x_Hog1_Sko1 >= 0 \select node if connected (incoming)
y_Slt2 - x_Hog1_Slt2 >= 0 \select node if connected (incoming)
y_Ste7 - x_Hog1_Ste7 >= 0 \select node if connected (incoming)
y_Ras2 - x_Ste7_Ras2 >= 0 \select node if connected (incoming)
y_Ras2 - x_Sch9_Ras2 >= 0 \select node if connected (incoming)
y_Ras2 - x_Tpk1_Ras2 >= 0 \select node if connected (incoming)
y_Pho85 - x_Rim15_Pho85 >= 0 \select node if connected (incoming)
y_Sch9 - x_Rim15_Sch9 >= 0 \select node if connected (incoming)
y_Tpk1 - x_Rim15_Tpk1 >= 0 \select node if connected (incoming)
y_Tif4631 - x_Tif2_Tif4631 >= 0 \select node if connected (incoming)
y_Yra1 - x_Tif4631_Yra1 >= 0 \select node if connected (incoming)
y_Yra1 - x_Tif4632_Yra1 >= 0 \select node if connected (incoming)
y_Tif4632 - x_Tif2_Tif4632 >= 0 \select node if connected (incoming)
y_Bud14 - x_Bem2_Bud14 >= 0 \select node if connected (incoming)
y_Bud14 - x_Glc7_Bud14 >= 0 \select node if connected (incoming)
y_Pbs2 - x_Hog1_Pbs2 >= 0 \select node if connected (incoming)
y_Pbs2 - x_Kss1_Pbs2 >= 0 \select node if connected (incoming)
y_Pbs2 - x_Fus3_Pbs2 >= 0 \select node if connected (incoming)
y_Pbs2 - x_Ste11_Pbs2 >= 0 \select node if connected (incoming)
y_Pbs2 - x_Mih1_Pbs2 >= 0 \select node if connected (incoming)
y_Pbs2 - x_Ptc1_Pbs2 >= 0 \select node if connected (incoming)
y_Pbs2 - x_Nbp2_Pbs2 >= 0 \select node if connected (incoming)
y_Pbs2 - x_Sho1_Pbs2 >= 0 \select node if connected (incoming)
y_Sky1 - x_Npl3_Sky1 >= 0 \select node if connected (incoming)
y_Gpa2 - x_Ras2_Gpa2 >= 0 \select node if connected (incoming)
y_Rck2 - x_Hog1_Rck2 >= 0 \select node if connected (incoming)
y_Sub2 - x_Yra1_Sub2 >= 0 \select node if connected (incoming)

binary
y_Bem2
y_Hog1
y_Kss1
y_Fus3
y_Ste11
y_Mih1
y_Ptc1
y_Nbp2
y_Ptp2
y_Rpo21
y_Glc7
y_Npl3
y_Sho1
y_Sko1
y_Slt2
y_Ste7
y_Ras2
y_Rim15
y_Pho85
y_Sch9
y_Tpk1
y_Tif2
y_Tif4631
y_Yra1
y_Tif4632
y_Bud14
y_Pbs2
y_Sky1
y_Gpa2
y_Rck2
y_Sub2
x_Bem2_Bud14
x_Hog1_Kss1
x_Hog1_Mih1
x_Hog1_Ptc1
x_Hog1_Ptp2
x_Hog1_Rpo21
x_Hog1_Sho1
x_Hog1_Sko1
x_Hog1_Slt2
x_Hog1_Ste7
x_Hog1_Pbs2
x_Hog1_Rck2
x_Kss1_Fus3
x_Kss1_Ste11
x_Kss1_Pbs2
x_Fus3_Pbs2
x_Ste11_Pbs2
x_Mih1_Pbs2
x_Ptc1_Nbp2
x_Ptc1_Pbs2
x_Nbp2_Pbs2
x_Ptp2_Ptc1
x_Rpo21_Glc7
x_Rpo21_Npl3
x_Glc7_Bud14
x_Npl3_Sky1
x_Sho1_Ste11
x_Sho1_Pbs2
x_Sko1_Kss1
x_Sko1_Fus3
x_Slt2_Ptc1
x_Ste7_Kss1
x_Ste7_Fus3
x_Ste7_Ste11
x_Ste7_Ras2
x_Ras2_Gpa2
x_Rim15_Pho85
x_Rim15_Sch9
x_Rim15_Tpk1
x_Pho85_Mih1
x_Sch9_Ras2
x_Tpk1_Ras2
x_Tif2_Tif4631
x_Tif2_Tif4632
x_Tif4631_Npl3
x_Tif4631_Yra1
x_Yra1_Sub2
x_Tif4632_Yra1

End