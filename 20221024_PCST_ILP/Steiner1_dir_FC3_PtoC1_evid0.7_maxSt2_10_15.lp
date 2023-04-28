\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 0.027 x_Ipl1_Cdc5 + 0.176 x_Cdc5_Mih1 + 0.035 x_Mih1_Hog1 + 0.175 x_Ipl1_Glc7 + 0.005 x_Glc7_Bud14 + 0.224 x_Bud14_Bem2 + 0.246 x_Glc7_Rpo21 + 0.214 x_Rpo21_Hog1 + 0.001 x_Snu66_Nam8 + 0.205 x_Nam8_Tif4632 + 0.124 x_Tif4632_Tif2 + 0.001 x_Snu66_Sto1 + 0.038 x_Sto1_Prt1 + 0.029 x_Prt1_Tif2 + 0.162 x_Sto1_Tif4631 + 0.00700000000000001 x_Tif4631_Tif2 + 0.162 x_Sto1_Tif4632 + 0.194 x_Spa2_Slt2 + 0.248 x_Slt2_Fps1 + 0.141 x_Fps1_Hog1 + 0.214 x_Slt2_Hog1 + 0.183 x_Slt2_Ptc1 + 0.186 x_Ptc1_Hog1 + 0.26 x_Slt2_Ptp2 + 0.095 x_Ptp2_Hog1 + 0.168 x_Slt2_Sko1 + 0.001 x_Sko1_Hog1 + 0.24 x_Slt2_Ste7 + 0.1 x_Ste7_Hog1 - 2.22879578185143 y_Ipl1 - 0.383515987743196 y_Cdc5 - 0.116628340365083 y_Mih1 - 0 y_Glc7 - 1.87884256040446 y_Bud14 - 1.8162863958549 y_Rpo21 - 1.12613426133923 y_Snu66 - 0 y_Nam8 - 0.517617071461765 y_Tif4632 - 0 y_Sto1 - 0.397976742486334 y_Prt1 - 0.367065514001302 y_Tif4631 - 4.20790455734154 y_Spa2 - 0 y_Slt2 - 0.795291672469519 y_Fps1 - 0.0854152687766444 y_Ptc1 - 0.173434899484957 y_Ptp2 - 1.50119873871931 y_Sko1 - 0 y_Ste7 - 1.81347223589972 y_Hog1 - 3.6868019743945 y_Bem2 - 2.893087867955 y_Tif2

subject to
x_Ipl1_Cdc5 + x_Ipl1_Glc7 >= 1 \edges from terminal
x_Snu66_Nam8 + x_Snu66_Sto1 >= 1 \edges from terminal
x_Spa2_Slt2 >= 1 \edges from terminal
x_Mih1_Hog1 + x_Rpo21_Hog1 + x_Slt2_Hog1 + x_Fps1_Hog1 + x_Ptc1_Hog1 + x_Ptp2_Hog1 + x_Sko1_Hog1 + x_Ste7_Hog1 >= 1 \edges to terminal
x_Bud14_Bem2 >= 1 \edges to terminal
x_Tif4632_Tif2 + x_Prt1_Tif2 + x_Tif4631_Tif2 >= 1 \edges to terminal
y_Ipl1 >= 1  \terminal
y_Snu66 >= 1  \terminal
y_Spa2 >= 1  \terminal
y_Hog1 >= 1  \terminal
y_Bem2 >= 1  \terminal
y_Tif2 >= 1  \terminal
y_Cdc5 - x_Cdc5_Mih1 <= 0 \Steiner node out
y_Mih1 - x_Mih1_Hog1 <= 0 \Steiner node out
y_Glc7 - x_Glc7_Bud14 - x_Glc7_Rpo21 <= 0 \Steiner node out
y_Bud14 - x_Bud14_Bem2 <= 0 \Steiner node out
y_Rpo21 - x_Rpo21_Hog1 <= 0 \Steiner node out
y_Nam8 - x_Nam8_Tif4632 <= 0 \Steiner node out
y_Tif4632 - x_Tif4632_Tif2 <= 0 \Steiner node out
y_Sto1 - x_Sto1_Tif4632 - x_Sto1_Prt1 - x_Sto1_Tif4631 <= 0 \Steiner node out
y_Prt1 - x_Prt1_Tif2 <= 0 \Steiner node out
y_Tif4631 - x_Tif4631_Tif2 <= 0 \Steiner node out
y_Slt2 - x_Slt2_Fps1 - x_Slt2_Ptc1 - x_Slt2_Ptp2 - x_Slt2_Sko1 - x_Slt2_Ste7 - x_Slt2_Hog1 <= 0 \Steiner node out
y_Fps1 - x_Fps1_Hog1 <= 0 \Steiner node out
y_Ptc1 - x_Ptc1_Hog1 <= 0 \Steiner node out
y_Ptp2 - x_Ptp2_Hog1 <= 0 \Steiner node out
y_Sko1 - x_Sko1_Hog1 <= 0 \Steiner node out
y_Ste7 - x_Ste7_Hog1 <= 0 \Steiner node out
y_Cdc5 - x_Ipl1_Cdc5 <= 0 \Steiner node in
y_Mih1 - x_Cdc5_Mih1 <= 0 \Steiner node in
y_Glc7 - x_Ipl1_Glc7 <= 0 \Steiner node in
y_Bud14 - x_Glc7_Bud14 <= 0 \Steiner node in
y_Rpo21 - x_Glc7_Rpo21 <= 0 \Steiner node in
y_Nam8 - x_Snu66_Nam8 <= 0 \Steiner node in
y_Tif4632 - x_Nam8_Tif4632 - x_Sto1_Tif4632 <= 0 \Steiner node in
y_Sto1 - x_Snu66_Sto1 <= 0 \Steiner node in
y_Prt1 - x_Sto1_Prt1 <= 0 \Steiner node in
y_Tif4631 - x_Sto1_Tif4631 <= 0 \Steiner node in
y_Slt2 - x_Spa2_Slt2 <= 0 \Steiner node in
y_Fps1 - x_Slt2_Fps1 <= 0 \Steiner node in
y_Ptc1 - x_Slt2_Ptc1 <= 0 \Steiner node in
y_Ptp2 - x_Slt2_Ptp2 <= 0 \Steiner node in
y_Sko1 - x_Slt2_Sko1 <= 0 \Steiner node in
y_Ste7 - x_Slt2_Ste7 <= 0 \Steiner node in
y_Ipl1 - x_Ipl1_Cdc5 >= 0 \select node if connected (outgoing)
y_Ipl1 - x_Ipl1_Glc7 >= 0 \select node if connected (outgoing)
y_Cdc5 - x_Cdc5_Mih1 >= 0 \select node if connected (outgoing)
y_Mih1 - x_Mih1_Hog1 >= 0 \select node if connected (outgoing)
y_Glc7 - x_Glc7_Bud14 >= 0 \select node if connected (outgoing)
y_Glc7 - x_Glc7_Rpo21 >= 0 \select node if connected (outgoing)
y_Bud14 - x_Bud14_Bem2 >= 0 \select node if connected (outgoing)
y_Rpo21 - x_Rpo21_Hog1 >= 0 \select node if connected (outgoing)
y_Snu66 - x_Snu66_Nam8 >= 0 \select node if connected (outgoing)
y_Snu66 - x_Snu66_Sto1 >= 0 \select node if connected (outgoing)
y_Nam8 - x_Nam8_Tif4632 >= 0 \select node if connected (outgoing)
y_Tif4632 - x_Tif4632_Tif2 >= 0 \select node if connected (outgoing)
y_Sto1 - x_Sto1_Tif4632 >= 0 \select node if connected (outgoing)
y_Sto1 - x_Sto1_Prt1 >= 0 \select node if connected (outgoing)
y_Sto1 - x_Sto1_Tif4631 >= 0 \select node if connected (outgoing)
y_Prt1 - x_Prt1_Tif2 >= 0 \select node if connected (outgoing)
y_Tif4631 - x_Tif4631_Tif2 >= 0 \select node if connected (outgoing)
y_Spa2 - x_Spa2_Slt2 >= 0 \select node if connected (outgoing)
y_Slt2 - x_Slt2_Fps1 >= 0 \select node if connected (outgoing)
y_Slt2 - x_Slt2_Ptc1 >= 0 \select node if connected (outgoing)
y_Slt2 - x_Slt2_Ptp2 >= 0 \select node if connected (outgoing)
y_Slt2 - x_Slt2_Sko1 >= 0 \select node if connected (outgoing)
y_Slt2 - x_Slt2_Ste7 >= 0 \select node if connected (outgoing)
y_Slt2 - x_Slt2_Hog1 >= 0 \select node if connected (outgoing)
y_Fps1 - x_Fps1_Hog1 >= 0 \select node if connected (outgoing)
y_Ptc1 - x_Ptc1_Hog1 >= 0 \select node if connected (outgoing)
y_Ptp2 - x_Ptp2_Hog1 >= 0 \select node if connected (outgoing)
y_Sko1 - x_Sko1_Hog1 >= 0 \select node if connected (outgoing)
y_Ste7 - x_Ste7_Hog1 >= 0 \select node if connected (outgoing)
y_Cdc5 - x_Ipl1_Cdc5 >= 0 \select node if connected (incoming)
y_Mih1 - x_Cdc5_Mih1 >= 0 \select node if connected (incoming)
y_Glc7 - x_Ipl1_Glc7 >= 0 \select node if connected (incoming)
y_Bud14 - x_Glc7_Bud14 >= 0 \select node if connected (incoming)
y_Rpo21 - x_Glc7_Rpo21 >= 0 \select node if connected (incoming)
y_Nam8 - x_Snu66_Nam8 >= 0 \select node if connected (incoming)
y_Tif4632 - x_Nam8_Tif4632 >= 0 \select node if connected (incoming)
y_Tif4632 - x_Sto1_Tif4632 >= 0 \select node if connected (incoming)
y_Sto1 - x_Snu66_Sto1 >= 0 \select node if connected (incoming)
y_Prt1 - x_Sto1_Prt1 >= 0 \select node if connected (incoming)
y_Tif4631 - x_Sto1_Tif4631 >= 0 \select node if connected (incoming)
y_Slt2 - x_Spa2_Slt2 >= 0 \select node if connected (incoming)
y_Fps1 - x_Slt2_Fps1 >= 0 \select node if connected (incoming)
y_Ptc1 - x_Slt2_Ptc1 >= 0 \select node if connected (incoming)
y_Ptp2 - x_Slt2_Ptp2 >= 0 \select node if connected (incoming)
y_Sko1 - x_Slt2_Sko1 >= 0 \select node if connected (incoming)
y_Ste7 - x_Slt2_Ste7 >= 0 \select node if connected (incoming)
y_Hog1 - x_Mih1_Hog1 >= 0 \select node if connected (incoming)
y_Hog1 - x_Rpo21_Hog1 >= 0 \select node if connected (incoming)
y_Hog1 - x_Slt2_Hog1 >= 0 \select node if connected (incoming)
y_Hog1 - x_Fps1_Hog1 >= 0 \select node if connected (incoming)
y_Hog1 - x_Ptc1_Hog1 >= 0 \select node if connected (incoming)
y_Hog1 - x_Ptp2_Hog1 >= 0 \select node if connected (incoming)
y_Hog1 - x_Sko1_Hog1 >= 0 \select node if connected (incoming)
y_Hog1 - x_Ste7_Hog1 >= 0 \select node if connected (incoming)
y_Bem2 - x_Bud14_Bem2 >= 0 \select node if connected (incoming)
y_Tif2 - x_Tif4632_Tif2 >= 0 \select node if connected (incoming)
y_Tif2 - x_Prt1_Tif2 >= 0 \select node if connected (incoming)
y_Tif2 - x_Tif4631_Tif2 >= 0 \select node if connected (incoming)

binary
y_Ipl1
y_Cdc5
y_Mih1
y_Glc7
y_Bud14
y_Rpo21
y_Snu66
y_Nam8
y_Tif4632
y_Sto1
y_Prt1
y_Tif4631
y_Spa2
y_Slt2
y_Fps1
y_Ptc1
y_Ptp2
y_Sko1
y_Ste7
y_Hog1
y_Bem2
y_Tif2
x_Ipl1_Cdc5
x_Ipl1_Glc7
x_Cdc5_Mih1
x_Mih1_Hog1
x_Glc7_Bud14
x_Glc7_Rpo21
x_Bud14_Bem2
x_Rpo21_Hog1
x_Snu66_Nam8
x_Snu66_Sto1
x_Nam8_Tif4632
x_Tif4632_Tif2
x_Sto1_Tif4632
x_Sto1_Prt1
x_Sto1_Tif4631
x_Prt1_Tif2
x_Tif4631_Tif2
x_Spa2_Slt2
x_Slt2_Fps1
x_Slt2_Ptc1
x_Slt2_Ptp2
x_Slt2_Sko1
x_Slt2_Ste7
x_Slt2_Hog1
x_Fps1_Hog1
x_Ptc1_Hog1
x_Ptp2_Hog1
x_Sko1_Hog1
x_Ste7_Hog1

End