\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 0.178 x_Gpa2_Ras2 + 0.115 x_Ras2_Tpk1 + 0.279 x_Tpk1_Pfk26 + 0.215 x_Ras2_Tpk2 + 0.279 x_Tpk2_Pfk26 + 0.023 x_Pbs2_Mih1 + 0.043 x_Mih1_Hsl1 + 0.229 x_Hsl1_Zds1 - 1.01894679922407 y_Gpa2 - 0.080694252350959 y_Ras2 - 0 y_Tpk1 - 0 y_Tpk2 - 1.7465565764114 y_Pbs2 - 0.116628340365083 y_Mih1 - 0.478166934331359 y_Hsl1 - 1.65431116623074 y_Pfk26 - 1.81054691116939 y_Zds1

subject to
x_Gpa2_Ras2 >= 1 \edges from terminal
x_Pbs2_Mih1 >= 1 \edges from terminal
x_Tpk1_Pfk26 + x_Tpk2_Pfk26 >= 1 \edges to terminal
x_Hsl1_Zds1 >= 1 \edges to terminal
y_Gpa2 >= 1  \terminal
y_Pbs2 >= 1  \terminal
y_Pfk26 >= 1  \terminal
y_Zds1 >= 1  \terminal
y_Ras2 - x_Ras2_Tpk1 - x_Ras2_Tpk2 <= 0 \Steiner node out
y_Tpk1 - x_Tpk1_Pfk26 <= 0 \Steiner node out
y_Tpk2 - x_Tpk2_Pfk26 <= 0 \Steiner node out
y_Mih1 - x_Mih1_Hsl1 <= 0 \Steiner node out
y_Hsl1 - x_Hsl1_Zds1 <= 0 \Steiner node out
y_Ras2 - x_Gpa2_Ras2 <= 0 \Steiner node in
y_Tpk1 - x_Ras2_Tpk1 <= 0 \Steiner node in
y_Tpk2 - x_Ras2_Tpk2 <= 0 \Steiner node in
y_Mih1 - x_Pbs2_Mih1 <= 0 \Steiner node in
y_Hsl1 - x_Mih1_Hsl1 <= 0 \Steiner node in
y_Gpa2 - x_Gpa2_Ras2 >= 0 \select node if connected (outgoing)
y_Ras2 - x_Ras2_Tpk1 >= 0 \select node if connected (outgoing)
y_Ras2 - x_Ras2_Tpk2 >= 0 \select node if connected (outgoing)
y_Tpk1 - x_Tpk1_Pfk26 >= 0 \select node if connected (outgoing)
y_Tpk2 - x_Tpk2_Pfk26 >= 0 \select node if connected (outgoing)
y_Pbs2 - x_Pbs2_Mih1 >= 0 \select node if connected (outgoing)
y_Mih1 - x_Mih1_Hsl1 >= 0 \select node if connected (outgoing)
y_Hsl1 - x_Hsl1_Zds1 >= 0 \select node if connected (outgoing)
y_Ras2 - x_Gpa2_Ras2 >= 0 \select node if connected (incoming)
y_Tpk1 - x_Ras2_Tpk1 >= 0 \select node if connected (incoming)
y_Tpk2 - x_Ras2_Tpk2 >= 0 \select node if connected (incoming)
y_Mih1 - x_Pbs2_Mih1 >= 0 \select node if connected (incoming)
y_Hsl1 - x_Mih1_Hsl1 >= 0 \select node if connected (incoming)
y_Pfk26 - x_Tpk1_Pfk26 >= 0 \select node if connected (incoming)
y_Pfk26 - x_Tpk2_Pfk26 >= 0 \select node if connected (incoming)
y_Zds1 - x_Hsl1_Zds1 >= 0 \select node if connected (incoming)

binary
y_Gpa2
y_Ras2
y_Tpk1
y_Tpk2
y_Pbs2
y_Mih1
y_Hsl1
y_Pfk26
y_Zds1
x_Gpa2_Ras2
x_Ras2_Tpk1
x_Ras2_Tpk2
x_Tpk1_Pfk26
x_Tpk2_Pfk26
x_Pbs2_Mih1
x_Mih1_Hsl1
x_Hsl1_Zds1

End