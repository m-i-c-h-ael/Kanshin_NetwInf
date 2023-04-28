\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 0.0659999999999999 x_Rad23_Rpn11 + 0.246 x_Rpn11_Nip1 + 0.15 x_Nip1_Ded1 + 0.037 x_Rad23_Rps31 + 0.043 x_Rps31_Nip1 + 0.078 x_Rps31_Prt1 + 0.09 x_Prt1_Ded1 + 0.039 x_Rps31_Rpg1 + 0.099 x_Rpg1_Ded1 + 0.078 x_Rps31_Tif35 + 0.105 x_Tif35_Ded1 + 0.038 x_Prm15_Cof1 + 0.004 x_Cof1_Act1 + 0.195 x_Act1_Crn1 + 0.075 x_Cof1_Arp9 + 0.04 x_Arp9_Rsc9 - 1.54767341365104 y_Rad23 - 0 y_Rpn11 - 0 y_Nip1 - 0.227967035078496 y_Rps31 - 0.397976742486334 y_Prt1 - 0.414174219075414 y_Rpg1 - 0.18472471674679 y_Tif35 - 0 y_Prm15 - 0 y_Cof1 - 0 y_Act1 - 0 y_Arp9 - 1.21987702051012 y_Ded1 - 2.56841741690909 y_Crn1 - 2.29148584420608 y_Rsc9

subject to
x_Rad23_Rpn11 + x_Rad23_Rps31 >= 1 \edges from terminal
x_Prm15_Cof1 >= 1 \edges from terminal
x_Nip1_Ded1 + x_Prt1_Ded1 + x_Rpg1_Ded1 + x_Tif35_Ded1 >= 1 \edges to terminal
x_Act1_Crn1 >= 1 \edges to terminal
x_Arp9_Rsc9 >= 1 \edges to terminal
y_Rad23 >= 1  \terminal
y_Prm15 >= 1  \terminal
y_Ded1 >= 1  \terminal
y_Crn1 >= 1  \terminal
y_Rsc9 >= 1  \terminal
y_Rpn11 - x_Rpn11_Nip1 <= 0 \Steiner node out
y_Nip1 - x_Nip1_Ded1 <= 0 \Steiner node out
y_Rps31 - x_Rps31_Nip1 - x_Rps31_Prt1 - x_Rps31_Rpg1 - x_Rps31_Tif35 <= 0 \Steiner node out
y_Prt1 - x_Prt1_Ded1 <= 0 \Steiner node out
y_Rpg1 - x_Rpg1_Ded1 <= 0 \Steiner node out
y_Tif35 - x_Tif35_Ded1 <= 0 \Steiner node out
y_Cof1 - x_Cof1_Act1 - x_Cof1_Arp9 <= 0 \Steiner node out
y_Act1 - x_Act1_Crn1 <= 0 \Steiner node out
y_Arp9 - x_Arp9_Rsc9 <= 0 \Steiner node out
y_Rpn11 - x_Rad23_Rpn11 <= 0 \Steiner node in
y_Nip1 - x_Rpn11_Nip1 - x_Rps31_Nip1 <= 0 \Steiner node in
y_Rps31 - x_Rad23_Rps31 <= 0 \Steiner node in
y_Prt1 - x_Rps31_Prt1 <= 0 \Steiner node in
y_Rpg1 - x_Rps31_Rpg1 <= 0 \Steiner node in
y_Tif35 - x_Rps31_Tif35 <= 0 \Steiner node in
y_Cof1 - x_Prm15_Cof1 <= 0 \Steiner node in
y_Act1 - x_Cof1_Act1 <= 0 \Steiner node in
y_Arp9 - x_Cof1_Arp9 <= 0 \Steiner node in
y_Rad23 - x_Rad23_Rpn11 >= 0 \select node if connected (outgoing)
y_Rad23 - x_Rad23_Rps31 >= 0 \select node if connected (outgoing)
y_Rpn11 - x_Rpn11_Nip1 >= 0 \select node if connected (outgoing)
y_Nip1 - x_Nip1_Ded1 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Nip1 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Prt1 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpg1 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Tif35 >= 0 \select node if connected (outgoing)
y_Prt1 - x_Prt1_Ded1 >= 0 \select node if connected (outgoing)
y_Rpg1 - x_Rpg1_Ded1 >= 0 \select node if connected (outgoing)
y_Tif35 - x_Tif35_Ded1 >= 0 \select node if connected (outgoing)
y_Prm15 - x_Prm15_Cof1 >= 0 \select node if connected (outgoing)
y_Cof1 - x_Cof1_Act1 >= 0 \select node if connected (outgoing)
y_Cof1 - x_Cof1_Arp9 >= 0 \select node if connected (outgoing)
y_Act1 - x_Act1_Crn1 >= 0 \select node if connected (outgoing)
y_Arp9 - x_Arp9_Rsc9 >= 0 \select node if connected (outgoing)
y_Rpn11 - x_Rad23_Rpn11 >= 0 \select node if connected (incoming)
y_Nip1 - x_Rpn11_Nip1 >= 0 \select node if connected (incoming)
y_Nip1 - x_Rps31_Nip1 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rad23_Rps31 >= 0 \select node if connected (incoming)
y_Prt1 - x_Rps31_Prt1 >= 0 \select node if connected (incoming)
y_Rpg1 - x_Rps31_Rpg1 >= 0 \select node if connected (incoming)
y_Tif35 - x_Rps31_Tif35 >= 0 \select node if connected (incoming)
y_Cof1 - x_Prm15_Cof1 >= 0 \select node if connected (incoming)
y_Act1 - x_Cof1_Act1 >= 0 \select node if connected (incoming)
y_Arp9 - x_Cof1_Arp9 >= 0 \select node if connected (incoming)
y_Ded1 - x_Nip1_Ded1 >= 0 \select node if connected (incoming)
y_Ded1 - x_Prt1_Ded1 >= 0 \select node if connected (incoming)
y_Ded1 - x_Rpg1_Ded1 >= 0 \select node if connected (incoming)
y_Ded1 - x_Tif35_Ded1 >= 0 \select node if connected (incoming)
y_Crn1 - x_Act1_Crn1 >= 0 \select node if connected (incoming)
y_Rsc9 - x_Arp9_Rsc9 >= 0 \select node if connected (incoming)

binary
y_Rad23
y_Rpn11
y_Nip1
y_Rps31
y_Prt1
y_Rpg1
y_Tif35
y_Prm15
y_Cof1
y_Act1
y_Arp9
y_Ded1
y_Crn1
y_Rsc9
x_Rad23_Rpn11
x_Rad23_Rps31
x_Rpn11_Nip1
x_Nip1_Ded1
x_Rps31_Nip1
x_Rps31_Prt1
x_Rps31_Rpg1
x_Rps31_Tif35
x_Prt1_Ded1
x_Rpg1_Ded1
x_Tif35_Ded1
x_Prm15_Cof1
x_Cof1_Act1
x_Cof1_Arp9
x_Act1_Crn1
x_Arp9_Rsc9

End