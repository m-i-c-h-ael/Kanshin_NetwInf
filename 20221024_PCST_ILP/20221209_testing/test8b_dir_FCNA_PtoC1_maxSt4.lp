\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 1 x_B_A + 2 x_A_D + 100 x_D_G + 200 x_B_E + 2 x_E_G + 100 x_B_H + 1 x_H_E + 2 x_B_I + 2 x_I_A - 1 y_B - 1 y_A - 1 y_D - 1 y_E - 1 y_H - 1 y_I - 1 y_G

subject to
x_B_A + x_B_E + x_B_H + x_B_I >= 1 \edges from terminal
x_D_G + x_E_G >= 1 \edges to terminal
y_B >= 1  \terminal
y_G >= 1  \terminal
y_A - x_A_D <= 0 \Steiner node out
y_D - x_D_G <= 0 \Steiner node out
y_E - x_E_G <= 0 \Steiner node out
y_H - x_H_E <= 0 \Steiner node out
y_I - x_I_A <= 0 \Steiner node out
y_A - x_B_A - x_I_A <= 0 \Steiner node in
y_D - x_A_D <= 0 \Steiner node in
y_E - x_B_E - x_H_E <= 0 \Steiner node in
y_H - x_B_H <= 0 \Steiner node in
y_I - x_B_I <= 0 \Steiner node in
y_B - x_B_A >= 0 \select node if connected (outgoing)
y_B - x_B_E >= 0 \select node if connected (outgoing)
y_B - x_B_H >= 0 \select node if connected (outgoing)
y_B - x_B_I >= 0 \select node if connected (outgoing)
y_A - x_A_D >= 0 \select node if connected (outgoing)
y_D - x_D_G >= 0 \select node if connected (outgoing)
y_E - x_E_G >= 0 \select node if connected (outgoing)
y_H - x_H_E >= 0 \select node if connected (outgoing)
y_I - x_I_A >= 0 \select node if connected (outgoing)
y_A - x_B_A >= 0 \select node if connected (incoming)
y_A - x_I_A >= 0 \select node if connected (incoming)
y_D - x_A_D >= 0 \select node if connected (incoming)
y_E - x_B_E >= 0 \select node if connected (incoming)
y_E - x_H_E >= 0 \select node if connected (incoming)
y_H - x_B_H >= 0 \select node if connected (incoming)
y_I - x_B_I >= 0 \select node if connected (incoming)
y_G - x_D_G >= 0 \select node if connected (incoming)
y_G - x_E_G >= 0 \select node if connected (incoming)

binary
y_B
y_A
y_D
y_E
y_H
y_I
y_G
x_B_A
x_B_E
x_B_H
x_B_I
x_A_D
x_D_G
x_E_G
x_H_E
x_I_A

End