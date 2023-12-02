\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 1 x_out_A + 1 x_A_B + 1 x_B_C + 1 x_C_D + 1 x_D_in + 1 x_C_in + 1 x_out_D + 1 x_D_E + 1 x_E_F + 1 x_F_in - 10 y_out - 10 y_A - 10 y_B - 10 y_C - 10 y_D - 10 y_E - 10 y_F - 10 y_in

subject to
x_out_A + x_out_D >= 1 \edges from terminal
x_C_in + x_D_in + x_F_in >= 1 \edges to terminal
y_out >= 1  \terminal
y_in >= 1  \terminal
y_A - x_A_B <= 0 \Steiner node out
y_B - x_B_C <= 0 \Steiner node out
y_C - x_C_D - x_C_in <= 0 \Steiner node out
y_D - x_D_E - x_D_in <= 0 \Steiner node out
y_E - x_E_F <= 0 \Steiner node out
y_F - x_F_in <= 0 \Steiner node out
y_A - x_out_A <= 0 \Steiner node in
y_B - x_A_B <= 0 \Steiner node in
y_C - x_B_C <= 0 \Steiner node in
y_D - x_out_D - x_C_D <= 0 \Steiner node in
y_E - x_D_E <= 0 \Steiner node in
y_F - x_E_F <= 0 \Steiner node in
y_out - x_out_A >= 0 \select node if connected (outgoing)
y_out - x_out_D >= 0 \select node if connected (outgoing)
y_A - x_A_B >= 0 \select node if connected (outgoing)
y_B - x_B_C >= 0 \select node if connected (outgoing)
y_C - x_C_D >= 0 \select node if connected (outgoing)
y_C - x_C_in >= 0 \select node if connected (outgoing)
y_D - x_D_E >= 0 \select node if connected (outgoing)
y_D - x_D_in >= 0 \select node if connected (outgoing)
y_E - x_E_F >= 0 \select node if connected (outgoing)
y_F - x_F_in >= 0 \select node if connected (outgoing)
y_A - x_out_A >= 0 \select node if connected (incoming)
y_B - x_A_B >= 0 \select node if connected (incoming)
y_C - x_B_C >= 0 \select node if connected (incoming)
y_D - x_out_D >= 0 \select node if connected (incoming)
y_D - x_C_D >= 0 \select node if connected (incoming)
y_E - x_D_E >= 0 \select node if connected (incoming)
y_F - x_E_F >= 0 \select node if connected (incoming)
y_in - x_C_in >= 0 \select node if connected (incoming)
y_in - x_D_in >= 0 \select node if connected (incoming)
y_in - x_F_in >= 0 \select node if connected (incoming)

binary
y_out
y_A
y_B
y_C
y_D
y_E
y_F
y_in
x_out_A
x_out_D
x_A_B
x_B_C
x_C_D
x_C_in
x_D_E
x_D_in
x_E_F
x_F_in

End