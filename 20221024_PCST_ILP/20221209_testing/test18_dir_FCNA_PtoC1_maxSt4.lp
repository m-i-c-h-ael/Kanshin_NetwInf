\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 1 x_a2_b2 + 1 x_b2_i + 1 x_i_a1 + 1 x_a1_b1 + 1 x_b1_d1 + 1 x_a1_c1 + 0.9 x_c1_d1 + 1 x_a1_d1 + 0.9 x_a2_c2 + 1 x_c2_b2 + 1 x_c2_d2 + 1 x_d2_b2 + 1 x_d2_e2 + 1 x_e2_b2 + 1 x_c2_e2 + 1 x_e2_d2 + 1 x_a2_d2 + 1 x_d2_c2 + 1 x_e2_c2 + 1 x_a2_e2 + 1 x_b1_c1 + 1 x_c1_b1 - 1 y_a2 - 1 y_b2 - 1 y_i - 1 y_a1 - 1 y_b1 - 1 y_c1 - 1 y_c2 - 1 y_d2 - 1 y_e2 - 1 y_d1

subject to
x_a2_b2 + x_a2_c2 + x_a2_d2 + x_a2_e2 >= 1 \edges from terminal
x_b2_i >= 1 \edges from terminal
x_b2_i >= 1 \edges to terminal
x_a1_d1 + x_b1_d1 + x_c1_d1 >= 1 \edges to terminal
y_a2 >= 1  \terminal
y_b2 >= 1  \terminal
y_i >= 1  \terminal
y_d1 >= 1  \terminal
y_a1 - x_a1_b1 - x_a1_c1 - x_a1_d1 <= 0 \Steiner node out
y_b1 - x_b1_c1 - x_b1_d1 <= 0 \Steiner node out
y_c1 - x_c1_b1 - x_c1_d1 <= 0 \Steiner node out
y_c2 - x_c2_b2 - x_c2_d2 - x_c2_e2 <= 0 \Steiner node out
y_d2 - x_d2_b2 - x_d2_c2 - x_d2_e2 <= 0 \Steiner node out
y_e2 - x_e2_b2 - x_e2_c2 - x_e2_d2 <= 0 \Steiner node out
y_a1 - x_i_a1 <= 0 \Steiner node in
y_b1 - x_a1_b1 - x_c1_b1 <= 0 \Steiner node in
y_c1 - x_a1_c1 - x_b1_c1 <= 0 \Steiner node in
y_c2 - x_a2_c2 - x_d2_c2 - x_e2_c2 <= 0 \Steiner node in
y_d2 - x_a2_d2 - x_c2_d2 - x_e2_d2 <= 0 \Steiner node in
y_e2 - x_a2_e2 - x_c2_e2 - x_d2_e2 <= 0 \Steiner node in
y_a2 - x_a2_b2 >= 0 \select node if connected (outgoing)
y_a2 - x_a2_c2 >= 0 \select node if connected (outgoing)
y_a2 - x_a2_d2 >= 0 \select node if connected (outgoing)
y_a2 - x_a2_e2 >= 0 \select node if connected (outgoing)
y_b2 - x_b2_i >= 0 \select node if connected (outgoing)
y_i - x_i_a1 >= 0 \select node if connected (outgoing)
y_a1 - x_a1_b1 >= 0 \select node if connected (outgoing)
y_a1 - x_a1_c1 >= 0 \select node if connected (outgoing)
y_a1 - x_a1_d1 >= 0 \select node if connected (outgoing)
y_b1 - x_b1_c1 >= 0 \select node if connected (outgoing)
y_b1 - x_b1_d1 >= 0 \select node if connected (outgoing)
y_c1 - x_c1_b1 >= 0 \select node if connected (outgoing)
y_c1 - x_c1_d1 >= 0 \select node if connected (outgoing)
y_c2 - x_c2_b2 >= 0 \select node if connected (outgoing)
y_c2 - x_c2_d2 >= 0 \select node if connected (outgoing)
y_c2 - x_c2_e2 >= 0 \select node if connected (outgoing)
y_d2 - x_d2_b2 >= 0 \select node if connected (outgoing)
y_d2 - x_d2_c2 >= 0 \select node if connected (outgoing)
y_d2 - x_d2_e2 >= 0 \select node if connected (outgoing)
y_e2 - x_e2_b2 >= 0 \select node if connected (outgoing)
y_e2 - x_e2_c2 >= 0 \select node if connected (outgoing)
y_e2 - x_e2_d2 >= 0 \select node if connected (outgoing)
y_b2 - x_a2_b2 >= 0 \select node if connected (incoming)
y_b2 - x_c2_b2 >= 0 \select node if connected (incoming)
y_b2 - x_d2_b2 >= 0 \select node if connected (incoming)
y_b2 - x_e2_b2 >= 0 \select node if connected (incoming)
y_i - x_b2_i >= 0 \select node if connected (incoming)
y_a1 - x_i_a1 >= 0 \select node if connected (incoming)
y_b1 - x_a1_b1 >= 0 \select node if connected (incoming)
y_b1 - x_c1_b1 >= 0 \select node if connected (incoming)
y_c1 - x_a1_c1 >= 0 \select node if connected (incoming)
y_c1 - x_b1_c1 >= 0 \select node if connected (incoming)
y_c2 - x_a2_c2 >= 0 \select node if connected (incoming)
y_c2 - x_d2_c2 >= 0 \select node if connected (incoming)
y_c2 - x_e2_c2 >= 0 \select node if connected (incoming)
y_d2 - x_a2_d2 >= 0 \select node if connected (incoming)
y_d2 - x_c2_d2 >= 0 \select node if connected (incoming)
y_d2 - x_e2_d2 >= 0 \select node if connected (incoming)
y_e2 - x_a2_e2 >= 0 \select node if connected (incoming)
y_e2 - x_c2_e2 >= 0 \select node if connected (incoming)
y_e2 - x_d2_e2 >= 0 \select node if connected (incoming)
y_d1 - x_a1_d1 >= 0 \select node if connected (incoming)
y_d1 - x_b1_d1 >= 0 \select node if connected (incoming)
y_d1 - x_c1_d1 >= 0 \select node if connected (incoming)
x_b1_c1 + x_c1_b1 <= 1 \avoid 2-node cycles
x_c2_d2 + x_d2_c2 <= 1 \avoid 2-node cycles
x_c2_e2 + x_e2_c2 <= 1 \avoid 2-node cycles
x_c2_d2 + x_d2_e2 + x_e2_c2 <= 2 \avoid 3-node cycles
x_c2_d2 + x_d2_e2 + x_e2_d2 + x_d2_c2 <= 3 \avoid 4-node cycles
x_c2_e2 + x_e2_d2 + x_d2_c2 <= 2 \avoid 3-node cycles
x_c2_e2 + x_e2_d2 + x_d2_e2 + x_e2_c2 <= 3 \avoid 4-node cycles
x_d2_e2 + x_e2_d2 <= 1 \avoid 2-node cycles

binary
y_a2
y_b2
y_i
y_a1
y_b1
y_c1
y_c2
y_d2
y_e2
y_d1
x_a2_b2
x_a2_c2
x_a2_d2
x_a2_e2
x_b2_i
x_i_a1
x_a1_b1
x_a1_c1
x_a1_d1
x_b1_c1
x_b1_d1
x_c1_b1
x_c1_d1
x_c2_b2
x_c2_d2
x_c2_e2
x_d2_b2
x_d2_c2
x_d2_e2
x_e2_b2
x_e2_c2
x_e2_d2

End