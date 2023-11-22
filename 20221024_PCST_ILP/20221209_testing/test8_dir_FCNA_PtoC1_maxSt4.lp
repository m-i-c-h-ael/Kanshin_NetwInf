\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 2 x_B_E + 2 x_E_G - 2 y_B - 5 y_E - 7 y_G

subject to
x_B_E >= 1 \edges from terminal
x_E_G >= 1 \edges to terminal
y_B >= 1  \terminal
y_G >= 1  \terminal
y_E - x_E_G <= 0 \Steiner node out
y_E - x_B_E <= 0 \Steiner node in
y_B - x_B_E >= 0 \select node if connected (outgoing)
y_E - x_E_G >= 0 \select node if connected (outgoing)
y_E - x_B_E >= 0 \select node if connected (incoming)
y_G - x_E_G >= 0 \select node if connected (incoming)

binary
y_B
y_E
y_G
x_B_E
x_E_G

End