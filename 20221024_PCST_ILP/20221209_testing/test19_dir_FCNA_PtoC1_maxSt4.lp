\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 1 x_A_B + 1 x_B_C - 1 y_A - 1 y_B - 1 y_C

subject to
x_A_B >= 1 \edges from terminal
x_B_C >= 1 \edges from terminal
x_B_C >= 1 \edges to terminal
y_A >= 1  \terminal
y_B >= 1  \terminal
y_C >= 1  \terminal
y_A - x_A_B >= 0 \select node if connected (outgoing)
y_B - x_B_C >= 0 \select node if connected (outgoing)
y_B - x_A_B >= 0 \select node if connected (incoming)
y_C - x_B_C >= 0 \select node if connected (incoming)

binary
y_A
y_B
y_C
x_A_B
x_B_C

End