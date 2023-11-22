\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 2 x_2_5 + 2 x_5_7 - 2 y_2 - 5 y_5 - 7 y_7

subject to
x_2_5 >= 1 \edges from terminal
x_5_7 >= 1 \edges to terminal
y_2 >= 1  \terminal
y_7 >= 1  \terminal
y_5 - x_5_7 <= 0 \Steiner node out
y_5 - x_2_5 <= 0 \Steiner node in
y_2 - x_2_5 >= 0 \select node if connected (outgoing)
y_5 - x_5_7 >= 0 \select node if connected (outgoing)
y_5 - x_2_5 >= 0 \select node if connected (incoming)
y_7 - x_5_7 >= 0 \select node if connected (incoming)

binary
y_2
y_5
y_7
x_2_5
x_5_7

End