\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 3 x_a17_a4 + 2 x_a4_a1 + 3 x_a1_a3 + 7 x_a1_a7 + 1 x_a1_a33 + 3 x_a12_a11 + 10 x_a11_a4 - 5 y_a17 - 8 y_a4 - 6 y_a1 - 4 y_a12 - 5 y_a11 - 6 y_a3 - 3 y_a7 - 7 y_a33

subject to
x_a17_a4 >= 1 \edges from terminal
x_a12_a11 >= 1 \edges from terminal
x_a11_a4 >= 1 \edges from terminal
x_a1_a3 >= 1 \edges to terminal
x_a1_a7 >= 1 \edges to terminal
x_a1_a33 >= 1 \edges to terminal
y_a17 >= 1  \terminal
y_a12 >= 1  \terminal
y_a11 >= 1  \terminal
y_a3 >= 1  \terminal
y_a7 >= 1  \terminal
y_a33 >= 1  \terminal
y_a4 - x_a4_a1 <= 0 \Steiner node out
y_a1 - x_a1_a3 - x_a1_a7 - x_a1_a33 <= 0 \Steiner node out
y_a4 - x_a17_a4 - x_a11_a4 <= 0 \Steiner node in
y_a1 - x_a4_a1 <= 0 \Steiner node in
y_a17 - x_a17_a4 >= 0 \select node if connected (outgoing)
y_a4 - x_a4_a1 >= 0 \select node if connected (outgoing)
y_a1 - x_a1_a3 >= 0 \select node if connected (outgoing)
y_a1 - x_a1_a7 >= 0 \select node if connected (outgoing)
y_a1 - x_a1_a33 >= 0 \select node if connected (outgoing)
y_a12 - x_a12_a11 >= 0 \select node if connected (outgoing)
y_a11 - x_a11_a4 >= 0 \select node if connected (outgoing)
y_a4 - x_a17_a4 >= 0 \select node if connected (incoming)
y_a4 - x_a11_a4 >= 0 \select node if connected (incoming)
y_a1 - x_a4_a1 >= 0 \select node if connected (incoming)
y_a11 - x_a12_a11 >= 0 \select node if connected (incoming)
y_a3 - x_a1_a3 >= 0 \select node if connected (incoming)
y_a7 - x_a1_a7 >= 0 \select node if connected (incoming)
y_a33 - x_a1_a33 >= 0 \select node if connected (incoming)

binary
y_a17
y_a4
y_a1
y_a12
y_a11
y_a3
y_a7
y_a33
x_a17_a4
x_a4_a1
x_a1_a3
x_a1_a7
x_a1_a33
x_a12_a11
x_a11_a4

End