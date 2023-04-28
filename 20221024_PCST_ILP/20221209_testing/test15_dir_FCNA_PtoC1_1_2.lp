\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 1 x_START1_MFP1 + 100 x_MFP1_MFP2 + 100 x_MFP2_MFP3 + 100 x_MFP3_END1 - 1 y_START1 - 1 y_MFP1 - 1 y_MFP2 - 1 y_MFP3 - 1 y_END1

subject to
x_START1_MFP1 >= 1 \edges from terminal
x_MFP3_END1 >= 1 \edges to terminal
y_START1 >= 1  \terminal
y_END1 >= 1  \terminal
y_MFP1 - x_MFP1_MFP2 <= 0 \Steiner node out
y_MFP2 - x_MFP2_MFP3 <= 0 \Steiner node out
y_MFP3 - x_MFP3_END1 <= 0 \Steiner node out
y_MFP1 - x_START1_MFP1 <= 0 \Steiner node in
y_MFP2 - x_MFP1_MFP2 <= 0 \Steiner node in
y_MFP3 - x_MFP2_MFP3 <= 0 \Steiner node in
y_START1 - x_START1_MFP1 >= 0 \select node if connected (outgoing)
y_MFP1 - x_MFP1_MFP2 >= 0 \select node if connected (outgoing)
y_MFP2 - x_MFP2_MFP3 >= 0 \select node if connected (outgoing)
y_MFP3 - x_MFP3_END1 >= 0 \select node if connected (outgoing)
y_MFP1 - x_START1_MFP1 >= 0 \select node if connected (incoming)
y_MFP2 - x_MFP1_MFP2 >= 0 \select node if connected (incoming)
y_MFP3 - x_MFP2_MFP3 >= 0 \select node if connected (incoming)
y_END1 - x_MFP3_END1 >= 0 \select node if connected (incoming)

binary
y_START1
y_MFP1
y_MFP2
y_MFP3
y_END1
x_START1_MFP1
x_MFP1_MFP2
x_MFP2_MFP3
x_MFP3_END1

End