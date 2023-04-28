\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 3 x_START1_SFP1 + 6 x_SFP1_SKP1 + 10 x_SKP1_SNU66 + 10 x_SNU66_SPE4 + 6 x_SPE4_END1 + 4 x_START2_MMS1 + 4 x_MMS1_MRS6 + 10 x_MRS6_MTR2 + 9 x_MTR2_END2 - 4 y_START1 - 10 y_SFP1 - 5 y_SKP1 - 2 y_SNU66 - 7 y_SPE4 - 9 y_START2 - 5 y_MMS1 - 9 y_MRS6 - 10 y_MTR2 - 9 y_END1 - 7 y_END2

subject to
x_START1_SFP1 >= 1 \edges from terminal
x_START2_MMS1 >= 1 \edges from terminal
x_SPE4_END1 >= 1 \edges to terminal
x_MTR2_END2 >= 1 \edges to terminal
y_START1 >= 1  \terminal
y_START2 >= 1  \terminal
y_END1 >= 1  \terminal
y_END2 >= 1  \terminal
y_SFP1 - x_SFP1_SKP1 <= 0 \Steiner node out
y_SKP1 - x_SKP1_SNU66 <= 0 \Steiner node out
y_SNU66 - x_SNU66_SPE4 <= 0 \Steiner node out
y_SPE4 - x_SPE4_END1 <= 0 \Steiner node out
y_MMS1 - x_MMS1_MRS6 <= 0 \Steiner node out
y_MRS6 - x_MRS6_MTR2 <= 0 \Steiner node out
y_MTR2 - x_MTR2_END2 <= 0 \Steiner node out
y_SFP1 - x_START1_SFP1 <= 0 \Steiner node in
y_SKP1 - x_SFP1_SKP1 <= 0 \Steiner node in
y_SNU66 - x_SKP1_SNU66 <= 0 \Steiner node in
y_SPE4 - x_SNU66_SPE4 <= 0 \Steiner node in
y_MMS1 - x_START2_MMS1 <= 0 \Steiner node in
y_MRS6 - x_MMS1_MRS6 <= 0 \Steiner node in
y_MTR2 - x_MRS6_MTR2 <= 0 \Steiner node in
y_START1 - x_START1_SFP1 >= 0 \select node if connected (outgoing)
y_SFP1 - x_SFP1_SKP1 >= 0 \select node if connected (outgoing)
y_SKP1 - x_SKP1_SNU66 >= 0 \select node if connected (outgoing)
y_SNU66 - x_SNU66_SPE4 >= 0 \select node if connected (outgoing)
y_SPE4 - x_SPE4_END1 >= 0 \select node if connected (outgoing)
y_START2 - x_START2_MMS1 >= 0 \select node if connected (outgoing)
y_MMS1 - x_MMS1_MRS6 >= 0 \select node if connected (outgoing)
y_MRS6 - x_MRS6_MTR2 >= 0 \select node if connected (outgoing)
y_MTR2 - x_MTR2_END2 >= 0 \select node if connected (outgoing)
y_SFP1 - x_START1_SFP1 >= 0 \select node if connected (incoming)
y_SKP1 - x_SFP1_SKP1 >= 0 \select node if connected (incoming)
y_SNU66 - x_SKP1_SNU66 >= 0 \select node if connected (incoming)
y_SPE4 - x_SNU66_SPE4 >= 0 \select node if connected (incoming)
y_MMS1 - x_START2_MMS1 >= 0 \select node if connected (incoming)
y_MRS6 - x_MMS1_MRS6 >= 0 \select node if connected (incoming)
y_MTR2 - x_MRS6_MTR2 >= 0 \select node if connected (incoming)
y_END1 - x_SPE4_END1 >= 0 \select node if connected (incoming)
y_END2 - x_MTR2_END2 >= 0 \select node if connected (incoming)

binary
y_START1
y_SFP1
y_SKP1
y_SNU66
y_SPE4
y_START2
y_MMS1
y_MRS6
y_MTR2
y_END1
y_END2
x_START1_SFP1
x_SFP1_SKP1
x_SKP1_SNU66
x_SNU66_SPE4
x_SPE4_END1
x_START2_MMS1
x_MMS1_MRS6
x_MRS6_MTR2
x_MTR2_END2

End