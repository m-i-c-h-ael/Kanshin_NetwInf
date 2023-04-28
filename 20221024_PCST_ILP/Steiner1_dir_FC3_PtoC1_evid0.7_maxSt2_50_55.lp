\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 0.251 x_Rtn1_Vps54 + 0.254 x_Vps54_Vps51 + 0.244 x_Vps51_Vps53 + 0.001 x_Vps54_Vps52 + 0.004 x_Vps52_Vps53 + 0.001 x_Vps54_Vps53 + 0.001 x_Zrc1_Cot1 + 0.226 x_Cot1_Fmp42 - 1.3606164549753 y_Rtn1 - 0.322061639797625 y_Vps54 - 0 y_Vps51 - 0 y_Vps52 - 1.57931769417956 y_Zrc1 - 1.84374380445618 y_Cot1 - 1.12224309419982 y_Vps53 - 0 y_Fmp42

subject to
x_Rtn1_Vps54 >= 1 \edges from terminal
x_Zrc1_Cot1 >= 1 \edges from terminal
x_Vps54_Vps53 + x_Vps51_Vps53 + x_Vps52_Vps53 >= 1 \edges to terminal
x_Cot1_Fmp42 >= 1 \edges to terminal
y_Rtn1 >= 1  \terminal
y_Zrc1 >= 1  \terminal
y_Vps53 >= 1  \terminal
y_Fmp42 >= 1  \terminal
y_Vps54 - x_Vps54_Vps51 - x_Vps54_Vps52 - x_Vps54_Vps53 <= 0 \Steiner node out
y_Vps51 - x_Vps51_Vps53 <= 0 \Steiner node out
y_Vps52 - x_Vps52_Vps53 <= 0 \Steiner node out
y_Cot1 - x_Cot1_Fmp42 <= 0 \Steiner node out
y_Vps54 - x_Rtn1_Vps54 <= 0 \Steiner node in
y_Vps51 - x_Vps54_Vps51 <= 0 \Steiner node in
y_Vps52 - x_Vps54_Vps52 <= 0 \Steiner node in
y_Cot1 - x_Zrc1_Cot1 <= 0 \Steiner node in
y_Rtn1 - x_Rtn1_Vps54 >= 0 \select node if connected (outgoing)
y_Vps54 - x_Vps54_Vps51 >= 0 \select node if connected (outgoing)
y_Vps54 - x_Vps54_Vps52 >= 0 \select node if connected (outgoing)
y_Vps54 - x_Vps54_Vps53 >= 0 \select node if connected (outgoing)
y_Vps51 - x_Vps51_Vps53 >= 0 \select node if connected (outgoing)
y_Vps52 - x_Vps52_Vps53 >= 0 \select node if connected (outgoing)
y_Zrc1 - x_Zrc1_Cot1 >= 0 \select node if connected (outgoing)
y_Cot1 - x_Cot1_Fmp42 >= 0 \select node if connected (outgoing)
y_Vps54 - x_Rtn1_Vps54 >= 0 \select node if connected (incoming)
y_Vps51 - x_Vps54_Vps51 >= 0 \select node if connected (incoming)
y_Vps52 - x_Vps54_Vps52 >= 0 \select node if connected (incoming)
y_Cot1 - x_Zrc1_Cot1 >= 0 \select node if connected (incoming)
y_Vps53 - x_Vps54_Vps53 >= 0 \select node if connected (incoming)
y_Vps53 - x_Vps51_Vps53 >= 0 \select node if connected (incoming)
y_Vps53 - x_Vps52_Vps53 >= 0 \select node if connected (incoming)
y_Fmp42 - x_Cot1_Fmp42 >= 0 \select node if connected (incoming)

binary
y_Rtn1
y_Vps54
y_Vps51
y_Vps52
y_Zrc1
y_Cot1
y_Vps53
y_Fmp42
x_Rtn1_Vps54
x_Vps54_Vps51
x_Vps54_Vps52
x_Vps54_Vps53
x_Vps51_Vps53
x_Vps52_Vps53
x_Zrc1_Cot1
x_Cot1_Fmp42

End