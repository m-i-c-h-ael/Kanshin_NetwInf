\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 0.214 x_Bul2_Rsp5 + 0.173 x_Rsp5_Rpn10 + 0.203 x_Rpn10_Sec18 + 0.193 x_Rsp5_Yhr131c - 2.51030846041447 y_Bul2 - 0 y_Rsp5 - 0 y_Rpn10 - 1.00307268675229 y_Sec18 - 2.33925737568746 y_Yhr131c

subject to
x_Bul2_Rsp5 >= 1 \edges from terminal
x_Rpn10_Sec18 >= 1 \edges to terminal
x_Rsp5_Yhr131c >= 1 \edges to terminal
y_Bul2 >= 1  \terminal
y_Sec18 >= 1  \terminal
y_Yhr131c >= 1  \terminal
y_Rsp5 - x_Rsp5_Rpn10 - x_Rsp5_Yhr131c <= 0 \Steiner node out
y_Rpn10 - x_Rpn10_Sec18 <= 0 \Steiner node out
y_Rsp5 - x_Bul2_Rsp5 <= 0 \Steiner node in
y_Rpn10 - x_Rsp5_Rpn10 <= 0 \Steiner node in
y_Bul2 - x_Bul2_Rsp5 >= 0 \select node if connected (outgoing)
y_Rsp5 - x_Rsp5_Rpn10 >= 0 \select node if connected (outgoing)
y_Rsp5 - x_Rsp5_Yhr131c >= 0 \select node if connected (outgoing)
y_Rpn10 - x_Rpn10_Sec18 >= 0 \select node if connected (outgoing)
y_Rsp5 - x_Bul2_Rsp5 >= 0 \select node if connected (incoming)
y_Rpn10 - x_Rsp5_Rpn10 >= 0 \select node if connected (incoming)
y_Sec18 - x_Rpn10_Sec18 >= 0 \select node if connected (incoming)
y_Yhr131c - x_Rsp5_Yhr131c >= 0 \select node if connected (incoming)

binary
y_Bul2
y_Rsp5
y_Rpn10
y_Sec18
y_Yhr131c
x_Bul2_Rsp5
x_Rsp5_Rpn10
x_Rsp5_Yhr131c
x_Rpn10_Sec18

End