\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 0.26 x_Akl1_Apm4 + 0.00900000000000001 x_Apm4_Apl1 + 0.297 x_Apl1_Chc1 + 0.037 x_Apm4_Apl2 + 0.081 x_Apl2_Chc1 + 0.005 x_Apm4_Apl3 + 0.037 x_Apl3_Sec16 + 0.038 x_Apm4_Sec16 + 0.209 x_Sec18_Rpn1 + 0.052 x_Rpn1_Rps31 + 0.2 x_Rps31_Sla1 + 0.00800000000000001 x_Rpn1_Ubi4 + 0.052 x_Ubi4_Sla1 + 0.203 x_Sec18_Rpn10 + 0.005 x_Rpn10_Rps31 + 0.173 x_Rpn10_Rsp5 + 0.181 x_Rsp5_Art5 + 0.221 x_Rsp5_Ecm21 + 0.093 x_Rsp5_Ldb19 + 0.001 x_Rpn10_Ubi4 + 0.22 x_Sec18_Rpn11 + 0.052 x_Rpn11_Rps31 + 0.00900000000000001 x_Rpn11_Ubi4 + 0.00700000000000001 x_Sec18_Sec17 + 0.038 x_Sec17_Sar1 + 0.281 x_Sar1_Sec16 + 0.099 x_Sec18_Sly1 + 0.29 x_Sly1_Sec23 + 0.0679999999999999 x_Sec23_Sec16 + 0.04 x_Ssn2_Med6 + 0.069 x_Med6_Rpb10 + 0.2 x_Rpb10_Bdp1 + 0.0649999999999999 x_Med6_Rpb5 + 0.2 x_Rpb5_Bdp1 + 0.00900000000000001 x_Med6_Rpb8 + 0.2 x_Rpb8_Bdp1 + 0.00900000000000001 x_Med6_Rpc10 + 0.2 x_Rpc10_Bdp1 + 0.069 x_Med6_Rpo26 + 0.196 x_Rpo26_Bdp1 + 0.2 x_Med6_Spt15 + 0.00900000000000001 x_Spt15_Bdp1 + 0.056 x_Ssn2_Med7 + 0.00800000000000001 x_Med7_Rpb10 + 0.00800000000000001 x_Med7_Rpb5 + 0.0580000000000001 x_Med7_Rpb8 + 0.069 x_Med7_Rpc10 + 0.00800000000000001 x_Med7_Rpo26 + 0.2 x_Med7_Spt15 + 0.022 x_Ssn2_Rgr1 + 0.069 x_Rgr1_Rpb10 + 0.048 x_Rgr1_Rpb5 + 0.069 x_Rgr1_Rpb8 + 0.069 x_Rgr1_Rpc10 + 0.049 x_Rgr1_Rpo26 + 0.2 x_Rgr1_Spt15 + 0.249 x_Ssn2_Rox3 + 0.2 x_Rox3_Rpb10 + 0.2 x_Rox3_Rpb5 + 0.2 x_Rox3_Rpb8 + 0.2 x_Rox3_Rpc10 + 0.2 x_Rox3_Rpo26 + 0.2 x_Rox3_Spt15 + 0.253 x_Ssn2_Soh1 + 0.084 x_Soh1_Rpb10 + 0.0570000000000001 x_Soh1_Rpb5 + 0.0679999999999999 x_Soh1_Rpb8 + 0.084 x_Soh1_Rpc10 + 0.084 x_Soh1_Rpo26 + 0.2 x_Soh1_Spt15 + 0.227 x_Ssn2_Srb4 + 0.0600000000000001 x_Srb4_Rpb10 + 0.0600000000000001 x_Srb4_Rpb5 + 0.0600000000000001 x_Srb4_Rpb8 + 0.0600000000000001 x_Srb4_Rpc10 + 0.0600000000000001 x_Srb4_Rpo26 + 0.112 x_Srb4_Spt15 + 0.219 x_Ssn2_Srb5 + 0.2 x_Srb5_Rpb10 + 0.143 x_Srb5_Rpb5 + 0.143 x_Srb5_Rpb8 + 0.143 x_Srb5_Rpc10 + 0.143 x_Srb5_Rpo26 + 0.112 x_Srb5_Spt15 + 0.233 x_Ssn2_Srb6 + 0.2 x_Srb6_Rpb10 + 0.2 x_Srb6_Rpb5 + 0.179 x_Srb6_Rpb8 + 0.2 x_Srb6_Rpc10 + 0.2 x_Srb6_Rpo26 + 0.2 x_Srb6_Spt15 + 0.249 x_Ssn2_Srb7 + 0.077 x_Srb7_Rpb10 + 0.0620000000000001 x_Srb7_Rpb5 + 0.077 x_Srb7_Rpb8 + 0.077 x_Srb7_Rpc10 + 0.077 x_Srb7_Rpo26 + 0.2 x_Srb7_Spt15 + 0.193 x_Yhr131c_Rsp5 + 0.00600000000000001 x_Rsp5_Rps31 + 0.195 x_Rsp5_Rvs167 + 0.229 x_Rvs167_Sla1 + 0.013 x_Rsp5_Ubi4 - 1.20478731460912 y_Akl1 - 0 y_Apm4 - 0.368956473565731 y_Apl1 - 0.428048790864964 y_Apl2 - 0.425537292960004 y_Apl3 - 1.00307268675229 y_Sec18 - 0.134940759104828 y_Rpn1 - 0.227967035078496 y_Rps31 - 0 y_Ubi4 - 0 y_Rpn10 - 0 y_Rsp5 - 0 y_Rpn11 - 0 y_Sec17 - 0 y_Sar1 - 0 y_Sly1 - 0.482891555889656 y_Sec23 - 1.39922909531988 y_Ssn2 - 0.106871738557176 y_Med6 - 0 y_Rpb10 - 0 y_Rpb5 - 0 y_Rpb8 - 0 y_Rpc10 - 0 y_Rpo26 - 0 y_Spt15 - 0 y_Med7 - 0.177114703895059 y_Rgr1 - 0 y_Rox3 - 0 y_Soh1 - 0 y_Srb4 - 0 y_Srb5 - 0 y_Srb6 - 0 y_Srb7 - 2.33925737568746 y_Yhr131c - 0.2916138569775 y_Rvs167 - 1.02568176866035 y_Chc1 - 1.36340347864178 y_Sec16 - 1.22235198568956 y_Sla1 - 1.85306782505179 y_Art5 - 1.69786774814898 y_Ecm21 - 1.71435498238203 y_Ldb19 - 1.36560429096284 y_Bdp1

subject to
x_Akl1_Apm4 >= 1 \edges from terminal
x_Sec18_Rpn1 + x_Sec18_Rpn10 + x_Sec18_Rpn11 + x_Sec18_Sec17 + x_Sec18_Sly1 >= 1 \edges from terminal
x_Ssn2_Med6 + x_Ssn2_Med7 + x_Ssn2_Rgr1 + x_Ssn2_Rox3 + x_Ssn2_Soh1 + x_Ssn2_Srb4 + x_Ssn2_Srb5 + x_Ssn2_Srb6 + x_Ssn2_Srb7 >= 1 \edges from terminal
x_Yhr131c_Rsp5 >= 1 \edges from terminal
x_Apl1_Chc1 + x_Apl2_Chc1 >= 1 \edges to terminal
x_Apm4_Sec16 + x_Apl3_Sec16 + x_Sar1_Sec16 + x_Sec23_Sec16 >= 1 \edges to terminal
x_Rps31_Sla1 + x_Ubi4_Sla1 + x_Rvs167_Sla1 >= 1 \edges to terminal
x_Rsp5_Art5 >= 1 \edges to terminal
x_Rsp5_Ecm21 >= 1 \edges to terminal
x_Rsp5_Ldb19 >= 1 \edges to terminal
x_Rpb10_Bdp1 + x_Rpb5_Bdp1 + x_Rpb8_Bdp1 + x_Rpc10_Bdp1 + x_Rpo26_Bdp1 + x_Spt15_Bdp1 >= 1 \edges to terminal
y_Akl1 >= 1  \terminal
y_Sec18 >= 1  \terminal
y_Ssn2 >= 1  \terminal
y_Yhr131c >= 1  \terminal
y_Chc1 >= 1  \terminal
y_Sec16 >= 1  \terminal
y_Sla1 >= 1  \terminal
y_Art5 >= 1  \terminal
y_Ecm21 >= 1  \terminal
y_Ldb19 >= 1  \terminal
y_Bdp1 >= 1  \terminal
y_Apm4 - x_Apm4_Apl1 - x_Apm4_Apl2 - x_Apm4_Apl3 - x_Apm4_Sec16 <= 0 \Steiner node out
y_Apl1 - x_Apl1_Chc1 <= 0 \Steiner node out
y_Apl2 - x_Apl2_Chc1 <= 0 \Steiner node out
y_Apl3 - x_Apl3_Sec16 <= 0 \Steiner node out
y_Rpn1 - x_Rpn1_Rps31 - x_Rpn1_Ubi4 <= 0 \Steiner node out
y_Rps31 - x_Rps31_Sla1 <= 0 \Steiner node out
y_Ubi4 - x_Ubi4_Sla1 <= 0 \Steiner node out
y_Rpn10 - x_Rpn10_Rps31 - x_Rpn10_Ubi4 - x_Rpn10_Rsp5 <= 0 \Steiner node out
y_Rsp5 - x_Rsp5_Rps31 - x_Rsp5_Ubi4 - x_Rsp5_Rvs167 - x_Rsp5_Art5 - x_Rsp5_Ecm21 - x_Rsp5_Ldb19 <= 0 \Steiner node out
y_Rpn11 - x_Rpn11_Rps31 - x_Rpn11_Ubi4 <= 0 \Steiner node out
y_Sec17 - x_Sec17_Sar1 <= 0 \Steiner node out
y_Sar1 - x_Sar1_Sec16 <= 0 \Steiner node out
y_Sly1 - x_Sly1_Sec23 <= 0 \Steiner node out
y_Sec23 - x_Sec23_Sec16 <= 0 \Steiner node out
y_Med6 - x_Med6_Rpb10 - x_Med6_Rpb5 - x_Med6_Rpb8 - x_Med6_Rpc10 - x_Med6_Rpo26 - x_Med6_Spt15 <= 0 \Steiner node out
y_Rpb10 - x_Rpb10_Bdp1 <= 0 \Steiner node out
y_Rpb5 - x_Rpb5_Bdp1 <= 0 \Steiner node out
y_Rpb8 - x_Rpb8_Bdp1 <= 0 \Steiner node out
y_Rpc10 - x_Rpc10_Bdp1 <= 0 \Steiner node out
y_Rpo26 - x_Rpo26_Bdp1 <= 0 \Steiner node out
y_Spt15 - x_Spt15_Bdp1 <= 0 \Steiner node out
y_Med7 - x_Med7_Rpb10 - x_Med7_Rpb5 - x_Med7_Rpb8 - x_Med7_Rpc10 - x_Med7_Rpo26 - x_Med7_Spt15 <= 0 \Steiner node out
y_Rgr1 - x_Rgr1_Rpb10 - x_Rgr1_Rpb5 - x_Rgr1_Rpb8 - x_Rgr1_Rpc10 - x_Rgr1_Rpo26 - x_Rgr1_Spt15 <= 0 \Steiner node out
y_Rox3 - x_Rox3_Rpb10 - x_Rox3_Rpb5 - x_Rox3_Rpb8 - x_Rox3_Rpc10 - x_Rox3_Rpo26 - x_Rox3_Spt15 <= 0 \Steiner node out
y_Soh1 - x_Soh1_Rpb10 - x_Soh1_Rpb5 - x_Soh1_Rpb8 - x_Soh1_Rpc10 - x_Soh1_Rpo26 - x_Soh1_Spt15 <= 0 \Steiner node out
y_Srb4 - x_Srb4_Rpb10 - x_Srb4_Rpb5 - x_Srb4_Rpb8 - x_Srb4_Rpc10 - x_Srb4_Rpo26 - x_Srb4_Spt15 <= 0 \Steiner node out
y_Srb5 - x_Srb5_Rpb10 - x_Srb5_Rpb5 - x_Srb5_Rpb8 - x_Srb5_Rpc10 - x_Srb5_Rpo26 - x_Srb5_Spt15 <= 0 \Steiner node out
y_Srb6 - x_Srb6_Rpb10 - x_Srb6_Rpb5 - x_Srb6_Rpb8 - x_Srb6_Rpc10 - x_Srb6_Rpo26 - x_Srb6_Spt15 <= 0 \Steiner node out
y_Srb7 - x_Srb7_Rpb10 - x_Srb7_Rpb5 - x_Srb7_Rpb8 - x_Srb7_Rpc10 - x_Srb7_Rpo26 - x_Srb7_Spt15 <= 0 \Steiner node out
y_Rvs167 - x_Rvs167_Sla1 <= 0 \Steiner node out
y_Apm4 - x_Akl1_Apm4 <= 0 \Steiner node in
y_Apl1 - x_Apm4_Apl1 <= 0 \Steiner node in
y_Apl2 - x_Apm4_Apl2 <= 0 \Steiner node in
y_Apl3 - x_Apm4_Apl3 <= 0 \Steiner node in
y_Rpn1 - x_Sec18_Rpn1 <= 0 \Steiner node in
y_Rps31 - x_Rpn1_Rps31 - x_Rpn10_Rps31 - x_Rsp5_Rps31 - x_Rpn11_Rps31 <= 0 \Steiner node in
y_Ubi4 - x_Rpn1_Ubi4 - x_Rpn10_Ubi4 - x_Rsp5_Ubi4 - x_Rpn11_Ubi4 <= 0 \Steiner node in
y_Rpn10 - x_Sec18_Rpn10 <= 0 \Steiner node in
y_Rsp5 - x_Rpn10_Rsp5 - x_Yhr131c_Rsp5 <= 0 \Steiner node in
y_Rpn11 - x_Sec18_Rpn11 <= 0 \Steiner node in
y_Sec17 - x_Sec18_Sec17 <= 0 \Steiner node in
y_Sar1 - x_Sec17_Sar1 <= 0 \Steiner node in
y_Sly1 - x_Sec18_Sly1 <= 0 \Steiner node in
y_Sec23 - x_Sly1_Sec23 <= 0 \Steiner node in
y_Med6 - x_Ssn2_Med6 <= 0 \Steiner node in
y_Rpb10 - x_Med6_Rpb10 - x_Med7_Rpb10 - x_Rgr1_Rpb10 - x_Rox3_Rpb10 - x_Soh1_Rpb10 - x_Srb4_Rpb10 - x_Srb5_Rpb10 - x_Srb6_Rpb10 - x_Srb7_Rpb10 <= 0 \Steiner node in
y_Rpb5 - x_Med6_Rpb5 - x_Med7_Rpb5 - x_Rgr1_Rpb5 - x_Rox3_Rpb5 - x_Soh1_Rpb5 - x_Srb4_Rpb5 - x_Srb5_Rpb5 - x_Srb6_Rpb5 - x_Srb7_Rpb5 <= 0 \Steiner node in
y_Rpb8 - x_Med6_Rpb8 - x_Med7_Rpb8 - x_Rgr1_Rpb8 - x_Rox3_Rpb8 - x_Soh1_Rpb8 - x_Srb4_Rpb8 - x_Srb5_Rpb8 - x_Srb6_Rpb8 - x_Srb7_Rpb8 <= 0 \Steiner node in
y_Rpc10 - x_Med6_Rpc10 - x_Med7_Rpc10 - x_Rgr1_Rpc10 - x_Rox3_Rpc10 - x_Soh1_Rpc10 - x_Srb4_Rpc10 - x_Srb5_Rpc10 - x_Srb6_Rpc10 - x_Srb7_Rpc10 <= 0 \Steiner node in
y_Rpo26 - x_Med6_Rpo26 - x_Med7_Rpo26 - x_Rgr1_Rpo26 - x_Rox3_Rpo26 - x_Soh1_Rpo26 - x_Srb4_Rpo26 - x_Srb5_Rpo26 - x_Srb6_Rpo26 - x_Srb7_Rpo26 <= 0 \Steiner node in
y_Spt15 - x_Med6_Spt15 - x_Med7_Spt15 - x_Rgr1_Spt15 - x_Rox3_Spt15 - x_Soh1_Spt15 - x_Srb4_Spt15 - x_Srb5_Spt15 - x_Srb6_Spt15 - x_Srb7_Spt15 <= 0 \Steiner node in
y_Med7 - x_Ssn2_Med7 <= 0 \Steiner node in
y_Rgr1 - x_Ssn2_Rgr1 <= 0 \Steiner node in
y_Rox3 - x_Ssn2_Rox3 <= 0 \Steiner node in
y_Soh1 - x_Ssn2_Soh1 <= 0 \Steiner node in
y_Srb4 - x_Ssn2_Srb4 <= 0 \Steiner node in
y_Srb5 - x_Ssn2_Srb5 <= 0 \Steiner node in
y_Srb6 - x_Ssn2_Srb6 <= 0 \Steiner node in
y_Srb7 - x_Ssn2_Srb7 <= 0 \Steiner node in
y_Rvs167 - x_Rsp5_Rvs167 <= 0 \Steiner node in
y_Akl1 - x_Akl1_Apm4 >= 0 \select node if connected (outgoing)
y_Apm4 - x_Apm4_Apl1 >= 0 \select node if connected (outgoing)
y_Apm4 - x_Apm4_Apl2 >= 0 \select node if connected (outgoing)
y_Apm4 - x_Apm4_Apl3 >= 0 \select node if connected (outgoing)
y_Apm4 - x_Apm4_Sec16 >= 0 \select node if connected (outgoing)
y_Apl1 - x_Apl1_Chc1 >= 0 \select node if connected (outgoing)
y_Apl2 - x_Apl2_Chc1 >= 0 \select node if connected (outgoing)
y_Apl3 - x_Apl3_Sec16 >= 0 \select node if connected (outgoing)
y_Sec18 - x_Sec18_Rpn1 >= 0 \select node if connected (outgoing)
y_Sec18 - x_Sec18_Rpn10 >= 0 \select node if connected (outgoing)
y_Sec18 - x_Sec18_Rpn11 >= 0 \select node if connected (outgoing)
y_Sec18 - x_Sec18_Sec17 >= 0 \select node if connected (outgoing)
y_Sec18 - x_Sec18_Sly1 >= 0 \select node if connected (outgoing)
y_Rpn1 - x_Rpn1_Rps31 >= 0 \select node if connected (outgoing)
y_Rpn1 - x_Rpn1_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Sla1 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Sla1 >= 0 \select node if connected (outgoing)
y_Rpn10 - x_Rpn10_Rps31 >= 0 \select node if connected (outgoing)
y_Rpn10 - x_Rpn10_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpn10 - x_Rpn10_Rsp5 >= 0 \select node if connected (outgoing)
y_Rsp5 - x_Rsp5_Rps31 >= 0 \select node if connected (outgoing)
y_Rsp5 - x_Rsp5_Ubi4 >= 0 \select node if connected (outgoing)
y_Rsp5 - x_Rsp5_Rvs167 >= 0 \select node if connected (outgoing)
y_Rsp5 - x_Rsp5_Art5 >= 0 \select node if connected (outgoing)
y_Rsp5 - x_Rsp5_Ecm21 >= 0 \select node if connected (outgoing)
y_Rsp5 - x_Rsp5_Ldb19 >= 0 \select node if connected (outgoing)
y_Rpn11 - x_Rpn11_Rps31 >= 0 \select node if connected (outgoing)
y_Rpn11 - x_Rpn11_Ubi4 >= 0 \select node if connected (outgoing)
y_Sec17 - x_Sec17_Sar1 >= 0 \select node if connected (outgoing)
y_Sar1 - x_Sar1_Sec16 >= 0 \select node if connected (outgoing)
y_Sly1 - x_Sly1_Sec23 >= 0 \select node if connected (outgoing)
y_Sec23 - x_Sec23_Sec16 >= 0 \select node if connected (outgoing)
y_Ssn2 - x_Ssn2_Med6 >= 0 \select node if connected (outgoing)
y_Ssn2 - x_Ssn2_Med7 >= 0 \select node if connected (outgoing)
y_Ssn2 - x_Ssn2_Rgr1 >= 0 \select node if connected (outgoing)
y_Ssn2 - x_Ssn2_Rox3 >= 0 \select node if connected (outgoing)
y_Ssn2 - x_Ssn2_Soh1 >= 0 \select node if connected (outgoing)
y_Ssn2 - x_Ssn2_Srb4 >= 0 \select node if connected (outgoing)
y_Ssn2 - x_Ssn2_Srb5 >= 0 \select node if connected (outgoing)
y_Ssn2 - x_Ssn2_Srb6 >= 0 \select node if connected (outgoing)
y_Ssn2 - x_Ssn2_Srb7 >= 0 \select node if connected (outgoing)
y_Med6 - x_Med6_Rpb10 >= 0 \select node if connected (outgoing)
y_Med6 - x_Med6_Rpb5 >= 0 \select node if connected (outgoing)
y_Med6 - x_Med6_Rpb8 >= 0 \select node if connected (outgoing)
y_Med6 - x_Med6_Rpc10 >= 0 \select node if connected (outgoing)
y_Med6 - x_Med6_Rpo26 >= 0 \select node if connected (outgoing)
y_Med6 - x_Med6_Spt15 >= 0 \select node if connected (outgoing)
y_Rpb10 - x_Rpb10_Bdp1 >= 0 \select node if connected (outgoing)
y_Rpb5 - x_Rpb5_Bdp1 >= 0 \select node if connected (outgoing)
y_Rpb8 - x_Rpb8_Bdp1 >= 0 \select node if connected (outgoing)
y_Rpc10 - x_Rpc10_Bdp1 >= 0 \select node if connected (outgoing)
y_Rpo26 - x_Rpo26_Bdp1 >= 0 \select node if connected (outgoing)
y_Spt15 - x_Spt15_Bdp1 >= 0 \select node if connected (outgoing)
y_Med7 - x_Med7_Rpb10 >= 0 \select node if connected (outgoing)
y_Med7 - x_Med7_Rpb5 >= 0 \select node if connected (outgoing)
y_Med7 - x_Med7_Rpb8 >= 0 \select node if connected (outgoing)
y_Med7 - x_Med7_Rpc10 >= 0 \select node if connected (outgoing)
y_Med7 - x_Med7_Rpo26 >= 0 \select node if connected (outgoing)
y_Med7 - x_Med7_Spt15 >= 0 \select node if connected (outgoing)
y_Rgr1 - x_Rgr1_Rpb10 >= 0 \select node if connected (outgoing)
y_Rgr1 - x_Rgr1_Rpb5 >= 0 \select node if connected (outgoing)
y_Rgr1 - x_Rgr1_Rpb8 >= 0 \select node if connected (outgoing)
y_Rgr1 - x_Rgr1_Rpc10 >= 0 \select node if connected (outgoing)
y_Rgr1 - x_Rgr1_Rpo26 >= 0 \select node if connected (outgoing)
y_Rgr1 - x_Rgr1_Spt15 >= 0 \select node if connected (outgoing)
y_Rox3 - x_Rox3_Rpb10 >= 0 \select node if connected (outgoing)
y_Rox3 - x_Rox3_Rpb5 >= 0 \select node if connected (outgoing)
y_Rox3 - x_Rox3_Rpb8 >= 0 \select node if connected (outgoing)
y_Rox3 - x_Rox3_Rpc10 >= 0 \select node if connected (outgoing)
y_Rox3 - x_Rox3_Rpo26 >= 0 \select node if connected (outgoing)
y_Rox3 - x_Rox3_Spt15 >= 0 \select node if connected (outgoing)
y_Soh1 - x_Soh1_Rpb10 >= 0 \select node if connected (outgoing)
y_Soh1 - x_Soh1_Rpb5 >= 0 \select node if connected (outgoing)
y_Soh1 - x_Soh1_Rpb8 >= 0 \select node if connected (outgoing)
y_Soh1 - x_Soh1_Rpc10 >= 0 \select node if connected (outgoing)
y_Soh1 - x_Soh1_Rpo26 >= 0 \select node if connected (outgoing)
y_Soh1 - x_Soh1_Spt15 >= 0 \select node if connected (outgoing)
y_Srb4 - x_Srb4_Rpb10 >= 0 \select node if connected (outgoing)
y_Srb4 - x_Srb4_Rpb5 >= 0 \select node if connected (outgoing)
y_Srb4 - x_Srb4_Rpb8 >= 0 \select node if connected (outgoing)
y_Srb4 - x_Srb4_Rpc10 >= 0 \select node if connected (outgoing)
y_Srb4 - x_Srb4_Rpo26 >= 0 \select node if connected (outgoing)
y_Srb4 - x_Srb4_Spt15 >= 0 \select node if connected (outgoing)
y_Srb5 - x_Srb5_Rpb10 >= 0 \select node if connected (outgoing)
y_Srb5 - x_Srb5_Rpb5 >= 0 \select node if connected (outgoing)
y_Srb5 - x_Srb5_Rpb8 >= 0 \select node if connected (outgoing)
y_Srb5 - x_Srb5_Rpc10 >= 0 \select node if connected (outgoing)
y_Srb5 - x_Srb5_Rpo26 >= 0 \select node if connected (outgoing)
y_Srb5 - x_Srb5_Spt15 >= 0 \select node if connected (outgoing)
y_Srb6 - x_Srb6_Rpb10 >= 0 \select node if connected (outgoing)
y_Srb6 - x_Srb6_Rpb5 >= 0 \select node if connected (outgoing)
y_Srb6 - x_Srb6_Rpb8 >= 0 \select node if connected (outgoing)
y_Srb6 - x_Srb6_Rpc10 >= 0 \select node if connected (outgoing)
y_Srb6 - x_Srb6_Rpo26 >= 0 \select node if connected (outgoing)
y_Srb6 - x_Srb6_Spt15 >= 0 \select node if connected (outgoing)
y_Srb7 - x_Srb7_Rpb10 >= 0 \select node if connected (outgoing)
y_Srb7 - x_Srb7_Rpb5 >= 0 \select node if connected (outgoing)
y_Srb7 - x_Srb7_Rpb8 >= 0 \select node if connected (outgoing)
y_Srb7 - x_Srb7_Rpc10 >= 0 \select node if connected (outgoing)
y_Srb7 - x_Srb7_Rpo26 >= 0 \select node if connected (outgoing)
y_Srb7 - x_Srb7_Spt15 >= 0 \select node if connected (outgoing)
y_Yhr131c - x_Yhr131c_Rsp5 >= 0 \select node if connected (outgoing)
y_Rvs167 - x_Rvs167_Sla1 >= 0 \select node if connected (outgoing)
y_Apm4 - x_Akl1_Apm4 >= 0 \select node if connected (incoming)
y_Apl1 - x_Apm4_Apl1 >= 0 \select node if connected (incoming)
y_Apl2 - x_Apm4_Apl2 >= 0 \select node if connected (incoming)
y_Apl3 - x_Apm4_Apl3 >= 0 \select node if connected (incoming)
y_Rpn1 - x_Sec18_Rpn1 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpn1_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpn10_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rsp5_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpn11_Rps31 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpn1_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpn10_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rsp5_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpn11_Ubi4 >= 0 \select node if connected (incoming)
y_Rpn10 - x_Sec18_Rpn10 >= 0 \select node if connected (incoming)
y_Rsp5 - x_Rpn10_Rsp5 >= 0 \select node if connected (incoming)
y_Rsp5 - x_Yhr131c_Rsp5 >= 0 \select node if connected (incoming)
y_Rpn11 - x_Sec18_Rpn11 >= 0 \select node if connected (incoming)
y_Sec17 - x_Sec18_Sec17 >= 0 \select node if connected (incoming)
y_Sar1 - x_Sec17_Sar1 >= 0 \select node if connected (incoming)
y_Sly1 - x_Sec18_Sly1 >= 0 \select node if connected (incoming)
y_Sec23 - x_Sly1_Sec23 >= 0 \select node if connected (incoming)
y_Med6 - x_Ssn2_Med6 >= 0 \select node if connected (incoming)
y_Rpb10 - x_Med6_Rpb10 >= 0 \select node if connected (incoming)
y_Rpb10 - x_Med7_Rpb10 >= 0 \select node if connected (incoming)
y_Rpb10 - x_Rgr1_Rpb10 >= 0 \select node if connected (incoming)
y_Rpb10 - x_Rox3_Rpb10 >= 0 \select node if connected (incoming)
y_Rpb10 - x_Soh1_Rpb10 >= 0 \select node if connected (incoming)
y_Rpb10 - x_Srb4_Rpb10 >= 0 \select node if connected (incoming)
y_Rpb10 - x_Srb5_Rpb10 >= 0 \select node if connected (incoming)
y_Rpb10 - x_Srb6_Rpb10 >= 0 \select node if connected (incoming)
y_Rpb10 - x_Srb7_Rpb10 >= 0 \select node if connected (incoming)
y_Rpb5 - x_Med6_Rpb5 >= 0 \select node if connected (incoming)
y_Rpb5 - x_Med7_Rpb5 >= 0 \select node if connected (incoming)
y_Rpb5 - x_Rgr1_Rpb5 >= 0 \select node if connected (incoming)
y_Rpb5 - x_Rox3_Rpb5 >= 0 \select node if connected (incoming)
y_Rpb5 - x_Soh1_Rpb5 >= 0 \select node if connected (incoming)
y_Rpb5 - x_Srb4_Rpb5 >= 0 \select node if connected (incoming)
y_Rpb5 - x_Srb5_Rpb5 >= 0 \select node if connected (incoming)
y_Rpb5 - x_Srb6_Rpb5 >= 0 \select node if connected (incoming)
y_Rpb5 - x_Srb7_Rpb5 >= 0 \select node if connected (incoming)
y_Rpb8 - x_Med6_Rpb8 >= 0 \select node if connected (incoming)
y_Rpb8 - x_Med7_Rpb8 >= 0 \select node if connected (incoming)
y_Rpb8 - x_Rgr1_Rpb8 >= 0 \select node if connected (incoming)
y_Rpb8 - x_Rox3_Rpb8 >= 0 \select node if connected (incoming)
y_Rpb8 - x_Soh1_Rpb8 >= 0 \select node if connected (incoming)
y_Rpb8 - x_Srb4_Rpb8 >= 0 \select node if connected (incoming)
y_Rpb8 - x_Srb5_Rpb8 >= 0 \select node if connected (incoming)
y_Rpb8 - x_Srb6_Rpb8 >= 0 \select node if connected (incoming)
y_Rpb8 - x_Srb7_Rpb8 >= 0 \select node if connected (incoming)
y_Rpc10 - x_Med6_Rpc10 >= 0 \select node if connected (incoming)
y_Rpc10 - x_Med7_Rpc10 >= 0 \select node if connected (incoming)
y_Rpc10 - x_Rgr1_Rpc10 >= 0 \select node if connected (incoming)
y_Rpc10 - x_Rox3_Rpc10 >= 0 \select node if connected (incoming)
y_Rpc10 - x_Soh1_Rpc10 >= 0 \select node if connected (incoming)
y_Rpc10 - x_Srb4_Rpc10 >= 0 \select node if connected (incoming)
y_Rpc10 - x_Srb5_Rpc10 >= 0 \select node if connected (incoming)
y_Rpc10 - x_Srb6_Rpc10 >= 0 \select node if connected (incoming)
y_Rpc10 - x_Srb7_Rpc10 >= 0 \select node if connected (incoming)
y_Rpo26 - x_Med6_Rpo26 >= 0 \select node if connected (incoming)
y_Rpo26 - x_Med7_Rpo26 >= 0 \select node if connected (incoming)
y_Rpo26 - x_Rgr1_Rpo26 >= 0 \select node if connected (incoming)
y_Rpo26 - x_Rox3_Rpo26 >= 0 \select node if connected (incoming)
y_Rpo26 - x_Soh1_Rpo26 >= 0 \select node if connected (incoming)
y_Rpo26 - x_Srb4_Rpo26 >= 0 \select node if connected (incoming)
y_Rpo26 - x_Srb5_Rpo26 >= 0 \select node if connected (incoming)
y_Rpo26 - x_Srb6_Rpo26 >= 0 \select node if connected (incoming)
y_Rpo26 - x_Srb7_Rpo26 >= 0 \select node if connected (incoming)
y_Spt15 - x_Med6_Spt15 >= 0 \select node if connected (incoming)
y_Spt15 - x_Med7_Spt15 >= 0 \select node if connected (incoming)
y_Spt15 - x_Rgr1_Spt15 >= 0 \select node if connected (incoming)
y_Spt15 - x_Rox3_Spt15 >= 0 \select node if connected (incoming)
y_Spt15 - x_Soh1_Spt15 >= 0 \select node if connected (incoming)
y_Spt15 - x_Srb4_Spt15 >= 0 \select node if connected (incoming)
y_Spt15 - x_Srb5_Spt15 >= 0 \select node if connected (incoming)
y_Spt15 - x_Srb6_Spt15 >= 0 \select node if connected (incoming)
y_Spt15 - x_Srb7_Spt15 >= 0 \select node if connected (incoming)
y_Med7 - x_Ssn2_Med7 >= 0 \select node if connected (incoming)
y_Rgr1 - x_Ssn2_Rgr1 >= 0 \select node if connected (incoming)
y_Rox3 - x_Ssn2_Rox3 >= 0 \select node if connected (incoming)
y_Soh1 - x_Ssn2_Soh1 >= 0 \select node if connected (incoming)
y_Srb4 - x_Ssn2_Srb4 >= 0 \select node if connected (incoming)
y_Srb5 - x_Ssn2_Srb5 >= 0 \select node if connected (incoming)
y_Srb6 - x_Ssn2_Srb6 >= 0 \select node if connected (incoming)
y_Srb7 - x_Ssn2_Srb7 >= 0 \select node if connected (incoming)
y_Rvs167 - x_Rsp5_Rvs167 >= 0 \select node if connected (incoming)
y_Chc1 - x_Apl1_Chc1 >= 0 \select node if connected (incoming)
y_Chc1 - x_Apl2_Chc1 >= 0 \select node if connected (incoming)
y_Sec16 - x_Apm4_Sec16 >= 0 \select node if connected (incoming)
y_Sec16 - x_Apl3_Sec16 >= 0 \select node if connected (incoming)
y_Sec16 - x_Sar1_Sec16 >= 0 \select node if connected (incoming)
y_Sec16 - x_Sec23_Sec16 >= 0 \select node if connected (incoming)
y_Sla1 - x_Rps31_Sla1 >= 0 \select node if connected (incoming)
y_Sla1 - x_Ubi4_Sla1 >= 0 \select node if connected (incoming)
y_Sla1 - x_Rvs167_Sla1 >= 0 \select node if connected (incoming)
y_Art5 - x_Rsp5_Art5 >= 0 \select node if connected (incoming)
y_Ecm21 - x_Rsp5_Ecm21 >= 0 \select node if connected (incoming)
y_Ldb19 - x_Rsp5_Ldb19 >= 0 \select node if connected (incoming)
y_Bdp1 - x_Rpb10_Bdp1 >= 0 \select node if connected (incoming)
y_Bdp1 - x_Rpb5_Bdp1 >= 0 \select node if connected (incoming)
y_Bdp1 - x_Rpb8_Bdp1 >= 0 \select node if connected (incoming)
y_Bdp1 - x_Rpc10_Bdp1 >= 0 \select node if connected (incoming)
y_Bdp1 - x_Rpo26_Bdp1 >= 0 \select node if connected (incoming)
y_Bdp1 - x_Spt15_Bdp1 >= 0 \select node if connected (incoming)

binary
y_Akl1
y_Apm4
y_Apl1
y_Apl2
y_Apl3
y_Sec18
y_Rpn1
y_Rps31
y_Ubi4
y_Rpn10
y_Rsp5
y_Rpn11
y_Sec17
y_Sar1
y_Sly1
y_Sec23
y_Ssn2
y_Med6
y_Rpb10
y_Rpb5
y_Rpb8
y_Rpc10
y_Rpo26
y_Spt15
y_Med7
y_Rgr1
y_Rox3
y_Soh1
y_Srb4
y_Srb5
y_Srb6
y_Srb7
y_Yhr131c
y_Rvs167
y_Chc1
y_Sec16
y_Sla1
y_Art5
y_Ecm21
y_Ldb19
y_Bdp1
x_Akl1_Apm4
x_Apm4_Apl1
x_Apm4_Apl2
x_Apm4_Apl3
x_Apm4_Sec16
x_Apl1_Chc1
x_Apl2_Chc1
x_Apl3_Sec16
x_Sec18_Rpn1
x_Sec18_Rpn10
x_Sec18_Rpn11
x_Sec18_Sec17
x_Sec18_Sly1
x_Rpn1_Rps31
x_Rpn1_Ubi4
x_Rps31_Sla1
x_Ubi4_Sla1
x_Rpn10_Rps31
x_Rpn10_Ubi4
x_Rpn10_Rsp5
x_Rsp5_Rps31
x_Rsp5_Ubi4
x_Rsp5_Rvs167
x_Rsp5_Art5
x_Rsp5_Ecm21
x_Rsp5_Ldb19
x_Rpn11_Rps31
x_Rpn11_Ubi4
x_Sec17_Sar1
x_Sar1_Sec16
x_Sly1_Sec23
x_Sec23_Sec16
x_Ssn2_Med6
x_Ssn2_Med7
x_Ssn2_Rgr1
x_Ssn2_Rox3
x_Ssn2_Soh1
x_Ssn2_Srb4
x_Ssn2_Srb5
x_Ssn2_Srb6
x_Ssn2_Srb7
x_Med6_Rpb10
x_Med6_Rpb5
x_Med6_Rpb8
x_Med6_Rpc10
x_Med6_Rpo26
x_Med6_Spt15
x_Rpb10_Bdp1
x_Rpb5_Bdp1
x_Rpb8_Bdp1
x_Rpc10_Bdp1
x_Rpo26_Bdp1
x_Spt15_Bdp1
x_Med7_Rpb10
x_Med7_Rpb5
x_Med7_Rpb8
x_Med7_Rpc10
x_Med7_Rpo26
x_Med7_Spt15
x_Rgr1_Rpb10
x_Rgr1_Rpb5
x_Rgr1_Rpb8
x_Rgr1_Rpc10
x_Rgr1_Rpo26
x_Rgr1_Spt15
x_Rox3_Rpb10
x_Rox3_Rpb5
x_Rox3_Rpb8
x_Rox3_Rpc10
x_Rox3_Rpo26
x_Rox3_Spt15
x_Soh1_Rpb10
x_Soh1_Rpb5
x_Soh1_Rpb8
x_Soh1_Rpc10
x_Soh1_Rpo26
x_Soh1_Spt15
x_Srb4_Rpb10
x_Srb4_Rpb5
x_Srb4_Rpb8
x_Srb4_Rpc10
x_Srb4_Rpo26
x_Srb4_Spt15
x_Srb5_Rpb10
x_Srb5_Rpb5
x_Srb5_Rpb8
x_Srb5_Rpc10
x_Srb5_Rpo26
x_Srb5_Spt15
x_Srb6_Rpb10
x_Srb6_Rpb5
x_Srb6_Rpb8
x_Srb6_Rpc10
x_Srb6_Rpo26
x_Srb6_Spt15
x_Srb7_Rpb10
x_Srb7_Rpb5
x_Srb7_Rpb8
x_Srb7_Rpc10
x_Srb7_Rpo26
x_Srb7_Spt15
x_Yhr131c_Rsp5
x_Rvs167_Sla1

End