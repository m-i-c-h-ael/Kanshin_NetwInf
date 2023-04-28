\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 0.214 x_Cyr1_Act1 + 0.004 x_Act1_Cof1 + 0.038 x_Cof1_Prm15 + 0.028 x_Cyr1_Srv2 + 0.011 x_Srv2_Cof1 + 0.131 x_Duf1_Ubi4 + 0.005 x_Ubi4_Cdc48 + 0.215 x_Cdc48_Rad23 + 0.026 x_Ubi4_Ddi1 + 0.00900000000000001 x_Ddi1_Rad23 + 0.005 x_Ubi4_Dsk2 + 0.034 x_Dsk2_Rad23 + 0.003 x_Ubi4_Rad23 + 0.00800000000000001 x_Ubi4_Rpn1 + 0.002 x_Rpn1_Rad23 + 0.001 x_Ubi4_Rpn10 + 0.003 x_Rpn10_Rad23 + 0.00900000000000001 x_Ubi4_Rpn11 + 0.0659999999999999 x_Rpn11_Rad23 + 0.002 x_Ubi4_Rps31 + 0.037 x_Rps31_Rad23 + 0.0679999999999999 x_Ubi4_Rpt1 + 0.029 x_Rpt1_Rad23 + 0.088 x_Ubi4_Rpt5 + 0.094 x_Rpt5_Rad23 + 0.118 x_Ubi4_Rpt6 + 0.071 x_Rpt6_Rad23 + 0.19 x_Ubi4_Ufd2 + 0.015 x_Ufd2_Rad23 + 0.024 x_Ede1_Ent1 + 0.056 x_Ent1_Chc1 + 0.155 x_Chc1_Ent5 + 0.177 x_Ent1_Ubi4 + 0.026 x_Ede1_Ent2 + 0.1 x_Ent2_Chc1 + 0.286 x_Ent2_Ubi4 + 0.172 x_Ede1_Rps31 + 0.169 x_Rps31_Cdc48 + 0.082 x_Rps31_Dsk2 + 0.052 x_Rps31_Rpn1 + 0.005 x_Rps31_Rpn10 + 0.052 x_Rps31_Rpn11 + 0.052 x_Rps31_Rpt1 + 0.00900000000000001 x_Rps31_Rpt5 + 0.078 x_Rps31_Rpt6 + 0.002 x_Rps31_Ubi4 + 0.028 x_Ede1_Ubi4 + 0.0659999999999999 x_Ede1_Vps27 + 0.0570000000000001 x_Vps27_Chc1 + 0.106 x_Vps27_Rps31 + 0.001 x_Vps27_Ubi4 + 0.00800000000000001 x_Fun12_Asc1 + 0.005 x_Asc1_Rps31 + 0.039 x_Asc1_Ubi4 + 0.271 x_Fun12_Gcd11 + 0.216 x_Gcd11_Rps31 + 0.00900000000000001 x_Fun12_Nip1 + 0.246 x_Nip1_Rpn11 + 0.043 x_Nip1_Rps31 + 0.00600000000000001 x_Fun12_Prt1 + 0.078 x_Prt1_Rps31 + 0.016 x_Fun12_Rpg1 + 0.039 x_Rpg1_Rps31 + 0.00800000000000001 x_Fun12_Rpl10 + 0.005 x_Rpl10_Rps31 + 0.028 x_Rpl10_Ubi4 + 0.032 x_Fun12_Rpl17a + 0.00800000000000001 x_Rpl17a_Rps31 + 0.221 x_Rpl17a_Ubi4 + 0.151 x_Fun12_Rpl17b + 0.037 x_Rpl17b_Rps31 + 0.221 x_Rpl17b_Ubi4 + 0.032 x_Fun12_Rpl19b + 0.00800000000000001 x_Rpl19b_Rps31 + 0.223 x_Rpl19b_Ubi4 + 0.00800000000000001 x_Fun12_Rpl25 + 0.005 x_Rpl25_Rps31 + 0.039 x_Rpl25_Ubi4 + 0.00800000000000001 x_Fun12_Rpl28 + 0.005 x_Rpl28_Rps31 + 0.028 x_Rpl28_Ubi4 + 0.024 x_Fun12_Rpl2a + 0.00800000000000001 x_Rpl2a_Rps31 + 0.172 x_Rpl2a_Ubi4 + 0.043 x_Fun12_Rpl2b + 0.00800000000000001 x_Rpl2b_Rps31 + 0.236 x_Rpl2b_Ubi4 + 0.13 x_Fun12_Rpl3 + 0.005 x_Rpl3_Rps31 + 0.266 x_Rpl3_Ubi4 + 0.00800000000000001 x_Fun12_Rpl30 + 0.005 x_Rpl30_Rps31 + 0.039 x_Rpl30_Ubi4 + 0.00800000000000001 x_Fun12_Rpl32 + 0.005 x_Rpl32_Rps31 + 0.028 x_Rpl32_Ubi4 + 0.032 x_Fun12_Rpl34a + 0.00800000000000001 x_Rpl34a_Rps31 + 0.16 x_Rpl34a_Ubi4 + 0.032 x_Fun12_Rpl35a + 0.00800000000000001 x_Rpl35a_Rps31 + 0.219 x_Rpl35a_Ubi4 + 0.032 x_Fun12_Rpl35b + 0.00800000000000001 x_Rpl35b_Rps31 + 0.219 x_Rpl35b_Ubi4 + 0.032 x_Fun12_Rpl37a + 0.00800000000000001 x_Rpl37a_Rps31 + 0.00800000000000001 x_Fun12_Rpl38 + 0.005 x_Rpl38_Rps31 + 0.039 x_Rpl38_Ubi4 + 0.00800000000000001 x_Fun12_Rpl5 + 0.005 x_Rpl5_Rps31 + 0.038 x_Rpl5_Ubi4 + 0.032 x_Fun12_Rpl7a + 0.00800000000000001 x_Rpl7a_Rps31 + 0.032 x_Fun12_Rps0a + 0.00800000000000001 x_Rps0a_Rps31 + 0.174 x_Rps0a_Ubi4 + 0.032 x_Fun12_Rps0b + 0.00800000000000001 x_Rps0b_Rps31 + 0.174 x_Rps0b_Ubi4 + 0.032 x_Fun12_Rps10a + 0.00800000000000001 x_Rps10a_Rps31 + 0.214 x_Rps10a_Ubi4 + 0.032 x_Fun12_Rps11a + 0.00800000000000001 x_Rps11a_Rps31 + 0.212 x_Rps11a_Ubi4 + 0.032 x_Fun12_Rps11b + 0.00800000000000001 x_Rps11b_Rps31 + 0.155 x_Rps11b_Ubi4 + 0.00900000000000001 x_Fun12_Rps12 + 0.005 x_Rps12_Rps31 + 0.028 x_Rps12_Ubi4 + 0.105 x_Fun12_Rps13 + 0.005 x_Rps13_Rps31 + 0.028 x_Rps13_Ubi4 + 0.051 x_Fun12_Rps14a + 0.00800000000000001 x_Rps14a_Rps31 + 0.18 x_Rps14a_Ubi4 + 0.051 x_Fun12_Rps14b + 0.00800000000000001 x_Rps14b_Rps31 + 0.18 x_Rps14b_Ubi4 + 0.00800000000000001 x_Fun12_Rps15 + 0.005 x_Rps15_Rps31 + 0.038 x_Rps15_Ubi4 + 0.032 x_Fun12_Rps16a + 0.00800000000000001 x_Rps16a_Rps31 + 0.172 x_Rps16a_Ubi4 + 0.032 x_Fun12_Rps16b + 0.00800000000000001 x_Rps16b_Rps31 + 0.172 x_Rps16b_Ubi4 + 0.032 x_Fun12_Rps17a + 0.00800000000000001 x_Rps17a_Rps31 + 0.151 x_Rps17a_Ubi4 + 0.032 x_Fun12_Rps17b + 0.00800000000000001 x_Rps17b_Rps31 + 0.211 x_Rps17b_Ubi4 + 0.032 x_Fun12_Rps18a + 0.00800000000000001 x_Rps18a_Rps31 + 0.032 x_Fun12_Rps18b + 0.00800000000000001 x_Rps18b_Rps31 + 0.228 x_Rps18b_Ubi4 + 0.032 x_Fun12_Rps19a + 0.00800000000000001 x_Rps19a_Rps31 + 0.185 x_Rps19a_Ubi4 + 0.032 x_Fun12_Rps1a + 0.00800000000000001 x_Rps1a_Rps31 + 0.125 x_Rps1a_Ubi4 + 0.032 x_Fun12_Rps1b + 0.00800000000000001 x_Rps1b_Rps31 + 0.174 x_Rps1b_Ubi4 + 0.104 x_Fun12_Rps2 + 0.005 x_Rps2_Rps31 + 0.038 x_Rps2_Ubi4 + 0.00800000000000001 x_Fun12_Rps20 + 0.005 x_Rps20_Rps31 + 0.029 x_Rps20_Ubi4 + 0.045 x_Fun12_Rps22a + 0.00800000000000001 x_Rps22a_Rps31 + 0.032 x_Fun12_Rps23a + 0.00800000000000001 x_Rps23a_Rps31 + 0.19 x_Rps23a_Ubi4 + 0.032 x_Fun12_Rps23b + 0.00800000000000001 x_Rps23b_Rps31 + 0.19 x_Rps23b_Ubi4 + 0.032 x_Fun12_Rps24a + 0.00800000000000001 x_Rps24a_Rps31 + 0.18 x_Rps24a_Ubi4 + 0.032 x_Fun12_Rps24b + 0.00800000000000001 x_Rps24b_Rps31 + 0.186 x_Rps24b_Ubi4 + 0.16 x_Fun12_Rps25a + 0.00800000000000001 x_Rps25a_Rps31 + 0.206 x_Rps25a_Ubi4 + 0.032 x_Fun12_Rps26a + 0.00800000000000001 x_Rps26a_Rps31 + 0.198 x_Rps26a_Ubi4 + 0.032 x_Fun12_Rps26b + 0.00800000000000001 x_Rps26b_Rps31 + 0.276 x_Rps26b_Ubi4 + 0.032 x_Fun12_Rps27a + 0.00800000000000001 x_Rps27a_Rps31 + 0.156 x_Rps27a_Ubi4 + 0.032 x_Fun12_Rps28a + 0.00800000000000001 x_Rps28a_Rps31 + 0.207 x_Rps28a_Ubi4 + 0.032 x_Fun12_Rps28b + 0.00800000000000001 x_Rps28b_Rps31 + 0.207 x_Rps28b_Ubi4 + 0.032 x_Fun12_Rps29a + 0.00800000000000001 x_Rps29a_Rps31 + 0.208 x_Rps29a_Ubi4 + 0.00800000000000001 x_Fun12_Rps3 + 0.005 x_Rps3_Rps31 + 0.028 x_Rps3_Ubi4 + 0.032 x_Fun12_Rps30a + 0.00800000000000001 x_Rps30a_Rps31 + 0.032 x_Fun12_Rps30b + 0.00800000000000001 x_Rps30b_Rps31 + 0.00800000000000001 x_Fun12_Rps31 + 0.00800000000000001 x_Fun12_Rps5 + 0.277 x_Rps5_Rpn11 + 0.005 x_Rps5_Rps31 + 0.108 x_Rps5_Ubi4 + 0.021 x_Fun12_Rps6a + 0.00800000000000001 x_Rps6a_Rps31 + 0.021 x_Fun12_Rps6b + 0.00800000000000001 x_Rps6b_Rps31 + 0.235 x_Rps6b_Ubi4 + 0.032 x_Fun12_Rps7a + 0.00800000000000001 x_Rps7a_Rps31 + 0.144 x_Rps7a_Ubi4 + 0.032 x_Fun12_Rps8a + 0.00800000000000001 x_Rps8a_Rps31 + 0.209 x_Rps8a_Ubi4 + 0.032 x_Fun12_Rps8b + 0.00800000000000001 x_Rps8b_Rps31 + 0.152 x_Rps8b_Ubi4 + 0.032 x_Fun12_Rps9a + 0.00800000000000001 x_Rps9a_Rps31 + 0.287 x_Rps9a_Ubi4 + 0.032 x_Fun12_Rps9b + 0.00800000000000001 x_Rps9b_Rps31 + 0.287 x_Rps9b_Ubi4 + 0.038 x_Fun12_Sui1 + 0.00800000000000001 x_Sui1_Rps31 + 0.017 x_Fun12_Tif11 + 0.00600000000000001 x_Tif11_Rps31 + 0.026 x_Fun12_Tif35 + 0.078 x_Tif35_Rps31 + 0.229 x_Fun12_Tif5 + 0.174 x_Tif5_Rps31 + 0.29 x_Sup35_Dom34 + 0.044 x_Dom34_Rps31 + 0.12 x_Sup35_Ssa1 + 0.255 x_Ssa1_Cdc48 + 0.039 x_Swe1_Chk1 + 0.033 x_Chk1_Cdc9 + 0.038 x_Cdc9_Rad23 - 2.06546962941424 y_Cyr1 - 0 y_Act1 - 0 y_Cof1 - 0 y_Srv2 - 0 y_Duf1 - 0 y_Ubi4 - 0.113085910813945 y_Cdc48 - 0 y_Ddi1 - 0 y_Dsk2 - 0.134940759104828 y_Rpn1 - 0 y_Rpn10 - 0 y_Rpn11 - 0.227967035078496 y_Rps31 - 0.17079794246036 y_Rpt1 - 0 y_Rpt5 - 0 y_Rpt6 - 0 y_Ufd2 - 1.39100988950645 y_Ede1 - 0.40488716780871 y_Ent1 - 1.02568176866035 y_Chc1 - 0.520447659565873 y_Ent2 - 0.0884911371317518 y_Vps27 - 1.20568332697494 y_Fun12 - 0.297909763660125 y_Asc1 - 0.246136646023169 y_Gcd11 - 0 y_Nip1 - 0.397976742486334 y_Prt1 - 0.414174219075414 y_Rpg1 - 0.183683848587923 y_Rpl10 - 0.260470332775798 y_Rpl17a - 0 y_Rpl17b - 0.136248251868258 y_Rpl19b - 0.270664526640099 y_Rpl25 - 0 y_Rpl28 - 0 y_Rpl2a - 0 y_Rpl2b - 0.37064662837935 y_Rpl3 - 0.0916242099948261 y_Rpl30 - 0.0780246825673995 y_Rpl32 - 0 y_Rpl34a - 0 y_Rpl35a - 0.216836386251174 y_Rpl35b - 0 y_Rpl37a - 0 y_Rpl38 - 0.120065169941506 y_Rpl5 - 0.190773283192771 y_Rpl7a - 0 y_Rps0a - 0 y_Rps0b - 0 y_Rps10a - 0 y_Rps11a - 0.276867820406057 y_Rps11b - 0 y_Rps12 - 0 y_Rps13 - 0 y_Rps14a - 0.21043919262963 y_Rps14b - 0 y_Rps15 - 0 y_Rps16a - 0.0917888101035867 y_Rps16b - 0 y_Rps17a - 0.20078081333483 y_Rps17b - 0.33165617846531 y_Rps18a - 0 y_Rps18b - 0 y_Rps19a - 0.181523970358348 y_Rps1a - 0.511183741969515 y_Rps1b - 0 y_Rps2 - 0 y_Rps20 - 0 y_Rps22a - 0 y_Rps23a - 0 y_Rps23b - 0.181971262328775 y_Rps24a - 0 y_Rps24b - 0 y_Rps25a - 0 y_Rps26a - 0.128719434628373 y_Rps26b - 0 y_Rps27a - 0 y_Rps28a - 0 y_Rps28b - 0.092884269459382 y_Rps29a - 0.118468258343378 y_Rps3 - 0 y_Rps30a - 0 y_Rps30b - 0 y_Rps5 - 0 y_Rps6a - 0.129832488722207 y_Rps6b - 0.149603196222064 y_Rps7a - 0.219566469769258 y_Rps8a - 0 y_Rps8b - 0 y_Rps9a - 0.22725818108385 y_Rps9b - 0.154199226094496 y_Sui1 - 0 y_Tif11 - 0.18472471674679 y_Tif35 - 0.216773820733091 y_Tif5 - 1.53949924297249 y_Sup35 - 0 y_Dom34 - 0.10510464059825 y_Ssa1 - 1.21946220627969 y_Swe1 - 0 y_Chk1 - 0.384891476437097 y_Cdc9 - 0 y_Prm15 - 1.54767341365104 y_Rad23 - 1.46860721495992 y_Ent5

subject to
x_Cyr1_Act1 + x_Cyr1_Srv2 >= 1 \edges from terminal
x_Duf1_Ubi4 >= 1 \edges from terminal
x_Ede1_Ubi4 + x_Ede1_Rps31 + x_Ede1_Ent1 + x_Ede1_Ent2 + x_Ede1_Vps27 >= 1 \edges from terminal
x_Fun12_Rps31 + x_Fun12_Asc1 + x_Fun12_Gcd11 + x_Fun12_Nip1 + x_Fun12_Prt1 + x_Fun12_Rpg1 + x_Fun12_Rpl10 + x_Fun12_Rpl17a + x_Fun12_Rpl17b + x_Fun12_Rpl19b + x_Fun12_Rpl25 + x_Fun12_Rpl28 + x_Fun12_Rpl2a + x_Fun12_Rpl2b + x_Fun12_Rpl3 + x_Fun12_Rpl30 + x_Fun12_Rpl32 + x_Fun12_Rpl34a + x_Fun12_Rpl35a + x_Fun12_Rpl35b + x_Fun12_Rpl37a + x_Fun12_Rpl38 + x_Fun12_Rpl5 + x_Fun12_Rpl7a + x_Fun12_Rps0a + x_Fun12_Rps0b + x_Fun12_Rps10a + x_Fun12_Rps11a + x_Fun12_Rps11b + x_Fun12_Rps12 + x_Fun12_Rps13 + x_Fun12_Rps14a + x_Fun12_Rps14b + x_Fun12_Rps15 + x_Fun12_Rps16a + x_Fun12_Rps16b + x_Fun12_Rps17a + x_Fun12_Rps17b + x_Fun12_Rps18a + x_Fun12_Rps18b + x_Fun12_Rps19a + x_Fun12_Rps1a + x_Fun12_Rps1b + x_Fun12_Rps2 + x_Fun12_Rps20 + x_Fun12_Rps22a + x_Fun12_Rps23a + x_Fun12_Rps23b + x_Fun12_Rps24a + x_Fun12_Rps24b + x_Fun12_Rps25a + x_Fun12_Rps26a + x_Fun12_Rps26b + x_Fun12_Rps27a + x_Fun12_Rps28a + x_Fun12_Rps28b + x_Fun12_Rps29a + x_Fun12_Rps3 + x_Fun12_Rps30a + x_Fun12_Rps30b + x_Fun12_Rps5 + x_Fun12_Rps6a + x_Fun12_Rps6b + x_Fun12_Rps7a + x_Fun12_Rps8a + x_Fun12_Rps8b + x_Fun12_Rps9a + x_Fun12_Rps9b + x_Fun12_Sui1 + x_Fun12_Tif11 + x_Fun12_Tif35 + x_Fun12_Tif5 >= 1 \edges from terminal
x_Sup35_Dom34 + x_Sup35_Ssa1 >= 1 \edges from terminal
x_Swe1_Chk1 >= 1 \edges from terminal
x_Cof1_Prm15 >= 1 \edges to terminal
x_Ubi4_Rad23 + x_Cdc48_Rad23 + x_Ddi1_Rad23 + x_Dsk2_Rad23 + x_Rpn1_Rad23 + x_Rpn10_Rad23 + x_Rpn11_Rad23 + x_Rps31_Rad23 + x_Rpt1_Rad23 + x_Rpt5_Rad23 + x_Rpt6_Rad23 + x_Ufd2_Rad23 + x_Cdc9_Rad23 >= 1 \edges to terminal
x_Chc1_Ent5 >= 1 \edges to terminal
y_Cyr1 >= 1  \terminal
y_Duf1 >= 1  \terminal
y_Ede1 >= 1  \terminal
y_Fun12 >= 1  \terminal
y_Sup35 >= 1  \terminal
y_Swe1 >= 1  \terminal
y_Prm15 >= 1  \terminal
y_Rad23 >= 1  \terminal
y_Ent5 >= 1  \terminal
y_Act1 - x_Act1_Cof1 <= 0 \Steiner node out
y_Cof1 - x_Cof1_Prm15 <= 0 \Steiner node out
y_Srv2 - x_Srv2_Cof1 <= 0 \Steiner node out
y_Ubi4 - x_Ubi4_Cdc48 - x_Ubi4_Ddi1 - x_Ubi4_Dsk2 - x_Ubi4_Rpn1 - x_Ubi4_Rpn10 - x_Ubi4_Rpn11 - x_Ubi4_Rps31 - x_Ubi4_Rpt1 - x_Ubi4_Rpt5 - x_Ubi4_Rpt6 - x_Ubi4_Ufd2 - x_Ubi4_Rad23 <= 0 \Steiner node out
y_Cdc48 - x_Cdc48_Rad23 <= 0 \Steiner node out
y_Ddi1 - x_Ddi1_Rad23 <= 0 \Steiner node out
y_Dsk2 - x_Dsk2_Rad23 <= 0 \Steiner node out
y_Rpn1 - x_Rpn1_Rad23 <= 0 \Steiner node out
y_Rpn10 - x_Rpn10_Rad23 <= 0 \Steiner node out
y_Rpn11 - x_Rpn11_Rad23 <= 0 \Steiner node out
y_Rps31 - x_Rps31_Ubi4 - x_Rps31_Cdc48 - x_Rps31_Dsk2 - x_Rps31_Rpn1 - x_Rps31_Rpn10 - x_Rps31_Rpn11 - x_Rps31_Rpt1 - x_Rps31_Rpt5 - x_Rps31_Rpt6 - x_Rps31_Rad23 <= 0 \Steiner node out
y_Rpt1 - x_Rpt1_Rad23 <= 0 \Steiner node out
y_Rpt5 - x_Rpt5_Rad23 <= 0 \Steiner node out
y_Rpt6 - x_Rpt6_Rad23 <= 0 \Steiner node out
y_Ufd2 - x_Ufd2_Rad23 <= 0 \Steiner node out
y_Ent1 - x_Ent1_Ubi4 - x_Ent1_Chc1 <= 0 \Steiner node out
y_Chc1 - x_Chc1_Ent5 <= 0 \Steiner node out
y_Ent2 - x_Ent2_Ubi4 - x_Ent2_Chc1 <= 0 \Steiner node out
y_Vps27 - x_Vps27_Ubi4 - x_Vps27_Rps31 - x_Vps27_Chc1 <= 0 \Steiner node out
y_Asc1 - x_Asc1_Ubi4 - x_Asc1_Rps31 <= 0 \Steiner node out
y_Gcd11 - x_Gcd11_Rps31 <= 0 \Steiner node out
y_Nip1 - x_Nip1_Rpn11 - x_Nip1_Rps31 <= 0 \Steiner node out
y_Prt1 - x_Prt1_Rps31 <= 0 \Steiner node out
y_Rpg1 - x_Rpg1_Rps31 <= 0 \Steiner node out
y_Rpl10 - x_Rpl10_Ubi4 - x_Rpl10_Rps31 <= 0 \Steiner node out
y_Rpl17a - x_Rpl17a_Ubi4 - x_Rpl17a_Rps31 <= 0 \Steiner node out
y_Rpl17b - x_Rpl17b_Ubi4 - x_Rpl17b_Rps31 <= 0 \Steiner node out
y_Rpl19b - x_Rpl19b_Ubi4 - x_Rpl19b_Rps31 <= 0 \Steiner node out
y_Rpl25 - x_Rpl25_Ubi4 - x_Rpl25_Rps31 <= 0 \Steiner node out
y_Rpl28 - x_Rpl28_Ubi4 - x_Rpl28_Rps31 <= 0 \Steiner node out
y_Rpl2a - x_Rpl2a_Ubi4 - x_Rpl2a_Rps31 <= 0 \Steiner node out
y_Rpl2b - x_Rpl2b_Ubi4 - x_Rpl2b_Rps31 <= 0 \Steiner node out
y_Rpl3 - x_Rpl3_Ubi4 - x_Rpl3_Rps31 <= 0 \Steiner node out
y_Rpl30 - x_Rpl30_Ubi4 - x_Rpl30_Rps31 <= 0 \Steiner node out
y_Rpl32 - x_Rpl32_Ubi4 - x_Rpl32_Rps31 <= 0 \Steiner node out
y_Rpl34a - x_Rpl34a_Ubi4 - x_Rpl34a_Rps31 <= 0 \Steiner node out
y_Rpl35a - x_Rpl35a_Ubi4 - x_Rpl35a_Rps31 <= 0 \Steiner node out
y_Rpl35b - x_Rpl35b_Ubi4 - x_Rpl35b_Rps31 <= 0 \Steiner node out
y_Rpl37a - x_Rpl37a_Rps31 <= 0 \Steiner node out
y_Rpl38 - x_Rpl38_Ubi4 - x_Rpl38_Rps31 <= 0 \Steiner node out
y_Rpl5 - x_Rpl5_Ubi4 - x_Rpl5_Rps31 <= 0 \Steiner node out
y_Rpl7a - x_Rpl7a_Rps31 <= 0 \Steiner node out
y_Rps0a - x_Rps0a_Ubi4 - x_Rps0a_Rps31 <= 0 \Steiner node out
y_Rps0b - x_Rps0b_Ubi4 - x_Rps0b_Rps31 <= 0 \Steiner node out
y_Rps10a - x_Rps10a_Ubi4 - x_Rps10a_Rps31 <= 0 \Steiner node out
y_Rps11a - x_Rps11a_Ubi4 - x_Rps11a_Rps31 <= 0 \Steiner node out
y_Rps11b - x_Rps11b_Ubi4 - x_Rps11b_Rps31 <= 0 \Steiner node out
y_Rps12 - x_Rps12_Ubi4 - x_Rps12_Rps31 <= 0 \Steiner node out
y_Rps13 - x_Rps13_Ubi4 - x_Rps13_Rps31 <= 0 \Steiner node out
y_Rps14a - x_Rps14a_Ubi4 - x_Rps14a_Rps31 <= 0 \Steiner node out
y_Rps14b - x_Rps14b_Ubi4 - x_Rps14b_Rps31 <= 0 \Steiner node out
y_Rps15 - x_Rps15_Ubi4 - x_Rps15_Rps31 <= 0 \Steiner node out
y_Rps16a - x_Rps16a_Ubi4 - x_Rps16a_Rps31 <= 0 \Steiner node out
y_Rps16b - x_Rps16b_Ubi4 - x_Rps16b_Rps31 <= 0 \Steiner node out
y_Rps17a - x_Rps17a_Ubi4 - x_Rps17a_Rps31 <= 0 \Steiner node out
y_Rps17b - x_Rps17b_Ubi4 - x_Rps17b_Rps31 <= 0 \Steiner node out
y_Rps18a - x_Rps18a_Rps31 <= 0 \Steiner node out
y_Rps18b - x_Rps18b_Ubi4 - x_Rps18b_Rps31 <= 0 \Steiner node out
y_Rps19a - x_Rps19a_Ubi4 - x_Rps19a_Rps31 <= 0 \Steiner node out
y_Rps1a - x_Rps1a_Ubi4 - x_Rps1a_Rps31 <= 0 \Steiner node out
y_Rps1b - x_Rps1b_Ubi4 - x_Rps1b_Rps31 <= 0 \Steiner node out
y_Rps2 - x_Rps2_Ubi4 - x_Rps2_Rps31 <= 0 \Steiner node out
y_Rps20 - x_Rps20_Ubi4 - x_Rps20_Rps31 <= 0 \Steiner node out
y_Rps22a - x_Rps22a_Rps31 <= 0 \Steiner node out
y_Rps23a - x_Rps23a_Ubi4 - x_Rps23a_Rps31 <= 0 \Steiner node out
y_Rps23b - x_Rps23b_Ubi4 - x_Rps23b_Rps31 <= 0 \Steiner node out
y_Rps24a - x_Rps24a_Ubi4 - x_Rps24a_Rps31 <= 0 \Steiner node out
y_Rps24b - x_Rps24b_Ubi4 - x_Rps24b_Rps31 <= 0 \Steiner node out
y_Rps25a - x_Rps25a_Ubi4 - x_Rps25a_Rps31 <= 0 \Steiner node out
y_Rps26a - x_Rps26a_Ubi4 - x_Rps26a_Rps31 <= 0 \Steiner node out
y_Rps26b - x_Rps26b_Ubi4 - x_Rps26b_Rps31 <= 0 \Steiner node out
y_Rps27a - x_Rps27a_Ubi4 - x_Rps27a_Rps31 <= 0 \Steiner node out
y_Rps28a - x_Rps28a_Ubi4 - x_Rps28a_Rps31 <= 0 \Steiner node out
y_Rps28b - x_Rps28b_Ubi4 - x_Rps28b_Rps31 <= 0 \Steiner node out
y_Rps29a - x_Rps29a_Ubi4 - x_Rps29a_Rps31 <= 0 \Steiner node out
y_Rps3 - x_Rps3_Ubi4 - x_Rps3_Rps31 <= 0 \Steiner node out
y_Rps30a - x_Rps30a_Rps31 <= 0 \Steiner node out
y_Rps30b - x_Rps30b_Rps31 <= 0 \Steiner node out
y_Rps5 - x_Rps5_Ubi4 - x_Rps5_Rpn11 - x_Rps5_Rps31 <= 0 \Steiner node out
y_Rps6a - x_Rps6a_Rps31 <= 0 \Steiner node out
y_Rps6b - x_Rps6b_Ubi4 - x_Rps6b_Rps31 <= 0 \Steiner node out
y_Rps7a - x_Rps7a_Ubi4 - x_Rps7a_Rps31 <= 0 \Steiner node out
y_Rps8a - x_Rps8a_Ubi4 - x_Rps8a_Rps31 <= 0 \Steiner node out
y_Rps8b - x_Rps8b_Ubi4 - x_Rps8b_Rps31 <= 0 \Steiner node out
y_Rps9a - x_Rps9a_Ubi4 - x_Rps9a_Rps31 <= 0 \Steiner node out
y_Rps9b - x_Rps9b_Ubi4 - x_Rps9b_Rps31 <= 0 \Steiner node out
y_Sui1 - x_Sui1_Rps31 <= 0 \Steiner node out
y_Tif11 - x_Tif11_Rps31 <= 0 \Steiner node out
y_Tif35 - x_Tif35_Rps31 <= 0 \Steiner node out
y_Tif5 - x_Tif5_Rps31 <= 0 \Steiner node out
y_Dom34 - x_Dom34_Rps31 <= 0 \Steiner node out
y_Ssa1 - x_Ssa1_Cdc48 <= 0 \Steiner node out
y_Chk1 - x_Chk1_Cdc9 <= 0 \Steiner node out
y_Cdc9 - x_Cdc9_Rad23 <= 0 \Steiner node out
y_Act1 - x_Cyr1_Act1 <= 0 \Steiner node in
y_Cof1 - x_Act1_Cof1 - x_Srv2_Cof1 <= 0 \Steiner node in
y_Srv2 - x_Cyr1_Srv2 <= 0 \Steiner node in
y_Ubi4 - x_Duf1_Ubi4 - x_Rps31_Ubi4 - x_Ede1_Ubi4 - x_Ent1_Ubi4 - x_Ent2_Ubi4 - x_Vps27_Ubi4 - x_Asc1_Ubi4 - x_Rpl10_Ubi4 - x_Rpl17a_Ubi4 - x_Rpl17b_Ubi4 - x_Rpl19b_Ubi4 - x_Rpl25_Ubi4 - x_Rpl28_Ubi4 - x_Rpl2a_Ubi4 - x_Rpl2b_Ubi4 - x_Rpl3_Ubi4 - x_Rpl30_Ubi4 - x_Rpl32_Ubi4 - x_Rpl34a_Ubi4 - x_Rpl35a_Ubi4 - x_Rpl35b_Ubi4 - x_Rpl38_Ubi4 - x_Rpl5_Ubi4 - x_Rps0a_Ubi4 - x_Rps0b_Ubi4 - x_Rps10a_Ubi4 - x_Rps11a_Ubi4 - x_Rps11b_Ubi4 - x_Rps12_Ubi4 - x_Rps13_Ubi4 - x_Rps14a_Ubi4 - x_Rps14b_Ubi4 - x_Rps15_Ubi4 - x_Rps16a_Ubi4 - x_Rps16b_Ubi4 - x_Rps17a_Ubi4 - x_Rps17b_Ubi4 - x_Rps18b_Ubi4 - x_Rps19a_Ubi4 - x_Rps1a_Ubi4 - x_Rps1b_Ubi4 - x_Rps2_Ubi4 - x_Rps20_Ubi4 - x_Rps23a_Ubi4 - x_Rps23b_Ubi4 - x_Rps24a_Ubi4 - x_Rps24b_Ubi4 - x_Rps25a_Ubi4 - x_Rps26a_Ubi4 - x_Rps26b_Ubi4 - x_Rps27a_Ubi4 - x_Rps28a_Ubi4 - x_Rps28b_Ubi4 - x_Rps29a_Ubi4 - x_Rps3_Ubi4 - x_Rps5_Ubi4 - x_Rps6b_Ubi4 - x_Rps7a_Ubi4 - x_Rps8a_Ubi4 - x_Rps8b_Ubi4 - x_Rps9a_Ubi4 - x_Rps9b_Ubi4 <= 0 \Steiner node in
y_Cdc48 - x_Ubi4_Cdc48 - x_Rps31_Cdc48 - x_Ssa1_Cdc48 <= 0 \Steiner node in
y_Ddi1 - x_Ubi4_Ddi1 <= 0 \Steiner node in
y_Dsk2 - x_Ubi4_Dsk2 - x_Rps31_Dsk2 <= 0 \Steiner node in
y_Rpn1 - x_Ubi4_Rpn1 - x_Rps31_Rpn1 <= 0 \Steiner node in
y_Rpn10 - x_Ubi4_Rpn10 - x_Rps31_Rpn10 <= 0 \Steiner node in
y_Rpn11 - x_Ubi4_Rpn11 - x_Rps31_Rpn11 - x_Nip1_Rpn11 - x_Rps5_Rpn11 <= 0 \Steiner node in
y_Rps31 - x_Ubi4_Rps31 - x_Ede1_Rps31 - x_Vps27_Rps31 - x_Fun12_Rps31 - x_Asc1_Rps31 - x_Gcd11_Rps31 - x_Nip1_Rps31 - x_Prt1_Rps31 - x_Rpg1_Rps31 - x_Rpl10_Rps31 - x_Rpl17a_Rps31 - x_Rpl17b_Rps31 - x_Rpl19b_Rps31 - x_Rpl25_Rps31 - x_Rpl28_Rps31 - x_Rpl2a_Rps31 - x_Rpl2b_Rps31 - x_Rpl3_Rps31 - x_Rpl30_Rps31 - x_Rpl32_Rps31 - x_Rpl34a_Rps31 - x_Rpl35a_Rps31 - x_Rpl35b_Rps31 - x_Rpl37a_Rps31 - x_Rpl38_Rps31 - x_Rpl5_Rps31 - x_Rpl7a_Rps31 - x_Rps0a_Rps31 - x_Rps0b_Rps31 - x_Rps10a_Rps31 - x_Rps11a_Rps31 - x_Rps11b_Rps31 - x_Rps12_Rps31 - x_Rps13_Rps31 - x_Rps14a_Rps31 - x_Rps14b_Rps31 - x_Rps15_Rps31 - x_Rps16a_Rps31 - x_Rps16b_Rps31 - x_Rps17a_Rps31 - x_Rps17b_Rps31 - x_Rps18a_Rps31 - x_Rps18b_Rps31 - x_Rps19a_Rps31 - x_Rps1a_Rps31 - x_Rps1b_Rps31 - x_Rps2_Rps31 - x_Rps20_Rps31 - x_Rps22a_Rps31 - x_Rps23a_Rps31 - x_Rps23b_Rps31 - x_Rps24a_Rps31 - x_Rps24b_Rps31 - x_Rps25a_Rps31 - x_Rps26a_Rps31 - x_Rps26b_Rps31 - x_Rps27a_Rps31 - x_Rps28a_Rps31 - x_Rps28b_Rps31 - x_Rps29a_Rps31 - x_Rps3_Rps31 - x_Rps30a_Rps31 - x_Rps30b_Rps31 - x_Rps5_Rps31 - x_Rps6a_Rps31 - x_Rps6b_Rps31 - x_Rps7a_Rps31 - x_Rps8a_Rps31 - x_Rps8b_Rps31 - x_Rps9a_Rps31 - x_Rps9b_Rps31 - x_Sui1_Rps31 - x_Tif11_Rps31 - x_Tif35_Rps31 - x_Tif5_Rps31 - x_Dom34_Rps31 <= 0 \Steiner node in
y_Rpt1 - x_Ubi4_Rpt1 - x_Rps31_Rpt1 <= 0 \Steiner node in
y_Rpt5 - x_Ubi4_Rpt5 - x_Rps31_Rpt5 <= 0 \Steiner node in
y_Rpt6 - x_Ubi4_Rpt6 - x_Rps31_Rpt6 <= 0 \Steiner node in
y_Ufd2 - x_Ubi4_Ufd2 <= 0 \Steiner node in
y_Ent1 - x_Ede1_Ent1 <= 0 \Steiner node in
y_Chc1 - x_Ent1_Chc1 - x_Ent2_Chc1 - x_Vps27_Chc1 <= 0 \Steiner node in
y_Ent2 - x_Ede1_Ent2 <= 0 \Steiner node in
y_Vps27 - x_Ede1_Vps27 <= 0 \Steiner node in
y_Asc1 - x_Fun12_Asc1 <= 0 \Steiner node in
y_Gcd11 - x_Fun12_Gcd11 <= 0 \Steiner node in
y_Nip1 - x_Fun12_Nip1 <= 0 \Steiner node in
y_Prt1 - x_Fun12_Prt1 <= 0 \Steiner node in
y_Rpg1 - x_Fun12_Rpg1 <= 0 \Steiner node in
y_Rpl10 - x_Fun12_Rpl10 <= 0 \Steiner node in
y_Rpl17a - x_Fun12_Rpl17a <= 0 \Steiner node in
y_Rpl17b - x_Fun12_Rpl17b <= 0 \Steiner node in
y_Rpl19b - x_Fun12_Rpl19b <= 0 \Steiner node in
y_Rpl25 - x_Fun12_Rpl25 <= 0 \Steiner node in
y_Rpl28 - x_Fun12_Rpl28 <= 0 \Steiner node in
y_Rpl2a - x_Fun12_Rpl2a <= 0 \Steiner node in
y_Rpl2b - x_Fun12_Rpl2b <= 0 \Steiner node in
y_Rpl3 - x_Fun12_Rpl3 <= 0 \Steiner node in
y_Rpl30 - x_Fun12_Rpl30 <= 0 \Steiner node in
y_Rpl32 - x_Fun12_Rpl32 <= 0 \Steiner node in
y_Rpl34a - x_Fun12_Rpl34a <= 0 \Steiner node in
y_Rpl35a - x_Fun12_Rpl35a <= 0 \Steiner node in
y_Rpl35b - x_Fun12_Rpl35b <= 0 \Steiner node in
y_Rpl37a - x_Fun12_Rpl37a <= 0 \Steiner node in
y_Rpl38 - x_Fun12_Rpl38 <= 0 \Steiner node in
y_Rpl5 - x_Fun12_Rpl5 <= 0 \Steiner node in
y_Rpl7a - x_Fun12_Rpl7a <= 0 \Steiner node in
y_Rps0a - x_Fun12_Rps0a <= 0 \Steiner node in
y_Rps0b - x_Fun12_Rps0b <= 0 \Steiner node in
y_Rps10a - x_Fun12_Rps10a <= 0 \Steiner node in
y_Rps11a - x_Fun12_Rps11a <= 0 \Steiner node in
y_Rps11b - x_Fun12_Rps11b <= 0 \Steiner node in
y_Rps12 - x_Fun12_Rps12 <= 0 \Steiner node in
y_Rps13 - x_Fun12_Rps13 <= 0 \Steiner node in
y_Rps14a - x_Fun12_Rps14a <= 0 \Steiner node in
y_Rps14b - x_Fun12_Rps14b <= 0 \Steiner node in
y_Rps15 - x_Fun12_Rps15 <= 0 \Steiner node in
y_Rps16a - x_Fun12_Rps16a <= 0 \Steiner node in
y_Rps16b - x_Fun12_Rps16b <= 0 \Steiner node in
y_Rps17a - x_Fun12_Rps17a <= 0 \Steiner node in
y_Rps17b - x_Fun12_Rps17b <= 0 \Steiner node in
y_Rps18a - x_Fun12_Rps18a <= 0 \Steiner node in
y_Rps18b - x_Fun12_Rps18b <= 0 \Steiner node in
y_Rps19a - x_Fun12_Rps19a <= 0 \Steiner node in
y_Rps1a - x_Fun12_Rps1a <= 0 \Steiner node in
y_Rps1b - x_Fun12_Rps1b <= 0 \Steiner node in
y_Rps2 - x_Fun12_Rps2 <= 0 \Steiner node in
y_Rps20 - x_Fun12_Rps20 <= 0 \Steiner node in
y_Rps22a - x_Fun12_Rps22a <= 0 \Steiner node in
y_Rps23a - x_Fun12_Rps23a <= 0 \Steiner node in
y_Rps23b - x_Fun12_Rps23b <= 0 \Steiner node in
y_Rps24a - x_Fun12_Rps24a <= 0 \Steiner node in
y_Rps24b - x_Fun12_Rps24b <= 0 \Steiner node in
y_Rps25a - x_Fun12_Rps25a <= 0 \Steiner node in
y_Rps26a - x_Fun12_Rps26a <= 0 \Steiner node in
y_Rps26b - x_Fun12_Rps26b <= 0 \Steiner node in
y_Rps27a - x_Fun12_Rps27a <= 0 \Steiner node in
y_Rps28a - x_Fun12_Rps28a <= 0 \Steiner node in
y_Rps28b - x_Fun12_Rps28b <= 0 \Steiner node in
y_Rps29a - x_Fun12_Rps29a <= 0 \Steiner node in
y_Rps3 - x_Fun12_Rps3 <= 0 \Steiner node in
y_Rps30a - x_Fun12_Rps30a <= 0 \Steiner node in
y_Rps30b - x_Fun12_Rps30b <= 0 \Steiner node in
y_Rps5 - x_Fun12_Rps5 <= 0 \Steiner node in
y_Rps6a - x_Fun12_Rps6a <= 0 \Steiner node in
y_Rps6b - x_Fun12_Rps6b <= 0 \Steiner node in
y_Rps7a - x_Fun12_Rps7a <= 0 \Steiner node in
y_Rps8a - x_Fun12_Rps8a <= 0 \Steiner node in
y_Rps8b - x_Fun12_Rps8b <= 0 \Steiner node in
y_Rps9a - x_Fun12_Rps9a <= 0 \Steiner node in
y_Rps9b - x_Fun12_Rps9b <= 0 \Steiner node in
y_Sui1 - x_Fun12_Sui1 <= 0 \Steiner node in
y_Tif11 - x_Fun12_Tif11 <= 0 \Steiner node in
y_Tif35 - x_Fun12_Tif35 <= 0 \Steiner node in
y_Tif5 - x_Fun12_Tif5 <= 0 \Steiner node in
y_Dom34 - x_Sup35_Dom34 <= 0 \Steiner node in
y_Ssa1 - x_Sup35_Ssa1 <= 0 \Steiner node in
y_Chk1 - x_Swe1_Chk1 <= 0 \Steiner node in
y_Cdc9 - x_Chk1_Cdc9 <= 0 \Steiner node in
y_Cyr1 - x_Cyr1_Act1 >= 0 \select node if connected (outgoing)
y_Cyr1 - x_Cyr1_Srv2 >= 0 \select node if connected (outgoing)
y_Act1 - x_Act1_Cof1 >= 0 \select node if connected (outgoing)
y_Cof1 - x_Cof1_Prm15 >= 0 \select node if connected (outgoing)
y_Srv2 - x_Srv2_Cof1 >= 0 \select node if connected (outgoing)
y_Duf1 - x_Duf1_Ubi4 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Cdc48 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Ddi1 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Dsk2 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpn1 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpn10 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpn11 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps31 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpt1 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpt5 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpt6 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Ufd2 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rad23 >= 0 \select node if connected (outgoing)
y_Cdc48 - x_Cdc48_Rad23 >= 0 \select node if connected (outgoing)
y_Ddi1 - x_Ddi1_Rad23 >= 0 \select node if connected (outgoing)
y_Dsk2 - x_Dsk2_Rad23 >= 0 \select node if connected (outgoing)
y_Rpn1 - x_Rpn1_Rad23 >= 0 \select node if connected (outgoing)
y_Rpn10 - x_Rpn10_Rad23 >= 0 \select node if connected (outgoing)
y_Rpn11 - x_Rpn11_Rad23 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Cdc48 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Dsk2 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpn1 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpn10 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpn11 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpt1 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpt5 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpt6 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rad23 >= 0 \select node if connected (outgoing)
y_Rpt1 - x_Rpt1_Rad23 >= 0 \select node if connected (outgoing)
y_Rpt5 - x_Rpt5_Rad23 >= 0 \select node if connected (outgoing)
y_Rpt6 - x_Rpt6_Rad23 >= 0 \select node if connected (outgoing)
y_Ufd2 - x_Ufd2_Rad23 >= 0 \select node if connected (outgoing)
y_Ede1 - x_Ede1_Ubi4 >= 0 \select node if connected (outgoing)
y_Ede1 - x_Ede1_Rps31 >= 0 \select node if connected (outgoing)
y_Ede1 - x_Ede1_Ent1 >= 0 \select node if connected (outgoing)
y_Ede1 - x_Ede1_Ent2 >= 0 \select node if connected (outgoing)
y_Ede1 - x_Ede1_Vps27 >= 0 \select node if connected (outgoing)
y_Ent1 - x_Ent1_Ubi4 >= 0 \select node if connected (outgoing)
y_Ent1 - x_Ent1_Chc1 >= 0 \select node if connected (outgoing)
y_Chc1 - x_Chc1_Ent5 >= 0 \select node if connected (outgoing)
y_Ent2 - x_Ent2_Ubi4 >= 0 \select node if connected (outgoing)
y_Ent2 - x_Ent2_Chc1 >= 0 \select node if connected (outgoing)
y_Vps27 - x_Vps27_Ubi4 >= 0 \select node if connected (outgoing)
y_Vps27 - x_Vps27_Rps31 >= 0 \select node if connected (outgoing)
y_Vps27 - x_Vps27_Chc1 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps31 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Asc1 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Gcd11 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Nip1 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Prt1 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpg1 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl10 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl17a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl17b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl19b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl25 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl28 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl2a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl2b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl3 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl30 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl32 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl34a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl35a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl35b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl37a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl38 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl5 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rpl7a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps0a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps0b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps10a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps11a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps11b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps12 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps13 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps14a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps14b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps15 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps16a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps16b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps17a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps17b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps18a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps18b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps19a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps1a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps1b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps2 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps20 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps22a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps23a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps23b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps24a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps24b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps25a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps26a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps26b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps27a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps28a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps28b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps29a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps3 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps30a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps30b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps5 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps6a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps6b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps7a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps8a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps8b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps9a >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Rps9b >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Sui1 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Tif11 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Tif35 >= 0 \select node if connected (outgoing)
y_Fun12 - x_Fun12_Tif5 >= 0 \select node if connected (outgoing)
y_Asc1 - x_Asc1_Ubi4 >= 0 \select node if connected (outgoing)
y_Asc1 - x_Asc1_Rps31 >= 0 \select node if connected (outgoing)
y_Gcd11 - x_Gcd11_Rps31 >= 0 \select node if connected (outgoing)
y_Nip1 - x_Nip1_Rpn11 >= 0 \select node if connected (outgoing)
y_Nip1 - x_Nip1_Rps31 >= 0 \select node if connected (outgoing)
y_Prt1 - x_Prt1_Rps31 >= 0 \select node if connected (outgoing)
y_Rpg1 - x_Rpg1_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl10 - x_Rpl10_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl10 - x_Rpl10_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl17a - x_Rpl17a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl17a - x_Rpl17a_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl17b - x_Rpl17b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl17b - x_Rpl17b_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl19b - x_Rpl19b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl19b - x_Rpl19b_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl25 - x_Rpl25_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl25 - x_Rpl25_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl28 - x_Rpl28_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl28 - x_Rpl28_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl2a - x_Rpl2a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl2a - x_Rpl2a_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl2b - x_Rpl2b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl2b - x_Rpl2b_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl3 - x_Rpl3_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl3 - x_Rpl3_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl30 - x_Rpl30_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl30 - x_Rpl30_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl32 - x_Rpl32_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl32 - x_Rpl32_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl34a - x_Rpl34a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl34a - x_Rpl34a_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl35a - x_Rpl35a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl35a - x_Rpl35a_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl35b - x_Rpl35b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl35b - x_Rpl35b_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl37a - x_Rpl37a_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl38 - x_Rpl38_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl38 - x_Rpl38_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl5 - x_Rpl5_Ubi4 >= 0 \select node if connected (outgoing)
y_Rpl5 - x_Rpl5_Rps31 >= 0 \select node if connected (outgoing)
y_Rpl7a - x_Rpl7a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps0a - x_Rps0a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps0a - x_Rps0a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps0b - x_Rps0b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps0b - x_Rps0b_Rps31 >= 0 \select node if connected (outgoing)
y_Rps10a - x_Rps10a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps10a - x_Rps10a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps11a - x_Rps11a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps11a - x_Rps11a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps11b - x_Rps11b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps11b - x_Rps11b_Rps31 >= 0 \select node if connected (outgoing)
y_Rps12 - x_Rps12_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps12 - x_Rps12_Rps31 >= 0 \select node if connected (outgoing)
y_Rps13 - x_Rps13_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps13 - x_Rps13_Rps31 >= 0 \select node if connected (outgoing)
y_Rps14a - x_Rps14a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps14a - x_Rps14a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps14b - x_Rps14b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps14b - x_Rps14b_Rps31 >= 0 \select node if connected (outgoing)
y_Rps15 - x_Rps15_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps15 - x_Rps15_Rps31 >= 0 \select node if connected (outgoing)
y_Rps16a - x_Rps16a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps16a - x_Rps16a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps16b - x_Rps16b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps16b - x_Rps16b_Rps31 >= 0 \select node if connected (outgoing)
y_Rps17a - x_Rps17a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps17a - x_Rps17a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps17b - x_Rps17b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps17b - x_Rps17b_Rps31 >= 0 \select node if connected (outgoing)
y_Rps18a - x_Rps18a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps18b - x_Rps18b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps18b - x_Rps18b_Rps31 >= 0 \select node if connected (outgoing)
y_Rps19a - x_Rps19a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps19a - x_Rps19a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps1a - x_Rps1a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps1a - x_Rps1a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps1b - x_Rps1b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps1b - x_Rps1b_Rps31 >= 0 \select node if connected (outgoing)
y_Rps2 - x_Rps2_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps2 - x_Rps2_Rps31 >= 0 \select node if connected (outgoing)
y_Rps20 - x_Rps20_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps20 - x_Rps20_Rps31 >= 0 \select node if connected (outgoing)
y_Rps22a - x_Rps22a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps23a - x_Rps23a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps23a - x_Rps23a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps23b - x_Rps23b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps23b - x_Rps23b_Rps31 >= 0 \select node if connected (outgoing)
y_Rps24a - x_Rps24a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps24a - x_Rps24a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps24b - x_Rps24b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps24b - x_Rps24b_Rps31 >= 0 \select node if connected (outgoing)
y_Rps25a - x_Rps25a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps25a - x_Rps25a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps26a - x_Rps26a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps26a - x_Rps26a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps26b - x_Rps26b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps26b - x_Rps26b_Rps31 >= 0 \select node if connected (outgoing)
y_Rps27a - x_Rps27a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps27a - x_Rps27a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps28a - x_Rps28a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps28a - x_Rps28a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps28b - x_Rps28b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps28b - x_Rps28b_Rps31 >= 0 \select node if connected (outgoing)
y_Rps29a - x_Rps29a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps29a - x_Rps29a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps3 - x_Rps3_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps3 - x_Rps3_Rps31 >= 0 \select node if connected (outgoing)
y_Rps30a - x_Rps30a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps30b - x_Rps30b_Rps31 >= 0 \select node if connected (outgoing)
y_Rps5 - x_Rps5_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps5 - x_Rps5_Rpn11 >= 0 \select node if connected (outgoing)
y_Rps5 - x_Rps5_Rps31 >= 0 \select node if connected (outgoing)
y_Rps6a - x_Rps6a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps6b - x_Rps6b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps6b - x_Rps6b_Rps31 >= 0 \select node if connected (outgoing)
y_Rps7a - x_Rps7a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps7a - x_Rps7a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps8a - x_Rps8a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps8a - x_Rps8a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps8b - x_Rps8b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps8b - x_Rps8b_Rps31 >= 0 \select node if connected (outgoing)
y_Rps9a - x_Rps9a_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps9a - x_Rps9a_Rps31 >= 0 \select node if connected (outgoing)
y_Rps9b - x_Rps9b_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps9b - x_Rps9b_Rps31 >= 0 \select node if connected (outgoing)
y_Sui1 - x_Sui1_Rps31 >= 0 \select node if connected (outgoing)
y_Tif11 - x_Tif11_Rps31 >= 0 \select node if connected (outgoing)
y_Tif35 - x_Tif35_Rps31 >= 0 \select node if connected (outgoing)
y_Tif5 - x_Tif5_Rps31 >= 0 \select node if connected (outgoing)
y_Sup35 - x_Sup35_Dom34 >= 0 \select node if connected (outgoing)
y_Sup35 - x_Sup35_Ssa1 >= 0 \select node if connected (outgoing)
y_Dom34 - x_Dom34_Rps31 >= 0 \select node if connected (outgoing)
y_Ssa1 - x_Ssa1_Cdc48 >= 0 \select node if connected (outgoing)
y_Swe1 - x_Swe1_Chk1 >= 0 \select node if connected (outgoing)
y_Chk1 - x_Chk1_Cdc9 >= 0 \select node if connected (outgoing)
y_Cdc9 - x_Cdc9_Rad23 >= 0 \select node if connected (outgoing)
y_Act1 - x_Cyr1_Act1 >= 0 \select node if connected (incoming)
y_Cof1 - x_Act1_Cof1 >= 0 \select node if connected (incoming)
y_Cof1 - x_Srv2_Cof1 >= 0 \select node if connected (incoming)
y_Srv2 - x_Cyr1_Srv2 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Duf1_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps31_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Ede1_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Ent1_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Ent2_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Vps27_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Asc1_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl10_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl17a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl17b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl19b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl25_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl28_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl2a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl2b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl3_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl30_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl32_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl34a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl35a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl35b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl38_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rpl5_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps0a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps0b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps10a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps11a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps11b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps12_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps13_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps14a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps14b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps15_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps16a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps16b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps17a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps17b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps18b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps19a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps1a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps1b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps2_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps20_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps23a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps23b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps24a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps24b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps25a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps26a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps26b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps27a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps28a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps28b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps29a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps3_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps5_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps6b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps7a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps8a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps8b_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps9a_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps9b_Ubi4 >= 0 \select node if connected (incoming)
y_Cdc48 - x_Ubi4_Cdc48 >= 0 \select node if connected (incoming)
y_Cdc48 - x_Rps31_Cdc48 >= 0 \select node if connected (incoming)
y_Cdc48 - x_Ssa1_Cdc48 >= 0 \select node if connected (incoming)
y_Ddi1 - x_Ubi4_Ddi1 >= 0 \select node if connected (incoming)
y_Dsk2 - x_Ubi4_Dsk2 >= 0 \select node if connected (incoming)
y_Dsk2 - x_Rps31_Dsk2 >= 0 \select node if connected (incoming)
y_Rpn1 - x_Ubi4_Rpn1 >= 0 \select node if connected (incoming)
y_Rpn1 - x_Rps31_Rpn1 >= 0 \select node if connected (incoming)
y_Rpn10 - x_Ubi4_Rpn10 >= 0 \select node if connected (incoming)
y_Rpn10 - x_Rps31_Rpn10 >= 0 \select node if connected (incoming)
y_Rpn11 - x_Ubi4_Rpn11 >= 0 \select node if connected (incoming)
y_Rpn11 - x_Rps31_Rpn11 >= 0 \select node if connected (incoming)
y_Rpn11 - x_Nip1_Rpn11 >= 0 \select node if connected (incoming)
y_Rpn11 - x_Rps5_Rpn11 >= 0 \select node if connected (incoming)
y_Rps31 - x_Ubi4_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Ede1_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Vps27_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Fun12_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Asc1_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Gcd11_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Nip1_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Prt1_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpg1_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl10_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl17a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl17b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl19b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl25_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl28_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl2a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl2b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl3_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl30_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl32_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl34a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl35a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl35b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl37a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl38_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl5_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rpl7a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps0a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps0b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps10a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps11a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps11b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps12_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps13_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps14a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps14b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps15_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps16a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps16b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps17a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps17b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps18a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps18b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps19a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps1a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps1b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps2_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps20_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps22a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps23a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps23b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps24a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps24b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps25a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps26a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps26b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps27a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps28a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps28b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps29a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps3_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps30a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps30b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps5_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps6a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps6b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps7a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps8a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps8b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps9a_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rps9b_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Sui1_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Tif11_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Tif35_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Tif5_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Dom34_Rps31 >= 0 \select node if connected (incoming)
y_Rpt1 - x_Ubi4_Rpt1 >= 0 \select node if connected (incoming)
y_Rpt1 - x_Rps31_Rpt1 >= 0 \select node if connected (incoming)
y_Rpt5 - x_Ubi4_Rpt5 >= 0 \select node if connected (incoming)
y_Rpt5 - x_Rps31_Rpt5 >= 0 \select node if connected (incoming)
y_Rpt6 - x_Ubi4_Rpt6 >= 0 \select node if connected (incoming)
y_Rpt6 - x_Rps31_Rpt6 >= 0 \select node if connected (incoming)
y_Ufd2 - x_Ubi4_Ufd2 >= 0 \select node if connected (incoming)
y_Ent1 - x_Ede1_Ent1 >= 0 \select node if connected (incoming)
y_Chc1 - x_Ent1_Chc1 >= 0 \select node if connected (incoming)
y_Chc1 - x_Ent2_Chc1 >= 0 \select node if connected (incoming)
y_Chc1 - x_Vps27_Chc1 >= 0 \select node if connected (incoming)
y_Ent2 - x_Ede1_Ent2 >= 0 \select node if connected (incoming)
y_Vps27 - x_Ede1_Vps27 >= 0 \select node if connected (incoming)
y_Asc1 - x_Fun12_Asc1 >= 0 \select node if connected (incoming)
y_Gcd11 - x_Fun12_Gcd11 >= 0 \select node if connected (incoming)
y_Nip1 - x_Fun12_Nip1 >= 0 \select node if connected (incoming)
y_Prt1 - x_Fun12_Prt1 >= 0 \select node if connected (incoming)
y_Rpg1 - x_Fun12_Rpg1 >= 0 \select node if connected (incoming)
y_Rpl10 - x_Fun12_Rpl10 >= 0 \select node if connected (incoming)
y_Rpl17a - x_Fun12_Rpl17a >= 0 \select node if connected (incoming)
y_Rpl17b - x_Fun12_Rpl17b >= 0 \select node if connected (incoming)
y_Rpl19b - x_Fun12_Rpl19b >= 0 \select node if connected (incoming)
y_Rpl25 - x_Fun12_Rpl25 >= 0 \select node if connected (incoming)
y_Rpl28 - x_Fun12_Rpl28 >= 0 \select node if connected (incoming)
y_Rpl2a - x_Fun12_Rpl2a >= 0 \select node if connected (incoming)
y_Rpl2b - x_Fun12_Rpl2b >= 0 \select node if connected (incoming)
y_Rpl3 - x_Fun12_Rpl3 >= 0 \select node if connected (incoming)
y_Rpl30 - x_Fun12_Rpl30 >= 0 \select node if connected (incoming)
y_Rpl32 - x_Fun12_Rpl32 >= 0 \select node if connected (incoming)
y_Rpl34a - x_Fun12_Rpl34a >= 0 \select node if connected (incoming)
y_Rpl35a - x_Fun12_Rpl35a >= 0 \select node if connected (incoming)
y_Rpl35b - x_Fun12_Rpl35b >= 0 \select node if connected (incoming)
y_Rpl37a - x_Fun12_Rpl37a >= 0 \select node if connected (incoming)
y_Rpl38 - x_Fun12_Rpl38 >= 0 \select node if connected (incoming)
y_Rpl5 - x_Fun12_Rpl5 >= 0 \select node if connected (incoming)
y_Rpl7a - x_Fun12_Rpl7a >= 0 \select node if connected (incoming)
y_Rps0a - x_Fun12_Rps0a >= 0 \select node if connected (incoming)
y_Rps0b - x_Fun12_Rps0b >= 0 \select node if connected (incoming)
y_Rps10a - x_Fun12_Rps10a >= 0 \select node if connected (incoming)
y_Rps11a - x_Fun12_Rps11a >= 0 \select node if connected (incoming)
y_Rps11b - x_Fun12_Rps11b >= 0 \select node if connected (incoming)
y_Rps12 - x_Fun12_Rps12 >= 0 \select node if connected (incoming)
y_Rps13 - x_Fun12_Rps13 >= 0 \select node if connected (incoming)
y_Rps14a - x_Fun12_Rps14a >= 0 \select node if connected (incoming)
y_Rps14b - x_Fun12_Rps14b >= 0 \select node if connected (incoming)
y_Rps15 - x_Fun12_Rps15 >= 0 \select node if connected (incoming)
y_Rps16a - x_Fun12_Rps16a >= 0 \select node if connected (incoming)
y_Rps16b - x_Fun12_Rps16b >= 0 \select node if connected (incoming)
y_Rps17a - x_Fun12_Rps17a >= 0 \select node if connected (incoming)
y_Rps17b - x_Fun12_Rps17b >= 0 \select node if connected (incoming)
y_Rps18a - x_Fun12_Rps18a >= 0 \select node if connected (incoming)
y_Rps18b - x_Fun12_Rps18b >= 0 \select node if connected (incoming)
y_Rps19a - x_Fun12_Rps19a >= 0 \select node if connected (incoming)
y_Rps1a - x_Fun12_Rps1a >= 0 \select node if connected (incoming)
y_Rps1b - x_Fun12_Rps1b >= 0 \select node if connected (incoming)
y_Rps2 - x_Fun12_Rps2 >= 0 \select node if connected (incoming)
y_Rps20 - x_Fun12_Rps20 >= 0 \select node if connected (incoming)
y_Rps22a - x_Fun12_Rps22a >= 0 \select node if connected (incoming)
y_Rps23a - x_Fun12_Rps23a >= 0 \select node if connected (incoming)
y_Rps23b - x_Fun12_Rps23b >= 0 \select node if connected (incoming)
y_Rps24a - x_Fun12_Rps24a >= 0 \select node if connected (incoming)
y_Rps24b - x_Fun12_Rps24b >= 0 \select node if connected (incoming)
y_Rps25a - x_Fun12_Rps25a >= 0 \select node if connected (incoming)
y_Rps26a - x_Fun12_Rps26a >= 0 \select node if connected (incoming)
y_Rps26b - x_Fun12_Rps26b >= 0 \select node if connected (incoming)
y_Rps27a - x_Fun12_Rps27a >= 0 \select node if connected (incoming)
y_Rps28a - x_Fun12_Rps28a >= 0 \select node if connected (incoming)
y_Rps28b - x_Fun12_Rps28b >= 0 \select node if connected (incoming)
y_Rps29a - x_Fun12_Rps29a >= 0 \select node if connected (incoming)
y_Rps3 - x_Fun12_Rps3 >= 0 \select node if connected (incoming)
y_Rps30a - x_Fun12_Rps30a >= 0 \select node if connected (incoming)
y_Rps30b - x_Fun12_Rps30b >= 0 \select node if connected (incoming)
y_Rps5 - x_Fun12_Rps5 >= 0 \select node if connected (incoming)
y_Rps6a - x_Fun12_Rps6a >= 0 \select node if connected (incoming)
y_Rps6b - x_Fun12_Rps6b >= 0 \select node if connected (incoming)
y_Rps7a - x_Fun12_Rps7a >= 0 \select node if connected (incoming)
y_Rps8a - x_Fun12_Rps8a >= 0 \select node if connected (incoming)
y_Rps8b - x_Fun12_Rps8b >= 0 \select node if connected (incoming)
y_Rps9a - x_Fun12_Rps9a >= 0 \select node if connected (incoming)
y_Rps9b - x_Fun12_Rps9b >= 0 \select node if connected (incoming)
y_Sui1 - x_Fun12_Sui1 >= 0 \select node if connected (incoming)
y_Tif11 - x_Fun12_Tif11 >= 0 \select node if connected (incoming)
y_Tif35 - x_Fun12_Tif35 >= 0 \select node if connected (incoming)
y_Tif5 - x_Fun12_Tif5 >= 0 \select node if connected (incoming)
y_Dom34 - x_Sup35_Dom34 >= 0 \select node if connected (incoming)
y_Ssa1 - x_Sup35_Ssa1 >= 0 \select node if connected (incoming)
y_Chk1 - x_Swe1_Chk1 >= 0 \select node if connected (incoming)
y_Cdc9 - x_Chk1_Cdc9 >= 0 \select node if connected (incoming)
y_Prm15 - x_Cof1_Prm15 >= 0 \select node if connected (incoming)
y_Rad23 - x_Ubi4_Rad23 >= 0 \select node if connected (incoming)
y_Rad23 - x_Cdc48_Rad23 >= 0 \select node if connected (incoming)
y_Rad23 - x_Ddi1_Rad23 >= 0 \select node if connected (incoming)
y_Rad23 - x_Dsk2_Rad23 >= 0 \select node if connected (incoming)
y_Rad23 - x_Rpn1_Rad23 >= 0 \select node if connected (incoming)
y_Rad23 - x_Rpn10_Rad23 >= 0 \select node if connected (incoming)
y_Rad23 - x_Rpn11_Rad23 >= 0 \select node if connected (incoming)
y_Rad23 - x_Rps31_Rad23 >= 0 \select node if connected (incoming)
y_Rad23 - x_Rpt1_Rad23 >= 0 \select node if connected (incoming)
y_Rad23 - x_Rpt5_Rad23 >= 0 \select node if connected (incoming)
y_Rad23 - x_Rpt6_Rad23 >= 0 \select node if connected (incoming)
y_Rad23 - x_Ufd2_Rad23 >= 0 \select node if connected (incoming)
y_Rad23 - x_Cdc9_Rad23 >= 0 \select node if connected (incoming)
y_Ent5 - x_Chc1_Ent5 >= 0 \select node if connected (incoming)
x_Ubi4_Rps31 + x_Rps31_Ubi4 <= 1 \avoid 2-node cycles

binary
y_Cyr1
y_Act1
y_Cof1
y_Srv2
y_Duf1
y_Ubi4
y_Cdc48
y_Ddi1
y_Dsk2
y_Rpn1
y_Rpn10
y_Rpn11
y_Rps31
y_Rpt1
y_Rpt5
y_Rpt6
y_Ufd2
y_Ede1
y_Ent1
y_Chc1
y_Ent2
y_Vps27
y_Fun12
y_Asc1
y_Gcd11
y_Nip1
y_Prt1
y_Rpg1
y_Rpl10
y_Rpl17a
y_Rpl17b
y_Rpl19b
y_Rpl25
y_Rpl28
y_Rpl2a
y_Rpl2b
y_Rpl3
y_Rpl30
y_Rpl32
y_Rpl34a
y_Rpl35a
y_Rpl35b
y_Rpl37a
y_Rpl38
y_Rpl5
y_Rpl7a
y_Rps0a
y_Rps0b
y_Rps10a
y_Rps11a
y_Rps11b
y_Rps12
y_Rps13
y_Rps14a
y_Rps14b
y_Rps15
y_Rps16a
y_Rps16b
y_Rps17a
y_Rps17b
y_Rps18a
y_Rps18b
y_Rps19a
y_Rps1a
y_Rps1b
y_Rps2
y_Rps20
y_Rps22a
y_Rps23a
y_Rps23b
y_Rps24a
y_Rps24b
y_Rps25a
y_Rps26a
y_Rps26b
y_Rps27a
y_Rps28a
y_Rps28b
y_Rps29a
y_Rps3
y_Rps30a
y_Rps30b
y_Rps5
y_Rps6a
y_Rps6b
y_Rps7a
y_Rps8a
y_Rps8b
y_Rps9a
y_Rps9b
y_Sui1
y_Tif11
y_Tif35
y_Tif5
y_Sup35
y_Dom34
y_Ssa1
y_Swe1
y_Chk1
y_Cdc9
y_Prm15
y_Rad23
y_Ent5
x_Cyr1_Act1
x_Cyr1_Srv2
x_Act1_Cof1
x_Cof1_Prm15
x_Srv2_Cof1
x_Duf1_Ubi4
x_Ubi4_Cdc48
x_Ubi4_Ddi1
x_Ubi4_Dsk2
x_Ubi4_Rpn1
x_Ubi4_Rpn10
x_Ubi4_Rpn11
x_Ubi4_Rps31
x_Ubi4_Rpt1
x_Ubi4_Rpt5
x_Ubi4_Rpt6
x_Ubi4_Ufd2
x_Ubi4_Rad23
x_Cdc48_Rad23
x_Ddi1_Rad23
x_Dsk2_Rad23
x_Rpn1_Rad23
x_Rpn10_Rad23
x_Rpn11_Rad23
x_Rps31_Ubi4
x_Rps31_Cdc48
x_Rps31_Dsk2
x_Rps31_Rpn1
x_Rps31_Rpn10
x_Rps31_Rpn11
x_Rps31_Rpt1
x_Rps31_Rpt5
x_Rps31_Rpt6
x_Rps31_Rad23
x_Rpt1_Rad23
x_Rpt5_Rad23
x_Rpt6_Rad23
x_Ufd2_Rad23
x_Ede1_Ubi4
x_Ede1_Rps31
x_Ede1_Ent1
x_Ede1_Ent2
x_Ede1_Vps27
x_Ent1_Ubi4
x_Ent1_Chc1
x_Chc1_Ent5
x_Ent2_Ubi4
x_Ent2_Chc1
x_Vps27_Ubi4
x_Vps27_Rps31
x_Vps27_Chc1
x_Fun12_Rps31
x_Fun12_Asc1
x_Fun12_Gcd11
x_Fun12_Nip1
x_Fun12_Prt1
x_Fun12_Rpg1
x_Fun12_Rpl10
x_Fun12_Rpl17a
x_Fun12_Rpl17b
x_Fun12_Rpl19b
x_Fun12_Rpl25
x_Fun12_Rpl28
x_Fun12_Rpl2a
x_Fun12_Rpl2b
x_Fun12_Rpl3
x_Fun12_Rpl30
x_Fun12_Rpl32
x_Fun12_Rpl34a
x_Fun12_Rpl35a
x_Fun12_Rpl35b
x_Fun12_Rpl37a
x_Fun12_Rpl38
x_Fun12_Rpl5
x_Fun12_Rpl7a
x_Fun12_Rps0a
x_Fun12_Rps0b
x_Fun12_Rps10a
x_Fun12_Rps11a
x_Fun12_Rps11b
x_Fun12_Rps12
x_Fun12_Rps13
x_Fun12_Rps14a
x_Fun12_Rps14b
x_Fun12_Rps15
x_Fun12_Rps16a
x_Fun12_Rps16b
x_Fun12_Rps17a
x_Fun12_Rps17b
x_Fun12_Rps18a
x_Fun12_Rps18b
x_Fun12_Rps19a
x_Fun12_Rps1a
x_Fun12_Rps1b
x_Fun12_Rps2
x_Fun12_Rps20
x_Fun12_Rps22a
x_Fun12_Rps23a
x_Fun12_Rps23b
x_Fun12_Rps24a
x_Fun12_Rps24b
x_Fun12_Rps25a
x_Fun12_Rps26a
x_Fun12_Rps26b
x_Fun12_Rps27a
x_Fun12_Rps28a
x_Fun12_Rps28b
x_Fun12_Rps29a
x_Fun12_Rps3
x_Fun12_Rps30a
x_Fun12_Rps30b
x_Fun12_Rps5
x_Fun12_Rps6a
x_Fun12_Rps6b
x_Fun12_Rps7a
x_Fun12_Rps8a
x_Fun12_Rps8b
x_Fun12_Rps9a
x_Fun12_Rps9b
x_Fun12_Sui1
x_Fun12_Tif11
x_Fun12_Tif35
x_Fun12_Tif5
x_Asc1_Ubi4
x_Asc1_Rps31
x_Gcd11_Rps31
x_Nip1_Rpn11
x_Nip1_Rps31
x_Prt1_Rps31
x_Rpg1_Rps31
x_Rpl10_Ubi4
x_Rpl10_Rps31
x_Rpl17a_Ubi4
x_Rpl17a_Rps31
x_Rpl17b_Ubi4
x_Rpl17b_Rps31
x_Rpl19b_Ubi4
x_Rpl19b_Rps31
x_Rpl25_Ubi4
x_Rpl25_Rps31
x_Rpl28_Ubi4
x_Rpl28_Rps31
x_Rpl2a_Ubi4
x_Rpl2a_Rps31
x_Rpl2b_Ubi4
x_Rpl2b_Rps31
x_Rpl3_Ubi4
x_Rpl3_Rps31
x_Rpl30_Ubi4
x_Rpl30_Rps31
x_Rpl32_Ubi4
x_Rpl32_Rps31
x_Rpl34a_Ubi4
x_Rpl34a_Rps31
x_Rpl35a_Ubi4
x_Rpl35a_Rps31
x_Rpl35b_Ubi4
x_Rpl35b_Rps31
x_Rpl37a_Rps31
x_Rpl38_Ubi4
x_Rpl38_Rps31
x_Rpl5_Ubi4
x_Rpl5_Rps31
x_Rpl7a_Rps31
x_Rps0a_Ubi4
x_Rps0a_Rps31
x_Rps0b_Ubi4
x_Rps0b_Rps31
x_Rps10a_Ubi4
x_Rps10a_Rps31
x_Rps11a_Ubi4
x_Rps11a_Rps31
x_Rps11b_Ubi4
x_Rps11b_Rps31
x_Rps12_Ubi4
x_Rps12_Rps31
x_Rps13_Ubi4
x_Rps13_Rps31
x_Rps14a_Ubi4
x_Rps14a_Rps31
x_Rps14b_Ubi4
x_Rps14b_Rps31
x_Rps15_Ubi4
x_Rps15_Rps31
x_Rps16a_Ubi4
x_Rps16a_Rps31
x_Rps16b_Ubi4
x_Rps16b_Rps31
x_Rps17a_Ubi4
x_Rps17a_Rps31
x_Rps17b_Ubi4
x_Rps17b_Rps31
x_Rps18a_Rps31
x_Rps18b_Ubi4
x_Rps18b_Rps31
x_Rps19a_Ubi4
x_Rps19a_Rps31
x_Rps1a_Ubi4
x_Rps1a_Rps31
x_Rps1b_Ubi4
x_Rps1b_Rps31
x_Rps2_Ubi4
x_Rps2_Rps31
x_Rps20_Ubi4
x_Rps20_Rps31
x_Rps22a_Rps31
x_Rps23a_Ubi4
x_Rps23a_Rps31
x_Rps23b_Ubi4
x_Rps23b_Rps31
x_Rps24a_Ubi4
x_Rps24a_Rps31
x_Rps24b_Ubi4
x_Rps24b_Rps31
x_Rps25a_Ubi4
x_Rps25a_Rps31
x_Rps26a_Ubi4
x_Rps26a_Rps31
x_Rps26b_Ubi4
x_Rps26b_Rps31
x_Rps27a_Ubi4
x_Rps27a_Rps31
x_Rps28a_Ubi4
x_Rps28a_Rps31
x_Rps28b_Ubi4
x_Rps28b_Rps31
x_Rps29a_Ubi4
x_Rps29a_Rps31
x_Rps3_Ubi4
x_Rps3_Rps31
x_Rps30a_Rps31
x_Rps30b_Rps31
x_Rps5_Ubi4
x_Rps5_Rpn11
x_Rps5_Rps31
x_Rps6a_Rps31
x_Rps6b_Ubi4
x_Rps6b_Rps31
x_Rps7a_Ubi4
x_Rps7a_Rps31
x_Rps8a_Ubi4
x_Rps8a_Rps31
x_Rps8b_Ubi4
x_Rps8b_Rps31
x_Rps9a_Ubi4
x_Rps9a_Rps31
x_Rps9b_Ubi4
x_Rps9b_Rps31
x_Sui1_Rps31
x_Tif11_Rps31
x_Tif35_Rps31
x_Tif5_Rps31
x_Sup35_Dom34
x_Sup35_Ssa1
x_Dom34_Rps31
x_Ssa1_Cdc48
x_Swe1_Chk1
x_Chk1_Cdc9
x_Cdc9_Rad23

End