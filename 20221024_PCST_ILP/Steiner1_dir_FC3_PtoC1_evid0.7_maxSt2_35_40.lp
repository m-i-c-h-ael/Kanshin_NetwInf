\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 0.084 x_Air2_Mtr4 + 0.2 x_Mtr4_Rps11a + 0.032 x_Rps11a_Fun12 + 0.2 x_Mtr4_Rps11b + 0.032 x_Rps11b_Fun12 + 0.196 x_Mtr4_Rps13 + 0.105 x_Rps13_Fun12 + 0.2 x_Mtr4_Rps14a + 0.051 x_Rps14a_Fun12 + 0.2 x_Mtr4_Rps16a + 0.032 x_Rps16a_Fun12 + 0.2 x_Mtr4_Rps16b + 0.032 x_Rps16b_Fun12 + 0.2 x_Mtr4_Rps1a + 0.032 x_Rps1a_Fun12 + 0.188 x_Mtr4_Rps23a + 0.032 x_Rps23a_Fun12 + 0.188 x_Mtr4_Rps23b + 0.032 x_Rps23b_Fun12 + 0.2 x_Mtr4_Rps24a + 0.032 x_Rps24a_Fun12 + 0.2 x_Mtr4_Rps24b + 0.032 x_Rps24b_Fun12 + 0.2 x_Mtr4_Rps27a + 0.032 x_Rps27a_Fun12 + 0.2 x_Mtr4_Rps28a + 0.032 x_Rps28a_Fun12 + 0.2 x_Mtr4_Rps5 + 0.00800000000000001 x_Rps5_Fun12 + 0.2 x_Mtr4_Rps6a + 0.021 x_Rps6a_Fun12 + 0.2 x_Mtr4_Rps6b + 0.021 x_Rps6b_Fun12 + 0.196 x_Mtr4_Rps7a + 0.032 x_Rps7a_Fun12 + 0.195 x_Mtr4_Rps8a + 0.032 x_Rps8a_Fun12 + 0.2 x_Mtr4_Rps8b + 0.032 x_Rps8b_Fun12 + 0.2 x_Mtr4_Rps9a + 0.032 x_Rps9a_Fun12 + 0.277 x_Air2_Npl3 + 0.298 x_Npl3_Rpo21 + 0.021 x_Rpo21_Rtf1 + 0.181 x_Art5_Rsp5 + 0.0580000000000001 x_Rsp5_Rpo21 + 0.00600000000000001 x_Rsp5_Rps31 + 0.172 x_Rps31_Ede1 + 0.00800000000000001 x_Rps31_Fun12 + 0.013 x_Rsp5_Ubi4 + 0.131 x_Ubi4_Duf1 + 0.028 x_Ubi4_Ede1 + 0.26 x_Ask10_Fps1 + 0.141 x_Fps1_Hog1 + 0.001 x_Hog1_Sko1 + 0.248 x_Fps1_Slt2 + 0.168 x_Slt2_Sko1 + 0.2 x_Bdp1_Rpb10 + 0.265 x_Rpb10_Cdc73 + 0.003 x_Cdc73_Rtf1 + 0.001 x_Rpb10_Rpb2 + 0.214 x_Rpb2_Rtf1 + 0.001 x_Rpb10_Rpb3 + 0.26 x_Rpb3_Rtf1 + 0.001 x_Rpb10_Rpo21 + 0.024 x_Rpb10_Spt4 + 0.218 x_Spt4_Rtf1 + 0.031 x_Rpb10_Spt5 + 0.013 x_Spt5_Rtf1 + 0.2 x_Bdp1_Rpb5 + 0.268 x_Rpb5_Cdc73 + 0.292 x_Rpb5_Leo1 + 0.003 x_Leo1_Rtf1 + 0.001 x_Rpb5_Rpb2 + 0.001 x_Rpb5_Rpb3 + 0.001 x_Rpb5_Rpo21 + 0.013 x_Rpb5_Spt4 + 0.00700000000000001 x_Rpb5_Spt5 + 0.038 x_Rpb5_Spt6 + 0.27 x_Spt6_Rtf1 + 0.2 x_Bdp1_Rpb8 + 0.001 x_Rpb8_Rpb2 + 0.001 x_Rpb8_Rpb3 + 0.001 x_Rpb8_Rpo21 + 0.226 x_Rpb8_Spt4 + 0.027 x_Rpb8_Spt5 + 0.2 x_Bdp1_Rpc10 + 0.261 x_Rpc10_Cdc73 + 0.003 x_Rpc10_Rpb2 + 0.001 x_Rpc10_Rpb3 + 0.001 x_Rpc10_Rpo21 + 0.024 x_Rpc10_Spt4 + 0.024 x_Rpc10_Spt5 + 0.196 x_Bdp1_Rpo26 + 0.261 x_Rpo26_Cdc73 + 0.253 x_Rpo26_Paf1 + 0.003 x_Paf1_Rtf1 + 0.001 x_Rpo26_Rpb2 + 0.001 x_Rpo26_Rpb3 + 0.001 x_Rpo26_Rpo21 + 0.011 x_Rpo26_Spt4 + 0.011 x_Rpo26_Spt5 + 0.117 x_Bdp1_Rpo31 + 0.037 x_Rpo31_Rpb3 + 0.038 x_Rpo31_Rpo21 + 0.00900000000000001 x_Bdp1_Spt15 + 0.124 x_Spt15_Rpb2 + 0.069 x_Spt15_Rpb3 + 0.002 x_Spt15_Rpo21 + 0.297 x_Chc1_Apl1 + 0.00900000000000001 x_Apl1_Apl3 + 0.031 x_Apl3_Ede1 + 0.081 x_Chc1_Apl2 + 0.038 x_Apl2_Apl3 + 0.056 x_Chc1_Ent1 + 0.152 x_Ent1_Apl3 + 0.024 x_Ent1_Ede1 + 0.177 x_Ent1_Ubi4 + 0.1 x_Chc1_Ent2 + 0.152 x_Ent2_Apl3 + 0.026 x_Ent2_Ede1 + 0.286 x_Ent2_Ubi4 + 0.278 x_Chc1_Las17 + 0.004 x_Las17_Act1 + 0.214 x_Act1_Cyr1 + 0.141 x_Chc1_Sla2 + 0.273 x_Sla2_Act1 + 0.102 x_Sla2_Ent1 + 0.00600000000000001 x_Chc1_Swa2 + 0.217 x_Swa2_Ubi4 + 0.0570000000000001 x_Chc1_Vps27 + 0.0659999999999999 x_Vps27_Ede1 + 0.106 x_Vps27_Rps31 + 0.001 x_Vps27_Ubi4 + 0.221 x_Ecm21_Rsp5 + 0.093 x_Ldb19_Rsp5 + 0.036 x_Msg5_Fus3 + 0.165 x_Fus3_Gpa1 + 0.254 x_Gpa1_Cyr1 + 0.152 x_Fus3_Kss1 + 0.268 x_Kss1_Sko1 + 0.235 x_Fus3_Sko1 + 0.215 x_Msg5_Kss1 + 0.152 x_Kss1_Fus3 + 0.232 x_Kss1_Hog1 + 0.077 x_Msg5_Slt2 + 0.208 x_Slt2_Cdc28 + 0.003 x_Cdc28_Swe1 + 0.214 x_Slt2_Hog1 + 0.037 x_Sec16_Apl3 + 0.152 x_Apl3_Ent1 + 0.152 x_Apl3_Ent2 + 0.038 x_Sec16_Apm4 + 0.005 x_Apm4_Apl3 + 0.198 x_Sla1_Abp1 + 0.102 x_Abp1_Act1 + 0.077 x_Abp1_Srv2 + 0.028 x_Srv2_Cyr1 + 0.084 x_Sla1_Las17 + 0.2 x_Sla1_Rps31 + 0.005 x_Rps31_Asc1 + 0.00800000000000001 x_Asc1_Fun12 + 0.044 x_Rps31_Dom34 + 0.29 x_Dom34_Sup35 + 0.216 x_Rps31_Gcd11 + 0.271 x_Gcd11_Fun12 + 0.043 x_Rps31_Nip1 + 0.00900000000000001 x_Nip1_Fun12 + 0.078 x_Rps31_Prt1 + 0.00600000000000001 x_Prt1_Fun12 + 0.039 x_Rps31_Rpg1 + 0.016 x_Rpg1_Fun12 + 0.005 x_Rps31_Rpl10 + 0.00800000000000001 x_Rpl10_Fun12 + 0.00800000000000001 x_Rps31_Rpl17a + 0.032 x_Rpl17a_Fun12 + 0.037 x_Rps31_Rpl17b + 0.151 x_Rpl17b_Fun12 + 0.00800000000000001 x_Rps31_Rpl19b + 0.032 x_Rpl19b_Fun12 + 0.005 x_Rps31_Rpl25 + 0.00800000000000001 x_Rpl25_Fun12 + 0.005 x_Rps31_Rpl28 + 0.00800000000000001 x_Rpl28_Fun12 + 0.00800000000000001 x_Rps31_Rpl2a + 0.024 x_Rpl2a_Fun12 + 0.00800000000000001 x_Rps31_Rpl2b + 0.043 x_Rpl2b_Fun12 + 0.005 x_Rps31_Rpl3 + 0.13 x_Rpl3_Fun12 + 0.005 x_Rps31_Rpl30 + 0.00800000000000001 x_Rpl30_Fun12 + 0.005 x_Rps31_Rpl32 + 0.00800000000000001 x_Rpl32_Fun12 + 0.00800000000000001 x_Rps31_Rpl34a + 0.032 x_Rpl34a_Fun12 + 0.00800000000000001 x_Rps31_Rpl35a + 0.032 x_Rpl35a_Fun12 + 0.00800000000000001 x_Rps31_Rpl35b + 0.032 x_Rpl35b_Fun12 + 0.00800000000000001 x_Rps31_Rpl37a + 0.032 x_Rpl37a_Fun12 + 0.005 x_Rps31_Rpl38 + 0.00800000000000001 x_Rpl38_Fun12 + 0.005 x_Rps31_Rpl5 + 0.00800000000000001 x_Rpl5_Fun12 + 0.00800000000000001 x_Rps31_Rpl7a + 0.032 x_Rpl7a_Fun12 + 0.00800000000000001 x_Rps31_Rps0a + 0.032 x_Rps0a_Fun12 + 0.00800000000000001 x_Rps31_Rps0b + 0.032 x_Rps0b_Fun12 + 0.00800000000000001 x_Rps31_Rps10a + 0.032 x_Rps10a_Fun12 + 0.00800000000000001 x_Rps31_Rps11a + 0.00800000000000001 x_Rps31_Rps11b + 0.005 x_Rps31_Rps12 + 0.00900000000000001 x_Rps12_Fun12 + 0.005 x_Rps31_Rps13 + 0.00800000000000001 x_Rps31_Rps14a + 0.00800000000000001 x_Rps31_Rps14b + 0.051 x_Rps14b_Fun12 + 0.005 x_Rps31_Rps15 + 0.00800000000000001 x_Rps15_Fun12 + 0.00800000000000001 x_Rps31_Rps16a + 0.00800000000000001 x_Rps31_Rps16b + 0.00800000000000001 x_Rps31_Rps17a + 0.032 x_Rps17a_Fun12 + 0.00800000000000001 x_Rps31_Rps17b + 0.032 x_Rps17b_Fun12 + 0.00800000000000001 x_Rps31_Rps18a + 0.032 x_Rps18a_Fun12 + 0.00800000000000001 x_Rps31_Rps18b + 0.032 x_Rps18b_Fun12 + 0.00800000000000001 x_Rps31_Rps19a + 0.032 x_Rps19a_Fun12 + 0.00800000000000001 x_Rps31_Rps1a + 0.00800000000000001 x_Rps31_Rps1b + 0.032 x_Rps1b_Fun12 + 0.005 x_Rps31_Rps2 + 0.104 x_Rps2_Fun12 + 0.005 x_Rps31_Rps20 + 0.00800000000000001 x_Rps20_Fun12 + 0.00800000000000001 x_Rps31_Rps22a + 0.045 x_Rps22a_Fun12 + 0.00800000000000001 x_Rps31_Rps23a + 0.00800000000000001 x_Rps31_Rps23b + 0.00800000000000001 x_Rps31_Rps24a + 0.00800000000000001 x_Rps31_Rps24b + 0.00800000000000001 x_Rps31_Rps25a + 0.16 x_Rps25a_Fun12 + 0.00800000000000001 x_Rps31_Rps26a + 0.032 x_Rps26a_Fun12 + 0.00800000000000001 x_Rps31_Rps26b + 0.032 x_Rps26b_Fun12 + 0.00800000000000001 x_Rps31_Rps27a + 0.00800000000000001 x_Rps31_Rps28a + 0.00800000000000001 x_Rps31_Rps28b + 0.032 x_Rps28b_Fun12 + 0.00800000000000001 x_Rps31_Rps29a + 0.032 x_Rps29a_Fun12 + 0.005 x_Rps31_Rps3 + 0.00800000000000001 x_Rps3_Fun12 + 0.00800000000000001 x_Rps31_Rps30a + 0.032 x_Rps30a_Fun12 + 0.00800000000000001 x_Rps31_Rps30b + 0.032 x_Rps30b_Fun12 + 0.005 x_Rps31_Rps5 + 0.00800000000000001 x_Rps31_Rps6a + 0.00800000000000001 x_Rps31_Rps6b + 0.00800000000000001 x_Rps31_Rps7a + 0.00800000000000001 x_Rps31_Rps8a + 0.00800000000000001 x_Rps31_Rps8b + 0.00800000000000001 x_Rps31_Rps9a + 0.00800000000000001 x_Rps31_Rps9b + 0.032 x_Rps9b_Fun12 + 0.00800000000000001 x_Rps31_Sui1 + 0.038 x_Sui1_Fun12 + 0.00600000000000001 x_Rps31_Tif11 + 0.017 x_Tif11_Fun12 + 0.078 x_Rps31_Tif35 + 0.026 x_Tif35_Fun12 + 0.002 x_Rps31_Ubi4 + 0.174 x_Rps31_Tif5 + 0.229 x_Tif5_Fun12 + 0.106 x_Rps31_Vps27 + 0.229 x_Sla1_Rvs167 + 0.191 x_Rvs167_Act1 + 0.225 x_Sla1_Sla2 + 0.052 x_Sla1_Ubi4 + 0.039 x_Ubi4_Asc1 + 0.177 x_Ubi4_Ent1 + 0.286 x_Ubi4_Ent2 + 0.028 x_Ubi4_Rpl10 + 0.221 x_Ubi4_Rpl17a + 0.221 x_Ubi4_Rpl17b + 0.223 x_Ubi4_Rpl19b + 0.039 x_Ubi4_Rpl25 + 0.028 x_Ubi4_Rpl28 + 0.172 x_Ubi4_Rpl2a + 0.236 x_Ubi4_Rpl2b + 0.266 x_Ubi4_Rpl3 + 0.039 x_Ubi4_Rpl30 + 0.028 x_Ubi4_Rpl32 + 0.16 x_Ubi4_Rpl34a + 0.219 x_Ubi4_Rpl35a + 0.219 x_Ubi4_Rpl35b + 0.039 x_Ubi4_Rpl38 + 0.038 x_Ubi4_Rpl5 + 0.174 x_Ubi4_Rps0a + 0.174 x_Ubi4_Rps0b + 0.214 x_Ubi4_Rps10a + 0.212 x_Ubi4_Rps11a + 0.155 x_Ubi4_Rps11b + 0.028 x_Ubi4_Rps12 + 0.028 x_Ubi4_Rps13 + 0.18 x_Ubi4_Rps14a + 0.18 x_Ubi4_Rps14b + 0.038 x_Ubi4_Rps15 + 0.172 x_Ubi4_Rps16a + 0.172 x_Ubi4_Rps16b + 0.151 x_Ubi4_Rps17a + 0.211 x_Ubi4_Rps17b + 0.228 x_Ubi4_Rps18b + 0.185 x_Ubi4_Rps19a + 0.125 x_Ubi4_Rps1a + 0.174 x_Ubi4_Rps1b + 0.038 x_Ubi4_Rps2 + 0.029 x_Ubi4_Rps20 + 0.19 x_Ubi4_Rps23a + 0.19 x_Ubi4_Rps23b + 0.18 x_Ubi4_Rps24a + 0.186 x_Ubi4_Rps24b + 0.206 x_Ubi4_Rps25a + 0.198 x_Ubi4_Rps26a + 0.276 x_Ubi4_Rps26b + 0.156 x_Ubi4_Rps27a + 0.207 x_Ubi4_Rps28a + 0.207 x_Ubi4_Rps28b + 0.208 x_Ubi4_Rps29a + 0.028 x_Ubi4_Rps3 + 0.002 x_Ubi4_Rps31 + 0.108 x_Ubi4_Rps5 + 0.235 x_Ubi4_Rps6b + 0.144 x_Ubi4_Rps7a + 0.209 x_Ubi4_Rps8a + 0.152 x_Ubi4_Rps8b + 0.287 x_Ubi4_Rps9a + 0.287 x_Ubi4_Rps9b + 0.001 x_Ubi4_Vps27 + 0.194 x_Sla1_Vrp1 + 0.224 x_Vrp1_Act1 + 0.00600000000000001 x_Ssk2_Pbs2 + 0.221 x_Pbs2_Fus3 + 0.001 x_Pbs2_Hog1 + 0.212 x_Pbs2_Kss1 + 0.023 x_Pbs2_Mih1 + 0.00700000000000001 x_Mih1_Swe1 + 0.25 x_Bud4_Bud3 + 0.224 x_Bud3_Clb2 + 0.142 x_Clb2_Swe1 + 0.159 x_Igo1_Cdc55 + 0.214 x_Cdc55_Swe1 + 0.181 x_Cdc55_Igo2 + 0.119 x_Syp1_Ede1 + 0.00700000000000001 x_Ypl247c_Yak1 + 0.283 x_Yak1_Rim15 + 0.286 x_Rim15_Msn2 + 0.182 x_Rim15_Sch9 + 0.191 x_Yak1_Tpk1 + 0.223 x_Tpk1_Sch9 - 1.19335864926735 y_Air2 - 0.236003220158082 y_Mtr4 - 0 y_Rps11a - 0.276867820406057 y_Rps11b - 0 y_Rps13 - 0 y_Rps14a - 0 y_Rps16a - 0.0917888101035867 y_Rps16b - 0.181523970358348 y_Rps1a - 0 y_Rps23a - 0 y_Rps23b - 0.181971262328775 y_Rps24a - 0 y_Rps24b - 0 y_Rps27a - 0 y_Rps28a - 0 y_Rps5 - 0 y_Rps6a - 0.129832488722207 y_Rps6b - 0.149603196222064 y_Rps7a - 0.219566469769258 y_Rps8a - 0 y_Rps8b - 0 y_Rps9a - 0.063200435312774 y_Npl3 - 1.8162863958549 y_Rpo21 - 1.85306782505179 y_Art5 - 0 y_Rsp5 - 0.227967035078496 y_Rps31 - 0 y_Ubi4 - 1.43527559798099 y_Ask10 - 0.795291672469519 y_Fps1 - 1.81347223589972 y_Hog1 - 0 y_Slt2 - 1.36560429096284 y_Bdp1 - 0 y_Rpb10 - 0 y_Cdc73 - 0 y_Rpb2 - 0 y_Rpb3 - 0 y_Spt4 - 0.456111625777312 y_Spt5 - 0 y_Rpb5 - 0.157230763432329 y_Leo1 - 0.216056193029292 y_Spt6 - 0 y_Rpb8 - 0 y_Rpc10 - 0 y_Rpo26 - 0.113116830788378 y_Paf1 - 0 y_Rpo31 - 0 y_Spt15 - 1.02568176866035 y_Chc1 - 0.368956473565731 y_Apl1 - 0.425537292960004 y_Apl3 - 0.428048790864964 y_Apl2 - 0.40488716780871 y_Ent1 - 0.520447659565873 y_Ent2 - 1.98110487445054 y_Las17 - 0 y_Act1 - 0.201670154249431 y_Sla2 - 0.234647618370013 y_Swa2 - 0.0884911371317518 y_Vps27 - 1.69786774814898 y_Ecm21 - 1.71435498238203 y_Ldb19 - 2.11713535523598 y_Msg5 - 0 y_Fus3 - 0 y_Gpa1 - 0 y_Kss1 - 0.121462293216658 y_Cdc28 - 1.36340347864178 y_Sec16 - 0 y_Apm4 - 1.22235198568956 y_Sla1 - 0.639220569280991 y_Abp1 - 0 y_Srv2 - 0.297909763660125 y_Asc1 - 0 y_Dom34 - 0.246136646023169 y_Gcd11 - 0 y_Nip1 - 0.397976742486334 y_Prt1 - 0.414174219075414 y_Rpg1 - 0.183683848587923 y_Rpl10 - 0.260470332775798 y_Rpl17a - 0 y_Rpl17b - 0.136248251868258 y_Rpl19b - 0.270664526640099 y_Rpl25 - 0 y_Rpl28 - 0 y_Rpl2a - 0 y_Rpl2b - 0.37064662837935 y_Rpl3 - 0.0916242099948261 y_Rpl30 - 0.0780246825673995 y_Rpl32 - 0 y_Rpl34a - 0 y_Rpl35a - 0.216836386251174 y_Rpl35b - 0 y_Rpl37a - 0 y_Rpl38 - 0.120065169941506 y_Rpl5 - 0.190773283192771 y_Rpl7a - 0 y_Rps0a - 0 y_Rps0b - 0 y_Rps10a - 0 y_Rps12 - 0.21043919262963 y_Rps14b - 0 y_Rps15 - 0 y_Rps17a - 0.20078081333483 y_Rps17b - 0.33165617846531 y_Rps18a - 0 y_Rps18b - 0 y_Rps19a - 0.511183741969515 y_Rps1b - 0 y_Rps2 - 0 y_Rps20 - 0 y_Rps22a - 0 y_Rps25a - 0 y_Rps26a - 0.128719434628373 y_Rps26b - 0 y_Rps28b - 0.092884269459382 y_Rps29a - 0.118468258343378 y_Rps3 - 0 y_Rps30a - 0 y_Rps30b - 0.22725818108385 y_Rps9b - 0.154199226094496 y_Sui1 - 0 y_Tif11 - 0.18472471674679 y_Tif35 - 0.216773820733091 y_Tif5 - 0.2916138569775 y_Rvs167 - 0.890339699745531 y_Vrp1 - 1.15805335341679 y_Ssk2 - 1.7465565764114 y_Pbs2 - 0.116628340365083 y_Mih1 - 1.14836639805159 y_Bud4 - 0.521359498804557 y_Bud3 - 0 y_Clb2 - 0 y_Igo1 - 0.150676667675947 y_Cdc55 - 1.35449142678341 y_Syp1 - 1.26346293940005 y_Ypl247c - 0.33456987011284 y_Yak1 - 1.00036416841318 y_Rim15 - 0 y_Tpk1 - 1.20568332697494 y_Fun12 - 1.26884043634153 y_Rtf1 - 1.39100988950645 y_Ede1 - 0 y_Duf1 - 1.50119873871931 y_Sko1 - 2.06546962941424 y_Cyr1 - 1.21946220627969 y_Swe1 - 1.53949924297249 y_Sup35 - 0 y_Igo2 - 1.22973675788113 y_Msn2 - 1.51300508365901 y_Sch9

subject to
x_Air2_Mtr4 + x_Air2_Npl3 >= 1 \edges from terminal
x_Art5_Rsp5 >= 1 \edges from terminal
x_Ask10_Fps1 >= 1 \edges from terminal
x_Bdp1_Rpb10 + x_Bdp1_Rpb5 + x_Bdp1_Rpb8 + x_Bdp1_Rpc10 + x_Bdp1_Rpo26 + x_Bdp1_Rpo31 + x_Bdp1_Spt15 >= 1 \edges from terminal
x_Chc1_Apl1 + x_Chc1_Apl2 + x_Chc1_Ent1 + x_Chc1_Ent2 + x_Chc1_Las17 + x_Chc1_Sla2 + x_Chc1_Swa2 + x_Chc1_Vps27 >= 1 \edges from terminal
x_Ecm21_Rsp5 >= 1 \edges from terminal
x_Ldb19_Rsp5 >= 1 \edges from terminal
x_Msg5_Slt2 + x_Msg5_Fus3 + x_Msg5_Kss1 >= 1 \edges from terminal
x_Sec16_Apl3 + x_Sec16_Apm4 >= 1 \edges from terminal
x_Sla1_Rps31 + x_Sla1_Ubi4 + x_Sla1_Las17 + x_Sla1_Sla2 + x_Sla1_Abp1 + x_Sla1_Rvs167 + x_Sla1_Vrp1 >= 1 \edges from terminal
x_Ssk2_Pbs2 >= 1 \edges from terminal
x_Bud4_Bud3 >= 1 \edges from terminal
x_Igo1_Cdc55 >= 1 \edges from terminal
x_Syp1_Ede1 >= 1 \edges from terminal
x_Ypl247c_Yak1 >= 1 \edges from terminal
x_Rps11a_Fun12 + x_Rps11b_Fun12 + x_Rps13_Fun12 + x_Rps14a_Fun12 + x_Rps16a_Fun12 + x_Rps16b_Fun12 + x_Rps1a_Fun12 + x_Rps23a_Fun12 + x_Rps23b_Fun12 + x_Rps24a_Fun12 + x_Rps24b_Fun12 + x_Rps27a_Fun12 + x_Rps28a_Fun12 + x_Rps5_Fun12 + x_Rps6a_Fun12 + x_Rps6b_Fun12 + x_Rps7a_Fun12 + x_Rps8a_Fun12 + x_Rps8b_Fun12 + x_Rps9a_Fun12 + x_Rps31_Fun12 + x_Asc1_Fun12 + x_Gcd11_Fun12 + x_Nip1_Fun12 + x_Prt1_Fun12 + x_Rpg1_Fun12 + x_Rpl10_Fun12 + x_Rpl17a_Fun12 + x_Rpl17b_Fun12 + x_Rpl19b_Fun12 + x_Rpl25_Fun12 + x_Rpl28_Fun12 + x_Rpl2a_Fun12 + x_Rpl2b_Fun12 + x_Rpl3_Fun12 + x_Rpl30_Fun12 + x_Rpl32_Fun12 + x_Rpl34a_Fun12 + x_Rpl35a_Fun12 + x_Rpl35b_Fun12 + x_Rpl37a_Fun12 + x_Rpl38_Fun12 + x_Rpl5_Fun12 + x_Rpl7a_Fun12 + x_Rps0a_Fun12 + x_Rps0b_Fun12 + x_Rps10a_Fun12 + x_Rps12_Fun12 + x_Rps14b_Fun12 + x_Rps15_Fun12 + x_Rps17a_Fun12 + x_Rps17b_Fun12 + x_Rps18a_Fun12 + x_Rps18b_Fun12 + x_Rps19a_Fun12 + x_Rps1b_Fun12 + x_Rps2_Fun12 + x_Rps20_Fun12 + x_Rps22a_Fun12 + x_Rps25a_Fun12 + x_Rps26a_Fun12 + x_Rps26b_Fun12 + x_Rps28b_Fun12 + x_Rps29a_Fun12 + x_Rps3_Fun12 + x_Rps30a_Fun12 + x_Rps30b_Fun12 + x_Rps9b_Fun12 + x_Sui1_Fun12 + x_Tif11_Fun12 + x_Tif35_Fun12 + x_Tif5_Fun12 >= 1 \edges to terminal
x_Rpo21_Rtf1 + x_Cdc73_Rtf1 + x_Rpb2_Rtf1 + x_Rpb3_Rtf1 + x_Spt4_Rtf1 + x_Spt5_Rtf1 + x_Leo1_Rtf1 + x_Spt6_Rtf1 + x_Paf1_Rtf1 >= 1 \edges to terminal
x_Rps31_Ede1 + x_Ubi4_Ede1 + x_Apl3_Ede1 + x_Ent1_Ede1 + x_Ent2_Ede1 + x_Vps27_Ede1 + x_Syp1_Ede1 >= 1 \edges to terminal
x_Ubi4_Duf1 >= 1 \edges to terminal
x_Hog1_Sko1 + x_Slt2_Sko1 + x_Fus3_Sko1 + x_Kss1_Sko1 >= 1 \edges to terminal
x_Act1_Cyr1 + x_Gpa1_Cyr1 + x_Srv2_Cyr1 >= 1 \edges to terminal
x_Cdc28_Swe1 + x_Mih1_Swe1 + x_Clb2_Swe1 + x_Cdc55_Swe1 >= 1 \edges to terminal
x_Dom34_Sup35 >= 1 \edges to terminal
x_Cdc55_Igo2 >= 1 \edges to terminal
x_Rim15_Msn2 >= 1 \edges to terminal
x_Rim15_Sch9 + x_Tpk1_Sch9 >= 1 \edges to terminal
y_Air2 >= 1  \terminal
y_Art5 >= 1  \terminal
y_Ask10 >= 1  \terminal
y_Bdp1 >= 1  \terminal
y_Chc1 >= 1  \terminal
y_Ecm21 >= 1  \terminal
y_Ldb19 >= 1  \terminal
y_Msg5 >= 1  \terminal
y_Sec16 >= 1  \terminal
y_Sla1 >= 1  \terminal
y_Ssk2 >= 1  \terminal
y_Bud4 >= 1  \terminal
y_Igo1 >= 1  \terminal
y_Syp1 >= 1  \terminal
y_Ypl247c >= 1  \terminal
y_Fun12 >= 1  \terminal
y_Rtf1 >= 1  \terminal
y_Ede1 >= 1  \terminal
y_Duf1 >= 1  \terminal
y_Sko1 >= 1  \terminal
y_Cyr1 >= 1  \terminal
y_Swe1 >= 1  \terminal
y_Sup35 >= 1  \terminal
y_Igo2 >= 1  \terminal
y_Msn2 >= 1  \terminal
y_Sch9 >= 1  \terminal
y_Mtr4 - x_Mtr4_Rps11a - x_Mtr4_Rps11b - x_Mtr4_Rps13 - x_Mtr4_Rps14a - x_Mtr4_Rps16a - x_Mtr4_Rps16b - x_Mtr4_Rps1a - x_Mtr4_Rps23a - x_Mtr4_Rps23b - x_Mtr4_Rps24a - x_Mtr4_Rps24b - x_Mtr4_Rps27a - x_Mtr4_Rps28a - x_Mtr4_Rps5 - x_Mtr4_Rps6a - x_Mtr4_Rps6b - x_Mtr4_Rps7a - x_Mtr4_Rps8a - x_Mtr4_Rps8b - x_Mtr4_Rps9a <= 0 \Steiner node out
y_Rps11a - x_Rps11a_Fun12 <= 0 \Steiner node out
y_Rps11b - x_Rps11b_Fun12 <= 0 \Steiner node out
y_Rps13 - x_Rps13_Fun12 <= 0 \Steiner node out
y_Rps14a - x_Rps14a_Fun12 <= 0 \Steiner node out
y_Rps16a - x_Rps16a_Fun12 <= 0 \Steiner node out
y_Rps16b - x_Rps16b_Fun12 <= 0 \Steiner node out
y_Rps1a - x_Rps1a_Fun12 <= 0 \Steiner node out
y_Rps23a - x_Rps23a_Fun12 <= 0 \Steiner node out
y_Rps23b - x_Rps23b_Fun12 <= 0 \Steiner node out
y_Rps24a - x_Rps24a_Fun12 <= 0 \Steiner node out
y_Rps24b - x_Rps24b_Fun12 <= 0 \Steiner node out
y_Rps27a - x_Rps27a_Fun12 <= 0 \Steiner node out
y_Rps28a - x_Rps28a_Fun12 <= 0 \Steiner node out
y_Rps5 - x_Rps5_Fun12 <= 0 \Steiner node out
y_Rps6a - x_Rps6a_Fun12 <= 0 \Steiner node out
y_Rps6b - x_Rps6b_Fun12 <= 0 \Steiner node out
y_Rps7a - x_Rps7a_Fun12 <= 0 \Steiner node out
y_Rps8a - x_Rps8a_Fun12 <= 0 \Steiner node out
y_Rps8b - x_Rps8b_Fun12 <= 0 \Steiner node out
y_Rps9a - x_Rps9a_Fun12 <= 0 \Steiner node out
y_Npl3 - x_Npl3_Rpo21 <= 0 \Steiner node out
y_Rpo21 - x_Rpo21_Rtf1 <= 0 \Steiner node out
y_Rsp5 - x_Rsp5_Rpo21 - x_Rsp5_Rps31 - x_Rsp5_Ubi4 <= 0 \Steiner node out
y_Rps31 - x_Rps31_Rps11a - x_Rps31_Rps11b - x_Rps31_Rps13 - x_Rps31_Rps14a - x_Rps31_Rps16a - x_Rps31_Rps16b - x_Rps31_Rps1a - x_Rps31_Rps23a - x_Rps31_Rps23b - x_Rps31_Rps24a - x_Rps31_Rps24b - x_Rps31_Rps27a - x_Rps31_Rps28a - x_Rps31_Rps5 - x_Rps31_Rps6a - x_Rps31_Rps6b - x_Rps31_Rps7a - x_Rps31_Rps8a - x_Rps31_Rps8b - x_Rps31_Rps9a - x_Rps31_Ubi4 - x_Rps31_Vps27 - x_Rps31_Asc1 - x_Rps31_Dom34 - x_Rps31_Gcd11 - x_Rps31_Nip1 - x_Rps31_Prt1 - x_Rps31_Rpg1 - x_Rps31_Rpl10 - x_Rps31_Rpl17a - x_Rps31_Rpl17b - x_Rps31_Rpl19b - x_Rps31_Rpl25 - x_Rps31_Rpl28 - x_Rps31_Rpl2a - x_Rps31_Rpl2b - x_Rps31_Rpl3 - x_Rps31_Rpl30 - x_Rps31_Rpl32 - x_Rps31_Rpl34a - x_Rps31_Rpl35a - x_Rps31_Rpl35b - x_Rps31_Rpl37a - x_Rps31_Rpl38 - x_Rps31_Rpl5 - x_Rps31_Rpl7a - x_Rps31_Rps0a - x_Rps31_Rps0b - x_Rps31_Rps10a - x_Rps31_Rps12 - x_Rps31_Rps14b - x_Rps31_Rps15 - x_Rps31_Rps17a - x_Rps31_Rps17b - x_Rps31_Rps18a - x_Rps31_Rps18b - x_Rps31_Rps19a - x_Rps31_Rps1b - x_Rps31_Rps2 - x_Rps31_Rps20 - x_Rps31_Rps22a - x_Rps31_Rps25a - x_Rps31_Rps26a - x_Rps31_Rps26b - x_Rps31_Rps28b - x_Rps31_Rps29a - x_Rps31_Rps3 - x_Rps31_Rps30a - x_Rps31_Rps30b - x_Rps31_Rps9b - x_Rps31_Sui1 - x_Rps31_Tif11 - x_Rps31_Tif35 - x_Rps31_Tif5 - x_Rps31_Fun12 - x_Rps31_Ede1 <= 0 \Steiner node out
y_Ubi4 - x_Ubi4_Rps11a - x_Ubi4_Rps11b - x_Ubi4_Rps13 - x_Ubi4_Rps14a - x_Ubi4_Rps16a - x_Ubi4_Rps16b - x_Ubi4_Rps1a - x_Ubi4_Rps23a - x_Ubi4_Rps23b - x_Ubi4_Rps24a - x_Ubi4_Rps24b - x_Ubi4_Rps27a - x_Ubi4_Rps28a - x_Ubi4_Rps5 - x_Ubi4_Rps6b - x_Ubi4_Rps7a - x_Ubi4_Rps8a - x_Ubi4_Rps8b - x_Ubi4_Rps9a - x_Ubi4_Rps31 - x_Ubi4_Ent1 - x_Ubi4_Ent2 - x_Ubi4_Vps27 - x_Ubi4_Asc1 - x_Ubi4_Rpl10 - x_Ubi4_Rpl17a - x_Ubi4_Rpl17b - x_Ubi4_Rpl19b - x_Ubi4_Rpl25 - x_Ubi4_Rpl28 - x_Ubi4_Rpl2a - x_Ubi4_Rpl2b - x_Ubi4_Rpl3 - x_Ubi4_Rpl30 - x_Ubi4_Rpl32 - x_Ubi4_Rpl34a - x_Ubi4_Rpl35a - x_Ubi4_Rpl35b - x_Ubi4_Rpl38 - x_Ubi4_Rpl5 - x_Ubi4_Rps0a - x_Ubi4_Rps0b - x_Ubi4_Rps10a - x_Ubi4_Rps12 - x_Ubi4_Rps14b - x_Ubi4_Rps15 - x_Ubi4_Rps17a - x_Ubi4_Rps17b - x_Ubi4_Rps18b - x_Ubi4_Rps19a - x_Ubi4_Rps1b - x_Ubi4_Rps2 - x_Ubi4_Rps20 - x_Ubi4_Rps25a - x_Ubi4_Rps26a - x_Ubi4_Rps26b - x_Ubi4_Rps28b - x_Ubi4_Rps29a - x_Ubi4_Rps3 - x_Ubi4_Rps9b - x_Ubi4_Ede1 - x_Ubi4_Duf1 <= 0 \Steiner node out
y_Fps1 - x_Fps1_Hog1 - x_Fps1_Slt2 <= 0 \Steiner node out
y_Hog1 - x_Hog1_Sko1 <= 0 \Steiner node out
y_Slt2 - x_Slt2_Hog1 - x_Slt2_Cdc28 - x_Slt2_Sko1 <= 0 \Steiner node out
y_Rpb10 - x_Rpb10_Rpo21 - x_Rpb10_Cdc73 - x_Rpb10_Rpb2 - x_Rpb10_Rpb3 - x_Rpb10_Spt4 - x_Rpb10_Spt5 <= 0 \Steiner node out
y_Cdc73 - x_Cdc73_Rtf1 <= 0 \Steiner node out
y_Rpb2 - x_Rpb2_Rtf1 <= 0 \Steiner node out
y_Rpb3 - x_Rpb3_Rtf1 <= 0 \Steiner node out
y_Spt4 - x_Spt4_Rtf1 <= 0 \Steiner node out
y_Spt5 - x_Spt5_Rtf1 <= 0 \Steiner node out
y_Rpb5 - x_Rpb5_Rpo21 - x_Rpb5_Cdc73 - x_Rpb5_Rpb2 - x_Rpb5_Rpb3 - x_Rpb5_Spt4 - x_Rpb5_Spt5 - x_Rpb5_Leo1 - x_Rpb5_Spt6 <= 0 \Steiner node out
y_Leo1 - x_Leo1_Rtf1 <= 0 \Steiner node out
y_Spt6 - x_Spt6_Rtf1 <= 0 \Steiner node out
y_Rpb8 - x_Rpb8_Rpo21 - x_Rpb8_Rpb2 - x_Rpb8_Rpb3 - x_Rpb8_Spt4 - x_Rpb8_Spt5 <= 0 \Steiner node out
y_Rpc10 - x_Rpc10_Rpo21 - x_Rpc10_Cdc73 - x_Rpc10_Rpb2 - x_Rpc10_Rpb3 - x_Rpc10_Spt4 - x_Rpc10_Spt5 <= 0 \Steiner node out
y_Rpo26 - x_Rpo26_Rpo21 - x_Rpo26_Cdc73 - x_Rpo26_Rpb2 - x_Rpo26_Rpb3 - x_Rpo26_Spt4 - x_Rpo26_Spt5 - x_Rpo26_Paf1 <= 0 \Steiner node out
y_Paf1 - x_Paf1_Rtf1 <= 0 \Steiner node out
y_Rpo31 - x_Rpo31_Rpo21 - x_Rpo31_Rpb3 <= 0 \Steiner node out
y_Spt15 - x_Spt15_Rpo21 - x_Spt15_Rpb2 - x_Spt15_Rpb3 <= 0 \Steiner node out
y_Apl1 - x_Apl1_Apl3 <= 0 \Steiner node out
y_Apl3 - x_Apl3_Ent1 - x_Apl3_Ent2 - x_Apl3_Ede1 <= 0 \Steiner node out
y_Apl2 - x_Apl2_Apl3 <= 0 \Steiner node out
y_Ent1 - x_Ent1_Ubi4 - x_Ent1_Apl3 - x_Ent1_Ede1 <= 0 \Steiner node out
y_Ent2 - x_Ent2_Ubi4 - x_Ent2_Apl3 - x_Ent2_Ede1 <= 0 \Steiner node out
y_Las17 - x_Las17_Act1 <= 0 \Steiner node out
y_Act1 - x_Act1_Cyr1 <= 0 \Steiner node out
y_Sla2 - x_Sla2_Ent1 - x_Sla2_Act1 <= 0 \Steiner node out
y_Swa2 - x_Swa2_Ubi4 <= 0 \Steiner node out
y_Vps27 - x_Vps27_Rps31 - x_Vps27_Ubi4 - x_Vps27_Ede1 <= 0 \Steiner node out
y_Fus3 - x_Fus3_Gpa1 - x_Fus3_Kss1 - x_Fus3_Sko1 <= 0 \Steiner node out
y_Gpa1 - x_Gpa1_Cyr1 <= 0 \Steiner node out
y_Kss1 - x_Kss1_Hog1 - x_Kss1_Fus3 - x_Kss1_Sko1 <= 0 \Steiner node out
y_Cdc28 - x_Cdc28_Swe1 <= 0 \Steiner node out
y_Apm4 - x_Apm4_Apl3 <= 0 \Steiner node out
y_Abp1 - x_Abp1_Act1 - x_Abp1_Srv2 <= 0 \Steiner node out
y_Srv2 - x_Srv2_Cyr1 <= 0 \Steiner node out
y_Asc1 - x_Asc1_Fun12 <= 0 \Steiner node out
y_Dom34 - x_Dom34_Sup35 <= 0 \Steiner node out
y_Gcd11 - x_Gcd11_Fun12 <= 0 \Steiner node out
y_Nip1 - x_Nip1_Fun12 <= 0 \Steiner node out
y_Prt1 - x_Prt1_Fun12 <= 0 \Steiner node out
y_Rpg1 - x_Rpg1_Fun12 <= 0 \Steiner node out
y_Rpl10 - x_Rpl10_Fun12 <= 0 \Steiner node out
y_Rpl17a - x_Rpl17a_Fun12 <= 0 \Steiner node out
y_Rpl17b - x_Rpl17b_Fun12 <= 0 \Steiner node out
y_Rpl19b - x_Rpl19b_Fun12 <= 0 \Steiner node out
y_Rpl25 - x_Rpl25_Fun12 <= 0 \Steiner node out
y_Rpl28 - x_Rpl28_Fun12 <= 0 \Steiner node out
y_Rpl2a - x_Rpl2a_Fun12 <= 0 \Steiner node out
y_Rpl2b - x_Rpl2b_Fun12 <= 0 \Steiner node out
y_Rpl3 - x_Rpl3_Fun12 <= 0 \Steiner node out
y_Rpl30 - x_Rpl30_Fun12 <= 0 \Steiner node out
y_Rpl32 - x_Rpl32_Fun12 <= 0 \Steiner node out
y_Rpl34a - x_Rpl34a_Fun12 <= 0 \Steiner node out
y_Rpl35a - x_Rpl35a_Fun12 <= 0 \Steiner node out
y_Rpl35b - x_Rpl35b_Fun12 <= 0 \Steiner node out
y_Rpl37a - x_Rpl37a_Fun12 <= 0 \Steiner node out
y_Rpl38 - x_Rpl38_Fun12 <= 0 \Steiner node out
y_Rpl5 - x_Rpl5_Fun12 <= 0 \Steiner node out
y_Rpl7a - x_Rpl7a_Fun12 <= 0 \Steiner node out
y_Rps0a - x_Rps0a_Fun12 <= 0 \Steiner node out
y_Rps0b - x_Rps0b_Fun12 <= 0 \Steiner node out
y_Rps10a - x_Rps10a_Fun12 <= 0 \Steiner node out
y_Rps12 - x_Rps12_Fun12 <= 0 \Steiner node out
y_Rps14b - x_Rps14b_Fun12 <= 0 \Steiner node out
y_Rps15 - x_Rps15_Fun12 <= 0 \Steiner node out
y_Rps17a - x_Rps17a_Fun12 <= 0 \Steiner node out
y_Rps17b - x_Rps17b_Fun12 <= 0 \Steiner node out
y_Rps18a - x_Rps18a_Fun12 <= 0 \Steiner node out
y_Rps18b - x_Rps18b_Fun12 <= 0 \Steiner node out
y_Rps19a - x_Rps19a_Fun12 <= 0 \Steiner node out
y_Rps1b - x_Rps1b_Fun12 <= 0 \Steiner node out
y_Rps2 - x_Rps2_Fun12 <= 0 \Steiner node out
y_Rps20 - x_Rps20_Fun12 <= 0 \Steiner node out
y_Rps22a - x_Rps22a_Fun12 <= 0 \Steiner node out
y_Rps25a - x_Rps25a_Fun12 <= 0 \Steiner node out
y_Rps26a - x_Rps26a_Fun12 <= 0 \Steiner node out
y_Rps26b - x_Rps26b_Fun12 <= 0 \Steiner node out
y_Rps28b - x_Rps28b_Fun12 <= 0 \Steiner node out
y_Rps29a - x_Rps29a_Fun12 <= 0 \Steiner node out
y_Rps3 - x_Rps3_Fun12 <= 0 \Steiner node out
y_Rps30a - x_Rps30a_Fun12 <= 0 \Steiner node out
y_Rps30b - x_Rps30b_Fun12 <= 0 \Steiner node out
y_Rps9b - x_Rps9b_Fun12 <= 0 \Steiner node out
y_Sui1 - x_Sui1_Fun12 <= 0 \Steiner node out
y_Tif11 - x_Tif11_Fun12 <= 0 \Steiner node out
y_Tif35 - x_Tif35_Fun12 <= 0 \Steiner node out
y_Tif5 - x_Tif5_Fun12 <= 0 \Steiner node out
y_Rvs167 - x_Rvs167_Act1 <= 0 \Steiner node out
y_Vrp1 - x_Vrp1_Act1 <= 0 \Steiner node out
y_Pbs2 - x_Pbs2_Hog1 - x_Pbs2_Fus3 - x_Pbs2_Kss1 - x_Pbs2_Mih1 <= 0 \Steiner node out
y_Mih1 - x_Mih1_Swe1 <= 0 \Steiner node out
y_Bud3 - x_Bud3_Clb2 <= 0 \Steiner node out
y_Clb2 - x_Clb2_Swe1 <= 0 \Steiner node out
y_Cdc55 - x_Cdc55_Swe1 - x_Cdc55_Igo2 <= 0 \Steiner node out
y_Yak1 - x_Yak1_Rim15 - x_Yak1_Tpk1 <= 0 \Steiner node out
y_Rim15 - x_Rim15_Msn2 - x_Rim15_Sch9 <= 0 \Steiner node out
y_Tpk1 - x_Tpk1_Sch9 <= 0 \Steiner node out
y_Mtr4 - x_Air2_Mtr4 <= 0 \Steiner node in
y_Rps11a - x_Mtr4_Rps11a - x_Rps31_Rps11a - x_Ubi4_Rps11a <= 0 \Steiner node in
y_Rps11b - x_Mtr4_Rps11b - x_Rps31_Rps11b - x_Ubi4_Rps11b <= 0 \Steiner node in
y_Rps13 - x_Mtr4_Rps13 - x_Rps31_Rps13 - x_Ubi4_Rps13 <= 0 \Steiner node in
y_Rps14a - x_Mtr4_Rps14a - x_Rps31_Rps14a - x_Ubi4_Rps14a <= 0 \Steiner node in
y_Rps16a - x_Mtr4_Rps16a - x_Rps31_Rps16a - x_Ubi4_Rps16a <= 0 \Steiner node in
y_Rps16b - x_Mtr4_Rps16b - x_Rps31_Rps16b - x_Ubi4_Rps16b <= 0 \Steiner node in
y_Rps1a - x_Mtr4_Rps1a - x_Rps31_Rps1a - x_Ubi4_Rps1a <= 0 \Steiner node in
y_Rps23a - x_Mtr4_Rps23a - x_Rps31_Rps23a - x_Ubi4_Rps23a <= 0 \Steiner node in
y_Rps23b - x_Mtr4_Rps23b - x_Rps31_Rps23b - x_Ubi4_Rps23b <= 0 \Steiner node in
y_Rps24a - x_Mtr4_Rps24a - x_Rps31_Rps24a - x_Ubi4_Rps24a <= 0 \Steiner node in
y_Rps24b - x_Mtr4_Rps24b - x_Rps31_Rps24b - x_Ubi4_Rps24b <= 0 \Steiner node in
y_Rps27a - x_Mtr4_Rps27a - x_Rps31_Rps27a - x_Ubi4_Rps27a <= 0 \Steiner node in
y_Rps28a - x_Mtr4_Rps28a - x_Rps31_Rps28a - x_Ubi4_Rps28a <= 0 \Steiner node in
y_Rps5 - x_Mtr4_Rps5 - x_Rps31_Rps5 - x_Ubi4_Rps5 <= 0 \Steiner node in
y_Rps6a - x_Mtr4_Rps6a - x_Rps31_Rps6a <= 0 \Steiner node in
y_Rps6b - x_Mtr4_Rps6b - x_Rps31_Rps6b - x_Ubi4_Rps6b <= 0 \Steiner node in
y_Rps7a - x_Mtr4_Rps7a - x_Rps31_Rps7a - x_Ubi4_Rps7a <= 0 \Steiner node in
y_Rps8a - x_Mtr4_Rps8a - x_Rps31_Rps8a - x_Ubi4_Rps8a <= 0 \Steiner node in
y_Rps8b - x_Mtr4_Rps8b - x_Rps31_Rps8b - x_Ubi4_Rps8b <= 0 \Steiner node in
y_Rps9a - x_Mtr4_Rps9a - x_Rps31_Rps9a - x_Ubi4_Rps9a <= 0 \Steiner node in
y_Npl3 - x_Air2_Npl3 <= 0 \Steiner node in
y_Rpo21 - x_Npl3_Rpo21 - x_Rsp5_Rpo21 - x_Rpb10_Rpo21 - x_Rpb5_Rpo21 - x_Rpb8_Rpo21 - x_Rpc10_Rpo21 - x_Rpo26_Rpo21 - x_Rpo31_Rpo21 - x_Spt15_Rpo21 <= 0 \Steiner node in
y_Rsp5 - x_Art5_Rsp5 - x_Ecm21_Rsp5 - x_Ldb19_Rsp5 <= 0 \Steiner node in
y_Rps31 - x_Rsp5_Rps31 - x_Ubi4_Rps31 - x_Vps27_Rps31 - x_Sla1_Rps31 <= 0 \Steiner node in
y_Ubi4 - x_Rsp5_Ubi4 - x_Rps31_Ubi4 - x_Ent1_Ubi4 - x_Ent2_Ubi4 - x_Swa2_Ubi4 - x_Vps27_Ubi4 - x_Sla1_Ubi4 <= 0 \Steiner node in
y_Fps1 - x_Ask10_Fps1 <= 0 \Steiner node in
y_Hog1 - x_Fps1_Hog1 - x_Slt2_Hog1 - x_Kss1_Hog1 - x_Pbs2_Hog1 <= 0 \Steiner node in
y_Slt2 - x_Fps1_Slt2 - x_Msg5_Slt2 <= 0 \Steiner node in
y_Rpb10 - x_Bdp1_Rpb10 <= 0 \Steiner node in
y_Cdc73 - x_Rpb10_Cdc73 - x_Rpb5_Cdc73 - x_Rpc10_Cdc73 - x_Rpo26_Cdc73 <= 0 \Steiner node in
y_Rpb2 - x_Rpb10_Rpb2 - x_Rpb5_Rpb2 - x_Rpb8_Rpb2 - x_Rpc10_Rpb2 - x_Rpo26_Rpb2 - x_Spt15_Rpb2 <= 0 \Steiner node in
y_Rpb3 - x_Rpb10_Rpb3 - x_Rpb5_Rpb3 - x_Rpb8_Rpb3 - x_Rpc10_Rpb3 - x_Rpo26_Rpb3 - x_Rpo31_Rpb3 - x_Spt15_Rpb3 <= 0 \Steiner node in
y_Spt4 - x_Rpb10_Spt4 - x_Rpb5_Spt4 - x_Rpb8_Spt4 - x_Rpc10_Spt4 - x_Rpo26_Spt4 <= 0 \Steiner node in
y_Spt5 - x_Rpb10_Spt5 - x_Rpb5_Spt5 - x_Rpb8_Spt5 - x_Rpc10_Spt5 - x_Rpo26_Spt5 <= 0 \Steiner node in
y_Rpb5 - x_Bdp1_Rpb5 <= 0 \Steiner node in
y_Leo1 - x_Rpb5_Leo1 <= 0 \Steiner node in
y_Spt6 - x_Rpb5_Spt6 <= 0 \Steiner node in
y_Rpb8 - x_Bdp1_Rpb8 <= 0 \Steiner node in
y_Rpc10 - x_Bdp1_Rpc10 <= 0 \Steiner node in
y_Rpo26 - x_Bdp1_Rpo26 <= 0 \Steiner node in
y_Paf1 - x_Rpo26_Paf1 <= 0 \Steiner node in
y_Rpo31 - x_Bdp1_Rpo31 <= 0 \Steiner node in
y_Spt15 - x_Bdp1_Spt15 <= 0 \Steiner node in
y_Apl1 - x_Chc1_Apl1 <= 0 \Steiner node in
y_Apl3 - x_Apl1_Apl3 - x_Apl2_Apl3 - x_Ent1_Apl3 - x_Ent2_Apl3 - x_Sec16_Apl3 - x_Apm4_Apl3 <= 0 \Steiner node in
y_Apl2 - x_Chc1_Apl2 <= 0 \Steiner node in
y_Ent1 - x_Ubi4_Ent1 - x_Chc1_Ent1 - x_Apl3_Ent1 - x_Sla2_Ent1 <= 0 \Steiner node in
y_Ent2 - x_Ubi4_Ent2 - x_Chc1_Ent2 - x_Apl3_Ent2 <= 0 \Steiner node in
y_Las17 - x_Chc1_Las17 - x_Sla1_Las17 <= 0 \Steiner node in
y_Act1 - x_Las17_Act1 - x_Sla2_Act1 - x_Abp1_Act1 - x_Rvs167_Act1 - x_Vrp1_Act1 <= 0 \Steiner node in
y_Sla2 - x_Chc1_Sla2 - x_Sla1_Sla2 <= 0 \Steiner node in
y_Swa2 - x_Chc1_Swa2 <= 0 \Steiner node in
y_Vps27 - x_Rps31_Vps27 - x_Ubi4_Vps27 - x_Chc1_Vps27 <= 0 \Steiner node in
y_Fus3 - x_Msg5_Fus3 - x_Kss1_Fus3 - x_Pbs2_Fus3 <= 0 \Steiner node in
y_Gpa1 - x_Fus3_Gpa1 <= 0 \Steiner node in
y_Kss1 - x_Msg5_Kss1 - x_Fus3_Kss1 - x_Pbs2_Kss1 <= 0 \Steiner node in
y_Cdc28 - x_Slt2_Cdc28 <= 0 \Steiner node in
y_Apm4 - x_Sec16_Apm4 <= 0 \Steiner node in
y_Abp1 - x_Sla1_Abp1 <= 0 \Steiner node in
y_Srv2 - x_Abp1_Srv2 <= 0 \Steiner node in
y_Asc1 - x_Rps31_Asc1 - x_Ubi4_Asc1 <= 0 \Steiner node in
y_Dom34 - x_Rps31_Dom34 <= 0 \Steiner node in
y_Gcd11 - x_Rps31_Gcd11 <= 0 \Steiner node in
y_Nip1 - x_Rps31_Nip1 <= 0 \Steiner node in
y_Prt1 - x_Rps31_Prt1 <= 0 \Steiner node in
y_Rpg1 - x_Rps31_Rpg1 <= 0 \Steiner node in
y_Rpl10 - x_Rps31_Rpl10 - x_Ubi4_Rpl10 <= 0 \Steiner node in
y_Rpl17a - x_Rps31_Rpl17a - x_Ubi4_Rpl17a <= 0 \Steiner node in
y_Rpl17b - x_Rps31_Rpl17b - x_Ubi4_Rpl17b <= 0 \Steiner node in
y_Rpl19b - x_Rps31_Rpl19b - x_Ubi4_Rpl19b <= 0 \Steiner node in
y_Rpl25 - x_Rps31_Rpl25 - x_Ubi4_Rpl25 <= 0 \Steiner node in
y_Rpl28 - x_Rps31_Rpl28 - x_Ubi4_Rpl28 <= 0 \Steiner node in
y_Rpl2a - x_Rps31_Rpl2a - x_Ubi4_Rpl2a <= 0 \Steiner node in
y_Rpl2b - x_Rps31_Rpl2b - x_Ubi4_Rpl2b <= 0 \Steiner node in
y_Rpl3 - x_Rps31_Rpl3 - x_Ubi4_Rpl3 <= 0 \Steiner node in
y_Rpl30 - x_Rps31_Rpl30 - x_Ubi4_Rpl30 <= 0 \Steiner node in
y_Rpl32 - x_Rps31_Rpl32 - x_Ubi4_Rpl32 <= 0 \Steiner node in
y_Rpl34a - x_Rps31_Rpl34a - x_Ubi4_Rpl34a <= 0 \Steiner node in
y_Rpl35a - x_Rps31_Rpl35a - x_Ubi4_Rpl35a <= 0 \Steiner node in
y_Rpl35b - x_Rps31_Rpl35b - x_Ubi4_Rpl35b <= 0 \Steiner node in
y_Rpl37a - x_Rps31_Rpl37a <= 0 \Steiner node in
y_Rpl38 - x_Rps31_Rpl38 - x_Ubi4_Rpl38 <= 0 \Steiner node in
y_Rpl5 - x_Rps31_Rpl5 - x_Ubi4_Rpl5 <= 0 \Steiner node in
y_Rpl7a - x_Rps31_Rpl7a <= 0 \Steiner node in
y_Rps0a - x_Rps31_Rps0a - x_Ubi4_Rps0a <= 0 \Steiner node in
y_Rps0b - x_Rps31_Rps0b - x_Ubi4_Rps0b <= 0 \Steiner node in
y_Rps10a - x_Rps31_Rps10a - x_Ubi4_Rps10a <= 0 \Steiner node in
y_Rps12 - x_Rps31_Rps12 - x_Ubi4_Rps12 <= 0 \Steiner node in
y_Rps14b - x_Rps31_Rps14b - x_Ubi4_Rps14b <= 0 \Steiner node in
y_Rps15 - x_Rps31_Rps15 - x_Ubi4_Rps15 <= 0 \Steiner node in
y_Rps17a - x_Rps31_Rps17a - x_Ubi4_Rps17a <= 0 \Steiner node in
y_Rps17b - x_Rps31_Rps17b - x_Ubi4_Rps17b <= 0 \Steiner node in
y_Rps18a - x_Rps31_Rps18a <= 0 \Steiner node in
y_Rps18b - x_Rps31_Rps18b - x_Ubi4_Rps18b <= 0 \Steiner node in
y_Rps19a - x_Rps31_Rps19a - x_Ubi4_Rps19a <= 0 \Steiner node in
y_Rps1b - x_Rps31_Rps1b - x_Ubi4_Rps1b <= 0 \Steiner node in
y_Rps2 - x_Rps31_Rps2 - x_Ubi4_Rps2 <= 0 \Steiner node in
y_Rps20 - x_Rps31_Rps20 - x_Ubi4_Rps20 <= 0 \Steiner node in
y_Rps22a - x_Rps31_Rps22a <= 0 \Steiner node in
y_Rps25a - x_Rps31_Rps25a - x_Ubi4_Rps25a <= 0 \Steiner node in
y_Rps26a - x_Rps31_Rps26a - x_Ubi4_Rps26a <= 0 \Steiner node in
y_Rps26b - x_Rps31_Rps26b - x_Ubi4_Rps26b <= 0 \Steiner node in
y_Rps28b - x_Rps31_Rps28b - x_Ubi4_Rps28b <= 0 \Steiner node in
y_Rps29a - x_Rps31_Rps29a - x_Ubi4_Rps29a <= 0 \Steiner node in
y_Rps3 - x_Rps31_Rps3 - x_Ubi4_Rps3 <= 0 \Steiner node in
y_Rps30a - x_Rps31_Rps30a <= 0 \Steiner node in
y_Rps30b - x_Rps31_Rps30b <= 0 \Steiner node in
y_Rps9b - x_Rps31_Rps9b - x_Ubi4_Rps9b <= 0 \Steiner node in
y_Sui1 - x_Rps31_Sui1 <= 0 \Steiner node in
y_Tif11 - x_Rps31_Tif11 <= 0 \Steiner node in
y_Tif35 - x_Rps31_Tif35 <= 0 \Steiner node in
y_Tif5 - x_Rps31_Tif5 <= 0 \Steiner node in
y_Rvs167 - x_Sla1_Rvs167 <= 0 \Steiner node in
y_Vrp1 - x_Sla1_Vrp1 <= 0 \Steiner node in
y_Pbs2 - x_Ssk2_Pbs2 <= 0 \Steiner node in
y_Mih1 - x_Pbs2_Mih1 <= 0 \Steiner node in
y_Bud3 - x_Bud4_Bud3 <= 0 \Steiner node in
y_Clb2 - x_Bud3_Clb2 <= 0 \Steiner node in
y_Cdc55 - x_Igo1_Cdc55 <= 0 \Steiner node in
y_Yak1 - x_Ypl247c_Yak1 <= 0 \Steiner node in
y_Rim15 - x_Yak1_Rim15 <= 0 \Steiner node in
y_Tpk1 - x_Yak1_Tpk1 <= 0 \Steiner node in
y_Air2 - x_Air2_Mtr4 >= 0 \select node if connected (outgoing)
y_Air2 - x_Air2_Npl3 >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps11a >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps11b >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps13 >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps14a >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps16a >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps16b >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps1a >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps23a >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps23b >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps24a >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps24b >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps27a >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps28a >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps5 >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps6a >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps6b >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps7a >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps8a >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps8b >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Mtr4_Rps9a >= 0 \select node if connected (outgoing)
y_Rps11a - x_Rps11a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps11b - x_Rps11b_Fun12 >= 0 \select node if connected (outgoing)
y_Rps13 - x_Rps13_Fun12 >= 0 \select node if connected (outgoing)
y_Rps14a - x_Rps14a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps16a - x_Rps16a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps16b - x_Rps16b_Fun12 >= 0 \select node if connected (outgoing)
y_Rps1a - x_Rps1a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps23a - x_Rps23a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps23b - x_Rps23b_Fun12 >= 0 \select node if connected (outgoing)
y_Rps24a - x_Rps24a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps24b - x_Rps24b_Fun12 >= 0 \select node if connected (outgoing)
y_Rps27a - x_Rps27a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps28a - x_Rps28a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps5 - x_Rps5_Fun12 >= 0 \select node if connected (outgoing)
y_Rps6a - x_Rps6a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps6b - x_Rps6b_Fun12 >= 0 \select node if connected (outgoing)
y_Rps7a - x_Rps7a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps8a - x_Rps8a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps8b - x_Rps8b_Fun12 >= 0 \select node if connected (outgoing)
y_Rps9a - x_Rps9a_Fun12 >= 0 \select node if connected (outgoing)
y_Npl3 - x_Npl3_Rpo21 >= 0 \select node if connected (outgoing)
y_Rpo21 - x_Rpo21_Rtf1 >= 0 \select node if connected (outgoing)
y_Art5 - x_Art5_Rsp5 >= 0 \select node if connected (outgoing)
y_Rsp5 - x_Rsp5_Rpo21 >= 0 \select node if connected (outgoing)
y_Rsp5 - x_Rsp5_Rps31 >= 0 \select node if connected (outgoing)
y_Rsp5 - x_Rsp5_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps11a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps11b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps13 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps14a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps16a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps16b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps1a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps23a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps23b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps24a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps24b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps27a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps28a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps5 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps6a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps6b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps7a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps8a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps8b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps9a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Ubi4 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Vps27 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Asc1 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Dom34 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Gcd11 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Nip1 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Prt1 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpg1 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl10 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl17a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl17b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl19b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl25 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl28 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl2a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl2b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl3 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl30 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl32 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl34a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl35a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl35b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl37a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl38 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl5 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rpl7a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps0a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps0b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps10a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps12 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps14b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps15 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps17a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps17b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps18a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps18b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps19a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps1b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps2 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps20 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps22a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps25a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps26a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps26b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps28b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps29a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps3 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps30a >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps30b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Rps9b >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Sui1 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Tif11 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Tif35 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Tif5 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Fun12 >= 0 \select node if connected (outgoing)
y_Rps31 - x_Rps31_Ede1 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps11a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps11b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps13 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps14a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps16a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps16b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps1a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps23a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps23b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps24a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps24b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps27a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps28a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps5 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps6b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps7a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps8a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps8b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps9a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps31 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Ent1 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Ent2 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Vps27 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Asc1 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl10 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl17a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl17b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl19b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl25 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl28 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl2a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl2b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl3 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl30 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl32 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl34a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl35a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl35b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl38 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rpl5 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps0a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps0b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps10a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps12 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps14b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps15 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps17a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps17b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps18b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps19a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps1b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps2 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps20 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps25a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps26a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps26b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps28b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps29a >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps3 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Rps9b >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Ede1 >= 0 \select node if connected (outgoing)
y_Ubi4 - x_Ubi4_Duf1 >= 0 \select node if connected (outgoing)
y_Ask10 - x_Ask10_Fps1 >= 0 \select node if connected (outgoing)
y_Fps1 - x_Fps1_Hog1 >= 0 \select node if connected (outgoing)
y_Fps1 - x_Fps1_Slt2 >= 0 \select node if connected (outgoing)
y_Hog1 - x_Hog1_Sko1 >= 0 \select node if connected (outgoing)
y_Slt2 - x_Slt2_Hog1 >= 0 \select node if connected (outgoing)
y_Slt2 - x_Slt2_Cdc28 >= 0 \select node if connected (outgoing)
y_Slt2 - x_Slt2_Sko1 >= 0 \select node if connected (outgoing)
y_Bdp1 - x_Bdp1_Rpb10 >= 0 \select node if connected (outgoing)
y_Bdp1 - x_Bdp1_Rpb5 >= 0 \select node if connected (outgoing)
y_Bdp1 - x_Bdp1_Rpb8 >= 0 \select node if connected (outgoing)
y_Bdp1 - x_Bdp1_Rpc10 >= 0 \select node if connected (outgoing)
y_Bdp1 - x_Bdp1_Rpo26 >= 0 \select node if connected (outgoing)
y_Bdp1 - x_Bdp1_Rpo31 >= 0 \select node if connected (outgoing)
y_Bdp1 - x_Bdp1_Spt15 >= 0 \select node if connected (outgoing)
y_Rpb10 - x_Rpb10_Rpo21 >= 0 \select node if connected (outgoing)
y_Rpb10 - x_Rpb10_Cdc73 >= 0 \select node if connected (outgoing)
y_Rpb10 - x_Rpb10_Rpb2 >= 0 \select node if connected (outgoing)
y_Rpb10 - x_Rpb10_Rpb3 >= 0 \select node if connected (outgoing)
y_Rpb10 - x_Rpb10_Spt4 >= 0 \select node if connected (outgoing)
y_Rpb10 - x_Rpb10_Spt5 >= 0 \select node if connected (outgoing)
y_Cdc73 - x_Cdc73_Rtf1 >= 0 \select node if connected (outgoing)
y_Rpb2 - x_Rpb2_Rtf1 >= 0 \select node if connected (outgoing)
y_Rpb3 - x_Rpb3_Rtf1 >= 0 \select node if connected (outgoing)
y_Spt4 - x_Spt4_Rtf1 >= 0 \select node if connected (outgoing)
y_Spt5 - x_Spt5_Rtf1 >= 0 \select node if connected (outgoing)
y_Rpb5 - x_Rpb5_Rpo21 >= 0 \select node if connected (outgoing)
y_Rpb5 - x_Rpb5_Cdc73 >= 0 \select node if connected (outgoing)
y_Rpb5 - x_Rpb5_Rpb2 >= 0 \select node if connected (outgoing)
y_Rpb5 - x_Rpb5_Rpb3 >= 0 \select node if connected (outgoing)
y_Rpb5 - x_Rpb5_Spt4 >= 0 \select node if connected (outgoing)
y_Rpb5 - x_Rpb5_Spt5 >= 0 \select node if connected (outgoing)
y_Rpb5 - x_Rpb5_Leo1 >= 0 \select node if connected (outgoing)
y_Rpb5 - x_Rpb5_Spt6 >= 0 \select node if connected (outgoing)
y_Leo1 - x_Leo1_Rtf1 >= 0 \select node if connected (outgoing)
y_Spt6 - x_Spt6_Rtf1 >= 0 \select node if connected (outgoing)
y_Rpb8 - x_Rpb8_Rpo21 >= 0 \select node if connected (outgoing)
y_Rpb8 - x_Rpb8_Rpb2 >= 0 \select node if connected (outgoing)
y_Rpb8 - x_Rpb8_Rpb3 >= 0 \select node if connected (outgoing)
y_Rpb8 - x_Rpb8_Spt4 >= 0 \select node if connected (outgoing)
y_Rpb8 - x_Rpb8_Spt5 >= 0 \select node if connected (outgoing)
y_Rpc10 - x_Rpc10_Rpo21 >= 0 \select node if connected (outgoing)
y_Rpc10 - x_Rpc10_Cdc73 >= 0 \select node if connected (outgoing)
y_Rpc10 - x_Rpc10_Rpb2 >= 0 \select node if connected (outgoing)
y_Rpc10 - x_Rpc10_Rpb3 >= 0 \select node if connected (outgoing)
y_Rpc10 - x_Rpc10_Spt4 >= 0 \select node if connected (outgoing)
y_Rpc10 - x_Rpc10_Spt5 >= 0 \select node if connected (outgoing)
y_Rpo26 - x_Rpo26_Rpo21 >= 0 \select node if connected (outgoing)
y_Rpo26 - x_Rpo26_Cdc73 >= 0 \select node if connected (outgoing)
y_Rpo26 - x_Rpo26_Rpb2 >= 0 \select node if connected (outgoing)
y_Rpo26 - x_Rpo26_Rpb3 >= 0 \select node if connected (outgoing)
y_Rpo26 - x_Rpo26_Spt4 >= 0 \select node if connected (outgoing)
y_Rpo26 - x_Rpo26_Spt5 >= 0 \select node if connected (outgoing)
y_Rpo26 - x_Rpo26_Paf1 >= 0 \select node if connected (outgoing)
y_Paf1 - x_Paf1_Rtf1 >= 0 \select node if connected (outgoing)
y_Rpo31 - x_Rpo31_Rpo21 >= 0 \select node if connected (outgoing)
y_Rpo31 - x_Rpo31_Rpb3 >= 0 \select node if connected (outgoing)
y_Spt15 - x_Spt15_Rpo21 >= 0 \select node if connected (outgoing)
y_Spt15 - x_Spt15_Rpb2 >= 0 \select node if connected (outgoing)
y_Spt15 - x_Spt15_Rpb3 >= 0 \select node if connected (outgoing)
y_Chc1 - x_Chc1_Apl1 >= 0 \select node if connected (outgoing)
y_Chc1 - x_Chc1_Apl2 >= 0 \select node if connected (outgoing)
y_Chc1 - x_Chc1_Ent1 >= 0 \select node if connected (outgoing)
y_Chc1 - x_Chc1_Ent2 >= 0 \select node if connected (outgoing)
y_Chc1 - x_Chc1_Las17 >= 0 \select node if connected (outgoing)
y_Chc1 - x_Chc1_Sla2 >= 0 \select node if connected (outgoing)
y_Chc1 - x_Chc1_Swa2 >= 0 \select node if connected (outgoing)
y_Chc1 - x_Chc1_Vps27 >= 0 \select node if connected (outgoing)
y_Apl1 - x_Apl1_Apl3 >= 0 \select node if connected (outgoing)
y_Apl3 - x_Apl3_Ent1 >= 0 \select node if connected (outgoing)
y_Apl3 - x_Apl3_Ent2 >= 0 \select node if connected (outgoing)
y_Apl3 - x_Apl3_Ede1 >= 0 \select node if connected (outgoing)
y_Apl2 - x_Apl2_Apl3 >= 0 \select node if connected (outgoing)
y_Ent1 - x_Ent1_Ubi4 >= 0 \select node if connected (outgoing)
y_Ent1 - x_Ent1_Apl3 >= 0 \select node if connected (outgoing)
y_Ent1 - x_Ent1_Ede1 >= 0 \select node if connected (outgoing)
y_Ent2 - x_Ent2_Ubi4 >= 0 \select node if connected (outgoing)
y_Ent2 - x_Ent2_Apl3 >= 0 \select node if connected (outgoing)
y_Ent2 - x_Ent2_Ede1 >= 0 \select node if connected (outgoing)
y_Las17 - x_Las17_Act1 >= 0 \select node if connected (outgoing)
y_Act1 - x_Act1_Cyr1 >= 0 \select node if connected (outgoing)
y_Sla2 - x_Sla2_Ent1 >= 0 \select node if connected (outgoing)
y_Sla2 - x_Sla2_Act1 >= 0 \select node if connected (outgoing)
y_Swa2 - x_Swa2_Ubi4 >= 0 \select node if connected (outgoing)
y_Vps27 - x_Vps27_Rps31 >= 0 \select node if connected (outgoing)
y_Vps27 - x_Vps27_Ubi4 >= 0 \select node if connected (outgoing)
y_Vps27 - x_Vps27_Ede1 >= 0 \select node if connected (outgoing)
y_Ecm21 - x_Ecm21_Rsp5 >= 0 \select node if connected (outgoing)
y_Ldb19 - x_Ldb19_Rsp5 >= 0 \select node if connected (outgoing)
y_Msg5 - x_Msg5_Slt2 >= 0 \select node if connected (outgoing)
y_Msg5 - x_Msg5_Fus3 >= 0 \select node if connected (outgoing)
y_Msg5 - x_Msg5_Kss1 >= 0 \select node if connected (outgoing)
y_Fus3 - x_Fus3_Gpa1 >= 0 \select node if connected (outgoing)
y_Fus3 - x_Fus3_Kss1 >= 0 \select node if connected (outgoing)
y_Fus3 - x_Fus3_Sko1 >= 0 \select node if connected (outgoing)
y_Gpa1 - x_Gpa1_Cyr1 >= 0 \select node if connected (outgoing)
y_Kss1 - x_Kss1_Hog1 >= 0 \select node if connected (outgoing)
y_Kss1 - x_Kss1_Fus3 >= 0 \select node if connected (outgoing)
y_Kss1 - x_Kss1_Sko1 >= 0 \select node if connected (outgoing)
y_Cdc28 - x_Cdc28_Swe1 >= 0 \select node if connected (outgoing)
y_Sec16 - x_Sec16_Apl3 >= 0 \select node if connected (outgoing)
y_Sec16 - x_Sec16_Apm4 >= 0 \select node if connected (outgoing)
y_Apm4 - x_Apm4_Apl3 >= 0 \select node if connected (outgoing)
y_Sla1 - x_Sla1_Rps31 >= 0 \select node if connected (outgoing)
y_Sla1 - x_Sla1_Ubi4 >= 0 \select node if connected (outgoing)
y_Sla1 - x_Sla1_Las17 >= 0 \select node if connected (outgoing)
y_Sla1 - x_Sla1_Sla2 >= 0 \select node if connected (outgoing)
y_Sla1 - x_Sla1_Abp1 >= 0 \select node if connected (outgoing)
y_Sla1 - x_Sla1_Rvs167 >= 0 \select node if connected (outgoing)
y_Sla1 - x_Sla1_Vrp1 >= 0 \select node if connected (outgoing)
y_Abp1 - x_Abp1_Act1 >= 0 \select node if connected (outgoing)
y_Abp1 - x_Abp1_Srv2 >= 0 \select node if connected (outgoing)
y_Srv2 - x_Srv2_Cyr1 >= 0 \select node if connected (outgoing)
y_Asc1 - x_Asc1_Fun12 >= 0 \select node if connected (outgoing)
y_Dom34 - x_Dom34_Sup35 >= 0 \select node if connected (outgoing)
y_Gcd11 - x_Gcd11_Fun12 >= 0 \select node if connected (outgoing)
y_Nip1 - x_Nip1_Fun12 >= 0 \select node if connected (outgoing)
y_Prt1 - x_Prt1_Fun12 >= 0 \select node if connected (outgoing)
y_Rpg1 - x_Rpg1_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl10 - x_Rpl10_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl17a - x_Rpl17a_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl17b - x_Rpl17b_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl19b - x_Rpl19b_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl25 - x_Rpl25_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl28 - x_Rpl28_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl2a - x_Rpl2a_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl2b - x_Rpl2b_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl3 - x_Rpl3_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl30 - x_Rpl30_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl32 - x_Rpl32_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl34a - x_Rpl34a_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl35a - x_Rpl35a_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl35b - x_Rpl35b_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl37a - x_Rpl37a_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl38 - x_Rpl38_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl5 - x_Rpl5_Fun12 >= 0 \select node if connected (outgoing)
y_Rpl7a - x_Rpl7a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps0a - x_Rps0a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps0b - x_Rps0b_Fun12 >= 0 \select node if connected (outgoing)
y_Rps10a - x_Rps10a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps12 - x_Rps12_Fun12 >= 0 \select node if connected (outgoing)
y_Rps14b - x_Rps14b_Fun12 >= 0 \select node if connected (outgoing)
y_Rps15 - x_Rps15_Fun12 >= 0 \select node if connected (outgoing)
y_Rps17a - x_Rps17a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps17b - x_Rps17b_Fun12 >= 0 \select node if connected (outgoing)
y_Rps18a - x_Rps18a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps18b - x_Rps18b_Fun12 >= 0 \select node if connected (outgoing)
y_Rps19a - x_Rps19a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps1b - x_Rps1b_Fun12 >= 0 \select node if connected (outgoing)
y_Rps2 - x_Rps2_Fun12 >= 0 \select node if connected (outgoing)
y_Rps20 - x_Rps20_Fun12 >= 0 \select node if connected (outgoing)
y_Rps22a - x_Rps22a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps25a - x_Rps25a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps26a - x_Rps26a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps26b - x_Rps26b_Fun12 >= 0 \select node if connected (outgoing)
y_Rps28b - x_Rps28b_Fun12 >= 0 \select node if connected (outgoing)
y_Rps29a - x_Rps29a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps3 - x_Rps3_Fun12 >= 0 \select node if connected (outgoing)
y_Rps30a - x_Rps30a_Fun12 >= 0 \select node if connected (outgoing)
y_Rps30b - x_Rps30b_Fun12 >= 0 \select node if connected (outgoing)
y_Rps9b - x_Rps9b_Fun12 >= 0 \select node if connected (outgoing)
y_Sui1 - x_Sui1_Fun12 >= 0 \select node if connected (outgoing)
y_Tif11 - x_Tif11_Fun12 >= 0 \select node if connected (outgoing)
y_Tif35 - x_Tif35_Fun12 >= 0 \select node if connected (outgoing)
y_Tif5 - x_Tif5_Fun12 >= 0 \select node if connected (outgoing)
y_Rvs167 - x_Rvs167_Act1 >= 0 \select node if connected (outgoing)
y_Vrp1 - x_Vrp1_Act1 >= 0 \select node if connected (outgoing)
y_Ssk2 - x_Ssk2_Pbs2 >= 0 \select node if connected (outgoing)
y_Pbs2 - x_Pbs2_Hog1 >= 0 \select node if connected (outgoing)
y_Pbs2 - x_Pbs2_Fus3 >= 0 \select node if connected (outgoing)
y_Pbs2 - x_Pbs2_Kss1 >= 0 \select node if connected (outgoing)
y_Pbs2 - x_Pbs2_Mih1 >= 0 \select node if connected (outgoing)
y_Mih1 - x_Mih1_Swe1 >= 0 \select node if connected (outgoing)
y_Bud4 - x_Bud4_Bud3 >= 0 \select node if connected (outgoing)
y_Bud3 - x_Bud3_Clb2 >= 0 \select node if connected (outgoing)
y_Clb2 - x_Clb2_Swe1 >= 0 \select node if connected (outgoing)
y_Igo1 - x_Igo1_Cdc55 >= 0 \select node if connected (outgoing)
y_Cdc55 - x_Cdc55_Swe1 >= 0 \select node if connected (outgoing)
y_Cdc55 - x_Cdc55_Igo2 >= 0 \select node if connected (outgoing)
y_Syp1 - x_Syp1_Ede1 >= 0 \select node if connected (outgoing)
y_Ypl247c - x_Ypl247c_Yak1 >= 0 \select node if connected (outgoing)
y_Yak1 - x_Yak1_Rim15 >= 0 \select node if connected (outgoing)
y_Yak1 - x_Yak1_Tpk1 >= 0 \select node if connected (outgoing)
y_Rim15 - x_Rim15_Msn2 >= 0 \select node if connected (outgoing)
y_Rim15 - x_Rim15_Sch9 >= 0 \select node if connected (outgoing)
y_Tpk1 - x_Tpk1_Sch9 >= 0 \select node if connected (outgoing)
y_Mtr4 - x_Air2_Mtr4 >= 0 \select node if connected (incoming)
y_Rps11a - x_Mtr4_Rps11a >= 0 \select node if connected (incoming)
y_Rps11a - x_Rps31_Rps11a >= 0 \select node if connected (incoming)
y_Rps11a - x_Ubi4_Rps11a >= 0 \select node if connected (incoming)
y_Rps11b - x_Mtr4_Rps11b >= 0 \select node if connected (incoming)
y_Rps11b - x_Rps31_Rps11b >= 0 \select node if connected (incoming)
y_Rps11b - x_Ubi4_Rps11b >= 0 \select node if connected (incoming)
y_Rps13 - x_Mtr4_Rps13 >= 0 \select node if connected (incoming)
y_Rps13 - x_Rps31_Rps13 >= 0 \select node if connected (incoming)
y_Rps13 - x_Ubi4_Rps13 >= 0 \select node if connected (incoming)
y_Rps14a - x_Mtr4_Rps14a >= 0 \select node if connected (incoming)
y_Rps14a - x_Rps31_Rps14a >= 0 \select node if connected (incoming)
y_Rps14a - x_Ubi4_Rps14a >= 0 \select node if connected (incoming)
y_Rps16a - x_Mtr4_Rps16a >= 0 \select node if connected (incoming)
y_Rps16a - x_Rps31_Rps16a >= 0 \select node if connected (incoming)
y_Rps16a - x_Ubi4_Rps16a >= 0 \select node if connected (incoming)
y_Rps16b - x_Mtr4_Rps16b >= 0 \select node if connected (incoming)
y_Rps16b - x_Rps31_Rps16b >= 0 \select node if connected (incoming)
y_Rps16b - x_Ubi4_Rps16b >= 0 \select node if connected (incoming)
y_Rps1a - x_Mtr4_Rps1a >= 0 \select node if connected (incoming)
y_Rps1a - x_Rps31_Rps1a >= 0 \select node if connected (incoming)
y_Rps1a - x_Ubi4_Rps1a >= 0 \select node if connected (incoming)
y_Rps23a - x_Mtr4_Rps23a >= 0 \select node if connected (incoming)
y_Rps23a - x_Rps31_Rps23a >= 0 \select node if connected (incoming)
y_Rps23a - x_Ubi4_Rps23a >= 0 \select node if connected (incoming)
y_Rps23b - x_Mtr4_Rps23b >= 0 \select node if connected (incoming)
y_Rps23b - x_Rps31_Rps23b >= 0 \select node if connected (incoming)
y_Rps23b - x_Ubi4_Rps23b >= 0 \select node if connected (incoming)
y_Rps24a - x_Mtr4_Rps24a >= 0 \select node if connected (incoming)
y_Rps24a - x_Rps31_Rps24a >= 0 \select node if connected (incoming)
y_Rps24a - x_Ubi4_Rps24a >= 0 \select node if connected (incoming)
y_Rps24b - x_Mtr4_Rps24b >= 0 \select node if connected (incoming)
y_Rps24b - x_Rps31_Rps24b >= 0 \select node if connected (incoming)
y_Rps24b - x_Ubi4_Rps24b >= 0 \select node if connected (incoming)
y_Rps27a - x_Mtr4_Rps27a >= 0 \select node if connected (incoming)
y_Rps27a - x_Rps31_Rps27a >= 0 \select node if connected (incoming)
y_Rps27a - x_Ubi4_Rps27a >= 0 \select node if connected (incoming)
y_Rps28a - x_Mtr4_Rps28a >= 0 \select node if connected (incoming)
y_Rps28a - x_Rps31_Rps28a >= 0 \select node if connected (incoming)
y_Rps28a - x_Ubi4_Rps28a >= 0 \select node if connected (incoming)
y_Rps5 - x_Mtr4_Rps5 >= 0 \select node if connected (incoming)
y_Rps5 - x_Rps31_Rps5 >= 0 \select node if connected (incoming)
y_Rps5 - x_Ubi4_Rps5 >= 0 \select node if connected (incoming)
y_Rps6a - x_Mtr4_Rps6a >= 0 \select node if connected (incoming)
y_Rps6a - x_Rps31_Rps6a >= 0 \select node if connected (incoming)
y_Rps6b - x_Mtr4_Rps6b >= 0 \select node if connected (incoming)
y_Rps6b - x_Rps31_Rps6b >= 0 \select node if connected (incoming)
y_Rps6b - x_Ubi4_Rps6b >= 0 \select node if connected (incoming)
y_Rps7a - x_Mtr4_Rps7a >= 0 \select node if connected (incoming)
y_Rps7a - x_Rps31_Rps7a >= 0 \select node if connected (incoming)
y_Rps7a - x_Ubi4_Rps7a >= 0 \select node if connected (incoming)
y_Rps8a - x_Mtr4_Rps8a >= 0 \select node if connected (incoming)
y_Rps8a - x_Rps31_Rps8a >= 0 \select node if connected (incoming)
y_Rps8a - x_Ubi4_Rps8a >= 0 \select node if connected (incoming)
y_Rps8b - x_Mtr4_Rps8b >= 0 \select node if connected (incoming)
y_Rps8b - x_Rps31_Rps8b >= 0 \select node if connected (incoming)
y_Rps8b - x_Ubi4_Rps8b >= 0 \select node if connected (incoming)
y_Rps9a - x_Mtr4_Rps9a >= 0 \select node if connected (incoming)
y_Rps9a - x_Rps31_Rps9a >= 0 \select node if connected (incoming)
y_Rps9a - x_Ubi4_Rps9a >= 0 \select node if connected (incoming)
y_Npl3 - x_Air2_Npl3 >= 0 \select node if connected (incoming)
y_Rpo21 - x_Npl3_Rpo21 >= 0 \select node if connected (incoming)
y_Rpo21 - x_Rsp5_Rpo21 >= 0 \select node if connected (incoming)
y_Rpo21 - x_Rpb10_Rpo21 >= 0 \select node if connected (incoming)
y_Rpo21 - x_Rpb5_Rpo21 >= 0 \select node if connected (incoming)
y_Rpo21 - x_Rpb8_Rpo21 >= 0 \select node if connected (incoming)
y_Rpo21 - x_Rpc10_Rpo21 >= 0 \select node if connected (incoming)
y_Rpo21 - x_Rpo26_Rpo21 >= 0 \select node if connected (incoming)
y_Rpo21 - x_Rpo31_Rpo21 >= 0 \select node if connected (incoming)
y_Rpo21 - x_Spt15_Rpo21 >= 0 \select node if connected (incoming)
y_Rsp5 - x_Art5_Rsp5 >= 0 \select node if connected (incoming)
y_Rsp5 - x_Ecm21_Rsp5 >= 0 \select node if connected (incoming)
y_Rsp5 - x_Ldb19_Rsp5 >= 0 \select node if connected (incoming)
y_Rps31 - x_Rsp5_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Ubi4_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Vps27_Rps31 >= 0 \select node if connected (incoming)
y_Rps31 - x_Sla1_Rps31 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rsp5_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Rps31_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Ent1_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Ent2_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Swa2_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Vps27_Ubi4 >= 0 \select node if connected (incoming)
y_Ubi4 - x_Sla1_Ubi4 >= 0 \select node if connected (incoming)
y_Fps1 - x_Ask10_Fps1 >= 0 \select node if connected (incoming)
y_Hog1 - x_Fps1_Hog1 >= 0 \select node if connected (incoming)
y_Hog1 - x_Slt2_Hog1 >= 0 \select node if connected (incoming)
y_Hog1 - x_Kss1_Hog1 >= 0 \select node if connected (incoming)
y_Hog1 - x_Pbs2_Hog1 >= 0 \select node if connected (incoming)
y_Slt2 - x_Fps1_Slt2 >= 0 \select node if connected (incoming)
y_Slt2 - x_Msg5_Slt2 >= 0 \select node if connected (incoming)
y_Rpb10 - x_Bdp1_Rpb10 >= 0 \select node if connected (incoming)
y_Cdc73 - x_Rpb10_Cdc73 >= 0 \select node if connected (incoming)
y_Cdc73 - x_Rpb5_Cdc73 >= 0 \select node if connected (incoming)
y_Cdc73 - x_Rpc10_Cdc73 >= 0 \select node if connected (incoming)
y_Cdc73 - x_Rpo26_Cdc73 >= 0 \select node if connected (incoming)
y_Rpb2 - x_Rpb10_Rpb2 >= 0 \select node if connected (incoming)
y_Rpb2 - x_Rpb5_Rpb2 >= 0 \select node if connected (incoming)
y_Rpb2 - x_Rpb8_Rpb2 >= 0 \select node if connected (incoming)
y_Rpb2 - x_Rpc10_Rpb2 >= 0 \select node if connected (incoming)
y_Rpb2 - x_Rpo26_Rpb2 >= 0 \select node if connected (incoming)
y_Rpb2 - x_Spt15_Rpb2 >= 0 \select node if connected (incoming)
y_Rpb3 - x_Rpb10_Rpb3 >= 0 \select node if connected (incoming)
y_Rpb3 - x_Rpb5_Rpb3 >= 0 \select node if connected (incoming)
y_Rpb3 - x_Rpb8_Rpb3 >= 0 \select node if connected (incoming)
y_Rpb3 - x_Rpc10_Rpb3 >= 0 \select node if connected (incoming)
y_Rpb3 - x_Rpo26_Rpb3 >= 0 \select node if connected (incoming)
y_Rpb3 - x_Rpo31_Rpb3 >= 0 \select node if connected (incoming)
y_Rpb3 - x_Spt15_Rpb3 >= 0 \select node if connected (incoming)
y_Spt4 - x_Rpb10_Spt4 >= 0 \select node if connected (incoming)
y_Spt4 - x_Rpb5_Spt4 >= 0 \select node if connected (incoming)
y_Spt4 - x_Rpb8_Spt4 >= 0 \select node if connected (incoming)
y_Spt4 - x_Rpc10_Spt4 >= 0 \select node if connected (incoming)
y_Spt4 - x_Rpo26_Spt4 >= 0 \select node if connected (incoming)
y_Spt5 - x_Rpb10_Spt5 >= 0 \select node if connected (incoming)
y_Spt5 - x_Rpb5_Spt5 >= 0 \select node if connected (incoming)
y_Spt5 - x_Rpb8_Spt5 >= 0 \select node if connected (incoming)
y_Spt5 - x_Rpc10_Spt5 >= 0 \select node if connected (incoming)
y_Spt5 - x_Rpo26_Spt5 >= 0 \select node if connected (incoming)
y_Rpb5 - x_Bdp1_Rpb5 >= 0 \select node if connected (incoming)
y_Leo1 - x_Rpb5_Leo1 >= 0 \select node if connected (incoming)
y_Spt6 - x_Rpb5_Spt6 >= 0 \select node if connected (incoming)
y_Rpb8 - x_Bdp1_Rpb8 >= 0 \select node if connected (incoming)
y_Rpc10 - x_Bdp1_Rpc10 >= 0 \select node if connected (incoming)
y_Rpo26 - x_Bdp1_Rpo26 >= 0 \select node if connected (incoming)
y_Paf1 - x_Rpo26_Paf1 >= 0 \select node if connected (incoming)
y_Rpo31 - x_Bdp1_Rpo31 >= 0 \select node if connected (incoming)
y_Spt15 - x_Bdp1_Spt15 >= 0 \select node if connected (incoming)
y_Apl1 - x_Chc1_Apl1 >= 0 \select node if connected (incoming)
y_Apl3 - x_Apl1_Apl3 >= 0 \select node if connected (incoming)
y_Apl3 - x_Apl2_Apl3 >= 0 \select node if connected (incoming)
y_Apl3 - x_Ent1_Apl3 >= 0 \select node if connected (incoming)
y_Apl3 - x_Ent2_Apl3 >= 0 \select node if connected (incoming)
y_Apl3 - x_Sec16_Apl3 >= 0 \select node if connected (incoming)
y_Apl3 - x_Apm4_Apl3 >= 0 \select node if connected (incoming)
y_Apl2 - x_Chc1_Apl2 >= 0 \select node if connected (incoming)
y_Ent1 - x_Ubi4_Ent1 >= 0 \select node if connected (incoming)
y_Ent1 - x_Chc1_Ent1 >= 0 \select node if connected (incoming)
y_Ent1 - x_Apl3_Ent1 >= 0 \select node if connected (incoming)
y_Ent1 - x_Sla2_Ent1 >= 0 \select node if connected (incoming)
y_Ent2 - x_Ubi4_Ent2 >= 0 \select node if connected (incoming)
y_Ent2 - x_Chc1_Ent2 >= 0 \select node if connected (incoming)
y_Ent2 - x_Apl3_Ent2 >= 0 \select node if connected (incoming)
y_Las17 - x_Chc1_Las17 >= 0 \select node if connected (incoming)
y_Las17 - x_Sla1_Las17 >= 0 \select node if connected (incoming)
y_Act1 - x_Las17_Act1 >= 0 \select node if connected (incoming)
y_Act1 - x_Sla2_Act1 >= 0 \select node if connected (incoming)
y_Act1 - x_Abp1_Act1 >= 0 \select node if connected (incoming)
y_Act1 - x_Rvs167_Act1 >= 0 \select node if connected (incoming)
y_Act1 - x_Vrp1_Act1 >= 0 \select node if connected (incoming)
y_Sla2 - x_Chc1_Sla2 >= 0 \select node if connected (incoming)
y_Sla2 - x_Sla1_Sla2 >= 0 \select node if connected (incoming)
y_Swa2 - x_Chc1_Swa2 >= 0 \select node if connected (incoming)
y_Vps27 - x_Rps31_Vps27 >= 0 \select node if connected (incoming)
y_Vps27 - x_Ubi4_Vps27 >= 0 \select node if connected (incoming)
y_Vps27 - x_Chc1_Vps27 >= 0 \select node if connected (incoming)
y_Fus3 - x_Msg5_Fus3 >= 0 \select node if connected (incoming)
y_Fus3 - x_Kss1_Fus3 >= 0 \select node if connected (incoming)
y_Fus3 - x_Pbs2_Fus3 >= 0 \select node if connected (incoming)
y_Gpa1 - x_Fus3_Gpa1 >= 0 \select node if connected (incoming)
y_Kss1 - x_Msg5_Kss1 >= 0 \select node if connected (incoming)
y_Kss1 - x_Fus3_Kss1 >= 0 \select node if connected (incoming)
y_Kss1 - x_Pbs2_Kss1 >= 0 \select node if connected (incoming)
y_Cdc28 - x_Slt2_Cdc28 >= 0 \select node if connected (incoming)
y_Apm4 - x_Sec16_Apm4 >= 0 \select node if connected (incoming)
y_Abp1 - x_Sla1_Abp1 >= 0 \select node if connected (incoming)
y_Srv2 - x_Abp1_Srv2 >= 0 \select node if connected (incoming)
y_Asc1 - x_Rps31_Asc1 >= 0 \select node if connected (incoming)
y_Asc1 - x_Ubi4_Asc1 >= 0 \select node if connected (incoming)
y_Dom34 - x_Rps31_Dom34 >= 0 \select node if connected (incoming)
y_Gcd11 - x_Rps31_Gcd11 >= 0 \select node if connected (incoming)
y_Nip1 - x_Rps31_Nip1 >= 0 \select node if connected (incoming)
y_Prt1 - x_Rps31_Prt1 >= 0 \select node if connected (incoming)
y_Rpg1 - x_Rps31_Rpg1 >= 0 \select node if connected (incoming)
y_Rpl10 - x_Rps31_Rpl10 >= 0 \select node if connected (incoming)
y_Rpl10 - x_Ubi4_Rpl10 >= 0 \select node if connected (incoming)
y_Rpl17a - x_Rps31_Rpl17a >= 0 \select node if connected (incoming)
y_Rpl17a - x_Ubi4_Rpl17a >= 0 \select node if connected (incoming)
y_Rpl17b - x_Rps31_Rpl17b >= 0 \select node if connected (incoming)
y_Rpl17b - x_Ubi4_Rpl17b >= 0 \select node if connected (incoming)
y_Rpl19b - x_Rps31_Rpl19b >= 0 \select node if connected (incoming)
y_Rpl19b - x_Ubi4_Rpl19b >= 0 \select node if connected (incoming)
y_Rpl25 - x_Rps31_Rpl25 >= 0 \select node if connected (incoming)
y_Rpl25 - x_Ubi4_Rpl25 >= 0 \select node if connected (incoming)
y_Rpl28 - x_Rps31_Rpl28 >= 0 \select node if connected (incoming)
y_Rpl28 - x_Ubi4_Rpl28 >= 0 \select node if connected (incoming)
y_Rpl2a - x_Rps31_Rpl2a >= 0 \select node if connected (incoming)
y_Rpl2a - x_Ubi4_Rpl2a >= 0 \select node if connected (incoming)
y_Rpl2b - x_Rps31_Rpl2b >= 0 \select node if connected (incoming)
y_Rpl2b - x_Ubi4_Rpl2b >= 0 \select node if connected (incoming)
y_Rpl3 - x_Rps31_Rpl3 >= 0 \select node if connected (incoming)
y_Rpl3 - x_Ubi4_Rpl3 >= 0 \select node if connected (incoming)
y_Rpl30 - x_Rps31_Rpl30 >= 0 \select node if connected (incoming)
y_Rpl30 - x_Ubi4_Rpl30 >= 0 \select node if connected (incoming)
y_Rpl32 - x_Rps31_Rpl32 >= 0 \select node if connected (incoming)
y_Rpl32 - x_Ubi4_Rpl32 >= 0 \select node if connected (incoming)
y_Rpl34a - x_Rps31_Rpl34a >= 0 \select node if connected (incoming)
y_Rpl34a - x_Ubi4_Rpl34a >= 0 \select node if connected (incoming)
y_Rpl35a - x_Rps31_Rpl35a >= 0 \select node if connected (incoming)
y_Rpl35a - x_Ubi4_Rpl35a >= 0 \select node if connected (incoming)
y_Rpl35b - x_Rps31_Rpl35b >= 0 \select node if connected (incoming)
y_Rpl35b - x_Ubi4_Rpl35b >= 0 \select node if connected (incoming)
y_Rpl37a - x_Rps31_Rpl37a >= 0 \select node if connected (incoming)
y_Rpl38 - x_Rps31_Rpl38 >= 0 \select node if connected (incoming)
y_Rpl38 - x_Ubi4_Rpl38 >= 0 \select node if connected (incoming)
y_Rpl5 - x_Rps31_Rpl5 >= 0 \select node if connected (incoming)
y_Rpl5 - x_Ubi4_Rpl5 >= 0 \select node if connected (incoming)
y_Rpl7a - x_Rps31_Rpl7a >= 0 \select node if connected (incoming)
y_Rps0a - x_Rps31_Rps0a >= 0 \select node if connected (incoming)
y_Rps0a - x_Ubi4_Rps0a >= 0 \select node if connected (incoming)
y_Rps0b - x_Rps31_Rps0b >= 0 \select node if connected (incoming)
y_Rps0b - x_Ubi4_Rps0b >= 0 \select node if connected (incoming)
y_Rps10a - x_Rps31_Rps10a >= 0 \select node if connected (incoming)
y_Rps10a - x_Ubi4_Rps10a >= 0 \select node if connected (incoming)
y_Rps12 - x_Rps31_Rps12 >= 0 \select node if connected (incoming)
y_Rps12 - x_Ubi4_Rps12 >= 0 \select node if connected (incoming)
y_Rps14b - x_Rps31_Rps14b >= 0 \select node if connected (incoming)
y_Rps14b - x_Ubi4_Rps14b >= 0 \select node if connected (incoming)
y_Rps15 - x_Rps31_Rps15 >= 0 \select node if connected (incoming)
y_Rps15 - x_Ubi4_Rps15 >= 0 \select node if connected (incoming)
y_Rps17a - x_Rps31_Rps17a >= 0 \select node if connected (incoming)
y_Rps17a - x_Ubi4_Rps17a >= 0 \select node if connected (incoming)
y_Rps17b - x_Rps31_Rps17b >= 0 \select node if connected (incoming)
y_Rps17b - x_Ubi4_Rps17b >= 0 \select node if connected (incoming)
y_Rps18a - x_Rps31_Rps18a >= 0 \select node if connected (incoming)
y_Rps18b - x_Rps31_Rps18b >= 0 \select node if connected (incoming)
y_Rps18b - x_Ubi4_Rps18b >= 0 \select node if connected (incoming)
y_Rps19a - x_Rps31_Rps19a >= 0 \select node if connected (incoming)
y_Rps19a - x_Ubi4_Rps19a >= 0 \select node if connected (incoming)
y_Rps1b - x_Rps31_Rps1b >= 0 \select node if connected (incoming)
y_Rps1b - x_Ubi4_Rps1b >= 0 \select node if connected (incoming)
y_Rps2 - x_Rps31_Rps2 >= 0 \select node if connected (incoming)
y_Rps2 - x_Ubi4_Rps2 >= 0 \select node if connected (incoming)
y_Rps20 - x_Rps31_Rps20 >= 0 \select node if connected (incoming)
y_Rps20 - x_Ubi4_Rps20 >= 0 \select node if connected (incoming)
y_Rps22a - x_Rps31_Rps22a >= 0 \select node if connected (incoming)
y_Rps25a - x_Rps31_Rps25a >= 0 \select node if connected (incoming)
y_Rps25a - x_Ubi4_Rps25a >= 0 \select node if connected (incoming)
y_Rps26a - x_Rps31_Rps26a >= 0 \select node if connected (incoming)
y_Rps26a - x_Ubi4_Rps26a >= 0 \select node if connected (incoming)
y_Rps26b - x_Rps31_Rps26b >= 0 \select node if connected (incoming)
y_Rps26b - x_Ubi4_Rps26b >= 0 \select node if connected (incoming)
y_Rps28b - x_Rps31_Rps28b >= 0 \select node if connected (incoming)
y_Rps28b - x_Ubi4_Rps28b >= 0 \select node if connected (incoming)
y_Rps29a - x_Rps31_Rps29a >= 0 \select node if connected (incoming)
y_Rps29a - x_Ubi4_Rps29a >= 0 \select node if connected (incoming)
y_Rps3 - x_Rps31_Rps3 >= 0 \select node if connected (incoming)
y_Rps3 - x_Ubi4_Rps3 >= 0 \select node if connected (incoming)
y_Rps30a - x_Rps31_Rps30a >= 0 \select node if connected (incoming)
y_Rps30b - x_Rps31_Rps30b >= 0 \select node if connected (incoming)
y_Rps9b - x_Rps31_Rps9b >= 0 \select node if connected (incoming)
y_Rps9b - x_Ubi4_Rps9b >= 0 \select node if connected (incoming)
y_Sui1 - x_Rps31_Sui1 >= 0 \select node if connected (incoming)
y_Tif11 - x_Rps31_Tif11 >= 0 \select node if connected (incoming)
y_Tif35 - x_Rps31_Tif35 >= 0 \select node if connected (incoming)
y_Tif5 - x_Rps31_Tif5 >= 0 \select node if connected (incoming)
y_Rvs167 - x_Sla1_Rvs167 >= 0 \select node if connected (incoming)
y_Vrp1 - x_Sla1_Vrp1 >= 0 \select node if connected (incoming)
y_Pbs2 - x_Ssk2_Pbs2 >= 0 \select node if connected (incoming)
y_Mih1 - x_Pbs2_Mih1 >= 0 \select node if connected (incoming)
y_Bud3 - x_Bud4_Bud3 >= 0 \select node if connected (incoming)
y_Clb2 - x_Bud3_Clb2 >= 0 \select node if connected (incoming)
y_Cdc55 - x_Igo1_Cdc55 >= 0 \select node if connected (incoming)
y_Yak1 - x_Ypl247c_Yak1 >= 0 \select node if connected (incoming)
y_Rim15 - x_Yak1_Rim15 >= 0 \select node if connected (incoming)
y_Tpk1 - x_Yak1_Tpk1 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps11a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps11b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps13_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps14a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps16a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps16b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps1a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps23a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps23b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps24a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps24b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps27a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps28a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps5_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps6a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps6b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps7a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps8a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps8b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps9a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps31_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Asc1_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Gcd11_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Nip1_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Prt1_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpg1_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl10_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl17a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl17b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl19b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl25_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl28_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl2a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl2b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl3_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl30_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl32_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl34a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl35a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl35b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl37a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl38_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl5_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rpl7a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps0a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps0b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps10a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps12_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps14b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps15_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps17a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps17b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps18a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps18b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps19a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps1b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps2_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps20_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps22a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps25a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps26a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps26b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps28b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps29a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps3_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps30a_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps30b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Rps9b_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Sui1_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Tif11_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Tif35_Fun12 >= 0 \select node if connected (incoming)
y_Fun12 - x_Tif5_Fun12 >= 0 \select node if connected (incoming)
y_Rtf1 - x_Rpo21_Rtf1 >= 0 \select node if connected (incoming)
y_Rtf1 - x_Cdc73_Rtf1 >= 0 \select node if connected (incoming)
y_Rtf1 - x_Rpb2_Rtf1 >= 0 \select node if connected (incoming)
y_Rtf1 - x_Rpb3_Rtf1 >= 0 \select node if connected (incoming)
y_Rtf1 - x_Spt4_Rtf1 >= 0 \select node if connected (incoming)
y_Rtf1 - x_Spt5_Rtf1 >= 0 \select node if connected (incoming)
y_Rtf1 - x_Leo1_Rtf1 >= 0 \select node if connected (incoming)
y_Rtf1 - x_Spt6_Rtf1 >= 0 \select node if connected (incoming)
y_Rtf1 - x_Paf1_Rtf1 >= 0 \select node if connected (incoming)
y_Ede1 - x_Rps31_Ede1 >= 0 \select node if connected (incoming)
y_Ede1 - x_Ubi4_Ede1 >= 0 \select node if connected (incoming)
y_Ede1 - x_Apl3_Ede1 >= 0 \select node if connected (incoming)
y_Ede1 - x_Ent1_Ede1 >= 0 \select node if connected (incoming)
y_Ede1 - x_Ent2_Ede1 >= 0 \select node if connected (incoming)
y_Ede1 - x_Vps27_Ede1 >= 0 \select node if connected (incoming)
y_Ede1 - x_Syp1_Ede1 >= 0 \select node if connected (incoming)
y_Duf1 - x_Ubi4_Duf1 >= 0 \select node if connected (incoming)
y_Sko1 - x_Hog1_Sko1 >= 0 \select node if connected (incoming)
y_Sko1 - x_Slt2_Sko1 >= 0 \select node if connected (incoming)
y_Sko1 - x_Fus3_Sko1 >= 0 \select node if connected (incoming)
y_Sko1 - x_Kss1_Sko1 >= 0 \select node if connected (incoming)
y_Cyr1 - x_Act1_Cyr1 >= 0 \select node if connected (incoming)
y_Cyr1 - x_Gpa1_Cyr1 >= 0 \select node if connected (incoming)
y_Cyr1 - x_Srv2_Cyr1 >= 0 \select node if connected (incoming)
y_Swe1 - x_Cdc28_Swe1 >= 0 \select node if connected (incoming)
y_Swe1 - x_Mih1_Swe1 >= 0 \select node if connected (incoming)
y_Swe1 - x_Clb2_Swe1 >= 0 \select node if connected (incoming)
y_Swe1 - x_Cdc55_Swe1 >= 0 \select node if connected (incoming)
y_Sup35 - x_Dom34_Sup35 >= 0 \select node if connected (incoming)
y_Igo2 - x_Cdc55_Igo2 >= 0 \select node if connected (incoming)
y_Msn2 - x_Rim15_Msn2 >= 0 \select node if connected (incoming)
y_Sch9 - x_Rim15_Sch9 >= 0 \select node if connected (incoming)
y_Sch9 - x_Tpk1_Sch9 >= 0 \select node if connected (incoming)
x_Rps31_Ubi4 + x_Ubi4_Rps31 <= 1 \avoid 2-node cycles
x_Rps31_Vps27 + x_Vps27_Rps31 <= 1 \avoid 2-node cycles
x_Rps31_Ubi4 + x_Ubi4_Vps27 + x_Vps27_Rps31 <= 2 \avoid 3-node cycles
x_Rps31_Ubi4 + x_Ubi4_Ent1 + x_Ent1_Ubi4 + x_Ubi4_Rps31 <= 3 \avoid 4-node cycles
x_Rps31_Ubi4 + x_Ubi4_Ent1 + x_Ent1_Apl3 + x_Apl3_Ent2 + x_Ent2_Ubi4 + x_Ubi4_Rps31 <= 5 \avoid 6-node cycles
x_Rps31_Ubi4 + x_Ubi4_Ent2 + x_Ent2_Ubi4 + x_Ubi4_Rps31 <= 3 \avoid 4-node cycles
x_Rps31_Ubi4 + x_Ubi4_Ent2 + x_Ent2_Apl3 + x_Apl3_Ent1 + x_Ent1_Ubi4 + x_Ubi4_Rps31 <= 5 \avoid 6-node cycles
x_Rps31_Ubi4 + x_Ubi4_Vps27 + x_Vps27_Ubi4 + x_Ubi4_Rps31 <= 3 \avoid 4-node cycles
x_Rps31_Vps27 + x_Vps27_Ubi4 + x_Ubi4_Rps31 <= 2 \avoid 3-node cycles
x_Rps31_Vps27 + x_Vps27_Ubi4 + x_Ubi4_Vps27 + x_Vps27_Rps31 <= 3 \avoid 4-node cycles
x_Rps31_Vps27 + x_Vps27_Ubi4 + x_Ubi4_Ent1 + x_Ent1_Ubi4 + x_Ubi4_Rps31 <= 4 \avoid 5-node cycles
x_Rps31_Vps27 + x_Vps27_Ubi4 + x_Ubi4_Ent2 + x_Ent2_Ubi4 + x_Ubi4_Rps31 <= 4 \avoid 5-node cycles
x_Ubi4_Ent1 + x_Ent1_Ubi4 <= 1 \avoid 2-node cycles
x_Ubi4_Ent2 + x_Ent2_Ubi4 <= 1 \avoid 2-node cycles
x_Ubi4_Vps27 + x_Vps27_Ubi4 <= 1 \avoid 2-node cycles
x_Ubi4_Ent1 + x_Ent1_Apl3 + x_Apl3_Ent1 + x_Ent1_Ubi4 <= 3 \avoid 4-node cycles
x_Ubi4_Ent1 + x_Ent1_Apl3 + x_Apl3_Ent2 + x_Ent2_Ubi4 <= 3 \avoid 4-node cycles
x_Ubi4_Ent2 + x_Ent2_Apl3 + x_Apl3_Ent1 + x_Ent1_Ubi4 <= 3 \avoid 4-node cycles
x_Ubi4_Ent2 + x_Ent2_Apl3 + x_Apl3_Ent2 + x_Ent2_Ubi4 <= 3 \avoid 4-node cycles
x_Apl3_Ent1 + x_Ent1_Apl3 <= 1 \avoid 2-node cycles
x_Apl3_Ent2 + x_Ent2_Apl3 <= 1 \avoid 2-node cycles
x_Fus3_Kss1 + x_Kss1_Fus3 <= 1 \avoid 2-node cycles

binary
y_Air2
y_Mtr4
y_Rps11a
y_Rps11b
y_Rps13
y_Rps14a
y_Rps16a
y_Rps16b
y_Rps1a
y_Rps23a
y_Rps23b
y_Rps24a
y_Rps24b
y_Rps27a
y_Rps28a
y_Rps5
y_Rps6a
y_Rps6b
y_Rps7a
y_Rps8a
y_Rps8b
y_Rps9a
y_Npl3
y_Rpo21
y_Art5
y_Rsp5
y_Rps31
y_Ubi4
y_Ask10
y_Fps1
y_Hog1
y_Slt2
y_Bdp1
y_Rpb10
y_Cdc73
y_Rpb2
y_Rpb3
y_Spt4
y_Spt5
y_Rpb5
y_Leo1
y_Spt6
y_Rpb8
y_Rpc10
y_Rpo26
y_Paf1
y_Rpo31
y_Spt15
y_Chc1
y_Apl1
y_Apl3
y_Apl2
y_Ent1
y_Ent2
y_Las17
y_Act1
y_Sla2
y_Swa2
y_Vps27
y_Ecm21
y_Ldb19
y_Msg5
y_Fus3
y_Gpa1
y_Kss1
y_Cdc28
y_Sec16
y_Apm4
y_Sla1
y_Abp1
y_Srv2
y_Asc1
y_Dom34
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
y_Rps12
y_Rps14b
y_Rps15
y_Rps17a
y_Rps17b
y_Rps18a
y_Rps18b
y_Rps19a
y_Rps1b
y_Rps2
y_Rps20
y_Rps22a
y_Rps25a
y_Rps26a
y_Rps26b
y_Rps28b
y_Rps29a
y_Rps3
y_Rps30a
y_Rps30b
y_Rps9b
y_Sui1
y_Tif11
y_Tif35
y_Tif5
y_Rvs167
y_Vrp1
y_Ssk2
y_Pbs2
y_Mih1
y_Bud4
y_Bud3
y_Clb2
y_Igo1
y_Cdc55
y_Syp1
y_Ypl247c
y_Yak1
y_Rim15
y_Tpk1
y_Fun12
y_Rtf1
y_Ede1
y_Duf1
y_Sko1
y_Cyr1
y_Swe1
y_Sup35
y_Igo2
y_Msn2
y_Sch9
x_Air2_Mtr4
x_Air2_Npl3
x_Mtr4_Rps11a
x_Mtr4_Rps11b
x_Mtr4_Rps13
x_Mtr4_Rps14a
x_Mtr4_Rps16a
x_Mtr4_Rps16b
x_Mtr4_Rps1a
x_Mtr4_Rps23a
x_Mtr4_Rps23b
x_Mtr4_Rps24a
x_Mtr4_Rps24b
x_Mtr4_Rps27a
x_Mtr4_Rps28a
x_Mtr4_Rps5
x_Mtr4_Rps6a
x_Mtr4_Rps6b
x_Mtr4_Rps7a
x_Mtr4_Rps8a
x_Mtr4_Rps8b
x_Mtr4_Rps9a
x_Rps11a_Fun12
x_Rps11b_Fun12
x_Rps13_Fun12
x_Rps14a_Fun12
x_Rps16a_Fun12
x_Rps16b_Fun12
x_Rps1a_Fun12
x_Rps23a_Fun12
x_Rps23b_Fun12
x_Rps24a_Fun12
x_Rps24b_Fun12
x_Rps27a_Fun12
x_Rps28a_Fun12
x_Rps5_Fun12
x_Rps6a_Fun12
x_Rps6b_Fun12
x_Rps7a_Fun12
x_Rps8a_Fun12
x_Rps8b_Fun12
x_Rps9a_Fun12
x_Npl3_Rpo21
x_Rpo21_Rtf1
x_Art5_Rsp5
x_Rsp5_Rpo21
x_Rsp5_Rps31
x_Rsp5_Ubi4
x_Rps31_Rps11a
x_Rps31_Rps11b
x_Rps31_Rps13
x_Rps31_Rps14a
x_Rps31_Rps16a
x_Rps31_Rps16b
x_Rps31_Rps1a
x_Rps31_Rps23a
x_Rps31_Rps23b
x_Rps31_Rps24a
x_Rps31_Rps24b
x_Rps31_Rps27a
x_Rps31_Rps28a
x_Rps31_Rps5
x_Rps31_Rps6a
x_Rps31_Rps6b
x_Rps31_Rps7a
x_Rps31_Rps8a
x_Rps31_Rps8b
x_Rps31_Rps9a
x_Rps31_Ubi4
x_Rps31_Vps27
x_Rps31_Asc1
x_Rps31_Dom34
x_Rps31_Gcd11
x_Rps31_Nip1
x_Rps31_Prt1
x_Rps31_Rpg1
x_Rps31_Rpl10
x_Rps31_Rpl17a
x_Rps31_Rpl17b
x_Rps31_Rpl19b
x_Rps31_Rpl25
x_Rps31_Rpl28
x_Rps31_Rpl2a
x_Rps31_Rpl2b
x_Rps31_Rpl3
x_Rps31_Rpl30
x_Rps31_Rpl32
x_Rps31_Rpl34a
x_Rps31_Rpl35a
x_Rps31_Rpl35b
x_Rps31_Rpl37a
x_Rps31_Rpl38
x_Rps31_Rpl5
x_Rps31_Rpl7a
x_Rps31_Rps0a
x_Rps31_Rps0b
x_Rps31_Rps10a
x_Rps31_Rps12
x_Rps31_Rps14b
x_Rps31_Rps15
x_Rps31_Rps17a
x_Rps31_Rps17b
x_Rps31_Rps18a
x_Rps31_Rps18b
x_Rps31_Rps19a
x_Rps31_Rps1b
x_Rps31_Rps2
x_Rps31_Rps20
x_Rps31_Rps22a
x_Rps31_Rps25a
x_Rps31_Rps26a
x_Rps31_Rps26b
x_Rps31_Rps28b
x_Rps31_Rps29a
x_Rps31_Rps3
x_Rps31_Rps30a
x_Rps31_Rps30b
x_Rps31_Rps9b
x_Rps31_Sui1
x_Rps31_Tif11
x_Rps31_Tif35
x_Rps31_Tif5
x_Rps31_Fun12
x_Rps31_Ede1
x_Ubi4_Rps11a
x_Ubi4_Rps11b
x_Ubi4_Rps13
x_Ubi4_Rps14a
x_Ubi4_Rps16a
x_Ubi4_Rps16b
x_Ubi4_Rps1a
x_Ubi4_Rps23a
x_Ubi4_Rps23b
x_Ubi4_Rps24a
x_Ubi4_Rps24b
x_Ubi4_Rps27a
x_Ubi4_Rps28a
x_Ubi4_Rps5
x_Ubi4_Rps6b
x_Ubi4_Rps7a
x_Ubi4_Rps8a
x_Ubi4_Rps8b
x_Ubi4_Rps9a
x_Ubi4_Rps31
x_Ubi4_Ent1
x_Ubi4_Ent2
x_Ubi4_Vps27
x_Ubi4_Asc1
x_Ubi4_Rpl10
x_Ubi4_Rpl17a
x_Ubi4_Rpl17b
x_Ubi4_Rpl19b
x_Ubi4_Rpl25
x_Ubi4_Rpl28
x_Ubi4_Rpl2a
x_Ubi4_Rpl2b
x_Ubi4_Rpl3
x_Ubi4_Rpl30
x_Ubi4_Rpl32
x_Ubi4_Rpl34a
x_Ubi4_Rpl35a
x_Ubi4_Rpl35b
x_Ubi4_Rpl38
x_Ubi4_Rpl5
x_Ubi4_Rps0a
x_Ubi4_Rps0b
x_Ubi4_Rps10a
x_Ubi4_Rps12
x_Ubi4_Rps14b
x_Ubi4_Rps15
x_Ubi4_Rps17a
x_Ubi4_Rps17b
x_Ubi4_Rps18b
x_Ubi4_Rps19a
x_Ubi4_Rps1b
x_Ubi4_Rps2
x_Ubi4_Rps20
x_Ubi4_Rps25a
x_Ubi4_Rps26a
x_Ubi4_Rps26b
x_Ubi4_Rps28b
x_Ubi4_Rps29a
x_Ubi4_Rps3
x_Ubi4_Rps9b
x_Ubi4_Ede1
x_Ubi4_Duf1
x_Ask10_Fps1
x_Fps1_Hog1
x_Fps1_Slt2
x_Hog1_Sko1
x_Slt2_Hog1
x_Slt2_Cdc28
x_Slt2_Sko1
x_Bdp1_Rpb10
x_Bdp1_Rpb5
x_Bdp1_Rpb8
x_Bdp1_Rpc10
x_Bdp1_Rpo26
x_Bdp1_Rpo31
x_Bdp1_Spt15
x_Rpb10_Rpo21
x_Rpb10_Cdc73
x_Rpb10_Rpb2
x_Rpb10_Rpb3
x_Rpb10_Spt4
x_Rpb10_Spt5
x_Cdc73_Rtf1
x_Rpb2_Rtf1
x_Rpb3_Rtf1
x_Spt4_Rtf1
x_Spt5_Rtf1
x_Rpb5_Rpo21
x_Rpb5_Cdc73
x_Rpb5_Rpb2
x_Rpb5_Rpb3
x_Rpb5_Spt4
x_Rpb5_Spt5
x_Rpb5_Leo1
x_Rpb5_Spt6
x_Leo1_Rtf1
x_Spt6_Rtf1
x_Rpb8_Rpo21
x_Rpb8_Rpb2
x_Rpb8_Rpb3
x_Rpb8_Spt4
x_Rpb8_Spt5
x_Rpc10_Rpo21
x_Rpc10_Cdc73
x_Rpc10_Rpb2
x_Rpc10_Rpb3
x_Rpc10_Spt4
x_Rpc10_Spt5
x_Rpo26_Rpo21
x_Rpo26_Cdc73
x_Rpo26_Rpb2
x_Rpo26_Rpb3
x_Rpo26_Spt4
x_Rpo26_Spt5
x_Rpo26_Paf1
x_Paf1_Rtf1
x_Rpo31_Rpo21
x_Rpo31_Rpb3
x_Spt15_Rpo21
x_Spt15_Rpb2
x_Spt15_Rpb3
x_Chc1_Apl1
x_Chc1_Apl2
x_Chc1_Ent1
x_Chc1_Ent2
x_Chc1_Las17
x_Chc1_Sla2
x_Chc1_Swa2
x_Chc1_Vps27
x_Apl1_Apl3
x_Apl3_Ent1
x_Apl3_Ent2
x_Apl3_Ede1
x_Apl2_Apl3
x_Ent1_Ubi4
x_Ent1_Apl3
x_Ent1_Ede1
x_Ent2_Ubi4
x_Ent2_Apl3
x_Ent2_Ede1
x_Las17_Act1
x_Act1_Cyr1
x_Sla2_Ent1
x_Sla2_Act1
x_Swa2_Ubi4
x_Vps27_Rps31
x_Vps27_Ubi4
x_Vps27_Ede1
x_Ecm21_Rsp5
x_Ldb19_Rsp5
x_Msg5_Slt2
x_Msg5_Fus3
x_Msg5_Kss1
x_Fus3_Gpa1
x_Fus3_Kss1
x_Fus3_Sko1
x_Gpa1_Cyr1
x_Kss1_Hog1
x_Kss1_Fus3
x_Kss1_Sko1
x_Cdc28_Swe1
x_Sec16_Apl3
x_Sec16_Apm4
x_Apm4_Apl3
x_Sla1_Rps31
x_Sla1_Ubi4
x_Sla1_Las17
x_Sla1_Sla2
x_Sla1_Abp1
x_Sla1_Rvs167
x_Sla1_Vrp1
x_Abp1_Act1
x_Abp1_Srv2
x_Srv2_Cyr1
x_Asc1_Fun12
x_Dom34_Sup35
x_Gcd11_Fun12
x_Nip1_Fun12
x_Prt1_Fun12
x_Rpg1_Fun12
x_Rpl10_Fun12
x_Rpl17a_Fun12
x_Rpl17b_Fun12
x_Rpl19b_Fun12
x_Rpl25_Fun12
x_Rpl28_Fun12
x_Rpl2a_Fun12
x_Rpl2b_Fun12
x_Rpl3_Fun12
x_Rpl30_Fun12
x_Rpl32_Fun12
x_Rpl34a_Fun12
x_Rpl35a_Fun12
x_Rpl35b_Fun12
x_Rpl37a_Fun12
x_Rpl38_Fun12
x_Rpl5_Fun12
x_Rpl7a_Fun12
x_Rps0a_Fun12
x_Rps0b_Fun12
x_Rps10a_Fun12
x_Rps12_Fun12
x_Rps14b_Fun12
x_Rps15_Fun12
x_Rps17a_Fun12
x_Rps17b_Fun12
x_Rps18a_Fun12
x_Rps18b_Fun12
x_Rps19a_Fun12
x_Rps1b_Fun12
x_Rps2_Fun12
x_Rps20_Fun12
x_Rps22a_Fun12
x_Rps25a_Fun12
x_Rps26a_Fun12
x_Rps26b_Fun12
x_Rps28b_Fun12
x_Rps29a_Fun12
x_Rps3_Fun12
x_Rps30a_Fun12
x_Rps30b_Fun12
x_Rps9b_Fun12
x_Sui1_Fun12
x_Tif11_Fun12
x_Tif35_Fun12
x_Tif5_Fun12
x_Rvs167_Act1
x_Vrp1_Act1
x_Ssk2_Pbs2
x_Pbs2_Hog1
x_Pbs2_Fus3
x_Pbs2_Kss1
x_Pbs2_Mih1
x_Mih1_Swe1
x_Bud4_Bud3
x_Bud3_Clb2
x_Clb2_Swe1
x_Igo1_Cdc55
x_Cdc55_Swe1
x_Cdc55_Igo2
x_Syp1_Ede1
x_Ypl247c_Yak1
x_Yak1_Rim15
x_Yak1_Tpk1
x_Rim15_Msn2
x_Rim15_Sch9
x_Tpk1_Sch9

End