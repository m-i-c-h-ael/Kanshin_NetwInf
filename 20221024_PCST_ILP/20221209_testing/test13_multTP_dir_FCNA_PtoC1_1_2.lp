\This implementation avoids singly-connected nodes that therefore do not
   \contribute to connecting terminal nodes, even if they contribute positive prize
   \it can produce unconnected subgraphs as long as each subgraph has at least two terminals
\Decision variables: edge-selected: x_ij; node-selected: y_j

Minimize 5 x_START1_MMS1 + 9 x_MMS1_MRS6 + 5 x_MRS6_MTR2 + 5 x_MTR2_PSD2 + 10 x_PSD2_END1 + 8 x_MTR2_VID30 + 4 x_VID30_END1 + 7 x_MTR2_END2 + 10 x_END2_END1 + 8 x_MRS6_PGM3 + 8 x_PGM3_BIT2 + 6 x_BIT2_END1 + 10 x_MRS6_SNU66 + 7 x_SNU66_PSD2 + 10 x_SNU66_SPE4 + 4 x_SPE4_END1 + 5 x_SNU66_END2 + 2 x_START1_SFP1 + 5 x_SFP1_RSC9 + 9 x_RSC9_PGM3 + 3 x_RSC9_VPS53 + 5 x_VPS53_END2 + 10 x_SFP1_SKP1 + 4 x_SKP1_MTR2 + 9 x_SKP1_PGM3 + 8 x_SKP1_SNU66 + 5 x_START2_MMS1 + 4 x_START2_SFP1 + 4 x_START2_YGL140C + 8 x_YGL140C_SKP1 + 10 x_YGL140C_YOL087C + 5 x_YOL087C_MTR2 + 5 x_YOL087C_PGM3 + 9 x_START2_YNL157W + 5 x_YNL157W_ART5 + 4 x_ART5_PGM3 + 6 x_ART5_SNU66 + 2 x_ART5_VPS53 + 4 x_YNL157W_MRS6 - 9 y_START1 - 9 y_MMS1 - 1 y_MRS6 - 10 y_MTR2 - 2 y_PSD2 - 2 y_VID30 - 9 y_END2 - 4 y_PGM3 - 10 y_BIT2 - 4 y_SNU66 - 9 y_SPE4 - 7 y_SFP1 - 4 y_RSC9 - 8 y_VPS53 - 6 y_SKP1 - 7 y_START2 - 5 y_YGL140C - 10 y_YOL087C - 9 y_YNL157W - 8 y_ART5 - 3 y_END1

subject to
x_START1_MMS1 + x_START1_SFP1 >= 1 \edges from terminal
x_START2_MMS1 + x_START2_SFP1 + x_START2_YGL140C + x_START2_YNL157W >= 1 \edges from terminal
x_MTR2_END2 + x_SNU66_END2 + x_VPS53_END2 >= 1 \edges to terminal
x_PSD2_END1 + x_VID30_END1 + x_END2_END1 + x_BIT2_END1 + x_SPE4_END1 >= 1 \edges to terminal
y_START1 >= 1  \terminal
y_START2 >= 1  \terminal
y_END2 >= 1  \terminal
y_END1 >= 1  \terminal
y_MMS1 - x_MMS1_MRS6 <= 0 \Steiner node out
y_MRS6 - x_MRS6_MTR2 - x_MRS6_PGM3 - x_MRS6_SNU66 <= 0 \Steiner node out
y_MTR2 - x_MTR2_PSD2 - x_MTR2_VID30 - x_MTR2_END2 <= 0 \Steiner node out
y_PSD2 - x_PSD2_END1 <= 0 \Steiner node out
y_VID30 - x_VID30_END1 <= 0 \Steiner node out
y_PGM3 - x_PGM3_BIT2 <= 0 \Steiner node out
y_BIT2 - x_BIT2_END1 <= 0 \Steiner node out
y_SNU66 - x_SNU66_PSD2 - x_SNU66_END2 - x_SNU66_SPE4 <= 0 \Steiner node out
y_SPE4 - x_SPE4_END1 <= 0 \Steiner node out
y_SFP1 - x_SFP1_RSC9 - x_SFP1_SKP1 <= 0 \Steiner node out
y_RSC9 - x_RSC9_PGM3 - x_RSC9_VPS53 <= 0 \Steiner node out
y_VPS53 - x_VPS53_END2 <= 0 \Steiner node out
y_SKP1 - x_SKP1_MTR2 - x_SKP1_PGM3 - x_SKP1_SNU66 <= 0 \Steiner node out
y_YGL140C - x_YGL140C_SKP1 - x_YGL140C_YOL087C <= 0 \Steiner node out
y_YOL087C - x_YOL087C_MTR2 - x_YOL087C_PGM3 <= 0 \Steiner node out
y_YNL157W - x_YNL157W_MRS6 - x_YNL157W_ART5 <= 0 \Steiner node out
y_ART5 - x_ART5_PGM3 - x_ART5_SNU66 - x_ART5_VPS53 <= 0 \Steiner node out
y_MMS1 - x_START1_MMS1 - x_START2_MMS1 <= 0 \Steiner node in
y_MRS6 - x_MMS1_MRS6 - x_YNL157W_MRS6 <= 0 \Steiner node in
y_MTR2 - x_MRS6_MTR2 - x_SKP1_MTR2 - x_YOL087C_MTR2 <= 0 \Steiner node in
y_PSD2 - x_MTR2_PSD2 - x_SNU66_PSD2 <= 0 \Steiner node in
y_VID30 - x_MTR2_VID30 <= 0 \Steiner node in
y_PGM3 - x_MRS6_PGM3 - x_RSC9_PGM3 - x_SKP1_PGM3 - x_YOL087C_PGM3 - x_ART5_PGM3 <= 0 \Steiner node in
y_BIT2 - x_PGM3_BIT2 <= 0 \Steiner node in
y_SNU66 - x_MRS6_SNU66 - x_SKP1_SNU66 - x_ART5_SNU66 <= 0 \Steiner node in
y_SPE4 - x_SNU66_SPE4 <= 0 \Steiner node in
y_SFP1 - x_START1_SFP1 - x_START2_SFP1 <= 0 \Steiner node in
y_RSC9 - x_SFP1_RSC9 <= 0 \Steiner node in
y_VPS53 - x_RSC9_VPS53 - x_ART5_VPS53 <= 0 \Steiner node in
y_SKP1 - x_SFP1_SKP1 - x_YGL140C_SKP1 <= 0 \Steiner node in
y_YGL140C - x_START2_YGL140C <= 0 \Steiner node in
y_YOL087C - x_YGL140C_YOL087C <= 0 \Steiner node in
y_YNL157W - x_START2_YNL157W <= 0 \Steiner node in
y_ART5 - x_YNL157W_ART5 <= 0 \Steiner node in
y_START1 - x_START1_MMS1 >= 0 \select node if connected (outgoing)
y_START1 - x_START1_SFP1 >= 0 \select node if connected (outgoing)
y_MMS1 - x_MMS1_MRS6 >= 0 \select node if connected (outgoing)
y_MRS6 - x_MRS6_MTR2 >= 0 \select node if connected (outgoing)
y_MRS6 - x_MRS6_PGM3 >= 0 \select node if connected (outgoing)
y_MRS6 - x_MRS6_SNU66 >= 0 \select node if connected (outgoing)
y_MTR2 - x_MTR2_PSD2 >= 0 \select node if connected (outgoing)
y_MTR2 - x_MTR2_VID30 >= 0 \select node if connected (outgoing)
y_MTR2 - x_MTR2_END2 >= 0 \select node if connected (outgoing)
y_PSD2 - x_PSD2_END1 >= 0 \select node if connected (outgoing)
y_VID30 - x_VID30_END1 >= 0 \select node if connected (outgoing)
y_END2 - x_END2_END1 >= 0 \select node if connected (outgoing)
y_PGM3 - x_PGM3_BIT2 >= 0 \select node if connected (outgoing)
y_BIT2 - x_BIT2_END1 >= 0 \select node if connected (outgoing)
y_SNU66 - x_SNU66_PSD2 >= 0 \select node if connected (outgoing)
y_SNU66 - x_SNU66_END2 >= 0 \select node if connected (outgoing)
y_SNU66 - x_SNU66_SPE4 >= 0 \select node if connected (outgoing)
y_SPE4 - x_SPE4_END1 >= 0 \select node if connected (outgoing)
y_SFP1 - x_SFP1_RSC9 >= 0 \select node if connected (outgoing)
y_SFP1 - x_SFP1_SKP1 >= 0 \select node if connected (outgoing)
y_RSC9 - x_RSC9_PGM3 >= 0 \select node if connected (outgoing)
y_RSC9 - x_RSC9_VPS53 >= 0 \select node if connected (outgoing)
y_VPS53 - x_VPS53_END2 >= 0 \select node if connected (outgoing)
y_SKP1 - x_SKP1_MTR2 >= 0 \select node if connected (outgoing)
y_SKP1 - x_SKP1_PGM3 >= 0 \select node if connected (outgoing)
y_SKP1 - x_SKP1_SNU66 >= 0 \select node if connected (outgoing)
y_START2 - x_START2_MMS1 >= 0 \select node if connected (outgoing)
y_START2 - x_START2_SFP1 >= 0 \select node if connected (outgoing)
y_START2 - x_START2_YGL140C >= 0 \select node if connected (outgoing)
y_START2 - x_START2_YNL157W >= 0 \select node if connected (outgoing)
y_YGL140C - x_YGL140C_SKP1 >= 0 \select node if connected (outgoing)
y_YGL140C - x_YGL140C_YOL087C >= 0 \select node if connected (outgoing)
y_YOL087C - x_YOL087C_MTR2 >= 0 \select node if connected (outgoing)
y_YOL087C - x_YOL087C_PGM3 >= 0 \select node if connected (outgoing)
y_YNL157W - x_YNL157W_MRS6 >= 0 \select node if connected (outgoing)
y_YNL157W - x_YNL157W_ART5 >= 0 \select node if connected (outgoing)
y_ART5 - x_ART5_PGM3 >= 0 \select node if connected (outgoing)
y_ART5 - x_ART5_SNU66 >= 0 \select node if connected (outgoing)
y_ART5 - x_ART5_VPS53 >= 0 \select node if connected (outgoing)
y_MMS1 - x_START1_MMS1 >= 0 \select node if connected (incoming)
y_MMS1 - x_START2_MMS1 >= 0 \select node if connected (incoming)
y_MRS6 - x_MMS1_MRS6 >= 0 \select node if connected (incoming)
y_MRS6 - x_YNL157W_MRS6 >= 0 \select node if connected (incoming)
y_MTR2 - x_MRS6_MTR2 >= 0 \select node if connected (incoming)
y_MTR2 - x_SKP1_MTR2 >= 0 \select node if connected (incoming)
y_MTR2 - x_YOL087C_MTR2 >= 0 \select node if connected (incoming)
y_PSD2 - x_MTR2_PSD2 >= 0 \select node if connected (incoming)
y_PSD2 - x_SNU66_PSD2 >= 0 \select node if connected (incoming)
y_VID30 - x_MTR2_VID30 >= 0 \select node if connected (incoming)
y_END2 - x_MTR2_END2 >= 0 \select node if connected (incoming)
y_END2 - x_SNU66_END2 >= 0 \select node if connected (incoming)
y_END2 - x_VPS53_END2 >= 0 \select node if connected (incoming)
y_PGM3 - x_MRS6_PGM3 >= 0 \select node if connected (incoming)
y_PGM3 - x_RSC9_PGM3 >= 0 \select node if connected (incoming)
y_PGM3 - x_SKP1_PGM3 >= 0 \select node if connected (incoming)
y_PGM3 - x_YOL087C_PGM3 >= 0 \select node if connected (incoming)
y_PGM3 - x_ART5_PGM3 >= 0 \select node if connected (incoming)
y_BIT2 - x_PGM3_BIT2 >= 0 \select node if connected (incoming)
y_SNU66 - x_MRS6_SNU66 >= 0 \select node if connected (incoming)
y_SNU66 - x_SKP1_SNU66 >= 0 \select node if connected (incoming)
y_SNU66 - x_ART5_SNU66 >= 0 \select node if connected (incoming)
y_SPE4 - x_SNU66_SPE4 >= 0 \select node if connected (incoming)
y_SFP1 - x_START1_SFP1 >= 0 \select node if connected (incoming)
y_SFP1 - x_START2_SFP1 >= 0 \select node if connected (incoming)
y_RSC9 - x_SFP1_RSC9 >= 0 \select node if connected (incoming)
y_VPS53 - x_RSC9_VPS53 >= 0 \select node if connected (incoming)
y_VPS53 - x_ART5_VPS53 >= 0 \select node if connected (incoming)
y_SKP1 - x_SFP1_SKP1 >= 0 \select node if connected (incoming)
y_SKP1 - x_YGL140C_SKP1 >= 0 \select node if connected (incoming)
y_YGL140C - x_START2_YGL140C >= 0 \select node if connected (incoming)
y_YOL087C - x_YGL140C_YOL087C >= 0 \select node if connected (incoming)
y_YNL157W - x_START2_YNL157W >= 0 \select node if connected (incoming)
y_ART5 - x_YNL157W_ART5 >= 0 \select node if connected (incoming)
y_END1 - x_PSD2_END1 >= 0 \select node if connected (incoming)
y_END1 - x_VID30_END1 >= 0 \select node if connected (incoming)
y_END1 - x_END2_END1 >= 0 \select node if connected (incoming)
y_END1 - x_BIT2_END1 >= 0 \select node if connected (incoming)
y_END1 - x_SPE4_END1 >= 0 \select node if connected (incoming)

binary
y_START1
y_MMS1
y_MRS6
y_MTR2
y_PSD2
y_VID30
y_END2
y_PGM3
y_BIT2
y_SNU66
y_SPE4
y_SFP1
y_RSC9
y_VPS53
y_SKP1
y_START2
y_YGL140C
y_YOL087C
y_YNL157W
y_ART5
y_END1
x_START1_MMS1
x_START1_SFP1
x_MMS1_MRS6
x_MRS6_MTR2
x_MRS6_PGM3
x_MRS6_SNU66
x_MTR2_PSD2
x_MTR2_VID30
x_MTR2_END2
x_PSD2_END1
x_VID30_END1
x_END2_END1
x_PGM3_BIT2
x_BIT2_END1
x_SNU66_PSD2
x_SNU66_END2
x_SNU66_SPE4
x_SPE4_END1
x_SFP1_RSC9
x_SFP1_SKP1
x_RSC9_PGM3
x_RSC9_VPS53
x_VPS53_END2
x_SKP1_MTR2
x_SKP1_PGM3
x_SKP1_SNU66
x_START2_MMS1
x_START2_SFP1
x_START2_YGL140C
x_START2_YNL157W
x_YGL140C_SKP1
x_YGL140C_YOL087C
x_YOL087C_MTR2
x_YOL087C_PGM3
x_YNL157W_MRS6
x_YNL157W_ART5
x_ART5_PGM3
x_ART5_SNU66
x_ART5_VPS53

End