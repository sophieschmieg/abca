package fields.tests;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

class MysteryStringTest {

	@Test
	void test() {
		String[] strings = { "F Ci 27mta Kod 38ztm Gndn-y 8N Nzrlt 6 Sndq Gztz ztr-ngq Qm 9a Mo 9e On Me Jntu Nntu",
				"D2ei 2mgq Jz 9b-mm Gm Pq Ry LN Nnwm Olt 7ywi Gmt-Kndr 9otq Rywv 8o 9e Pmti Nmd 2 Sn 92 Tma",
				"Xl E Dum Mz 7n M7-m 9i Gogm Rm LN Nyt 8lwi Kz 9e Gm 9-Pmv 7zti Lztz o 9e On Me Qnd-Sodm",
				"F Ci 27zgf 8mdq Mm Me Gn Mm My 8N Nz 9e Qlwe Nye Gm Mi Mm 96 Qmgr 9n Mb-o 9e Pmtu Rmt 6 Jotm Tma",
				"D2U Udywi Jmtz 7zt-Gm 9e Qmc N Nzt 2 Qlwf 7m 9u Gzd 7zdf 7owr 9y Mq Oo 9e Pmta Kn M2 Nmdu Tma",
				"Tpr Dsn Me Jn 9i Onwe Gn Pu Qns N Nze Mlt-Qmm Gotzyti Nzy Rmd-Mo 9e On M6 Pot 2 Om M6 Tma",
				"D2ei 2n Mb zwm Sowy Gzwv 8m LN Notj 8lt-My 9y Gmta Moda Nm 92 Ryty So 9e Pmta Kn 92 Qmt 2 Tma",
				"D2U Udot  owr-y 9m Godq Loc N Nne Olwm Pmta Gmgj 7ndn n Mi Mndi No 9e Pmdi Ln Mm Potm Tmq",
				"H Ionintn-z 9u Pogm Gn Me Szs N Nogf-lwjzq Gmgq Sn 9y Pndf 7mdm Lo 9e Ootu Lm 9a Nodq Tma",
				"Tag 6m I Uqngi Nzdu Nn 9i Gmge Jns N Not 2 Rlt-Szgu Gzt 2 Oodf ne Nodzo 9e On 9m Qn Mq Om 9e",
				"Tag 6 U7vcngj-zt 2 Lnu Godr 8mc N Nmde Slwe Kmd 2 Gzdz 9n M3 mgf 7yt 2 Ro 9e Pmt 6 Sn 9q Lnty Tma",
				"Zj JB Nm Pnmdi Rnti Gzgm Pn LN Nm M2 Klt 6 Jm Mq Gy 9a Nz 9a Mmdv mwu No 9e Pm 96 Qm 9i Rndi Tma" };
		for (String input : strings) {
			String[] whitespaceSplit = input.split(" ");
			String whiteSpaceStripped = input.replaceAll(" *", "");
			System.out.println(input.length() + "|" + whitespaceSplit.length + "|" + whiteSpaceStripped.length() + "|" + input.replaceAll("[^ ]*", "").length());
		}
		String string1 = "66 59 23 b2 6c 1e 2b 3b 5e 26 87 6e 1a 7b 67 ff 2f 0d 37 3b 6b ee 5c 22 3a 7b 6e 1b 38 32 3a 7b 5e 3e 77 7a d3 23 1e 43 4f 5e 3e 77 66 2a 77 7e 2a 63 32 4e 66";
		String string2 = "2d 72 df f2 6b 7e 4b 2c 26 3a 7f 72 1b 2f 66 2f 37 0d 36 78 23 f6 5c 22 46 7b 72 1b 3b 73 ee 88 32 3a 6c 2e 2e 77 6a 3f 4f 5e 3e 77 66 2a 77 7e 2a 7f 6e 4e 66";
		String[] asArray1 = string1.split(" ");
		String[] asArray2 = string2.split(" ");
		int[] asBytes1 = new int[asArray1.length];
		int[] asBytes2 = new int[asArray2.length];
		for (int i = 0; i < asArray1.length; i++) {
			asBytes1[i] = Integer.parseInt(asArray1[i], 16);
		}
		for (int i = 0; i < asArray2.length; i++) {
			asBytes2[i] = Integer.parseInt(asArray2[i], 16);
		}
		int length = asBytes1.length;
		assertEquals(length, asBytes2.length);
		int[] xor = new int[length];
		for (int i = 0; i < length; i++) {
			xor[i] = asBytes1[i] ^ asBytes2[i];
			String hex = Integer.toHexString(xor[i]);
			hex = hex.length() == 1 ? "0" + hex : hex;
			System.out.print(hex + " ");
		}

	}
}
