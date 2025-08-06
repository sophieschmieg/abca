package fields.tests;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.TreeMap;

import org.junit.jupiter.api.Test;

class Aoc23Test {

	class Trie {
		private Map<String, Trie> children;
		private boolean leaf;
		private int leafValue;

		public Trie() {
			children = new TreeMap<>();
			leaf = false;
			leafValue = -1;
		}

		public void add(String search, int value) {
			if (search.length() == 0) {
				leaf = true;
				leafValue = value;
				return;
			}
			String first = search.substring(0, 1);
			String last = search.substring(1);
			if (!children.containsKey(first)) {
				children.put(first, new Trie());
			}
			children.get(first).add(last, value);
		}

		public Optional<Trie> lookup(String letter) {
			if (!children.containsKey(letter)) {
				return Optional.empty();
			}
			return Optional.of(children.get(letter));
		}

		public Optional<Integer> value() {
			if (!leaf) {
				return Optional.empty();
			}
			return Optional.of(leafValue);
		}
		
		@Override
		public String toString() {
			if (leaf) {
				return Integer.toString(leafValue);
			}
			return children.toString();
		}
	}

	@Test
	void test() {
		Trie root = new Trie();
		root.add("one", 1);
		root.add("two", 2);
		root.add("three", 3);
		root.add("four", 4);
		root.add("five", 5);
		root.add("six", 6);
		root.add("seven", 7);
		root.add("eight", 8);
		root.add("nine", 9);
		root.add("1", 1);
		root.add("2", 2);
		root.add("3", 3);
		root.add("4", 4);
		root.add("5", 5);
		root.add("6", 6);
		root.add("7", 7);
		root.add("8", 8);
		root.add("9", 9);
		int sum = 0;
		for (String word : words) {
			List<Trie> active = new ArrayList<>();
			int digit = 0;
			boolean first = true;
			for (int i = 0; i < word.length(); i++) {
				active.add(root);
				String letter = word.substring(i, i + 1);
				List<Trie> newActive = new ArrayList<>();
				for (Trie trie : active) {
					Optional<Trie> child = trie.lookup(letter);
					if (child.isPresent()) {
						Optional<Integer> value = child.get().value();
						if (value.isPresent()) {
							digit = value.get();
							if (first) {
								sum += 10 * digit;
								first = false;
							}
						} else {
							newActive.add(child.get());
						}
					}
				}
				active = newActive;
			}
			sum += digit;
		}
		System.out.println(sum);
	}

	private String[] words = { "9sixsevenz3", "seven1cvdvnhpgthfhfljmnq", "6tvxlgrsevenjvbxbfqrsk4seven", "9zml",
			"52sevenone", "41onevfsgvssxnpsix38four", "15ninedzhkpfstrscggbqhktwo", "rxbfsvhpnjvsixmxfhhmvdvg26rgrfj43",
			"gcbq2sghsv4fiveeightrlhchsfs2hsrjknfz", "tworgqpdjzrzf7one", "fivefive18",
			"sixfztrcxmbzktwofour3fiveeight", "2g4stjrjjmbngvljfvg24", "5jkcdxkltwo71", "vhpttjh2",
			"3threegmsppjrngfive7", "ktkgsvkthreevone2xxrxzgdqpnone2xnf", "llxzczt3seventwotwosix5",
			"eighteighthvllljmgg82eightseven2", "9hcjzphpktwo76xfpdvnhklzthreejrzkdknhrktwo",
			"fourtwojrvzzctzs5eight2vlm", "zoneightthree2hvhfsprqxmlsix7onevn4", "56sixssevenfour5twonine",
			"58fourfour", "683", "8fourninenhfpz9eightptsnnb", "gddmrzlpn9", "g3qsbqmbhqqp1eightkjggcxpmzgg",
			"6bzjqs7nxpvgtwoseven", "271lonepxp2flbmbz", "5sgshnrfn9qxt8xsnhtdtx6", "cnhflx4sevensixonethree",
			"595six93", "6eight8jbtmxdqpj96mqkrdxt9tpbpppl", "68qvmgth371lqcrglmvqxtwonfk", "mbqbhpmhspsbbxgflg98zn",
			"8oneqzvdcrh", "kb6", "mhgvpmfgjfourfourtwo4", "81mvcgqxlktbpkkrsgbdgeightbqn65", "hlbhthree1seven9",
			"vqczmtvqsnd3sh63", "qrzngvmk8", "chpldqtkhqbvdvmnqd5", "6three6two", "sevensix4tjvkcfgpone98",
			"96zgpcq9four", "cbvksjgvvklhnksixfive7kmmxg", "threeljfr1vzskggfvjsccc", "8pfourthreetwosix",
			"l2sbdcvngvmtzrtq", "hjdxfj94", "2hbhfbjfteightnr5sixfpfjkn", "2sevenjqlpprggjlkddqv9oneightpj",
			"fivebfmkdrrrdkclkscqmtwo6five7", "fourrp8fnbp9d1fl7", "six8pffqfpjl2nine", "91sdjlghq5",
			"fivemrmcnjmfcnck2", "v1", "78htfvqsztlsixeight", "5zpjtbgpkvkxbgpsp3cgklflkhdteightwortv",
			"mcxqfgxt49sevenone9", "rt2dsvpf", "9dzxmfour76six", "tsbmbdrgmzsjeightonezqhqb5qn", "seven8zllmz",
			"373onefourthree", "1fourfivezthdcxfr", "996seven2four", "8smnghninezczdlzxxgf", "7cbqfzrmhxdmrdr94tz",
			"6vstr65tfivelmhjshhj", "248fourlckvrtnzz4sxdqpgkvh", "zfivejfgfgdhfrhr6", "one2tworsmtnzfjhvvqjnp",
			"5threefourjsnfzninethreejrknzfn", "seven9nine", "gthree15fivemzpnqgttcq1", "fivevvvkgtvs3four4",
			"five121three", "7fourn", "sixdnineqvfqglmn98516", "ninesixrppx4", "threesevenxbskhlqbone276",
			"fivesixfour82", "tmdnvgrkjxfour9", "61nine", "vtxpkcgb48jzx8rgeight7", "twoonetwo3",
			"four167bnbdtxq4xvdgckkpb2", "xgbxvx7ninehmvqlldtxls88qhztfldr4", "5qxsfdchthree41",
			"fournb5gfqzfrlnc5fivethree1", "t6fourh3", "ninefgzksevenseven5sbcpnczgflqptjhk", "6three4",
			"5tnplthreeeightvsk2fivesix9", "three8hcvnmvbhbtwo", "6jnfourdbcgzgfzczbc", "z53",
			"fourxzhgjfrrbmkcheightfive7seven8oneightb", "2fivehkhnbtm", "rptgktsnzvsslcvfsevensixfive72", "95chcvtxv4",
			"ninenldnphcpn93fivetwo7", "threeeight675vmhvgtxxc2bp", "hfdmpv6", "xpffsljjcm27", "8sevensix8",
			"7hmszcvzpf9eightqqxgjdnhpfour68", "2czrnzdpxsjsdldcxq2vzgfrthsssix", "four5bdxdnvtgcbdzzsxjdlbckninetwo5",
			"35fourfive6", "z5eight3r84hvptgskcclbsgh", "2eightljldmconeeightgvhm", "4sixsixsixfourthree",
			"tnnthreelfmhone1onepfivevdg", "8nine9lmpptxstrrbqbjrjtm4seven", "2ltqd9mnrbhcsprz",
			"two5ptcjzbqnk2spvlcxqvglndgvrlm5fiveeight", "52sqstnqpdeight55", "2nvskvthreeqhkqvjzqrk1", "gbkqz127four",
			"8zmngr5381rb", "klrlqvkhslvbxdtwofour4cmlpbdm6", "eightninenjnmc8zrnrltgfmfour", "38sv7",
			"eightvtkkmhhninebmgdvconenine7sixjkdzxcr", "h14bdgvknbnjq4fv", "74svtwo", "7three5sixone4",
			"zseightwo6five", "shkznseven5ninefivefour8", "mcmztmrgctqszzpqs8tworfvkseven", "onejjmdfour15nine4",
			"7bxqxxxvlthree4tb9five", "twoqzpx6hj", "68hbdjdd45two", "68", "4hxlfbgvtdmbhnpfzxpnltfivemtmczjcbffour4",
			"fourseven4hfffmtqoneglgnnsrqsmbnxr", "dhbbmcxtv7five9ffjjcmshv31seven",
			"foursbgvhgbggnine6three4zznzcvmxxqktpdhff", "681cvdh", "3ninefivefive", "jg4dnnxzdtwoqdbvfsix",
			"g4dbeight6rkc", "236bjcppxmlnine", "5crthreefive8", "6twoknqfbmeight5c",
			"64eightrkvvtjtsonefour9eightwoqxq", "sm7sixfivexbrkvhtl9", "threeninetwotwo1one",
			"tjjeightwoj3gsdmdseventwo1", "63sevenk4blssxzkcq29ztjlnx", "ckpleightxgngnhkj2threem",
			"8knckxzrfbbpxldzninennnmkxkzzstsflv9two", "rmldrlcvccgdbvdnxvqr3ninefznpqrrtd1cfpjrvjzs",
			"5twojdflbxnlslnsjtklqkrbklvkthreejgxkfzeight", "sevenxfpzvss6", "one61", "2seven7eight",
			"fivefive5fourtwo9vjmxmtpv", "4pxfour", "twofive1eight9ldjbfbzhhffdb", "onefive7four93pjchnkrzvz",
			"4194hgfslhdczj1ztjtnccpjseven8", "sixtkjmtcthreezjfcgbdhvqlmbpzmc7",
			"nineeightvnd55sftvbf4sevenrkgbjskxbj", "kqoneighttpqcnsphsfour7eight9sixthree",
			"5foursvfsztnfivekkkftdfbptmrgcone", "sevensix8one5zjlfcrqzvp", "1hmvpbhzcbgbgvdkblnhgklltxx",
			"55xszsjgbd2fk", "hdrhmddj6ninebcfhbd", "f91six", "3tltgktntqpkeighthbqhlxcqnffgbfvxrkmr2", "3jthree",
			"2eightplclghthree", "771", "5gccqqmzmx81twob", "6pgnqbsqglk44mvslnxghckcbtffxvkq2", "twooneone5",
			"twohbxsxj2six826hlkdjnnz", "17tbdmntksvn89fivethreeseven", "sixplznhmbgzmlh3six1", "h2kvpbnkcxrssevenfive",
			"jhpmnksix7fgnbhgv8fztxfrpfvmnrvhzr", "vskpbnine3vjqdlmv9", "vmg86ht36sbnzphxfive", "klqxeightnine6",
			"4nfninelzshl26", "zqsb243hfdtvgx", "two9tzcrmrsix81xsgvldl", "ggqeightwo5bdpqtjkgzdclxsqptwo1eight1",
			"pvfzftdxqzfourmmbzlbzcfsevenfivesix5", "fivelt5fivegshvshpvkpqlddmfj6", "cjpzgsvgsz7kzcpr1sixfour",
			"8two4onesevenseven7v", "ztwonefouronezdvggfive8eightzflpvlnxfspgrtwo", "one7tphmqvfltwo", "49seven2three",
			"9rntqtxseveneightgflqzkrxhglrnvrvz1", "four6nine7", "5nine3", "4ddpkhft12hrhjvbtdxthree1tpknk", "7k",
			"fmlrpkxqktzdldrs84onefptgdjspldljvnine8", "oneseven4nineqbqlbxf2kqnpr", "ffivetwo97mp", "47h92ffcd6",
			"two7tcjvxdvctfivepcqxjnzfourfcndmdmnjseven", "cjdjxmnmpqr6cqvcscdpg", "onefour95", "seven9onenine",
			"trpvv6six64kfmxqrdvbsevenhqzk9", "tsgpcfxsgsgmhd7rvfz", "66788", "five26b", "chcbxone1",
			"6mrlfkqjssdxbtzkshvc8hjltjjxthkl", "four5sdsevenvztqmzkm21threenine", "snxrlctcztlbntnxkpmbs2",
			"652nnhvz6vcgx", "2zrxljdgsnnfour5fourp4", "85lzkjln55239", "51fourrnbkb4", "xbvsjmgsjdvhjpxdjhl1",
			"6sevencfhtbm4svgpzeight8three4", "827cclvxpdgqlhrjvrgxbxthree", "8pgqr", "three2ninesjztqjdhdv",
			"538one59one", "6lcmxone4ninehzzb13", "phckzkgtxdcnine8onefour744", "five61", "skbknb9qknrlszdt",
			"fivenine6jnxfsseven7", "lfj2onethree4198", "bb9", "kfljcb5zssbrlzml", "twohmspvflmxnine2eight8",
			"4onegfnxxqpqnfour7", "qvg7seven4", "85bheightghccknine1three", "fivedpthreeseven1fiveseven2",
			"6nine3gzlbtvtvlnine6eightmdn", "seven3seven58gzmdjmchdrzxgkgbbfrf4", "2lxbnnsixgdcv6hrqjnfhdmz",
			"8threefiveflbfmthreegtvckvpxkd856", "jbnkkbvfivexqzm238zxvztfl", "68kglf7", "9fivetzjsixkvfrngxbfbbjpd69",
			"ninejgz82fivenltntmfs6xpxc", "fourlseight3qrrrrddzlone", "seventwoqfkj2qxzddcgtb348", "xgjddt5",
			"q8bfhspkgmsevenninevdqmlzxznhmdlg", "qc2sqqlkfrkj652xpgzjskr3four", "two3tpl1cvmldrrghr",
			"qhoneoneseven9zfivevrkkjhtdf", "ddjczzpcvkksjzdcmxkhmbds2", "qnvgskzdb3nine77sdvfhfqsnv2kjffgsvz",
			"six8three3vbpnkb", "rgjmvgtgfour36qqbqznkjv78cbpdqb", "fournshtzbqfourthreefive8hsbpflngrvzdhone",
			"toneightone9four", "fourkkxtzpfivetwovnvxmtkeight4t9", "tbxlvkc5vgcmdckzv", "7threethree",
			"9ntqnzpldshfqlc2six", "8zs7five", "3threethreentqthree", "jmxjl4four5", "jnmrhzpdvbvvg9qcxjjmv2msrszndl",
			"onezjrkeight4two2seven", "32jlkhszgnkklbrsgpplphh", "1gsdqlbgt5eight1", "219", "93gcgx2twothree9xggt",
			"5fivenine1fiveplgmlffsvzbtqpb", "qv7lthszlgxeightnine", "7sevenmhmkcrkkq", "mxmmgqp4rptkbhfourvjh",
			"3mrmhgjmdv6pvfkbmconelmlckqkxjfiveznxg", "7one1kbjnmtpskgsix99", "bsvjlfgrcvvhmnfjzktdeight9914",
			"csfrjtwofivebcrcmggfpfsevenlmhncfb1two", "3five1vphfournine1kvfvzrtm", "nineeight3vxrdvttwo",
			"snncrseventwo5sevencjgl", "ngmgxnlsjjhvqpcxjshninesix63", "4threeone2twotzsseven4prqdrnjln", "fivegfmn391",
			"three5xflqlnrjgfpzt5", "28two", "six729twoeight7", "eight11", "sevenseven1fqcsevenonekrgxmone", "8rvlj62",
			"one3ninethree", "89s6two2twofive", "pnhhshxqb312", "vlnj5eight9", "5seven5eights3eight1",
			"6kqvkbjrfveighttwo6", "5twofour", "ztnkthreefive1fivetwovqjpx3three", "qs2", "2vft3", "zjr1zjpxndcsc",
			"76threesix", "1four88n7", "9mrrkjzlxdc", "1ll89", "fgvcseven81", "4789vtvtcseven", "2lcfbmqcqt",
			"5sixonefive", "cdoneight1onetwo96lpllgksff4hrzjkxng", "gntghdtwo1", "eightninefpdttgflvr2",
			"8zcnfm4krxhfive67seven6", "2hxpgvxgsmc", "trlvltteighteightrzkxntpdtpl24", "hkttzcmnck7s1seven",
			"sevensnmhgdxpbksngnflnthreemlqgdvphzk5tvmzjvdzbcseven", "eightpsbcshqcbppkgfxcnrgtwoeightfour6four",
			"53xktsrztnxninehpjjjktqnsixfivefive7", "sevenonebkdseven1seven", "11sixzdqbbppninehhkpxdbmlv2chddf",
			"3one3fivesevenlljjmxvzbcnqtszvzspsevenb", "1txlrnsb1vq28kpvv", "1two855", "2lxzhcjl1qfslvldkpdcxf",
			"9seven3fcqtzfive", "eight4phtznrb69xqbmxmdvxnsstnine", "eightfive2sixnhxcffkq95", "five28cppfive1two6",
			"p8", "kn8", "6hvzkkr4nine2seventwo4", "39134", "dztwone3kqlbbbknfive", "xvvdhddonekcgqqqzgxhlg2seven",
			"237rhppmlcmhsevenvnjxhzfnzbzrsdl", "1kgq1qsphhpcdeightjnfsdggnlnckgfbj2", "bnvpzxtnt16two",
			"mkhvcrfqdtwo3qfmhs", "1hkjncflcz5four7bbgpgcmnv94seven", "7seventmmhfgngfxt4", "98238",
			"3zmnhxqjqhjtptwoqtptmxfdp3seven", "bkqxlrtf4", "onefive8tvb6fiveone6", "1four3seventhree5mvsbsdjz",
			"jtjsjflv5gxvhdgvrbgbdcjjtchkvmf", "sixthreefcrbqzqs7psczkdc", "dc572twonejgl",
			"2xkrpz9dfslbjvmbdkgsixhkgcvdgktq", "five114", "8two85", "9fbh79mqbfsxcnn1two6", "6one9fivetwo",
			"gsjljnkhhvqlmmvcthreemfcbnjsbkvzzsnkb7", "nine49twojdqjsfbzsixrhbjhn4", "3threesix1", "52six1twoseven",
			"eight9dbb6hhjnt", "gsnqmninefivefbqcrlneight1", "jxczccqqpxbcq8", "1ckvkgqtvqrpvdrfivejjhhfkxvvhfm",
			"sixdxgvgglfh4qzczn8jpgqmgbzkmcdjfnhmh9", "eight7nine9", "tzqv6fivefourhfz",
			"sixnine2vkmdnkgtgnbbkcxvvsc2bf22", "78two41", "srpeightwovrhmbnkpnsix1",
			"fourninemfcvlstbmflzqf9lgvlvrlff2srxpzvrp2", "kcgtwone9eight", "7three1onesix1m3",
			"rgk2fiveeightthreegckdbd9dml", "sevencx9onenineeight6", "9ninetwocb4thttbkqj", "five26vrc2krtfivejpgdmtjg",
			"2ninethree7cnxbkpvthreefiveqclhc", "98seven28", "onetwosevensix9three8", "sevensix2",
			"ninesixlrdgpbrzs7onedrtlqpfour1", "pmvxzronexxxvbdrjr7", "15eight", "crksmfive64",
			"onetwo6fivemqkddjfxndjfpzmeight5xzk", "gq3ninetwo9", "8mgkvbpmbzpd7", "xfbvzlbpvb79", "kkbrvppqcg5",
			"94dkngltfzs98", "rchmjsrh7onejcknbl", "9threesixxrbzjt", "3963seven48sixeight", "98992", "7zqd",
			"sevenseven814htznfour1", "ninefivedfrtwoone5", "qzgrrng8six", "3742zsgbqgfvzlgsgfmxql86eight",
			"seven7zdfrhonemfhcfmclxj25three", "8pnhpnsqxh", "62seven", "eightpfgmmdg945ksctbstnh6cdxvgqbl",
			"932zsnvmcone1ktfqfmbnsfive", "dglzhqjthreemzpm78", "one7sk13ghnmvsrprg",
			"jdvzccvczspscxj5rzzdqdd44txvprhqx", "sevengmtcflgvpzonecvmbtgknine7kngpspbg",
			"xmqk6jmhmdtvh2kbchmsgpjrv2four", "tntsjnine3jksrrvone", "xfzgbzfive3ninekmnjrlqkzq96foursix", "3xbcth723",
			"zbvtsrxh94s", "threesix8mjchcpvmdgfive17xb", "8eightnbbzhfhf4fiveeightmnzqldztsnfour",
			"tkglxsb6one62rhzggrtgqxqnvjzfmb", "gtj2onempqp34qfjnlxtztvjvsv", "2eightzxdgc3", "32dgskjkh2",
			"sevendnvrcm5166", "5knxkkjh5", "fssjcvvckqjrcghcmlrkcvxl22", "2jbhdlfjtbfivetkjbjmtrgxrdmxsix",
			"fivedm8pzjfngzfsk35", "two31ldnvx", "5mfive9ggkq", "vqdhfsfkrpp44trqpnkqsdxvvdxc", "lkdqckshmn1rgrvtjqj",
			"ttwoneqzmsxzskbbnkfh8onespphhmsix31", "jknbnfvbf1snjscz5", "sevenvjn9", "ninetjzhbdjsffthreecfctlsfrz2",
			"6ninedqmxfqxssbhrrffpdhjvhqtkxfour9", "6tlrxcvhtllkrhjxqt3two3", "eighteightfive79",
			"5sixfive7bjlkninefour", "xrgrxsvm5bhcmzggbkrljnssxgpgdlv2four8", "1mxlnrtsjgdlcsncktrsfour",
			"8xtkgjjbtjfnc62vsshkjp", "91twoninetwo", "ninefour5nbnnzhtfiveggrjf7zqzblbml", "sixfour42rtbvlpcnv4",
			"oneczchdtwoglj9279", "jfvv5rdkpzldsxh", "ninexghzhdqk67", "gnllbntksevenseven7sixeightvgnfd",
			"qzqjxcqrs8fivesevenrnvnq9nqnchpjpmfkgtjqcsvtv", "6sixsixspzppcstlhqlssvt", "4hhfndc17rjrrzvbjj",
			"xm6sevenseven2fivendnn", "threesixd6nrxmxgcbjfrtlmpkjnoneshmrn", "ncnpmsixfour4two",
			"pqvz9fourkkmfvzbjqbfour", "szspfourhsqkfkfcndcnineone5khdb", "six29sixlkfbphnrzcjl", "cdtwonenine73",
			"mtlqpjcqcseven29", "8onekdvdmjbmbjtdngxhjjchdv7bvsbjqszlhvht", "fivevqnjsvnhvnine1kxcsjmzx",
			"6six32twonine", "nnsix5fcsrdvoneightcn", "11rzzpnrtnsevenmvnhgrsgngthree", "3fourjfonefoursevenbbfour",
			"fiveclfour7seventwoeight", "4zvzlfive7hkzhbqrleight8", "4eightthreemthreeq", "66nine", "dlnm59eight",
			"fourthreexvgpp7", "85sixfnrjqvmzmtbpxttrn", "one2seven", "63fourztwofdssrbjbcvhltg7",
			"fdstgbg1qdznxnvftfvfnr6djj", "5dcmbkrlvsrfvbdfqfour", "xp6fourrcfxdnktdctwofour", "9three4threeone6",
			"8pjzglttwofive9", "2gsix9sixthree", "ptwofive3t", "4hqkqzjkqddnnhxkrfnhgbkthreethree2",
			"threedthcktqkcthree22", "2gdlntwoseven527tzxbzkdjbv", "6hbxrgxzcnlbnz6", "ds7mvjbvfkn",
			"sixrpd7eightfour6", "qfoursspgghsflcrvqeight6", "94hxsj5", "bhhhnfnnkninexnjtjxrphrkc9mdmjp",
			"pbttrcplcsrldftsgk4991threepcbhxxfrgjddpz", "6four1fgjmjcnj3nbxxnnxhjhv", "6onefourshbzqgxjnhpmz",
			"9onedvscbrdj5", "jgjrpgcjztvkqseven8sznzl16dzrhmhnq7", "ctlbzmcctmkzpqtsdztbmllqnoneoneonenine28",
			"znkkbdsix23msxkcs", "seventwo365one2nlmlbgh1", "mtzthdjtonezfsixms5mxjpftkd",
			"onenine8ksxnslf16njqldnnkjx", "pzjl78nztcgnj211zxnmhxzrjjh", "twozkf3fiveshpt", "seventwo387bnfsix",
			"2xvvtzthreexjfxxg4seven", "8vgjxlvjc4twoglhmpgtsbfonethreeqbpcfbmz", "mhlgjlbm68rmeightnineone",
			"threezhldnzscqfour9seventfkj", "cnprmbrprlhmonesixlksvgzseventhreehdzh2", "smp215344two",
			"5nineonefive8sevenzvlfrprdq", "rvf583", "four8nine4crzpff", "rhnsqgnncnjzhpx297qmj18", "oneqlx716452",
			"mg3zptnqbm", "4ssgqpsevenznzlnnsbjpxkfvlkrjckxrhqsevenzrgjlkcffb",
			"six9fourshjkhspkcffiveeightqqkvktvjsmmnfrthree", "sixdljbtwoninetwo6ninesevenn", "kcjqfvldclrmsd52qn",
			"352one", "7nine7", "9fourjcveighth", "5twoz73", "bbsvvbjc1one5gtts3nine", "four2six", "7nine5",
			"nxjfthslhc2sspz", "69onethreefour8nq4eightwoqh", "8mfpmmthree6eightbksq7", "98bzjlfjm",
			"six541threepklngxmfthreevvkdgxpfour", "seventhree4ninesix", "7ninejkcxsixsixlvffhbkjkfive2rdnd",
			"9flxnpfddx68mfgnkkdqlonenblbszdvvbfour", "threemfk22vf", "threeninemh6", "88six", "4eightrtthree",
			"4rkjjgppc2lkvnhbnhffnrvqsj8", "65gdpxxz", "5fourfour6", "six9m9eightfour", "21szgzbcg", "krzfvtsgp3",
			"six2threethreebskteight4twoeight", "58bjrnzbqg", "8eight4four3", "four6qcrnthree2",
			"j8nine7gqlvtqxlzdtcbonetwo", "sixvlthree55seven37", "48sevengfgzczbg82", "eightsix2qsqfgmsrscfive",
			"zpbxnxbnmjthreebmksffpkd7seight8nine", "vfvhcvnsevenkbnhxxhhc3csk79", "147threesixtwo",
			"hldlvjninehnqxd7fiveeight7sixkkbmjxnjf", "six3gjbkbpgzdszl3f", "donethree3seven7xvfgnthree",
			"3lgpzhzmfftwoeightjhlngfn6five57", "5fnfl6vkncfsrrrsmxsl", "7ss", "five6gsncdxsixnhphl85",
			"3rklhfbkhfour884", "xsjseven7", "58mxlmsevennine7cssix", "1jnbbfourv", "r21kstwo",
			"zrqlghpcsixeightseventlnpnxpzv3", "887", "zeightwo3xzjvmqxpjtbvj", "1brtqrcvbrrqsqbpsl1sevenseven22",
			"4qcddseven", "mg3xxlgdkblgg1htdjtfthpv", "9eightwoslp", "ggqvtwonsfsevenfour73ccmzkz1", "7rfjfqh9nine",
			"gnpnrzzkdtfcseven45hvv", "34one1", "xzkfl9fivebxhcxrmfzszonekgzcvk1six",
			"8pzscngsixsrtjqqlk25ldczspzfdxrzhh", "zktjxcb72twosix", "9mvtk5jqhtwo", "two82klbfdf",
			"xstrdlsmx2onesevenrmgpjlrtfsfourfive4", "four1djgssbq8sixfive", "j1jkgzrjmnfxonejsvgzznsevenprbqbjn",
			"two65ffglkhgqb8xdxmldhtllj3", "lqqqsrjf27four9dgpjcgz", "fptvlnfzpfm99eight8five5xtgqt", "sgrdlpplhb2six6",
			"one4811eightrvkbchcngmone", "bcsv46dgxglfxsrn9two", "spshpxtwo6c", "1sevenrtbxc23",
			"sxppftvvfx4qvhbrrpcgfb", "rmcvzcl8dj4ninenpgjhlblbeight", "threeseven7lljggng5two9",
			"four47bnnine2chscshsrone", "kxtknjfourjr1pxxkqlxjsix9", "sixtwoseveneight8sixthree", "8djjznt4",
			"pj48hjlqfour7", "onedvpqsevensncd5fivenhftwo", "twonls1fivednlbmttzszjczdjfqcspxpt5", "11hbdpbbrjddrmlghs",
			"7seventhreeonemsszrjxd", "4smtgkzzrxlcslkhlnrft14fbtznvq5six", "rdv9four39", "four8fourvthreebrnmqmhb",
			"msgrzsdxfxfivecpjqxtphffm3zvgqxkrtwo", "1sevenseven3vlcmzlgseight3", "three1bkmpjqxsjmjsmlvhfrnine",
			"213sevenfkslf4ldg", "fivegflgdlqrfkmh81", "jrkqklq5sevensevenddlzbljvf", "5zjxsxqpthreeone", "fvh6",
			"154xqxjfvt", "32fnpjclfgjplmkf83", "6nine785six", "5one9fiventcxnine37", "24seveneightrp",
			"hbd6two1knqnfgcftq839six", "6eightzzqghd", "fiveptonesixnkslhvpkpsbfltfqnnqjfqjlhthree1",
			"ztdjmvrspdtbqsffive8onetwoblslmlxssjks", "876fivenhzmftccrcgvqnssixthreejhg", "jpvxbhbrv55ninefive1three",
			"ndpdvz2sevenn", "dkzpptsmbfhnnrqqtwodllstjtd87", "eight1fivejftckmgxmzbsq416one", "2q16ninesevenfour",
			"thdst7fourfivetwosixgfjsnvftone", "ztwone8tjlsrsqpnqsxfll9five39eight", "threefive757nzddvxh6three",
			"71two", "nbgsnbdvqkqr6nvsqhr", "four487eight", "mgqtnlqxpvbjrgqffourzbtvcxj4xht3twoone",
			"fourninefivefive97hzfxr", "4dhrgdrsr4onejzgs5x1", "cvnmzrdmtwotlvk4five", "7sixzghv5rxlbslsninefrsxht",
			"6lcmjnckbjrtwosevenvmtxcnfvkfqvgzp7", "dvndcdmtfxbtznqjrprksseven6twosevennine3",
			"threethreeone11ninesixeightone", "dvvcbsix96cxhgjpmqdsixone", "five4sevenmnzgktwo",
			"9rlsseventhree6xhtbxdnn", "bxlbjtmone15jsixeighteight", "1fivejlsevenmbbfksks",
			"eightsix35fourgcptzrfjhkrnfbbznine", "seven8sixfive72eightqjmcfjx", "seven8616xszttt3seven", "53twoone22",
			"onefour9ninetsrvczrnbsbrfstwogkkzs", "zpsplntwoeighttwo49three", "twoeightseventhree4",
			"eight5twofourgpnlsmppt1", "fourccknqs7dndlqhbsdhfmqgr1lvnpmxjtnkshgm", "five9ptrnrfdfgdkgxzlr6three",
			"pjxgrrgmfk57hdgclbftr", "gdcxzqldf3z2prsmfivenlklnrtbhfst", "999", "85lctfptljktwobtfnfrttlxrvlfvdnbsm",
			"172", "tmrcrhkvnfsixkvzhjxmngcrfmkfpzqcfivermnkxlfive7vtvvmmnfsz", "sixseven6two8kfrpjksixeight",
			"75pzmmlvsjn9987dftrvbf", "4z1fivethreefive", "foureightrdnnpxlnn32knstrmxg71", "svzljxbj113zpnpshmnf",
			"hj2hmvp3eight38gxngvkdmnzzzcxjkl", "gnrhzjzvrzsix33vqnnrmgtdvvkbsmglckd8",
			"onelgjfczsevenninehkhkxcskvcvnncbpj5", "onesevenninethree1", "dqqnkzfv9bnine", "8sevensjtppblkhh4seven",
			"ninekggzmfsfpbfnvtv3sevenczrhzztlsfour", "eight59hkthhk", "blskkshczzone7", "6one4pffxsgmc", "4lrbncnn8",
			"314ninefivesevennine6seven", "threeone3ccrfthz2seventwoctg3", "bfxknbtwortcfgnrcqsvqfcrxzmlmk83krb5",
			"fiveone8", "nine6661", "eight8884zbdzcsseventwonexgg", "16ninetnzftqzlpvgd2", "jxfccks1",
			"sixtfvgpjrvpr6pmsseveneightvbrxnq", "6eight5", "tdbs3cbfhglpfdnxlt", "nbcgkzchlj6pbcx4", "661",
			"jcjggnnn1", "eightdbtfgpfivecnlmnkrpgf1nine", "33bngjpkhgfqp4tl97one",
			"6rkmbdjztnfninecdlhnbnf6ninekmvxrqzbl", "465", "nfklpzppbq65threeeightsix9", "2gddbjlcdkx",
			"kclqnmpsixsix4fivepb", "stvlnrfdgcslqmveightbmbgmnzlrq8gfjndq6rsv", "sjzcfxn7hs",
			"oneeight32nine2ghx2nine", "2lhx9", "6eightseven9", "three9nineone7five", "1nrmk",
			"dccbdjqhvfoneeighthjsmfp89hheight9", "sixrfjshf9nhzngkgeight", "zsixone5fiveeightsix7",
			"vkthreesevenq42982", "2h", "8five9jkvqtwo9", "lkxvfm7qdhvnkt", "8jtttvnmxt", "1nine91sixsixfour",
			"8threenmffourone", "two2zjj39seventhree", "c4eighteight7hssvhvlm2six", "9eightfivezmnknpl5eight8seven7",
			"threenine4eight14vzmmhczfhxppqs", "seven22fivehgtttqthreeseven", "5nineeightthree",
			"hfffcgvnkrp5threevccpjmnfn3", "19seventwonelj", "4nine8", "2sixfggckcdt91three6", "9fpztvd", "pbm384",
			"rqxzzqtlsx8one91tjmqtcmkxhplcmns5three", "four13sixgreight64", "seven22four", "2six946", "nine1767three",
			"mdpvkhvbqstqpskhdxgbzt2zsdvsclhlzbcskckz", "1775", "pgdbfmvffninezhthree6qlrdkbvqthree", "29fourtwo",
			"fiveone14bsnrd", "loneight1fourpvcgxjsscssftbfxtkq", "fhtwobfnmvjxqzbzctxseven8lhxv", "ninefive67rqvgnbt",
			"2fivesevenfive4rhpvklfjz9ninezszc", "616vvbxfjplsppgpx", "rttqfddgone3rcvdljn88jqrlbdmxgv4",
			"8282eightseven", "mkhxlkksgsjrczffqmzzsevenv7seven6zv", "5hvkqdnpgtzfjbqrtzx5tncqbmxjpqmmzcf9",
			"llpkjvonesixlnf7one8oneseven", "five65qzjtwo", "xhcjxj1jghxktnmbxml5bbpklmdcthreedzt",
			"h3sevenonevdnjp5zpzfmch", "seveneight1mrrkcpbqd", "8one5rn", "qbffrljhl48qtg1jhngrrbsdhxl6",
			"ndc5sixcxlcgxpbstwoqfffive5", "fpfsqrzfjthreehzbcmhss4fivegbtwo8", "nine2two",
			"mddbqdmtcjrkqhv2dxfvdg8eight", "jbrkj2llgmg36twocvhmxnb", "mhdbsnine1", "four1mcvjdkmthmhcsz4",
			"59fbsnx7qrtclvrkfoursdpmhdz", "6d6four27zeight", "1sevenvjbqrd", "9fourlmjqn9rd",
			"gfxsrconexrgdzrhzcsh4six", "4one5jhsztrspthree", "four5eight2", "onebtwo4eightfourhkrsgeight",
			"51kjcqqxrjcnnine5", "kprdj25twotwovsdhzgmc", "htwone4344five1", "dvmkvcfcpsqrh1",
			"9twogzkc572sixhktmslseven", "fivefive7qnll7seven", "svbpx64n31onevzjhhhl", "8fbjkdcttwofourtwokj",
			"6twotwo5zkcnxczszfive32", "mbfcmsjmg9hmqngl", "ninesixsgnfzsmbgrlxbxjstkmmfxc5",
			"two8fourthree9sdxzvpgseventwonez", "td3two", "onespktrhrktzrcvdgqvdxgbgctdhjmm7shqcbzvfxhzlt",
			"5zfpfksszthtzxznxgkrpc", "8xvmsseven592ssmzjdmz", "64fkcmhmqdxnseven", "vdhkbktf1seven5", "tzfvfour3three",
			"fourl7four", "dkbbtpd5qbqgb", "qzdlttqfhn8chxxbnplt4", "qmd78hqdqxtx2rrdvkvvfourtwo",
			"91mmlbnbs5peightmznzhskfjv", "2rftqscv", "4oneonenineddktjvjlhone", "ninetwo5four",
			"1f36xndmtmmbpx1qzqmdkpbp3", "hxsevenjg6fiveeightwodps", "4mrndsix18", "r6klzqlz", "ndbrrsvp9",
			"sixsnkh1gvcnine5", "bsstxkninethree5ktwo7five", "fiveeight9five7", "3nine5sixkqlfrpdpcfive3",
			"lrrqkznlrcmbvdr6", "8xklphsevenonetworjgpjlrllgqcrxhlskfhpq", "fourmsthcgcxjsixcvnvninebdhttzm85",
			"8gnqnhptgkfivesix", "8onenine", "fivesixone3kpzvnbrjf", "9four54", "6twothreefthreefivetx",
			"4nxqzkkbgvthree7qxdhtpjv", "six5three2ninesevenfive", "3cqkmxnbkkh6tnszgzxqk", "82tzncrpvjts7",
			"39one5fivenine", "51jrdstpqnjdfbbtjz9three6", "nqrmg8", "7threesix", "one6fsxsflbnfivesixthree5",
			"bxbonethree55one", "7onextpttrflql6snmbdtbnvvfive", "threeklkjkvqzone2vhzsqdg", "eight231eightsix5",
			"cflpngxndfivefiveeight4rjrsfrmmtwonen", "three8sdone", "nine4nine", "two2dxzjxkbb2knvg", "bplttc53",
			"69hqtkfivesixtwopffgltlsj1rhhslz", "6fbg1scccrkjjsnnhpqmphksevendt5pcdl", "9onesixsevenfmmxtkdzone",
			"1onetwocrhcqhrxt7sevenr", "sixeightxjjqndfqtwo35", "9cqvtmfsqrfqhhbkjgbdk", "foursixfrtkpcbxgxx5one",
			"429three", "84mxhzmbdk8", "nine79", "xgntrzninemhxtqnine4ltvx", "four3seven6qhrzznzctwofour", "5h",
			"16eightninejmddlknrxfone", "3gjpbjdone", "654twocbczrzjnhkgdpqdd", "2xkptbxsixnprkhfslj",
			"19mlqcgbfpdonegdvzghjjb", "dnslxvmdlpmlsggq1one18", "lbxnz487vjlhhsxvcl",
			"sevengvgkk8mfbplfshlhqnrvbtwofour", "djfptvqgmkqgnzdvstwomzpcxfthvzpfsglc7", "1nlsztzzcbmqseven4fourdqq2",
			"cqcvkcthreeflhbcsbddg8", "sevensvtbtdkfkxzbfqznlh1", "sevenksn54", "4dzv6", "lsixninesevenrxn2seven2",
			"bsix2hqsvvvxvkpbg", "hplfzmghmbddz2htfkcfblqcdzfrvqpssbxdone", "19sqgxkn8four", "64ksmvcseven1748",
			"oneoneonedlmsdc1mjn", "7218", "3xdskntkmlcldqjgxbgx38rkbddntz13", "tnqnvshmhrkxbjvxcvdhmx1",
			"twocrzgfourbvtkcthdkrqpbsevenkfv8zczzszpf", "fdkdqfgbgnttlpnjrvnine6dvpdhtchfourbv",
			"qqrznptxjseven9twofourtwo", "hpdrqkonetvgfour5onepdk", "6eightsixninetclcsllxknspxfgxmlxqddvone",
			"twosixmcbbjthreekclp3", "kklgzxnk2eight", "zgmnine8oneseven", "gmcgzggsixvjzzgrs3gbzmxninezrlcfsphzhseven",
			"7sevensix5fivefourh", "rh1qzxvcmqjmtspknine", "12threecjltwozchdsfnkmchhgv",
			"vktfhngfb391ghtnrqfourfiveone", "9eight46xdxkqtkflqdv59four", "839sonesix2btrctxfm",
			"hzpjkvqdfg6three2twocsq8tskmdnvdl", "5f", "dgmvxqbpbjpbronefivehlf8ls9four", "sixm7m8three",
			"vgppvrgdlb26623csvkhsd1", "2hjfccgbjnhl8176xkpftwo", "seven12vgdnrvmmp1", "oneoneone7eighttworvpvsjzl",
			"738one99six", "sixrjm3", "15six44qndpslhnine8twonehkb", "912", "3blcn", "eight81fivexsbkzcthree",
			"kdqlzbnbnkh2mrpz82six3six", "sckfhxxjxfivejgtlmdhc3threeoneightrc", "znjgdjd8six6onesevenfour4qpnmvtdnnf",
			"7gnbonesixninehreightlmjone", "six3fhhlfgmdlgvhvqctcrxxh", "982vbjgptnc",
			"sgfnlppbvfzrmntwo9ggqzsixsixxgqvjvffour", "ninenine8", "kzfmgls5seveneight2", "8ninefourzrgjgrqkxrmjlzqb5",
			"eight7five", "fbjt4eight7rnhvfkl5knpvjhqdhvmvczxbvrx", "35five", "7gzjqslr13qsqxltsninetwoxmhgzhl",
			"threefour2three", "onethree45", "zvt8nvxctwo6", "sevenkfournine1drrmrmljsclgbgsd",
			"jfcnrxjjnbsrlblzpvxc84seven3six", "sevenonemt7eightseven6", "sevenrddndpj85fzzn4zhvthzp",
			"78zfdbmrfgeightgjtqnx", "ninerdqndgffive2dlsblldpfthree9", "lvzlfqzsixeightkqbnlv5njjsc4plh",
			"jsqhmbt1xvmkgfbghzdplkxdmgvcrkbngrjlpfj", "lfnvhdxcx7twoq", "five85", "1two3two22bdpbskrlph",
			"ninefc31fnhsnhf5", "2cxnrgtlvfvmvvmnfjllshmdvvfc48sevennjfk", "xbglffkvrzsmz1",
			"czlrrchbkhmz5qkdbtcjlffd5", "eight65sn4", "7seven2four99mlpskrgoneighthm",
			"9kfpfgzdjdgxjkltdkbkeightmxteightthree", "9bdsbeightjvkrmhdkghfive73four3",
			"xeightwoninehcrsdbnvtwovtbkhtxktjslsix3", "15fourlgrsk", "5xjqd9", "four8ttpzxpnrqnkz1",
			"fourvbfhg1rbrngbgfj6nineldqfxvrx", "sdpnkkfive9twodz23", "sixone4twoktcx",
			"8rjgbnxsixfivebsnthree2fivenmjxx", "fourddtxngtd4jvlttthhmz", "2nine95four1six9nine",
			"2gmxtrrkftjfnknknineqjqnscctfourzrqdrgs", "flghzhfgmn9tckbpmkgsix9", "jg9svtdrmlzm31rsrqvc4mggcj",
			"eightdpvfplptwofdgrkstvh8qseven87", "eightfour2fourvzksqhxmlkpkfktmdzpmthreetwonehv",
			"nine86kzqvkjqtjfourmpcggd8", "8nstjmtmstcnffnksqh", "bvgcmbcrgqfourpvs5xs", "8three12",
			"5398db9sixvnvcrztrqz", "7one62fourlndnshczz522", "qxrhp5eight183tfour", "fhpzgkt81two57", "ktlfdnbone6",
			"stsninecqxpfmdhk41vlpq", "eighttqcc5fqnfour84", "25gmh12threeltfnfdrxhh5", "57four", "mqgdhfour67",
			"37ninetxkddhfive", "rzrsskzrlzjbcgthreeghbqhdpxfvgjfqclcf4", "fourvone2vbpltlrj", "xz5four3nineseven",
			"1szrhcmzkftwo9eight2ltjmgjzcblzone", "zlnkddtgsb1sixsxvjxgxp2", "26sixpzpsixtwozqff",
			"seven99fzqxfmttfgxm", "9twonineonefourpttbgkxt8two", "fv9", "5qcmjsfk6zxjld1", "fkjstnvmchsr9q699",
			"nine78three", "4rcs6bhbbgzhsstwomnineksbxfzj8", "4fmblhqninefive6qbkm", "zsgjbfrjfour1sp3",
			"zbfeightfive1oneonernfd", "5bxtfvzczbhtzfourqglqdxsc", "f9five7five8ddvseven", "23bszpdxfjmzg",
			"fivegctmd3vlcgfgnine", "63hbdkxljlq", "64eight6eight6gxdpmtnbfone", "28xcbtt1", "1six5", "four289" };
}
