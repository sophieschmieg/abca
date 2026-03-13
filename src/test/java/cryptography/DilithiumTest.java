package cryptography;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.Charset;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

import cryptography.Dilithium.DilithiumPrivateKey;
import cryptography.Dilithium.DilithiumPublicKey;
import cryptography.Dilithium.DilithiumSignature;
import cryptography.interfaces.ExtendedOutputFunction;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.helper.NTTRing;
import fields.helper.NTTRing.NTT;
import fields.integers.Integers;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;

class DilithiumTest {

	@Test
	void test() {
		Dilithium d3 = Dilithium.DILITHIUM3;
		DilithiumPrivateKey privateKey = d3.createPrivateKey();
		DilithiumPublicKey publicKey = d3.createPublicKey(privateKey);
		for (int tc = 0; tc < 10; tc++) {
			byte[] message = new byte[100];
			new SecureRandom().nextBytes(message);
			DilithiumSignature signature = d3.sign(message, privateKey);
			assertTrue(d3.verify(message, signature, publicKey));
			assertFalse(d3.verify(new byte[100], signature, publicKey));
		}
	}

	@Test
	void testVectors() {
		Dilithium d3 = Dilithium.DILITHIUM3;
		DilithiumPrivateKey privateKey = d3.createPrivateKeyFromSeed(
				ByteArray.fromString("2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a").array());
		DilithiumPublicKey publicKey = d3.createPublicKey(privateKey);
		ByteArray expected = ByteArray.fromString(
				"f5408337d0fee65c28851226a5fa81b58464632c78e2a9bef70d330f2e3a5f74d9cf676aedd1067c91a5dd5d4edc46f868a93ffec9f44e254e44f682a153aeadf228e8db7c5fcfed30cc3408e261ab896876bee56660d2a7c1d7eac20c5754255206a178f7156295065ce7876f90c48f44bc37f3a00e32eefd3a4bb1e298fe283d106eaef92a33a594253a2a0790976a1d04636f8672d28c06c852ea8bb43b84bff512996e7616963d5b9a2906466a152c7ea9be178be35405683b44367af85d2daad87630c1e21ba5490154f0141780f5ed0407cb0b975dd56d5930f9b26413b843b83f3693304b0038bd3e4bb398868060ea18c9c67099376470a50deb052e4056743fbcdf0341b192663bd1c21ba3b3d5666e0d0e29c4e1ed0759ab0bd9d1d355011b94e0ff0c049b03ddb7138640667144fcacd7265f55a07e5387f1abd30c037cf14d436aa855f827049215440d8007f61460500d943f57ffb6bfee6fedd2fcec52882d7d8da1aab29e892c8beac3df3234b4a7d2eca3a45c6623c52bbdd07c1c94314b706988a52029f8f8b06e874b741d72926652c78c6ace2cfd8864eadb2e4b39cafe6e03e4edbafa2747db9bc42f92af8b031e3e380846b1bfd15ade88c285d6a6fffe91eafc8b17de6cbc68575f323cc09fc20e49e8efd76f9568bec486b78df4245428d8d0d5f53873e11de65fda4c770b521a8c67f5c51d48cc26358954514447881fd9a42e5891dac7e1db5249d7861b322111e5fb929bee9ff5e9d5a2667ba93e63fc03040d2e82648f89e89dec1d1d2dfb9efeceb7940f7dcbebeb5a239cc1c54d8f7d52cba220d0634e15df46a58280bc5a48840bd39274cfde150f9ad9a40f6398d715350925f0e0501944409f32331a362bdaaafb3d8ce71c964332d6afb7e684f99951246d88081c86744ae68133f22c53a4b5ae258f230a98491d2d43a79a6d0f4d54a3b62013965ac7c82d0507125a38a0277f81cbc1d46cef2a131c6f51b88ec0baae0c82a6a0e72831cb06f9116cff5111d597e01057d32805a008f52c9aec3311139bfb35982789ff83bdd0c31e9f1080e8ed8eb99fde66bafb29e3357389fe3785b60c78e229ef073e1b65e34d848bd4d8a4f251551e2d38d2546afbc205d3c6dab34d2b962b1afb44f1d22fc10c6744fcd6b636afd3cb414b16c2e0d708fe9f51ff19120bde693b028b6d1e6dbe37b4b8b3bc7c6f7a842701603869d3ded572500f085502efc8d3cc62b30e5cdbcb5e86d9c0d42973bf755df539cc0aea58f9148386db67bd2bf70cd12ccd96d5c66fb271416b772465228dc44b079178f9b766370b66a79b871faca246ca6f8f63be9f0668297ac446cad5cf4a83318b1b00ecbd283f0eecee60a9a37a27abdbdbe382e307970002837dfc0bd3934ebd008918fd4bd383c02c9d37f694996e989a49075767ebc4a2981ef5275455e026cb0bd70946cdd1fadaf251381d324f9efbb860d1b280c29685bab97d010676273b45cca12ac3966aae342c84e2357eccf252577743b8787967b40b07ef2d3d9e6c1a3bcb059cba0fdb7f0d4f815c242b8e14acd3375e608e9230ba3cf8718f43882a3e1e661a2bbe81830d34741f33473e263b3790abe67acf29f5df44865b2ffbc96975fd62738a64112deda5a2534fb0a23b3b3024df986391badf9041c593c313a7ca1e1fcffcb65b07b9a99337b4a4acf616cbe1553eb9541f38aa6247342905995233a28172ca13396b2a9662970120f82b92a213f43de7a232ccca3268265c9ce042d50915430a6c455f32277da42f9962fb9163b623231ebc080fa7b8e9f9021fcf85b98f9c483e4d2226b9326a5bcb2e7449ef029ae142d3a0f0c28bd4f7e9c51a12e1336f24dfacbc3f808a8f7dd683027bc948763b808fb0037394b8b41bc9b2ec7887e67584e03d11b15ca203b2bcb43f8881638c4e4eee7f846d09c7f89b7739df22b2c3acc235032ba8f7ae27b5b9d25733143e80a4cdde6770719c1e66ec2ce683612233e88fafff84c0745a98aa1254c8219c6c556348c2b5d1beeb61532d6bf7bde153271dc647460beb65fe0055b33fd6480dcbb9d7d471952cfa5be260c39721a8c5c89b9e966ae2dc9036451ec9f2c49433b2225e13f23e20c2bfba81a7b3a555883449238f7d48213e9f10ce19e76f1bdcfc73ee5524bd7d8be0a4b46784e238233c04fb99383ec7726f9717e1179dd14fba9ad6c2ebd1699f0ab0e57e6cad23875b029e89cfda06f51266ecd2eed4edafb51e82f2a506d57ba74da611774ca5fa2fff4a976519de425885e7d09219cf815b1767d4fc5a72c18918991a285086a6a766614a4d245387da50f28dd778fb33ab88c0918feba3768c55bb1f07aec33cfeed33d6faa4d34fd7227b365533c1e67dbc89f0b20195cf1cbd480d333ade1c9bb28308085b72ced430268c1492a27050c43668adc9cf8b8509447cfcd3c8f8d8eb554f704101786aa9ebca86991d250776a37a1f56fbf7d08e591f978da49c3870625879f70e2418aec5cba32fa8c346fa9038baebc35ad0068a4d03537aee14c2e71570a87490377fa8dd66f995aa044a522f0c7025a7ab2dd5ad30a64268dc112b7f9fa156df64d631f55f1d6edc55cec570a9c7372e29e02c8d4867bae249431dcf6ed2794a0183f0f7501201feca4a81d334c642fc8d38e9a90fa77429665e09e214797dfa455ff47c4f219d3a2cb0176bc2236455123c1c5da714ad29d580fb194f87173a18dc");
		ByteArray actual = new ByteArray(d3.serializePublicKey(publicKey));
		assertEquals(expected, actual);
		ByteArray msg = ByteArray.fromString("0700000000000000000000000000000000000000000000000000000000000000");
		DilithiumSignature signature = d3.sign(msg.array(), privateKey);
		ByteArray expectedSignature = ByteArray.fromString(
				"424ac14fd2cf3f7eb28b104a7418909e84eff42d2f6e7c4ed3874838f8580d918a5d229f8ad72a0e093959da15f05e9942fad12ac1fe73c9de377b85a0f9d7203e6f8d3d4e83e5ac52a82b7a77c0096d92d40f2c1b7e194673faeaf16d7bf11342956d14d7bcd45f98ac9003ec86da5a1b2775d0723ba7350b81ecabc464be480d834b3965fd6554d340b9a11add210c6228f29a8932b8b6a00263777039533e1580faca62a7885372fe15b4e441a7595a13966fc05a95926fec9eed7563b1347d2b8d55e15212148e9f7e6e6b3e8dc30c6a1d213f66f574588ddf5d559438b1bfe01e851b15108baa97b64273d2e850b8ad504d2defb158e1ef3d4af5e1cbebfdf8a0e8037ba4def90bfdde7f2eadb3da96becd9e6a95e286e120871217aca6a8756f3ed7a11c2f75e1a82f9ab0234b1cc2c8178ad690e8c2798f1eb8b83178f39d7fdd09f44d7df268e7e062566a4516b1a1ceb97049df4646e8156e678c5f1fe2705b0c261db71a9932b8f41950e326a2200908c4c0bf3681384ae1501f5ddef8ab98133519452c725eb84f36814577108508809015af4d8758c23afea3815a09537477e8f6c7694a81c1d08c822843b9e67d95592eff819f6c2927e403f9962f36364ec2dbfe494049884d61695331e32abe6272d5a44961b8679879a67b142aa40729630d8d24bc69caf17774c6e2b3021b9181a8a498586a0f7d4109c5c82841454785873670dd2f60a4c278a373986333120c1b17b4297620bf72c7fc126205c1b240c6d60399b715579adfd96492e350a4145ede8d075ae544d051f58bb9991575a6ddab224c9c83a39b6fa9b49444fb0d2b8e32464c216e42d02091ff25122f85f9454f6d818034a30e6252cb0d25b572de3ab1c669709e04fbdaeffa26d0fc67d250da0864c71d74f8c08b5ad15cb95a33ca5b03d9573e83b35cc173afd97e88b9ea859861384cf8ee6338bbe5768cdf9c84cabdb8b332aaa7574ba200d496aa1d6be562f59be022d7d9cc63fe6a35100b664694b5bd1f5a610cc3d7a670f45b47befb8cf705339a8a5adec9fbb9f1fc4196fd3ad2a5ec9f3c9bd48e96bd6ac7276bb00705aee192106388f2fef2cd3c7cf00a159459314260dc49d9630d9bdfd0245db95fff41703f3236ed1bce77582dd6de015ed6d8dec0f2f2952fafc0e11479487a21c1cd61f1a9410dc730dfdf01527c6ae4072eb1591fb5f9b6117dabd9458331d9fdf582acde4b9001c6b6c6e048120cf15534db40e7b8f77fdff60b08a0ed6a802095fb76c71062dc34a612fe34f6598529c2028487457530c154d18709381cdb36b4ec697a2b921ec9056c86c5e73dd74b63f1844b1240c703c6f6373018c2daebb7beb55ec0365c76725063a3fd212ba80bc951baee1b58ceab532e58754ea44b0cf815a0624e2a0cf5dea9bfcae56b7ac74e4627a05bd8a2539ea4d25a3aceb003c8e7b750f0e8485eebc84b130a54ce8926b54cbb54b94db057f993ab426fe1a3be7a997782fab08cd704b1628b9d40af916730014d71d71d82bc09a7a17855c314b9febd079e371598103a4fc6923fa47a1979fea54c0bf7542a6e54d5cb92814c824d94753095d2fab232153f204d16aaa6acba4c72ae358eab79a8b51ef06db93cd5acb606ddf4dfbbc031f697486c284ef2627e3e65a5998f746d6ed934cedc1204c72694f6738ef0de7cf25bd147c938f3e7efd6def80c01b690ed200c21f56df02f7f31e90c0f05e61c60e2279b307664f491fac8f3315e48040ed63f6b68cd7821bd9447a8cd9d95aca3b16d721ef8e90d668636ebb4974553f41e7e9b075eee4e83e118fa74f2d7197fa135dc53146e7c194645e4f07412693320117cc3cdb1bca814e0e30fcf5590ce237e0479bf088cb0f68c2dc0b53203c7c007a4bc6361a32e3e153cbef1163eba0a3d7e5e133fb30c6f37b3f92e94abe270e7c3e8c97056bd12a1d588f919d65db51accf624e0a356d6e6e0da3c67dc298746a39b6313bb658e3451a2b1abf7efb6c1e63b83cfaeea9a17ac8f22587254227f5adb2f41f367bb9e4530c7fb93b78cbf60b7619a25f2f7a00aff9d935c294e82102eb0e96da7c3d4c87ccb896a5e5ee3ee13da3d94655a66c0c3d781a67c23bbd53cd7eea563a8d63100144ede22e3af96b69ca18d5d5be5849d27c1db28298d485e2fba3c7f68297c7f883049d73ffbb3a300c5aacc1d8a805205602e410e76fc48adb109cd47624c2822e6dcc0db0a28be5401a79659d19b8bbdf131ce04fd6278130aa2f203b9e0f52134a3dec5a97c4753da3127288ea0969fcee232837ea067e8a34b6dd4e1ecd7f37b8e956d137be78be955854a75d84a339f3e80bebe8d295ec1149bbaf115bfc84e8cce1cf994fe067dbcca2ce1c378180bbceda95444ae535a9ffe1819caac5cf95b86219e95c62581d1110daefeb4830043a7a9552de5a55f250e7f1aa5993f87647c92c27256753d6daf32da18c1cef9bd5ff5a7a8e182d071f1e75d10be89f02c861a0d05257c08c9e12f877a22c8382ec883c592a34722f1138a1e8b5ded41e970b97ba2a69c386f94b508b38a59b464dee983741e61d0096f499f3c552f43ca317ddf073e0450e800585053830681d4a8a8ecd5b86c901285c84278c06b86375a04fa7002b259cd61d5ecca426aa334659284945ea6c5adedb7d40835c6134271d641bdb19cc3617bfdd53618252b4a684a259b82982b0a13cc90ca931c4e90915d2d7a0abdee3f78764b41e3236e2ac21e138c0cbce20740cd17f34eea236720acff25161ecd4c903e33395f51d2c51583407b35a044d54caee2b66dcaba84359218c6803040b19a39c0e4ef570fad7330041e024cd194326891547985b10e505c6e6b1abdaf3fcd63c3f2f08b2d93aad41e83167d969e54b12745b328b5cfb93266afa793bf00c86270ab98e2ab44b2be16b5e522ca280cfb0a9a4f744900f0f44f6126b787274b2a36232c84b563dc3a5879c83c34cc6d644fedc3912b10e9806f944440b00883bb76d51782f928139117c7358cd38a460ff73af4d417278b7bd463da86ba11bc88cff14340a24f675905c7edacb121fd625187b90a3152115b3ce1f49432409fa629384509cf51426948c71e76da8efbde2eeb7d3f3e576f18d4fb61f90d99de2bfea2cec28468368610670890ef26e05bb108dde48095251574c20611ce68890144750c62a9792815a89ffc50ff0d76fa0f708d917f3bd9f840b11085445784d2017e54266c69f3fbef4c8fa82c0723ba58736aa8466a81b00d563dfa074a9a06016ffb04f029bc12c53585183f237d84c343a13b08b484e1625be1b0c661f697ad5097040fed2d3d46c6b1229157729b306694bcdfeb4eece5b5dfa8043a6fcd607e2f82aeda46e88f8d59f2032fa439365bfda8ac7b9e09ed8020c045757d510d5de5dd0b2520b0f74a99d8fe7b1915d82491a270d1a1c554d0c05d14a1419e71dbbe5df225fea5c49b73b84fbf079dbef25ae8861e58ab8bab6ec4ad07dc10072c5627029d4a74700193d12f439a5bdc41a284b92499b3ab0f3973529ea148053da54125c2969603fa55d3dc30170568ac9642aedcdc817b0f4b8d5d31925428afbedd21fc7f1132fc7e85284723b4f26ba1b50bd5fea2d81b12f4451f97331c6af0c82aa9d962ce034df061c32a1113990a9f320291993c9d614c42e2aadb5dc52c58daf8d22bb302e1e64b32e16394cc0c44cf93bcf811e7ad772754b6e9c2dda630449e005a0303996a69f5ad647f98908af4f857af804294d416db784c34d2d5aad7d840e135ecf6b636dd33f4ba29d0f12fe3967260013662d8bc23a34a0f8a5df93e111e8792140487eb04854bfd9bf679e1b38de028ea5fc457cbc8e2706fced55f5363425b0c2c13cbad1ef96d18f11b0fb3749a561504b8fee8cfb4ef9c948e3efaa2de6ca38fb445382af66855f3c0603ab175aca17e1b005deec73b03eb952cefcb3cd5860c40578e8a366da62a431916ee212d9e5740fbfc42f6d0e0951fd2f213a8f6f21984e01832c0dadf56b3de4833b849b7d0c9125d7b9ae5c6f245703d77c0298c9b6abd563529b3f2986ce54462f050d18dfc85a840b1633b2b0373e585f77c106d70ff711d03362af483f7812965405e586d6988318effc81e9fa918ab278e97e0d6d3932e24ac17ba537f568f8b16bdf3d45488cb34f696295f7321b63a76101ec283f0427a533544f752cc9d3c277c59997e01604192746353421884142fa23dcccb3acdac509a5688089ee919230a9e0ec2804254a153f544977a634bb13880f3241daab08d1b73bdd8fe4a187e9b4c3758eba88cdcd9fc104360b2512788f4d31a23c80b78c150af10e1a41dd279e895ae9e66ecaf2bf32786c11573a74cbbb93b82255f9ed190dafce6803f07271889828696f5f8173772a7142ab46899c609136713cf7a3f4c5e9718c84ee8bce8b4b5ba83a09e421c830eafda214c0e0b3dc28d0f1adde9b5fd9fee6a07a328807d0e983cef6383ab94c415dd1e649d9804f571581ae0c7ab2896e86a470a9eab94f2fe9e138e91ac103cd770d27af9648ef83b5254696d72d7174188bafb05114a788d9726c9d7eafb182348618cacbec6224186909c00000000000000000000000000000000000000070c12171f24");
		ByteArray actualSignature = new ByteArray(d3.serializeSignature(signature));
		System.out.println(expectedSignature);
		System.out.println(actualSignature);
		assertEquals(expectedSignature, actualSignature);
	}

	@Test
	void createTestVectors() throws IOException {
		Dilithium[] dilithium = new Dilithium[] { Dilithium.DILITHIUM2, Dilithium.DILITHIUM3, Dilithium.DILITHIUM5 };
		ExtendedOutputFunction xof = Sha3.SHAKE_128;
		InputStream stream = xof.stream("This is a seed".getBytes(Charset.defaultCharset()));
		for (Dilithium d : dilithium) {
			DilithiumPrivateKey privateKey = d.createPrivateKeyFromSeed(
					ByteArray.fromString("2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a2a").array());
			DilithiumPublicKey publicKey = d.createPublicKey(privateKey);
			PrintStream out = new PrintStream("ML-DSA" + d.getRows() + "" + d.getCols() + "-violations");
			out.println("ML-DSA" + d.getRows() + "" + d.getCols());
			out.println("Public key");
			out.println(new ByteArray(d.serializePublicKey(publicKey)));
			byte[] msg = "If you can read this you know how to hex decode".getBytes(Charset.defaultCharset());
			out.println("Message");
			out.println(new ByteArray(msg));
			out.println();
			for (int i = 0; i < d.getCols(); i++) {
				for (int j = 0; j < 10; j++) {
					int index = stream.read();
					boolean positive = (stream.read() & 0x01) == 0;
					DilithiumSignature signature = d.signViolateZ(msg, privateKey, i * d.getDegree() + index, positive);
					assertTrue(d.verifyNoCheck(msg, signature, publicKey));
					assertFalse(d.verify(msg, signature, publicKey));
					out.println("Violate index " + (i * d.getDegree() +  index) + (positive ? " in the positive direction" : " in the negative direction"));
					out.println(new ByteArray(d.serializeSignature(signature)));
				}
			}
			out.close();
		}
	}

	private class XorKeyInfo {
		private int[][] xorS;
		private int[][] xorE;
		private DilithiumPublicKey publicKey1;
		private Vector<NTT<PFE>> uncompressedPublicKey1;
		private DilithiumPublicKey publicKey2;
		private Vector<NTT<PFE>> uncompressedPublicKey2;

		public XorKeyInfo(Dilithium dilithium, DilithiumPrivateKey key1, DilithiumPrivateKey key2) {
			Integers z = Integers.z();
			NTTRing<PFE, PrimeField> ntt = dilithium.getNTTRing();
			this.publicKey1 = dilithium.createPublicKey(key1);
			this.uncompressedPublicKey1 = dilithium.getFreeModule().add(
					dilithium.getMatrixModule().multiply(dilithium.expandMatrix(publicKey1.getRho()), key1.getS()),
					key1.getE());
			this.publicKey2 = dilithium.createPublicKey(key2);
			this.uncompressedPublicKey2 = dilithium.getFreeModule().add(
					dilithium.getMatrixModule().multiply(dilithium.expandMatrix(publicKey2.getRho()), key2.getS()),
					key2.getE());
			int eta = dilithium.getEta();
			this.xorS = new int[dilithium.getCols()][dilithium.getDegree()];
			for (int i = 0; i < dilithium.getCols(); i++) {
				UnivariatePolynomial<PFE> coeff1 = ntt.asPolynomial(key1.getS().get(i + 1));
				UnivariatePolynomial<PFE> coeff2 = ntt.asPolynomial(key2.getS().get(i + 1));
				for (int j = 0; j < dilithium.getDegree(); j++) {
					this.xorS[i][j] = (eta
							- z.centeredLift(coeff1.univariateCoefficient(j), Dilithium.PRIME).intValueExact())
							^ (eta - z.centeredLift(coeff2.univariateCoefficient(j), Dilithium.PRIME).intValueExact());
				}
			}
			this.xorE = new int[dilithium.getRows()][dilithium.getDegree()];
			for (int i = 0; i < dilithium.getRows(); i++) {
				UnivariatePolynomial<PFE> coeff1 = ntt.asPolynomial(key1.getE().get(i + 1));
				UnivariatePolynomial<PFE> coeff2 = ntt.asPolynomial(key2.getE().get(i + 1));
				for (int j = 0; j < dilithium.getDegree(); j++) {
					this.xorE[i][j] = (eta
							- z.centeredLift(coeff1.univariateCoefficient(j), Dilithium.PRIME).intValueExact())
							^ (eta - z.centeredLift(coeff2.univariateCoefficient(j), Dilithium.PRIME).intValueExact());
				}
			}
		}
	}

	private Matrix<PFE> toIntMatrix(Dilithium dilithium, Matrix<NTT<PFE>> matrix) {
		PrimeField fp = PrimeField.getPrimeField(Dilithium.PRIME);
		UnivariatePolynomialRing<PFE> polynomials = fp.getUnivariatePolynomialRing();
		NTTRing<PFE, PrimeField> ntt = dilithium.getNTTRing();
		MatrixModule<NTT<PFE>> mm = matrix.getModule(ntt);
		FreeModule<NTT<PFE>> free = mm.domain();
		List<Vector<PFE>> result = new ArrayList<>();
		for (int i = 0; i < matrix.columns(); i++) {
			for (int j = 0; j < dilithium.getDegree(); j++) {
				Vector<NTT<PFE>> basisVector = free.scalarMultiply(ntt.fromPolynomial(polynomials.getVarPower(j)),
						free.getUnitVector(i + 1));
				Vector<NTT<PFE>> mapped = mm.multiply(matrix, basisVector);
				List<PFE> resultVector = new ArrayList<>();
				for (int k = 0; k < matrix.rows(); k++) {
					UnivariatePolynomial<PFE> mappedPolynomial = ntt.asPolynomial(mapped.get(k + 1));
					for (int l = 0; l < dilithium.getDegree(); l++) {
						resultVector.add(mappedPolynomial.univariateCoefficient(l));
					}
				}
				result.add(new Vector<>(resultVector));
			}
		}
		return Matrix.fromColumns(result);
	}

	private Vector<PFE> toIntVector(Dilithium dilithium, Vector<NTT<PFE>> t) {
		List<PFE> result = new ArrayList<>();
		for (NTT<PFE> coeff : t.asList()) {
			UnivariatePolynomial<PFE> poly = dilithium.getNTTRing().asPolynomial(coeff);
			for (int i = 0; i < dilithium.getDegree(); i++) {
				result.add(poly.univariateCoefficient(i));
			}
		}
		return new Vector<>(result);
	}

	private DilithiumPrivateKey xorKeyRecovery(Dilithium dilithium, List<XorKeyInfo> xors) {
		PrimeField fp = PrimeField.getPrimeField(Dilithium.PRIME);
		int bits = (int) Math.ceil(Math.log(2 * dilithium.getEta() + 1) / Math.log(2));
		List<List<Matrix<PFE>>> adjustedBlocks = new ArrayList<>();
		List<List<Matrix<PFE>>> blocks = new ArrayList<>();
		List<PFE> rhs = new ArrayList<>();
		blocks.add(new ArrayList<>());
		blocks.get(0).add(toIntMatrix(dilithium, dilithium.expandMatrix(xors.get(0).publicKey1.getRho())));
		blocks.get(0).add(toIntMatrix(dilithium, dilithium.getFreeModule().matrixAlgebra().one()));
		rhs.addAll(toIntVector(dilithium, xors.get(0).uncompressedPublicKey1).asList());
		List<Matrix<PFE>> diagonalBlocks = new ArrayList<>();
		MatrixModule<PFE> mm = new MatrixModule<>(fp, dilithium.getRows() * dilithium.getDegree() * (xors.size() + 1),
				2 * dilithium.getCols() * dilithium.getDegree() * bits);
		List<PFE> constant = new ArrayList<>();
		List<PFE> rhsAdjustment = new ArrayList<>();
		for (int i = 0; i < dilithium.getRows() + dilithium.getCols(); i++) {
			for (int j = 0; j < dilithium.getDegree(); j++) {
				List<PFE> row = new ArrayList<>();
				constant.add(fp.getElement(dilithium.getEta()));
				for (int k = 0; k < bits; k++) {
					row.add(fp.getElement(-(1 << k)));
				}
				diagonalBlocks.add(Matrix.fromRows(Collections.singletonList(new Vector<>(row))));
			}
		}
		Matrix<PFE> blockMatrix = Matrix.fromBlockMatrix(blocks);
		Matrix<PFE> adjustment = Matrix.fromBlockDiagonalMatrix(diagonalBlocks, fp.zero());
		adjustedBlocks.add(Collections.singletonList(mm.multiply(blockMatrix, adjustment)));
		rhsAdjustment.addAll(mm.multiply(blockMatrix, new Vector<>(constant)).asList());
		for (XorKeyInfo info : xors) {
			blocks.clear();
			diagonalBlocks.clear();
			constant.clear();
			blocks.add(new ArrayList<>());
			blocks.get(0).add(toIntMatrix(dilithium, dilithium.expandMatrix(info.publicKey2.getRho())));
			blocks.get(0).add(toIntMatrix(dilithium, dilithium.getFreeModule().matrixAlgebra().one()));
			rhs.addAll(toIntVector(dilithium, info.uncompressedPublicKey2).asList());
			for (int i = 0; i < dilithium.getRows() + dilithium.getCols(); i++) {
				for (int j = 0; j < dilithium.getDegree(); j++) {
					int c = dilithium.getEta();
					List<PFE> row = new ArrayList<>();
					for (int k = 0; k < bits; k++) {
						boolean xor = ((1 << k) & (i < dilithium.getCols() ? info.xorS[i][j]
								: info.xorE[i - dilithium.getRows()][j])) != 0;
						row.add(fp.getElement((xor ? 1 : -1) * (1 << k)));
						if (xor) {
							c -= 1 << k;
						}
					}
					constant.add(fp.getElement(c));
					diagonalBlocks.add(Matrix.fromRows(Collections.singletonList(new Vector<>(row))));
				}
			}
			blockMatrix = Matrix.fromBlockMatrix(blocks);
			adjustment = Matrix.fromBlockDiagonalMatrix(diagonalBlocks, fp.zero());
			adjustedBlocks.add(Collections.singletonList(mm.multiply(blockMatrix, adjustment)));
			rhsAdjustment.addAll(mm.multiply(blockMatrix, new Vector<>(constant)).asList());
		}
		Matrix<PFE> full = Matrix.fromBlockMatrix(adjustedBlocks);
		Vector<PFE> rhsVector = new Vector<>(rhs);
		rhsVector = mm.codomain().subtract(rhsVector, new Vector<>(rhsAdjustment));
		Vector<PFE> solved = mm.solve(full, rhsVector);
		System.out.println(solved);
		List<NTT<PFE>> s = new ArrayList<>();
		List<NTT<PFE>> e = new ArrayList<>();
		for (int i = 0; i < dilithium.getCols(); i++) {
			List<PFE> coeffs = new ArrayList<>();
			for (int j = 0; j < dilithium.getDegree(); j++) {
				int coeff = 0;
				for (int k = 0; k < bits; k++) {
					coeff += (1 << k) * solved.get(i * dilithium.getDegree() * bits + j * bits + k + 1).getValue()
							.intValueExact();
				}
				coeffs.add(fp.getInteger(dilithium.getEta() - coeff));
			}
			s.add(dilithium.getNTTRing().fromPolynomial(fp.getUnivariatePolynomialRing().getPolynomial(coeffs)));
		}
		for (int i = 0; i < dilithium.getRows(); i++) {
			List<PFE> coeffs = new ArrayList<>();
			for (int j = 0; j < dilithium.getDegree(); j++) {
				int coeff = 0;
				for (int k = 0; k < bits; k++) {
					coeff += (1 << k) * solved.get(dilithium.getCols() * dilithium.getDegree() * bits
							+ i * dilithium.getDegree() * bits + j * bits + k + 1).getValue().intValueExact();
				}
				coeffs.add(fp.getInteger(dilithium.getEta() - coeff));
			}
			e.add(dilithium.getNTTRing().fromPolynomial(fp.getUnivariatePolynomialRing().getPolynomial(coeffs)));
		}
		return dilithium.createPrivateKeyFromParameters(xors.get(0).publicKey1.getRho(), new byte[32], new Vector<>(s),
				new Vector<>(e));
	}

	@Test
	void testXorKeyRecovery() {
		// 29.615.817 s
		PrimeField fp = PrimeField.getPrimeField(Dilithium.PRIME);
		Dilithium dilithium = Dilithium.DILITHIUM2; /*
													 * new Dilithium(Dilithium.PRIME, 8, fp.power(fp.getElement(1753),
													 * 256 / 8).getValue().intValueExact(), 4, 4, 2, 13, 39, 17,
													 * (Dilithium.PRIME - 1) / 88, 80, 4, Sha3.SHAKE_128,
													 * Sha3.SHAKE_256, Sha3.SHAKE_256);
													 */
		List<XorKeyInfo> xors = new ArrayList<>();
		DilithiumPrivateKey key = dilithium.createPrivateKey();
		int bits = (int) Math.ceil(Math.log(2 * dilithium.getEta() + 1) / Math.log(2));
		for (int i = 0; i < 2 * bits; i++) {
			DilithiumPrivateKey key2 = dilithium.createPrivateKey();
			xors.add(new XorKeyInfo(dilithium, key, key2));
		}
		DilithiumPrivateKey recovered = xorKeyRecovery(dilithium, xors);
		System.out.println(key.getS());
		System.out.println(recovered.getS());
		System.out.println(key.getE());
		System.out.println(recovered.getE());
		// DilithiumSignature signature = dilithium.sign("This is a test
		// message".getBytes(), recovered);
		// assertTrue(dilithium.verify("This is a test message".getBytes(), signature,
		// dilithium.createPublicKey(key)));
	}
}
