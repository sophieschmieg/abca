package cryptography;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

import cryptography.Kyber.ParsedKem;
import cryptography.Kyber.ParsedPrivateKey;
import cryptography.Kyber.ParsedPublicKey;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.NTTRing.NTT;
import fields.integers.Integers;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.numberfields.NumberFieldIntegers.SmallestNumberFieldIntegerSolutionPreparation;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.RealLattice;
import fields.vectors.Vector;
import util.FunctionMathMap;
import util.Pair;

class KyberTest {

	@Test
	void testVectorTest() {
		Kyber kyber = new Kyber(2, 3, 2, 10, 4, Sha3.SHA3_512, Sha3.SHA3_256, Sha3.SHAKE_256_PRF, Sha3.SHAKE_128,
				Sha3.SHAKE_256);
		ByteArray publicKey = ByteArray.fromString(
				"bb46cb622988f1d70ebbdb760f498fece49012b73d71e49866ec780cc83005cc252da512e3c5b43997b29470c6827a729e6acced8a214832915bfb629b83a1f477c835c29ce0e175847a63a1a69870c4adc264663b896b0e7b984871c2c8802f7380cd3c754c8371c33ec82d2d7518ec2638326aca9e0a9bec2abf68ca70390bc2dd2aa691a502f9c0ba9a829ba89925b69a13a6ba5d2390117809887b51b7fcc6c5af31c464cc17e81221813a837562a15a0ac9c81c289689890ccb6635f15d23bbb4ec248a3ebc86bc9b283b5a3460d3a06ca326b519091095836eca075fd6392402955da84fe7538cc18763a91a2df2790039d464b7750250d9b93e774cc656bddcf3c17c7b8820db883e8331a0003de25825bd0176cb3c4051a384b687309421acd0bb8bd2ec6a422a5d6509958eb25a5bc1ba3ee43b5453c53ad7776021b19e2a5b1f553b36134cf3421951cc5645d74813565f3d8b6e83d3607ef5751d09ce8fba792865231d60c82cd78ccd26001ed1b3b0378c87131b0cd3a389b963cc01b8b4255a31522a4c83b8706a7c1981b555730804836a730b45cee7853f771869306896d8c94fd261d057b48a17327fd8634eb77823239f623ac3d4d38a1b1b0d36736aa7f28250912ff5a710ed4776be53c264d687a9a1bd3c0b99453742bc4554fba93841ea2fdae6214f49090a24b93cc287e1c322a1298bfeb21313fa4367b84d928c1d6c1923991341429b3d6a59c3d0721f474c753c3093809bc8894b63f747b2d34a219923467464698e595131e3606a5812ef08351ef931679bb5a36a6fe547a103401c7634a7790283f4db6d15bac3cf81389cc23b74302d7b61389198ca7241870fc351ef159b894c50f93cb24c945bbb70306f1a7a0c5063d642433f852b3ef6c227d2b49c9a4c72ab21aae567c9bc6d93817a6161c12d6200d6a08afc1325a7da9a3782b8191aa5d6f30751d01c6af308b2c06c92084e4086c4f3e306e9d5cd82fc151a700bb066c4f97191a46833b8ac55fce41200772450e826b8e5985d830610fb44f216aae0e889de6857db4a0d4b8479a1800f92fba629a6cf40d7000000000000000000000000000000000000000000000000000000");
		kyber.parsePublicKey(publicKey.array());
		publicKey = ByteArray.fromString(
				"55533942a315d6053fedc19de4f70c2cf7378e8b3ab28605d0f2c7e6f5892f4365e8eb087ec0c00d879ce1c87655e9b7b5292374f8bce4a75422d64d7efac0e8e62df4d952c972727d3c3199627d6d012e1c22a09a1498d423b2d9cc84fe168484b4c8f53198ba39322ae9c20f30601b59949b846708025fade360f8c98617bba23e777e51bb20f23c095bf06ce4b7aa9ce03c80a63f48e83fbe783f718355cdf67b019807eb129932e441db7a5210b85d92b370676685b66220e1cc5f084a55bf073d9562370cfc70250c4c9e2150eeb942b11c09450a42a7e708968786d359af48327693010a21f79106734f30f3220d38afc9b12fe316bd818c2cbeb3be1cf8af8aa119f6e5b939a561634875028ccca426382635273ee3c724d63db832cacf9acee749cd1bb101c5ebaf7688a37595a56f73125be4ab734b5c9b451dfe764f88d524bae85dad3c21573632c2043482fa19e80740275031f62804062535c1b691304696b305674cb9370fcb71c14b45b301ce8793c766f282071123babb74fd63ccadba777168bf2be438a96668d5e3be4304833233bebb282dc01085a906c0662676f41c1d8a5c1fc43b71979183c33b4a2f80b059d9be802647c987167dbc0a3e40067bd4793d093556c53f2684656d52af3e3109d61144cabc8cc124b5be20281db0c9f89b2b54a4933194162f1c62cd453f5bf26b42c85f572b2a961c23b5a1bdba638bb3268575b0cdeed6076a7931dd587df2eb06a73648c085843c94641ed2790b1ca64eda335a6a91b5311282f2a5e86282ef7986d667ce7c9c335dc733c890c0b7234253255876fba9b511a3989766427aacabb83f9bc9468f07323aba6c0588c6e1653db6c7c9b5eaa41d5262e367a2f3d47cfdb404027a1719f92f45735237a4c2e4215b47a59b605a35477887f80091f2a8a17a466251e12ff071a68a485171d5446b5b925fd4a12a7c70eda625933515c97c3cde8436ec461e77c247e29a4a4e114281313e7e86bcd5acc5201328418901f49bacb0b572e6bb5703e312b0eba18337bf34f9b09be2c77a403017c3a434071940eb1047fb8d23c5ab9094000000000000000000000000000000000000000000000000000000");
		Kyber.ParsedPublicKey parsed = kyber.parsePublicKey(publicKey.array());
		System.out.println(parsed.getMatrix());
		publicKey = ByteArray.fromString(
				"1fbca0f5a93566588d65880f44d4c1ec8a16e5e508f7726107761cbf026b89c580a62abc82137d9b8b9edbf082980416c43444072c9008e2184b557dbc1a79598c5d2c95c2bd7c12aee084c150b90bfa832859c841b235eb377ebe533469660ecf784faf862c6b026feb36508130c3c09196a146423823a210b4250e2b3135f1686a4c65d3f878f6233a3d48264b3302858774be5b328fd801164b21fd3c7e547c65d83894580536ea919fe00609a8f018f5aace9534c3d98c2508841e448440750747d3e6ba2907928129843825b428f650e8917e26787883276ad1116ccbb551de12755a4b62dccc5ff6ac180672735e354dd1d018cad55bc997bb93ba31d74c49dc99cdbbe359d834a65eb65c044b2b77754cd380abb6442d67440adbca3716a36b321766c6f40051617e716a49475ba2995b73b7fa45d4878517c45358ec0510e849a181be28c2000b376783c238c6c80300d25e4e1724783c182a668f9164c8e31112c2723f88087576a4c58e4a8d345c88fa4674ab9a2c3c0ba849847a64cb00d4d3b512758ea9d3892cd014370940cfcc4fad27c4ca09a9cf9a2fecdc11eeb93d88361914762c95e844d80529f586259ca750d72269b9a28336668855460222bb27af073260297910f1250336169de469e93093596405123398c861cb081c6833a9af25911739e9c1ab1a2817548009ec8e6654521e8b6efa3abf96c9354b7002380ab76f251d7a0436d69345ad46cafbab6edaa66123f2389397af5226c99401291b9722ee634d286619319bb5ee1773f799bbe54b6ec04aa6dde72570e8983e244b8bb1a5b0ecae0ae23618fb307c322d21017a241808cea14c2ad2cd9a474a11a2ba15cb99e29464e724322763914b1cb3a2859e716292fac42367120cefc7a9e66a63ef891ebf21739286083e194baea6bc4f84932443adfcc37eef5452cab5c3e67425fd751189e87cdfc3220e041bc1b5399fb28c5da2093993c66b8c5736190446d698be94abad3a7ce14b11454a8b1f2421c2022a5cf00a3329582fccb104c46f99ecad91d1c6caa2845f4202ca5520ba32376f93bce1e8128518b5e578000000000000000000000000000000000000000000000000000000");
		parsed = kyber.parsePublicKey(publicKey.array());
		FreeModule<NTT<FFE>> free = new FreeModule<>(kyber.getNTTRing(), 2);
	System.out.println(free.matrixAlgebra().multiply(parsed.getMatrix(), 	free.matrixAlgebra().inverse(parsed.getMatrix())));
	}

	@Test
	void kyberTest() {
		Kyber kyber = new Kyber(2, 3, 2, 10, 4, Sha3.SHA3_512, Sha3.SHA3_256, Sha3.SHAKE_256_PRF, Sha3.SHAKE_128,
				Sha3.SHAKE_256);
		ByteArray privateKey = kyber.createPrivateKey();
		ByteArray publicKey = kyber.createPublicKey(privateKey);
		System.out.println(privateKey);
		System.out.println(publicKey);
		Pair<VariableLengthKey, ByteArray> encapsulated = kyber.encapsulate(publicKey);
		System.out.println(encapsulated.getSecond());
		byte[] aliceKey = encapsulated.getFirst().key(32);
		System.out.println(new ByteArray(aliceKey));
		VariableLengthKey decrypted = kyber.decapsulate(encapsulated.getSecond(), privateKey);
		byte[] bobKey = decrypted.key(32);
		System.out.println(new ByteArray(bobKey));
		assertEquals(new ByteArray(aliceKey), new ByteArray(bobKey));
	}

	@Test
	void funSizeKyberTest() throws IOException {
		/*PrimeField f65537 = PrimeField.getPrimeField(65537);
		PFE fiveByThree = f65537.divide(f65537.getInteger(5), f65537.getInteger(3));
		System.out.println(fiveByThree);
		System.out.println(Integers.z().lift(fiveByThree));
		Reals r = Reals.r(128);
		List<Vector<Real>> cs = new ArrayList<>();
		cs.add(new Vector<>(r.getInteger(Integers.z().lift(fiveByThree)), r.one()));
		cs.add(new Vector<>(r.getInteger(65537), r.zero()));
		Matrix<Real> m = Matrix.fromColumns(cs);
		Vector<Real> rhsExact = new Vector<>(r.getInteger(5), r.getInteger(3));
		Vector<Real> rhsOff = new Vector<>(r.getDouble(5.8), r.getDouble(2.4));
		MatrixModule<Real> mm = m.getModule(r);
		Vector<Real> solvedExact = mm.solve(m, rhsExact);
		Vector<Real> solvedOff = mm.solve(m, rhsOff);
		System.out.println(rhsExact);
		System.out.println(rhsOff);
		System.out.println(solvedExact);
		System.out.println(solvedOff);
		Vector<Real> approxExact = mm.multiply(m,
				Vector.mapVector(new FunctionMathMap<>((Real x) -> r.getInteger(x.round())), solvedExact));
		Vector<Real> approxOff = mm.multiply(m,
				Vector.mapVector(new FunctionMathMap<>((Real x) -> r.getInteger(x.round())), solvedOff));
		System.out.println(approxExact);
		System.out.println(approxOff);
		FiniteRealVectorSpace sp = new FiniteRealVectorSpace(r, 2);
		System.out.println(sp.matrixAlgebra().determinant(m));
		System.out.println(sp.conditionNumber(m));
		RealLattice lat = new RealLattice(sp, cs);
		System.out.println(lat.getModuleGenerators());
		System.out.println(sp.matrixAlgebra().determinant(lat.generatorsAsMatrix()));
		System.out.println(sp.conditionNumber(lat.generatorsAsMatrix()));
		fail();*/
		int degree = 8;
		int prime = 41;
		int k = 2;
		Kyber kyber = new Kyber(degree, prime, 3, k, 1, 1, 4, 3, new TruncatedHashFunction(Sha2.SHA2_512, 2),
				new TruncatedHashFunction(Sha2.SHA2_256, 1), new Hmac(Sha2.SHA2_512), new Hkdf(new Hmac(Sha2.SHA2_512)),
				new Hkdf(new Hmac(Sha2.SHA2_512), new byte[0x01]));
		ByteArray privateKey = kyber.createPrivateKey();
		ByteArray publicKey = kyber.createPublicKey(privateKey);
		System.out.println("Private Key: " + privateKey);
		System.out.println("Public Key:  " + publicKey);
		Pair<VariableLengthKey, ByteArray> encapsulated = kyber.encapsulate(publicKey);
		System.out.println("Key Encapsulation Message: " + encapsulated.getSecond());
		byte[] aliceKey = encapsulated.getFirst().key(1);
		System.out.println("Alice's shared key: " + new ByteArray(aliceKey));
		VariableLengthKey decrypted = kyber.decapsulate(encapsulated.getSecond(), privateKey);
		byte[] bobKey = decrypted.key(1);
		System.out.println("Bob's shared key:   " + new ByteArray(bobKey));
		assertEquals(new ByteArray(aliceKey), new ByteArray(bobKey));
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		UnivariatePolynomialRing<Fraction> rationalPolynomials = q.getUnivariatePolynomialRing();
		NumberField nf = NumberField.getNumberField(rationalPolynomials.parse("X^" + degree + " + 1"));
		NumberFieldIntegers order = nf.maximalOrder();
		System.out.println(nf);
		System.out.println(order.idealsOver(prime));
		FreeModule<NFE> space = new FreeModule<>(order, k);
		ParsedPublicKey parsedPublicKey = kyber.parsePublicKey(publicKey.array());
		Matrix<NFE> matrix = Matrix
				.mapMatrix(new FunctionMathMap<>((NTT<FFE> t) -> nf.fromPolynomial(rationalPolynomials.getEmbedding(
						z.centeredLiftUnivariatePolynomial(kyber.asPolynomial(t), BigInteger.valueOf(prime)),
						q.getEmbeddingMap()))), parsedPublicKey.getMatrix());
		matrix = space.matrixAlgebra().transpose(matrix);
		Vector<NFE> t = Vector.mapVector(new FunctionMathMap<>((NTT<FFE> s) -> nf.fromPolynomial(rationalPolynomials
				.getEmbedding(z.centeredLiftUnivariatePolynomial(kyber.asPolynomial(s), BigInteger.valueOf(prime)),
						q.getEmbeddingMap()))),
				parsedPublicKey.getT());
		ParsedKem parsedKem = kyber.parseKem(encapsulated.getSecond().array());
		Vector<NFE> u = Vector.mapVector(
				new FunctionMathMap<>((Polynomial<PFE> s) -> nf.fromPolynomial(rationalPolynomials.getEmbedding(
						z.centeredLiftUnivariatePolynomial(s, BigInteger.valueOf(prime)), q.getEmbeddingMap()))),
				parsedKem.getU());
		NFE v = nf.fromPolynomial(rationalPolynomials.getEmbedding(
				z.centeredLiftUnivariatePolynomial(parsedKem.getV(), BigInteger.valueOf(prime)), q.getEmbeddingMap()));
		ParsedPrivateKey parsedPrivateKey = kyber.parsePrivateKey(privateKey.array());
		Vector<NFE> s = Vector.mapVector(new FunctionMathMap<>((NTT<FFE> g) -> nf.fromPolynomial(rationalPolynomials
				.getEmbedding(z.centeredLiftUnivariatePolynomial(kyber.asPolynomial(g), BigInteger.valueOf(prime)),
						q.getEmbeddingMap()))),
				parsedPrivateKey.getS());
		Vector<NFE> e = Vector
				.mapVector(
						new FunctionMathMap<NFE, NFE>(
								(NFE x) -> nf
										.fromPolynomial(rationalPolynomials.getEmbedding(
												z.centeredLiftUnivariatePolynomial(
														z.reduceUnivariatePolynomial(
																z.getUnivariatePolynomialRing().getEmbedding(
																		x.asPolynomial(), q.getAsIntegerMap()),
																z.getInteger(prime)),
														BigInteger.valueOf(prime)),
												q.getEmbeddingMap()))),
						space.subtract(t, space.matrixAlgebra().multiply(matrix, s)));
		System.out.println("s: " + s);
		System.out.println("e: " + e);
		MathMap<Vector<NFE>, Vector<Real>> asVector = new FunctionMathMap<>((Vector<NFE> x) -> {
			List<Real> result = new ArrayList<>();
			for (NFE coeff : x.asList()) {
				result.addAll(order.embedding(coeff).asList());
			}
			return new Vector<>(result);
		});
		FiniteRealVectorSpace realSpace = new FiniteRealVectorSpace(Reals.r(128), k * degree);
		System.out.println("Size s: " + realSpace.valueNorm(asVector.evaluate(s)));
		System.out.println("Size e: " + realSpace.valueNorm(asVector.evaluate(e)));
		NumberFieldIdeal ideal = order.getIdeal(Collections.singletonList(order.getInteger(prime)));
		MathMap<NFE, NFE> idealReduction = new FunctionMathMap<>((NFE a) -> ideal.residue(a));
		Vector<NFE> expected = space.add(space.matrixAlgebra().multiply(matrix, s), e);
		assertEquals(Vector.mapVector(idealReduction, t), Vector.mapVector(idealReduction, expected));
		List<Vector<NFE>> columns = new ArrayList<>();
		columns.addAll(matrix.asColumnList());
		for (int i = 0; i < k; i++) {
			columns.add(space.getUnitVector(i + 1));
		}
		SmallestNumberFieldIntegerSolutionPreparation preparation = order.prepareSmallestIntegerSolution(columns,
				ideal);
		System.out.println("Preparation complete!");
		Vector<NFE> result = order.smallestIntegerSolution(t, preparation);
		System.out.println(result);
		Vector<NFE> reconstructedS = new Vector<>(result.asList().subList(0, k));
		Vector<NFE> reconstructedE = new Vector<>(result.asList().subList(k, 2 * k));
		Vector<NFE> reconstructed = space.add(space.matrixAlgebra().multiply(matrix, reconstructedS), reconstructedE);
		System.out.println("Size reconstructed s: " + realSpace.valueNorm(asVector.evaluate(reconstructedS)));
		System.out.println("Size reconstructed e: " + realSpace.valueNorm(asVector.evaluate(reconstructedE)));
		assertEquals(Vector.mapVector(idealReduction, t), Vector.mapVector(idealReduction, reconstructed));
	}
}
