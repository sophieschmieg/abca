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
		PrimeField f65537 = PrimeField.getPrimeField(65537);
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
		System.out.println(sp.matrixAlgebra().determinant( lat.generatorsAsMatrix()));
		System.out.println(sp.conditionNumber(lat.generatorsAsMatrix()));
			fail();
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
