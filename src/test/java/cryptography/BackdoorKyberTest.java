package cryptography;

import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

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
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.MatrixModule;
import fields.vectors.MatrixModule.SmithNormalFormResult;
import fields.vectors.RealIntegerLattice;
import fields.vectors.Vector;
import util.FunctionMathMap;

class BackdoorKyberTest {

	private Vector<IntE> toIntVector(Vector<NTT<FFE>> vector, Kyber kyber) {
		Integers z = Integers.z();
		List<IntE> intList = new ArrayList<>();
		for (int i = 0; i < kyber.getRank(); i++) {
			UnivariatePolynomial<PFE> polynomial = kyber.asPolynomial(vector.get(i + 1));
			for (int j = 0; j < kyber.getDegree(); j++) {
				intList.add(z.centeredLift(polynomial.univariateCoefficient(j), kyber.getPrime()));
			}
		}
		return new Vector<>(intList);
	}

	@Test
	void test() throws IOException {
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		int n = 8;
		int k = 2;
		int prime;
		int multiplier = 11;
		do {
			multiplier += 2;
			prime = multiplier * n + 1;
		} while (!z.isPrime(z.getInteger(multiplier)) || !z.isPrime(z.getInteger(prime)));
		System.out.println(prime);
		System.out.println(multiplier);
		// int prime = 137;
		PrimeField fp = PrimeField.getPrimeField(prime);
		int zeta = 0;
		do {
			zeta++;
		} while (!fp.power(fp.getInteger(zeta), n / 2).equals(fp.negative(fp.one())));
		System.out.println(zeta);
		Kyber kyber = new Kyber(n, prime, zeta, k, 1, 1, 4, 3, new TruncatedHashFunction(Sha2.SHA2_512, n / 4),
				new TruncatedHashFunction(Sha2.SHA2_256, n / 8), new Hmac(Sha2.SHA2_512),
				new Hkdf(new Hmac(Sha2.SHA2_512)), new Hkdf(new Hmac(Sha2.SHA2_512), new byte[0x01]));
		ByteArray privateKey = kyber.createPrivateKey();
		ByteArray publicKey = kyber.createPublicKey(privateKey);
		ParsedPublicKey parsedPublicKey = kyber.parsePublicKey(publicKey.array());
		ParsedPrivateKey parsedPrivateKey = kyber.parsePrivateKey(privateKey.array());
		UnivariatePolynomialRing<Fraction> rationalPolynomials = q.getUnivariatePolynomialRing();
		NumberField nf = NumberField.getNumberField(rationalPolynomials.parse("X^" + n + " + 1"));
		NumberFieldIntegers order = nf.maximalOrder();
		System.out.println(nf);
		final int p = prime;
		Matrix<NFE> matrix = Matrix
				.mapMatrix(new FunctionMathMap<>((NTT<FFE> t) -> nf.fromPolynomial(rationalPolynomials.getEmbedding(
						z.centeredLiftUnivariatePolynomial(kyber.asPolynomial(t), BigInteger.valueOf(p)),
						q.getEmbeddingMap()))), parsedPublicKey.getMatrix());
		MatrixModule<NFE> mm = matrix.getModule(order);
		matrix = mm.transpose(matrix);
//		Vector<NFE> t = Vector.mapVector(new FunctionMathMap<>((NTT<FFE> s) -> nf.fromPolynomial(rationalPolynomials
//				.getEmbedding(z.centeredLiftUnivariatePolynomial(kyber.asPolynomial(s), BigInteger.valueOf(prime)),
//						q.getEmbeddingMap()))),
//				parsed.getT());
		Vector<IntE> t = toIntVector(parsedPublicKey.getT(), kyber);
		Vector<IntE> s = toIntVector(parsedPrivateKey.getS(), kyber);
		FreeModule<IntE> intSpace = new FreeModule<>(z, n * k);
		List<List<Matrix<IntE>>> blocks = new ArrayList<>();
		for (int i = 0; i < k; i++) {
			List<Matrix<IntE>> row = new ArrayList<>();
			for (int j = 0; j < k; j++) {
				row.add(Matrix.mapMatrix(q.getAsIntegerMap(), nf.asMatrix(matrix.entry(i + 1, j + 1))));
			}
			blocks.add(row);
		}
		Matrix<IntE> intMatrix = Matrix.fromBlockMatrix(blocks);
		MatrixAlgebra<IntE> intAlgebra = intSpace.matrixAlgebra();
		Vector<IntE> e = Vector.mapVector(new FunctionMathMap<>((IntE a) -> z.centeredLift(fp.reduce(a), p)),
				intSpace.subtract(t, intAlgebra.multiply(intMatrix, s)));
		FreeModule<IntE> free = new FreeModule<>(z, 2 * n * k);
		MatrixAlgebra<IntE> m = free.matrixAlgebra();
		intMatrix = intAlgebra.transpose(intMatrix);
		List<Vector<IntE>> basis = new ArrayList<>();
		// A^t.s + e = t
		// A.x = y
		// t^t.x = (A^t.s + e)^t.x
		// = (s^t.A + e^t).x
		// = s^t.A.x + e^t.x
		// = s^t.y + e^t.x
		for (int i = 0; i < n * k; i++) {
			Vector<IntE> x = intSpace.getUnitVector(i + 1);
			Vector<IntE> y = intAlgebra.multiply(intMatrix, x);
			List<IntE> vectorList = new ArrayList<>();
			vectorList.addAll(x.asList());
			vectorList.addAll(y.asList());
			basis.add(new Vector<>(vectorList));
		}
		for (int i = 0; i < n * k; i++) {
			List<IntE> vectorList = new ArrayList<>();
			vectorList.addAll(intSpace.zero().asList());
			vectorList.addAll(intSpace.scalarMultiply(prime, intSpace.getUnitVector(i + 1)).asList());
			basis.add(new Vector<>(vectorList));
		}
		System.out.println(basis);
		Reals r = Reals.r(64);
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(r, 2 * n * k);
		RealIntegerLattice<Real, Vector<Real>> lattice = new RealIntegerLattice<>(space, free, basis);
		List<Vector<IntE>> originalBasis = basis;
		basis = lattice.getModuleGenerators();
		System.out.println(basis);
		for (Vector<IntE> basisVector : originalBasis) {
			System.out.println(lattice.asVector(basisVector));
		}
		List<IntE> results = new ArrayList<>();
		for (int i = 0; i < 2 * n * k; i++) {
			Vector<IntE> firstHalf = new Vector<>(basis.get(i).asList().subList(0, n * k));
			IntE product = z.centeredLift(fp.reduce(intSpace.innerProduct(t, firstHalf)), z.getInteger(prime));
			results.add(product);
		}
		Vector<IntE> resultVector = new Vector<>(results);
		System.out.println(resultVector);
		Matrix<IntE> reducedMatrix = Matrix.fromColumns(basis);
		Matrix<Real> realMatrix = Matrix.mapMatrix(new FunctionMathMap<>((IntE a) -> r.getInteger(a)),
				m.transpose(reducedMatrix));
		Vector<Real> realResults = Vector.mapVector(new FunctionMathMap<>((IntE a) -> r.getInteger(a)), resultVector);
		MatrixModule<Real> realModule = realMatrix.getModule(r);
		Vector<Real> solved = realModule.solve(realMatrix, realResults);
		List<IntE> reconstructedErrorList = new ArrayList<>();
		List<IntE> reconstructedSecretList = new ArrayList<>();
		for (int i = 0; i < n * k; i++) {
			reconstructedErrorList.add(solved.get(i + 1).round());
		}
		for (int i = 0; i < n * k; i++) {
			reconstructedSecretList.add(solved.get(n * k + i + 1).round());
		}
		Vector<IntE> reconstructedSecret = new Vector<>(reconstructedSecretList);
		Vector<IntE> reconstructedError = new Vector<>(reconstructedErrorList);
		System.out.println();
		System.out.println(s);
		System.out.println(reconstructedSecret);
		System.out.println();
		System.out.println(e);
		System.out.println(reconstructedError);
	}

}
