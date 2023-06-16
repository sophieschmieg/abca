package fields.helper;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

import java.io.IOException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.helper.GenericAlgebraicRingExtension.GenericAlgebraicExtensionElement;
import fields.helper.NTTRing.NTT;
import fields.interfaces.Ring.QuotientAndRemainderResult;
import fields.interfaces.UnivariatePolynomialRing;

class NTTRingTest {

	@Test
	void test() {
		PrimeField fp = PrimeField.getPrimeField(17);
		PFE root = fp.primitiveRoot();
		NTTRing<PFE, PrimeField> ntt = new NTTRing<>(fp, 8, root);
		System.out.println(ntt.fromPolynomial(fp.getUnivariatePolynomialRing().one()));
		System.out.println(ntt.fromPolynomial(fp.getUnivariatePolynomialRing().getVar()));
		GenericAlgebraicRingExtension<PFE> ext = new GenericAlgebraicRingExtension<>(fp.getUnivariatePolynomialRing()
				.add(fp.getUnivariatePolynomialRing().getVarPower(8), fp.getUnivariatePolynomialRing().one()), fp);
		for (int tc = 0; tc < 10; tc++) {
			GenericAlgebraicExtensionElement<PFE> ext1 = ext.getRandomElement();
			NTT<PFE> ntt1 = ntt.fromPolynomial(ext1.asPolynomial());
			NTT<PFE> ntt2 = ntt.getRandomElement();
			GenericAlgebraicExtensionElement<PFE> ext2 = ext.fromPolynomial(ntt.asPolynomial(ntt2));
			assertEquals(ext1, ext.fromPolynomial(ntt.asPolynomial(ntt1)));
			assertEquals(ntt2, ntt.fromPolynomial(ext2.asPolynomial()));
			assertEquals(ext.add(ext1, ext2), ext.fromPolynomial(ntt.asPolynomial(ntt.add(ntt1, ntt2))));
			assertEquals(ext.subtract(ext1, ext2), ext.fromPolynomial(ntt.asPolynomial(ntt.subtract(ntt1, ntt2))));
			assertEquals(ext.multiply(ext1, ext2), ext.fromPolynomial(ntt.asPolynomial(ntt.multiply(ntt1, ntt2))));
			assertEquals(ext.isDivisible(ext1, ext2), ntt.isDivisible(ntt1, ntt2));
			QuotientAndRemainderResult<NTT<PFE>> qr = ntt.quotientAndRemainder(ntt1, ntt2);
			assertEquals(ntt1, ntt.add(ntt.multiply(qr.getQuotient(), ntt2), qr.getRemainder()));
		}
	}

	@Test
	void testKyberNtt() throws IOException {
		int numbers = 1 << 20;
		int eta = 2;
		Random rand = new SecureRandom();
		PrimeField fp = PrimeField.getPrimeField(3329);
		int degree = 256;
//		PFE root = fp.getElement(17);
//		System.out.println(fp.getMultiplicativeGroup().getOrder(root));
		FiniteField fq = FiniteField.getFiniteField(fp.getUnivariatePolynomialRing().parse("X^2 + -17"), fp);
//		System.out.println(fq.getMultiplicativeGroup().getOrder(fq.getEmbedding(root)));
//		System.out.println(fq.getMultiplicativeGroup().getOrder(fq.alpha()));
		NTTRing<FFE, FiniteField> ntt = new NTTRing<>(fq, degree, fq.alpha());
		UnivariatePolynomialRing<FFE> polynomials = fq.getUnivariatePolynomialRing();
		Map<PFE, Integer> counts = new TreeMap<>();
		for (int i = 0; i < numbers; i++) {
			List<FFE> coeffs = new ArrayList<>();
			for (int j = 0; j < degree; j++) {
				int coeff = 0;
				for (int k = 0; k < eta; k++) {
					if (rand.nextBoolean()) {
						coeff++;
					}
				}
				for (int k = 0; k < eta; k++) {
					if (rand.nextBoolean()) {
						coeff--;
					}
				}
				coeffs.add(fq.getInteger(coeff));
			}
			NTT<FFE> element = ntt.fromPolynomial(polynomials.getPolynomial(coeffs));
			PFE first = ntt.asList(element).get(0).asPolynomial().univariateCoefficient(0);
			counts.put(first, counts.getOrDefault(first, 0) + 1);
		}
		int sum = 0;
		for (PFE fieldElement : fp) {
			int actual = counts.getOrDefault(fieldElement, 0);
			sum += actual * actual;
			// System.out.println(fieldElement + ": " + counts.getOrDefault(fieldElement,
			// 0));
		}
		double testValue = (double) sum * 3329.0 / (double) numbers - (double) numbers;
		if (testValue < 3521.76) {
			System.out.println(testValue + " < 3521.76");
		} else {
			System.err.println(testValue + " >= 3521.76");
			fail();
		}
//		NTT<FFE> one = ntt.fromPolynomial(polynomials.one());
//		System.out.println(one);
//		assertEquals(polynomials.one(), ntt.asPolynomial(one));
//		NTT<FFE> var = ntt.fromPolynomial(polynomials.getVar());
//		System.out.println(var);
//		System.out.println(ntt.power(var, degree));
//		assertEquals(polynomials.getVar(), ntt.asPolynomial(var));
//		NTT<FFE> var2 = ntt.fromPolynomial(polynomials.getVarPower(2));
//		System.out.println(var2);
//		System.out.println(ntt.asPolynomial(var2));
//		int[] values = new int[256];
//		UnivariatePolynomial<FFE> testPolynomial = polynomials.zero();
//		for (int i = 0; i < degree; i++) {
//			testPolynomial = polynomials.add(testPolynomial, polynomials.getEmbedding(fq.getInteger(values[i]), i));
//		}
//		System.out.println(testPolynomial);
//		UnivariatePolynomial<FFE> testPolynomial2 = polynomials.zero();
//		for (int i = 0; i < degree; i++) {
//			testPolynomial2 = polynomials.add(testPolynomial2, polynomials.getEmbedding(fq.getInteger(3329 - i), i));
//		}
//		System.out.println(testPolynomial2);
//		NTT<FFE> nttTestPolynomial = ntt.fromPolynomial(testPolynomial);
//		NTT<FFE> nttTestPolynomial2 = ntt.fromPolynomial(testPolynomial2);
//		NTT<FFE> productNtt = ntt.multiply(nttTestPolynomial, nttTestPolynomial2);
//		UnivariatePolynomial<FFE> product = ntt.asPolynomial(productNtt);
//		System.out.println(product);
	}

	@Test
	void testKyberLikeNtt() throws IOException {
		PrimeField fp = PrimeField.getPrimeField(17);
		int degree = 16;
		PFE root = fp.primitiveRoot();
		System.out.println(fp.getMultiplicativeGroup().getOrder(root));
		FiniteField fq = FiniteField.getFiniteField(
				fp.getUnivariatePolynomialRing().subtract(fp.getUnivariatePolynomialRing().getVarPower(2),
						fp.getUnivariatePolynomialRing().getEmbedding(root)),
				fp);
		System.out.println(fq.getMultiplicativeGroup().getOrder(fq.getEmbedding(root)));
		System.out.println(fq.getMultiplicativeGroup().getOrder(fq.alpha()));
		NTTRing<FFE, FiniteField> ntt = new NTTRing<>(fq, degree, fq.alpha());
		GenericAlgebraicRingExtension<PFE> ext = new GenericAlgebraicRingExtension<>(fp.getUnivariatePolynomialRing()
				.add(fp.getUnivariatePolynomialRing().getVarPower(16), fp.getUnivariatePolynomialRing().one()), fp);
		UnivariatePolynomialRing<FFE> polynomials = fq.getUnivariatePolynomialRing();
		NTT<FFE> one = ntt.fromPolynomial(polynomials.one());
		System.out.println(one);
		assertEquals(polynomials.one(), ntt.asPolynomial(one));
		NTT<FFE> var = ntt.fromPolynomial(polynomials.getVar());
		System.out.println(var);
		System.out.println(ntt.power(var, degree));
		assertEquals(polynomials.getVar(), ntt.asPolynomial(var));
		for (int tc = 0; tc < 10; tc++) {
			GenericAlgebraicExtensionElement<PFE> ext1 = ext.getRandomElement();
			NTT<FFE> ntt1 = ntt.fromPolynomial(polynomials.getEmbedding(ext1.asPolynomial(), fq.getEmbeddingMap()));
			GenericAlgebraicExtensionElement<PFE> ext2 = ext.getRandomElement();
			NTT<FFE> ntt2 = ntt.fromPolynomial(polynomials.getEmbedding(ext2.asPolynomial(), fq.getEmbeddingMap()));
			assertEquals(ext1, ext.fromPolynomial(
					fp.getUnivariatePolynomialRing().getEmbedding(ntt.asPolynomial(ntt1), fq.asBaseFieldElementMap())));
			assertEquals(ext2, ext.fromPolynomial(
					fp.getUnivariatePolynomialRing().getEmbedding(ntt.asPolynomial(ntt2), fq.asBaseFieldElementMap())));
			assertEquals(ext.add(ext1, ext2), ext.fromPolynomial(fp.getUnivariatePolynomialRing()
					.getEmbedding(ntt.asPolynomial(ntt.add(ntt1, ntt2)), fq.asBaseFieldElementMap())));
			assertEquals(ext.subtract(ext1, ext2), ext.fromPolynomial(fp.getUnivariatePolynomialRing()
					.getEmbedding(ntt.asPolynomial(ntt.subtract(ntt1, ntt2)), fq.asBaseFieldElementMap())));
			assertEquals(ext.multiply(ext1, ext2), ext.fromPolynomial(fp.getUnivariatePolynomialRing()
					.getEmbedding(ntt.asPolynomial(ntt.multiply(ntt1, ntt2)), fq.asBaseFieldElementMap())));
			// assertEquals(ext.isDivisible(ext1, ext2), ntt.isDivisible(ntt1, ntt2));
			QuotientAndRemainderResult<NTT<FFE>> qr = ntt.quotientAndRemainder(ntt1, ntt2);
			assertEquals(ntt1, ntt.add(ntt.multiply(qr.getQuotient(), ntt2), qr.getRemainder()));
		}
	}
}
