package fields.numberfields;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Complex;
import fields.floatingpoint.Reals;
import fields.helper.FieldEmbedding;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.MathMap;
import fields.interfaces.Ring.ChineseRemainderPreparation;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.IdealClassGroup.IdealClass;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.vectors.Vector;

class RelativeNumberFieldTest {

	// @Test
	void test5() {
		Rationals q = Rationals.q();
		NumberField nf = new NumberField(
				q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(5), q.zero(), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		UnivariatePolynomialRing<NFE> polynomials = nf.getUnivariatePolynomialRing();
		List<UnivariatePolynomial<NFE>> options = new ArrayList<>();
		/// options.add(polynomials.getPolynomial(nf.getInteger(-5), nf.zero(),
		/// nf.one()));
		// options.add(polynomials.getPolynomial(nf.getInteger(-2), nf.zero(),
		// nf.one()));
		options.add(polynomials.getPolynomial(nf.getInteger(1), nf.zero(), nf.one()));
		for (UnivariatePolynomial<NFE> option : options) {
			System.out.println();
			System.out.println(option);
			FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding = nf.getEmbeddedExtension(option);
			NumberFieldIntegers upperOrder = fieldEmbedding.getField().maximalOrder();
			System.out.println(order.primitiveRootsOfUnity());
			System.out.println(upperOrder.primitiveRootsOfUnity());
			System.out.println(fieldEmbedding.getField() + " integral basis: " + upperOrder.getModuleGenerators());
			for (UnivariatePolynomial<NFE> otherOption : options) {
				UnivariatePolynomial<NFE> embedded = fieldEmbedding.getField().getUnivariatePolynomialRing()
						.getEmbedding(otherOption, fieldEmbedding.getEmbeddingMap());
				System.out.println(embedded + " = " + fieldEmbedding.getField().factorization(embedded));
			}
			NumberFieldIdeal upperDifferent = upperOrder.different();
			IntE upperDiscriminant = upperOrder.discriminant();
			NumberFieldIdeal lowerDifferent = order.different();
			IntE lowerDiscriminant = order.discriminant();
			NumberFieldIdeal differentOver = upperOrder.differentOver(fieldEmbedding);
			NumberFieldIdeal discriminantOver = upperOrder.discriminantOver(fieldEmbedding);
			System.out.println("Different(" + fieldEmbedding.getField() + ") = " + upperDifferent + " = "
					+ upperOrder.idealFactorization(upperDifferent));
			System.out.println("Discriminant(" + fieldEmbedding.getField() + ") = " + upperDiscriminant + " = "
					+ Integers.z().uniqueFactorization(upperDiscriminant));
			System.out.println("Norm (different/" + fieldEmbedding.getField() + ") = " + upperDifferent.norm());
			System.out.println(
					"Different(" + nf + ") = " + lowerDifferent + " = " + order.idealFactorization(lowerDifferent));
			System.out.println("Discriminant(" + nf + ") = " + lowerDiscriminant + " = "
					+ Integers.z().uniqueFactorization(lowerDiscriminant));
			System.out.println("Different(" + fieldEmbedding + ") = " + differentOver + " = "
					+ upperOrder.idealFactorization(differentOver));
			System.out.println("Discriminant(" + fieldEmbedding + ") = " + discriminantOver + " = "
					+ order.idealFactorization(discriminantOver));
			NumberFieldIdeal ideal = order.idealsOver(2).get(0);
			NumberFieldIdeal ideal2 = order.idealsOver(3).get(0);
			System.out.println(
					ideal + "*" + fieldEmbedding.getField() + " = " + upperOrder.extend(ideal, fieldEmbedding));
			System.out.println(
					ideal2 + "*" + fieldEmbedding.getField() + " = " + upperOrder.extend(ideal2, fieldEmbedding));
			NumberFieldIdeal ideal3 = order.idealsOver(7).get(0);
			System.out.println(ideal3);
			NumberFieldIdeal extended = upperOrder.extend(ideal3, fieldEmbedding);
			System.out.println(extended);
			NFE g = extended.uniformizer();
			System.out.println(fieldEmbedding.asVector(g));
			System.out.println(
					fieldEmbedding.getField().factorization(fieldEmbedding.getField().getUnivariatePolynomialRing()
							.getPolynomial(upperOrder.one(), upperOrder.zero(), upperOrder.one())));
			System.out.println(
					fieldEmbedding.getField().factorization(fieldEmbedding.getField().getUnivariatePolynomialRing()
							.getPolynomial(upperOrder.getInteger(5), upperOrder.zero(), upperOrder.one())));
		}
	}

	// @Test
	void test23() {
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		NumberField nf = new NumberField(
				q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(23), q.zero(), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		nf.idealClassGroup();
		List<IntE> primes = new ArrayList<>();
		primes.add(z.getInteger(2));
		primes.add(z.getInteger(3));
		primes.add(z.getInteger(23));
		ChineseRemainderPreparation<IntE> prepare = z.prepareChineseRemainderTheoremModuli(primes);
		UnivariatePolynomialRing<IntE> integerPolynomials = z.getUnivariatePolynomialRing();
		List<UnivariatePolynomial<IntE>> polynomialList = new ArrayList<>();
		for (IntE prime : primes) {
			FiniteField fq = FiniteField.getFiniteField(z.power(prime, 3).getValue());
			polynomialList.add(integerPolynomials.getEmbedding(fq.minimalPolynomial(), new MathMap<>() {
				@Override
				public IntE evaluate(PFE t) {
					return fq.getBaseField().liftToInteger(t);
				}
			}));
		}
		List<IntE> c = new ArrayList<>();
		for (int i = 0; i <= 3; i++) {
			List<IntE> coefficients = new ArrayList<>();
			for (UnivariatePolynomial<IntE> p : polynomialList) {
				coefficients.add(p.univariateCoefficient(i));
			}
			c.add(z.chineseRemainderTheorem(coefficients, prepare));
		}
		UnivariatePolynomial<IntE> result = integerPolynomials.getPolynomial(c);
		IdealClass ic = null;
		for (IdealClass ict : nf.idealClassGroup()) {
			System.out.println(ict);
			if (ic == null && !ict.isPrincipal()) {
				ic = ict;
			}
		}
		NumberFieldIdeal generator = ic.representative();
		NumberFieldIdeal generator2 = order.multiply(generator, generator);
		NumberFieldIdeal generator3 = order.multiply(generator2, generator);
		System.out.println(generator);
		System.out.println(generator2);
		System.out.println(generator3);
		UnivariatePolynomialRing<NFE> polynomials = nf.getUnivariatePolynomialRing();
		FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding = nf.getEmbeddedExtension(
				polynomials.getPolynomial(nf.negative(generator3.generators().get(0)), nf.zero(), nf.zero(), nf.one()));
		NumberFieldIntegers upperOrder = fieldEmbedding.getField().maximalOrder();
		NumberFieldIdeal upperDifferent = upperOrder.different();
		IntE upperDiscriminant = upperOrder.discriminant();
		NumberFieldIdeal lowerDifferent = order.different();
		IntE lowerDiscriminant = order.discriminant();
		NumberFieldIdeal differentOver = upperOrder.differentOver(fieldEmbedding);
		NumberFieldIdeal discriminantOver = upperOrder.discriminantOver(fieldEmbedding);
		System.out.println("Different(" + fieldEmbedding.getField() + ") = " + upperDifferent + " = "
				+ upperOrder.idealFactorization(upperDifferent));
		System.out.println("Discriminant(" + fieldEmbedding.getField() + ") = " + upperDiscriminant + " = "
				+ Integers.z().uniqueFactorization(upperDiscriminant));
		System.out.println(
				"Different(" + nf + ") = " + lowerDifferent + " = " + order.idealFactorization(lowerDifferent));
		System.out.println("Discriminant(" + nf + ") = " + lowerDiscriminant + " = "
				+ Integers.z().uniqueFactorization(lowerDiscriminant));
		System.out.println("Different(" + fieldEmbedding + ") = " + differentOver + " = "
				+ upperOrder.idealFactorization(differentOver));
		System.out.println("Discriminant(" + fieldEmbedding + ") = " + discriminantOver + " = "
				+ order.idealFactorization(discriminantOver));
		NumberFieldIdeal ideal = order.idealsOver(2).get(0);
		NumberFieldIdeal ideal2 = order.idealsOver(3).get(0);
		System.out.println(ideal + "*" + fieldEmbedding.getField() + " = " + upperOrder.extend(ideal, fieldEmbedding));
		System.out
				.println(ideal2 + "*" + fieldEmbedding.getField() + " = " + upperOrder.extend(ideal2, fieldEmbedding));
		System.out.println(order.multiply(ideal, ideal) + "*" + fieldEmbedding.getField() + " = "
				+ upperOrder.extend(ideal, fieldEmbedding));
		System.out.println(order.multiply(ideal2, ideal2) + "*" + fieldEmbedding.getField() + " = "
				+ upperOrder.extend(ideal2, fieldEmbedding));
		System.out.println(upperOrder.idealsOver(ideal2, fieldEmbedding));
	}

	//@Test
	void test17() {
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		NumberField nf = new NumberField(
				q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(17), q.zero(), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		for (IdealClass ic : nf.idealClassGroup()) {
			System.out.println(nf.idealClassGroup().getOrder(ic));
		}
	}
	
	//@Test
	void testReal5() {
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		NumberField nf = new NumberField(
				q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(-5), q.zero(), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		for (IdealClass ic : nf.idealClassGroup()) {
			System.out.println(nf.idealClassGroup().getOrder(ic));
		}
		System.out.println(order.freeUnitGroupGenerators());
			for (NFE g : order.getModuleGenerators()) {
			System.out.println(g);
			System.out.println(nf.minkowskiEmbedding(g));
			System.out.println(nf.logRepresentation(g));
			System.out.println();
				}
	}

	//@Test
	void testReal3() {
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		NumberField nf = new NumberField(
				q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(-3), q.zero(), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		for (IdealClass ic : nf.idealClassGroup()) {
			System.out.println(nf.idealClassGroup().getOrder(ic));
		}
		System.out.println(order.freeUnitGroupGenerators());
			for (NFE g : order.getModuleGenerators()) {
			System.out.println(g);
			System.out.println(nf.minkowskiEmbedding(g));
			System.out.println(nf.logRepresentation(g));
			System.out.println();
				}
	}
	
	//@Test
	void testDirichletZetaFunction () {
		NumberField nf = new NumberField();
		System.out.println(nf.dirichletZetaFunction(10).evaluate(Complex.c(128).getEmbedding(Reals.r(128).getFraction(Rationals.q().getFraction(1025, 1024)))));
	}

	//@Test
	void testReal7() {
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		NumberField nf = new NumberField(
				q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(-7), q.zero(), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		Reals r = (Reals)nf.minkowskiEmbeddingSpace().getField();
		nf.minkowskiEmbeddingSpace().latticePointsInParallelotope(new Vector<>(r.getInteger(10), r.getInteger(1)), order);
		for (IdealClass ic : nf.idealClassGroup()) {
			System.out.println(nf.idealClassGroup().getOrder(ic));
		}
		System.out.println(order.freeUnitGroupGenerators());
			for (NFE g : order.getModuleGenerators()) {
			System.out.println(g);
			System.out.println(nf.minkowskiEmbedding(g));
			System.out.println(nf.logRepresentation(g));
			System.out.println();
				}
	}
	
	// @Test
	void test21() {
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		NumberField nf = new NumberField(
				q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(21), q.zero(), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		nf.idealClassGroup();
	}
}
